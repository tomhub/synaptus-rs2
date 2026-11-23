//! **A. NUFFT-PSM**: the same discrete sum as [`crate::psm::psm`], evaluated with a
//! type-1 NUFFT along depth instead of an O(N_ω·N_z) loop.
//!
//! For fixed (k⊥, layer l) PSM computes I(j) = Σ_m a_m e^{i j x_m} with
//! x_m = k_z(ω_m)·Δz (dossier v2 §5.1). With j = k + N/2 this is a type-1 sum with
//! strengths a_m e^{i(N/2)x_m}. No interpolation, no Jacobian, no kz grid: the
//! ω-samples are used where they are. Cost per column O(N_ω·w + N log N).
//!
//! Kernels: cancellation-free (factored) kz, as `PsmOptions::high_precision_kernels`.
//! The result equals v1 PSM up to the NUFFT tolerance (test N4).

use crate::error::{invalid, Result};
use crate::fftutil::{self, Fft2, C64};
use crate::grid;
use crate::kernels::{cis, KzFormula};
use crate::nufft::Nufft1d;
use crate::psm::{build_spectrum, kz_table, layer_depths, spatial_fft_sizes, DepthStep};
use crate::types::{validate, Acquisition, FocusedImage, LayerImage, Medium, PulseEchoData};
use rayon::prelude::*;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct PsmNufftOptions {
    pub depth_step: DepthStep,
    /// Time FFT zero-padding multiplier (finer ω grid): same meaning as v1.
    pub t_fft_mult: usize,
    pub x_fft_mult: usize,
    pub y_fft_mult: usize,
    /// NUFFT relative tolerance.
    pub tol: f64,
}

impl Default for PsmNufftOptions {
    fn default() -> Self {
        Self { depth_step: DepthStep::HighFrequency, t_fft_mult: 1, x_fft_mult: 1, y_fft_mult: 1, tol: 1e-12 }
    }
}

pub(crate) fn even_at_least(n: usize) -> usize {
    (n + 1) & !1
}

/// Shared tail: per-layer columns `cols[ik][j]` (nz each) → cropped spatial image.
pub(crate) fn columns_to_layer(cols: &[C64], nz: usize, nfx: usize, nfy: usize, data: &PulseEchoData, z: Vec<f64>) -> LayerImage {
    let nk = nfx * nfy;
    let fft2_inv = Fft2::new(nfx, nfy, true);
    let planes: Vec<Vec<C64>> = (0..nz)
        .into_par_iter()
        .map(|j| {
            let mut plane: Vec<C64> = (0..nk).map(|ik| cols[ik * nz + j]).collect();
            let mut out = Vec::with_capacity(data.nx * data.ny);
            fftutil::inverse_plane_crop(&fft2_inv, &mut plane, data.nx, data.ny, &mut out);
            out
        })
        .collect();
    LayerImage { nz, nx: data.nx, ny: data.ny, z, data: planes.concat() }
}

pub fn psm_nufft(data: &PulseEchoData, acq: &Acquisition, medium: &Medium, opts: &PsmNufftOptions) -> Result<FocusedImage> {
    let z_offset = validate(data, acq, medium)?;
    let nfft_t = grid::fft_len(opts.t_fft_mult, data.nt as f64)?;
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    let (wl, wh) = (2.0 * PI * acq.f_low, 2.0 * PI * acq.f_high);
    let keep: Vec<usize> = (0..nfft_t).filter(|&k| omega_full[k] >= wl && omega_full[k] <= wh).collect();
    if keep.is_empty() {
        return invalid("no FFT bins inside [f_low, f_high]");
    }
    let mut spec = build_spectrum(data, acq, nfft_t, &keep, None, nfx, nfy);
    let (nw, nk) = (spec.nw(), spec.nk());
    let mut layers = Vec::with_capacity(medium.n_layers());

    for l in 0..medium.n_layers() {
        let c = medium.velocities[l];
        let dz = match opts.depth_step {
            DepthStep::HighFrequency => (c / 2.0) / acq.f_high,
            DepthStep::Bandwidth => (c / 2.0) / (acq.f_high - acq.f_low),
        };
        let kz = kz_table(&spec, c, KzFormula::Factored);
        // cumulative evanescent mask, as v1
        spec.w.par_chunks_mut(nw).zip(kz.par_chunks(nw)).for_each(|(col, kzc)| {
            for (v, k) in col.iter_mut().zip(kzc) {
                if k.is_nan() {
                    *v = C64::new(0.0, 0.0);
                }
            }
        });
        let kz: Vec<f64> = kz.into_iter().map(|k| if k.is_nan() { 0.0 } else { k }).collect();
        let depths = layer_depths(medium, l, dz, z_offset);
        let nz = depths.len();
        let n_modes = even_at_least(nz.max(2));
        let plan = Nufft1d::new(n_modes, opts.tol)?;
        let start = if l == 0 { z_offset } else { 0.0 };
        let half = (n_modes / 2) as f64;

        let mut cols = vec![C64::new(0.0, 0.0); nk * nz];
        cols.par_chunks_mut(nz).enumerate().for_each(|(ik, out)| {
            let base = ik * nw;
            let mut x = Vec::with_capacity(nw);
            let mut a = Vec::with_capacity(nw);
            for iw in 0..nw {
                let v = spec.w[base + iw];
                if v == C64::new(0.0, 0.0) {
                    continue; // masked (evanescent) or empty: contributes nothing
                }
                let k = kz[base + iw];
                let xm = k * dz;
                x.push(xm);
                a.push(v * cis(k * start + half * xm));
            }
            if x.is_empty() {
                return;
            }
            let f = plan.type1(&x, &a);
            out.copy_from_slice(&f[..nz]);
        });
        layers.push(columns_to_layer(&cols, nz, nfx, nfy, data, depths));

        if l + 1 < medium.n_layers() {
            let d = medium.thicknesses[l];
            spec.w.par_iter_mut().zip(kz.par_iter()).for_each(|(v, &k)| *v *= cis(k * d));
        }
    }
    let x = (0..data.nx).map(|i| i as f64 * acq.x_step).collect();
    let y = (0..data.ny).map(|i| i as f64 * acq.y_step).collect();
    Ok(FocusedImage { layers, x, y })
}
