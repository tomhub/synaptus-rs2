//! Phase shift migration (PSM), single- and multi-layer, 2D and 3D.
//!
//! Algorithm: thesis (Skjelvareid 2012) eqs. 3.17–3.21, 4.12, 4.14 and §4.A.4;
//! structure identical to `matlab/core/psm.m`. The numerical kernels are
//! selectable through [`PsmOptions`]; [`PsmOptions::legacy`] reproduces `psm.m`'s
//! arithmetic (MATLAB kz formula, recursive phase shift, plain summation).

use crate::error::Result;
use crate::fftutil::{self, Fft2, C64};
use crate::grid;
use crate::kernels::{cis, kz_two_way, sum_naive, sum_neumaier, KzFormula, Summation};
use crate::types::{validate, Acquisition, FocusedImage, LayerImage, Medium, PulseEchoData};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Depth sampling step per layer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepthStep {
    /// dz = (c/2)/f_high — `psm.m`; Nyquist for the complex image (dossier §5.4).
    HighFrequency,
    /// dz = (c/2)/(f_high − f_low) — Python `PhaseShiftMigration` and `mulok.m`.
    /// Under-samples steep-angle content; pixel values themselves are unaffected.
    Bandwidth,
}

/// How the wavefield is advanced from one image depth to the next.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseUpdate {
    /// Multiply by e^{i kz dz} repeatedly (thesis §3.5.3; `psm.m`). Rounding error
    /// in the accumulated phase/magnitude grows with the number of steps.
    Recursive,
    /// Evaluate e^{i kz (j·dz)} directly for every depth j. No accumulation;
    /// one sin/cos per element per depth.
    Direct,
    /// Recursive, but re-anchored by a direct evaluation every `n` steps; the
    /// accumulated error is bounded by ~n steps' worth.
    Reanchored(usize),
}

#[derive(Debug, Clone)]
pub struct PsmOptions {
    pub kz: KzFormula,
    pub phase: PhaseUpdate,
    pub sum: Summation,
    pub depth_step: DepthStep,
    /// FFT size multipliers (MATLAB `tFftMult`, `xFftMult`, `yFftMult`).
    pub t_fft_mult: usize,
    pub x_fft_mult: usize,
    pub y_fft_mult: usize,
    /// Two-way impulse response for matched filtering (MATLAB `hh`); `None` = no filter.
    pub hh: Option<Vec<f64>>,
}

impl PsmOptions {
    /// Arithmetic of the original `psm.m`.
    pub fn legacy() -> Self {
        Self {
            kz: KzFormula::Legacy,
            phase: PhaseUpdate::Recursive,
            sum: Summation::Naive,
            depth_step: DepthStep::HighFrequency,
            t_fft_mult: 1,
            x_fft_mult: 1,
            y_fft_mult: 1,
            hh: None,
        }
    }
    /// Cancellation-free kz, direct phase evaluation and compensated summation.
    /// NOTE: the rounding study (T5, `runs/*-psm-rounding`) found **no measurable
    /// image-level difference** from [`PsmOptions::legacy`] in double precision
    /// (all variants within 4e-13 of a 40-digit reference over 842 depth steps),
    /// while `Direct` costs ~50% more time. Kept for experimentation only.
    pub fn high_precision_kernels() -> Self {
        Self { kz: KzFormula::Factored, phase: PhaseUpdate::Direct, sum: Summation::Neumaier, ..Self::legacy() }
    }
}

impl Default for PsmOptions {
    /// `psm.m` arithmetic — measured to be accurate to ~1e-13 relative in double.
    fn default() -> Self {
        Self::legacy()
    }
}

/// Band-limited (ω, kx, ky) spectrum of the data, layout `[ik][iw]`.
pub(crate) struct Spectrum {
    pub omega: Vec<f64>,
    pub kx: Vec<f64>,
    pub ky: Vec<f64>,
    pub nfx: usize,
    pub nfy: usize,
    pub w: Vec<C64>,
}

impl Spectrum {
    pub fn nw(&self) -> usize {
        self.omega.len()
    }
    pub fn nk(&self) -> usize {
        self.nfx * self.nfy
    }
    pub fn k_of(&self, ik: usize) -> (f64, f64) {
        (self.kx[ik / self.nfy], self.ky[ik % self.nfy])
    }
}

/// FFT sizes shared by PSM and MULOK. `nfy == 1` for 2D data.
pub(crate) fn spatial_fft_sizes(data: &PulseEchoData, xm: usize, ym: usize) -> Result<(usize, usize)> {
    let nfx = grid::fft_len(xm, data.nx as f64)?;
    let nfy = if data.is_3d() { grid::fft_len(ym, data.ny as f64)? } else { 1 };
    Ok((nfx, nfy))
}

/// Time FFT (length `nfft_t`), matched filter + t_delay compensation, selection of
/// bins `keep`, spatial FFT. Mirrors `psm.m`: the per-ω factor
/// conj(HH)·e^{-iω t_delay} is formed first and then multiplied onto the data.
pub(crate) fn build_spectrum(
    data: &PulseEchoData,
    acq: &Acquisition,
    nfft_t: usize,
    keep: &[usize],
    hh: Option<&[f64]>,
    nfx: usize,
    nfy: usize,
) -> Spectrum {
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    let hh_f = hh.map(|h| fftutil::fft_real_padded(h, nfft_t));
    let factor: Vec<C64> = keep
        .iter()
        .map(|&k| {
            let delay = C64::new(0.0, -omega_full[k] * acq.t_delay).exp();
            match &hh_f {
                Some(hf) => hf[k].conj() * delay,
                None => delay,
            }
        })
        .collect();
    let mut cols = fftutil::time_fft_select(data, nfft_t, keep);
    let nw = keep.len();
    cols.par_chunks_mut(nw).for_each(|c| {
        for (v, f) in c.iter_mut().zip(&factor) {
            *v *= f;
        }
    });
    let w = fftutil::spatial_fft(&cols, nw, data.nx, data.ny, nfx, nfy);
    let kx = grid::k_vec(nfx, acq.x_step);
    let ky = if nfy > 1 { grid::k_vec(nfy, acq.y_step) } else { vec![0.0] };
    Spectrum { omega: keep.iter().map(|&k| omega_full[k]).collect(), kx, ky, nfx, nfy, w }
}

/// Number of image planes and their depths for layer `l` (as `psm.m`/`mulok.m`):
/// layer 0 starts at the measurement z offset, others at their top interface;
/// `ceil(extent/dz) + 1` planes.
pub(crate) fn layer_depths(medium: &Medium, l: usize, dz: f64, z_offset: f64) -> Vec<f64> {
    let z_if = medium.interfaces();
    let (start, extent) = if l == 0 {
        (z_offset, medium.thicknesses[0] - z_offset)
    } else {
        (z_if[l], medium.thicknesses[l])
    };
    let n = (extent / dz).ceil() as usize + 1;
    (0..n).map(|j| j as f64 * dz + start).collect()
}

/// Phase shift migration. Returns one complex image per layer; depth step in
/// layer l is dz_l = (c_l/2)/f_high (as `psm.m`).
pub fn psm(data: &PulseEchoData, acq: &Acquisition, medium: &Medium, opts: &PsmOptions) -> Result<FocusedImage> {
    let z_offset = validate(data, acq, medium)?;
    let nh = opts.hh.as_ref().map_or(1, |h| h.len().max(1));
    let nfft_t = grid::fft_len(opts.t_fft_mult, (data.nt + nh - 1) as f64)?;
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    let (wl, wh) = (2.0 * PI * acq.f_low, 2.0 * PI * acq.f_high);
    let keep: Vec<usize> = (0..nfft_t).filter(|&k| omega_full[k] >= wl && omega_full[k] <= wh).collect();
    if keep.is_empty() {
        return crate::error::invalid("no FFT bins inside [f_low, f_high]");
    }
    let spec = build_spectrum(data, acq, nfft_t, &keep, opts.hh.as_deref(), nfx, nfy);
    migrate_layers(spec, data, acq, medium, z_offset, opts)
}

/// Diagnostics: the band-limited spectrum PSM migrates (after time FFT, t_delay
/// compensation, band selection and spatial FFT). Returns (ω, kx, ky, W) with W in
/// layout `[ik][iw]`, ik = ikx·nfy + iky. Used by the rounding study (T5).
#[doc(hidden)]
pub fn debug_spectrum(
    data: &PulseEchoData,
    acq: &Acquisition,
    medium: &Medium,
    opts: &PsmOptions,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<C64>)> {
    validate(data, acq, medium)?;
    let nh = opts.hh.as_ref().map_or(1, |h| h.len().max(1));
    let nfft_t = grid::fft_len(opts.t_fft_mult, (data.nt + nh - 1) as f64)?;
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    let (wl, wh) = (2.0 * PI * acq.f_low, 2.0 * PI * acq.f_high);
    let keep: Vec<usize> = (0..nfft_t).filter(|&k| omega_full[k] >= wl && omega_full[k] <= wh).collect();
    let s = build_spectrum(data, acq, nfft_t, &keep, opts.hh.as_deref(), nfx, nfy);
    Ok((s.omega, s.kx, s.ky, s.w))
}

pub(crate) fn migrate_layers(
    mut spec: Spectrum,
    data: &PulseEchoData,
    acq: &Acquisition,
    medium: &Medium,
    z_offset: f64,
    opts: &PsmOptions,
) -> Result<FocusedImage> {
    let (nw, nk) = (spec.nw(), spec.nk());
    let fft2_inv = Fft2::new(spec.nfx, spec.nfy, true);
    let mut layers = Vec::with_capacity(medium.n_layers());

    for l in 0..medium.n_layers() {
        let c = medium.velocities[l];
        let dz = match opts.depth_step {
            DepthStep::HighFrequency => (c / 2.0) / acq.f_high,
            DepthStep::Bandwidth => (c / 2.0) / (acq.f_high - acq.f_low),
        };
        // kz per (ik, iw); evanescent components are masked out of the reference
        // wavefield cumulatively (as psm.m) and get kz = 0.
        let kz = kz_table(&spec, c, opts.kz);
        spec.w.par_chunks_mut(nw).zip(kz.par_chunks(nw)).for_each(|(col, kzc)| {
            for (v, k) in col.iter_mut().zip(kzc) {
                if k.is_nan() {
                    *v = C64::new(0.0, 0.0);
                }
            }
        });
        let kz: Vec<f64> = kz.into_iter().map(|k| if k.is_nan() { 0.0 } else { k }).collect();

        let depths = layer_depths(medium, l, dz, z_offset);
        let start_shift = if l == 0 { z_offset } else { 0.0 };
        // Wavefield at the first image depth of this layer.
        let w0: Vec<C64> = spec
            .w
            .par_iter()
            .zip(kz.par_iter())
            .map(|(v, &k)| if l == 0 { v * cis(k * start_shift) } else { *v })
            .collect();
        let ph: Vec<C64> = match opts.phase {
            PhaseUpdate::Direct => Vec::new(),
            _ => kz.par_iter().map(|&k| cis(k * dz)).collect(),
        };
        let mut cur = w0.clone();
        let mut img = Vec::with_capacity(depths.len() * data.nx * data.ny);
        let mut plane = vec![C64::new(0.0, 0.0); nk];

        for j in 0..depths.len() {
            let use_direct = match opts.phase {
                PhaseUpdate::Direct => true,
                PhaseUpdate::Recursive => false,
                PhaseUpdate::Reanchored(n) => n > 0 && j % n == 0 && j > 0,
            };
            let zj = j as f64 * dz;
            plane
                .par_iter_mut()
                .zip(cur.par_chunks_mut(nw))
                .enumerate()
                .for_each(|(ik, (p, col))| {
                    let base = ik * nw;
                    if matches!(opts.phase, PhaseUpdate::Direct) {
                        let it = (0..nw).map(|iw| w0[base + iw] * cis(kz[base + iw] * zj));
                        *p = match opts.sum {
                            Summation::Naive => sum_naive(it),
                            Summation::Neumaier => sum_neumaier(it),
                        };
                        return;
                    }
                    if use_direct {
                        for iw in 0..nw {
                            col[iw] = w0[base + iw] * cis(kz[base + iw] * zj);
                        }
                    }
                    *p = match opts.sum {
                        Summation::Naive => sum_naive(col.iter().copied()),
                        Summation::Neumaier => sum_neumaier(col.iter().copied()),
                    };
                    for iw in 0..nw {
                        col[iw] *= ph[base + iw];
                    }
                });
            fftutil::inverse_plane_crop(&fft2_inv, &mut plane, data.nx, data.ny, &mut img);
        }
        layers.push(LayerImage { nz: depths.len(), nx: data.nx, ny: data.ny, z: depths, data: img });

        // Extrapolate the reference wavefield to the next interface (eq. 4.14).
        if l + 1 < medium.n_layers() {
            let d = medium.thicknesses[l];
            spec.w.par_iter_mut().zip(kz.par_iter()).for_each(|(v, &k)| *v *= cis(k * d));
        }
    }
    let x = (0..data.nx).map(|i| i as f64 * acq.x_step).collect();
    let y = (0..data.ny).map(|i| i as f64 * acq.y_step).collect();
    Ok(FocusedImage { layers, x, y })
}

/// Table `kz[ik*nw+iw]` of two-way kz for velocity `c`; NaN marks evanescent
/// components (kz² < 0). Internal sentinel only — never returned to callers.
pub(crate) fn kz_table(s: &Spectrum, c: f64, f: KzFormula) -> Vec<f64> {
    let nw = s.nw();
    let mut kz = vec![0.0f64; s.nk() * nw];
    kz.par_chunks_mut(nw).enumerate().for_each(|(ik, col)| {
        let (kx, ky) = s.k_of(ik);
        for (iw, k) in col.iter_mut().enumerate() {
            *k = kz_two_way(f, s.omega[iw], kx, ky, c).unwrap_or(f64::NAN);
        }
    });
    kz
}
