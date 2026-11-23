//! Multilayer omega-k (Stolt) migration as specified in the thesis
//! (Skjelvareid 2012, eqs. 4.15–4.17 and §4.A.3, eqs. 4.25–4.29).
//!
//! Differences from `matlab/core/mulok.m` (each can be re-enabled through
//! [`MulokOptions`] for attribution studies):
//! * Jacobian A = ĉ·kz/√(kz²+k⊥²) (thesis 4.16) instead of kz²/(kz²+k⊥²).
//! * kz grid starts at 0 (thesis 4.27) with dz = ĉ/f_high (same depth grid as PSM),
//!   instead of starting at 2πf_low/ĉ with dz = ĉ/(f_high − f_low).
//! * Frequency oversampling by zero-padding in time (exact DFT samples) instead
//!   of `fft(ifft(band), up·n)` with a mislabelled ω axis.
//! * Evanescent components are masked when migrating between interfaces (4.A.2).
//! * Normalisation: the image approximates the PSM image (the ω-sum at the base
//!   frequency spacing), so PSM and MULOK are directly comparable.
//!
//! Accuracy options beyond the thesis: higher-order interpolation (cubic Lagrange,
//! Lanczos-windowed sinc) and *demodulated* interpolation (interpolate
//! W(ω)·e^{iωt_c}, which varies slowly in ω, then restore e^{−iωt_c} exactly).

use crate::error::{invalid, Error, Result};
use crate::fftutil::{Fft2, C64};
use crate::grid;
use crate::kernels::{cis, kz_two_way, KzFormula};
use crate::psm::{build_spectrum, layer_depths, spatial_fft_sizes};
use crate::types::{validate, Acquisition, FocusedImage, LayerImage, Medium, PulseEchoData};
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Interp {
    Linear,
    /// 4-point Lagrange (local cubic) interpolation.
    Cubic,
    /// Lanczos-windowed sinc with half-width `a` samples (2a-point stencil).
    Lanczos(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Jacobian {
    /// ĉ·kz/√(kz²+k⊥²) — thesis eqs. 3.25, 4.16, 4.26.
    Thesis,
    /// ĉ·kz²/(kz²+k⊥²) — the factor used in `mulok.m`/Python (times ĉ, so that only
    /// the k-dependence differs from `Thesis`).
    LegacySquared,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KzSupport {
    /// kz ∈ [0, kz_max) — thesis eq. 4.27.
    Full,
    /// Zero all kz < 2πf_low/ĉ, emulating the truncation of `mulok.m`'s kz grid.
    TruncatedBelowKLow,
}

#[derive(Debug, Clone)]
pub struct MulokOptions {
    pub interp: Interp,
    /// Frequency oversampling factor via zero-padding in time (≥ 1).
    pub upsample: usize,
    pub jacobian: Jacobian,
    pub kz_support: KzSupport,
    /// Shift of the interpolation query, in units of the base bin spacing
    /// 2πfs/nfft_base. 0 = correct. `mulok.m` with upSamp=4 corresponds to −0.75.
    pub omega_shift_bins: f64,
    /// Interpolate the demodulated spectrum W·e^{iωt_c} (t_c = centre of the
    /// remaining record at vertical incidence) instead of W itself.
    pub demodulate: bool,
    /// kz formula for interface migration.
    pub kz: KzFormula,
    pub t_fft_mult: usize,
    pub x_fft_mult: usize,
    pub y_fft_mult: usize,
    pub z_fft_mult: usize,
    pub hh: Option<Vec<f64>>,
}

impl MulokOptions {
    /// Thesis algorithm with the interpolation used in `mulok.m` (linear, 4× oversampling).
    pub fn thesis() -> Self {
        Self {
            interp: Interp::Linear,
            upsample: 4,
            jacobian: Jacobian::Thesis,
            kz_support: KzSupport::Full,
            omega_shift_bins: 0.0,
            demodulate: false,
            kz: KzFormula::Factored,
            t_fft_mult: 1,
            x_fft_mult: 1,
            y_fft_mult: 1,
            z_fft_mult: 1,
            hh: None,
        }
    }
    /// Thesis algorithm with the accuracy improvements: demodulated cubic
    /// interpolation, 4x frequency oversampling and a 4x longer depth period.
    /// Chosen from the measurements in `runs/*-accuracy-mulok*`, `*-mulok-floor`
    /// (see ACCURACY.md). The longer depth period removes wrap-around of band-edge
    /// ringing, which the thesis criterion (eq. 4.28) does not cover.
    pub fn accurate() -> Self {
        Self { interp: Interp::Cubic, upsample: 4, demodulate: true, z_fft_mult: 4, ..Self::thesis() }
    }
}

impl Default for MulokOptions {
    fn default() -> Self {
        Self::accurate()
    }
}

const GUARD: usize = 8; // extra bins on each side of the band for interpolation stencils

fn interp_at(y: &[C64], p: f64, m: Interp) -> Option<C64> {
    let n = y.len() as isize;
    let i0 = p.floor() as isize;
    let t = p - i0 as f64;
    let get = |i: isize| if i >= 0 && i < n { Some(y[i as usize]) } else { None };
    match m {
        Interp::Linear => Some(get(i0)? * (1.0 - t) + get(i0 + 1)? * t),
        Interp::Cubic => {
            let (ym, y0, y1, y2) = (get(i0 - 1)?, get(i0)?, get(i0 + 1)?, get(i0 + 2)?);
            let wm = -t * (t - 1.0) * (t - 2.0) / 6.0;
            let w0 = (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0;
            let w1 = -(t + 1.0) * t * (t - 2.0) / 2.0;
            let w2 = (t + 1.0) * t * (t - 1.0) / 6.0;
            Some(ym * wm + y0 * w0 + y1 * w1 + y2 * w2)
        }
        Interp::Lanczos(a) => {
            let a = a as isize;
            let mut s = C64::new(0.0, 0.0);
            for i in (i0 - a + 1)..=(i0 + a) {
                let x = p - i as f64;
                let w = if x == 0.0 {
                    1.0
                } else {
                    let px = PI * x;
                    (a as f64) * px.sin() * (px / a as f64).sin() / (px * px)
                };
                s += get(i)? * w;
            }
            Some(s)
        }
    }
}

pub fn mulok(data: &PulseEchoData, acq: &Acquisition, medium: &Medium, opts: &MulokOptions) -> Result<FocusedImage> {
    let z_offset = validate(data, acq, medium)?;
    if opts.upsample == 0 || opts.z_fft_mult == 0 {
        return invalid("upsample and z_fft_mult must be >= 1");
    }
    if let Interp::Lanczos(a) = opts.interp {
        if a == 0 || a > GUARD {
            return invalid(format!("Lanczos half-width must be in 1..={GUARD}"));
        }
    }
    let nh = opts.hh.as_ref().map_or(1, |h| h.len().max(1));
    let n_rec = data.nt + nh - 1;
    let nfft_base = grid::fft_len(opts.t_fft_mult, n_rec as f64)?;
    let nfft_t = nfft_base.checked_mul(opts.upsample).ok_or_else(|| Error::SizeOverflow("nfft_t".into()))?;
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let d_omega = (2.0 * PI * acq.fs) / nfft_t as f64;
    let d_omega_base = (2.0 * PI * acq.fs) / nfft_base as f64;
    let (wl, wh) = (2.0 * PI * acq.f_low, 2.0 * PI * acq.f_high);
    // Stored bins: band plus guard, all non-negative frequencies below Nyquist.
    let k_lo = ((wl / d_omega).ceil() as usize).saturating_sub(GUARD);
    let k_hi = (((wh / d_omega).floor() as usize) + GUARD).min(nfft_t / 2);
    let keep: Vec<usize> = (k_lo..=k_hi).collect();
    let mut spec = build_spectrum(data, acq, nfft_t, &keep, opts.hh.as_deref(), nfx, nfy);
    let nw = spec.nw();
    let nk = spec.nk();
    let t_rec = n_rec as f64 / acq.fs;
    let shift = opts.omega_shift_bins * d_omega_base;

    let mut layers = Vec::with_capacity(medium.n_layers());
    for l in 0..medium.n_layers() {
        let c = medium.velocities[l];
        let ch = c / 2.0;
        let dz = ch / acq.f_high;
        let depths = layer_depths(medium, l, dz, z_offset);
        // z period must exceed the remaining record (thesis 4.28) and the layer.
        let period = if l == 0 { ch * t_rec } else { ch * (acq.t_delay + t_rec) };
        let need = period.max(depths.len() as f64 * dz) / dz;
        let n_kz = grid::fft_len(opts.z_fft_mult, need)?;
        let d_kz = 2.0 * PI / (n_kz as f64 * dz);
        let k_low = wl / ch;
        // Two-way travel time to this layer's top at vertical incidence.
        let t_top: f64 = (0..l).map(|m| 2.0 * medium.thicknesses[m] / medium.velocities[m]).sum();
        let t_c = if opts.demodulate { acq.t_delay + 0.5 * t_rec - t_top } else { 0.0 };
        let scale = d_kz / d_omega_base;

        let mut pkz = vec![C64::new(0.0, 0.0); nk * n_kz];
        let spec_ref = &spec;
        pkz.par_chunks_mut(n_kz).enumerate().for_each(|(ik, out)| {
            let (kx, ky) = spec_ref.k_of(ik);
            let kp = kx.hypot(ky);
            let col = &spec_ref.w[ik * nw..(ik + 1) * nw];
            // demodulated column: q(ω_i) = W(ω_i)·e^{iω_i t_c}
            let q: Vec<C64> = if t_c != 0.0 {
                col.iter().zip(&spec_ref.omega).map(|(v, &w)| v * cis(w * t_c)).collect()
            } else {
                col.to_vec()
            };
            for (j, o) in out.iter_mut().enumerate() {
                let kz = j as f64 * d_kz;
                if opts.kz_support == KzSupport::TruncatedBelowKLow && kz < k_low {
                    continue;
                }
                let kk = kz.hypot(kp);
                let om = ch * kk;
                if !(om >= wl && om <= wh) || kk == 0.0 {
                    continue;
                }
                let a = match opts.jacobian {
                    Jacobian::Thesis => ch * kz / kk,
                    Jacobian::LegacySquared => ch * (kz / kk) * (kz / kk),
                };
                let wq = om + shift;
                let p = wq / d_omega - k_lo as f64;
                let v = match interp_at(&q, p, opts.interp) {
                    Some(v) => v,
                    None => continue, // outside stored range; cannot happen inside band (GUARD)
                };
                let v = if t_c != 0.0 { v * cis(-wq * t_c) } else { v };
                let mut v = v * (a * scale);
                if l == 0 {
                    v *= cis(kz * z_offset);
                }
                *o = v;
            }
        });
        // Σ_k P(kz_k) e^{i kz_k z_j}, z_j = j·dz  (unnormalised inverse DFT).
        let inv_z = FftPlanner::new().plan_fft_inverse(n_kz);
        pkz.par_chunks_mut(n_kz).for_each(|col| inv_z.process(col));
        let fft2_inv = Fft2::new(nfx, nfy, true);
        let mut img = Vec::with_capacity(depths.len() * data.nx * data.ny);
        let mut plane = vec![C64::new(0.0, 0.0); nk];
        for iz in 0..depths.len() {
            for ik in 0..nk {
                plane[ik] = pkz[ik * n_kz + iz];
            }
            crate::fftutil::inverse_plane_crop(&fft2_inv, &mut plane, data.nx, data.ny, &mut img);
        }
        layers.push(LayerImage { nz: depths.len(), nx: data.nx, ny: data.ny, z: depths, data: img });

        if l + 1 < medium.n_layers() {
            let d = medium.thicknesses[l];
            let omega = spec.omega.clone();
            let (kxv, kyv) = (spec.kx.clone(), spec.ky.clone());
            spec.w.par_chunks_mut(nw).enumerate().for_each(|(ik, col)| {
                let (kx, ky) = (kxv[ik / nfy], kyv[ik % nfy]);
                for (iw, v) in col.iter_mut().enumerate() {
                    *v = match kz_two_way(opts.kz, omega[iw], kx, ky, c) {
                        Some(kz) => *v * cis(kz * d),
                        None => C64::new(0.0, 0.0),
                    };
                }
            });
        }
    }
    let x = (0..data.nx).map(|i| i as f64 * acq.x_step).collect();
    let y = (0..data.ny).map(|i| i as f64 * acq.y_step).collect();
    Ok(FocusedImage { layers, x, y })
}
