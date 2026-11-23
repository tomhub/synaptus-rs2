//! Cylindrical phase shift migration (CPSM), thesis (Skjelvareid 2012) §6.1,
//! eqs. 6.4–6.10. Geometry: transducer on a cylinder of radius r0 pointing
//! outward; data p(t, φ, z) (2D: p(t, φ)).
//!
//! [`cpsm`] is the corrected version and the default (TEST_PLAN amendment A3):
//! integer angular orders (φ zero-padded to the full circle), evanescent
//! components removed, no DC bin, Hankel-ratio transfer function. The literal
//! port of `cpsm.m`, defects included, is [`crate::legacy::cpsm`].

use crate::error::{invalid, Result};
use crate::fftutil::{self, Fft2, C64};
use crate::grid;
use crate::types::PulseEchoData;
use rayon::prelude::*;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransferFunction {
    /// Haun et al. (2002) approximation (thesis eq. 6.9); `cpsm.m`'s default. Its
    /// amplitude √(r/r0) is the reciprocal of the Hankel asymptote (dossier C8, open).
    Haun,
    /// Hankel ratio H⁽²⁾ₙ(k_r r)/H⁽²⁾ₙ(k_r r0), Gardner–Collin closed form. Amplitude
    /// consistent with thesis eq. 6.7. Default of [`cpsm`].
    Gardner,
}

#[derive(Debug, Clone, Copy)]
pub struct CylAcquisition {
    pub fs: f64,
    pub t_delay: f64,
    /// Lower band edge [Hz]; 0 is allowed (as in `cpsm.m`; includes the DC bin).
    pub f_low: f64,
    pub f_high: f64,
    /// Angular step [rad].
    pub phi_step: f64,
    /// Axial step [m] (ignored for 2D data).
    pub z_step: f64,
    /// Scanning radius [m].
    pub r0: f64,
    /// Wave velocity [m/s].
    pub c: f64,
}

/// CPSM image: `data[(ir * nphi + iphi) * nz + iz]`.
#[derive(Debug, Clone)]
pub struct CylImage {
    pub r: Vec<f64>,
    /// φ coordinates, centred: ((0:nφ−1) − (nφ−1)/2)·Δφ.
    pub phi: Vec<f64>,
    pub z: Vec<f64>,
    pub nphi: usize,
    pub nz: usize,
    pub data: Vec<C64>,
}

impl CylImage {
    pub fn at(&self, ir: usize, iphi: usize, iz: usize) -> C64 {
        self.data[(ir * self.nphi + iphi) * self.nz + iz]
    }
}

/// `approxHankel2` from `cpsm.m` (Gardner & Collin 2000); 0 where |z| <= 1.1|n|.
pub(crate) fn approx_hankel2(z: f64, n: f64) -> C64 {
    if !(z.abs() > 1.1 * n.abs()) {
        return C64::new(0.0, 0.0);
    }
    let sq = (z * z - n * n).sqrt();
    let amp = (2.0 / (PI * sq)).sqrt();
    let phase = sq - PI / 4.0 - n * (PI / 2.0) + n * (n / sq).atan();
    amp * C64::new(0.0, -phase).exp()
}

/// Range vector `r_start:r_step:r_end` (MATLAB colon semantics up to rounding of
/// the last element; pass an explicit vector to [`cpsm`] if exact agreement matters).
pub fn range_vector(r_start: f64, r_step: f64, r_end: f64) -> Vec<f64> {
    let n = ((r_end - r_start) / r_step + 1e-10).floor() as usize + 1;
    (0..n).map(|i| r_start + i as f64 * r_step).collect()
}

/// Options for the corrected [`cpsm`]. `Default` is the corrected algorithm.
#[derive(Debug, Clone, Copy)]
pub struct CpsmOptions {
    pub transfer: TransferFunction,
    /// Time FFT zero-padding multiplier.
    pub t_fft_mult: usize,
}

impl Default for CpsmOptions {
    fn default() -> Self {
        Self { transfer: TransferFunction::Gardner, t_fft_mult: 1 }
    }
}

/// Corrected CPSM (TEST_PLAN amendment A3). Differences from `cpsm.m`:
/// (i) φ is zero-padded to the full circle N = round(2π/Δφ), so the angular
/// orders are integers (thesis eq. 6.4) and the scanned arc is not wrapped onto
/// itself; (ii) evanescent components are removed instead of entering with zero
/// phase; (iii) the DC bin is excluded; (iv) Hankel-ratio transfer by default.
pub fn cpsm(data: &PulseEchoData, acq: &CylAcquisition, ranges: &[f64], opts: &CpsmOptions) -> Result<CylImage> {
    let pos = |v: f64| v.is_finite() && v > 0.0;
    if !pos(acq.fs) || !pos(acq.c) || !pos(acq.r0) || !pos(acq.phi_step) {
        return invalid("fs, c, r0 and phi_step must be finite and > 0");
    }
    if data.is_3d() && !pos(acq.z_step) {
        return invalid("z_step must be finite and > 0 for 3D data");
    }
    if !(acq.f_low.is_finite() && acq.f_low >= 0.0 && acq.f_high > acq.f_low && acq.f_high <= acq.fs / 2.0) {
        return invalid("need 0 <= f_low < f_high <= fs/2");
    }
    if !(acq.t_delay.is_finite() && acq.t_delay >= 0.0) {
        return invalid("t_delay must be finite and >= 0");
    }
    if ranges.is_empty() || ranges.iter().any(|r| !pos(*r)) {
        return invalid("ranges must be non-empty, finite and > 0");
    }
    let full = (2.0 * PI / acq.phi_step).round();
    if !(full >= data.nx as f64) {
        return invalid(format!("the scanned arc ({} × {} rad) exceeds the full circle", data.nx, acq.phi_step));
    }
    let nfn = full as usize;
    let nfft_t = grid::fft_len(opts.t_fft_mult, data.nt as f64)?;
    let nfz = if data.is_3d() { grid::fft_len(1, data.ny as f64)? } else { 1 };
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    // negative-ω band as cpsm.m (sign convention of H⁽²⁾), DC bin excluded
    let keep: Vec<usize> = (0..nfft_t)
        .filter(|&k| omega_full[k] >= -(2.0 * PI * acq.f_high) && omega_full[k] <= -(2.0 * PI * acq.f_low) && omega_full[k] != 0.0)
        .collect();
    if keep.is_empty() {
        return invalid("no FFT bins inside the band");
    }
    let omega: Vec<f64> = keep.iter().map(|&k| omega_full[k]).collect();
    let nw = omega.len();
    let cols = fftutil::time_fft_select(data, nfft_t, &keep);
    let mut w = fftutil::spatial_fft(&cols, nw, data.nx, data.ny, nfn, nfz); // [ik][iw]
    drop(cols);
    // integer angular orders: n = signed index · 2π/(N·Δφ) (= the integer itself when N·Δφ = 2π)
    let order_scale = 2.0 * PI / (nfn as f64 * acq.phi_step);
    let nn: Vec<f64> = (0..nfn).map(|k| grid::signed_index(k, nfn) as f64 * order_scale).collect();
    let kzv = if nfz > 1 { grid::k_vec(nfz, acq.z_step) } else { vec![0.0] };
    let ch = acq.c / 2.0;
    let nk = nfn * nfz;
    let is3d = data.is_3d();

    // k_r per (ik, iw); NaN marks evanescent (removed below)
    let mut kr = vec![0.0f64; nk * nw];
    kr.par_chunks_mut(nw).zip(w.par_chunks_mut(nw)).enumerate().for_each(|(ik, (krc, wc))| {
        let kz = kzv[ik % nfz];
        for iw in 0..nw {
            let a = omega[iw].abs() / ch;
            krc[iw] = if is3d {
                if a > kz.abs() { ((a - kz.abs()) * (a + kz.abs())).sqrt() } else { f64::NAN }
            } else {
                a
            };
            wc[iw] = if krc[iw].is_nan() { C64::new(0.0, 0.0) } else { wc[iw] * C64::new(0.0, -omega[iw] * acq.t_delay).exp() };
        }
    });
    if opts.transfer == TransferFunction::Gardner {
        w.par_chunks_mut(nw).zip(kr.par_chunks(nw)).enumerate().for_each(|(ik, (wc, krc))| {
            let n = nn[ik / nfz];
            for iw in 0..nw {
                if krc[iw].is_nan() { continue; }
                let h = approx_hankel2(krc[iw] * acq.r0, n);
                wc[iw] = if h == C64::new(0.0, 0.0) { C64::new(0.0, 0.0) } else { wc[iw] / h };
            }
        });
    }

    let fft2_inv = Fft2::new(nfn, nfz, true);
    let mut img = Vec::with_capacity(ranges.len() * data.nx * data.ny);
    let mut plane = vec![C64::new(0.0, 0.0); nk];
    for &r in ranges {
        plane.par_iter_mut().enumerate().for_each(|(ik, p)| {
            let n = nn[ik / nfz];
            let (wc, krc) = (&w[ik * nw..(ik + 1) * nw], &kr[ik * nw..(ik + 1) * nw]);
            let mut s = C64::new(0.0, 0.0);
            for iw in 0..nw {
                if krc[iw].is_nan() || wc[iw] == C64::new(0.0, 0.0) { continue; }
                let g = match opts.transfer {
                    TransferFunction::Haun => {
                        let m = krc[iw] * krc[iw] - (n * n) / (acq.r0 * r);
                        if m <= 0.0 { continue; } // evanescent in the Haun approximation: removed
                        (r / acq.r0).sqrt() * C64::new(0.0, -(r - acq.r0) * m.sqrt()).exp()
                    }
                    TransferFunction::Gardner => approx_hankel2(krc[iw] * r, n),
                };
                s += wc[iw] * g;
            }
            *p = s;
        });
        fftutil::inverse_plane_crop(&fft2_inv, &mut plane, data.nx, data.ny, &mut img);
    }
    let phi = (0..data.nx).map(|i| (i as f64 - (data.nx as f64 - 1.0) / 2.0) * acq.phi_step).collect();
    let z = (0..data.ny).map(|i| i as f64 * acq.z_step).collect();
    Ok(CylImage { r: ranges.to_vec(), phi, z, nphi: data.nx, nz: data.ny, data: img })
}
