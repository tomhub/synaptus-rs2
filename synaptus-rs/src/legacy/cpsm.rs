//! Literal port of `matlab/core/cpsm.m` ('haun' and 'gardner'), kept to reproduce
//! original results (test T7). It reproduces the original's behaviour on purpose:
//! negative-ω band including the DC bin when f_low = 0; evanescent components
//! clamped to zero phase instead of removed; Haun amplitude √(r/r0) (dossier C8);
//! angular FFT zero-padded to a power of two, giving non-integer orders (C10).
//! Use [`crate::cpsm::cpsm`] unless you need exactly what `cpsm.m` computes.

use crate::cpsm::{approx_hankel2, CylAcquisition, CylImage, TransferFunction};
use crate::error::{invalid, Result};
use crate::fftutil::{self, Fft2, C64};
use crate::grid;
use crate::types::PulseEchoData;
use rayon::prelude::*;
use std::f64::consts::PI;

pub fn cpsm(data: &PulseEchoData, acq: &CylAcquisition, ranges: &[f64], tf: TransferFunction, t_fft_mult: usize) -> Result<CylImage> {
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
    let nfft_t = grid::fft_len(t_fft_mult, data.nt as f64)?;
    let nfn = grid::fft_len(1, data.nx as f64)?;
    let nfz = if data.is_3d() { grid::fft_len(1, data.ny as f64)? } else { 1 };
    let omega_full = grid::omega_vec(nfft_t, acq.fs);
    let keep: Vec<usize> = (0..nfft_t)
        .filter(|&k| omega_full[k] >= -(2.0 * PI * acq.f_high) && omega_full[k] <= -(2.0 * PI * acq.f_low))
        .collect();
    if keep.is_empty() {
        return invalid("no FFT bins inside the band");
    }
    let omega: Vec<f64> = keep.iter().map(|&k| omega_full[k]).collect();
    let nw = omega.len();
    let cols = fftutil::time_fft_select(data, nfft_t, &keep);
    let mut w = fftutil::spatial_fft(&cols, nw, data.nx, data.ny, nfn, nfz); // [ik][iw]
    drop(cols);
    let nn = grid::fft_freqs(nfn, (2.0 * PI) / acq.phi_step);
    let kzv = if nfz > 1 { grid::k_vec(nfz, acq.z_step) } else { vec![0.0] };
    let ch = acq.c / 2.0;
    let nk = nfn * nfz;

    // KR per (ik, iw), and t_delay compensation (applied after KR as in cpsm.m).
    let mut kr = vec![0.0f64; nk * nw];
    let is3d = data.is_3d();
    kr.par_chunks_mut(nw).zip(w.par_chunks_mut(nw)).enumerate().for_each(|(ik, (krc, wc))| {
        let kz = kzv[ik % nfz];
        for iw in 0..nw {
            krc[iw] = if is3d {
                let kr2 = (omega[iw] / ch) * (omega[iw] / ch) - kz * kz;
                (kr2 * if kr2 >= 0.0 { 1.0 } else { 0.0 }).sqrt()
            } else {
                omega[iw].abs() / ch
            };
            wc[iw] *= C64::new(0.0, -omega[iw] * acq.t_delay).exp();
        }
    });
    if tf == TransferFunction::Gardner {
        w.par_chunks_mut(nw).zip(kr.par_chunks(nw)).enumerate().for_each(|(ik, (wc, krc))| {
            let n = nn[ik / nfz];
            for iw in 0..nw {
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
                let g = match tf {
                    TransferFunction::Haun => {
                        let m = krc[iw] * krc[iw] - (n * n) / (acq.r0 * r);
                        let m = (m * if m > 0.0 { 1.0 } else { 0.0 }).sqrt();
                        (r / acq.r0).sqrt() * C64::new(0.0, -(r - acq.r0) * m).exp()
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
