//! Literal port of `matlab/core/mulok.m`, `interpol = 'linear'` path.
//! The `'chirpz'` path is not ported.

use crate::error::{invalid, Error, Result};
use crate::fftutil::{self, Fft2, C64};
use crate::grid;
use crate::psm::{layer_depths, spatial_fft_sizes};
use crate::types::{validate, Acquisition, FocusedImage, LayerImage, Medium, PulseEchoData};
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct LegacyMulokOptions {
    pub t_fft_mult: usize,
    pub x_fft_mult: usize,
    pub y_fft_mult: usize,
    pub z_fft_mult: usize,
    /// Band upsampling factor (`upSamp`, default 4).
    pub up_samp: usize,
}

impl Default for LegacyMulokOptions {
    fn default() -> Self {
        Self { t_fft_mult: 1, x_fft_mult: 1, y_fft_mult: 1, z_fft_mult: 1, up_samp: 4 }
    }
}

/// Octave `interp1(x, y, xi, 'linear')` for ascending `x`; `None` outside [x₀, x_{n−1}].
/// Octave evaluates the piecewise polynomial as slope·(xi − x_k) + y_k.
fn interp1_linear(x: &[f64], y: &[C64], xi: f64) -> Option<C64> {
    let n = x.len();
    if n < 2 || xi.is_nan() || xi < x[0] || xi > x[n - 1] {
        return None;
    }
    // largest k with x[k] <= xi, clamped to the last interval
    let k = (x.partition_point(|&v| v <= xi).max(1) - 1).min(n - 2);
    let dx = x[k + 1] - x[k];
    let dy = y[k + 1] - y[k];
    let slope = C64::new(dy.re / dx, dy.im / dx);
    Some(slope * (xi - x[k]) + y[k])
}

pub fn mulok(data: &PulseEchoData, acq: &Acquisition, medium: &Medium, opts: &LegacyMulokOptions) -> Result<FocusedImage> {
    let z_offset = validate(data, acq, medium)?;
    if opts.up_samp == 0 {
        return invalid("up_samp must be >= 1");
    }
    let (fs, f_low, f_high) = (acq.fs, acq.f_low, acq.f_high);
    let nl = medium.n_layers();
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let nfft_t = grid::fft_len(opts.t_fft_mult, data.nt as f64)?; // hh = 1
    let d_omega = (2.0 * PI * fs) / nfft_t as f64;
    let omega = grid::omega_vec(nfft_t, fs);
    let keep: Vec<usize> =
        (0..nfft_t).filter(|&k| omega[k] >= -(2.0 * PI * f_high) && omega[k] <= -(2.0 * PI * f_low)).collect();
    let n_band = keep.len();
    if n_band < 2 {
        return invalid("fewer than 2 FFT bins inside the band");
    }
    let mut omega_band: Vec<f64> = keep.iter().map(|&k| omega[k]).collect();
    let mut cols = fftutil::time_fft_select(data, nfft_t, &keep); // [c][iw]
    let mut nw = n_band;

    // Band upsampling: Poxy = fft(ifft(Poxy), up*n), with the original omega grid.
    if opts.up_samp != 1 {
        let m = opts.up_samp * n_band;
        let mut p = FftPlanner::new();
        let (inv, fwd) = (p.plan_fft_inverse(n_band), p.plan_fft_forward(m));
        let mut out = vec![C64::new(0.0, 0.0); cols.len() / n_band * m];
        out.par_chunks_mut(m).zip(cols.par_chunks(n_band)).for_each(|(o, c)| {
            let mut b = c.to_vec();
            inv.process(&mut b);
            let s = 1.0 / n_band as f64;
            for (oi, bi) in o.iter_mut().zip(&b) {
                *oi = bi * s;
            }
            fwd.process(o);
        });
        cols = out;
        let step = d_omega / opts.up_samp as f64;
        let end = *omega_band.last().unwrap();
        omega_band = (0..m).map(|j| -((m - 1 - j) as f64 * step) + end).collect();
        nw = m;
    }
    // tDelay compensation on the (possibly mislabelled) grid, as in mulok.m.
    let delay: Vec<C64> = omega_band.iter().map(|&w| C64::new(0.0, -w * acq.t_delay).exp()).collect();
    cols.par_chunks_mut(nw).for_each(|c| c.iter_mut().zip(&delay).for_each(|(v, d)| *v *= d));
    let mut w = fftutil::spatial_fft(&cols, nw, data.nx, data.ny, nfx, nfy); // [ik][iw]
    drop(cols);
    let kxv = grid::k_vec(nfx, acq.x_step);
    let kyv = if nfy > 1 { grid::k_vec(nfy, acq.y_step) } else { vec![0.0] };
    let is3d = nfy > 1 || data.is_3d();
    let nk = nfx * nfy;

    let mut layers = Vec::with_capacity(nl);
    for l in 0..nl {
        let c = medium.velocities[l];
        let dz = (c / 2.0) / (f_high - f_low);
        let n_fft_z = grid::fft_len(opts.z_fft_mult, (data.nt as f64 * (f_high - f_low)) / (fs / 2.0))?;
        let kzs = (2.0 * PI) / dz;
        let kz_pos: Vec<f64> = (0..n_fft_z)
            .map(|j| (2.0 * PI * f_low) / (c / 2.0) + j as f64 * (kzs / n_fft_z as f64))
            .collect();
        let kz: Vec<f64> = (0..n_fft_z).map(|j| -kz_pos[n_fft_z - 1 - j]).collect();
        let depths = layer_depths(medium, l, dz, z_offset);
        if depths.len() > n_fft_z {
            return Err(Error::SizeOverflow(format!(
                "layer {l}: {} image planes exceed the z FFT length {n_fft_z} (MATLAB would fail with an index error)",
                depths.len()
            )));
        }

        // Stolt resampling per (kx,ky) column; result layout [ik][ikz].
        let mut pkz = vec![C64::new(0.0, 0.0); nk * n_fft_z];
        pkz.par_chunks_mut(n_fft_z).enumerate().for_each(|(ik, out)| {
            let (kx, ky) = (kxv[ik / nfy], kyv[ik % nfy]);
            let col = &w[ik * nw..(ik + 1) * nw];
            for (j, o) in out.iter_mut().enumerate() {
                let kzj = kz[j];
                let (kk, a) = if is3d {
                    (kzj * kzj + kx * kx + ky * ky, 1.0 / (1.0 + (kx * kx + ky * ky) / (kzj * kzj)))
                } else {
                    (kzj * kzj + kx * kx, 1.0 / (1.0 + (kx * kx) / (kzj * kzj)))
                };
                let a = if a.is_nan() { 0.0 } else { a };
                let kk_st = -(kk * if kk > 0.0 { 1.0 } else { 0.0 }).sqrt();
                let om = (c / 2.0) * kk_st;
                let mut v = match interp1_linear(&omega_band, col, om) {
                    Some(v) if !(v.re.is_nan() || v.im.is_nan()) => v,
                    _ => C64::new(0.0, 0.0),
                };
                v *= a;
                if l == 0 {
                    v *= C64::new(0.0, kzj * z_offset).exp();
                }
                *o = v;
            }
        });
        // ifftn over (kz, kx[, ky]) with 1/N normalisation.
        let inv_z = FftPlanner::new().plan_fft_inverse(n_fft_z);
        pkz.par_chunks_mut(n_fft_z).for_each(|col| inv_z.process(col));
        let fft2_inv = Fft2::new(nfx, nfy, true);
        let zn = 1.0 / n_fft_z as f64;
        let mut img = Vec::with_capacity(depths.len() * data.nx * data.ny);
        let mut plane = vec![C64::new(0.0, 0.0); nk];
        for iz in 0..depths.len() {
            for ik in 0..nk {
                plane[ik] = pkz[ik * n_fft_z + iz] * zn;
            }
            fftutil::inverse_plane_crop(&fft2_inv, &mut plane, data.nx, data.ny, &mut img);
        }
        layers.push(LayerImage { nz: depths.len(), nx: data.nx, ny: data.ny, z: depths, data: img });

        // Migrate to next interface. Note: mulok.m does NOT zero evanescent
        // components here (they get phase factor 1); reproduced as is.
        if l + 1 < nl {
            let d = medium.thicknesses[l];
            w.par_chunks_mut(nw).enumerate().for_each(|(ik, col)| {
                let (kx, ky) = (kxv[ik / nfy], kyv[ik % nfy]);
                for (iw, v) in col.iter_mut().enumerate() {
                    let kk = omega_band[iw] / (c / 2.0);
                    let kz2 = if is3d { kk * kk - kx * kx - ky * ky } else { kk * kk - kx * kx };
                    let kzp = (kz2 * if kz2 > 0.0 { 1.0 } else { 0.0 }).sqrt();
                    *v *= C64::new(0.0, -kzp * d).exp();
                }
            });
        }
    }
    let x = (0..data.nx).map(|i| i as f64 * acq.x_step).collect();
    let y = (0..data.ny).map(|i| i as f64 * acq.y_step).collect();
    Ok(FocusedImage { layers, x, y })
}
