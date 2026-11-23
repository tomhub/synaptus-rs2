//! **B. SQ-PSM** (spectral-quadrature phase shift migration).
//!
//! Target (dossier v2 §5.2): the m → ∞ limit of PSM with a time FFT zero-padded
//! m-fold (sum scaled by 1/m), i.e. the continuous integral
//!
//!   J_l(z; k⊥) = (1/Δω₁) ∫_{max(ω_lo, ω*_l)}^{ω_hi} W(ω,k⊥) e^{−iω t_d} Φ_l(ω) e^{i k_{z,l}(ω) z} dω,
//!
//! where W is the DTFT of the recorded samples, Φ_l the phase to the top of layer
//! l, ω*_l = max_{j≤l} c_j k⊥/2 the cumulative evanescent cut-off and
//! Δω₁ = 2π f_s / N_fft(m=1) (so J ≈ v1 PSM in amplitude).
//!
//! Method: ω = ω*_l + s² removes the √ branch point at the cut-off (§5.3), then
//! n-point Gauss–Legendre in s (geometric convergence for analytic integrands,
//! Trefethen 2008 Thm 4.5); W at the nodes by a type-2 NUFFT from the time
//! samples (§5.5, exact DTFT up to ε); the depth column by a type-1 NUFFT (§5.1).
//! Node count n = ⌈κ·N_b⌉ + 16 rounded up to a multiple of 32, with
//! N_b = (interval length)·T_eff/2π and T_eff = max(2π/Δω₁, t_d + n_t/f_s)
//! (§5.4 as corrected by TEST_PLAN_V2 amendment A2; κ fixed by test N5-A2).
//!
//! Not supported: matched filtering (`hh`).

use crate::error::{invalid, Result};
use crate::fftutil::{self, C64};
use crate::grid;
use crate::kernels::cis;
use crate::nufft::Nufft1d;
use crate::psm::{layer_depths, spatial_fft_sizes, DepthStep};
use crate::psm_nufft::{columns_to_layer, even_at_least};
use crate::quad::GaussCache;
use crate::types::{validate, Acquisition, FocusedImage, Medium, PulseEchoData};
use rayon::prelude::*;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct PsmSqOptions {
    pub depth_step: DepthStep,
    pub x_fft_mult: usize,
    pub y_fft_mult: usize,
    /// Gauss–Legendre nodes per Δω₁ of integration interval.
    pub kappa: f64,
    /// NUFFT relative tolerance (both transforms).
    pub tol: f64,
    /// Matched filter: unsupported, must be `None`.
    pub hh: Option<Vec<f64>>,
}

/// Default node factor. Fixed by test N5 (TEST_PLAN_V2) — see ACCURACY_V2.md.
pub const DEFAULT_KAPPA: f64 = 2.0;

impl Default for PsmSqOptions {
    fn default() -> Self {
        Self { depth_step: DepthStep::HighFrequency, x_fft_mult: 1, y_fft_mult: 1, kappa: DEFAULT_KAPPA, tol: 1e-12, hh: None }
    }
}

pub fn psm_sq(data: &PulseEchoData, acq: &Acquisition, medium: &Medium, opts: &PsmSqOptions) -> Result<FocusedImage> {
    let z_offset = validate(data, acq, medium)?;
    if opts.hh.is_some() {
        return invalid("psm_sq does not support matched filtering (hh)");
    }
    if !(opts.kappa.is_finite() && opts.kappa > 0.0) {
        return invalid(format!("kappa must be finite and > 0 (got {})", opts.kappa));
    }
    let nfft1 = grid::fft_len(1, data.nt as f64)?;
    let dw1 = 2.0 * PI * acq.fs / nfft1 as f64;
    let t_eff = (2.0 * PI / dw1).max(acq.t_delay + data.nt as f64 / acq.fs);
    let (nfx, nfy) = spatial_fft_sizes(data, opts.x_fft_mult, opts.y_fft_mult)?;
    let nk = nfx * nfy;
    let (wl, wh) = (2.0 * PI * acq.f_low, 2.0 * PI * acq.f_high);
    let kx = grid::k_vec(nfx, acq.x_step);
    let ky = if nfy > 1 { grid::k_vec(nfy, acq.y_step) } else { vec![0.0] };

    // Spatial FFT of the time samples: layout [ik][n], n = 0..nt-1.
    let nt = data.nt;
    let mut cols_t = vec![C64::new(0.0, 0.0); nt * data.nx * data.ny];
    for iy in 0..data.ny {
        for ix in 0..data.nx {
            let c = ix + data.nx * iy;
            for (n, &v) in data.column(ix, iy).iter().enumerate() {
                cols_t[c * nt + n] = C64::new(v, 0.0);
            }
        }
    }
    let p_k = fftutil::spatial_fft(&cols_t, nt, data.nx, data.ny, nfx, nfy);
    drop(cols_t);

    // Type-2 plan over time samples: W(ω) = Σ_n P_n e^{−iωn/fs} = e^{−i(N/2)x} Σ_k P_{k+N/2} e^{−ikx}, x = ω/fs.
    let n2 = even_at_least(nt.max(2));
    let plan_w = Nufft1d::new(n2, opts.tol)?;
    let half_t = (n2 / 2) as f64;

    let nl = medium.n_layers();
    let dzs: Vec<f64> = medium
        .velocities
        .iter()
        .map(|&c| match opts.depth_step {
            DepthStep::HighFrequency => (c / 2.0) / acq.f_high,
            DepthStep::Bandwidth => (c / 2.0) / (acq.f_high - acq.f_low),
        })
        .collect();
    let depths: Vec<Vec<f64>> = (0..nl).map(|l| layer_depths(medium, l, dzs[l], z_offset)).collect();
    let plans_z: Vec<Nufft1d> = depths.iter().map(|d| Nufft1d::new(even_at_least(d.len().max(2)), opts.tol)).collect::<Result<_>>()?;
    let gl = GaussCache::default();

    // per column, per layer: depth column of length nz_l
    let col_results: Vec<Vec<Vec<C64>>> = (0..nk)
        .into_par_iter()
        .map(|ik| {
            let kp = kx[ik / nfy].hypot(ky[ik % nfy]);
            let mut f = p_k[ik * nt..(ik + 1) * nt].to_vec();
            f.resize(n2, C64::new(0.0, 0.0));
            // rotate so that index k+N/2 holds P_{k+N/2}: identity (f[i] = P_i), as type2 expects f at k+N/2
            let mut out = Vec::with_capacity(nl);
            for l in 0..nl {
                let nz = depths[l].len();
                // cumulative cut-off and per-layer offsets d_j = ω* − ω_{c,j} ≥ 0
                let wc: Vec<f64> = (0..=l).map(|j| medium.velocities[j] * kp / 2.0).collect();
                let wstar = wc.iter().cloned().fold(0.0, f64::max);
                let a = wl.max(wstar);
                if a >= wh {
                    out.push(vec![C64::new(0.0, 0.0); nz]);
                    continue;
                }
                let (s_lo, s_hi) = ((a - wstar).sqrt(), (wh - wstar).sqrt());
                // Amendment A2: the phase rate in ω reaches t_d + record length in shallow layers.
                let nb = (wh - a) * t_eff / (2.0 * PI);
                let n = (((opts.kappa * nb).ceil() as usize + 16 + 31) / 32) * 32;
                let rule = gl.get(n);
                let (t, tw) = (&rule.0, &rule.1);
                let half_ds = 0.5 * (s_hi - s_lo);
                let mut omega = Vec::with_capacity(n);
                let mut weight = Vec::with_capacity(n);
                let mut ssq = Vec::with_capacity(n);
                for (ti, wi) in t.iter().zip(tw) {
                    let s = s_lo + half_ds * (ti + 1.0);
                    omega.push(wstar + s * s);
                    ssq.push(s * s);
                    weight.push(wi * half_ds * 2.0 * s / dw1);
                }
                let xw: Vec<f64> = omega.iter().map(|w| w / acq.fs).collect();
                let wv = plan_w.type2(&f, &xw);
                // kz_j(ω) = (2/c_j) √((ω−ω_cj)(ω+ω_cj)), with ω−ω_cj = s² + (ω*−ω_cj) (no cancellation)
                let kz_of = |j: usize, q: usize| -> f64 {
                    let dm = ssq[q] + (wstar - wc[j]);
                    let dp = omega[q] + wc[j];
                    (2.0 / medium.velocities[j]) * (dm * dp).sqrt()
                };
                let dz = dzs[l];
                let start = if l == 0 { z_offset } else { 0.0 };
                let plan = &plans_z[l];
                let half_z = (plan.n_modes() / 2) as f64;
                let mut xs = Vec::with_capacity(n);
                let mut cs = Vec::with_capacity(n);
                for q in 0..n {
                    let kzl = kz_of(l, q);
                    let mut phase = -omega[q] * acq.t_delay - half_t * xw[q] + kzl * start + half_z * kzl * dz;
                    for j in 0..l {
                        phase += kz_of(j, q) * medium.thicknesses[j];
                    }
                    xs.push(kzl * dz);
                    cs.push(wv[q] * cis(phase) * weight[q]);
                }
                let col = plan.type1(&xs, &cs);
                out.push(col[..nz].to_vec());
            }
            out
        })
        .collect();

    let mut layers = Vec::with_capacity(nl);
    for l in 0..nl {
        let nz = depths[l].len();
        let mut cols = vec![C64::new(0.0, 0.0); nk * nz];
        for ik in 0..nk {
            cols[ik * nz..(ik + 1) * nz].copy_from_slice(&col_results[ik][l]);
        }
        layers.push(columns_to_layer(&cols, nz, nfx, nfy, data, depths[l].clone()));
    }
    let x = (0..data.nx).map(|i| i as f64 * acq.x_step).collect();
    let y = (0..data.ny).map(|i| i as f64 * acq.y_step).collect();
    Ok(FocusedImage { layers, x, y })
}
