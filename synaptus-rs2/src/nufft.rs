//! 1D nonuniform FFT, types 1 and 2, with the "exponential of semicircle" (ES)
//! kernel of Barnett, Magland & af Klinteberg, SIAM J. Sci. Comput. 41(5):C479–C504
//! (2019), §3.1.1. Conventions (their eqs. 1.1, 1.3), modes k ∈ I_N = {−N/2 … N/2−1}:
//!
//!   type 1:  f_k = Σ_j c_j e^{+i k x_j}      (nonuniform → uniform)
//!   type 2:  c_j = Σ_k f_k e^{−i k x_j}      (uniform → nonuniform)
//!
//! Outputs/inputs over modes are stored at index k + N/2. Any real x_j is
//! accepted (reduced mod 2π; exact for integer k).
//!
//! Parameters (paper eq. 3.2): σ = 2, w = ⌈log₁₀(1/ε)⌉ + 1, β = 2.30·w.
//! Deconvolution p_k = h/ψ̂(k) (eq. 3.9) with ψ̂ by Gauss–Legendre quadrature
//! (eq. 3.10) using p = 2w + 2 positive nodes (paper: p ≥ 1.5w + 2 suffices).
//!
//! Error: rigorously max_k |f̃_k − f_k| ≤ ε∞‖c‖₁ (paper eq. 4.4) with ε∞ ≈ ε by
//! design; relative ℓ2 error ≈ ε is a heuristic (paper Remark 10).
//! Single-threaded per call; callers parallelise over independent transforms.

use crate::error::{invalid, Result};
use crate::fftutil::C64;
use crate::quad::gauss_legendre;
use rustfft::{Fft, FftPlanner};
use std::f64::consts::PI;
use std::sync::Arc;

pub struct Nufft1d {
    n_modes: usize,
    n_fine: usize,
    w: usize,
    beta: f64,
    alpha: f64,
    h: f64,
    corr: Vec<f64>,
    fwd: Arc<dyn Fft<f64>>,
    inv: Arc<dyn Fft<f64>>,
}

fn next_smooth_even(mut n: usize) -> usize {
    if n % 2 == 1 {
        n += 1;
    }
    loop {
        let mut m = n;
        for p in [2, 3, 5] {
            while m % p == 0 {
                m /= p;
            }
        }
        if m == 1 {
            return n;
        }
        n += 2;
    }
}

impl Nufft1d {
    /// Plan for `n_modes` (must be even, ≥ 2) output/input modes at relative
    /// tolerance `tol` ∈ [1e-15, 1e-1].
    pub fn new(n_modes: usize, tol: f64) -> Result<Self> {
        if n_modes < 2 || n_modes % 2 != 0 {
            return invalid(format!("NUFFT mode count must be even and >= 2 (got {n_modes})"));
        }
        if !(1e-15..=1e-1).contains(&tol) {
            return invalid(format!("NUFFT tolerance must lie in [1e-15, 1e-1] (got {tol})"));
        }
        let w = (1.0 / tol).log10().ceil() as usize + 1;
        let beta = 2.30 * w as f64;
        let n_fine = next_smooth_even((2 * n_modes).max(2 * w));
        let h = 2.0 * PI / n_fine as f64;
        let alpha = PI * w as f64 / n_fine as f64;
        // ψ̂(k) = 2α Σ_{q_j > 0} w_j φ(q_j) cos(α k q_j)   (2p-node rule, positive half)
        let p = 2 * w + 2;
        let (q, qw) = gauss_legendre(2 * p);
        let phi = |z: f64| (beta * ((1.0 - z * z).max(0.0).sqrt() - 1.0)).exp();
        let pos: Vec<(f64, f64)> = q.iter().zip(&qw).filter(|(x, _)| **x > 0.0).map(|(x, w)| (*x, *w * phi(*x))).collect();
        let half = n_modes as i64 / 2;
        let corr = (-half..half)
            .map(|k| {
                let psi_hat: f64 = 2.0 * alpha * pos.iter().map(|(x, wf)| wf * (alpha * k as f64 * x).cos()).sum::<f64>();
                h / psi_hat
            })
            .collect();
        let mut planner = FftPlanner::new();
        Ok(Self {
            n_modes,
            n_fine,
            w,
            beta,
            alpha,
            h,
            corr,
            fwd: planner.plan_fft_forward(n_fine),
            inv: planner.plan_fft_inverse(n_fine),
        })
    }

    pub fn n_modes(&self) -> usize {
        self.n_modes
    }

    #[inline]
    fn phi(&self, z: f64) -> f64 {
        (self.beta * ((1.0 - z * z).max(0.0).sqrt() - 1.0)).exp()
    }

    /// Fine-grid indices l (wrapped) and kernel values ψ(l h − x) for one point.
    #[inline]
    fn footprint(&self, x: f64, mut f: impl FnMut(usize, f64)) {
        let xr = x.rem_euclid(2.0 * PI);
        let lo = ((xr - self.alpha) / self.h).ceil() as i64;
        let hi = ((xr + self.alpha) / self.h).floor() as i64;
        let n = self.n_fine as i64;
        for l in lo..=hi {
            let z = (l as f64 * self.h - xr) / self.alpha;
            if z.abs() <= 1.0 {
                f(l.rem_euclid(n) as usize, self.phi(z));
            }
        }
    }

    /// Type 1: f_k = Σ_j c_j e^{+i k x_j}, returned at index k + N/2.
    pub fn type1(&self, x: &[f64], c: &[C64]) -> Vec<C64> {
        assert_eq!(x.len(), c.len());
        let mut b = vec![C64::new(0.0, 0.0); self.n_fine];
        for (&xj, &cj) in x.iter().zip(c) {
            self.footprint(xj, |l, v| b[l] += cj * v);
        }
        // b̂_k = Σ_l e^{+2πilk/n} b_l  (eq. 3.7): unnormalised inverse FFT.
        self.inv.process(&mut b);
        let half = self.n_modes as i64 / 2;
        let n = self.n_fine as i64;
        (-half..half).zip(&self.corr).map(|(k, p)| b[k.rem_euclid(n) as usize] * *p).collect()
    }

    /// Type 2: c_j = Σ_k f_k e^{−i k x_j}, with f given at index k + N/2.
    pub fn type2(&self, f: &[C64], x: &[f64]) -> Vec<C64> {
        assert_eq!(f.len(), self.n_modes);
        let half = self.n_modes as i64 / 2;
        let n = self.n_fine as i64;
        let mut g = vec![C64::new(0.0, 0.0); self.n_fine];
        for ((k, fk), p) in (-half..half).zip(f).zip(&self.corr) {
            g[k.rem_euclid(n) as usize] = fk * *p;
        }
        // b_l = Σ_k g_k e^{−2πilk/n}: forward FFT.
        self.fwd.process(&mut g);
        x.iter()
            .map(|&xj| {
                let mut acc = C64::new(0.0, 0.0);
                self.footprint(xj, |l, v| acc += g[l] * v);
                acc
            })
            .collect()
    }

    /// Kernel width w (grid points) — for diagnostics.
    pub fn width(&self) -> usize {
        self.w
    }
}
