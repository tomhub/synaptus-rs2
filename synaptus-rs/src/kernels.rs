//! Numerical kernels: vertical wavenumber, phase factors, summation.

use num_complex::Complex64 as C64;

/// How kz² = (2ω/c)² − k⊥² is evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KzFormula {
    /// MATLAB form `(2/c)^2*ω.^2 - kx.^2 - ky.^2`. Suffers cancellation near the
    /// evanescent cut-off (relative error of kz grows like u·(2ω/c)²/kz²).
    Legacy,
    /// Factored form (a − k⊥)(a + k⊥) with a = 2|ω|/c and k⊥ = hypot(kx, ky).
    /// Each factor is exact or nearly so (Sterbenz) near the cut-off, giving
    /// O(u) relative error in kz² everywhere.
    Factored,
}

/// kz for two-way propagation (ĉ = c/2). Returns `None` for evanescent components
/// (kz² < 0). `two_over_c_sq` must be `(2/c)^2` computed by the caller (MATLAB order).
#[inline]
pub fn kz_two_way(formula: KzFormula, omega: f64, kx: f64, ky: f64, c: f64) -> Option<f64> {
    match formula {
        KzFormula::Legacy => {
            let a2 = (2.0 / c) * (2.0 / c);
            let kz2 = a2 * (omega * omega) - kx * kx - ky * ky;
            if kz2 >= 0.0 { Some(kz2.sqrt()) } else { None }
        }
        KzFormula::Factored => {
            let a = 2.0 * omega.abs() / c;
            let kp = kx.hypot(ky);
            if a >= kp { Some(((a - kp) * (a + kp)).sqrt()) } else { None }
        }
    }
}

/// e^{iφ}.
#[inline]
pub fn cis(phi: f64) -> C64 {
    let (s, c) = phi.sin_cos();
    C64::new(c, s)
}

/// Summation strategy for the ω integral (imaging condition t = 0).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Summation {
    /// Plain left-to-right accumulation (what MATLAB/Octave `sum` does for a column).
    Naive,
    /// Neumaier compensated summation applied to real and imaginary parts.
    Neumaier,
}

#[inline]
pub fn sum_naive(v: impl Iterator<Item = C64>) -> C64 {
    let mut s = C64::new(0.0, 0.0);
    for x in v {
        s += x;
    }
    s
}

/// Neumaier (improved Kahan–Babuška) summation, componentwise. Error bound
/// |ŝ − s| ≤ 2u|s| + O(n u²)·Σ|x_i| instead of (n−1)u·Σ|x_i| for naive summation
/// (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §4.3).
#[inline]
pub fn sum_neumaier(v: impl Iterator<Item = C64>) -> C64 {
    let (mut sr, mut cr, mut si, mut ci) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for x in v {
        let t = sr + x.re;
        cr += if sr.abs() >= x.re.abs() { (sr - t) + x.re } else { (x.re - t) + sr };
        sr = t;
        let t = si + x.im;
        ci += if si.abs() >= x.im.abs() { (si - t) + x.im } else { (x.im - t) + si };
        si = t;
    }
    C64::new(sr + cr, si + ci)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn neumaier_beats_naive_on_cancellation() {
        // 1 + 1e100 + 1 - 1e100 = 2 exactly; naive gives 0.
        let v = [1.0, 1e100, 1.0, -1e100].map(|r| C64::new(r, -r));
        assert_eq!(sum_naive(v.into_iter()).re, 0.0);
        let s = sum_neumaier(v.into_iter());
        assert_eq!(s, C64::new(2.0, -2.0));
    }
    #[test]
    fn kz_formulas_agree_away_from_cutoff_and_mask_evanescent() {
        let (w, c) = (2.0 * std::f64::consts::PI * 2e6, 1480.0);
        for &kx in &[0.0, 1000.0, 5000.0] {
            let a = kz_two_way(KzFormula::Legacy, w, kx, 0.0, c).unwrap();
            let b = kz_two_way(KzFormula::Factored, w, kx, 0.0, c).unwrap();
            assert!((a - b).abs() <= 1e-12 * a);
        }
        let kmax = 2.0 * w / c;
        assert!(kz_two_way(KzFormula::Factored, w, 1.01 * kmax, 0.0, c).is_none());
        assert!(kz_two_way(KzFormula::Legacy, w, 1.01 * kmax, 0.0, c).is_none());
    }
}
