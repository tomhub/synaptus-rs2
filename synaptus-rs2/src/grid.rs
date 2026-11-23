//! Frequency / wavenumber grids with MATLAB-identical arithmetic.

use crate::error::{invalid, Error, Result};
use std::f64::consts::PI;

/// MATLAB `nextpow2` for x > 0, restricted to non-negative exponents:
/// smallest p >= 0 such that 2^p >= x.
pub fn nextpow2(x: f64) -> Result<u32> {
    if !(x.is_finite() && x > 0.0) {
        return invalid(format!("nextpow2 argument must be finite and > 0 (got {x})"));
    }
    let mut p = 0u32;
    while (2f64).powi(p as i32) < x {
        p += 1;
        if p > 62 {
            return Err(Error::SizeOverflow(format!("nextpow2({x}) too large")));
        }
    }
    Ok(p)
}

/// FFT length `mult * 2^nextpow2(x)` with overflow checks.
pub fn fft_len(mult: usize, x: f64) -> Result<usize> {
    if mult == 0 {
        return invalid("FFT size multiplier must be >= 1");
    }
    let p = nextpow2(x)?;
    1usize
        .checked_shl(p)
        .and_then(|v| v.checked_mul(mult))
        .ok_or_else(|| Error::SizeOverflow(format!("FFT length {mult}*2^{p}")))
}

/// Signed frequency index of FFT bin `k` of an `n`-point FFT (as `ifftshift` of
/// `(0:n-1) - floor(n/2)`, equivalently NumPy `fftfreq * n`).
pub fn signed_index(k: usize, n: usize) -> i64 {
    if k < n.div_ceil(2) {
        k as i64
    } else {
        k as i64 - n as i64
    }
}

/// FFT-ordered frequency vector: `m * (scale / n)` with m the signed index.
/// The operation order matches MATLAB's `((0:n-1)-floor(n/2)) * (scale/n)`
/// so that band masks are bitwise identical.
pub fn fft_freqs(n: usize, scale: f64) -> Vec<f64> {
    let step = scale / n as f64;
    (0..n).map(|k| signed_index(k, n) as f64 * step).collect()
}

/// Angular frequency vector ω for an `n`-point time FFT at sampling rate `fs`.
pub fn omega_vec(n: usize, fs: f64) -> Vec<f64> {
    fft_freqs(n, 2.0 * PI * fs)
}

/// Wavenumber vector k for an `n`-point spatial FFT with step `dx`.
pub fn k_vec(n: usize, dx: f64) -> Vec<f64> {
    fft_freqs(n, (2.0 * PI) / dx)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nextpow2_matches_matlab_semantics() {
        assert_eq!(nextpow2(1.0).unwrap(), 0);
        assert_eq!(nextpow2(2.0).unwrap(), 1);
        assert_eq!(nextpow2(3.0).unwrap(), 2);
        assert_eq!(nextpow2(1024.0).unwrap(), 10);
        assert_eq!(nextpow2(1024.5).unwrap(), 11);
        assert!(nextpow2(0.0).is_err());
        assert!(nextpow2(f64::NAN).is_err());
    }
    #[test]
    fn signed_index_even_odd() {
        let e: Vec<i64> = (0..4).map(|k| signed_index(k, 4)).collect();
        assert_eq!(e, vec![0, 1, -2, -1]);
        let o: Vec<i64> = (0..5).map(|k| signed_index(k, 5)).collect();
        assert_eq!(o, vec![0, 1, 2, -2, -1]);
        assert_eq!(signed_index(0, 1), 0);
    }
}
