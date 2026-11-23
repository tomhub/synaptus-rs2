//! Envelope via the analytic signal (as MATLAB/Octave `hilbert` and
//! `scipy.signal.hilbert`: FFT of length N, no padding).

use crate::fftutil::C64;
use crate::types::PulseEchoData;
use rayon::prelude::*;
use rustfft::FftPlanner;

/// |analytic signal| of every t-column; same layout as the input samples.
pub fn envelope(data: &PulseEchoData) -> Vec<f64> {
    let n = data.nt;
    let mut p = FftPlanner::new();
    let (f, fi) = (p.plan_fft_forward(n), p.plan_fft_inverse(n));
    // h: 1 at DC, 2 for positive frequencies, 1 at Nyquist (even n), 0 for negative.
    let h: Vec<f64> = (0..n)
        .map(|k| {
            if k == 0 || (n % 2 == 0 && k == n / 2) {
                1.0
            } else if k < n.div_ceil(2) {
                2.0
            } else {
                0.0
            }
        })
        .collect();
    let mut out = vec![0.0; data.samples.len()];
    out.par_chunks_mut(n).zip(data.samples.par_chunks(n)).for_each(|(o, x)| {
        let mut b: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
        f.process(&mut b);
        for (bi, hi) in b.iter_mut().zip(&h) {
            *bi *= *hi;
        }
        fi.process(&mut b);
        let s = 1.0 / n as f64;
        for (oi, bi) in o.iter_mut().zip(&b) {
            *oi = bi.norm() * s;
        }
    });
    out
}
