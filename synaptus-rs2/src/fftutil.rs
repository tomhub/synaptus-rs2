//! FFT helpers on top of `rustfft` (unnormalised in both directions; callers
//! apply 1/N for inverse transforms explicitly, matching MATLAB `ifft`).

use num_complex::Complex64;
use rayon::prelude::*;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

pub type C64 = Complex64;

/// Plans for a 2D transform over a row-major plane `[ix][iy]` of size nfx × nfy.
pub struct Fft2 {
    pub nfx: usize,
    pub nfy: usize,
    fx: Arc<dyn Fft<f64>>,
    fy: Arc<dyn Fft<f64>>,
}

impl Fft2 {
    pub fn new(nfx: usize, nfy: usize, inverse: bool) -> Self {
        let mut p = FftPlanner::new();
        let (fx, fy) = if inverse {
            (p.plan_fft_inverse(nfx), p.plan_fft_inverse(nfy))
        } else {
            (p.plan_fft_forward(nfx), p.plan_fft_forward(nfy))
        };
        Self { nfx, nfy, fx, fy }
    }

    /// In-place unnormalised 2D transform (y then x). `col` is scratch of length nfx.
    pub fn process(&self, plane: &mut [C64], col: &mut Vec<C64>) {
        debug_assert_eq!(plane.len(), self.nfx * self.nfy);
        if self.nfy > 1 {
            for row in plane.chunks_exact_mut(self.nfy) {
                self.fy.process(row);
            }
        }
        if self.nfx > 1 {
            col.resize(self.nfx, C64::new(0.0, 0.0));
            for iy in 0..self.nfy {
                for ix in 0..self.nfx {
                    col[ix] = plane[ix * self.nfy + iy];
                }
                self.fx.process(col);
                for ix in 0..self.nfx {
                    plane[ix * self.nfy + iy] = col[ix];
                }
            }
        }
    }
}

/// Transform columns `cols[c][iw]` (c = ix + nx*iy, nw frequencies each) to the
/// spatial-frequency domain with zero padding to nfx × nfy. Output layout is
/// `[ik][iw]` with ik = ikx*nfy + iky (ω contiguous per wavenumber column).
pub fn spatial_fft(cols: &[C64], nw: usize, nx: usize, ny: usize, nfx: usize, nfy: usize) -> Vec<C64> {
    let nk = nfx * nfy;
    let fft2 = Fft2::new(nfx, nfy, false);
    // planes[iw][ik], computed in parallel, then transposed.
    let planes: Vec<Vec<C64>> = (0..nw)
        .into_par_iter()
        .map(|iw| {
            let mut plane = vec![C64::new(0.0, 0.0); nk];
            for iy in 0..ny {
                for ix in 0..nx {
                    plane[ix * nfy + iy] = cols[(ix + nx * iy) * nw + iw];
                }
            }
            let mut scratch = Vec::new();
            fft2.process(&mut plane, &mut scratch);
            plane
        })
        .collect();
    let mut out = vec![C64::new(0.0, 0.0); nk * nw];
    out.par_chunks_mut(nw).enumerate().for_each(|(ik, colw)| {
        for (iw, v) in colw.iter_mut().enumerate() {
            *v = planes[iw][ik];
        }
    });
    out
}

/// Inverse 2D transform of a wavenumber plane (1/(nfx·nfy) normalisation, as
/// MATLAB `ifft` along x then y), cropped to nx × ny, appended to `dst` as `[ix][iy]`.
pub fn inverse_plane_crop(fft2_inv: &Fft2, plane: &mut [C64], nx: usize, ny: usize, dst: &mut Vec<C64>) {
    let mut scratch = Vec::new();
    fft2_inv.process(plane, &mut scratch);
    let norm = 1.0 / (fft2_inv.nfx * fft2_inv.nfy) as f64;
    for ix in 0..nx {
        for iy in 0..ny {
            dst.push(plane[ix * fft2_inv.nfy + iy] * norm);
        }
    }
}

/// Zero-padded forward FFT of every real column of the data along t, keeping only
/// the bins selected by `keep` (in the order given). Output `[c][j]`.
pub fn time_fft_select(
    data: &crate::PulseEchoData,
    nfft: usize,
    keep: &[usize],
) -> Vec<C64> {
    let plan = FftPlanner::new().plan_fft_forward(nfft);
    let ncol = data.nx * data.ny;
    let nk = keep.len();
    let mut out = vec![C64::new(0.0, 0.0); ncol * nk];
    out.par_chunks_mut(nk).enumerate().for_each(|(c, dst)| {
        let (ix, iy) = (c % data.nx, c / data.nx);
        let mut buf = vec![C64::new(0.0, 0.0); nfft];
        for (b, &v) in buf.iter_mut().zip(data.column(ix, iy)) {
            *b = C64::new(v, 0.0);
        }
        plan.process(&mut buf);
        for (d, &k) in dst.iter_mut().zip(keep) {
            *d = buf[k];
        }
    });
    out
}

/// Full forward FFT of a real vector zero-padded to n.
pub fn fft_real_padded(x: &[f64], n: usize) -> Vec<C64> {
    let mut buf = vec![C64::new(0.0, 0.0); n];
    for (b, &v) in buf.iter_mut().zip(x) {
        *b = C64::new(v, 0.0);
    }
    FftPlanner::new().plan_fft_forward(n).process(&mut buf);
    buf
}
