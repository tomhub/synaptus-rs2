//! Synthetic pulse-echo data with analytically known scatterer positions.
//! Independent of the library: no imports from `synaptus` except plain data types.
#![allow(dead_code)]

use std::f64::consts::PI;

pub struct Pulse {
    pub fc: f64,
    pub sigma_t: f64,
}

impl Pulse {
    /// Gaussian-modulated cosine g(t) = exp(−t²/2σ²)·cos(2πf_c t).
    pub fn eval(&self, t: f64) -> f64 {
        (-(t * t) / (2.0 * self.sigma_t * self.sigma_t)).exp() * (2.0 * PI * self.fc * t).cos()
    }
}

/// Frequency-independent directivity weight exp(−(sinθ/s0)²).
pub fn directivity(sin_theta: f64, s0: f64) -> f64 {
    (-(sin_theta / s0).powi(2)).exp()
}

/// One-way travel time and launch-angle sine from (x, 0) through a horizontal
/// interface at depth d1 (velocity c1 above, c2 below) to (x0, z0), z0 > d1.
/// Fermat's principle: minimise t(xc) over the crossing point by golden-section
/// search (t is convex in xc).
pub fn two_layer_time(x: f64, x0: f64, z0: f64, d1: f64, c1: f64, c2: f64) -> (f64, f64) {
    let t = |xc: f64| ((xc - x).hypot(d1)) / c1 + ((x0 - xc).hypot(z0 - d1)) / c2;
    let (mut a, mut b) = (x.min(x0), x.max(x0));
    let g = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = b - g * (b - a);
    let mut d = a + g * (b - a);
    for _ in 0..200 {
        if t(c) < t(d) { b = d } else { a = c }
        c = b - g * (b - a);
        d = a + g * (b - a);
        if (b - a).abs() < 1e-15 { break; }
    }
    let xc = 0.5 * (a + b);
    (t(xc), (xc - x) / (xc - x).hypot(d1))
}

pub struct Scene {
    pub fs: f64,
    pub nt: usize,
    pub nx: usize,
    pub x_step: f64,
    pub t_delay: f64,
    pub pulse: Pulse,
    pub s0: f64,
}

impl Scene {
    /// Column-major (t fastest) samples for a single-layer point scatterer.
    pub fn single_layer(&self, c: f64, x0: f64, z0: f64) -> Vec<f64> {
        let mut v = vec![0.0; self.nt * self.nx];
        for ix in 0..self.nx {
            let x = ix as f64 * self.x_step;
            let r = (x - x0).hypot(z0);
            let tau = 2.0 * r / c;
            let w = directivity((x - x0) / r, self.s0);
            for it in 0..self.nt {
                let t = self.t_delay + it as f64 / self.fs;
                v[it + self.nt * ix] = w * self.pulse.eval(t - tau);
            }
        }
        v
    }

    /// Two-layer point scatterer (interface at depth d1, scatterer in layer 2).
    pub fn two_layer(&self, c1: f64, c2: f64, d1: f64, x0: f64, z0: f64) -> Vec<f64> {
        let mut v = vec![0.0; self.nt * self.nx];
        for ix in 0..self.nx {
            let x = ix as f64 * self.x_step;
            let (t1, s) = two_layer_time(x, x0, z0, d1, c1, c2);
            let tau = 2.0 * t1;
            let w = directivity(s, self.s0);
            for it in 0..self.nt {
                let t = self.t_delay + it as f64 / self.fs;
                v[it + self.nt * ix] = w * self.pulse.eval(t - tau);
            }
        }
        v
    }
}
