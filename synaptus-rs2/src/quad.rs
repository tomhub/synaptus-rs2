//! Gauss–Legendre quadrature nodes and weights on [−1, 1].
//!
//! Newton iteration on the three-term recurrence for P_n (classical; e.g.
//! Golub & Welsch 1969 give the eigenvalue alternative). O(n²) per rule, so
//! callers cache rules by n (see [`GaussCache`]).

use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

/// (P_n(x), P_n'(x)) by the recurrence (k+1)P_{k+1} = (2k+1)xP_k − kP_{k−1}.
fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, x);
    if n == 0 {
        return (1.0, 0.0);
    }
    for k in 1..n {
        let p2 = ((2 * k + 1) as f64 * x * p1 - k as f64 * p0) / (k + 1) as f64;
        p0 = p1;
        p1 = p2;
    }
    // P_n'(x) = n (x P_n − P_{n−1}) / (x² − 1)
    let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
    (p1, dp)
}

/// n-point Gauss–Legendre rule on [−1, 1], nodes ascending.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(n >= 1);
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let m = (n + 1) / 2;
    for i in 0..m {
        // Tricomi-type initial guess for the i-th largest root.
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (p, dp) = legendre(n, z);
            let dz = p / dp;
            z -= dz;
            if dz.abs() <= 1e-16 * z.abs().max(1e-300) {
                break;
            }
        }
        let (_, dp) = legendre(n, z);
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        x[n - 1 - i] = z;
        x[i] = -z;
        w[n - 1 - i] = wi;
        w[i] = wi;
    }
    if n % 2 == 1 {
        x[n / 2] = 0.0;
    }
    (x, w)
}

/// Thread-safe cache of Gauss–Legendre rules keyed by n.
#[derive(Default)]
pub struct GaussCache(Mutex<HashMap<usize, Arc<(Vec<f64>, Vec<f64>)>>>);

impl GaussCache {
    pub fn get(&self, n: usize) -> Arc<(Vec<f64>, Vec<f64>)> {
        if let Some(r) = self.0.lock().unwrap().get(&n) {
            return r.clone();
        }
        let r = Arc::new(gauss_legendre(n));
        self.0.lock().unwrap().entry(n).or_insert(r).clone()
    }
}
