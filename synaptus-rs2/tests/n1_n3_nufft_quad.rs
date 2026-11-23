//! N1–N3 (TEST_PLAN_V2.md): NUFFT types 1/2 vs direct sums; Gauss–Legendre rules.
//! Oracles are direct O(MN) sums with Neumaier summation, written here.
use num_complex::Complex64 as C64;
use synaptus::nufft::Nufft1d;
use synaptus::quad::gauss_legendre;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        // splitmix64 → uniform [0,1)
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn sym(&mut self) -> f64 { 2.0 * self.next() - 1.0 }
}

fn neumaier(it: impl Iterator<Item = C64>) -> C64 {
    let (mut s, mut c) = (C64::new(0.0, 0.0), C64::new(0.0, 0.0));
    for v in it {
        for (sp, cp, vp) in [(&mut s.re, &mut c.re, v.re), (&mut s.im, &mut c.im, v.im)] {
            let t = *sp + vp;
            if sp.abs() >= vp.abs() { *cp += (*sp - t) + vp } else { *cp += (vp - t) + *sp }
            *sp = t;
        }
    }
    s + c
}

const CASES: [(usize, usize); 3] = [(100, 100), (1000, 200), (5000, 1000)];
const TOLS: [f64; 3] = [1e-6, 1e-9, 1e-12];

#[test]
fn n1_type1_vs_direct() {
    let mut fails = vec![];
    for (ci, &(m, n)) in CASES.iter().enumerate() {
        let mut rng = Rng(1000 + ci as u64);
        let x: Vec<f64> = (0..m).map(|_| std::f64::consts::PI * rng.sym()).collect();
        let c: Vec<C64> = (0..m).map(|_| C64::new(rng.sym(), rng.sym())).collect();
        let half = n as i64 / 2;
        let exact: Vec<C64> = (-half..half).map(|k| neumaier(x.iter().zip(&c).map(|(&xj, &cj)| cj * C64::from_polar(1.0, k as f64 * xj)))).collect();
        let c1: f64 = c.iter().map(|v| v.norm()).sum();
        for &eps in &TOLS {
            let f = Nufft1d::new(n, eps).unwrap().type1(&x, &c);
            let maxe = f.iter().zip(&exact).map(|(a, b)| (a - b).norm()).fold(0.0, f64::max);
            let rel2 = (f.iter().zip(&exact).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>() / exact.iter().map(|b| b.norm_sqr()).sum::<f64>()).sqrt();
            let bound = 10.0 * eps * c1 + 1e-13 * c1;
            println!("N1 M={m:5} N={n:5} eps={eps:.0e}: max|err|/||c||1 = {:.2e} (limit {:.2e}), rel L2 = {rel2:.2e}", maxe / c1, bound / c1);
            if !(maxe <= bound) { fails.push(format!("M={m} N={n} eps={eps:e}: {:.2e}", maxe / c1)); }
        }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn n2_type2_vs_direct() {
    let mut fails = vec![];
    for (ci, &(m, n)) in CASES.iter().enumerate() {
        let mut rng = Rng(2000 + ci as u64);
        let x: Vec<f64> = (0..m).map(|_| std::f64::consts::PI * rng.sym()).collect();
        let f: Vec<C64> = (0..n).map(|_| C64::new(rng.sym(), rng.sym())).collect();
        let half = n as i64 / 2;
        let exact: Vec<C64> = x.iter().map(|&xj| neumaier((-half..half).zip(&f).map(|(k, &fk)| fk * C64::from_polar(1.0, -(k as f64) * xj)))).collect();
        let f1: f64 = f.iter().map(|v| v.norm()).sum();
        for &eps in &TOLS {
            let c = Nufft1d::new(n, eps).unwrap().type2(&f, &x);
            let maxe = c.iter().zip(&exact).map(|(a, b)| (a - b).norm()).fold(0.0, f64::max);
            let rel2 = (c.iter().zip(&exact).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>() / exact.iter().map(|b| b.norm_sqr()).sum::<f64>()).sqrt();
            let bound = 10.0 * eps * f1 + 1e-13 * f1;
            println!("N2 M={m:5} N={n:5} eps={eps:.0e}: max|err|/||f||1 = {:.2e} (limit {:.2e}), rel L2 = {rel2:.2e}", maxe / f1, bound / f1);
            if !(maxe <= bound) { fails.push(format!("M={m} N={n} eps={eps:e}: {:.2e}", maxe / f1)); }
        }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn n3_gauss_legendre() {
    let mut worst = (0.0f64, 0.0f64);
    for n in 1..=200usize {
        let (x, w) = gauss_legendre(n);
        let sw: f64 = w.iter().sum();
        assert!((sw - 2.0).abs() <= 1e-14 * n as f64, "n={n}: sum w = {sw}");
        assert!(x.windows(2).all(|p| p[0] < p[1]) && x[0] > -1.0 && x[n - 1] < 1.0, "n={n}: nodes not increasing in (-1,1)");
        let mut k2 = 0;
        while k2 <= 2 * n - 1 {
            let exact = 2.0 / (k2 + 1) as f64;
            let q: f64 = x.iter().zip(&w).map(|(xi, wi)| wi * xi.powi(k2 as i32)).sum();
            let rel = (q - exact).abs() / exact;
            worst.0 = worst.0.max(rel / n as f64);
            assert!(rel <= 1e-13 * n as f64, "n={n}, x^{k2}: rel err {rel:.2e}");
            k2 += 2;
        }
        worst.1 = worst.1.max((sw - 2.0).abs());
    }
    println!("N3 worst rel moment error / n = {:.2e}, worst |sum w - 2| = {:.2e}", worst.0, worst.1);
}
