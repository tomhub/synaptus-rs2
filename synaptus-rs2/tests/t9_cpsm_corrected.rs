//! T9 (TEST_PLAN amendment A3): corrected CPSM on synthetic cylindrical point
//! scatterers with analytically known positions. Data generator written here,
//! independent of the library (plain geometry + pulse).
use std::f64::consts::PI;
use synaptus::cpsm::{cpsm, CpsmOptions, CylAcquisition, CylImage, TransferFunction};
use synaptus::{legacy, PulseEchoData};

const R0: f64 = 0.058;
const C: f64 = 1480.0;
const FS: f64 = 25e6;
const NT: usize = 1024;
const NPHI: usize = 96;
const NZ: usize = 64;
const DZ: f64 = 0.8e-3;

fn dphi() -> f64 { 2.0 * PI / 720.0 }
fn phi_of(i: usize) -> f64 { (i as f64 - (NPHI as f64 - 1.0) / 2.0) * dphi() }

fn pulse(t: f64) -> f64 {
    let (fc, s) = (2e6, 0.6e-6);
    (-(t * t) / (2.0 * s * s)).exp() * (2.0 * PI * fc * t).cos()
}

/// Scatterers (r, φ, z); returns data (t fastest, then φ, then z) and t_delay.
fn synth(sc: &[(f64, f64, f64)]) -> (Vec<f64>, f64) {
    let t_delay = 2.0 * (0.065 - R0) / C;
    let mut v = vec![0.0; NT * NPHI * NZ];
    for iz in 0..NZ { for ip in 0..NPHI {
        let (pt, zt) = (phi_of(ip), iz as f64 * DZ);
        let (xt, yt) = (R0 * pt.cos(), R0 * pt.sin());
        for &(rs, ps, zs) in sc {
            let (dx, dy, dzz) = (rs * ps.cos() - xt, rs * ps.sin() - yt, zs - zt);
            let d = (dx * dx + dy * dy + dzz * dzz).sqrt();
            // angle to the outward normal (cos pt, sin pt, 0)
            let cos_th = (dx * pt.cos() + dy * pt.sin()) / d;
            let w = (-((1.0 - cos_th * cos_th).max(0.0) / 0.04)).exp(); // exp(-(sinθ/0.2)²)
            let tau = 2.0 * d / C;
            for it in 0..NT {
                v[it + NT * (ip + NPHI * iz)] += w * pulse(t_delay + it as f64 / FS - tau);
            }
        }
    }}
    (v, t_delay)
}

fn peak_near(img: &CylImage, rr: &[f64], target: (f64, f64, f64)) -> ((f64, f64, f64), f64) {
    let ir0 = rr.iter().enumerate().min_by(|a, b| (a.1 - target.0).abs().partial_cmp(&(b.1 - target.0).abs()).unwrap()).unwrap().0;
    let ip0 = img.phi.iter().enumerate().min_by(|a, b| (a.1 - target.1).abs().partial_cmp(&(b.1 - target.1).abs()).unwrap()).unwrap().0;
    let iz0 = (target.2 / DZ).round() as usize;
    let (mut best, mut at) = (-1.0, (0, 0, 0));
    let rng = |c: usize, n: usize| c.saturating_sub(3)..(c + 4).min(n);
    for ir in rng(ir0, rr.len()) { for ip in rng(ip0, NPHI) { for iz in rng(iz0, NZ) {
        let v = img.at(ir, ip, iz).norm();
        if v > best { best = v; at = (ir, ip, iz); }
    }}}
    let dr = rr[1] - rr[0];
    let err = ((rr[at.0] - target.0).abs() / dr, (img.phi[at.1] - target.1).abs() / dphi(), (img.z[at.2] - target.2).abs() / DZ);
    // focus: fraction of total energy within ±2 samples of the true position
    let tot: f64 = img.data.iter().map(|c| c.norm_sqr()).sum();
    let mut near = 0.0;
    let rng2 = |c: usize, n: usize| c.saturating_sub(2)..(c + 3).min(n);
    for ir in rng2(ir0, rr.len()) { for ip in rng2(ip0, NPHI) { for iz in rng2(iz0, NZ) { near += img.at(ir, ip, iz).norm_sqr(); }}}
    (err, near / tot)
}

#[test]
fn t9_corrected_cpsm_locates_cylindrical_scatterers() {
    let sc = [(0.075, 3.0f64.to_radians(), 24.0e-3), (0.082, (-5.0f64).to_radians(), 30.4e-3)];
    let (v, t_delay) = synth(&sc);
    let data = PulseEchoData::new(NT, NPHI, NZ, v).unwrap();
    let acq = CylAcquisition { fs: FS, t_delay, f_low: 0.0, f_high: 4e6, phi_step: dphi(), z_step: DZ, r0: R0, c: C };
    let r_start = R0 + t_delay * C / 2.0;
    let r_step = (C / 2.0) / 4e6;
    let rr: Vec<f64> = (0..).map(|i| r_start + i as f64 * r_step).take_while(|&r| r <= 0.092).collect();
    let fixed = cpsm(&data, &acq, &rr, &CpsmOptions::default()).unwrap();
    let orig = legacy::cpsm(&data, &acq, &rr, TransferFunction::Haun, 1).unwrap();
    let mut fails = vec![];
    for (k, &s) in sc.iter().enumerate() {
        let (e, f) = peak_near(&fixed, &rr, s);
        let (eo, fo) = peak_near(&orig, &rr, s);
        println!("T9 scatterer {k}: corrected  err (r, phi, z) = ({:.2}, {:.2}, {:.2}) samples, energy within ±2 = {:.3}", e.0, e.1, e.2, f);
        println!("T9 scatterer {k}: legacy     err (r, phi, z) = ({:.2}, {:.2}, {:.2}) samples, energy within ±2 = {:.3}  [report only]", eo.0, eo.1, eo.2, fo);
        if !(e.0 <= 1.0 && e.1 <= 1.0 && e.2 <= 1.0) { fails.push(format!("scatterer {k}: {e:?}")); }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn t9_invalid_cpsm_inputs_are_errors() {
    let d = PulseEchoData::new(64, 8, 1, vec![0.0; 512]).unwrap();
    let ok = CylAcquisition { fs: FS, t_delay: 0.0, f_low: 0.0, f_high: 4e6, phi_step: dphi(), z_step: DZ, r0: R0, c: C };
    let o = CpsmOptions::default();
    assert!(cpsm(&d, &ok, &[0.07], &o).is_ok());
    assert!(cpsm(&d, &CylAcquisition { phi_step: 1.0, ..ok }, &[0.07], &o).is_err(), "arc > 2π accepted");
    assert!(cpsm(&d, &CylAcquisition { f_high: 13e6, ..ok }, &[0.07], &o).is_err(), "f_high > Nyquist accepted");
    assert!(cpsm(&d, &ok, &[], &o).is_err(), "empty ranges accepted");
}
