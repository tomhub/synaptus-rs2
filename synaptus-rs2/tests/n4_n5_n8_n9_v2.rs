//! N4, N5, N8, N9 (TEST_PLAN_V2.md).
mod common;
use common::synth::{Pulse, Scene};
use std::path::Path;
use synaptus::psm::{psm, PsmOptions};
use synaptus::psm_nufft::{psm_nufft, PsmNufftOptions};
use synaptus::psm_sq::{psm_sq, PsmSqOptions};
use synaptus::{io, Acquisition, FocusedImage, Medium, PulseEchoData};

fn scene(nt: usize, t_delay: f64) -> Scene {
    Scene { fs: 20e6, nt, nx: 256, x_step: 0.25e-3, t_delay, pulse: Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2 }
}
fn acq(s: &Scene) -> Acquisition {
    Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step }
}
/// T2 / T3 synthetic sets, exactly as v1 tests.
fn synth_t2() -> (PulseEchoData, Acquisition, Medium, usize, f64, f64) {
    let (c, x0, z0) = (1480.0, 32.1e-3, 40.3e-3);
    let s = scene(1024, 30e-6);
    (PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(c, x0, z0)).unwrap(), acq(&s), Medium::single(c, 60e-3), 0, z0, x0)
}
fn synth_t3() -> (PulseEchoData, Acquisition, Medium, usize, f64, f64) {
    let (c1, c2, d1, x0, z0) = (1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3);
    let s = scene(1000, 20e-6);
    (PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(c1, c2, d1, x0, z0)).unwrap(), acq(&s),
     Medium { velocities: vec![c1, c2], thicknesses: vec![d1, 25e-3] }, 1, z0, x0)
}

fn rel_l2_layers(a: &FocusedImage, b: &FocusedImage) -> Vec<f64> {
    a.layers.iter().zip(&b.layers).map(|(la, lb)| {
        assert_eq!(la.z, lb.z, "depth vectors differ");
        let num: f64 = la.data.iter().zip(&lb.data).map(|(x, y)| (x - y).norm_sqr()).sum();
        let den: f64 = lb.data.iter().map(|y| y.norm_sqr()).sum();
        (num / den).sqrt()
    }).collect()
}

#[test]
fn n4_nufft_psm_equals_v1_psm() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");
    let mut fails = vec![];
    for (name, m) in [("wire", 1), ("pins", 1), ("plane", 1), ("pins", 4)] {
        let (d, a, med) = io::load_dataset(&dir, name).unwrap();
        let v1 = psm(&d, &a, &med, &PsmOptions { t_fft_mult: m, ..PsmOptions::high_precision_kernels() }).unwrap();
        let v2 = psm_nufft(&d, &a, &med, &PsmNufftOptions { t_fft_mult: m, ..Default::default() }).unwrap();
        for (l, r) in rel_l2_layers(&v2, &v1).iter().enumerate() {
            println!("N4 {name:5} m={m} layer {}: rel L2 = {r:.3e}", l + 1);
            if !(*r <= 1e-10) { fails.push(format!("{name} m={m} L{}: {r:.3e}", l + 1)); }
        }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn n5_sq_self_convergence_and_kappa_choice() {
    let kappas = [1.0, 1.5, 2.0, 3.0, 4.0];
    let mut ok_per_kappa = vec![true; kappas.len()];
    for (label, (d, a, m, _, _, _)) in [("T2", synth_t2()), ("T3", synth_t3())] {
        let reference = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: 6.0, ..Default::default() }).unwrap();
        for (i, &k) in kappas.iter().enumerate() {
            let img = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: k, ..Default::default() }).unwrap();
            let r = rel_l2_layers(&img, &reference);
            let worst = r.iter().cloned().fold(0.0, f64::max);
            println!("N5 {label} kappa={k:3.1}: rel L2 to kappa=6 per layer {:?} (worst {worst:.3e})", r.iter().map(|v| format!("{v:.2e}")).collect::<Vec<_>>());
            if !(worst <= 1e-9) { ok_per_kappa[i] = false; }
        }
    }
    let chosen = kappas.iter().zip(&ok_per_kappa).find(|(_, ok)| **ok).map(|(k, _)| *k);
    println!("N5 chosen default kappa (smallest passing): {chosen:?}");
    assert!(chosen.is_some(), "no kappa <= 4 reached 1e-9 against kappa = 6");
}

fn peak(img: &FocusedImage, layer: usize) -> (f64, f64, f64, f64) {
    let li = &img.layers[layer];
    let (mut best, mut bz, mut bx) = (-1.0, 0, 0);
    for iz in 0..li.nz { for ix in 0..li.nx {
        let v = li.at(iz, ix, 0).norm();
        if v > best { (best, bz, bx) = (v, iz, ix); }
    }}
    (li.z[bz], img.x[bx], li.z[1] - li.z[0], img.x[1] - img.x[0])
}

#[test]
fn n8_point_location() {
    let mut fails = vec![];
    for (label, (d, a, m, layer, z0, x0)) in [("T2", synth_t2()), ("T3", synth_t3())] {
        for (name, img) in [
            ("nufft-psm", psm_nufft(&d, &a, &m, &PsmNufftOptions::default()).unwrap()),
            ("sq-psm", psm_sq(&d, &a, &m, &PsmSqOptions::default()).unwrap()),
        ] {
            let (z, x, dz, dx) = peak(&img, layer);
            let (ez, ex) = ((z - z0).abs() / dz, (x - x0).abs() / dx);
            println!("N8 {label} {name:9}: z err {ez:.2} px, x err {ex:.2} px");
            if !(ez <= 1.0 && ex <= 1.0) { fails.push(format!("{label} {name}: {ez:.2}, {ex:.2}")); }
        }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn n9_invalid_inputs_are_errors() {
    let ok_data = || PulseEchoData::new(64, 4, 1, vec![0.0; 256]).unwrap();
    let ok_acq = Acquisition { fs: 20e6, t_delay: 0.0, f_low: 1e6, f_high: 3e6, x_step: 1e-3, y_step: 1e-3 };
    let ok_med = Medium::single(1480.0, 0.01);
    let run = |d: &PulseEchoData, a: &Acquisition, m: &Medium| {
        [psm_nufft(d, a, m, &PsmNufftOptions::default()).is_err(), psm_sq(d, a, m, &PsmSqOptions::default()).is_err()]
    };
    let cases: Vec<(&str, Acquisition, Medium)> = vec![
        ("f_low>=f_high", Acquisition { f_low: 3e6, ..ok_acq }, ok_med.clone()),
        ("f_high>Nyquist", Acquisition { f_high: 11e6, ..ok_acq }, ok_med.clone()),
        ("x_step<=0", Acquisition { x_step: 0.0, ..ok_acq }, ok_med.clone()),
        ("fs NaN", Acquisition { fs: f64::NAN, ..ok_acq }, ok_med.clone()),
        ("velocity<=0", ok_acq, Medium::single(-1.0, 0.01)),
        ("mismatched layers", ok_acq, Medium { velocities: vec![1480.0, 2000.0], thicknesses: vec![0.01] }),
        ("offset beyond layer", Acquisition { t_delay: 1e-3, ..ok_acq }, ok_med.clone()),
    ];
    for (name, a, m) in cases {
        let r = run(&ok_data(), &a, &m);
        assert!(r.iter().all(|&e| e), "{name}: not all v2 algorithms returned Err: {r:?}");
    }
    assert!(psm_sq(&ok_data(), &ok_acq, &ok_med, &PsmSqOptions { hh: Some(vec![1.0]), ..Default::default() }).is_err(), "hh accepted");
    assert!(psm_sq(&ok_data(), &ok_acq, &ok_med, &PsmSqOptions { kappa: 0.0, ..Default::default() }).is_err(), "kappa 0 accepted");
    let r = run(&ok_data(), &ok_acq, &ok_med);
    assert!(r.iter().all(|&e| !e), "valid input rejected: {r:?}");
}
