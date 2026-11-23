//! T2, T3 (analytic focus location) and T6 (input validation) from TEST_PLAN.md.

mod common;
use common::synth::{Pulse, Scene};
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::psm::{psm, PsmOptions};
use synaptus::{legacy, Acquisition, FocusedImage, Medium, PulseEchoData};

fn scene(nt: usize, t_delay: f64) -> Scene {
    Scene { fs: 20e6, nt, nx: 256, x_step: 0.25e-3, t_delay, pulse: Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2 }
}

fn acq(s: &Scene) -> Acquisition {
    Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step }
}

/// Location (z, x) of max |image| in a given layer.
fn peak(img: &FocusedImage, layer: usize) -> (f64, f64, f64, f64) {
    let li = &img.layers[layer];
    let (mut best, mut bz, mut bx) = (-1.0, 0, 0);
    for iz in 0..li.nz {
        for ix in 0..li.nx {
            let v = li.at(iz, ix, 0).norm();
            if v > best { (best, bz, bx) = (v, iz, ix); }
        }
    }
    let dz = li.z[1] - li.z[0];
    (li.z[bz], img.x[bx], dz, img.x[1] - img.x[0])
}

fn check(name: &str, img: &FocusedImage, layer: usize, z0: f64, x0: f64) -> Option<String> {
    let (z, x, dz, dx) = peak(img, layer);
    let (ez, ex) = ((z - z0).abs() / dz, (x - x0).abs() / dx);
    println!("{name:28} peak z={:.4} mm (err {ez:.2} px)  x={:.4} mm (err {ex:.2} px)", z * 1e3, x * 1e3);
    if ez <= 1.0 && ex <= 1.0 { None } else { Some(format!("{name}: z err {ez:.2} px, x err {ex:.2} px")) }
}

#[test]
fn t2_single_layer_point_scatterer() {
    let (c, x0, z0) = (1480.0, 32.1e-3, 40.3e-3);
    let s = scene(1024, 30e-6);
    let data = PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(c, x0, z0)).unwrap();
    let a = acq(&s);
    let medium = Medium::single(c, 60e-3);
    let mut fails = vec![];
    fails.extend(check("legacy psm", &legacy::psm(&data, &a, &medium).unwrap(), 0, z0, x0));
    fails.extend(check("psm accurate", &psm(&data, &a, &medium, &PsmOptions::default()).unwrap(), 0, z0, x0));
    fails.extend(check("legacy mulok", &legacy::mulok(&data, &a, &medium, &Default::default()).unwrap(), 0, z0, x0));
    fails.extend(check("mulok thesis", &mulok(&data, &a, &medium, &MulokOptions::thesis()).unwrap(), 0, z0, x0));
    fails.extend(check("mulok accurate", &mulok(&data, &a, &medium, &MulokOptions::accurate()).unwrap(), 0, z0, x0));
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn t3_two_layer_point_scatterer() {
    let (c1, c2, d1, x0, z0) = (1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3);
    let s = scene(1000, 20e-6);
    let data = PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(c1, c2, d1, x0, z0)).unwrap();
    let a = acq(&s);
    let medium = Medium { velocities: vec![c1, c2], thicknesses: vec![d1, 25e-3] };
    let mut fails = vec![];
    fails.extend(check("legacy psm (2 layers)", &legacy::psm(&data, &a, &medium).unwrap(), 1, z0, x0));
    fails.extend(check("psm accurate (2 layers)", &psm(&data, &a, &medium, &PsmOptions::default()).unwrap(), 1, z0, x0));
    fails.extend(check("mulok thesis (2 layers)", &mulok(&data, &a, &medium, &MulokOptions::thesis()).unwrap(), 1, z0, x0));
    fails.extend(check("mulok accurate (2 layers)", &mulok(&data, &a, &medium, &MulokOptions::accurate()).unwrap(), 1, z0, x0));
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn t6_invalid_inputs_are_errors() {
    let ok_data = || PulseEchoData::new(8, 4, 1, vec![0.0; 32]).unwrap();
    let ok_acq = Acquisition { fs: 20e6, t_delay: 0.0, f_low: 1e6, f_high: 3e6, x_step: 1e-3, y_step: 1e-3 };
    let ok_med = Medium::single(1480.0, 0.01);
    let run = |d: &PulseEchoData, a: &Acquisition, m: &Medium| {
        [
            psm(d, a, m, &PsmOptions::default()).is_err(),
            legacy::psm(d, a, m).is_err(),
            legacy::mulok(d, a, m, &Default::default()).is_err(),
            mulok(d, a, m, &MulokOptions::accurate()).is_err(),
        ]
    };
    assert!(PulseEchoData::new(0, 4, 1, vec![]).is_err(), "empty");
    assert!(PulseEchoData::new(8, 4, 1, vec![0.0; 31]).is_err(), "length");
    let mut nanv = vec![0.0; 32];
    nanv[5] = f64::NAN;
    assert!(PulseEchoData::new(8, 4, 1, nanv).is_err(), "NaN");
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
        assert!(r.iter().all(|&e| e), "{name}: not all algorithms returned Err: {r:?}");
    }
    // sanity: a realistically sized valid input succeeds for every algorithm.
    // (The 8-sample input above has a single in-band FFT bin, which legacy MULOK
    // correctly rejects — MATLAB's interp1 needs >= 2 points.)
    let big = PulseEchoData::new(256, 4, 1, vec![0.0; 1024]).unwrap();
    let r = run(&big, &ok_acq, &Medium::single(1480.0, 0.005));
    assert!(r.iter().all(|&e| !e), "valid input rejected: {r:?}");
}
