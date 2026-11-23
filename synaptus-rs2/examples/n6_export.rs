//! N6 export (TEST_PLAN_V2.md): inputs + SQ-PSM values at the pre-declared pixels,
//! for the independent Python oracle `verify/sq_oracle.py`.
//!   cargo run --release --example n6_export   ->  verify/n6/<case>/
#[path = "../tests/common/synth.rs"]
#[allow(dead_code)]
mod synth;

use std::io::Write;
use std::path::Path;
use synaptus::psm::{psm, PsmOptions};
use synaptus::psm_sq::{psm_sq, PsmSqOptions};
use synaptus::{io, Acquisition, Medium, PulseEchoData};
use synth::{Pulse, Scene};

fn export(dir: &Path, d: &PulseEchoData, a: &Acquisition, m: &Medium) {
    std::fs::create_dir_all(dir).unwrap();
    let mut f = std::fs::File::create(dir.join("data.f64")).unwrap();
    f.write_all(&d.samples.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
    let nfx = 1usize << (d.nx as f64).log2().ceil() as u32;
    let sq = psm_sq(d, a, m, &PsmSqOptions::default()).unwrap();
    let v1 = psm(d, a, m, &PsmOptions::default()).unwrap();
    let mut px = vec![];
    for l in 0..m.velocities.len() {
        let (li, lv) = (&sq.layers[l], &v1.layers[l]);
        let (nz, nx) = (li.nz, li.nx);
        let (mut best, mut bz, mut bx) = (-1.0, 0, 0);
        for iz in 0..nz { for ix in 0..nx { let v = lv.at(iz, ix, 0).norm(); if v > best { (best, bz, bx) = (v, iz, ix); } } }
        let lmax = li.data.iter().map(|c| c.norm()).fold(0.0, f64::max);
        for (iz, ix) in [(bz, bx), (nz / 4, nx / 4), (nz / 2, nx / 2), (3 * nz / 4, 3 * nx / 4)] {
            let v = li.at(iz, ix, 0);
            px.push(format!(
                "{{\"layer\": {l}, \"iz\": {iz}, \"ix\": {ix}, \"z\": {:.17e}, \"x\": {:.17e}, \"sq_re\": {:.17e}, \"sq_im\": {:.17e}, \"layer_max\": {lmax:.17e}}}",
                li.z[iz], sq.x[ix], v.re, v.im
            ));
        }
    }
    let fmt = |v: &[f64]| v.iter().map(|x| format!("{x:.17e}")).collect::<Vec<_>>().join(", ");
    let params = format!(
        "{{\"nt\": {}, \"nx\": {}, \"fs\": {:.17e}, \"t_delay\": {:.17e}, \"f_low\": {:.17e}, \"f_high\": {:.17e}, \"x_step\": {:.17e}, \"nfx\": {nfx}, \"velocities\": [{}], \"thicknesses\": [{}], \"kappa\": {}, \"pixels\": [\n  {}\n]}}\n",
        d.nt, d.nx, a.fs, a.t_delay, a.f_low, a.f_high, a.x_step, fmt(&m.velocities), fmt(&m.thicknesses),
        synaptus::psm_sq::DEFAULT_KAPPA, px.join(",\n  ")
    );
    std::fs::write(dir.join("case.json"), params).unwrap();
    println!("wrote {}", dir.display());
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let out = root.join("verify/n6");
    let scene = |nt, t_delay| Scene { fs: 20e6, nt, nx: 256, x_step: 0.25e-3, t_delay, pulse: Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2 };
    let acq = |s: &Scene| Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step };
    // T2 and T3 synthetic sets, identical to the v1/v2 tests
    let s = scene(1024, 30e-6);
    let d = PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(1480.0, 32.1e-3, 40.3e-3)).unwrap();
    export(&out.join("t2"), &d, &acq(&s), &Medium::single(1480.0, 60e-3));
    let s = scene(1000, 20e-6);
    let d = PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3)).unwrap();
    export(&out.join("t3"), &d, &acq(&s), &Medium { velocities: vec![1480.0, 2730.0], thicknesses: vec![30e-3, 25e-3] });
    let (d, a, m) = io::load_dataset(&root.join("data"), "pins").unwrap();
    export(&out.join("pins"), &d, &a, &m);
}
