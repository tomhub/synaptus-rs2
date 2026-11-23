//! T5 (TEST_PLAN.md): dump inputs and PSM outputs for the high-precision
//! reference in `verify/psm_rounding_ref.py`.
//! Usage: cargo run --release --example psm_rounding -- <out_dir>

#[path = "../tests/common/synth.rs"]
mod synth;

use std::io::Write;
use synaptus::kernels::{KzFormula, Summation};
use synaptus::psm::{debug_spectrum, psm, PhaseUpdate, PsmOptions};
use synaptus::{Acquisition, Medium, PulseEchoData};

fn write_f64(path: &std::path::Path, v: &[f64]) {
    let mut f = std::fs::File::create(path).unwrap();
    for x in v {
        f.write_all(&x.to_le_bytes()).unwrap();
    }
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).expect("out dir"));
    std::fs::create_dir_all(&out).unwrap();
    let s = synth::Scene {
        fs: 20e6, nt: 512, nx: 48, x_step: 0.25e-3, t_delay: 30e-6,
        pulse: synth::Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2,
    };
    let c = 1480.0;
    let data = PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(c, 6.1e-3, 35.3e-3)).unwrap();
    let acq = Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step };
    let medium = Medium::single(c, 0.2); // deep layer => ~840 recursive steps
    let (omega, kx, _ky, w) = debug_spectrum(&data, &acq, &medium, &PsmOptions::legacy()).unwrap();
    write_f64(&out.join("omega.f64"), &omega);
    write_f64(&out.join("kx.f64"), &kx);
    let wv: Vec<f64> = w.iter().flat_map(|z| [z.re, z.im]).collect();
    write_f64(&out.join("W.f64"), &wv);
    let dz = (c / 2.0) / acq.f_high;
    let z_offset = acq.t_delay * (c / 2.0);
    let mut meta = std::fs::File::create(out.join("meta.txt")).unwrap();
    writeln!(meta, "c {c:.17e}\ndz {dz:.17e}\nz_offset {z_offset:.17e}\nnx {}\nnfx {}\nnw {}", s.nx, kx.len(), omega.len()).unwrap();

    let mut variants = vec![];
    for kz in [KzFormula::Legacy, KzFormula::Factored] {
        for phase in [PhaseUpdate::Recursive, PhaseUpdate::Direct, PhaseUpdate::Reanchored(16)] {
            for sum in [Summation::Naive, Summation::Neumaier] {
                variants.push((kz, phase, sum));
            }
        }
    }
    let mut names = std::fs::File::create(out.join("variants.txt")).unwrap();
    for (i, (kz, phase, sum)) in variants.iter().enumerate() {
        let o = PsmOptions { kz: *kz, phase: *phase, sum: *sum, ..PsmOptions::legacy() };
        let t = std::time::Instant::now();
        let img = psm(&data, &acq, &medium, &o).unwrap();
        let secs = t.elapsed().as_secs_f64();
        let l = &img.layers[0];
        let v: Vec<f64> = l.data.iter().flat_map(|z| [z.re, z.im]).collect();
        write_f64(&out.join(format!("img{i}.f64")), &v);
        writeln!(names, "{i} {kz:?} {phase:?} {sum:?} nz={} time={secs:.3}", l.nz).unwrap();
    }
    println!("wrote {} variants to {}", variants.len(), out.display());
}
