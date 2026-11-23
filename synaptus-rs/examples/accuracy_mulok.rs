//! T4 (TEST_PLAN.md): MULOK vs PSM on the same depth grid.
//! PSM evaluates eq. 3.21 exactly on the sampled ω grid (no interpolation), and
//! thesis-MULOK is an exact change of variables of the same integral, so
//! ‖MULOK − PSM‖/‖PSM‖ measures MULOK's implementation error.
//! Usage: cargo run --release --example accuracy_mulok [dataset ...]

#[path = "../tests/common/synth.rs"]
mod synth;

use std::path::Path;
use std::time::Instant;
use synaptus::fftutil::C64;
use synaptus::mulok::{mulok, Interp, Jacobian, KzSupport, MulokOptions};
use synaptus::psm::{psm, DepthStep, PsmOptions};
use synaptus::{io, legacy, Acquisition, FocusedImage, Medium, PulseEchoData};

fn flatten(img: &FocusedImage) -> Vec<C64> {
    img.layers.iter().flat_map(|l| l.data.iter().copied()).collect()
}

/// (raw relative error, relative error after least-squares complex scale fit)
fn errors(reference: &[C64], x: &[C64]) -> (f64, f64) {
    assert_eq!(reference.len(), x.len(), "grid mismatch");
    let nr: f64 = reference.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
    let raw = reference.iter().zip(x).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>().sqrt() / nr;
    let num: C64 = x.iter().zip(reference).map(|(m, r)| m.conj() * r).sum();
    let den: f64 = x.iter().map(|m| m.norm_sqr()).sum();
    let alpha = num / den;
    let fit = reference.iter().zip(x).map(|(a, b)| (a - alpha * b).norm_sqr()).sum::<f64>().sqrt() / nr;
    (raw, fit)
}

/// Envelope metric (Amendment A2): (raw, real-scale-fitted) ‖|x|−α|ref|‖/‖|ref|‖.
fn env_errors(reference: &[C64], x: &[C64]) -> (f64, f64) {
    let r: Vec<f64> = reference.iter().map(|v| v.norm()).collect();
    let m: Vec<f64> = x.iter().map(|v| v.norm()).collect();
    let nr = r.iter().map(|v| v * v).sum::<f64>().sqrt();
    let raw = r.iter().zip(&m).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt() / nr;
    let alpha = m.iter().zip(&r).map(|(a, b)| a * b).sum::<f64>() / m.iter().map(|a| a * a).sum::<f64>();
    let fit = r.iter().zip(&m).map(|(a, b)| (a - alpha * b).powi(2)).sum::<f64>().sqrt() / nr;
    (raw, fit)
}

fn load(name: &str) -> (PulseEchoData, Acquisition, Medium) {
    let mk_scene = |nt, td| synth::Scene {
        fs: 20e6, nt, nx: 256, x_step: 0.25e-3, t_delay: td,
        pulse: synth::Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2,
    };
    let acq = |s: &synth::Scene| Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step };
    match name {
        "synth1" => {
            let s = mk_scene(1024, 30e-6);
            let d = PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(1480.0, 32.1e-3, 40.3e-3)).unwrap();
            (d, acq(&s), Medium::single(1480.0, 60e-3))
        }
        "synth2" => {
            let s = mk_scene(1000, 20e-6);
            let d = PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3)).unwrap();
            (d, acq(&s), Medium { velocities: vec![1480.0, 2730.0], thicknesses: vec![30e-3, 25e-3] })
        }
        other => io::load_dataset(&Path::new(env!("CARGO_MANIFEST_DIR")).join("data"), other).unwrap(),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(|s| s.as_str()) == Some("--amendment-a2") {
        return amendment_a2(&args[1..]);
    }
    let sets: Vec<String> = if args.is_empty() {
        ["synth1", "synth2", "wire", "pins", "plane"].iter().map(|s| s.to_string()).collect()
    } else { args };
    println!("dataset | configuration | raw rel.err | scale-fitted rel.err | time [s]");
    println!("---|---|---|---|---");
    for name in &sets {
        let (d, a, m) = load(name);
        let ref_hf = flatten(&psm(&d, &a, &m, &PsmOptions::default()).unwrap());
        let ref_bw = flatten(&psm(&d, &a, &m, &PsmOptions { depth_step: DepthStep::Bandwidth, ..PsmOptions::default() }).unwrap());
        let row = |label: &str, img: FocusedImage, reference: &[C64], t: f64| {
            let (raw, fit) = errors(reference, &flatten(&img));
            println!("{name} | {label} | {raw:.3e} | {fit:.3e} | {t:.2}");
        };
        let t = Instant::now();
        let lg = legacy::mulok(&d, &a, &m, &Default::default()).unwrap();
        row("legacy mulok.m (vs PSM on its dz grid)", lg, &ref_bw, t.elapsed().as_secs_f64());

        let base = MulokOptions::thesis();
        let attrib: Vec<(&str, MulokOptions)> = vec![
            ("thesis: linear, up=4", base.clone()),
            ("thesis + legacy squared Jacobian", MulokOptions { jacobian: Jacobian::LegacySquared, ..base.clone() }),
            ("thesis + legacy kz truncation", MulokOptions { kz_support: KzSupport::TruncatedBelowKLow, ..base.clone() }),
            ("thesis + legacy -0.75 bin ω offset", MulokOptions { omega_shift_bins: -0.75, ..base.clone() }),
            ("thesis + all three legacy defects", MulokOptions {
                jacobian: Jacobian::LegacySquared, kz_support: KzSupport::TruncatedBelowKLow, omega_shift_bins: -0.75, ..base.clone() }),
        ];
        for (label, o) in attrib {
            let t = Instant::now();
            let img = mulok(&d, &a, &m, &o).unwrap();
            row(label, img, &ref_hf, t.elapsed().as_secs_f64());
        }
        for demod in [false, true] {
            for interp in [Interp::Linear, Interp::Cubic, Interp::Lanczos(6)] {
                for up in [1usize, 2, 4, 8] {
                    let o = MulokOptions { interp, upsample: up, demodulate: demod, ..base.clone() };
                    let t = Instant::now();
                    let img = mulok(&d, &a, &m, &o).unwrap();
                    row(&format!("interp={interp:?} up={up} demod={demod}"), img, &ref_hf, t.elapsed().as_secs_f64());
                }
            }
        }
    }
}

/// Amendment A2: envelope metric, attribution on an accurate baseline, PSM self-check.
fn amendment_a2(args: &[String]) {
    let sets: Vec<String> = if args.is_empty() {
        ["synth1", "synth2", "wire", "pins", "plane"].iter().map(|s| s.to_string()).collect()
    } else { args.to_vec() };
    println!("dataset | configuration | complex raw | complex fit | envelope raw | envelope fit");
    println!("---|---|---|---|---|---");
    for name in &sets {
        let (d, a, m) = load(name);
        let ref_hf = flatten(&psm(&d, &a, &m, &PsmOptions::default()).unwrap());
        let ref_bw = flatten(&psm(&d, &a, &m, &PsmOptions { depth_step: DepthStep::Bandwidth, ..PsmOptions::default() }).unwrap());
        let row = |label: &str, x: &[C64], reference: &[C64]| {
            let (cr, cf) = errors(reference, x);
            let (er, ef) = env_errors(reference, x);
            println!("{name} | {label} | {cr:.3e} | {cf:.3e} | {er:.3e} | {ef:.3e}");
        };
        // PSM self-check: quadrature (ω-grid) discretisation of the reference itself
        let psm2: Vec<C64> = flatten(&psm(&d, &a, &m, &PsmOptions { t_fft_mult: 2, ..PsmOptions::default() }).unwrap())
            .into_iter().map(|v| v * 0.5).collect();
        row("PSM self-check: nfft_t x2 (scaled 1/2) vs x1", &psm2, &ref_hf);
        row("legacy mulok.m (vs PSM, its dz grid)", &flatten(&legacy::mulok(&d, &a, &m, &Default::default()).unwrap()), &ref_bw);
        let base = MulokOptions { interp: Interp::Cubic, upsample: 4, demodulate: true, ..MulokOptions::thesis() };
        let variants: Vec<(&str, MulokOptions)> = vec![
            ("accurate baseline (cubic, up=4, demod)", base.clone()),
            ("baseline + legacy squared Jacobian", MulokOptions { jacobian: Jacobian::LegacySquared, ..base.clone() }),
            ("baseline + legacy kz truncation", MulokOptions { kz_support: KzSupport::TruncatedBelowKLow, ..base.clone() }),
            ("baseline + legacy -0.75 bin ω offset", MulokOptions { omega_shift_bins: -0.75, ..base.clone() }),
            ("baseline + all three legacy defects", MulokOptions {
                jacobian: Jacobian::LegacySquared, kz_support: KzSupport::TruncatedBelowKLow, omega_shift_bins: -0.75, ..base.clone() }),
            ("thesis as coded in mulok.m interp (linear, up=4, no demod)", MulokOptions::thesis()),
        ];
        for (label, o) in variants {
            row(label, &flatten(&mulok(&d, &a, &m, &o).unwrap()), &ref_hf);
        }
    }
}
