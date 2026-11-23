//! PSM ω-quadrature convergence: images with t_fft_mult m = 1, 2, 4, 8 (scaled 1/m).
//! Reports ‖I_m − I_8‖/‖I_8‖ (complex and envelope). Discretisation, not rounding.
#[path = "../tests/common/synth.rs"]
mod synth;
use synaptus::fftutil::C64;
use synaptus::psm::{psm, PsmOptions};
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::{io, Acquisition, Medium, PulseEchoData};

fn load(name: &str) -> (PulseEchoData, Acquisition, Medium) {
    let sc = |nt, td| synth::Scene { fs: 20e6, nt, nx: 256, x_step: 0.25e-3, t_delay: td,
        pulse: synth::Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2 };
    let ac = |s: &synth::Scene| Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step };
    match name {
        "synth1" => { let s = sc(1024, 30e-6); (PulseEchoData::new(s.nt, s.nx, 1, s.single_layer(1480.0, 32.1e-3, 40.3e-3)).unwrap(), ac(&s), Medium::single(1480.0, 60e-3)) }
        "synth2" => { let s = sc(1000, 20e-6); (PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3)).unwrap(), ac(&s),
            Medium { velocities: vec![1480.0, 2730.0], thicknesses: vec![30e-3, 25e-3] }) }
        o => io::load_dataset(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data"), o).unwrap(),
    }
}
fn main() {
    println!("dataset | method (vs PSM m=8 limit) | complex ‖I_m−I_8‖/‖I_8‖ | envelope | time [s]\n---|---|---|---|---");
    for name in ["synth1", "synth2", "wire", "pins", "plane"] {
        let (d, a, m) = load(name);
        let mut imgs = vec![];
        for mult in [1usize, 2, 4, 8] {
            let t = std::time::Instant::now();
            let img = psm(&d, &a, &m, &PsmOptions { t_fft_mult: mult, ..PsmOptions::default() }).unwrap();
            let v: Vec<C64> = img.layers.iter().flat_map(|l| l.data.iter().map(|z| z / mult as f64)).collect();
            imgs.push((mult, v, t.elapsed().as_secs_f64()));
        }
        let r = imgs[3].1.clone();
        let nr = r.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        let ne = r.iter().map(|z| z.norm().powi(2)).sum::<f64>().sqrt();
        for (mult, v, t) in &imgs[..3] {
            let ec = v.iter().zip(&r).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>().sqrt() / nr;
            let ee = v.iter().zip(&r).map(|(a, b)| (a.norm() - b.norm()).powi(2)).sum::<f64>().sqrt() / ne;
            println!("{name} | PSM m={mult} | {ec:.3e} | {ee:.3e} | {t:.2}");
        }
        for (label, o) in [("MULOK accurate (cubic, up=4, demod)", MulokOptions::accurate()),
                           ("MULOK accurate, up=8", MulokOptions { upsample: 8, ..MulokOptions::accurate() }),
                           ("MULOK thesis, linear up=4 no demod", MulokOptions::thesis())] {
            let t = std::time::Instant::now();
            let img = mulok(&d, &a, &m, &o).unwrap();
            let t = t.elapsed().as_secs_f64();
            let v: Vec<C64> = img.layers.iter().flat_map(|l| l.data.iter().copied()).collect();
            let ec = v.iter().zip(&r).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>().sqrt() / nr;
            let ee = v.iter().zip(&r).map(|(a, b)| (a.norm() - b.norm()).powi(2)).sum::<f64>().sqrt() / ne;
            println!("{name} | {label} | {ec:.3e} | {ee:.3e} | {t:.2}");
        }
    }
}
