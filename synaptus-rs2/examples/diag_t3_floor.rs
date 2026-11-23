//! Diagnostic (not a pre-registered test): is the T3 layer-1 floor in N5 a NUFFT
//! tolerance effect or a quadrature effect? Vary tol at fixed kappa set.
#[path = "../tests/common/synth.rs"]
#[allow(dead_code)]
mod synth;
use synaptus::psm_sq::{psm_sq, PsmSqOptions};
use synaptus::{Acquisition, Medium, PulseEchoData};
use synth::{Pulse, Scene};
fn main() {
    let s = Scene { fs: 20e6, nt: 1000, nx: 256, x_step: 0.25e-3, t_delay: 20e-6, pulse: Pulse { fc: 2e6, sigma_t: 0.6e-6 }, s0: 0.2 };
    let a = Acquisition { fs: s.fs, t_delay: s.t_delay, f_low: 0.5e6, f_high: 3.5e6, x_step: s.x_step, y_step: s.x_step };
    let d = PulseEchoData::new(s.nt, s.nx, 1, s.two_layer(1480.0, 2730.0, 30e-3, 32.1e-3, 45.3e-3)).unwrap();
    let m = Medium { velocities: vec![1480.0, 2730.0], thicknesses: vec![30e-3, 25e-3] };
    let data_norm: f64 = d.samples.iter().map(|v| v * v).sum::<f64>().sqrt();
    for tol in [1e-12, 1e-14] {
        let r = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: 6.0, tol, ..Default::default() }).unwrap();
        let rn: f64 = r.layers[0].data.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        let r2: f64 = r.layers[1].data.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        println!("tol={tol:.0e}: ||layer1|| = {rn:.3e}, ||layer2|| = {r2:.3e} (ratio {:.1e}); ||data|| = {data_norm:.3e}", rn / r2);
        for k in [2.0, 3.0, 4.0] {
            let i = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: k, tol, ..Default::default() }).unwrap();
            let num: f64 = i.layers[0].data.iter().zip(&r.layers[0].data).map(|(x, y)| (x - y).norm_sqr()).sum();
            println!("  kappa={k}: layer-1 rel L2 to kappa=6 (same tol) = {:.3e}", (num.sqrt()) / rn);
        }
    }
}
