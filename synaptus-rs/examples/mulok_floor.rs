//! Diagnose the MULOK floor on real data: z-period (wrap-around) hypothesis.
use synaptus::fftutil::C64;
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::psm::{psm, PsmOptions};
use synaptus::io;
fn main() {
    println!("dataset | z_fft_mult | envelope err vs PSM m=8 | time [s]\n---|---|---|---");
    for name in ["wire", "pins", "plane"] {
        let (d, a, m) = io::load_dataset(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data"), name).unwrap();
        let r: Vec<C64> = psm(&d, &a, &m, &PsmOptions { t_fft_mult: 8, ..PsmOptions::default() }).unwrap()
            .layers.iter().flat_map(|l| l.data.iter().map(|z| z / 8.0)).collect();
        let ne = r.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        for zm in [1usize, 2, 4] {
            let t = std::time::Instant::now();
            let v: Vec<C64> = mulok(&d, &a, &m, &MulokOptions { z_fft_mult: zm, ..MulokOptions::accurate() }).unwrap()
                .layers.iter().flat_map(|l| l.data.iter().copied()).collect();
            let e = v.iter().zip(&r).map(|(a, b)| (a.norm() - b.norm()).powi(2)).sum::<f64>().sqrt() / ne;
            println!("{name} | {zm} | {e:.3e} | {:.2}", t.elapsed().as_secs_f64());
        }
    }
}
