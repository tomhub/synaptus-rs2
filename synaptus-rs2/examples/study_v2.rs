//! N7 (pass/fail) and S1 (study) from TEST_PLAN_V2.md, plus arrays for the v2 comparison page.
//!   cargo run --release --example study_v2 -- <out_dir>
//! Writes <out_dir>/n7.tsv, s1.tsv and figure arrays under <out_dir>/fig/.
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::psm::{psm, DepthStep, PsmOptions};
use synaptus::psm_nufft::{psm_nufft, PsmNufftOptions};
use synaptus::psm_sq::{psm_sq, PsmSqOptions, DEFAULT_KAPPA};
use synaptus::{io, legacy, FocusedImage};

/// v1 ACCURACY E-metric on envelopes, all layers concatenated: ‖|M| − α|P|‖/‖α|P|‖.
fn e_metric(m: &FocusedImage, p: &FocusedImage) -> f64 {
    let (mut mp, mut pp) = (0.0, 0.0);
    let pairs: Vec<(f64, f64)> = m.layers.iter().zip(&p.layers).flat_map(|(a, b)| {
        assert_eq!(a.z.len(), b.z.len());
        a.data.iter().zip(&b.data).map(|(x, y)| (x.norm(), y.norm())).collect::<Vec<_>>()
    }).collect();
    for (a, b) in &pairs { mp += a * b; pp += b * b; }
    let alpha = mp / pp;
    let num: f64 = pairs.iter().map(|(a, b)| (a - alpha * b).powi(2)).sum();
    num.sqrt() / (alpha * pp.sqrt())
}

fn bench<T>(f: impl Fn() -> T) -> (T, f64, f64, f64) {
    let _ = f(); // warm-up
    let mut ts = vec![];
    let mut last = None;
    for _ in 0..5 {
        let t = Instant::now();
        last = Some(f());
        ts.push(t.elapsed().as_secs_f64());
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (last.unwrap(), ts[2], ts[0], ts[4])
}

fn save(dir: &Path, name: &str, dims: &[usize], v: &[f64]) {
    let mut f = std::fs::File::create(dir.join(format!("{name}.f64"))).unwrap();
    f.write_all(&v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
    std::fs::write(dir.join(format!("{name}.txt")), format!("dims {}\n", dims.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(" "))).unwrap();
}

/// Envelope (layers concatenated along z, (z, x) column-major for y = 0) + z axis.
fn save_env(dir: &Path, name: &str, img: &FocusedImage) {
    let nx = img.layers[0].nx;
    let mut v = vec![];
    let mut z = vec![];
    let nzt: usize = img.layers.iter().map(|l| l.nz).sum();
    v.resize(nzt * nx, 0.0);
    let mut off = 0;
    for l in &img.layers {
        for iz in 0..l.nz { for ix in 0..nx { v[off + iz + nzt * ix] = l.at(iz, ix, 0).norm(); } }
        z.extend_from_slice(&l.z);
        off += l.nz;
    }
    save(dir, name, &[nzt, nx], &v);
    save(dir, &format!("{name}_z"), &[nzt], &z);
}

fn main() {
    let out = PathBuf::from(std::env::args().nth(1).expect("out dir"));
    let fig = out.join("fig");
    std::fs::create_dir_all(&fig).unwrap();
    let ddir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");
    let mut n7 = vec!["dataset\tm\tE_vs_SQ".to_string()];
    let mut s1 = vec!["dataset\tmethod\tE_vs_ref\tt_median_s\tt_min_s\tt_max_s".to_string()];
    let mut n7_pass = true;
    for name in ["wire", "pins", "plane"] {
        let (d, a, m) = io::load_dataset(&ddir, name).unwrap();
        let sq_opts = PsmSqOptions::default();
        let (sq, ts, ts0, ts1) = bench(|| psm_sq(&d, &a, &m, &sq_opts).unwrap());
        let reference = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: 2.0 * DEFAULT_KAPPA, ..Default::default() }).unwrap();
        let ref_bw = psm_sq(&d, &a, &m, &PsmSqOptions { kappa: 2.0 * DEFAULT_KAPPA, depth_step: DepthStep::Bandwidth, ..Default::default() }).unwrap();

        // ---- N7: PSM_m (via A, tied to v1 by N4) converges to SQ ----
        let ms: Vec<usize> = vec![1, 2, 4, 8, 16];
        let mut e_m = vec![];
        for &mm in &ms {
            let img = psm_nufft(&d, &a, &m, &PsmNufftOptions { t_fft_mult: mm, ..Default::default() }).unwrap();
            let e = e_metric(&img, &sq);
            println!("N7 {name:5} m={mm:2}: E(PSM_m, SQ) = {e:.3e}");
            n7.push(format!("{name}\t{mm}\t{e:.6e}"));
            e_m.push(e);
        }
        let ok = e_m[4] <= e_m[0] / 4.0;
        println!("N7 {name}: E(16)/E(1) = {:.3e} -> {}", e_m[4] / e_m[0], if ok { "PASS" } else { "FAIL" });
        n7_pass &= ok;

        // ---- S1: accuracy and cost vs reference (SQ, kappa = 2 x default) ----
        let mut row = |label: &str, img: &FocusedImage, r: &FocusedImage, t: (f64, f64, f64)| {
            let e = e_metric(img, r);
            println!("S1 {name:5} {label:28} E={e:.3e}  t={:.3}s [{:.3}, {:.3}]", t.0, t.1, t.2);
            s1.push(format!("{name}\t{label}\t{e:.6e}\t{:.4}\t{:.4}\t{:.4}", t.0, t.1, t.2));
        };
        let (p1, t, t0, t1) = bench(|| psm(&d, &a, &m, &PsmOptions::default()).unwrap());
        row("v1 PSM (m=1)", &p1, &reference, (t, t0, t1));
        let (p8, t, t0, t1) = bench(|| psm(&d, &a, &m, &PsmOptions { t_fft_mult: 8, ..PsmOptions::default() }).unwrap());
        row("v1 PSM (m=8)", &p8, &reference, (t, t0, t1));
        let (lm, t, t0, t1) = bench(|| legacy::mulok(&d, &a, &m, &Default::default()).unwrap());
        row("v1 legacy MULOK", &lm, &ref_bw, (t, t0, t1));
        let (am, t, t0, t1) = bench(|| mulok(&d, &a, &m, &MulokOptions::accurate()).unwrap());
        row("v1 accurate MULOK", &am, &reference, (t, t0, t1));
        let mut nufft_ms = vec![1usize, 8, 32];
        if name == "plane" { nufft_ms = vec![1, 8, 16]; } // memory cap (~2 GB at m = 32 in 3D)
        let mut a_imgs = vec![];
        for mm in nufft_ms {
            let (img, t, t0, t1) = bench(|| psm_nufft(&d, &a, &m, &PsmNufftOptions { t_fft_mult: mm, ..Default::default() }).unwrap());
            row(&format!("v2 NUFFT-PSM (m={mm})"), &img, &reference, (t, t0, t1));
            a_imgs.push((mm, img));
        }
        row(&format!("v2 SQ-PSM (kappa={DEFAULT_KAPPA})"), &sq, &reference, (ts, ts0, ts1));

        if name != "plane" {
            save_env(&fig, &format!("{name}_ref"), &reference);
            save_env(&fig, &format!("{name}_refbw"), &ref_bw);
            save_env(&fig, &format!("{name}_v1psm"), &p1);
            save_env(&fig, &format!("{name}_v1psm8"), &p8);
            save_env(&fig, &format!("{name}_legmulok"), &lm);
            save_env(&fig, &format!("{name}_accmulok"), &am);
            save_env(&fig, &format!("{name}_sq"), &sq);
            for (mm, img) in &a_imgs { save_env(&fig, &format!("{name}_nufft_m{mm}"), img); }
        }
    }
    std::fs::write(out.join("n7.tsv"), n7.join("\n") + "\n").unwrap();
    std::fs::write(out.join("s1.tsv"), s1.join("\n") + "\n").unwrap();
    println!("N7 overall: {}", if n7_pass { "PASS" } else { "FAIL" });
}
