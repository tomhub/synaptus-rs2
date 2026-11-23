//! Compute every array needed to re-create the Synaptus figures (and the
//! comparison/accuracy plots) with the Rust implementation. Plot scripts in
//! `figures/` only render these arrays (log scaling, thresholds, colour maps).
//!
//! Output: figures/data/<name>.f64 (raw f64 LE, column-major / first index
//! fastest) + <name>.txt ("dims ..."). Run from the crate root:
//!   cargo run --release --example figures

use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;
use synaptus::cpsm::{CylAcquisition, TransferFunction};
use synaptus::legacy::cpsm;
use synaptus::envelope::envelope;
use synaptus::fftutil::C64;
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::psm::{debug_spectrum, psm, DepthStep, PsmOptions};
use synaptus::{io, legacy, FocusedImage, Medium, PulseEchoData};

struct Out(PathBuf);
impl Out {
    fn save(&self, name: &str, dims: &[usize], v: &[f64]) {
        assert_eq!(dims.iter().product::<usize>(), v.len(), "{name}");
        let mut f = std::fs::File::create(self.0.join(format!("{name}.f64"))).unwrap();
        let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
        f.write_all(&bytes).unwrap();
        let d: Vec<String> = dims.iter().map(|x| x.to_string()).collect();
        std::fs::write(self.0.join(format!("{name}.txt")), format!("dims {}\n", d.join(" "))).unwrap();
    }
    /// |layer| as (z, x, y) column-major, plus its z axis.
    fn save_layer(&self, name: &str, img: &FocusedImage, l: usize) {
        let li = &img.layers[l];
        let mut v = vec![0.0; li.data.len()];
        for iz in 0..li.nz { for ix in 0..li.nx { for iy in 0..li.ny {
            v[iz + li.nz * (ix + li.nx * iy)] = li.at(iz, ix, iy).norm();
        }}}
        self.save(name, &[li.nz, li.nx, li.ny], &v);
        self.save(&format!("{name}_z"), &[li.nz], &li.z);
    }
}

fn timed<T>(label: &str, log: &mut Vec<String>, f: impl FnOnce() -> T) -> T {
    let t = Instant::now();
    let r = f();
    let s = t.elapsed().as_secs_f64();
    println!("{label}: {s:.3} s");
    log.push(format!("{label}\t{s:.4}"));
    r
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let out = Out(root.join("figures/data"));
    std::fs::create_dir_all(&out.0).unwrap();
    let ddir = root.join("data");
    let mut log = vec![];

    // ---------- A. Python quick-start (WireTargets) ----------
    let (d, a, m) = io::load_dataset(&ddir, "wire").unwrap();
    let env = timed("wire envelope", &mut log, || envelope(&d));
    out.save("wire_raw_env", &[d.nt, d.nx], &env);
    out.save("wire_t", &[d.nt], &(0..d.nt).map(|i| a.t_delay + i as f64 / a.fs).collect::<Vec<_>>());
    out.save("wire_x", &[d.nx], &(0..d.nx).map(|i| i as f64 * a.x_step).collect::<Vec<_>>());
    // Python PhaseShiftMigration uses dz = (c/2)/(f_high - f_low); psm.m uses (c/2)/f_high.
    let bw = timed("wire PSM (dz bandwidth, as Python)", &mut log, || psm(&d, &a, &m, &PsmOptions { depth_step: DepthStep::Bandwidth, ..PsmOptions::default() }).unwrap());
    out.save_layer("wire_psm_bw", &bw, 0);
    let hf = timed("wire PSM (dz high-freq, as psm.m)", &mut log, || psm(&d, &a, &m, &PsmOptions::default()).unwrap());
    out.save_layer("wire_psm_hf", &hf, 0);
    let hf8 = timed("wire PSM t_fft_mult=8", &mut log, || psm(&d, &a, &m, &PsmOptions { t_fft_mult: 8, ..PsmOptions::default() }).unwrap());
    out.save_layer("wire_psm_hf8", &hf8, 0);
    // (ω, kx) spectrum magnitude in the band (fun plot)
    let (om, kx, _ky, w) = debug_spectrum(&d, &a, &m, &PsmOptions::default()).unwrap();
    let (nw, nk) = (om.len(), kx.len());
    let mut sp = vec![0.0; nw * nk];
    for ik in 0..nk { for iw in 0..nw { sp[iw + nw * ik] = w[ik * nw + iw].norm(); } }
    out.save("wire_spec", &[nw, nk], &sp);
    out.save("wire_spec_f", &[nw], &om.iter().map(|w| w / (2.0 * std::f64::consts::PI)).collect::<Vec<_>>());
    out.save("wire_spec_kx", &[nk], &kx);

    // ---------- B. Thesis Fig. 5.4 / 5.5 (PlaneScan3D) ----------
    let (d, a, m) = io::load_dataset(&ddir, "plane").unwrap();
    let env = timed("plane envelope", &mut log, || envelope(&d));
    out.save("plane_raw_env", &[d.nt, d.nx, d.ny], &env);
    out.save("plane_t", &[d.nt], &(0..d.nt).map(|i| a.t_delay + i as f64 / a.fs).collect::<Vec<_>>());
    let p3 = timed("plane 3D PSM", &mut log, || psm(&d, &a, &m, &PsmOptions::default()).unwrap());
    for l in 0..3 { out.save_layer(&format!("plane_psm3d_L{}", l + 1), &p3, l); }
    // 2D PSM applied to every x–z slice (Fig. 5.5a)
    let p2: Vec<FocusedImage> = timed("plane 2D PSM per y-slice", &mut log, || (0..d.ny).map(|iy| {
        let mut s = Vec::with_capacity(d.nt * d.nx);
        for ix in 0..d.nx { s.extend_from_slice(d.column(ix, iy)); }
        psm(&PulseEchoData::new(d.nt, d.nx, 1, s).unwrap(), &a, &m, &PsmOptions::default()).unwrap()
    }).collect());
    for l in 0..3 {
        let (nz, nx, ny) = (p2[0].layers[l].nz, d.nx, d.ny);
        let mut v = vec![0.0; nz * nx * ny];
        for (iy, img) in p2.iter().enumerate() { for iz in 0..nz { for ix in 0..nx {
            v[iz + nz * (ix + nx * iy)] = img.layers[l].at(iz, ix, 0).norm();
        }}}
        out.save(&format!("plane_psm2d_L{}", l + 1), &[nz, nx, ny], &v);
    }
    out.save("plane_interfaces", &[3], &m.interfaces().iter().zip(&m.thicknesses).map(|(z, _)| *z).chain([m.thicknesses.iter().sum()]).skip(1).collect::<Vec<_>>());

    // ---------- C. Thesis Fig. 7.13 + example_focused_image (ObjectsInPipe, CPSM) ----------
    let p = io::read_params(&ddir.join("pipe.txt")).unwrap();
    let dm = &p["dims"];
    let d = PulseEchoData::new(dm[0] as usize, dm[1] as usize, dm[2] as usize, io::read_f64(&ddir.join("pipe.f64")).unwrap()).unwrap();
    let g = |k: &str| p[k][0];
    let ca = CylAcquisition { fs: g("fs"), t_delay: g("tDelay"), f_low: g("fLow"), f_high: g("fHigh"),
        phi_step: g("phiStep"), z_step: g("zStep"), r0: g("r0"), c: g("cc") };
    let rr = p["rr"].clone();
    let img = timed("pipe CPSM (haun)", &mut log, || cpsm(&d, &ca, &rr, TransferFunction::Haun, 1).unwrap());
    let env = timed("pipe envelope", &mut log, || envelope(&d));
    let (np, nz) = (d.nx, d.ny);
    // raw: max envelope over t; range r = r0 + (tDelay + (1:nT)/fs)·c/2 (cpsm_test3D_ObjectsInPipe.m)
    let r_raw: Vec<f64> = (1..=d.nt).map(|i| ca.r0 + (ca.t_delay + i as f64 / ca.fs) * (ca.c / 2.0)).collect();
    let (mut a_raw, mut r_at_raw, mut a_im, mut r_at_im) = (vec![0.0; np * nz], vec![0.0; np * nz], vec![0.0; np * nz], vec![0.0; np * nz]);
    for iz in 0..nz { for ip in 0..np {
        let col = &env[d.nt * (ip + np * iz)..d.nt * (ip + np * iz + 1)];
        let (k, v) = col.iter().enumerate().fold((0, -1.0), |b, (k, &v)| if v > b.1 { (k, v) } else { b });
        a_raw[ip + np * iz] = v; r_at_raw[ip + np * iz] = r_raw[k];
        let (k, v) = (0..rr.len()).map(|ir| (ir, img.at(ir, ip, iz).norm())).fold((0, -1.0), |b, x| if x.1 > b.1 { x } else { b });
        a_im[ip + np * iz] = v; r_at_im[ip + np * iz] = rr[k];
    }}
    out.save("pipe_cscan_raw", &[np, nz], &a_raw);
    out.save("pipe_topo_raw", &[np, nz], &r_at_raw);
    out.save("pipe_cscan_im", &[np, nz], &a_im);
    out.save("pipe_topo_im", &[np, nz], &r_at_im);
    out.save("pipe_phi", &[np], &img.phi);
    out.save("pipe_z", &[nz], &img.z);

    // ---------- D. Accuracy comparison: legacy vs thesis/accurate (PinsPlexiAluSDH) ----------
    let (d, a, m): (PulseEchoData, _, Medium) = io::load_dataset(&ddir, "pins").unwrap();
    let ref8 = timed("pins PSM t_fft_mult=8 (reference)", &mut log, || psm(&d, &a, &m, &PsmOptions { t_fft_mult: 8, ..PsmOptions::default() }).unwrap());
    let ps = timed("pins PSM default", &mut log, || psm(&d, &a, &m, &PsmOptions::default()).unwrap());
    let lm = timed("pins legacy MULOK", &mut log, || legacy::mulok(&d, &a, &m, &Default::default()).unwrap());
    let am = timed("pins accurate MULOK", &mut log, || mulok(&d, &a, &m, &MulokOptions::accurate()).unwrap());
    let lp = PsmOptions { depth_step: DepthStep::Bandwidth, t_fft_mult: 8, ..PsmOptions::default() };
    let ref8bw = psm(&d, &a, &m, &lp).unwrap();
    for l in 0..3 {
        out.save_layer(&format!("pins_ref8_L{}", l + 1), &ref8, l);
        out.save_layer(&format!("pins_psm_L{}", l + 1), &ps, l);
        out.save_layer(&format!("pins_legmulok_L{}", l + 1), &lm, l);
        out.save_layer(&format!("pins_accmulok_L{}", l + 1), &am, l);
        out.save_layer(&format!("pins_ref8bw_L{}", l + 1), &ref8bw, l);
    }
    let _ = C64::new(0.0, 0.0);
    std::fs::write(out.0.join("timings_rust.tsv"), log.join("\n") + "\n").unwrap();
    println!("done -> {}", out.0.display());
}
