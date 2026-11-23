//! T7: legacy::cpsm (literal port) vs Octave cpsm.m; T8: analytic-signal envelope (TEST_PLAN.md).
use std::path::Path;
use synaptus::cpsm::{CylAcquisition, TransferFunction};
use synaptus::legacy::cpsm;
use synaptus::{envelope::envelope, io, PulseEchoData};

#[test]
fn t7_cpsm_matches_octave() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");
    let p = io::read_params(&dir.join("pipe.txt")).unwrap();
    let d = &p["dims"];
    let (nt, np, nz) = (d[0] as usize, d[1] as usize, d[2] as usize);
    let data = PulseEchoData::new(nt, np, nz, io::read_f64(&dir.join("pipe.f64")).unwrap()).unwrap();
    let g = |k: &str| p[k][0];
    let acq = CylAcquisition { fs: g("fs"), t_delay: g("tDelay"), f_low: g("fLow"), f_high: g("fHigh"),
        phi_step: g("phiStep"), z_step: g("zStep"), r0: g("r0"), c: g("cc") };
    let rr = p["rr"].clone();
    let mut fails = vec![];
    for (name, tf) in [("haun", TransferFunction::Haun), ("gardner", TransferFunction::Gardner)] {
        let img = cpsm(&data, &acq, &rr, tf, 1).unwrap();
        // Octave im(r, phi, z) column-major
        let raw = io::read_c64(&dir.join(format!("pipe_cpsm_{name}.c64"))).unwrap();
        let nr = rr.len();
        assert_eq!(raw.len(), nr * np * nz);
        let (mut num, mut den) = (0.0, 0.0);
        for iz in 0..nz { for ip in 0..np { for ir in 0..nr {
            let o = raw[ir + nr * (ip + np * iz)];
            num += (img.at(ir, ip, iz) - o).norm_sqr();
            den += o.norm_sqr();
        }}}
        let rel = (num / den).sqrt();
        println!("cpsm {name}: rel.err = {rel:.3e} (nr={nr})");
        if !(rel <= 1e-9) { fails.push(format!("{name}: {rel:.3e}")); }
    }
    assert!(fails.is_empty(), "{fails:?}");
}

#[test]
fn t8_envelope_of_gaussian_burst() {
    let (n, fs, f0, sig) = (1024usize, 25e6, 3e6, 1e-6);
    let t0 = n as f64 / 2.0 / fs;
    let g: Vec<f64> = (0..n).map(|i| { let t = i as f64 / fs - t0; (-(t * t) / (2.0 * sig * sig)).exp() }).collect();
    let x: Vec<f64> = (0..n).map(|i| g[i] * (2.0 * std::f64::consts::PI * f0 * (i as f64 / fs - t0)).cos()).collect();
    let e = envelope(&PulseEchoData::new(n, 1, 1, x).unwrap());
    let err = e.iter().zip(&g).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    println!("envelope max abs error = {err:.3e}");
    assert!(err <= 1e-12, "{err:.3e}");
}
