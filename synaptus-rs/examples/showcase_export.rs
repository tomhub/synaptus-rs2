//! Write the Rust legacy PSM/MULOK and CPSM outputs in exactly the layout of the
//! Octave reference files in data/ (complex f64 LE interleaved, column-major
//! (z, x, y) / (r, φ, z)), so original and port can be shown side by side.
//!
//!   cargo run --release --example showcase_export   # -> figures/data/showcase/
use std::io::Write;
use std::path::Path;
use synaptus::cpsm::{CylAcquisition, TransferFunction};
use synaptus::legacy::cpsm;
use synaptus::fftutil::C64;
use synaptus::{io, legacy, PulseEchoData};

fn write_c64(path: &Path, dims: &[usize], v: &[C64]) {
    let mut f = std::fs::File::create(path.with_extension("c64")).unwrap();
    let b: Vec<u8> = v.iter().flat_map(|c| [c.re.to_le_bytes(), c.im.to_le_bytes()]).flatten().collect();
    f.write_all(&b).unwrap();
    let d: Vec<String> = dims.iter().map(|x| x.to_string()).collect();
    std::fs::write(path.with_extension("txt"), format!("dims {}\n", d.join(" "))).unwrap();
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let (ddir, out) = (root.join("data"), root.join("figures/data/showcase"));
    std::fs::create_dir_all(&out).unwrap();
    for name in ["wire", "pins", "plane"] {
        let (d, a, m) = io::load_dataset(&ddir, name).unwrap();
        let imgs = [("psm", legacy::psm(&d, &a, &m).unwrap()),
                    ("mulok", legacy::mulok(&d, &a, &m, &Default::default()).unwrap())];
        for (alg, img) in imgs {
            for (l, li) in img.layers.iter().enumerate() {
                let mut v = vec![C64::new(0.0, 0.0); li.data.len()];
                for iy in 0..li.ny { for ix in 0..li.nx { for iz in 0..li.nz {
                    v[iz + li.nz * (ix + li.nx * iy)] = li.at(iz, ix, iy);
                }}}
                write_c64(&out.join(format!("{name}_{alg}_layer{}", l + 1)), &[li.nz, li.nx, li.ny], &v);
            }
        }
    }
    let p = io::read_params(&ddir.join("pipe.txt")).unwrap();
    let dm = &p["dims"];
    let (nt, np, nz) = (dm[0] as usize, dm[1] as usize, dm[2] as usize);
    let d = PulseEchoData::new(nt, np, nz, io::read_f64(&ddir.join("pipe.f64")).unwrap()).unwrap();
    let g = |k: &str| p[k][0];
    let ca = CylAcquisition { fs: g("fs"), t_delay: g("tDelay"), f_low: g("fLow"), f_high: g("fHigh"),
        phi_step: g("phiStep"), z_step: g("zStep"), r0: g("r0"), c: g("cc") };
    let rr = p["rr"].clone();
    let img = cpsm(&d, &ca, &rr, TransferFunction::Haun, 1).unwrap();
    let nr = rr.len();
    let mut v = vec![C64::new(0.0, 0.0); nr * np * nz];
    for iz in 0..nz { for ip in 0..np { for ir in 0..nr { v[ir + nr * (ip + np * iz)] = img.at(ir, ip, iz); }}}
    write_c64(&out.join("pipe_cpsm_haun"), &[nr, np, nz], &v);
    println!("done -> {}", out.display());
}
