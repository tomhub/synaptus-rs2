//! Arrays for the v2 gallery (figures/FIGURES.md): what the original code computes
//! (via the `legacy` ports, equal to Octave to ~1e-15) next to what v2 computes by
//! default. Only computes; `figures/plot_gallery.py` renders.
//!   cargo run --release --example gallery_v2      # -> figures/data/gallery/
use std::io::Write;
use std::path::Path;
use synaptus::cpsm::{cpsm, CpsmOptions, CylAcquisition, TransferFunction};
use synaptus::mulok::{mulok, MulokOptions};
use synaptus::psm::{psm, PsmOptions};
use synaptus::psm_sq::{psm_sq, PsmSqOptions};
use synaptus::{io, legacy, FocusedImage, PulseEchoData};

fn save(dir: &Path, name: &str, dims: &[usize], v: &[f64]) {
    let mut f = std::fs::File::create(dir.join(format!("{name}.f64"))).unwrap();
    f.write_all(&v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
    std::fs::write(dir.join(format!("{name}.txt")), format!("dims {}\n", dims.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(" "))).unwrap();
}

/// Complex image (2D, y = 0), layers concatenated: re/im as (z, x) column-major + z axis.
fn save_img(dir: &Path, name: &str, img: &FocusedImage) {
    let nx = img.layers[0].nx;
    let nz: usize = img.layers.iter().map(|l| l.nz).sum();
    let (mut re, mut im, mut z) = (vec![0.0; nz * nx], vec![0.0; nz * nx], vec![]);
    let mut off = 0;
    for l in &img.layers {
        for iz in 0..l.nz { for ix in 0..nx { let c = l.at(iz, ix, 0); re[off + iz + nz * ix] = c.re; im[off + iz + nz * ix] = c.im; } }
        z.extend_from_slice(&l.z);
        off += l.nz;
    }
    save(dir, &format!("{name}_re"), &[nz, nx], &re);
    save(dir, &format!("{name}_im"), &[nz, nx], &im);
    save(dir, &format!("{name}_z"), &[nz], &z);
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let out = root.join("figures/data/gallery");
    std::fs::create_dir_all(&out).unwrap();
    let ddir = root.join("data");

    // Wire targets: original PSM (psm.m arithmetic) vs SQ-PSM (the limit)
    let (d, a, m) = io::load_dataset(&ddir, "wire").unwrap();
    save_img(&out, "wire_orig", &legacy::psm(&d, &a, &m).unwrap());
    save_img(&out, "wire_v1", &psm(&d, &a, &m, &PsmOptions::default()).unwrap());
    save_img(&out, "wire_sq", &psm_sq(&d, &a, &m, &PsmSqOptions::default()).unwrap());

    // Pins: original MULOK vs default MULOK vs SQ-PSM
    let (d, a, m) = io::load_dataset(&ddir, "pins").unwrap();
    save_img(&out, "pins_origmulok", &legacy::mulok(&d, &a, &m, &Default::default()).unwrap());
    save_img(&out, "pins_mulok", &mulok(&d, &a, &m, &MulokOptions::default()).unwrap());
    save_img(&out, "pins_sq", &psm_sq(&d, &a, &m, &PsmSqOptions::default()).unwrap());
    // SQ-PSM on MULOK's original (bandwidth) depth grid, so the original can be compared pixel by pixel
    save_img(&out, "pins_sq_bw", &psm_sq(&d, &a, &m, &PsmSqOptions { depth_step: synaptus::psm::DepthStep::Bandwidth, ..Default::default() }).unwrap());

    // Flat-bottom holes (3D): original psm.m vs SQ-PSM, envelope max over y per layer
    let (d, a, m) = io::load_dataset(&ddir, "plane").unwrap();
    for (name, img) in [("plane_orig", legacy::psm(&d, &a, &m).unwrap()), ("plane_sq", psm_sq(&d, &a, &m, &PsmSqOptions::default()).unwrap())] {
        let nx = img.layers[0].nx;
        let nz: usize = img.layers.iter().map(|l| l.nz).sum();
        let (mut v, mut z) = (vec![0.0; nz * nx], vec![]);
        let mut off = 0;
        for l in &img.layers {
            for iz in 0..l.nz { for ix in 0..nx {
                v[off + iz + nz * ix] = (0..l.ny).map(|iy| l.at(iz, ix, iy).norm()).fold(0.0, f64::max);
            }}
            z.extend_from_slice(&l.z);
            off += l.nz;
        }
        save(&out, &format!("{name}_envmaxy"), &[nz, nx], &v);
        save(&out, &format!("{name}_z"), &[nz], &z);
    }

    // Objects in pipe: original CPSM (cpsm.m, 'haun') vs corrected CPSM default
    let p = io::read_params(&ddir.join("pipe.txt")).unwrap();
    let dm = &p["dims"];
    let d = PulseEchoData::new(dm[0] as usize, dm[1] as usize, dm[2] as usize, io::read_f64(&ddir.join("pipe.f64")).unwrap()).unwrap();
    let g = |k: &str| p[k][0];
    let ca = CylAcquisition { fs: g("fs"), t_delay: g("tDelay"), f_low: g("fLow"), f_high: g("fHigh"), phi_step: g("phiStep"), z_step: g("zStep"), r0: g("r0"), c: g("cc") };
    let rr = p["rr"].clone();
    for (name, img) in [
        ("pipe_orig", legacy::cpsm(&d, &ca, &rr, TransferFunction::Haun, 1).unwrap()),
        ("pipe_fixed", cpsm(&d, &ca, &rr, &CpsmOptions::default()).unwrap()),
        ("pipe_fixed_haun", cpsm(&d, &ca, &rr, &CpsmOptions { transfer: TransferFunction::Haun, ..Default::default() }).unwrap()),
    ] {
        let (np, nz, nr) = (img.nphi, img.nz, img.r.len());
        let (mut amp, mut rng) = (vec![0.0; np * nz], vec![0.0; np * nz]);
        for iz in 0..nz { for ip in 0..np {
            let (k, v) = (0..nr).map(|ir| (ir, img.at(ir, ip, iz).norm())).fold((0, -1.0), |b, x| if x.1 > b.1 { x } else { b });
            amp[ip + np * iz] = v;
            rng[ip + np * iz] = img.r[k];
        }}
        save(&out, &format!("{name}_cscan"), &[np, nz], &amp);
        save(&out, &format!("{name}_topo"), &[np, nz], &rng);
        if name == "pipe_orig" {
            save(&out, "pipe_phi", &[np], &img.phi);
            save(&out, "pipe_z", &[nz], &img.z);
        }
    }
    println!("done -> {}", out.display());
}
