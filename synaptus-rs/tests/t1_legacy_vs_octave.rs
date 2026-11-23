//! T1 (TEST_PLAN.md): legacy Rust port vs original Octave psm.m / mulok.m.
//! Consistency check only — agreement does not establish correctness.

use std::path::Path;
use synaptus::{io, legacy, FocusedImage};

const TOL: f64 = 1e-9; // pre-registered in TEST_PLAN.md (T1)

fn compare(name: &str, alg: &str, img: &FocusedImage) -> Vec<String> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");
    let mut failures = vec![];
    for (l, layer) in img.layers.iter().enumerate() {
        let (z, nz, nx, ny, refd) = io::load_reference_layer(&dir, &format!("{name}_{alg}"), l + 1).unwrap();
        assert_eq!((nz, nx, ny), (layer.nz, layer.nx, layer.ny), "{name} {alg} layer {l}: shape");
        let zerr = z.iter().zip(&layer.z).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let num: f64 = layer.data.iter().zip(&refd).map(|(a, b)| (a - b).norm_sqr()).sum::<f64>().sqrt();
        let den: f64 = refd.iter().map(|b| b.norm_sqr()).sum::<f64>().sqrt();
        let rel = num / den;
        println!("{name:6} {alg:6} layer {} : nz={nz:4} rel.err={rel:.3e}  max|dz|={zerr:.1e}", l + 1);
        if !(rel <= TOL) || !(zerr <= 1e-12) {
            failures.push(format!("{name} {alg} layer {}: rel={rel:.3e}, zerr={zerr:.1e}", l + 1));
        }
    }
    failures
}

#[test]
fn t1_legacy_port_matches_octave() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");
    let mut failures = vec![];
    for name in ["wire", "pins", "plane"] {
        let (data, acq, medium) = io::load_dataset(&dir, name).unwrap();
        let p = legacy::psm(&data, &acq, &medium).unwrap();
        failures.extend(compare(name, "psm", &p));
        let m = legacy::mulok(&data, &acq, &medium, &Default::default()).unwrap();
        failures.extend(compare(name, "mulok", &m));
    }
    assert!(failures.is_empty(), "T1 failures:\n{}", failures.join("\n"));
}
