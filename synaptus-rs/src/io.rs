//! Reader for datasets exported by `octave/export_and_reference.m`:
//! `<name>.f64` (raw f64 LE, MATLAB column-major) plus `<name>.txt` (key values...).

use crate::error::{Error, Result};
use crate::fftutil::C64;
use crate::{Acquisition, Medium, PulseEchoData};
use std::collections::HashMap;
use std::path::Path;

pub type Params = HashMap<String, Vec<f64>>;

pub fn read_params(path: &Path) -> Result<Params> {
    let s = std::fs::read_to_string(path).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
    let mut m = HashMap::new();
    for line in s.lines().filter(|l| !l.trim().is_empty()) {
        let mut it = line.split_whitespace();
        let key = it.next().unwrap().to_string();
        let vals = it
            .map(|t| t.parse::<f64>().map_err(|e| Error::Io(format!("{}: {key}: {e}", path.display()))))
            .collect::<Result<Vec<_>>>()?;
        m.insert(key, vals);
    }
    Ok(m)
}

pub fn read_f64(path: &Path) -> Result<Vec<f64>> {
    let b = std::fs::read(path).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
    if b.len() % 8 != 0 {
        return Err(Error::Io(format!("{}: length not a multiple of 8", path.display())));
    }
    Ok(b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect())
}

pub fn read_c64(path: &Path) -> Result<Vec<C64>> {
    Ok(read_f64(path)?.chunks_exact(2).map(|c| C64::new(c[0], c[1])).collect())
}

fn one(p: &Params, k: &str) -> Result<f64> {
    p.get(k).and_then(|v| v.first().copied()).ok_or_else(|| Error::Io(format!("missing parameter '{k}'")))
}

/// Load `<dir>/<name>.{f64,txt}` as data, acquisition and medium.
pub fn load_dataset(dir: &Path, name: &str) -> Result<(PulseEchoData, Acquisition, Medium)> {
    let p = read_params(&dir.join(format!("{name}.txt")))?;
    let dims = p.get("dims").ok_or_else(|| Error::Io("missing dims".into()))?;
    let (nt, nx, ny) = (dims[0] as usize, dims[1] as usize, dims[2] as usize);
    let data = PulseEchoData::new(nt, nx, ny, read_f64(&dir.join(format!("{name}.f64")))?)?;
    let acq = Acquisition {
        fs: one(&p, "fs")?,
        t_delay: one(&p, "tDelay")?,
        f_low: one(&p, "fLow")?,
        f_high: one(&p, "fHigh")?,
        x_step: one(&p, "xStep")?,
        y_step: one(&p, "yStep")?,
    };
    let medium = Medium {
        velocities: p.get("cc").cloned().unwrap_or_default(),
        thicknesses: p.get("thick").cloned().unwrap_or_default(),
    };
    Ok((data, acq, medium))
}

/// Reference layer image from Octave: complex values in MATLAB column-major order
/// over dims (nz, nx[, ny]); returned re-ordered to the crate layout [iz][ix][iy].
pub fn load_reference_layer(dir: &Path, prefix: &str, layer: usize) -> Result<(Vec<f64>, usize, usize, usize, Vec<C64>)> {
    let p = read_params(&dir.join(format!("{prefix}_layer{layer}.txt")))?;
    let dims = p.get("dims").ok_or_else(|| Error::Io("missing dims".into()))?;
    let (nz, nx, ny) = (dims[0] as usize, dims[1] as usize, dims.get(2).copied().unwrap_or(1.0) as usize);
    let z = p.get("z").cloned().unwrap_or_default();
    let raw = read_c64(&dir.join(format!("{prefix}_layer{layer}.c64")))?;
    if raw.len() != nz * nx * ny {
        return Err(Error::Io(format!("{prefix} layer {layer}: size mismatch")));
    }
    let mut out = vec![C64::new(0.0, 0.0); raw.len()];
    for iy in 0..ny {
        for ix in 0..nx {
            for iz in 0..nz {
                out[(iz * nx + ix) * ny + iy] = raw[iz + nz * (ix + nx * iy)];
            }
        }
    }
    Ok((z, nz, nx, ny, out))
}
