use crate::error::{invalid, Result};
use num_complex::Complex64;

/// Raw pulse-echo data p(t, x, y).
///
/// Storage is **column-major with t fastest** (identical to MATLAB's `ptxy(:)`):
/// `samples[t + nt * (x + nx * y)]`. 2D data has `ny == 1`.
#[derive(Debug, Clone)]
pub struct PulseEchoData {
    pub nt: usize,
    pub nx: usize,
    pub ny: usize,
    pub samples: Vec<f64>,
}

impl PulseEchoData {
    pub fn new(nt: usize, nx: usize, ny: usize, samples: Vec<f64>) -> Result<Self> {
        if nt == 0 || nx == 0 || ny == 0 {
            return invalid(format!("empty data: nt={nt}, nx={nx}, ny={ny}"));
        }
        let n = nt
            .checked_mul(nx)
            .and_then(|v| v.checked_mul(ny))
            .ok_or_else(|| crate::Error::SizeOverflow("nt*nx*ny".into()))?;
        if samples.len() != n {
            return invalid(format!("samples.len()={} but nt*nx*ny={n}", samples.len()));
        }
        if let Some(i) = samples.iter().position(|v| !v.is_finite()) {
            return invalid(format!("non-finite sample at linear index {i}"));
        }
        Ok(Self { nt, nx, ny, samples })
    }

    pub fn is_3d(&self) -> bool {
        self.ny > 1
    }

    /// Column (all t) for spatial position (x, y).
    pub fn column(&self, x: usize, y: usize) -> &[f64] {
        let start = self.nt * (x + self.nx * y);
        &self.samples[start..start + self.nt]
    }
}

/// Acquisition parameters (SI units).
#[derive(Debug, Clone, Copy)]
pub struct Acquisition {
    /// Sampling frequency [Hz].
    pub fs: f64,
    /// Delay between pulse transmission and first sample [s].
    pub t_delay: f64,
    /// Lower edge of the processed band [Hz].
    pub f_low: f64,
    /// Upper edge of the processed band [Hz].
    pub f_high: f64,
    /// Spatial sampling interval along x [m].
    pub x_step: f64,
    /// Spatial sampling interval along y [m] (ignored for 2D data).
    pub y_step: f64,
}

/// Layered medium. Layer `l` has velocity `velocities[l]` and thickness
/// `thicknesses[l]`; the last thickness sets the end of the region of interest.
#[derive(Debug, Clone)]
pub struct Medium {
    pub velocities: Vec<f64>,
    pub thicknesses: Vec<f64>,
}

impl Medium {
    pub fn single(velocity: f64, thickness: f64) -> Self {
        Self { velocities: vec![velocity], thicknesses: vec![thickness] }
    }
    pub fn n_layers(&self) -> usize {
        self.velocities.len()
    }
    /// z coordinate of the top of each layer (cumulative thickness), length L.
    pub fn interfaces(&self) -> Vec<f64> {
        let mut z = Vec::with_capacity(self.thicknesses.len());
        let mut acc = 0.0;
        for &d in &self.thicknesses {
            z.push(acc);
            acc += d;
        }
        z
    }
}

/// Focused image of one layer: complex values, `data[(iz * nx + ix) * ny + iy]`.
#[derive(Debug, Clone)]
pub struct LayerImage {
    pub nz: usize,
    pub nx: usize,
    pub ny: usize,
    /// Depth of each image plane [m].
    pub z: Vec<f64>,
    pub data: Vec<Complex64>,
}

impl LayerImage {
    pub fn at(&self, iz: usize, ix: usize, iy: usize) -> Complex64 {
        self.data[(iz * self.nx + ix) * self.ny + iy]
    }
}

/// Focused image: one [`LayerImage`] per layer, plus lateral pixel coordinates.
#[derive(Debug, Clone)]
pub struct FocusedImage {
    pub layers: Vec<LayerImage>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Validate the inputs shared by all algorithms. Returns the z offset c₁·t_delay/2.
pub(crate) fn validate(data: &PulseEchoData, acq: &Acquisition, medium: &Medium) -> Result<f64> {
    let finite_pos = |v: f64| v.is_finite() && v > 0.0;
    if !finite_pos(acq.fs) {
        return invalid(format!("fs must be finite and > 0 (got {})", acq.fs));
    }
    if !(acq.t_delay.is_finite() && acq.t_delay >= 0.0) {
        return invalid(format!("t_delay must be finite and >= 0 (got {})", acq.t_delay));
    }
    if !(acq.f_low.is_finite() && acq.f_low > 0.0) {
        return invalid(format!("f_low must be finite and > 0 (got {})", acq.f_low));
    }
    if !(acq.f_high.is_finite() && acq.f_high > acq.f_low) {
        return invalid(format!("need f_high > f_low (got {} <= {})", acq.f_high, acq.f_low));
    }
    if acq.f_high > acq.fs / 2.0 {
        return invalid(format!("f_high ({}) exceeds Nyquist fs/2 ({})", acq.f_high, acq.fs / 2.0));
    }
    if !finite_pos(acq.x_step) {
        return invalid(format!("x_step must be finite and > 0 (got {})", acq.x_step));
    }
    if data.is_3d() && !finite_pos(acq.y_step) {
        return invalid(format!("y_step must be finite and > 0 for 3D data (got {})", acq.y_step));
    }
    let l = medium.velocities.len();
    if l == 0 || medium.thicknesses.len() != l {
        return invalid(format!(
            "need equal, non-zero numbers of velocities ({l}) and thicknesses ({})",
            medium.thicknesses.len()
        ));
    }
    if let Some(v) = medium.velocities.iter().find(|v| !finite_pos(**v)) {
        return invalid(format!("velocities must be finite and > 0 (got {v})"));
    }
    if let Some(d) = medium.thicknesses.iter().find(|d| !finite_pos(**d)) {
        return invalid(format!("thicknesses must be finite and > 0 (got {d})"));
    }
    let z_offset = acq.t_delay * (medium.velocities[0] / 2.0);
    if z_offset >= medium.thicknesses[0] {
        return invalid(format!(
            "measurement window starts at z={z_offset} m, at or below the first layer's bottom ({} m)",
            medium.thicknesses[0]
        ));
    }
    Ok(z_offset)
}
