//! Literal ports of the original MATLAB/Octave implementations.
//!
//! These reproduce the original behaviour **including** the deviations from the
//! thesis documented in `research/synaptus-fourier-saft.md` (claims C4–C6):
//! squared Stolt amplitude factor, kz grid starting at 2πf_low/ĉ, and the
//! −(1−1/up)·Δω frequency-axis offset after band upsampling. Use them to
//! reproduce original results; use [`crate::mulok`] for thesis-faithful MULOK.
//! [`cpsm`] is the literal `cpsm.m` (defects C8/C10, evanescent clamping, DC bin);
//! use [`crate::cpsm::cpsm`] for the corrected version.

mod cpsm;
mod mulok;
pub use cpsm::cpsm;
pub use mulok::{mulok, LegacyMulokOptions};

use crate::psm::PsmOptions;
use crate::{Acquisition, FocusedImage, Medium, PulseEchoData, Result};

/// `psm.m` (2D/3D, multilayer) with its original arithmetic.
pub fn psm(data: &PulseEchoData, acq: &Acquisition, medium: &Medium) -> Result<FocusedImage> {
    crate::psm::psm(data, acq, medium, &PsmOptions::legacy())
}
