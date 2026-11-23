//! Rust port of the Synaptus Fourier-domain synthetic-aperture imaging algorithms.
//!
//! * [`legacy`] — literal ports of the original MATLAB/Octave `psm.m` and `mulok.m`
//!   (linear-interpolation mode), including their known deviations from the thesis.
//!   Use these to reproduce published/original results.
//! * [`psm`] — phase shift migration with selectable numerical kernels
//!   (cancellation-free kz, direct/re-anchored phase, compensated summation).
//! * [`mulok`] — multilayer omega-k (Stolt) migration as specified in the thesis
//!   (Skjelvareid 2012, eqs. 4.15–4.17, 4.25–4.29), with switches that re-introduce
//!   each legacy deviation individually for attribution studies.
//!
//! See `research/synaptus-fourier-saft.md` for the derivations and the list of
//! deviations between the original code and the thesis.
//!
//! Conventions (all modules): forward FFT uses e^{-iωt}; the *positive* frequency
//! band ω ∈ [2πf_low, 2πf_high] is processed, kz ≥ 0, and downward extrapolation by
//! Δz is multiplication with e^{+i kz Δz}. ĉ = c/2 (exploding reflector model).
//! Images are complex; take `.norm()` for the envelope.

pub mod cpsm;
pub mod envelope;
pub mod error;
pub mod fftutil;
pub mod grid;
pub mod io;
pub mod kernels;
pub mod legacy;
pub mod mulok;
pub mod psm;
pub mod types;

pub use error::{Error, Result};
pub use types::{Acquisition, FocusedImage, LayerImage, Medium, PulseEchoData};
