# Test plan (written 2026-11-23, before any Rust code)

These criteria are fixed in advance. Any later change is appended as a dated amendment, and the original result is kept.

## T1. Port consistency: legacy Rust vs original Octave (`data/*_psm_*`, `data/*_mulok_*`)
- Metric: per layer, ‖rust − octave‖₂ / ‖octave‖₂ over all complex pixels. Pixel grids and z vectors must match exactly in count, with z equal within 1e-12 m.
- Tolerance: **1e-9**.
  Justification: both run the same algorithm in IEEE double. The expected discrepancy comes from the FFT library (FFTW vs rustfft, about log₂N·u ≈ 2.5e-15 per transform), summation order, and recursion amplification (J ≤ ~1000 steps, giving ~J·u ≈ 2e-13). The normwise metric is dominated by bright pixels. 1e-9 leaves more than three orders of margin above the predicted ~1e-12 while still detecting any algorithmic difference (an algorithmic slip gives errors ≫ 1e-6).
- This is a **consistency** check only. Agreement does not mean the algorithm is correct.

## T2. Synthetic point scatterer, single layer (analytic oracle)
- Data: p(t,x) = g(t − 2r(x)/c), with r = √((x−x₀)² + z₀²) and g a Gaussian-modulated cosine (f_c = 2 MHz, σ_t = 0.4 µs), sampled at fs = 20 MHz. The band [0.5, 3.5] MHz contains the pulse spectrum down to < 1e-6 of peak.
- Expected: the maximum of |image| lies within **1 pixel** of (x₀, z₀) in x and within **1 depth sample** in z. These follow from the sampling grid: the true location is generally off-grid, so the nearest sample is within half a pixel, and one pixel allows for peak-shape asymmetry.
- Applies to: legacy PSM, improved PSM, legacy MULOK, thesis MULOK.

## T3. Synthetic point scatterer, two layers (independent ray-traced oracle)
- Travel time from Fermat's principle (the minimum over the interface crossing point, by golden-section search), for a scatterer in layer 2. Same pulse.
- Expected: the peak within 1 pixel in x and 1 depth sample in z of the true position. Applies to PSM and thesis MULOK.

## T4. Thesis-MULOK ≡ PSM (C3 in the dossier: an exact change of variables)
- Metric: ‖MULOK − PSM‖ / ‖PSM‖ on the same depth/x grid (both with dz = ĉ/f_high), for the synthetic data.
- Expected: decreases as interpolation improves (linear → cubic → windowed sinc) and with time zero-padding (up = 1, 2, 4, 8). No fixed numeric threshold is declared in advance. The *claim* tested is monotone convergence toward a floor set by the band-edge discretisation. Legacy MULOK is evaluated with the same metric after a least-squares complex scale fit (its normalisation is arbitrary by design). This scale fit is declared here, in advance.

## T5. PSM rounding (improvement claim)
- Reference: high-precision (mpmath, 40 digits) evaluation of eq. (3.21) from the same double-precision inputs (ω, kx, c, dz, wavefield).
- Metric: max relative error per depth line vs depth index j, for recursive, direct and re-anchored updates, and naive vs Neumaier summation.
- Claim to test: recursion error grows roughly linearly with j, while direct and re-anchored stay at O(u). No threshold. Report the measured curves.

## T6. Input validation
- Invalid inputs (f_low ≥ f_high, f_high > fs/2, non-positive steps/velocities, mismatched layer vectors, non-finite samples, empty data, z offset beyond the first layer) must return `Err`. They must not panic or produce output.

## Amendment A1 (2026-11-23, before any synthetic test was run)
T2 is internally inconsistent: with σ_t = 0.4 µs the Gaussian spectrum at the band edges (±1.5 MHz from f_c) is exp(−2π²σ²Δf²) = e^{−7.1} ≈ 8e-4, not < 1e-6 as required. **Change:** σ_t = 0.6 µs, giving e^{−16} ≈ 1.1e-7 < 1e-6. To keep the synthetic data free of spatial aliasing (not addressed originally), the synthetic scatterer also gets a frequency-independent directivity weight exp(−(sinθ/0.2)²), and the spatial step is set to 0.25 mm. No result had been seen when this amendment was made.
T4 addition: to compare legacy MULOK (depth step ĉ/(f_h−f_l)) pixel-for-pixel, PSM gets an option for the same depth step. Its PSM reference uses that grid.

## Amendment A2 (2026-11-23, AFTER seeing T4 run `runs/20261123-145346-accuracy-mulok`)
Observation: legacy MULOK gave complex-metric errors of 0.87–1.0 vs PSM even after a complex scale fit, although it focused correctly in T2.
Cause (DERIVED, to be checked): legacy processes the negative-ω band, which conjugates the image, and inverse-transforms a kz grid starting at k_low, which multiplies the image by e^{−i k_low z}. Both are phase conventions that leave |image| unchanged. The pre-declared complex metric (A1/T4) therefore cannot separate image-quality defects from phase conventions.
Additions (the original results are kept unchanged):
- **E-metric**: ‖ |M| − α|P| ‖ / ‖ |P| ‖ on envelopes, with α ≥ 0 the least-squares real scale (α = 1 reported separately as "raw").
- **Attribution re-run** with an accurate interpolation baseline (cubic, up = 4, demodulated), so each legacy defect is measured without linear-interpolation error masking it.
- **PSM quadrature self-check**: PSM with t_fft_mult = 1 vs PSM with t_fft_mult = 2 (×½). This measures the discretisation error of the PSM reference itself, i.e. how far "agreement with PSM" can be trusted.

## T7 (added 2026-11-23, before the CPSM port was written): legacy CPSM vs Octave `cpsm.m`
- Dataset `CylScan3D_ObjectsInPipe` cropped as in `cpsm_test3D_ObjectsInPipe.m` (r ≤ 0.147 m), 'haun' transfer function. Also 'gardner' on the same data.
- Metric and tolerance as T1: relative L2 error ≤ **1e-9**. Same justification: same algorithm in IEEE double, with FFT/summation rounding of order 1e-13. Range vectors are passed identically.

## T8 (added same time): analytic-signal envelope
- x(t) = g(t)·cos(2πf₀t), with g Gaussian (σ = 1 µs, f₀ = 3 MHz, fs = 25 MHz, N = 1024, centred). Expected |H{x}| = g to within the spectral overlap with DC/Nyquist, which is below e^{−(2πσf₀)²/2} ≈ e^{−177}, i.e. negligible. Tolerance **1e-12** absolute (rounding of a length-1024 FFT pair, ~log₂N·u·max|x|, with 3 orders of margin).

## Amendment A3 (2026-11-23, before the corrected CPSM was written): CPSM defaults must not reproduce known defects
**Change.** The literal port of `cpsm.m` moves to `legacy::cpsm`; **T7 now tests `legacy::cpsm`**, with criteria unchanged. `cpsm::cpsm` becomes a corrected implementation, used by default:
(i) the angular FFT is zero-padded to the full circle N = round(2π/Δφ), so orders are integers (thesis eq. 6.4; dossier C10);
(ii) evanescent components (k_r² ≤ 0, and the Haun radicand ≤ 0) are **removed**, not clamped to zero phase;
(iii) the DC bin is excluded;
(iv) the default transfer function is the Hankel ratio (Gardner–Collin closed form), whose amplitude agrees with thesis eq. 6.7. Haun stays available. C8 stays open.

## T9 (added with A3, before the code): corrected CPSM, synthetic cylindrical point scatterers (analytic oracle)
- Geometry: transducer at radius r₀ = 58 mm pointing radially outward, Δφ = 2π/720, n_φ = 96, Δz = 0.8 mm, n_z = 64, f_s = 25 MHz, n_t = 1024, c = 1480 m/s, band [0, 4] MHz (exercises the DC exclusion). Pulse: Gaussian-modulated cosine, f_c = 2 MHz, σ_t = 0.6 µs. Directivity exp(−(sin θ/0.2)²) about the outward normal.
- Data: p = Σ_s w·g(t − 2d_s/c), with d_s the straight-line distance from (r₀, φ_i, z_j) to scatterer s. Two scatterers: (r, φ, z) = (75 mm, +3.0°, 24.0 mm) and (82 mm, −5.0°, 30.4 mm). φ convention as `CylImage` (centred).
- Expected: for `cpsm::cpsm` with default options, the maximum of |image| within ±3 samples (in each of r, φ, z) of each true position lies within **1 sample** of it in r, φ and z.
- Reported, not pass/fail: the same location error for `legacy::cpsm` ('haun', its default in the original), and a focus measure: the fraction of the image energy within ±2 samples of each true position.
