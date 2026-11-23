# Numerical accuracy report: Synaptus Rust port

Date: 2026-11-23. Every number below is `OBSERVED` in the run directory cited, which also holds the environment record and the command used. Evidence labels follow the evidence labels in the root `README.md`. **Independent verification: PENDING.** These results were produced by the implementer.

## Summary

| Question | Answer | Evidence |
|---|---|---|
| Does the legacy Rust port reproduce the original Octave code? | Yes. Rel. error 4.9e-16 to 1.3e-15 on all 14 layer images (tolerance 1e-9, set in advance) | `runs/20261123-150211-final-tests`, T1 |
| Does the CPSM port reproduce `cpsm.m`? | Yes. Rel. error 5.8e-16 ('haun' and 'gardner') on ObjectsInPipe (tolerance 1e-9, set before the port was written). Envelope: max error 1.8e-14 against an analytic Gaussian burst (tolerance 1e-12) | `runs/20261123-152041-tests-with-cpsm`, T7/T8 |
| Do all algorithms focus correctly? | Yes. Peak within 0.4 px of the analytic (1 layer) or ray-traced (2 layers) scatterer position | same run, T2/T3 |
| Is PSM rounding a problem in double precision? | **No.** Every kernel variant is within 4e-13 of a 40-digit reference over 842 depth steps. The "improved" kernels give no measurable gain | `runs/20261123-145720-psm-rounding` |
| What limits PSM accuracy? | ω-quadrature discretisation (1.3e-2 on `wire` at default FFT size). It converges as the time FFT is zero-padded (m = 2: 4.0e-3) | `runs/20261123-145817-psm-quadrature` |
| How large are the legacy MULOK defects? | Envelope error vs PSM 1.4e-2 to 5.7e-2. Accurate MULOK reaches 3.6e-6 to 1.5e-2 on the same comparison | `runs/20261123-145508-accuracy-mulok-A2` |
| Best accuracy per unit cost | Accurate MULOK: within 1.2e-3 / 4.7e-4 / 9.4e-4 (wire/pins/plane) of the fine-quadrature limit, 7–21× faster than PSM at comparable accuracy (single timing runs) | `runs/20261123-150254-psm-quadrature-vs-mulok-final` |

## 1. Port consistency (T1)

The legacy ports (`legacy::psm`, `legacy::mulok`) were compared with the original `psm.m`/`mulok.m` run in GNU Octave 11.3.0 on three datasets (2D single layer, 2D three layers, 3D three layers). Relative L2 error was ≤ 1.3e-15 on every layer, and depth vectors were bitwise identical. This shows the port is **consistent** with the original. It says nothing about whether the original algorithm is correct; sections 3–4 address that.

## 2. Analytic oracles (T2, T3)

Synthetic point scatterers: an exact single-layer travel time, and a two-layer travel time computed by Fermat's principle with golden-section search. Legacy PSM, PSM, legacy MULOK (single layer only), thesis MULOK and accurate MULOK all placed the |image| maximum within 0.39–0.40 px in x and 0.23–0.39 px in z. These are the expected sub-pixel offsets of off-grid targets. This test checks *location* only. It is insensitive to amplitude or side-lobe errors, which sections 3–4 measure.

## 3. PSM

### 3.1 Rounding (T5): no improvement needed
The reference is mpmath 1.3.0 at 40 digits, evaluating eq. (3.21) from the same double-precision inputs. Over 842 recursive depth steps, the maximum error relative to the image maximum grows from about 1e-14 to at most 3.6e-13, **for every one of the 12 kernel combinations**: kz formula {MATLAB, factored} × phase update {recursive, direct, re-anchored(16)} × summation {naive, Neumaier}.
- Recursive phase update (`psm.m`) is not measurably worse than direct evaluation. This matches the derivation that forming φ = kz·z in double already costs u·|φ|, which also grows linearly with depth.
- Neumaier summation and the cancellation-free kz have no visible effect at the image level.
- Direct evaluation costs about 50% more time.

**Conclusion:** the PSM default stays at the `psm.m` arithmetic. The alternative kernels remain available as `PsmOptions::high_precision_kernels()`, documented as giving no measured benefit. Untested: single precision. Of the bundled datasets only `CylScan3D_LeadSpheres` is stored as `single`; the port was not run on it.

### 3.2 ω-quadrature: the real accuracy limit
Envelope error of PSM with the time FFT zero-padded m-fold (sum scaled by 1/m), relative to m = 8:

| dataset | m = 1 (default) | m = 2 | m = 4 |
|---|---|---|---|
| synth1 | 3.6e-6 | 4.9e-7 | 1.1e-7 |
| synth2 | 3.5e-4 | 1.4e-4 | 5.0e-5 |
| wire | 1.3e-2 | 4.0e-3 | 1.4e-3 |
| pins | 1.6e-3 | 3.5e-4 | 3.4e-4 |
| plane | 4.9e-3 | 1.9e-3 | 8.3e-4 |

The roughly 1/m convergence is consistent with integrating across the hard band and evanescent edges (`INFERRED`). The m = 8 reference is itself not exact. `t_fft_mult = 2` is a measured improvement, at twice the run time.

## 4. MULOK

### 4.1 Legacy defects: attribution (Amendment A2)
Envelope error vs PSM (m = 1), real scale fitted. Baseline: thesis MULOK with cubic interpolation, 4× oversampling, demodulated.

| dataset | baseline | + squared Jacobian | + kz truncation | + −0.75 bin ω offset | + all three | legacy `mulok.m` |
|---|---|---|---|---|---|---|
| synth1 | 3.6e-6 | 7.3e-3 | 3.6e-6 | 5.5e-3 | 8.0e-3 | 1.4e-2 |
| synth2 | 3.5e-4 | 2.3e-2 | 4.0e-4 | 1.6e-2 | 2.7e-2 | 5.7e-2 |
| wire | 1.5e-2 | 1.7e-2 | 1.7e-2 | 2.2e-2 | 2.2e-2 | 4.9e-2 |
| pins | 2.7e-3 | 5.9e-3 | 7.7e-3 | 6.2e-3 | 1.1e-2 | 3.0e-2 |
| plane | 4.4e-3 | 9.7e-3 | 1.2e-2 | 1.0e-2 | 1.6e-2 | 5.4e-2 |

Each of the three code-vs-thesis deviations (dossier C4–C6) degrades the envelope measurably on at least four of the five datasets. The ω offset also introduces a large *phase* error (complex metric about 0.4–1.0), which matters if complex images are post-processed. The remaining gap between "all three" and true legacy comes from the legacy linear interpolation without demodulation and the coarser depth grid. The pre-declared complex metric (`runs/20261123-145346-accuracy-mulok`) could not separate phase conventions from defects. It is kept as recorded; see TEST_PLAN Amendment A2.

### 4.2 Improvements (beyond the thesis, each measured)
1. **Demodulated interpolation:** interpolate W(ω)·e^{iωt_c} and restore e^{−iωt_c} exactly. This is the largest single gain. synth1 with linear interpolation at 4×: 2.4e-1 → 1.2e-4.
2. **Cubic Lagrange interpolation** instead of linear: with demodulation it reaches the reference floor at 2–4× oversampling. Lanczos(6) was worse on smooth data (5e-4 vs 3.7e-6 on synth1).
3. **Exact frequency oversampling** by zero-padding in time, instead of `fft(ifft(band))`, which removes the −0.75-bin offset and wrap-around samples.
4. **Longer depth period** (`z_fft_mult = 4`). The thesis criterion (4.28) prevents wrap-around of the record, but not of band-edge ringing. Measured floor vs PSM m = 8: wire 3.6e-3 → 1.2e-3; pins 2.6e-3 → 4.7e-4; plane 4.1e-3 → 9.4e-4.

### 4.3 Accuracy vs cost against the fine-quadrature limit (PSM m = 8)
Envelope error; all values from `runs/20261123-150254-psm-quadrature-vs-mulok-final` (final defaults), except the last column, which is from `runs/20261123-145930-psm-quadrature-vs-mulok`.

| dataset | PSM m=1 (default) | PSM m=4 (time) | accurate MULOK (time) | thesis MULOK, linear |
|---|---|---|---|---|
| synth1 | 3.6e-6 | 1.1e-7 (0.05 s) | 1.1e-7 (0.01 s) | 2.4e-1 |
| synth2 | 3.5e-4 | 5.0e-5 (0.04 s) | 3.9e-5 (0.01 s) | 3.7e-2 |
| wire | 1.3e-2 | 1.4e-3 (4.3 s) | 1.2e-3 (0.20 s) | 6.3e-2 |
| pins | 1.6e-3 | 3.4e-4 (0.21 s) | 4.7e-4 (0.03 s) | 3.1e-2 |
| plane | 4.9e-3 | 8.3e-4 (7.2 s) | 9.4e-4 (0.95 s) | 4.2e-2 |

Times are wall-clock on a 16-thread AMD Ryzen 7 5700U, from single runs (not repeated, so treat as indicative).

## 4b. CPSM defaults (TEST_PLAN amendment A3, T9)
The literal port moved to `legacy::cpsm` and still matches Octave `cpsm.m` to 5.8e-16 (T7). `cpsm::cpsm` is now corrected by default: integer angular orders, evanescent components removed, DC bin excluded, Hankel-ratio transfer function.

| Synthetic scatterer (T9) | corrected: location error (r, φ, z) [samples] | corrected: energy within ±2 | legacy: location error | legacy: energy within ±2 |
|---|---|---|---|---|
| (75 mm, +3°, 24 mm) | (0.05, 0.50, 0.00) | 0.315 | (0.05, 0.50, 0.00) | 0.280 |
| (82 mm, −5°, 30.4 mm) | (0.11, 0.50, 0.00) | 0.523 | (0.11, 0.50, 0.00) | 0.558 |

T9 **passes** (limit 1 sample). Both versions locate the scatterers equally well, and the focus measure is mixed. `OBSERVED`: on this test the corrections make the method consistent with its derivation, not visibly sharper. On the real ObjectsInPipe data the amplitude C-scan changes by 2.9% (rel. L2, each normalised): 1.9% from the three fixes and 2.3% from the transfer function. See `../synaptus-rs2/FIGURES.md`. Runs: `runs/20261123-171221-tests-A3-cpsm-defaults` (v1) and `../synaptus-rs2/runs/20261123-171236-tests-A3-cpsm-defaults` (v2).

## 5. Limitations and untested regimes
- The reference for sections 3.2 and 4 is PSM with a finer ω grid, i.e. the *same discretised model* (eq. 3.21). Agreement with it does not validate the physics: exploding-reflector model, angle-independent interface transmission (4.13), or scalar waves.
- Datasets cover 2D and 3D Cartesian, 1–3 layers, and one cylindrical scan (ObjectsInPipe). The literal CPSM is tested for consistency with `cpsm.m` (T7); the corrected default CPSM is tested for scatterer location only (T9), not for amplitude, so the dossier's amplitude question (C8) is avoided by choosing the Hankel ratio, not settled. Array PSM, the chirp-z MULOK path, CPSM's 'exact' transfer function and matched filtering (`hh`) on real data were **not** ported or tested. `hh` is implemented but has no test.
- Timing comes from single runs, not a benchmark with repeats.
- The synthetic data are idealised: noise-free, frequency-independent directivity, no attenuation.
