# Test plan v2 (written 2026-11-23, before any v2 code)

Scope: the two methods added in `synaptus-rs2` (see `../research/synaptus-v2-improvements-2026-11-23.md`):
- **A. NUFFT-PSM** (`psm_nufft`): the same discrete sum as v1 PSM, evaluated with a type-1 NUFFT.
- **B. SQ-PSM** (`psm_sq`): the continuous ω-integral J_l (dossier §5.2), by Gauss–Legendre after ω = ω* + s², with W at the nodes from a type-2 NUFFT.

All v1 tests (T1–T8) stay in this crate and must still pass. The criteria below are fixed now. Any later change is appended as a dated amendment, and the original result is kept.

Metric conventions: "rel L2" = ‖a − b‖₂/‖b‖₂ over all complex pixels of a layer. "E" = the envelope metric of v1 ACCURACY (least-squares real scale on |·|), used only in the study S1. NUFFT tolerance ε = 1e-12 inside A and B unless stated.

## N1. Type-1 NUFFT vs direct sum
- Inputs: seeded pseudo-random x ∈ [−π, π) and complex c (both components uniform in [−1, 1]); (M, N) ∈ {(100, 100), (1000, 200), (5000, 1000)}; ε ∈ {1e-6, 1e-9, 1e-12}. Oracle: direct double sum with Neumaier summation.
- **Pass (rigorous form, V6):** max_k |f̃_k − f_k| ≤ 10·ε·‖c‖₁ + 1e-13·‖c‖₁. The 10 covers ε∞ ≈ ε being a design heuristic; the second term is rounding (V9: N·u ≲ 2.2e-13).
- Reported, not pass/fail: rel L2 (heuristic V7).

## N2. Type-2 NUFFT vs direct sum
Same as N1 with ‖f‖₁ in place of ‖c‖₁.

## N3. Gauss–Legendre nodes
n = 1…200: |Σ w − 2| ≤ 1e-14·n; ∫₋₁¹ x^{2k} dx = 2/(2k+1) reproduced to relative 1e-13·n for all 2k ≤ 2n − 1; nodes strictly increasing in (−1, 1).

## N4. NUFFT-PSM ≡ v1 PSM (same discrete target)
- Datasets wire, pins, plane; t_fft_mult = 1; also pins with t_fft_mult = 4. Reference: v1 `psm` with `PsmOptions::high_precision_kernels()` (direct phase, Neumaier), same depth grid.
- **Pass:** rel L2 ≤ 1e-10 per layer, depth vectors identical. Justification: ε = 1e-12 requested; V7 is heuristic, so allow 100×; v1 is accurate to ~1e-13 (T5).

## N5. SQ-PSM self-convergence and choice of κ (node factor), synthetic data only
- Data: the v1 T2 (single-layer) and T3 (two-layer) synthetic scatterers.
- κ ∈ {1, 1.5, 2, 3, 4} against κ = 6. **Rule fixed now:** the default κ is the smallest value whose rel L2 to κ = 6 is ≤ 1e-9 on both synthetic sets. It is chosen **before** SQ-PSM is run on any real dataset. If no κ ≤ 4 qualifies, B fails N5 and no accuracy claims are made for it.

## N6. SQ-PSM vs independent quadrature oracle (Python, written separately)
- Oracle (`verify/sq_oracle.py`, NumPy + mpmath 1.3.0): reads the same input files; spatial FFT by NumPy; W(ω) by direct DTFT in NumPy; ∫ over ω by mpmath tanh–sinh (`mp.quad`) on [max(ω_lo, ω*), ω_hi] **directly in ω** (no substitution, no Gauss–Legendre, no NUFFT); then the inverse spatial DFT at the selected x.
- Pixels per layer, fixed now: the argmax of |v1 PSM| (t_fft_mult = 1) in that layer, plus (⌊nz/4⌋, ⌊nx/4⌋), (⌊nz/2⌋, ⌊nx/2⌋), (⌊3nz/4⌋, ⌊3nx/4⌋).
- Datasets: T2 and T3 synthetic, and `pins` (all 3 layers).
- **Pass:** |SQ − oracle| ≤ 1e-9 · max|SQ layer| at every selected pixel. (Oracle quadrature tolerance ~1e-12; NUFFT ε 1e-12; double rounding in the DTFT ~N_t·u·‖p‖.)

## N7. v1 PSM converges to SQ-PSM as t_fft_mult grows (falsification test)
- wire, pins, plane: E(PSM_m / m, SQ) for m = 1, 2, 4, 8, 16 (PSM_m by A, which N4 ties to v1).
- **Pass:** E(m = 16) ≤ E(m = 1)/4 on every dataset. (If SQ is the m → ∞ limit with O(1/m) error, the ratio is ≈ 1/16; 1/4 leaves room for the erratic constant noted in dossier §5.2.) Failure would mean SQ is not the limit of PSM.

## N8. Point location (analytic / ray-traced oracle, as v1 T2/T3)
A and B on the T2 and T3 synthetic data: |image| maximum within 1 pixel in x and 1 depth sample in z of the true position.

## N9. Input validation
A and B reject everything v1 T6 rejects. B also rejects `hh = Some(..)` (unsupported) with `Err`.

## S1. Study (no pass/fail): accuracy and cost
- Reference: SQ-PSM with κ = 2 × default (after N5–N7 pass).
- Methods: v1 PSM m = 1, v1 PSM m = 8, v1 legacy MULOK, v1 accurate MULOK, A with m = 1, 8, 32, B default κ.
- Metric: E per dataset (all layers concatenated, as v1 ACCURACY §4.3).
- Time: wall-clock, **median of 5 runs** after 1 warm-up, same machine, RAYON default threads. Min and max are reported as well. Labelled indicative.

## Amendment A1 (2026-11-23, before N6 was run; no N6 output had been seen)
**Change to the N6 oracle's quadrature only.** `mp.quad` evaluates the integrand one scalar at a time. For `pins` (3 layers × 128 wavenumbers × ~1000 sub-intervals) that is an estimated 10⁸–10⁹ calls, i.e. hours. Replacement, still independent of the method under test: a **composite tanh–sinh rule implemented in NumPy** (x = tanh(π/2·sinh t), t ∈ [−3.5, 3.5], step h; endpoint offsets computed as (b−a)/(1+e^{∓2u}) to keep the √ endpoint accurate), on sub-intervals of width Δω₁/2, integrating **directly in ω** (no substitution, no Gauss–Legendre, no NUFFT), with W(ω) by direct DTFT (NumPy matrix products). Convergence control: the oracle value is the level h = 1/16 result and is accepted only if it agrees with h = 1/8 to ≤ 1e-12 relative to the layer's maximum |oracle| over the checked pixels. Otherwise N6 is reported **INCONCLUSIVE** for that pixel, not passed. Pixels, datasets and the pass tolerance (1e-9 × max|SQ layer|) are unchanged. mpmath is no longer used by N6.

## N6b (added 2026-11-23 after N6 ran; separately identified follow-up, N6 verdicts unchanged)
N6 returned 21 PASS and 3 INCONCLUSIVE (oracle not converged under A1: h = 1/8 vs 1/16 differed by 1.3e-11, 1.1e-10, 6.0e-12 relative to the layer's max |oracle|). N6b re-runs **all** N6 pixels with oracle levels h = 1/16 (accepted only if it agrees with) and h = 1/32. Convergence bar (1e-12) and pass tolerance (1e-9 × max|SQ layer|) are unchanged. N6 is still reported as "not all pass".

## Amendment A2 (2026-11-23, after N5/N6/N6b; originals kept in runs/20261123-160351-v2-tests-pre-A2)
**Finding.** N5 passed only marginally (T3 layer 1: 9.8e-10 vs 1e-9, with a floor of ~3e-10 at κ = 3–4), and N6b showed SQ − oracle ≈ 5–15e-11 on T3 layer 1 while the oracle had converged to < 1e-12. Root cause (DERIVED): the node estimate (dossier §5.4) assumed the integrand's phase rate in ω is ≤ T_fft = 2π/Δω₁. The integrand contains e^{−iωt_d}, and in a shallow layer nothing compensates it, so the rate can reach t_d + n_t/f_s (T3: 70 µs vs T_fft = 51.2 µs).
**Change (v2 default from now on).** N_b := (ω_hi − a)·T_eff/(2π) with T_eff = max(T_fft, t_d + n_t/f_s). κ is re-selected with the **unchanged N5 rule**, reported as N5-A2. N6/N6b/N8 are re-run as N6-A2/N8-A2. All pre-A2 results stay recorded and are not replaced.
