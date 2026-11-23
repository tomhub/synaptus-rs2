# Accuracy report v2: NUFFT-PSM (A) and SQ-PSM (B)

Date: 2026-11-23. Labels follow the evidence labels in the root `README.md`. Every number is `OBSERVED` in the run directory cited. **Independent verification: PENDING.** The implementer wrote the methods, the tests and the Python oracle, so the oracle is independent in method (different quadrature, no substitution, no NUFFT, NumPy FFT/DTFT) but **not** in authorship.

Runs: `runs/20261123-160351-v2-tests-pre-A2` (before amendment A2), `runs/20261123-160703-v2-tests-A2` (all tests after A2), `runs/20261123-161018-v2-n7-s1` (N7, S1). `runs/20261123-160954-v2-n7-s1` is a failed start (`/usr/bin/time` missing), kept with a note.

## Summary

| Question | Answer | Evidence |
|---|---|---|
| Is A the same computation as v1 PSM? | Yes, to ≤ 2.9e-12 rel. L2 on 10 layer images (limit 1e-10) | N4 |
| Is the NUFFT correct? | Max error ≤ 8.1e-13·‖c‖₁ at ε = 1e-12 (limit 1.01e-11), types 1 and 2, 9 cases each | N1, N2 |
| Does B compute the integral it claims? | At 24 pixels of 3 cases: |SQ − oracle| ≤ 4.1e-11 of the layer maximum (most ≤ 1e-13), limit 1e-9. Pre-registered N6: 21 PASS, **3 INCONCLUSIVE** (the oracle had not converged). Follow-up N6b with a finer oracle: 24 PASS | N6, N6b (after A2) |
| Does v1 PSM converge to B as the time FFT is padded? | Yes on wire, pins, plane: E(m=16)/E(m=1) = 0.022, 0.040, 0.057 (limit 0.25). **Not monotone in m** | N7 |
| How far was v1 from the m → ∞ limit? | v1 PSM default: 1.3e-2 (wire), 1.6e-3 (pins), 4.8e-3 (plane). v1's own reference (PSM m = 8): 4.6e-4, 2.1e-4, 1.8e-4 | S1 |
| Do A and B still focus correctly? | Peak within 0.40 px of the analytic / ray-traced position (limit 1 px) | N8 |

## 1. What was planned
TEST_PLAN_V2.md: N1–N9 with pass criteria and study S1, written before any v2 code. Amendment A1 (before N6 ran) replaced `mp.quad` in the oracle by a NumPy tanh–sinh rule for run-time reasons. Amendment A2 (after N5/N6) changed B's node count; see §4.

## 2. What was executed and what passed
| Test | Result | Key numbers |
|---|---|---|
| N1 type-1 vs direct | PASS | worst max|err|/‖c‖₁ at ε = 1e-12: 8.1e-13; rel L2 ≈ 1.3–1.7·ε (heuristic, reported only) |
| N2 type-2 vs direct | PASS | worst 3.8e-13 at ε = 1e-12 |
| N3 Gauss–Legendre | PASS | moments exact to 4.7e-16·n; |Σw − 2| ≤ 1.8e-15 |
| N4 A ≡ v1 PSM | PASS | 3.6e-14 … 2.9e-12 (wire worst) |
| N5 κ choice (synthetic only) | PASS, κ = 2 | after A2: T2 5.6e-13, T3 6.1e-10 (layer 1), 1.5e-13 (layer 2); see §4 |
| N6 oracle (A1 rule) | 21 PASS, 3 INCONCLUSIVE, 0 FAIL | the inconclusive pixels are t2 (135,192), t3 L1 (72,117) and (54,192); oracle h = 1/8 vs 1/16 differed by 1.3e-11, 1.1e-10, 6.0e-12 (bar 1e-12) |
| N6b follow-up (finer oracle) | 24 PASS | worst |SQ − oracle| = 4.1e-11 of layer max (t3 layer 1); all other pixels ≤ 9.1e-14 |
| N7 PSM_m → SQ | PASS | wire 1.3e-2 → 3.0e-4; pins 1.6e-3 → 6.4e-5; plane 4.8e-3 → 2.7e-4 (m = 1 → 16) |
| N8 point location | PASS | T2 0.39/0.40 px; T3 0.23/0.40 px (A and B identical) |
| N9 input validation | PASS | — |
| v1 tests T1–T8 in this crate | PASS | unchanged |

## 3. Study S1: accuracy and cost (`OBSERVED`, run 20261123-161018)
E against the reference SQ-PSM with κ = 4 (2× default); legacy MULOK against the same reference on its own (bandwidth) depth grid. Time: median of 5 runs after a warm-up [min, max], wall-clock, 16 threads, AMD Ryzen 7 5700U.

| method | wire E | wire time | pins E | pins time | plane E | plane time |
|---|---|---|---|---|---|---|
| v1 PSM (m = 1) | 1.3e-2 | 0.70 s | 1.6e-3 | 0.028 s | 4.8e-3 | 1.80 s |
| v1 PSM (m = 8), v1's old reference | 4.6e-4 | 8.36 s | 2.1e-4 | 0.76 s | 1.8e-4 | 12.5 s |
| v1 legacy MULOK | 4.9e-2 | 0.050 s | 3.0e-2 | 0.011 s | 5.4e-2 | 0.37 s |
| v1 accurate MULOK | 8.0e-4 | 0.25 s | 6.1e-4 | 0.028 s | 8.4e-4 | 0.77 s |
| v2 A, m = 1 | 1.3e-2 | 0.026 s | 1.6e-3 | 0.009 s | 4.8e-3 | 0.26 s |
| v2 A, m = 8 | 4.6e-4 | 0.65 s | 2.1e-4 | 0.074 s | 1.8e-4 | 2.02 s |
| v2 A, m = 32 (plane: 16) | 9.0e-5 | 2.66 s | 4.8e-5 | 0.44 s | 2.7e-4 | 4.27 s |
| v2 B, κ = 2 | 7.4e-13\* | 0.15 s | 1.7e-13\* | 0.067 s | 8.7e-14\* | 1.25 s |

\* **Self-convergence, not accuracy**: the reference is B itself with twice the nodes. B's accuracy against something it does not define rests on N6/N6b (pixel oracle, ≤ 4.1e-11 of layer max) and N7 (v1 PSM converges toward it).

Observations (`OBSERVED`; interpretations `INFERRED`):
- A reproduces v1 PSM bit-for-bit in effect (N4) at 27× (wire), 3× (pins) and 7× (plane) less time for m = 1. A at m = 8 takes about the time v1 takes at m = 1 on wire.
- **v1's accuracy claims were measured against a reference that was itself 1.8–4.6e-4 from the limit.** Re-measured against the limit, v1 accurate MULOK is 8.0e-4 / 6.1e-4 / 8.4e-4 (v1 ACCURACY §4.3 said 1.2e-3 / 4.7e-4 / 9.4e-4 vs PSM m = 8). The v1 conclusion "accurate MULOK ≈ the fine-quadrature limit at a fraction of PSM's cost" survives; the decimals do not.
- Convergence in m is not monotone: pins m = 4 → 8 (1.6e-4 → 2.1e-4) and plane m = 8 → 16 (1.8e-4 → 2.7e-4). This is consistent with the band-edge offset changing erratically with m (dossier §5.2, `INFERRED`). So "just pad more" has no reliable rate, and Richardson extrapolation in m would be unjustified.
- Timings are indicative: one machine, threads not pinned, five repeats.

## 4. What failed, what was amended, what was inconclusive
- **A2 (B's node count).** N5 first passed only marginally (T3 layer 1: 9.8e-10 against 1e-9), and N6b showed SQ − oracle ≈ 5–15e-11 there while the oracle had converged. Root cause: §5.4's node estimate ignored the e^{−iωt_d} phase rate, which is uncompensated in shallow layers. Fixed by T_eff = max(T_fft, t_d + record length). Effect: at κ = 1.5, T2 1.1e-7 → 2.4e-12 and T3 layer 1 0.79 → 1.0e-9. The pre-A2 results are kept in `runs/20261123-160351-v2-tests-pre-A2`.
- **Remaining T3 layer-1 floor (~1e-10 relative to that layer).** A diagnostic (`examples/diag_t3_floor.rs`, run A2) shows it does not shrink with κ (1.2e-10 → 1.0e-10 for κ = 2 → 4 at ε = 1e-14) and barely with ε. Layer 1 holds only 3.2e-4 of layer 2's energy, so a rounding floor of ~3e-14 relative to the bright content appears as 1e-10 relative to the dark layer. `INFERRED`: conditioning, not a defect. The N5 metric (normalised per layer) cannot separate the two for near-empty layers.
- **N6: 3 INCONCLUSIVE** under the pre-registered A1 rule. N6 is reported as "not all pass". N6b is a separately identified follow-up.

## 5. What could not be tested / remains unknown
- **Oracle independence is partial.** Same author; same problem definition (dossier §5.2); same Δω₁ and nfx conventions (inherited from `psm.m`). A shared misreading of the target (e.g. the sign of t_d, or the layer-0 depth origin) would not be caught by N6. N4 (A ≡ v1) and N7 (v1 → B) argue against such a misreading: v1 PSM is tied to Octave `psm.m` by T1.
- 3D: B in 3D (plane) is covered only by N7 and S1, not by the pixel oracle.
- B's node rule is heuristic (κ chosen on two synthetic sets). Data with a long t_d, very long records or close layer velocities were not tested.
- **The model.** None of this touches the physics: angle-dependent transmission (≈6% across the beam in the MULOK paper's copper example), attenuation, elastic waves. B's 1e-13 self-consistency is precision on the model, not accuracy about the specimen. `UNKNOWN` how large the model error is on these datasets.
- Whether spectral quadrature of the migration ω-integral is published elsewhere: `UNKNOWN` (one search found nothing).
