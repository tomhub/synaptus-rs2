# Research dossier: Synaptus, Fourier-domain synthetic aperture ultrasound imaging (PSM, MULOK, CPSM, array PSM)

> Paths such as `matlab/core/psm.m`, `python/synaptus/…` and `docs/thesis.pdf` refer to the original repository, <https://github.com/mh-skjelvareid/synaptus>, at commit `3ffbe45a237cb3740001f972931fc4f07590ccb2`. It is not included here.


Date: 2026-11-23
Repository examined: `synaptus/` at commit `3ffbe45` (origin `github.com/mh-skjelvareid/synaptus`)
Evidence labels follow the evidence labels in the root `README.md`.

## 1. Research question and scope

Do the Synaptus implementations (MATLAB/Octave `matlab/core`, `matlab/learn`, Python `python/synaptus`) correctly implement the algorithms described in the JOSS paper (Skjelvareid 2022) and the PhD thesis it relies on (Skjelvareid 2012, `docs/thesis.pdf`)? Are the claims in the paper supported? This dossier is the basis for a Rust port and for numerical-accuracy improvements.

In scope: PSM (single- and multi-layer, 2D/3D), MULOK (multilayer Stolt / ω-k), CPSM (cylindrical PSM), and array PSM, checked against the thesis equations; sign and sampling conventions; numerical aspects.
Out of scope: re-deriving acoustic transmission physics at interfaces; checking the experimental results on the physical test blocks; the chirp-z MULOK variant beyond reading it.

## 2. Short answer

- The PSM implementation (`psm.m`) is **consistent with thesis eqs. (3.17)–(3.21), (4.12), (4.14)** after accounting for the FFT sign convention (`DERIVED`; independent numerical confirmation pending in the Rust phase).
- MULOK deviates from the thesis in three ways:
  1. **The Stolt amplitude factor is squared**: the code uses kz²/(kz²+kx²), but thesis eqs. (3.25), (4.16) and (4.26) give the Jacobian ∝ kz/√(kz²+kx²). This holds in MATLAB core, MATLAB learn and Python. `PUBLISHED-CLAIM` (thesis) vs. code; the thesis form is `DERIVED`-confirmed below.
  2. **The kz grid is truncated**: the code uses kz ∈ [2πf_low/ĉ, 2πf_high/ĉ), but thesis eq. (4.27) specifies kz starting at 0. This discards spectrum for kx ≠ 0 (`DERIVED`).
  3. **MATLAB only: the frequency axis after band upsampling is offset by −(1−1/up)·Δω = −0.75 Δω** (`OBSERVED`, `research/checks/runs/*/check_upsampling_grid.log`). The Python port does not have this defect.
- CPSM's "Haun" transfer function has amplitude √(r/R). The exact Hankel ratio it approximates has asymptotic amplitude √(R/r). Both the thesis (6.7 vs 6.9) and the code contain this inconsistency (`DERIVED`). It affects range gain only, not focusing.
- Python PSM uses a coarser depth step than MATLAB PSM ((c/2)/(f_h−f_l) vs (c/2)/f_h). This is a sampling difference, not an error in the computed pixel values (`DERIVED`).
- The image *values* produced by each MULOK variant have not yet been quantified against PSM. That comparison is planned in §11. `UNKNOWN` until then.

## 3. Definitions, conventions and mathematical formulation

- Exploding reflector model: effective velocity ĉ = c/2 (thesis §3.2; code `cc/2`). `DEFINITION`
- Thesis transform convention (eq. 3.15/3.16): inverse uses e^{+i kx x} e^{−iωt}; forward uses e^{+iωt}. MATLAB/NumPy `fft` uses e^{−iωt} in the **forward** transform, so FFT bin ω corresponds to thesis frequency −ω. `DEFINITION`
- Dispersion: kz = −sgn(ω)·√(ω²/ĉ² − kx² [− ky²]) (thesis 3.19). Evanescent components (negative radicand) are set to 0 (thesis 4.A.2). `DEFINITION`
- Extrapolation: P(ω,kx,z) = P(ω,kx,Z)·e^{i kz (z−Z)} (3.17). Imaging condition t = 0 gives i(x,z) = ∬ P(ω,kx,Z) e^{i kz (z−Z)} e^{i kx x} dkx dω (3.21). `PROVEN` (exact under the scalar-wave, source-free-half-space and exploding-reflector assumptions, which the thesis states)
- **Code convention (DERIVED).** `psm.m` keeps FFT bins ω ∈ [+2πf_l, +2πf_h], uses kz = +√(4ω²/c² − kx²) and multiplies by e^{+i kz Δz}. With the FFT's e^{−iωt}, an echo from depth z carries phase e^{−iω·2z/c} ≈ e^{−i kz z} in the (ω,kx) domain. Multiplying by e^{+i kz z} cancels it, so the code focuses correctly. The MATLAB `learn_*` scripts, `mulok.m` and `cpsm.m` use negative FFT bins with e^{−i|kz|Δz}, which is equivalent. Python uses positive bins like `psm.m`.
- Time-delay compensation: data sample n is at t = t_delay + n/fs, so the true spectrum is P_meas(ω)·e^{−iω t_delay} with the FFT convention. All implementations apply exactly this. `DERIVED`
- Multilayer: P(ω,kx,Z_l) ∝ P(ω,kx,Z_1)·exp(i Σ_{m<l} k_{z,m} d_m) (4.14). This assumes angle-independent interface transmission, as the thesis states explicitly in (4.13). `PUBLISHED-CLAIM` (an approximation with stated hypotheses)
- Stolt (3.22)–(3.25): ω(kz,kx) = −sgn(kz)·ĉ·√(kx²+kz²), and P(kz,kx) = A·P(ω(kz,kx),kx) with A = ∂ω/∂kz = ĉ/√(1+kx²/kz²).
  **DERIVED check:** ∂/∂kz [ĉ√(kx²+kz²)] = ĉ·kz/√(kx²+kz²) = ĉ/√(1+kx²/kz²) for kz > 0. The thesis Jacobian is correct. With this A, the Stolt image equals the PSM image (3.21) exactly, as a change of integration variable over the propagating region.
- Cylindrical (6.3)–(6.7): exact transfer G = H_n^{(2)}(k_r r)/H_n^{(2)}(k_r R), with k_r = √((ω/ĉ)² − kz²). `PROVEN` (separation of variables, inward waves only, sources at r > R)

## 4. Claim ledger

| # | Claim | Label | Hypotheses / conditions | Source (location) | Full text inspected? | Notes |
|---|---|---|---|---|---|---|
| C1 | PSM image (3.21) is the exact back-propagation focusing of the exploding-reflector field | PROVEN | scalar wave eq.; constant ĉ in the layer; sources only below the recording plane; exploding-reflector approximation (ĉ = c/2) | Thesis §3.5.2–3.5.3; Gazdag 1978 | Thesis yes; Gazdag no | Standard result |
| C2 | Recursive phase-shift (multiplying by e^{i kz Δz} repeatedly) is equivalent to direct evaluation | PROVEN (exact arithmetic) | Exact arithmetic only. In floating point, rounding accumulates (§8) | Thesis §3.5.3 citing Olofsson et al. 2010 | Thesis yes; Olofsson no | Accuracy cost is quantified in the Rust phase |
| C3 | Stolt Jacobian A = ĉ/√(1+kx²/kz²) | PROVEN / DERIVED | kz ≠ 0, constant ĉ | Thesis eqs. 3.25, 4.16, 4.26 | yes | Re-derived in §3 |
| C4 | `mulok.m`, `learn_stolt.m` and Python `_stolt_transform` implement C3 | **REFUTED by inspection** | — | `mulok.m` (`Akzkxky = 1./(1 + (KX.^2+KY.^2)./KZ.^2)`); `learn_stolt.m`; `multilayer_cartesian_migration.py` (`Akzkx = 1/(1 + ...)`) | yes | The code factor equals A² up to a constant |
| C5 | Stolt kz grid should start at kz = 0 and reach 2πf_max/ĉ (4.27–4.29) | PUBLISHED-CLAIM + DERIVED | Propagating components with kx ≠ 0 have kz ∈ (0, 2ωH/c] | Thesis §4.A.3 | yes | The code grid starts at 2πf_low/ĉ, discarding kz < 2πf_low/ĉ for kx ≠ 0 |
| C6 | `mulok.m` band upsampling assigns the correct ω grid | **REFUTED** | up = 4 (default) | `mulok.m` lines building `omegaBand` after `fft(ifft(Poxy),up*n)`; also `learn_stolt.m` | yes | `OBSERVED` offset −0.750000 Δω; interpolation error vs the analytic spectrum is 1.9e-1 on the code grid vs 2.7e-2 on the true grid |
| C7 | Python MULOK band upsampling is correct | DERIVED | — | Python zero-pads in time (`nfft_t *= omega_upsampling_factor`) | yes | Correct by construction |
| C8 | Haun approximation Ĝ ≈ √(r/R)·exp(i(r−R)√(ω²/ĉ² − n²/(Rr) − kz²)) | PUBLISHED-CLAIM | narrow beam, cos(φ−φ₀) ≈ 1 − (φ−φ₀)²/2 | Thesis eq. 6.9, citing Haun et al. 2002 | Thesis yes; **Haun 2002 not retrieved** | Its amplitude is the reciprocal of the asymptotic \|G\| ≈ √(R/r) from (6.7). DISPUTED internally. It affects only range-dependent gain |
| C9 | Gardner–Collin closed form used in `approxHankel2` | DERIVED (matches the Debye asymptotic) | \|z\| > 1.1\|n\| (code cut-off) | `cpsm.m` subfunction; Gardner & Collin 2000 | Code yes; **paper not retrieved** | Phase = √(z²−n²) − n·arccos(n/z) − π/4, written with atan. Also correct for n < 0 via H_{−ν}^{(2)} = e^{−iπν}H_ν^{(2)} |
| C10 | CPSM angular orders are integers (Fourier series over 2π, eq. 6.4) | DEFINITION (thesis) | full-circle scan | Thesis 6.4 | yes | The code zero-pads φ to 2^k over a **partial arc** (CylScan3D: 216 × 0.00503 rad ≈ 62°). The orders become non-integer (k·4.88). This is an approximation, not stated in the thesis |
| C11 | Lateral resolution ≈ D/2, independent of depth | PUBLISHED-CLAIM | planar uniformly excited piston; −6 dB beamwidth λ/D; small angle; δx < D/4; single layer, or multilayer with angle-independent transmission (§4.3) | JOSS paper "Example use"; thesis 3.28–3.29, §4.3 | yes | The paper states depth independence without the multilayer caveat that the thesis gives in §4.3 |
| C12 | The JOSS example satisfies the sampling criterion | DERIVED | δx = 1 mm, D = 6 mm → D/4 = 1.5 mm | Paper; dataset `PlaneScan3D_PlexiAluFBH.mat` (xStep = yStep = 1e-3, transDiam = 6e-3) | yes | Consistent |
| C13 | Paper velocities 1480/2730/6320 m/s match the data | OBSERVED | — | Paper; `.mat` `cc` | yes | Matches |
| C14 | Fourier-domain processing is "in many cases more computationally efficient" than delay-and-sum | PUBLISHED-CLAIM | the thesis assumes a constant synthetic-aperture size | Paper; thesis §4.4, Tables 4.1–4.2 | yes | PSM is O(MN²) vs O(NM) for MLDAS at constant aperture (thesis's own analysis), so PSM is **not** asymptotically faster. MULOK is O(LMN log MN). The paper's statement is qualified ("in many cases") and not contradicted, but the thesis shows it depends on the algorithm and the aperture assumption |
| C15 | `psm.m` accepts the layer-thickness vector in any orientation | **REFUTED** | column `thick` with ≥ 3 layers | `psm.m` `nZ = ceil(sum([thick(1)-zOffset; thick(2:end).']...` | yes | `OBSERVED` Octave error "vertical dimensions mismatch". It crashes rather than giving a wrong result |
| C16 | MULOK NaN in the amplitude factor ("TODO: Check NaN - why?") | DERIVED | — | `mulok.m` | yes | 0/0 occurs only where KZ = KX = 0, which requires f_low = 0. `mulok.m` sets it to 0 and `learn_stolt.m` sets it to 1 (inconsistent). It is harmless when f_low > 0 |

## 5. Theory and derivation

### 5.1 Why the code's Stolt weight is A² (DERIVED)
With u = kx²/kz²: the code uses 1/(1+u), and thesis (3.25) gives ĉ/√(1+u). Then 1/(1+u) = (√(1+u))^{−2} = (A/ĉ)². The code therefore weights each (kz,kx) component by an extra factor cos θ = kz/√(kx²+kz²), where θ is the propagation angle. The effect is to attenuate steep-angle components (large |kx|), which carry the lateral resolution. PSM applies no such weighting, so MULOK-as-coded cannot equal PSM even with perfect interpolation.

### 5.2 kz truncation (DERIVED)
For fixed kx, ω ∈ [ωL, ωH] maps to kz ∈ [√(max(ωL²/ĉ² − kx², 0)), √(ωH²/ĉ² − kx²)]. The lower end lies below ωL/ĉ whenever kx ≠ 0, and the derivative K/√(K²−kx²) > 1 means the mapped band is also wider than the kx = 0 band. The code's grid [ωL/ĉ, ωL/ĉ + 2π/dz) with dz = ĉ/(f_h−f_l) has exactly the kx = 0 width. So for |kx| > 0 the part kz < ωL/ĉ is dropped, while the top edge is preserved. Thesis (4.27)–(4.29) avoid this by starting the grid at 0.

### 5.3 Upsampling grid (DERIVED, then OBSERVED)
`fft(ifft(x), M)` with M = up·n returns samples X'_j at ω_1 + j·Δω/up (j = 0…M−1), with X'_{up·k} = x_k. The code assigns ω_code(j) = ω_n − (M−1−j)·Δω/up = ω_true(j) − (1 − 1/up)·Δω. In addition, the last up−1 samples interpolate across the circular wrap from x_n back to x_1, and the code places them inside the band.

### 5.4 Depth sampling (DERIVED)
For PSM, the image's kz support at a given kx is ⊂ [0, 2ωH/c], a band of width W = 2ωH/c = 4πf_h/c. Complex sampling without aliasing therefore needs dz ≤ 2π/W = (c/2)/f_h, which is `psm.m`'s choice. Python PSM's (c/2)/(f_h−f_l) under-samples steep-angle content. The pixel values it does compute are still exact; only interpolation between depth samples and peak capture suffer (scalloping of |image|).

## 6. Algorithms (faithful to source; deviations flagged)

**PSM (thesis Fig. 3.9, §4.A.4; `psm.m`)**
```
P ← FFT_t(p, nFFTt) · conj(H) · e^{-iω t_delay};  keep ω ∈ [2πf_l, 2πf_h]
P ← FFT_x[,y](P)
for layer l: KZ_l = sqrt(max(4ω²/c_l² − kx² − ky², 0)); P ← P·[KZ² ≥ 0]
   Ps ← P·e^{i KZ_l·z_start}           (z_start = c_1 t_delay/2 in layer 1, else 0)
   for each depth j: img_j = IFFT_x[,y]( Σ_ω Ps );  Ps ← Ps·e^{i KZ_l dz_l}   (recursive)
   P ← P·e^{i KZ_l d_l}
```
**MULOK (thesis eqs. 4.15–4.17, §4.A.3)**: per layer, interpolate P(ω,kx) at ω = ĉ_l√(kx²+kz²) on an equispaced kz grid, multiply by A_l, IFFT over (kz,kx), crop to the layer, then extrapolate to the next interface as in PSM. Deviations in the code: C4, C5, C6.

**CPSM**: per range r, multiply by G (exact Hankel ratio, Gardner approximation, or Haun approximation), sum over ω, IFFT over (n, kz).

**Array PSM**: one-way kz = √(ω²/c² − kx²). The up-going (received) and down-going (modelled transmit) fields are extrapolated jointly, and the imaging condition is Σ_ω U·D per depth. The thesis does not cover this algorithm (the code comment says so), so it has no published derivation to check against. `UNKNOWN` beyond internal consistency.

## 7. Competing methods or interpretations

- **Stolt Jacobian conventions in the literature.** Some seismic implementations deliberately omit or alter the Jacobian, or apply obliquity factors for amplitude-preserving migration. Using a different factor is a legitimate *design choice* only if it is documented. Here the code contradicts its own reference document (thesis 3.25/4.16/4.26) with no comment. `INFERRED`: this looks like an implementation slip, not a choice.
- **Stolt interpolation.** Linear interpolation (as coded) vs. sinc or band-limited interpolation, chirp-z (code option), or NUFFT-based approaches. Interpolation quality controls the artefact level (thesis §3.5.4). The NUFFT alternative is not discussed in the repository.
- **PSM recursion vs. direct exponentials.** The thesis chooses recursion for speed. Direct computation of e^{i kz z_j} per depth costs more transcendental evaluations but has no accumulated error.

## 8. Known limitations, failure modes and finite-precision behaviour

- Recursive phase shift: after J steps, the computed e^{iφ}^J carries a relative error of roughly J·u (u = unit round-off) in both phase and magnitude. For J ~ 10³ in double precision this is ~1e-13, which is negligible. In single precision it would be ~1e-4. All bundled datasets except CylScan3D are double; CylScan3D's `ptpz` is **single** precision, and MATLAB/Octave propagate `single` through the whole CPSM computation. `DERIVED`; effect size to be `OBSERVED`.
- kz² = 4ω²/c² − kx² suffers cancellation near the evanescent cut-off. The factored form (2ω/c − |kx|)(2ω/c + |kx|) has small relative error. The error in kz·z near the cut-off is amplified by 1/kz. `DERIVED`
- A hard rectangular band mask in ω and a hard evanescent mask cause Gibbs ringing (a design choice, not a defect).
- Circular convolution in x/y and in z (MULOK): zero-padding only to the next power of two. MULOK's z FFT length ignores t_delay. Wrap-around is possible if t_delay ≫ window length. `DERIVED`; checked not to occur for the shipped datasets.
- `psm.m`: the last image plane of each layer extends past the interface (ceil(...)+1 planes) and is computed with the upper layer's velocity. This creates slight overlap with the next layer.
- The multilayer models assume angle-independent transmission (4.13), i.e. narrow beams.

## 9. Reference problems and values for verification

| Problem | Parameters | Reference value | Precision | Source |
|---|---|---|---|---|
| Point scatterer, single layer, synthetic | p(t,x) = w(t − 2√((x−x₀)²+z₀²)/c), Gaussian-modulated pulse | Focus peak at (x₀, z₀), within one pixel | exact by construction | DERIVED (exploding-reflector model) |
| PSM ≡ direct (3.21) evaluation | any data | Identity up to rounding | ≲ J·u·‖·‖ | C2 |
| MULOK (exact Jacobian, full kz grid, accurate interpolation) ≈ PSM at the same (x,z) samples | synthetic point scatterer | Relative difference → 0 as interpolation improves | interpolation-limited | §3, C3 |
| Upsampling grid identity | analytic spectrum e^{−iωτ} | X'_{up·k} = x_k | 7.9e-16 observed | `research/checks` |
| WireTargets dataset | fs = 25 MHz, c = 1480 m/s, t_delay = 40 µs, x_step = 0.5 mm, band 0.4–2.5 MHz (Python test) | none published numerically | — | `datasets/` |
| PinsPlexiAluSDH | 3 layers, c = [1480 2730 6320], thick = [73.3 31 54] mm | thesis §4.5.1 figures only | qualitative | thesis |

No numerical reference *images* are published, so dataset comparisons can only check consistency between implementations, not correctness.

## 10. Conflicts, errata, evidence gaps and open questions

1. Code vs. thesis Stolt Jacobian (C4), kz grid (C5), upsampling grid (C6). I found no erratum for the thesis, and no issue in the repository that discusses these. A GitHub issue search was not performed (offline).
2. Haun amplitude √(r/R) vs. the exact |G| asymptote √(R/r) (C8). This cannot be resolved without Haun et al. 2002. **Open.**
3. Non-integer angular orders in CPSM for partial-arc scans (C10). This is not discussed in the thesis.
4. The practical image-quality impact of C4–C6 is not yet quantified. **Open; planned.**
5. The array-PSM imaging condition and normalisation have no published derivation in this repository.

## 11. Implications for implementation and verification

- The Rust port should implement the algorithms **as published in the thesis** (correct Jacobian, kz grid from 0, correct upsampling). It should also offer a `legacy` mode reproducing the MATLAB behaviour exactly, so that differences can be attributed and quantified rather than asserted.
- Planned numerical improvements, each to be measured rather than assumed:
  (a) direct per-depth phase computation, or periodically re-anchored recursion, instead of pure recursion;
  (b) the cancellation-free kz² form;
  (c) compensated (Neumaier) or pairwise summation over ω;
  (d) Stolt interpolation better than linear (e.g. band-limited or cubic), judged against PSM as the interpolation-free reference.
- Oracles for the independent verifier: synthetic point scatterers (analytic focus location), PSM-vs-direct (3.21), and MULOK-vs-PSM convergence. Cross-implementation agreement (Rust vs. Octave) is consistency only.

## 12. References

- Skjelvareid, M. H. (2022). Synaptus: A Matlab/Octave toolbox for synthetic aperture ultrasound imaging. *Journal of Open Source Software*, 7(76), 4185. doi:10.21105/joss.04185. **peer-reviewed.** Inspected: the repository source `paper/paper.md`, not the published PDF. DOI not resolved online in this session.
- Skjelvareid, M. H. (2012). *Synthetic aperture ultrasound imaging with application to interior pipe inspection.* PhD thesis, University of Tromsø. hdl:10037/4649. **thesis.** Full text inspected (`docs/thesis.pdf`, §§3.5, 3.6, 4.2, 4.3, 4.4, 4.A, 6.1).
- Stolt, R. H. (1978). Migration by Fourier transform. *Geophysics*, 43(1), 23–48. **peer-reviewed.** Not retrieved; cited as in the thesis.
- Gazdag, J. (1978). Wave equation migration with the phase-shift method. *Geophysics*, 43(7), 1342–1351. **peer-reviewed.** Not retrieved.
- Olofsson, T. (2010). Phase shift migration for imaging layered objects and objects immersed in water. *IEEE Trans. UFFC*, 57(11), 2522–2530. **peer-reviewed.** Not retrieved.
- Olofsson, T., Skjelvareid, M. H., Barkefors, A. (2010). Ultrasonic imaging of immersed objects using migration techniques. *EUSAR 2010*, 442–445. **conference.** Not retrieved.
- Haun, M. A., Jones, D. L., O'Brien Jr., W. D. (2002). Efficient three-dimensional imaging from a small cylindrical aperture. *IEEE Trans. UFFC*, 49(7), 861–870. **peer-reviewed.** Not retrieved. The thesis lists the third author as "Jr., W. D. O."; the expansion "O'Brien" is my reading of that abbreviation and is UNVERIFIED.
- Gardner, J., Collin, R. E. (2000). An accurate closed-form approximate representation for the Hankel function of the second kind. *IEEE Trans. Antennas Propag.*, 48(10), 1699–1700. **peer-reviewed.** Not retrieved.

## 13. Limits of this research

- No external source was retrieved online. All primary-literature statements come via the thesis's citations and are therefore only as reliable as the thesis.
- Haun 2002 and Gardner & Collin 2000 were not read. C8 and the accuracy claims for C9 remain open.
- The image-level impact of the MULOK defects is derived, not yet measured.
- The chirp-z MULOK path, `polFocusedImage`, `polRawDataImage` and the Jupyter notebooks were not examined in detail.

**Most important open uncertainty:** how much C4–C6 actually degrade MULOK images relative to PSM on realistic data. The defects are certain, but their practical significance is `UNKNOWN` until it is measured.

---

## Addendum (2026-11-23, after implementation and measurement)

The questions left open in §10 item 4 and §2 are now measured (`synaptus-rs/ACCURACY.md`, run directories under `synaptus-rs/runs/`). Labels are updated; the original text above is left unchanged.

| # | Update | New label |
|---|---|---|
| C2 | Recursive phase shift in double precision: error ≤ 3.6e-13 of the image maximum after 842 steps, the same as direct evaluation (40-digit mpmath reference). The recursion is numerically harmless in double | OBSERVED |
| C4 | The squared Jacobian degrades MULOK envelopes on all 5 test datasets (e.g. synth2: 3.5e-4 → 2.3e-2 vs PSM) | OBSERVED |
| C5 | kz truncation degrades envelopes on 4 of 5 datasets (pins 2.7e-3 → 7.7e-3; plane 4.4e-3 → 1.2e-2). No effect on synth1, which has no low-kz energy | OBSERVED |
| C6 | The −0.75 Δω offset degrades envelopes 1.5–1500× above the accurate baseline (5 of 5 datasets) and causes large phase errors | OBSERVED |
| new C17 | Thesis eq. 4.28 (δkz ≤ 2πfs/(ĉN_t)) prevents wrap-around of the recorded signal, but not of band-edge ringing. A 4× longer depth period lowers the MULOK error vs the fine-quadrature limit 3–6× on real data | OBSERVED |
| new C18 | PSM's dominant error is ω-quadrature discretisation (1.3e-2 envelope on `wire` at default FFT size), converging roughly as 1/m with time zero-padding | OBSERVED (convergence rate INFERRED) |
| new C19 | Thesis MULOK ≡ PSM (an exact change of variables) is confirmed numerically: with accurate interpolation, MULOK reaches the PSM fine-quadrature limit to 1.1e-7 (synthetic) | OBSERVED, consistent with C3 |

These remain open: C8 (Haun amplitude; paper not retrieved), C10 (CPSM non-integer orders), CPSM amplitude accuracy (CPSM was ported later the same day: the literal port, now `legacy::cpsm`, reproduces C8 and C10; the default `cpsm::cpsm` uses integer orders, which resolves C10, and the Hankel ratio, which sidesteps C8 without settling it; see TEST_PLAN amendment A3 and T9), and array PSM, which was not ported.
