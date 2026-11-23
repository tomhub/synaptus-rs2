# Research dossier: numerical improvements to Fourier-domain multilayer SAFT (PSM / MULOK) for a second Rust implementation
Date: 2026-11-23

Companion to `research/synaptus-fourier-saft.md` (the first dossier; its claim IDs C1–C10 are referenced as such). Claim IDs here are V1–V14.

## 1. Research question and scope
Given the first Rust port (`synaptus-rs`) and its measured accuracy limits (`synaptus-rs/ACCURACY.md`), which published methods (or derivations from published theory) can make the focused image **more accurate for the same model**, or equally accurate at lower cost, without changing the physical model? And which "improvements" would be illusory?

In scope: the ω→kz mapping of MULOK, the ω-sum ("imaging condition") of PSM, nonuniform FFT (NUFFT) algorithms, quadrature theory. Out of scope (noted only): physics model changes (angle-dependent transmission, attenuation, elastic waves), cylindrical CPSM, array (FMC) imaging.

What would answer it: a published or derivable algorithm with a stated error bound and hypotheses, plus an executable check against an oracle that does not share the algorithm. What would falsify a proposed improvement: no convergence in a self-convergence study; disagreement with an independent quadrature oracle beyond the pre-registered tolerance; no accuracy gain over PSM with a zero-padded time FFT at equal cost.

## 2. Short answer
1. `PUBLISHED-CLAIM` Replacing Stolt interpolation by a NUFFT is established practice in SAR (Subiza et al. 2003) and is suggested, but not implemented, for plane-wave ultrasound (Garcia et al. 2013, statement seen in a search snippet only). The MULOK paper itself lists an interpolation-free ω→kz mapping as future work (Skjelvareid et al. 2011, §V; full text read).
2. `DERIVED` For horizontally layered media, the per-(k⊥, layer) PSM image column is *exactly* a 1D type-1 nonuniform DFT over the ω samples. A NUFFT evaluates it to tolerance ε in O(N_ω w + N_z log N_z) per column instead of O(N_ω N_z). This gives **PSM's result at MULOK-like cost**, without interpolation or a Jacobian, which dissolves defects C4–C6 instead of fixing them.
3. `DERIVED` PSM's dominant error (ACCURACY §3.2: 1.3e-2 on `wire`) is a rectangle-rule quadrature of an ω-integral with a square-root branch point at the evanescent cut-off. Substituting ω = ω_c + s² removes the branch point. Gauss–Legendre then converges geometrically (Trefethen 2008, Thm 4.5, for analytic integrands). The spectrum at the nodes is an exact DTFT, a type-2 NUFFT. **No source applying this to migration was found** (`UNKNOWN` whether it is new).
4. `PUBLISHED-CLAIM` (full text read) Angle-dependent transmission across interfaces, which PSM/MULOK ignore, varies by about 6% within the beam in a copper example (Skjelvareid et al. 2011, p. 1041). The improvements here are numerical. They make the *model* evaluated more accurately and do not make the model more correct.

## 3. Definitions, conventions and mathematical formulation
`DEFINITION` (as in `synaptus-rs`, inherited from `psm.m`): data p(t, x, y), t = t_d + n/f_s; time FFT with MATLAB sign, W(ω, k) = Σ_n Σ_x p e^{-iω n/f_s} e^{-i k·x}, then multiplied by e^{-iω t_d}. Two-way vertical wavenumber in layer l: k_{z,l}(ω, k⊥) = √((2ω/c_l)² − k⊥²), real for ω ≥ ω_{c,l} := c_l k⊥ / 2. Band B = [ω_lo, ω_hi] = 2π[f_low, f_high]. Downward extrapolation e^{+i k_z Δz}. Evanescent components are masked cumulatively: column (ω, k⊥) contributes to layer l only if ω ≥ ω*_l := max_{j≤l} ω_{c,j}.

`DEFINITION` PSM image (layer l, depth index j, column k⊥), with z_j = s_l + jΔz_l (s_0 = the measurement z offset, s_l = 0 otherwise, as `psm.rs`):

  I_l(j; k⊥) = Σ_{ω_m ∈ B, ω_m ≥ ω*_l} W(ω_m, k⊥) e^{-iω_m t_d} Φ_l(ω_m) e^{i k_{z,l}(ω_m) z_j},   Φ_l = Π_{i<l} e^{i k_{z,i} d_i},

with ω_m = m Δω, Δω = 2π f_s / N_fft. The spatial image is the inverse 2D DFT over k⊥ of I_l, cropped.

`DEFINITION` (FINUFFT, Barnett et al. 2019, eqs. 1.1 and 1.3): type 1 f_k = Σ_j c_j e^{i k x_j}, type 2 c_j = Σ_k f_k e^{-i k x_j}, k ∈ I_N = {−N/2, …, N/2−1} (N even), x_j ∈ [−π, π) (any real x_j may be reduced mod 2π because k is an integer).

## 4. Claim ledger
| # | Claim | Label | Hypotheses / conditions | Source (exact location) | Full text inspected? | Notes |
|---|---|---|---|---|---|---|
| V1 | Stolt interpolation + final range IFFT in the ω-k algorithm can be replaced by one NUFFT | PUBLISHED-CLAIM | SAR range-migration algorithm | Subiza et al., IGARSS 2003, vol. 6, pp. 4089–4091, doi:10.1109/IGARSS.2003.1295371 | **No** (abstract/indexing only) | Accuracy figures ("<0.1% deviation") seen only in a search summary: not relied on |
| V2 | A NUFFT could replace the interpolation in f-k migration of plane-wave ultrasound; the authors used linear interpolation for simplicity | PUBLISHED-CLAIM | plane-wave transmit | Garcia et al., IEEE TUFFC 60:1853–1867, 2013, doi:10.1109/tuffc.2013.2771 | **No** (publisher blocks PMC full-text XML; statement from search snippet) | Metadata verified via Europe PMC |
| V3 | If MULOK's interpolation is not accurate enough, the image has visible artifacts; an interpolation-free ω→kz mapping (e.g. chirp-z) is future work | PUBLISHED-CLAIM | — | Skjelvareid et al., IEEE TUFFC 58(5):1037–1048, 2011, §V (Discussion), p. 1045–1046 | Yes | Also Tables I–III: PSM O(MN²), MULOK O(MN log N) in N time samples |
| V4 | Transmission factors depend on incidence angle; e.g. ≈6% variation within the beam for a 2.25 MHz, 10 mm transducer into copper | PUBLISHED-CLAIM | immersion, that transducer | Skjelvareid et al. 2011, §II, p. 1041 (citing their ref. [18]) | Yes | Model error scale; PSM/MULOK treat transmission as angle-independent |
| V5 | Type-1 NUFFT algorithm: spread with ES kernel φ_β(z)=e^{β(√(1−z²)−1)}, FFT on n ≥ σN grid, deconvolve by p_k = h/ψ̂(k); ψ̂ by Gauss–Legendre quadrature with p ≥ 1.5w+2 positive nodes | DEFINITION (algorithm) + PUBLISHED-CLAIM (quadrature sufficiency) | σ = 2; w = ⌈log₁₀ 1/ε⌉+1, β = 2.30w | Barnett, Magland, af Klinteberg, SIAM J. Sci. Comput. 41(5):C479–C504, 2019, eqs. (1.8), (3.2)–(3.10) | Yes (arXiv v2 text) | "p ≥ 1.5w+2" is from a convergence study, not proven |
| V6 | Rigorous type-1 error: max_k |f̃_k − f_k| ≤ ε∞ ‖c‖₁ | PROVEN (relied on) | exact arithmetic; ε∞ from (4.3) | ibid., eq. (4.4) and refs therein | Yes | ε∞ → requested ε is a design heuristic (Remark 10) |
| V7 | Relative ℓ2 error ≈ requested ε | PUBLISHED-CLAIM (heuristic) | "typical" data, ‖f‖₂ ≈ √N ‖c‖₂ | ibid., Remark 10, Fig. 4.1 | Yes | **No general bound when M > N** (ibid. after (4.10)). Our M (ω nodes) can exceed N (depths) |
| V8 | ES kernel aliasing error converges nearly at the Kaiser–Bessel rate in w | PROVEN (relied on) | 1D, types 1–2, fixed σ, γ<1 | Barnett, Appl. Comput. Harmon. Anal. 51:1–16, 2021 (arXiv 2001.09405); stated as Thm 7 in V5 source | Statement only (in V5 text) | |
| V9 | Rounding error of 1D NUFFT ≈ N·u relative when M ≈ N | OBSERVED (by the authors) | double precision | V5 source, Remark 9 | Yes | Our N ≤ ~10³ → ≲1e-13 |
| V10 | (n+1)-point Gauss–Legendre: if f analytic in the Bernstein ellipse E_ρ with |f| ≤ M, |I − I_n| ≤ 64M / (15(1−ρ⁻²)ρ^{2n+2}) | PROVEN (relied on; paper gives a "partial proof" and defers sharpening of constants) | f analytic in E_ρ, ρ > 1 | Trefethen, SIAM Rev. 50(1):67–87, 2008, Thm 4.5, eq. (4.14) | Yes | For an oscillatory integrand M grows with ρ: geometric convergence begins only once n exceeds the oscillation count |
| V11 | Same theorem, finite smoothness: |I − I_n| ≤ 32V / (15πk(2n+1−k)^k) | PROVEN (relied on) | f^(k) of bounded variation V | ibid., eq. (4.13) | Yes | Applies **before** substitution: √ branch point ⇒ algebraic rate only |
| V12 | Phase-shift/Fourier migration is periodic; enough zero padding in time and space prevents wraparound | PUBLISHED-CLAIM | — | S. Levin in Claerbout, *Imaging the Earth's Interior*, SEP web edition, "Stopping phase-shift migration wraparound" (1997 page) | Yes (web page) | Explains why ACCURACY §3.2 sees ~1/m convergence with t_fft_mult |
| V13 | The DFT-grid ω-sum in PSM is a rectangle-rule quadrature of ∫_B … dω/Δω with endpoints not on the grid; error O(Δω) generically | DERIVED (§5.2) | integrand bounded, piecewise smooth | — | — | Consistent with OBSERVED 1/m convergence (ACCURACY §3.2) |
| V14 | After ω = ω*_l + s², the integrand is analytic in s on the integration interval (and in a neighbourhood), so V10 applies | DERIVED (§5.3) | k⊥ fixed; ω*_l attained by a single layer; other layers' branch points lie at s² < 0 | — | — | Ellipse size limited by the nearest other branch point: an empirical convergence check is mandatory |

## 5. Theory and derivation

### 5.1 PSM column = type-1 nonuniform DFT (DERIVED)
For fixed (k⊥, l), set a_m := W(ω_m) e^{-iω_m t_d} Φ_l(ω_m) e^{i k_{z,l}(ω_m) s_l} and x_m := k_{z,l}(ω_m) Δz_l. Then I_l(j) = Σ_m a_m e^{i j x_m}, j = 0,…,N_z−1. Choose N even ≥ N_z, write j = k + N/2 with k ∈ I_N: I_l(k+N/2) = Σ_m (a_m e^{i(N/2) x_m}) e^{i k x_m}, a type-1 sum (V5 eq. 1.1) with points x_m mod 2π. ∎. The identity is exact: no interpolation, no Jacobian (the Jacobian of MULOK only arises from changing variables to a uniform kz grid). Cost per column O(M w + N log N) vs O(M N).

### 5.2 The ω-sum is a rectangle rule (DERIVED)
With t_fft_mult = m, Δω_m = Δω_1/m and the sum (scaled by 1/m, as ACCURACY §3.2) is (1/Δω_1) Σ_{ω ∈ B∩grid} F(ω) Δω_m with F(ω) = W(ω) e^{-iωt_d} Φ_l e^{i k_z z}. W(ω) here is the DTFT Σ_n p_n e^{-iωn/f_s}, independent of m (zero padding only samples it more finely). As m → ∞ this tends to (1/Δω_1)∫_{B, ω ≥ ω*} F(ω) dω. For a bounded integrand with endpoints not on the grid, the rectangle-rule error is O(Δω_m) = O(1/m). The constant also depends erratically on the fractional position of the band edges on the grid, so Richardson extrapolation in m is **not** justified. `DEFINITION` (v2 target): the continuous integral **J_l(z; k⊥) := (1/Δω_1) ∫_{max(ω_lo, ω*_l)}^{ω_hi} F(ω) dω**. This is the m → ∞ limit of the PSM sum, and V12 explains why it is the wraparound-free answer.

### 5.3 Branch point and substitution (DERIVED)
k_{z,l}² = (2/c_l)²(ω − ω_{c,l})(ω + ω_{c,l}). If ω* = ω_{c,l*} for the layer l* ≤ l with the largest velocity, then near ω* the factor e^{i k_{z,l*} d} (or e^{i k_{z,l} z} if l* = l) behaves like exp(i a √(ω−ω*)). It is continuous with an unbounded derivative, so the integrand has bounded variation but no bounded derivative, and V11 gives only an O(n⁻¹) rate. Substitute ω = ω* + s², dω = 2s ds, s ∈ [√(max(ω_lo, ω*) − ω*), √(ω_hi − ω*)]. Then (2/c)√((ω−ω*)(ω+ω*)) = (2/c) s √(2ω* + s²), which is analytic for |s|² < 2ω*. Layers i with c_i < c_{l*} have branch points at s² = ω_{c,i} − ω* < 0, i.e. on the imaginary s-axis, off the real interval. The integrand times 2s is analytic in a neighbourhood of the interval, so V10 gives geometric convergence. Caveats:
- The ellipse parameter ρ is limited by the nearest imaginary branch point, which is close to the interval when another layer has a velocity close to c_{l*}, or when ω* ≈ 0 (k⊥ ≈ 0, where √(2ω* + s²) → |s|, which is analytic anyway because it equals s for s > 0; DERIVED: for k⊥ = 0, k_z = 2ω/c is entire).
- W(ω) is entire (a finite exponential sum) but grows like e^{|Im ω| T} off the real axis, as does e^{-iωt_d}. M in V10 therefore grows with the ellipse, and convergence starts only once n exceeds roughly the number of oscillations of F over the interval, about the number of DFT bins in it (§5.4).

### 5.4 Node count (DERIVED estimate; to be confirmed empirically)
The total phase variation of F over an interval of N_b DFT bins (spacing Δω_1 = 2π/T_fft) is about 2πN_b·(T_eff/T_fft) ≤ 2πN_b. The standard requirement that Gauss–Legendre on [−1,1] integrates e^{iλx} to high accuracy once n ≳ λ/2 + O(λ^{1/3}) (heuristic; `UNVERIFIED REFERENCE` for the exact form) gives n ≳ πN_b/2 ≈ 1.6 N_b. The s-substitution makes the phase rate non-uniform (largest at the top of the band), which may require more. Default: n = ⌈κ N_b⌉ + 16, with κ fixed by a convergence study on **synthetic** data before any real dataset is evaluated.

### 5.5 Spectrum at arbitrary ω (DERIVED)
After the spatial FFT (exact, uniform), P(n, k⊥) is known for n = 0…N_t−1. Then W(ω_q, k⊥) = Σ_n P_n e^{-iω_q n/f_s}. With N = N_t (even; zero-pad by one if odd) and n = k + N/2, this is e^{-i(N/2)x_q} Σ_{k∈I_N} P_{k+N/2} e^{-ikx_q} with x_q = ω_q/f_s, a type-2 NUFFT (V5 eq. 1.3). Exact up to the NUFFT tolerance, with no interpolation of the spectrum.

## 6. Algorithm(s)
**A. NUFFT-PSM** (same discrete target as PSM with a given t_fft_mult): build the band spectrum as `psm.rs` does; for each column and layer form a_m, x_m (§5.1), call type-1 with N = N_z (rounded up to even), and extrapolate Φ to the next layer as `psm.rs` does; then run the inverse spatial FFT per depth plane. Faithful to V5 §3.1.1 except that (i) the fine grid size is the next even 5-smooth integer ≥ max(2N, 2w), and (ii) p = 2w+2 quadrature nodes for ψ̂ (≥ the V5 recommendation).

**B. SQ-PSM** (spectral-quadrature PSM; target J_l of §5.2): for each column: spatial FFT of time samples; for each layer: compute ω*_l, the s-interval, and n Gauss–Legendre nodes (§5.3–5.4); W at the nodes by type-2 NUFFT (§5.5); strengths c_q = (w_q 2s_q / Δω_1) W(ω_q) e^{-iω_q t_d} Φ_l(ω_q) e^{i k_{z,l} s_l}; image column by type-1 NUFFT (§5.1). Matched filtering (`hh`) is not supported (error).

## 7. Competing methods or interpretations
- **More time zero padding (t_fft_mult = m)** in plain PSM: converges O(1/m) (§5.2, OBSERVED in ACCURACY §3.2). Cost O(m N_ω N_z) per column. A is the same target at lower cost; B is a different (limit) target.
- **Chirp-z transform ω→kz mapping** (Lanari 1995, cited in V3 source as [25]): approximate per V3 ("approximate, but efficient"). Not read: `UNVERIFIED` beyond the citation.
- **Better MULOK interpolation** (cubic + demodulation, already in v1 `MulokOptions::accurate()`): error floor 4.7e-4 to 1.2e-3 vs PSM m = 8 (ACCURACY §4.3). Still interpolation, and still needs the Jacobian.
- **Sinc/tapered-sinc interpolation, centring by e^{iωT/2}** (Claerbout/Harlan, V12 source "Controlling Stolt migration wraparound", only the search summary was read): the demodulation idea is the same one v1 already uses.
- **Band tapering (smooth window on [ω_lo, ω_hi])**: would make the trapezoid/rectangle rule converge fast (cf. Trefethen & Weideman, SIAM Rev. 56(3):385–458, 2014, doi:10.1137/130932132: geometric convergence for smooth, effectively periodic integrands; read via abstract and secondary summary only). But it **changes the target image** (a different filter), so it is not an accuracy improvement for the same model. Not adopted.

## 8. Known limitations, failure modes and finite-precision behaviour
- V7: relative ℓ2 error of the type-1 output is a heuristic. With M > N, dark regions can carry relative errors far above ε. The rigorous statement is absolute: ≤ ε∞‖c‖₁ (V6). Metrics must be relative to the image norm or peak, not per pixel.
- Rounding: V9. N ≤ ~10³ here, so ≲1e-13.
- B's node count is heuristic (§5.4). If too few nodes are used, convergence is silent unless checked, so a mandatory self-convergence study is needed (κ vs 2κ).
- B integrates the hard band [ω_lo, ω_hi], whose edges stay non-smooth in ω. That is the model's choice. Ringing in the image is part of the target, not an error.
- **Spatial wraparound** (k⊥ sum on a DFT grid, x_fft_mult) remains in both A and B, exactly as in PSM. Neither method addresses it.
- **Model error** (V4 and the first dossier's open items: angle-independent transmission, scalar waves, exploding reflector) is untouched and plausibly percent-level, i.e. far above the numerical errors targeted here.

## 9. Reference problems and values for verification
| Problem | Parameters | Reference value | Precision quoted | Source |
|---|---|---|---|---|
| Type-1/2 NUFFT vs direct sum | random x ∈ [−π,π), complex Gaussian c, M, N ∈ {100, 1000}, ε ∈ {1e-6, 1e-9, 1e-12} | direct O(MN) sum | double | DEFINITION |
| Gauss–Legendre nodes | n = 1…200 | Σw = 2; exact for x^k, k ≤ 2n−1 | ~u | V10 (exactness property) |
| A vs PSM (v1, Direct phase kernel) | wire, pins, plane; m = 1 | v1 output | v1 accurate to ~1e-13 (ACCURACY §3.1) | first dossier / T5 |
| B vs independent quadrature | synthetic point scatterer (v1 T2 generator) and `pins`; fixed pixels | adaptive tanh–sinh quadrature (mpmath 1.3.0) of ∫ over ω directly, no substitution, W by direct DTFT in NumPy | ~1e-12 relative to the image peak | independent oracle, written separately |
| B point location | v1 T2/T3 synthetic scatterers | analytic / ray-traced position | sub-pixel | first dossier T2/T3 |

## 10. Conflicts, errata, evidence gaps and open questions
- No erratum found for V5 or V10. For V10 the paper itself says the proof is partial (constants sharpened in its §5); the geometric rate is classical (Bernstein, cited there).
- `UNKNOWN`: whether spectral quadrature of the PSM ω-integral (B) has been published. One targeted search found nothing; absence is not established.
- `UNKNOWN`: whether J_l (the m → ∞ limit) is a *better image* physically. It is the wraparound-free evaluation of the same model (V12), and nothing more is claimed.
- V1/V2 accuracy figures were not read in full text.

## 11. Implications for implementation and verification
- Keep v1 untouched; implement A and B in a separate crate (`synaptus-rs2`) starting from a copy of v1.
- Pre-register tolerances (TEST_PLAN) before any v2 run: NUFFT vs direct sum, A vs v1 PSM, B vs an independent Python quadrature oracle, B self-convergence in κ, B point-location tests.
- Report accuracy relative to the image norm (E-metric of v1 ACCURACY), never per pixel (§8).
- The claim "B is more accurate" may only be made against a target that B does not define by itself. Use the independent oracle (§9) and the monotone convergence of v1 PSM with m toward B.

## 12. References
- Barnett A. H., Magland J., af Klinteberg L. A parallel nonuniform fast Fourier transform library based on an "exponential of semicircle" kernel. SIAM J. Sci. Comput. 41(5):C479–C504, 2019. doi:10.1137/18M120885X; arXiv:1808.06736. *peer-reviewed* (full text of arXiv v2 read).
- Barnett A. H. Aliasing error of the exp(β√(1−z²)) kernel in the nonuniform fast Fourier transform. Appl. Comput. Harmon. Anal. 51:1–16, 2021. arXiv:2001.09405. *peer-reviewed* (not read; result as stated in the 2019 paper).
- Garcia D., Le Tarnec L., Muth S., Montagnon E., Porée J., Cloutier G. Stolt's f-k migration for plane wave ultrasound imaging. IEEE TUFFC 60:1853–1867, 2013. doi:10.1109/tuffc.2013.2771. *peer-reviewed* (full text not accessible; snippet only).
- Levin S. Stopping phase-shift migration wraparound. In: Claerbout J., *Imaging the Earth's Interior*, Stanford Exploration Project web edition (page dated 1997). *book (web)* (page read).
- Skjelvareid M. H., Olofsson T., Birkelund Y., Larsen Y. Synthetic aperture focusing of ultrasonic data from multilayered media using an omega-k algorithm. IEEE TUFFC 58(5):1037–1048, 2011. *peer-reviewed* (full text read, author copy at Uppsala University).
- Subiza B., Gimeno-Nieves E., Lopez-Sanchez J. M., Fortuny-Guasch J. An approach to SAR imaging by means of non-uniform FFTs. Proc. IGARSS 2003, vol. 6, pp. 4089–4091. doi:10.1109/IGARSS.2003.1295371. *peer-reviewed conference* (abstract only).
- Trefethen L. N. Is Gauss quadrature better than Clenshaw–Curtis? SIAM Rev. 50(1):67–87, 2008. doi:10.1137/060659831. *peer-reviewed* (full text read).
- Trefethen L. N., Weideman J. A. C. The exponentially convergent trapezoidal rule. SIAM Rev. 56(3):385–458, 2014. doi:10.1137/130932132. *peer-reviewed* (abstract and secondary summary only).
- Lanari R. A new method for the compensation of the SAR range cell migration based on the chirp z-transform. IEEE TGRS 33(5):1296–1299, 1995. *peer-reviewed* (not read; cited via the 2011 MULOK paper).

## 13. Limits of this research
- Three of the nine sources were read in full (FINUFFT 2019, MULOK 2011, Trefethen 2008), plus one web page (Levin). Subiza 2003, Garcia 2013, Barnett 2021, Trefethen–Weideman 2014 and Lanari 1995 were not read in full.
- The literature search for spectral quadrature of the migration ω-integral was a single query. Novelty of B is `UNKNOWN`.
- The thesis (Skjelvareid 2012) was not re-read for this dossier; conventions come from the code and the first dossier.

**Most important uncertainty:** B's accuracy rests on the derived analyticity after substitution (V14) and on a heuristic node count (§5.4). Both must be confirmed by self-convergence and by the independent quadrature oracle before any accuracy claim is made.
