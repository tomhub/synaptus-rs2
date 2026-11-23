# synaptus-rs2

The second Rust implementation of [Synaptus](https://github.com/mh-skjelvareid/synaptus). It started as a copy of `../synaptus-rs` (v1), which stays untouched next to it. Everything v1 had is still here, with the same tests (T1–T9), and the same rule: defaults never reproduce a known bug. v2 adds two ways of evaluating phase shift migration. One is faster, one is more exact, and neither makes the physics any more right:

| Function | What it is |
|---|---|
| `psm_nufft::psm_nufft` | **A. NUFFT-PSM.** The same discrete sum as v1 PSM, but each depth column is a type-1 nonuniform FFT instead of an O(N_ω·N_z) loop. Same answer as v1 PSM (≤ 2.9e-12, test N4), so a finer ω grid (`t_fft_mult`) becomes affordable. |
| `psm_sq::psm_sq` | **B. SQ-PSM** (spectral-quadrature PSM). Computes the limit PSM converges to as the time FFT is padded without end: the continuous ω-integral. The substitution ω = ω* + s² removes the square-root branch point at the evanescent cut-off, Gauss–Legendre integrates in s, and a type-2 NUFFT gives the spectrum at the nodes (to 1e-12, no interpolation). |
| `nufft::Nufft1d` | 1D type-1/type-2 NUFFT with the "exponential of semicircle" kernel (Barnett, Magland & af Klinteberg 2019). No external dependency. |
| `quad::gauss_legendre` | Gauss–Legendre rules (Newton on the Legendre recurrence). |

**Pictures first?** [`FIGURES.md`](FIGURES.md) puts the original thesis figures next to what the original code computes and what v2 computes by default. For the interactive version, open [`figures/showcase_v2.html`](figures/showcase_v2.html): drag between the original and v2, then press "Show the bugs". The analytical deep dive is [`figures/v2_comparison.html`](figures/v2_comparison.html).

Read before relying on the output:
- `../research/synaptus-v2-improvements-2026-11-23.md`: sources, derivations and what is *not* established. B is our own derivation; we found no publication applying it to migration. That makes it either new or well hidden.
- `TEST_PLAN_V2.md`: tolerances fixed before testing, plus amendments A1 (oracle quadrature) and A2 (node count).
- `ACCURACY_V2.md`: what was measured, what passed, what was inconclusive.

## Build and test
```
cargo test --release                                   # v1 T1–T9 + v2 N1–N5, N8, N9
cargo run --release --example n6_export                # inputs for the independent oracle
python3 -I verify/sq_oracle.py verify/n6/{t2,t3,pins}         # N6 (NumPy only, tanh–sinh in ω)
python3 -I verify/sq_oracle.py --n6b verify/n6/{t2,t3,pins}   # N6b follow-up (finer oracle)
cargo run --release --example study_v2 -- <out_dir>   # N7 (PSM → SQ convergence) and S1 study
python3 -I figures/build_v2_page.py                    # figures/v2_comparison.html from the recorded runs
cargo run --release --example gallery_v2 && python3 -I figures/plot_gallery.py   # FIGURES.md images
python3 -I figures/build_showcase_v2.py                # figures/showcase_v2.html (needs gallery_v2 first)
```
v2 borrows v1's data. Build it in `../synaptus-rs` first (see the v1 README, which needs the original [Synaptus](https://github.com/mh-skjelvareid/synaptus) cloned into `../synaptus`), then link it here with `ln -s ../synaptus-rs/data data`. The link is not in the repository. The oracle needs only NumPy.

## Not addressed
- Spatial wraparound (the kₓ sum on a DFT grid, `x_fft_mult`) is the same as in v1.
- The physics: angle-independent transmission, scalar waves, the exploding-reflector model. A more exact evaluation of the same model is not a better model.
- Matched filtering in SQ-PSM; 3D is supported but was checked only through N4/N7, not by the independent oracle.

## Licence and authorship
GPL-3.0-or-later (see `../LICENSE`), because this is a derivative of [Synaptus](https://github.com/mh-skjelvareid/synaptus) © Martin Hansen Skjelvareid, which is GPL-3.0. Written by Claude Opus 5.5, guided by Tomas.

**Status:** tested by the implementer. Independent verification is pending.
