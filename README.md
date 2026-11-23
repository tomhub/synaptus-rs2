# Synaptus, in Rust, twice

A Rust port of [Synaptus](https://github.com/mh-skjelvareid/synaptus), the MATLAB/Octave toolbox for synthetic aperture ultrasound imaging by Martin Hansen Skjelvareid. Then a second Rust port, because matching the original to fifteen digits apparently wasn't enough.

**Who did this:** every line of code, every test, every report and both comparison pages were written by **Claude Opus 5.5** (Anthropic), **guided by Tomas**. Tomas asked the questions, set the direction, said "no" when needed, and is the reason this exists. Claude did the typing and the arithmetic. Bugs are ours. The good ideas in the physics belong to the original thesis.

## What's in the box

| Folder | What it is | One-line verdict |
|---|---|---|
| [`synaptus-rs/`](synaptus-rs/) | **v1.** PSM, MULOK and cylindrical PSM, each twice: a faithful copy of the original (bugs included, under `legacy::`) and a corrected version that is the default | The copies match the original Octave code to 1.3e-15. Boringly correct, on purpose |
| [`synaptus-rs2/`](synaptus-rs2/) | **v2.** A copy of v1 plus two new ways to compute phase shift migration: with a nonuniform FFT, and as an exact frequency integral | More exact than v1 in ways that matter less than you'd hope (see "The elephant") |
| [`research/`](research/) | Two research dossiers: what the papers claim, what the code does, where they disagree | Every claim is labelled: proven, derived, observed, or "we don't know" |

Three self-contained HTML pages tell the story with pictures. Open them in a browser:
- [`synaptus-rs/figures/showcase.html`](synaptus-rs/figures/showcase.html): original vs v1, drag a slider, find nothing.
- [`synaptus-rs2/figures/showcase_v2.html`](synaptus-rs2/figures/showcase_v2.html): original vs v2, drag a slider, then press "Show the bugs".
- [`synaptus-rs2/figures/v2_comparison.html`](synaptus-rs2/figures/v2_comparison.html): "Chasing the fourth digit", every error on one log ruler.

And one gallery that renders right here on GitHub: [`synaptus-rs2/FIGURES.md`](synaptus-rs2/FIGURES.md), with the original thesis figures next to the original code's output and v2's.

## The results, without the scroll

| Claim | Number | Where it's proven (or not) |
|---|---|---|
| v1 reproduces the original code | ≤ 1.3e-15 relative error on 14 image layers, cylindrical PSM at 5.8e-16 | [`synaptus-rs/ACCURACY.md`](synaptus-rs/ACCURACY.md), tests T1/T7 |
| The original MULOK has three bugs (below) | Each one measurably degrades the image. Together with its linear interpolation, the original MULOK lands 3.0e-2 – 5.4e-2 from the exact answer on the real datasets | [`synaptus-rs/ACCURACY.md`](synaptus-rs/ACCURACY.md) §4, [`synaptus-rs2/ACCURACY_V2.md`](synaptus-rs2/ACCURACY_V2.md) §3 |
| v2's NUFFT-PSM is the same computation as PSM | ≤ 2.9e-12 apart, and 3–27× faster | [`synaptus-rs2/ACCURACY_V2.md`](synaptus-rs2/ACCURACY_V2.md), test N4 |
| v2's SQ-PSM computes the frequency integral that PSM only approximates | Within 4.1e-11 (relative to the brightest pixel) of a separate Python oracle at 24 of 24 checked pixels. Different method, same author. It took two tries: the first was inconclusive at 3 pixels, and that's on record too | [`synaptus-rs2/ACCURACY_V2.md`](synaptus-rs2/ACCURACY_V2.md), tests N6/N6b |
| v1 graded itself against a slightly bent ruler | v1's reference was 1.8e-4 – 4.6e-4 from the limit it was meant to stand in for | [`synaptus-rs2/ACCURACY_V2.md`](synaptus-rs2/ACCURACY_V2.md) §3 |

## Things we found in the original code

Found while porting, measured, and **not yet confirmed by the original author**. None of them is reproduced by default here; the faithful copies under `legacy::` keep them on purpose, for comparison:

1. `mulok.m` applies the Stolt amplitude factor **squared** (thesis eqs. 3.25/4.16/4.26 have it once). The Python version does the same.
2. `mulok.m` starts the kz grid at 2π·f_low/ĉ instead of 0 (thesis §4.A.3), so it throws away some oblique echoes.
3. After band upsampling, `mulok.m` labels the frequency axis **0.75 bins** away from where the samples are. The Python version gets this right.
4. `psm.m` crashes on a column vector of three or more layer thicknesses. Crashing is the polite way to be wrong.
5. `cpsm.m` treats a partial scan arc as a full circle (non-integer angular orders), lets evanescent components through with zero phase, and includes the DC bin. Its amplitude factor also looks inverted compared with the thesis' own asymptotics; that one is still an open question. The corrected version changes the pipe image's brightness by about 3% and moves almost nothing.

🦝 Most of these were found by going through frequency bins, which is how the raccoon became the project's unofficial bug-finder.

Details, file locations and measured impact: [`research/synaptus-fourier-saft.md`](research/synaptus-fourier-saft.md).

## The elephant

v2 computes the model to about ten digits. The model ignores that transmission through an interface depends on angle, which the original MULOK paper estimates at about 6% across the beam in a copper example. It also assumes scalar waves, no attenuation and an exploding reflector. So we now have ten good digits about a model that may be wrong in the second. Nobody has measured how wrong it is on these datasets. That's the interesting next problem, and no quadrature rule will solve it. The badger and the raccoon both tried pushing the elephant. It did not move.

## Evidence labels

The reports tag every material statement with one of these, so you can tell a theorem from a hunch:

| Label | Meaning |
|---|---|
| `DEFINITION` | What a term or problem means here |
| `PROVEN` | A theorem with its hypotheses (cited or proved), not just a convincing plot |
| `DERIVED` | Our own derivation, steps shown, never attributed to a source |
| `PUBLISHED-CLAIM` | A source says so; not proven there or not checked here |
| `OBSERVED` | Came out of a recorded run (see `*/runs/`). True for that run, not a general fact |
| `REPRODUCED` | An earlier or published claim we re-obtained independently |
| `INFERRED` | An interpretation of the evidence |
| `HYPOTHESIS` | Plausible, unconfirmed |
| `DISPUTED` | Credible sources or results disagree |
| `UNKNOWN` | Not established. Said out loud instead of guessed |

"Two implementations agree" is never `PROVEN`. It is consistency.

## Run it yourself

You need Rust (stable), GNU Octave (for the reference outputs) and Python 3 with NumPy. Matplotlib only if you want the figures, mpmath 1.3.0 only for v1's rounding check.

```bash
# 1. The original code and data. Not included here; pinned to the commit we tested against.
git clone https://github.com/mh-skjelvareid/synaptus synaptus
git -C synaptus checkout 3ffbe45a237cb3740001f972931fc4f07590ccb2

# 2. Export the datasets and the original Octave outputs (writes synaptus-rs/data/)
cd synaptus-rs
octave --no-gui octave/export_and_reference.m
octave --no-gui octave/export_cpsm_reference.m

# 3. Test v1, then v2
cargo test --release
cd ../synaptus-rs2 && ln -sfn ../synaptus-rs/data data && cargo test --release   # v2 borrows v1's data
```

Each crate's README has the rest: accuracy studies, the independent oracle, and how to rebuild the figures and pages.

## How seriously to take this

- **Tested by the implementer.** Pass marks were written down before the tests ran, and every change after that is a dated amendment with the old result kept. Still, the same author wrote the code, the tests and the oracle. **Independent verification is pending.** Please be the independent verifier. 🦡 The badger would like a colleague.
- Every number traces to a run directory under `*/runs/` with the command, environment and file hashes. Home paths and the hostname in those records have been redacted to `<repo>`, `~` and `<host>`.
- Timings come from one laptop. They're indicative, not a benchmark.
- "Matches the original" means *consistent*, not *correct*. Faithfully reproducing a bug is still a bug. The analytic point-scatterer tests and the oracle are what speak to correctness.

## Licence

**GPL-3.0-or-later**, see [`LICENSE`](LICENSE). MIT was requested, but this is a port of Synaptus, which is GPL-3.0, and a derivative of GPL code stays GPL. Synaptus is © Martin Hansen Skjelvareid and contributors.

## Credits and references

- Skjelvareid, M. H. (2022). *Synaptus: A Matlab/Octave toolbox for synthetic aperture ultrasound imaging.* JOSS 7(76), 4185. doi:10.21105/joss.04185. Code: <https://github.com/mh-skjelvareid/synaptus>
- Skjelvareid, M. H., Olofsson, T., Birkelund, Y., Larsen, Y. (2011). Synthetic aperture focusing of ultrasonic data from multilayered media using an omega-k algorithm. *IEEE TUFFC* 58(5), 1037–1048.
- Skjelvareid, M. H. (2012). *Synthetic aperture ultrasound imaging with application to interior pipe inspection.* PhD thesis, University of Tromsø. hdl:10037/4649.
- Barnett, A. H., Magland, J., af Klinteberg, L. (2019). A parallel nonuniform fast Fourier transform library based on an "exponential of semicircle" kernel. *SIAM J. Sci. Comput.* 41(5), C479–C504.
- Trefethen, L. N. (2008). Is Gauss quadrature better than Clenshaw–Curtis? *SIAM Review* 50(1), 67–87.

If this helps your work, cite the original Synaptus paper first. It did the hard part.
