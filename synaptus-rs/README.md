# synaptus-rs

A Rust port of the Fourier-domain synthetic aperture imaging algorithms in [Synaptus](https://github.com/mh-skjelvareid/synaptus) (Skjelvareid, JOSS 2022, doi:10.21105/joss.04185). It covers PSM and MULOK for 2D/3D, single- and multi-layer media, plus cylindrical PSM (CPSM) and the analytic-signal envelope.

This is the careful one. It reproduces the original first, bugs included, and only then fixes things, in separate functions, so every difference can be measured instead of argued about. **The defaults never reproduce a known bug:** the faithful copies live under `legacy::`, and you have to ask for them by name.

| Function | What it is |
|---|---|
| `legacy::psm`, `legacy::mulok`, `legacy::cpsm` | Literal ports of `psm.m`, `mulok.m` (linear mode) and `cpsm.m`. They match Octave to ≤ 1.3e-15 relative, **bugs included**. Nothing calls them unless you do |
| `psm::psm` | Phase shift migration. The default is `psm.m` arithmetic, which measured accurate to ~1e-13. `t_fft_mult = 2` gives a measured reduction in quadrature error |
| `mulok::mulok` | MULOK as specified in the thesis (eqs. 4.15–4.17, 4.25–4.29). `MulokOptions::accurate()` adds demodulated cubic interpolation and a longer depth period; `MulokOptions::thesis()` is the plain thesis algorithm |
| `cpsm::cpsm` | Cylindrical PSM, corrected by default: integer angular orders, evanescent components removed, no DC bin, Hankel-ratio transfer function (`CpsmOptions`; Haun still available). Locates synthetic scatterers within half a sample (T9) |
| `envelope::envelope` | Envelope of the analytic signal along t (MATLAB `abs(hilbert(x))`) |

Read these before relying on the output:
- `../research/synaptus-fourier-saft.md`: what was checked in the paper, thesis and original code, including three MULOK defects and a CPSM amplitude inconsistency.
- `ACCURACY.md`: measured accuracy, attribution of the defects, and the improvements made.
- `TEST_PLAN.md`: tolerances fixed before testing, plus amendments A1 and A2.

**Later news from v2:** the accuracy reference used here (PSM with `t_fft_mult = 8`) turned out to sit 1.8e-4 – 4.6e-4 from the true limit. The conclusions in `ACCURACY.md` hold; some decimals moved. See `../synaptus-rs2/ACCURACY_V2.md` §3.

## Build and test
The tests compare against the original code, so you need it first: clone [Synaptus](https://github.com/mh-skjelvareid/synaptus) into `../synaptus` at commit `3ffbe45` (see `../README.md`).
```
cargo build --release
octave --no-gui octave/export_and_reference.m     # once: exports datasets + Octave reference outputs to data/
octave --no-gui octave/export_cpsm_reference.m    # once: ObjectsInPipe + Octave cpsm.m outputs (T7)
cargo test --release                              # T1/T7 (vs Octave), T2/T3/T8/T9 (analytic oracles), T6 (validation)
cargo run --release --example accuracy_mulok      # T4 study  (add `-- --amendment-a2` for the envelope/attribution study)
cargo run --release --example psm_quadrature      # PSM quadrature convergence and MULOK vs fine-quadrature limit
cargo run --release --example psm_rounding -- <dir> && .venv/bin/python -I verify/psm_rounding_ref.py <dir>   # T5, needs mpmath==1.3.0
```
For T5, make the environment once: `python3 -m venv .venv && .venv/bin/pip install mpmath==1.3.0`.

Run records (results, commands, environment, source hashes) are in `runs/`.

## Figures
```
cargo run --release --example figures             # computes every plotted array -> figures/data/ (+ Rust timings)
octave --no-gui octave/time_originals.m           # optional: timings of the original Octave code
python3 -I figures/plot_figures.py                # renders figures/*.png (numpy + matplotlib)
cargo run --release --example showcase_export     # Rust outputs in the Octave reference layout
python3 -I figures/build_showcase.py              # figures/showcase.html: original vs Rust, drag to compare
```
The plot script only applies display transforms (dB, thresholds, colour maps); all values come from the Rust code.

| File | Shows |
|---|---|
| `wire_targets.png`, `wire_spectrum.png` | WireTargets raw vs PSM with both depth-step conventions; the band-limited (ω, kₓ) spectrum |
| `plane_isosurface.png`, `plane_projections.png` | Thesis Fig. 5.4/5.5: 2D PSM per slice vs 3D PSM on PlaneScan3D (voxels above 1/5 of the layer maximum, interface echoes gated ±5 mm) |
| `pipe_cscans.png` | Thesis Fig. 7.13: amplitude and topography C-scans of ObjectsInPipe, raw vs CPSM |
| `pins_mulok_comparison.png` | Pins dataset: PSM, legacy MULOK and accurate MULOK against a PSM t_fft_mult = 8 reference, with the E-metric of ACCURACY.md |
| `timings.png` | Rust vs original Octave on the same tasks (single runs, indicative only) |

## Conventions
- Input: `PulseEchoData`, stored column-major with t fastest (identical to MATLAB `ptxy(:)`).
- Output: complex images, `data[(iz*nx + ix)*ny + iy]`; use `.norm()` for the envelope.
- Positive-ω band, kz ≥ 0, downward extrapolation e^{+i kz Δz}, ĉ = c/2.
- Errors: invalid input returns `Err`; nothing panics on user input.

## Not ported
Array PSM, the chirp-z MULOK path, and CPSM's 'exact' Hankel transfer function. Matched filtering (`hh`) is implemented but untested.

## Licence and authorship
GPL-3.0-or-later (see `../LICENSE`), because this is a derivative of [Synaptus](https://github.com/mh-skjelvareid/synaptus) © Martin Hansen Skjelvareid, which is GPL-3.0. Written by Claude Opus 5.5, guided by Tomas. Independent, unofficial project: not affiliated with, endorsed or reviewed by the original author(s) of Synaptus.

**Status:** tested by the implementer. Independent verification is pending.
