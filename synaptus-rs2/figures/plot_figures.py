"""Render the arrays written by `cargo run --release --example figures`.

Only display transforms happen here (dB scaling, thresholds, colour maps); every
number comes from the Rust implementation (or, for timings, the Octave script
octave/time_originals.m). Run from the crate root:

    python3 -I figures/plot_figures.py        # needs numpy + matplotlib

Writes figures/*.png.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
GRAY = "bone"  # sequential, as in the original Synaptus figures
SERIES = {"Rust": "#2a78d6", "Octave (original)": "#eb6834"}  # validated pair


def load(name):
    dims = [int(x) for x in (DATA / f"{name}.txt").read_text().split()[1:]]
    return np.fromfile(DATA / f"{name}.f64", dtype="<f8").reshape(dims, order="F")


def db(a, ref=None):
    ref = np.max(a) if ref is None else ref
    return 20 * np.log10(np.maximum(a, 1e-300) / ref)


def zaxis(prefix, l):
    # 2D-per-slice plane images share the 3D PSM depth grid (figures.rs saves it once)
    return load(f"{prefix.replace('psm2d', 'psm3d')}_L{l}_z")


def layers(prefix, n=3):
    """Concatenate per-layer |image| (z, x, y) along z, with the z axis. Each layer's
    grid runs slightly past its bottom interface; samples at or beyond the next
    layer's first depth are dropped so z is strictly increasing."""
    imgs = [load(f"{prefix}_L{l}") for l in range(1, n + 1)]
    zs = [zaxis(prefix, l) for l in range(1, n + 1)]
    for k in range(n - 1):
        keep = zs[k] < zs[k + 1][0]
        imgs[k], zs[k] = imgs[k][keep], zs[k][keep]
    return np.concatenate(imgs), np.concatenate(zs)


def save(fig, name):
    fig.savefig(ROOT / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", ROOT / name)


def bscan(ax, img, x, z, title, rng=40, zlabel="z [mm]", zscale=1e3):
    # pcolormesh on the true coordinates: concatenated layers have different Δz,
    # so a uniform-extent imshow would misplace every layer boundary.
    im = ax.pcolormesh(np.asarray(x) * 1e3, np.asarray(z) * zscale, db(img), cmap=GRAY,
                       vmin=-rng, vmax=0, shading="nearest", rasterized=True)
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set(title=title, xlabel="x [mm]", ylabel=zlabel)
    return im


# ---------------------------------------------------------------- A. wire targets
def fig_wire():
    t, x = load("wire_t"), load("wire_x")
    fig, axs = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    bscan(axs[0], load("wire_raw_env"), x, t, "Raw data (envelope)", zlabel="t [µs]", zscale=1e6)
    bw = load("wire_psm_bw")[:, :, 0]
    bscan(axs[1], bw, x, load("wire_psm_bw_z"), "PSM, Δz = (c/2)/(f_high − f_low)\n(Python default)")
    im = bscan(axs[2], load("wire_psm_hf")[:, :, 0], x, load("wire_psm_hf_z"),
               "PSM, Δz = (c/2)/f_high\n(psm.m default)")
    fig.colorbar(im, ax=axs, label="dB", shrink=0.8)
    fig.suptitle("LineScan2D_WireTargets: raw vs focused (Rust)")
    save(fig, "wire_targets.png")

    sp, f, kx = load("wire_spec"), load("wire_spec_f"), load("wire_spec_kx")
    o = np.argsort(kx)
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    im = ax.imshow(db(sp[:, o]), aspect="auto", origin="lower", cmap=GRAY, vmin=-50, vmax=0,
                   extent=[kx[o][0], kx[o][-1], f[0] / 1e6, f[-1] / 1e6])
    ax.set(xlabel="kₓ [rad/m]", ylabel="f [MHz]", title="Wire targets: |P(ω, kₓ)| in the processed band")
    fig.colorbar(im, ax=ax, label="dB")
    save(fig, "wire_spectrum.png")


# ------------------------------------------------- B. plane scan (thesis Fig. 5.4/5.5)
def fig_plane(gate_mm=5.0):  # interface ringing spans ≈5 mm; nearest scatterer is 7 mm away
    zi = load("plane_interfaces")
    # Layer 1 is the water path above the plexiglass; scatterers are in layers 2 and 3.
    fig = plt.figure(figsize=(14, 6), constrained_layout=True)
    for k, (prefix, title) in enumerate([("plane_psm2d", "(a) 2D PSM per x–z plane"),
                                         ("plane_psm3d", "(b) 3D PSM")]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        for l in (2, 3):
            a, z = load(f"{prefix}_L{l}"), zaxis(prefix, l)
            # The interface echoes are the layer maxima; gate them out (± gate_mm around
            # each interface) so the thesis threshold, 1/5 of the maximum amplitude in
            # each layer, applies to the scatterers as in Fig. 5.5.
            top, bot = zi[l - 2], zi[l - 1]
            keep = (z > top + gate_mm * 1e-3) & (z < bot - gate_mm * 1e-3)
            a = np.where(keep[:, None, None], a, 0.0)
            iz, ix, iy = np.nonzero(a > a.max() / 5)
            ax.scatter(ix * 1.0, iy * 1.0, z[iz] * 1e3, s=4, c="#b22222", depthshade=True)
        nx, ny = a.shape[1], a.shape[2]
        xx, yy = np.meshgrid([0, nx - 1], [0, ny - 1])
        for zz in zi:
            ax.plot_surface(xx, yy, np.full_like(xx, zz * 1e3, dtype=float), color="#7fd17f", alpha=0.3)
        ax.set(xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]", title=title)
        ax.invert_zaxis()
        ax.view_init(elev=20, azim=-60)
    fig.suptitle("PlaneScan3D_PlexiAluFBH: voxels above 1/5 of the layer maximum (Rust)\n"
                 f"interface echoes gated out (±{gate_mm:g} mm); planes are layer interfaces")
    save(fig, "plane_isosurface.png")

    t = load("plane_t")
    raw = load("plane_raw_env")
    x = np.arange(raw.shape[1]) * 1e-3
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    bscan(axs[0], raw.max(axis=2), x, t, "Raw envelope, max over y", zlabel="t [µs]", zscale=1e6)
    for ax, prefix, title in [(axs[1], "plane_psm2d", "2D PSM, max over y"),
                              (axs[2], "plane_psm3d", "3D PSM, max over y")]:
        a, z = layers(prefix)
        im = bscan(ax, a.max(axis=2), x, z, title)
        for zz in zi:
            ax.axhline(zz * 1e3, color="#7fd17f", lw=0.8, ls="--")
    fig.colorbar(im, ax=axs, label="dB", shrink=0.8)
    save(fig, "plane_projections.png")


# ----------------------------------------------------- C. pipe (thesis Fig. 7.13)
def fig_pipe(mask_db=-15):
    phi, z = np.degrees(load("pipe_phi")), load("pipe_z") * 1e3
    ext = [z[0], z[-1], phi[-1], phi[0]]
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    for row, (amp, topo, lab) in enumerate([("pipe_cscan_raw", "pipe_topo_raw", "raw data"),
                                            ("pipe_cscan_im", "pipe_topo_im", "CPSM image")]):
        a = db(load(amp))
        im = axs[row, 0].imshow(a, aspect="auto", extent=ext, cmap=GRAY, vmin=-25, vmax=0)
        fig.colorbar(im, ax=axs[row, 0], label="dB")
        axs[row, 0].set(title=f"Amplitude C-scan, {lab}")
        r = np.where(a > mask_db, load(topo) * 1e3, np.nan)
        cm = plt.get_cmap("Blues_r").copy()
        cm.set_bad("#f0f0f0")
        im = axs[row, 1].imshow(r, aspect="auto", extent=ext, cmap=cm, vmin=135, vmax=150)
        fig.colorbar(im, ax=axs[row, 1], label="r [mm]")
        axs[row, 1].set(title=f"Topography C-scan, {lab} (shown where amplitude > {mask_db} dB)")
        for ax in axs[row]:
            ax.set(xlabel="z [mm]", ylabel="φ [deg]")
    fig.suptitle("CylScan3D_ObjectsInPipe: max envelope over r and its range (Rust CPSM, 'haun')")
    save(fig, "pipe_cscans.png")


# ------------------------------------------- D. pins: legacy vs accurate MULOK
def env_err(m, p):
    """ACCURACY.md E-metric: ‖|M| − α|P|‖ / ‖α|P|‖, α the least-squares real scale.
    Normalising by the *scaled* reference makes it invariant to the reference's
    amplitude convention (the t_fft_mult = 8 reference is 8× the default)."""
    alpha = np.vdot(p, m).real / np.vdot(p, p).real
    return np.linalg.norm(m - alpha * p) / np.linalg.norm(alpha * p)


def fig_pins(x_step=1e-3):  # LineScan2D_PinsPlexiAluSDH xStep
    ref, z = layers("pins_ref8")
    refbw, zbw = layers("pins_ref8bw")
    panels = [
        ("pins_ref8", None, "Reference: PSM, t_fft_mult = 8"),
        ("pins_psm", ref, "PSM (default)"),
        ("pins_legmulok", refbw, "Legacy MULOK (mulok.m port)"),
        ("pins_accmulok", ref, "MULOK, accurate options"),
    ]
    fig, axs = plt.subplots(1, 4, figsize=(18, 6), constrained_layout=True)
    for ax, (prefix, r, title) in zip(axs, panels):
        a, za = layers(prefix)
        a = a[:, :, 0]
        if r is not None:
            title += f"\nenvelope error E vs reference: {env_err(a, r[:, :, 0]):.1e}"
        im = bscan(ax, a, np.arange(a.shape[1]) * x_step, za, title)
    fig.colorbar(im, ax=axs, label="dB", shrink=0.8)
    fig.suptitle("LineScan2D_PinsPlexiAluSDH, three layers (legacy MULOK is compared with a "
                 "reference on its own Δz grid)")
    save(fig, "pins_mulok_comparison.png")


# ------------------------------------------------------------------- E. timings
def fig_timings():
    def tsv(n):
        return dict(l.split("\t") for l in (DATA / n).read_text().splitlines() if l.strip())
    rs, oc = tsv("timings_rust.tsv"), tsv("timings_octave.tsv")
    pairs = [("Wire PSM (2D)", "wire PSM (dz high-freq, as psm.m)", "wire PSM (psm.m)"),
             ("Plane PSM (3D)", "plane 3D PSM", "plane 3D PSM (psm.m)"),
             ("Pins PSM (2D, 3 layers)", "pins PSM default", "pins PSM (psm.m)"),
             ("Pins MULOK (legacy)", "pins legacy MULOK", "pins MULOK (mulok.m)"),
             ("Pipe CPSM (3D)", "pipe CPSM (haun)", "pipe CPSM haun (cpsm.m)")]
    labels = [p[0] for p in pairs]
    vals = {"Rust": [float(rs[p[1]]) for p in pairs],
            "Octave (original)": [float(oc[p[2]]) for p in pairs]}
    y = np.arange(len(pairs))
    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    h = 0.38
    for k, (name, v) in enumerate(vals.items()):
        bars = ax.barh(y + (k - 0.5) * h, v, height=h - 0.04, color=SERIES[name], label=name)
        ax.bar_label(bars, labels=[f"{s:.3g} s" for s in v], padding=3, fontsize=8, color="#333")
    speedup = [o / r for o, r in zip(vals["Octave (original)"], vals["Rust"])]
    ax.set_xscale("log")
    ax.set_yticks(y, [f"{l}\nRust ×{s:.0f} faster" for l, s in zip(labels, speedup)])
    ax.invert_yaxis()
    ax.set_xlabel("wall-clock time [s], log scale (single run each, indicative only)")
    ax.legend(loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ddd", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Same tasks, Rust port vs original Octave code")
    save(fig, "timings.png")


if __name__ == "__main__":
    fig_wire()
    fig_plane()
    fig_pipe()
    fig_pins()
    fig_timings()
