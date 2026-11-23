"""Render the v2 gallery (figures/gallery/*.png) from figures/data/gallery/.

Only display transforms (dB, log error scales); all values come from
`cargo run --release --example gallery_v2`. Run from the crate root:
    python3 -I figures/plot_gallery.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
D = ROOT / "data/gallery"
OUT = ROOT / "gallery"
OUT.mkdir(exist_ok=True)
GRAY, ERR = "bone", "magma"


def L(n):
    dims = [int(x) for x in (D / f"{n}.txt").read_text().split()[1:]]
    return np.fromfile(D / f"{n}.f64", "<f8").reshape(dims, order="F")


def C(n):
    return L(n + "_re") + 1j * L(n + "_im")


def db(a):
    return 20 * np.log10(np.maximum(a, 1e-300) / a.max())


def save(fig, name):
    fig.savefig(OUT / name, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def mesh(ax, img, x_mm, z_mm, **kw):
    im = ax.pcolormesh(x_mm, z_mm, img, shading="nearest", rasterized=True, **kw)
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    return im


def monotone(z):
    """Indices keeping z strictly increasing (layer joins overlap by one sample)."""
    keep = [0]
    for i in range(1, len(z)):
        if z[i] > z[keep[-1]]:
            keep.append(i)
    return np.array(keep)


# ---------- 1. pipe: original CPSM vs corrected CPSM ----------
phi, z = np.degrees(L("pipe_phi")), L("pipe_z") * 1e3
o, f = L("pipe_orig_cscan"), L("pipe_fixed_cscan")
on, fn = o / o.max(), f / f.max()
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
for ax, a, t in [(axs[0], on, "Original cpsm.m (via legacy::cpsm)"), (axs[1], fn, "v2 default: corrected CPSM")]:
    im = ax.pcolormesh(z, phi, db(a), cmap=GRAY, vmin=-25, vmax=0, shading="nearest", rasterized=True)
    ax.invert_yaxis()
    ax.set(title=t, xlabel="z [mm]", ylabel="φ [deg]")
fig.colorbar(im, ax=axs[:2], label="dB", shrink=0.85)
dd = np.abs(fn - on)
im = axs[2].pcolormesh(z, phi, np.log10(np.maximum(dd, 1e-6)), cmap=ERR, vmin=-4, vmax=-0.5, shading="nearest", rasterized=True)
axs[2].invert_yaxis()
axs[2].set(title=f"|difference| (each normalised to its peak)\nrel. L2 = {np.linalg.norm(fn - on) / np.linalg.norm(on):.1%}", xlabel="z [mm]")
fig.colorbar(im, ax=axs[2], label="log₁₀ of peak", shrink=0.85)
fig.suptitle("Objects in a pipe (thesis Fig. 7.13): amplitude C-scan, original vs corrected CPSM")
save(fig, "pipe_original_vs_corrected.png")
pipe_rel = np.linalg.norm(fn - on) / np.linalg.norm(on)

# ---------- 2. wire: original PSM vs the exact integral ----------
o, s = C("wire_orig"), C("wire_sq")
zw = L("wire_orig_z") * 1e3
xw = np.arange(o.shape[1]) * 0.5
e = np.abs(np.abs(o) - np.abs(s)) / np.abs(s).max()
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
for ax, a, t in [(axs[0], np.abs(o), "Original psm.m (via legacy::psm)"), (axs[1], np.abs(s), "v2 SQ-PSM (the exact ω-integral)")]:
    im = mesh(ax, db(a), xw, zw, cmap=GRAY, vmin=-40, vmax=0)
    ax.set(title=t, xlabel="x [mm]", ylabel="z [mm]")
fig.colorbar(im, ax=axs[:2], label="dB", shrink=0.85)
im = mesh(axs[2], np.log10(np.maximum(e, 1e-8)), xw, zw, cmap=ERR, vmin=-6, vmax=-3)
axs[2].set(title=f"|envelope difference| / peak (max {e.max():.1e})", xlabel="x [mm]")
fig.colorbar(im, ax=axs[2], label="log₁₀ of peak", shrink=0.85)
fig.suptitle("Wire targets: spot the difference (you can't), then look at the right panel")
save(fig, "wire_original_vs_sq.png")
rows = e.max(axis=1)
wire_top, wire_bot = float(rows[: len(rows) // 10].mean()), float(rows[-len(rows) // 10:].mean())

# ---------- 3. pins: original MULOK vs default MULOK vs SQ-PSM ----------
def env_on_uniform(name, nrows=260):
    a = np.abs(C(name))
    zz = L(name + "_z")
    k = monotone(zz)
    a, zz = a[k], zz[k]
    zt = np.linspace(zz[0], zz[-1], nrows)
    idx = np.clip(np.searchsorted(zz, zt), 1, len(zz) - 1)
    idx -= (zt - zz[idx - 1]) < (zz[idx] - zt)
    return a[idx], zt * 1e3


def env_err(name, ref):
    m, p = np.abs(C(name)), np.abs(C(ref))
    alpha = (m * p).sum() / (p * p).sum()
    return np.linalg.norm(m - alpha * p) / np.linalg.norm(alpha * p), np.abs(m - alpha * p) / (alpha * p).max(), L(name + "_z")


fig, axs = plt.subplots(1, 5, figsize=(20, 6), constrained_layout=True)
xp = np.arange(np.abs(C("pins_sq")).shape[1]) * 1.0
for ax, n, t in [(axs[0], "pins_origmulok", "Original mulok.m"), (axs[1], "pins_mulok", "v1/v2 default MULOK"), (axs[2], "pins_sq", "v2 SQ-PSM")]:
    a, zt = env_on_uniform(n)
    im = mesh(ax, db(a), xp, zt, cmap=GRAY, vmin=-40, vmax=0)
    ax.set(title=t, xlabel="x [mm]", ylabel="z [mm]")
fig.colorbar(im, ax=axs[:3], label="dB", shrink=0.8, location="bottom")
for ax, n, ref, t in [(axs[3], "pins_origmulok", "pins_sq_bw", "Original MULOK: error"), (axs[4], "pins_mulok", "pins_sq", "Default MULOK: error")]:
    E, err, zz = env_err(n, ref)
    k = monotone(zz)
    im2 = mesh(ax, np.log10(np.maximum(err[k], 1e-8)), xp, zz[k] * 1e3, cmap=ERR, vmin=-5, vmax=-1)
    ax.set(title=f"{t}\nE = {E:.1e} vs SQ-PSM", xlabel="x [mm]")
fig.colorbar(im2, ax=axs[3:], label="log₁₀ |envelope error| / peak", shrink=0.8, location="bottom")
fig.suptitle("Pins in water, plexiglass and aluminium: the three MULOK bugs, before and after")
save(fig, "pins_mulok_before_after.png")
E_orig = env_err("pins_origmulok", "pins_sq_bw")[0]
E_def = env_err("pins_mulok", "pins_sq")[0]

(OUT / "numbers.txt").write_text(
    f"pipe_cscan_rel_l2 {pipe_rel:.4e}\nwire_err_max {e.max():.4e}\nwire_err_top10pct {wire_top:.4e}\nwire_err_bottom10pct {wire_bot:.4e}\n"
    f"pins_E_original_mulok {E_orig:.4e}\npins_E_default_mulok {E_def:.4e}\n")
print((OUT / "numbers.txt").read_text())
