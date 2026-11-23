"""Build figures/showcase_v2.html: the original code's images vs v2's, drag to compare.

Inputs: figures/data/gallery/ from `cargo run --release --example gallery_v2`.
"Original" images come from the legacy ports, which match the original Octave
code to ~1e-15 (tests T1/T7). Every number on the page is computed here.
Run from the crate root:  python3 -I figures/build_showcase_v2.py
"""
import base64
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "figures/data/gallery"
ECHO = LinearSegmentedColormap.from_list("echo", ["#0b1424", "#173f5f", "#20788a", "#3caea3", "#f6d55c", "#fffbe8"])
DIFF = LinearSegmentedColormap.from_list("diff", ["#0a1020", "#26245e", "#6b2a8a", "#c2477f", "#f08a5d", "#ffe3a3"])
DLO, DHI = -5.0, -1.0


def L(n):
    dims = [int(x) for x in (D / f"{n}.txt").read_text().split()[1:]]
    return np.fromfile(D / f"{n}.f64", "<f8").reshape(dims, order="F")


def env(n):
    return np.abs(L(n + "_re") + 1j * L(n + "_im"))


def on_grid(a, z, zt):
    """Nearest-sample resampling of (z, x) onto target depths zt; joins de-duplicated."""
    keep = [0]
    for i in range(1, len(z)):
        if z[i] > z[keep[-1]]:
            keep.append(i)
    a, z = a[keep], z[keep]
    idx = np.clip(np.searchsorted(z, zt), 1, len(z) - 1)
    idx -= (zt - z[idx - 1]) < (z[idx] - zt)
    return a[idx]


def pair_on_common_grid(a, za, b, zb, nrows):
    zt = np.linspace(max(za[0], zb[0]), min(za[-1], zb[-1]), nrows)
    return on_grid(a, za, zt), on_grid(b, zb, zt), (float(zt[0]) * 1e3, float(zt[-1]) * 1e3)


def png(a, cmap, lo, hi, scale):
    a = np.clip((a - lo) / (hi - lo), 0, 1)
    if scale > 1:
        a = np.repeat(np.repeat(a, scale, 0), scale, 1)
    buf = io.BytesIO()
    plt.imsave(buf, a, cmap=cmap, vmin=0, vmax=1, format="png")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def db(a):
    return 20 * np.log10(np.maximum(a, 1e-300) / a.max())


def e_metric(m, p):
    """v1 ACCURACY E-metric on envelopes: ‖m − αp‖/‖αp‖, α least squares."""
    alpha = (m * p).sum() / (p * p).sum()
    return float(np.linalg.norm(m - alpha * p) / np.linalg.norm(alpha * p))


def err_map(m, p):
    """|m − αp| / max(αp) on the method's own grid (α least squares)."""
    alpha = (m * p).sum() / (p * p).sum()
    return np.abs(m - alpha * p) / (alpha * p).max()


def lg(a):
    return np.log10(np.maximum(a, 1e-12))


def pack(o, v, eo, ev, rng, scale):
    return dict(orig=png(db(o), ECHO, -rng, 0, scale), v2=png(db(v), ECHO, -rng, 0, scale),
                err_o=png(lg(eo), DIFF, DLO, DHI, scale), err_v=png(lg(ev), DIFF, DLO, DHI, scale), rng=rng)


tabs = {}
# wire: original psm.m vs SQ-PSM, which is itself the reference (exact ω-integral)
o, s = env("wire_orig"), env("wire_sq")
t = pack(o, s, err_map(o, s), np.zeros_like(s), 40, 1)
t.update(e_orig=e_metric(o, s), e_v2=None)
tabs["wire"] = t

# pins: original mulok.m vs default MULOK, each against SQ-PSM on its own depth grid
oo, zo = env("pins_origmulok"), L("pins_origmulok_z")
vv, zv = env("pins_mulok"), L("pins_mulok_z")
eo, ev = err_map(oo, env("pins_sq_bw")), err_map(vv, env("pins_sq"))
og, vg, zr = pair_on_common_grid(oo, zo, vv, zv, 230)
eog, evg, _ = pair_on_common_grid(eo, zo, ev, zv, 230)
t = pack(og, vg, eog, evg, 40, 3)
t.update(e_orig=e_metric(oo, env("pins_sq_bw")), e_v2=e_metric(vv, env("pins_sq")), z=zr)
tabs["pins"] = t

# plane: original psm.m vs SQ-PSM (the reference), brightest value over y
po, zpo = L("plane_orig_envmaxy"), L("plane_orig_z")
ps, zps = L("plane_sq_envmaxy"), L("plane_sq_z")
og, sg, zr = pair_on_common_grid(po, zpo, ps, zps, 170)
t = pack(og, sg, err_map(og, sg), np.zeros_like(sg), 40, 3)
t.update(e_orig=e_metric(po, ps), e_v2=None, z=zr)
tabs["plane"] = t

# pipe: original cpsm.m vs corrected CPSM; no exact reference, so both sides show the difference
co, cf = L("pipe_orig_cscan"), L("pipe_fixed_cscan")
dd = np.abs(cf / cf.max() - co / co.max())
t = pack(co, cf, dd, dd, 25, 3)
t.update(rel=float(np.linalg.norm(cf / cf.max() - co / co.max()) / np.linalg.norm(co / co.max())))
tabs["pipe"] = t

ramp = lambda cm: "linear-gradient(90deg," + ",".join(
    matplotlib.colors.to_hex(cm(x)) + f" {x * 100:.0f}%" for x in np.linspace(0, 1, 9)) + ")"
page = (ROOT / "figures/showcase_v2.template.html").read_text()
page = page.replace("/*__DATA__*/", "const T = " + json.dumps(tabs) + ";\n"
                    + f"const RAMP_ECHO = {json.dumps(ramp(ECHO))}; const RAMP_DIFF = {json.dumps(ramp(DIFF))}; const DLO = {DLO}, DHI = {DHI};")
(ROOT / "figures/showcase_v2.html").write_text(page)
for k, v in tabs.items():
    print(k, {a: b for a, b in v.items() if a not in ("orig", "v2", "err_o", "err_v")})
print("wrote figures/showcase_v2.html", f"{len(page) / 1e6:.2f} MB")
