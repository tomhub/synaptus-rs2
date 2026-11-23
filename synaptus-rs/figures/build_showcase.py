"""Build figures/showcase.html: original Synaptus (Octave) vs the Rust port.

Inputs: Octave reference outputs in data/ (octave/export_*.m), Rust outputs in
figures/data/showcase/ (cargo run --release --example showcase_export) and
figures/data/ (cargo run --release --example figures). Every number on the page
is computed here; the template only lays it out. Run from the crate root:

    python3 -I figures/build_showcase.py
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
REF, RS, FIG = ROOT / "data", ROOT / "figures/data/showcase", ROOT / "figures/data"
ECHO = LinearSegmentedColormap.from_list("echo", ["#0b1424", "#173f5f", "#20788a", "#3caea3", "#f6d55c", "#fffbe8"])
DIFF = LinearSegmentedColormap.from_list("diff", ["#0b1424", "#3b1f6b", "#a2338a", "#ef6f6c", "#ffd6c9"])


def rc(path):
    dims = [int(x) for x in path.with_suffix(".txt").read_text().split("\n")[0].split()[1:]]
    v = np.fromfile(path.with_suffix(".c64"), "<f8")
    return (v[0::2] + 1j * v[1::2]).reshape(dims + [1] * (3 - len(dims)), order="F")


def zgrid(prefix, l):
    for line in (REF / f"{prefix}_layer{l}.txt").read_text().splitlines():
        if line.startswith("z "):
            return np.array([float(x) for x in line.split()[1:]])


def f64(name):
    dims = [int(x) for x in (FIG / f"{name}.txt").read_text().split()[1:]]
    return np.fromfile(FIG / f"{name}.f64", "<f8").reshape(dims, order="F")


def stack(layers, zs, nrows):
    """Place per-layer (z, x) images on one uniform depth axis (nearest sample),
    so each layer occupies its true physical thickness despite different Δz."""
    for k in range(len(zs) - 1):  # drop overlap past the next layer's first depth
        keep = zs[k] < zs[k + 1][0]
        layers[k], zs[k] = layers[k][keep], zs[k][keep]
    img, z = np.concatenate(layers), np.concatenate(zs)
    zt = np.linspace(z[0], z[-1], nrows)
    idx = np.clip(np.searchsorted(z, zt), 1, len(z) - 1)
    idx -= (zt - z[idx - 1]) < (z[idx] - zt)
    return img[idx], (z[0], z[-1])


def png(a, cmap, vmin, vmax, scale=1):
    a = np.clip((a - vmin) / (vmax - vmin), 0, 1)
    if scale > 1:
        a = np.repeat(np.repeat(a, scale, 0), scale, 1)
    buf = io.BytesIO()
    plt.imsave(buf, a, cmap=cmap, vmin=0, vmax=1, format="png")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def ramp(cmap):
    return "linear-gradient(90deg," + ",".join(
        matplotlib.colors.to_hex(cmap(t)) + f" {t * 100:.0f}%" for t in np.linspace(0, 1, 9)) + ")"


def db(a, ref):
    return 20 * np.log10(np.maximum(a, 1e-300) / ref)


DIFF_MAX = 8e-16
stats = {}


def pair(key, name, alg, layers_n, project, rng, title, scale=1, nrows=None):
    os_, rs_, zs = [], [], []
    rels, mxs, n_eq, n = [], [], 0, 0
    for l in range(1, layers_n + 1):
        o, r = rc(REF / f"{name}_{alg}_layer{l}"), rc(RS / f"{name}_{alg}_layer{l}")
        rels.append(np.linalg.norm(r - o) ** 2), mxs.append(np.linalg.norm(o) ** 2)
        n_eq += int(np.sum(r == o)); n += o.size
        os_.append(o), rs_.append(r), zs.append(zgrid(f"{name}_{alg}", l))
    peak = max(np.abs(o).max() for o in os_)
    proj = lambda arrs: [project(a) for a in arrs]
    ao, zr = stack(proj([np.abs(o) for o in os_]), list(zs), nrows)
    ar, _ = stack(proj([np.abs(r) for r in rs_]), list(zs), nrows)
    dd, _ = stack(proj([np.abs(r - o) / peak for o, r in zip(os_, rs_)]), list(zs), nrows)
    stats[key] = dict(title=title, rel=float(np.sqrt(sum(rels) / sum(mxs))),
                      worst=float(max(np.abs(r - o).max() for o, r in zip(os_, rs_)) / peak),
                      identical=n_eq / n, pixels=n, z=[zr[0] * 1e3, zr[1] * 1e3])
    return dict(orig=png(db(ao, ao.max()), ECHO, -rng, 0, scale), rust=png(db(ar, ar.max()), ECHO, -rng, 0, scale),
                diff=png(dd, DIFF, 0, DIFF_MAX, scale), rng=rng)


imgs = {}
yproj = lambda a: a.max(axis=2)  # (z, x, y) -> (z, x), max over y
imgs["wire"] = pair("wire", "wire", "psm", 1, yproj, 40, "Wire targets in water, PSM", nrows=682)
imgs["pins"] = pair("pins", "pins", "mulok", 3, yproj, 40, "Pins in plexiglass and aluminium, MULOK", scale=3, nrows=230)
imgs["plane"] = pair("plane", "plane", "psm", 3, yproj, 40, "Flat-bottom holes, 3D PSM (brightest along y)", scale=3, nrows=170)

# CPSM: C-scan (max over r) of the cylindrical image
o, r = rc(REF / "pipe_cpsm_haun"), rc(RS / "pipe_cpsm_haun")
peak = np.abs(o).max()
co, cr = np.abs(o).max(axis=0), np.abs(r).max(axis=0)
stats["pipe"] = dict(title="Objects inside a pipe, cylindrical PSM", rel=float(np.linalg.norm(r - o) / np.linalg.norm(o)),
                     worst=float(np.abs(r - o).max() / peak), identical=float(np.mean(r == o)), pixels=int(o.size))
imgs["pipe"] = dict(orig=png(db(co, co.max()), ECHO, -25, 0, 3), rust=png(db(cr, cr.max()), ECHO, -25, 0, 3),
                    diff=png((np.abs(r - o) / peak).max(axis=0), DIFF, 0, DIFF_MAX, 3), rng=25)

# Beyond the original: envelope error maps on the pins dataset (E-metric of ACCURACY.md)
def env_layers(prefix):
    return [f64(f"{prefix}_L{l}")[:, :, 0] for l in (1, 2, 3)], [f64(f"{prefix}_L{l}_z") for l in (1, 2, 3)]


def err_map(m_prefix, ref_prefix):
    m, zm = env_layers(m_prefix)
    p, _ = env_layers(ref_prefix)
    M, P = np.concatenate(m), np.concatenate(p)
    alpha = np.vdot(P, M).real / np.vdot(P, P).real
    E = np.linalg.norm(M - alpha * P) / np.linalg.norm(alpha * P)
    pk = (alpha * P).max()
    e, _ = stack([np.abs(a - alpha * b) / pk for a, b in zip(m, p)], zm, 120)
    return e, float(E)


ERR_MAX = 0.01
e_leg, E_leg = err_map("pins_legmulok", "pins_ref8bw")
e_acc, E_acc = err_map("pins_accmulok", "pins_ref8")
imgs["better"] = dict(legacy=png(e_leg, DIFF, 0, ERR_MAX, 3), accurate=png(e_acc, DIFF, 0, ERR_MAX, 3))
stats["better"] = dict(legacy=E_leg, accurate=E_acc, ratio=E_leg / E_acc)
stats["total_pixels"] = sum(stats[k]["pixels"] for k in ("pipe",)) + sum(
    rc(REF / p.stem).size for p in REF.glob("*_layer*.c64"))
stats["worst_all"] = max(stats[k]["worst"] for k in ("wire", "pins", "plane", "pipe"))
# worst over *all* T1 layers, not only the four shown
for p in REF.glob("*_layer*.c64"):
    o, r = rc(p.with_suffix("")), rc(RS / p.stem)
    stats["worst_all"] = max(stats["worst_all"], float(np.abs(r - o).max() / np.abs(o).max()))

page = (ROOT / "figures/showcase.template.html").read_text()
page = page.replace("/*__DATA__*/", "const IMGS = " + json.dumps(imgs) + ";\nconst STATS = " + json.dumps(stats) + ";\n"
                    + f"const RAMP_ECHO = {json.dumps(ramp(ECHO))}; const RAMP_DIFF = {json.dumps(ramp(DIFF))};\n"
                    + f"const DIFF_MAX = {DIFF_MAX}; const ERR_MAX = {ERR_MAX};")
(ROOT / "figures/showcase.html").write_text(page)
print(json.dumps(stats, indent=1))
print("wrote", ROOT / "figures/showcase.html", f"{len(page) / 1e6:.2f} MB")
