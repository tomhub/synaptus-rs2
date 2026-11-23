"""Build figures/v2_comparison.html from recorded runs (no recomputation of methods).

Inputs (all OBSERVED, see ACCURACY_V2.md):
  runs/20261123-161018-v2-n7-s1/{s1.tsv,n7.tsv,fig/}     N7, S1, envelope images
  runs/20261123-160703-v2-tests-A2/{test_output.txt,n6b_A2_output.txt}  N1, N4, N6b
  ../synaptus-rs/runs/20261123-152041-tests-with-cpsm/test_output.txt   v1 T1
Run from the crate root:  python3 -I figures/build_v2_page.py
"""
import base64
import io
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "runs/20261123-161018-v2-n7-s1"
TESTS = ROOT / "runs/20261123-160703-v2-tests-A2"
V1RUN = ROOT.parent / "synaptus-rs/runs/20261123-152041-tests-with-cpsm"
ERR = LinearSegmentedColormap.from_list("err", ["#0a1020", "#26245e", "#6b2a8a", "#c2477f", "#f08a5d", "#ffe3a3"])
ERR_LO, ERR_HI = -10.0, -1.0  # log10 colour range


def tsv(path):
    lines = path.read_text().strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def floats(pattern, text):
    return [float(x) for x in re.findall(pattern, text)]


def load(name):
    dims = [int(x) for x in (RUN / "fig" / f"{name}.txt").read_text().split()[1:]]
    return np.fromfile(RUN / "fig" / f"{name}.f64", "<f8").reshape(dims, order="F")


def uniform_rows(img, z, nrows):
    """Nearest-sample resampling onto a uniform depth axis (layers have different Δz);
    overlapping samples at layer joins are dropped so z is increasing."""
    keep = np.ones(len(z), bool)
    for i in range(1, len(z)):
        if z[i] <= z[:i][keep[:i]].max():
            keep[i] = False
    img, z = img[keep], z[keep]
    zt = np.linspace(z[0], z[-1], nrows)
    idx = np.clip(np.searchsorted(z, zt), 1, len(z) - 1)
    idx -= (zt - z[idx - 1]) < (z[idx] - zt)
    return img[idx], float(z[0]), float(z[-1])


def png(a, scale):
    a = np.clip((a - ERR_LO) / (ERR_HI - ERR_LO), 0, 1)
    a = np.repeat(np.repeat(a, scale, 0), scale, 1)
    buf = io.BytesIO()
    plt.imsave(buf, a, cmap=ERR, vmin=0, vmax=1, format="png")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def err_map(ds, meth, ref, nrows, scale):
    m, p, z = load(f"{ds}_{meth}"), load(f"{ds}_{ref}"), load(f"{ds}_{meth}_z")
    alpha = (m * p).sum() / (p * p).sum()
    e = np.abs(m - alpha * p) / (alpha * p).max()
    E = np.linalg.norm(m - alpha * p) / np.linalg.norm(alpha * p)
    img, z0, z1 = uniform_rows(e, z, nrows)
    return png(np.log10(np.maximum(img, 1e-300)), scale), float(E), [z0 * 1e3, z1 * 1e3]


data = {}
# ---- S1 / N7 ----
data["s1"] = [dict(dataset=r["dataset"], method=r["method"], E=float(r["E_vs_ref"]),
                   t=float(r["t_median_s"]), tmin=float(r["t_min_s"]), tmax=float(r["t_max_s"])) for r in tsv(RUN / "s1.tsv")]
data["n7"] = [dict(dataset=r["dataset"], m=int(r["m"]), E=float(r["E_vs_SQ"])) for r in tsv(RUN / "n7.tsv")]

# ---- headline numbers from test logs ----
t = (TESTS / "test_output.txt").read_text()
n1 = floats(r"N[12] M=\s*\d+ N=\s*\d+ eps=1e-12: max\|err\|/\|\|[cf]\|\|1 = ([0-9.e+-]+)", t)
n4 = floats(r"N4 \w+\s+m=\d layer \d: rel L2 = ([0-9.e+-]+)", t)
n6b = (TESTS / "n6b_A2_output.txt").read_text()
n6b_err = floats(r"\|SQ-oracle\|/max\|SQ\|=([0-9.e+-]+)", n6b)
n6b_pass = len(re.findall(r"^\S+\s+L\d.*PASS$", n6b, re.M))
v1 = (V1RUN / "test_output.txt").read_text()
t1 = floats(r"rel\.err=([0-9.e+-]+)", v1)
sq_self = max(r["E"] for r in data["s1"] if r["method"].startswith("v2 SQ"))
by = lambda ds, prefix: next(r for r in data["s1"] if r["dataset"] == ds and r["method"].startswith(prefix))
rng = lambda prefix: [min(by(d, prefix)["E"] for d in ("wire", "pins", "plane")), max(by(d, prefix)["E"] for d in ("wire", "pins", "plane"))]
data["facts"] = dict(
    nufft_worst=max(n1), n4_worst=max(n4), n6b_worst=max(n6b_err), n6b_pass=n6b_pass, n6b_total=len(n6b_err),
    t1_worst=max(t1[:14]), sq_self=sq_self,
    v1_ref=rng("v1 PSM (m=8)"), acc_mulok=rng("v1 accurate MULOK"), v1_psm=rng("v1 PSM (m=1)"), leg_mulok=rng("v1 legacy MULOK"),
)

# ---- error maps ----
maps = {}
for ds, nrows, scale in (("pins", 200, 3), ("wire", 340, 1)):
    maps[ds] = []
    for meth, ref, label, who in (("legmulok", "refbw", "Original MULOK (legacy port)", "v1"),
                                  ("v1psm", "ref", "v1 PSM, default", "v1"),
                                  ("accmulok", "ref", "v1 accurate MULOK", "v1"),
                                  ("nufft_m32", "ref", "v2 NUFFT-PSM, 32× padding", "v2"),
                                  ("sq", "ref", "v2 SQ-PSM (graded by itself)", "v2")):
        img, E, zr = err_map(ds, meth, ref, nrows, scale)
        maps[ds].append(dict(img=img, E=E, label=label, who=who, z=zr))
data["maps"] = maps
data["err_range"] = [ERR_LO, ERR_HI]
data["err_ramp"] = "linear-gradient(90deg," + ",".join(
    matplotlib.colors.to_hex(ERR(x)) + f" {x * 100:.0f}%" for x in np.linspace(0, 1, 9)) + ")"

# ---- gallery (figures/gallery, from examples/gallery_v2.rs + figures/plot_gallery.py) ----
GAL = ROOT / "figures/gallery"
gnum = dict(l.split() for l in (GAL / "numbers.txt").read_text().splitlines() if l.strip())
patt = (GAL / "pipe_attribution.txt").read_text()
pick = lambda key: float(re.search(re.escape(key) + r"\s*([0-9.e+-]+)", patt).group(1))
b64 = lambda f: "data:image/png;base64," + base64.b64encode((GAL / f).read_bytes()).decode()
data["gallery"] = dict(
    pins_png=b64("pins_mulok_before_after.png"), pipe_png=b64("pipe_original_vs_corrected.png"),
    pins_orig=float(gnum["pins_E_original_mulok"]), pins_def=float(gnum["pins_E_default_mulok"]),
    pipe_total=pick("total:"), pipe_fixes=pick("fixes only (orig vs corrected-Haun):"),
    pipe_transfer=pick("transfer only (corrected-Haun vs corrected-Hankel):"),
    pipe_moved=pick("topography: pixels above -15 dB with different range sample:"),
)

page = (ROOT / "figures/v2_comparison.template.html").read_text()
page = page.replace("/*__DATA__*/", "const D = " + json.dumps(data) + ";")
(ROOT / "figures/v2_comparison.html").write_text(page)
print(json.dumps(data["facts"], indent=1))
print("wrote figures/v2_comparison.html", f"{len(page) / 1e6:.2f} MB")
