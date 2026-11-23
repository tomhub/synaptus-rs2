"""T5: high-precision reference for the PSM migration step (eq. 3.21).

Independent of the Rust implementation: it takes only the double-precision
INPUTS that Rust migrates (omega, kx, W, c, dz, z_offset) and evaluates

    i(x_m, z_j) = (1/nfx) * sum_k e^{2 pi i k m / nfx} * sum_w W[k,w] * exp(i kz(w,k) z_j)

with kz = sqrt((2 w / c)^2 - kx^2) (0 and masked when evanescent), and
z_j = z_offset + j*dz, in 40-digit arithmetic (mpmath 1.3.0).
Usage: python -I verify/psm_rounding_ref.py <dir written by examples/psm_rounding.rs>
"""
import struct
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40


def rd(p):
    b = Path(p).read_bytes()
    return list(struct.unpack("<%dd" % (len(b) // 8), b))


d = Path(sys.argv[1])
meta = {k: float(v) for k, v in (l.split() for l in (d / "meta.txt").read_text().split("\n") if l)}
omega, kx = rd(d / "omega.f64"), rd(d / "kx.f64")
nw, nfx, nx = len(omega), len(kx), int(meta["nx"])
wr = rd(d / "W.f64")
W = [[mp.mpc(wr[2 * (k * nw + i)], wr[2 * (k * nw + i) + 1]) for i in range(nw)] for k in range(nfx)]
c, dz, z0 = mp.mpf(meta["c"]), mp.mpf(meta["dz"]), mp.mpf(meta["z_offset"])
kz = [[None] * nw for _ in range(nfx)]
n_evan = 0
for k in range(nfx):
    for i in range(nw):
        a = 2 * mp.mpf(omega[i]) / c
        kz2 = a * a - mp.mpf(kx[k]) ** 2
        if kz2 >= 0:
            kz[k][i] = mp.sqrt(kz2)
        else:
            n_evan += 1
variants = [l for l in (d / "variants.txt").read_text().split("\n") if l]
nz = int(variants[0].split("nz=")[1].split()[0])
lines = sorted(set([0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, nz - 1]) & set(range(nz)))
twiddle = [[mp.expjpi(mp.mpf(2 * k * m) / nfx) for m in range(nx)] for k in range(nfx)]

ref = {}
for j in lines:
    z = z0 + j * dz
    S = []
    for k in range(nfx):
        s = mp.mpc(0)
        for i in range(nw):
            if kz[k][i] is not None:
                s += W[k][i] * mp.expj(kz[k][i] * z)
        S.append(s)
    ref[j] = [sum(S[k] * twiddle[k][m] for k in range(nfx)) / nfx for m in range(nx)]
gmax = max(abs(v) for j in lines for v in ref[j])

print(f"# evanescent elements: {n_evan}; image planes: {nz}; lines checked: {lines}")
print(f"# error metric: max_m |img(m,j) - ref(m,j)| / max|ref| (over checked lines); unit roundoff u = 1.11e-16")
hdr = "variant | " + " | ".join(f"j={j}" for j in lines) + " | time"
print(hdr)
print("|".join(["---"] * (len(lines) + 2)))
for v in variants:
    idx, rest = v.split(" ", 1)
    img = rd(d / f"img{idx}.f64")
    errs = []
    for j in lines:
        e = max(abs(mp.mpc(img[2 * (j * nx + m)], img[2 * (j * nx + m) + 1]) - ref[j][m]) for m in range(nx))
        errs.append(float(e / gmax))
    name = rest.split(" nz=")[0]
    t = rest.split("time=")[1]
    print(name + " | " + " | ".join(f"{e:.1e}" for e in errs) + f" | {t}s")
