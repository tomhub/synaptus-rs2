"""Independent oracle for N6 (TEST_PLAN_V2.md, amendment A1).

Evaluates, at selected pixels, the target of SQ-PSM as *defined* in
research/synaptus-v2-improvements-2026-11-23.md §5.2:

  pixel(z, x) = (1/nfx) Σ_k e^{i k_x x} J_l(z; k_x),
  J_l = (1/Δω₁) ∫_{max(ω_lo, ω*_l)}^{ω_hi} W(ω,k) e^{-iωt_d} e^{i[Σ_{j<l} k_{z,j} d_j + k_{z,l}(z − z_top,l)]} dω

with W the DTFT of the samples after a spatial FFT (NumPy), ω*_l = max_{j≤l} c_j|k|/2,
k_{z,j} = (2/c_j)√((ω−ω_cj)(ω+ω_cj)), Δω₁ = 2π f_s / 2^⌈log₂ n_t⌉.

Independent of the Rust method: no Gauss–Legendre, no s-substitution, no NUFFT.
Quadrature: composite tanh–sinh in ω on sub-intervals of width Δω₁/2, levels
h = 1/8 and 1/16. Shares with the Rust code only the problem definition and
the data files.

Usage:  python3 -I verify/sq_oracle.py verify/n6/<case> [...]
"""
import json
import sys
from pathlib import Path

import numpy as np

T_MAX = 3.5


def ts_rule(h):
    """tanh–sinh on [0,1]: fraction from lower end, fraction from upper end, weight."""
    t = np.arange(-T_MAX, T_MAX + h / 2, h)
    u = 0.5 * np.pi * np.sinh(t)
    frac_lo = 1.0 / (1.0 + np.exp(-2.0 * u))   # (1+x)/2
    frac_hi = 1.0 / (1.0 + np.exp(2.0 * u))    # (1-x)/2
    w = h * 0.5 * np.pi * np.cosh(t) / np.cosh(u) ** 2 / 2.0  # weight on [0,1]
    return frac_lo, frac_hi, w


def oracle(case_dir, h):
    case = json.loads((case_dir / "case.json").read_text())
    nt, nx, nfx = case["nt"], case["nx"], case["nfx"]
    fs, td, dx = case["fs"], case["t_delay"], case["x_step"]
    wl, wh = 2 * np.pi * case["f_low"], 2 * np.pi * case["f_high"]
    c, d = np.array(case["velocities"]), np.array(case["thicknesses"])
    p = np.fromfile(case_dir / "data.f64", dtype="<f8").reshape((nx, nt)).T  # (nt, nx), t fastest
    P = np.fft.fft(p, n=nfx, axis=1)                                          # (nt, nfx)
    kx = 2 * np.pi * np.fft.fftfreq(nfx, d=dx)
    dw1 = 2 * np.pi * fs / (1 << int(np.ceil(np.log2(nt))))
    n_idx = np.arange(nt) / fs
    fl, fh, fw = ts_rule(h)

    # global sub-interval edges
    edges = np.arange(wl, wh, dw1 / 2)
    edges = np.append(edges, wh)
    lo, hi = edges[:-1], edges[1:]
    om = (lo[:, None] + (hi - lo)[:, None] * fl[None, :])                     # nodes (nsub, nts)
    # use upper-end form where closer to the upper end
    om_hi = (hi[:, None] - (hi - lo)[:, None] * fh[None, :])
    om = np.where(fl[None, :] > 0.5, om_hi, om)
    wq = (hi - lo)[:, None] * fw[None, :]
    om_flat, wq_flat = om.ravel(), wq.ravel()
    sub_of = np.repeat(np.arange(len(lo)), len(fl))

    # W at all global nodes for all k: W = E @ P, E = exp(-i ω n/fs)
    W = np.empty((om_flat.size, nfx), dtype=complex)
    blk = 8192
    for s in range(0, om_flat.size, blk):
        E = np.exp(-1j * np.outer(om_flat[s:s + blk], n_idx))
        W[s:s + blk] = E @ P

    out = {}
    pixels = case["pixels"]
    nl = len(c)
    for l in range(nl):
        pl = [q for q in pixels if q["layer"] == l]
        zs = np.array([q["z"] for q in pl])
        xs = np.array([q["x"] for q in pl])
        ztop = 0.0 if l == 0 else float(np.sum(d[:l]))
        J = np.zeros((len(pl), nfx), dtype=complex)
        for k in range(nfx):
            wc = c[: l + 1] * abs(kx[k]) / 2.0
            wstar = wc.max()
            a = max(wl, wstar)
            if a >= wh:
                continue
            # full global sub-intervals inside [a, wh]
            full = lo >= a
            sel = full[sub_of]
            omg, wgt, Wk = om_flat[sel], wq_flat[sel], W[sel, k]
            off = omg - a  # ≥ 0
            # partial first sub-interval [a, e*] (if a is not on an edge)
            e_star = lo[full][0] if full.any() else wh
            if e_star > a:
                L = e_star - a
                off_p = np.where(fl > 0.5, L - L * fh, L * fl)
                om_p = a + off_p
                E = np.exp(-1j * np.outer(om_p, n_idx))
                omg = np.concatenate([om_p, omg])
                off = np.concatenate([off_p, off])
                wgt = np.concatenate([L * fw, wgt])
                Wk = np.concatenate([E @ P[:, k], Wk])
            # vertical wavenumbers; ω − ω_cj = (ω − a) + (a − ω_cj) avoids cancellation at the cut-off
            kz = [(2.0 / c[j]) * np.sqrt((off + (a - wc[j])) * (omg + wc[j])) for j in range(l + 1)]
            phi0 = -omg * td + sum(kz[j] * d[j] for j in range(l))
            base = wgt * Wk / dw1
            for ip, z in enumerate(zs):
                J[ip, k] = np.sum(base * np.exp(1j * (phi0 + kz[l] * (z - ztop))))
        vals = (J * np.exp(1j * np.outer(xs, kx))).sum(axis=1) / nfx
        for q, v in zip(pl, vals):
            out[(q["layer"], q["iz"], q["ix"])] = v
    return case, out


def main():
    all_ok = True
    args = sys.argv[1:]
    fine, coarse, tag = (1 / 16, 1 / 8, "") 
    if args and args[0] == "--n6b":  # follow-up N6b: levels 1/32 (value) vs 1/16 (check)
        fine, coarse, tag, args = 1 / 32, 1 / 16, "_n6b", args[1:]
    for arg in args:
        case_dir = Path(arg)
        case, o16 = oracle(case_dir, fine)
        _, o8 = oracle(case_dir, coarse)
        rows = []
        for l in sorted({q["layer"] for q in case["pixels"]}):
            pl = [q for q in case["pixels"] if q["layer"] == l]
            omax = max(abs(o16[(q["layer"], q["iz"], q["ix"])]) for q in pl)
            for q in pl:
                key = (q["layer"], q["iz"], q["ix"])
                sq = complex(q["sq_re"], q["sq_im"])
                conv = abs(o16[key] - o8[key]) / omax
                err = abs(sq - o16[key]) / q["layer_max"]
                if conv > 1e-12:
                    verdict = "INCONCLUSIVE"
                elif err <= 1e-9:
                    verdict = "PASS"
                else:
                    verdict = "FAIL"
                all_ok &= verdict == "PASS"
                rows.append(dict(layer=l + 1, iz=q["iz"], ix=q["ix"], oracle_abs=abs(o16[key]), sq_abs=abs(sq),
                                 err_rel_layer_max=err, oracle_conv=conv, verdict=verdict))
                print(f"{case_dir.name:5} L{l + 1} (iz={q['iz']:4}, ix={q['ix']:3}): |oracle|={abs(o16[key]):.6e} "
                      f"|SQ-oracle|/max|SQ|={err:.2e}  oracle conv={conv:.1e}  {verdict}")
        (case_dir / f"oracle_result{tag}.json").write_text(json.dumps(rows, indent=1))
    print(f"N6{tag.replace('_n6b','b')} overall:", "PASS" if all_ok else "NOT ALL PASS")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
