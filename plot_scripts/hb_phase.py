"""Condon-domain coexistence (B1, B2) on every dHvA loop, from `nsrs study hb`.

At fixed mu_n the field is in equilibrium where the potential
    g(B) = 4 pi eps_vac(B) - 4 pi P_m(mu_n, B)        [G^2, Gaussian]
is convex, since H = dg/dB (H = B - 4 pi M for Maxwell). Where g is not convex,
a uniform B is replaced by coexisting domains B1 < B2 given by the common
tangent of g, i.e. the 2x2 system
    H(B1) = H(B2)                                     (same H in both phases)
    g(B2) - g(B1) = H(B1) (B2 - B1)                   (equal areas)
with 4 pi eps_vac = B^2/2 (Maxwell) or xi^2 ln(1 + B^2/2 xi^2) (Log). g comes
straight from the P_m column, so the area condition needs no integration of M.

Every coexistence pair is a bridging edge of the lower convex hull of g(B):
the hull gives the starting point (and handles overlapping loops, where the
global tangent spans several of them), then Newton refines (B1, B2) on the
interpolated data. With Log above sqrt(2) xi, g stays concave, so the bridge
runs to the end of the scan: no closing B2 (flagged open = true).

Usage:
  nsrs study hb --densities 2,4 --bmin 1e16 --bmax 1e19 --per-decade 10000
  python plot_scripts/hb_phase.py --hb results/study_hb/hb.csv
Writes <out>/hb_domains_phase.csv and fig_hb_phase.{pdf,png}.
"""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

STYLE = Path(__file__).with_name("paper.mplstyle")
if STYLE.exists():
    plt.style.use(STYLE)

ERG_PER_MEV_FM3 = 1.602176634e33  # MeV/fm^3 -> erg/cm^3 = G^2
PAGE_W = 7.0
MODEL_COLOR = {"GM1": "#0072B2", "GM3": "#D55E00", "FSU2": "#009E73"}


def load(path: Path):
    """{(model, n_ref): (B, P_m [G^2], 4 pi M [G], n/n0)} sorted in B."""
    data = defaultdict(list)
    for r in csv.DictReader(open(path, encoding="utf-8")):
        data[(r["model"], r["n_ref_over_n0"])].append(
            (float(r["B_G"]), float(r["P_m_MeV_fm3"]), float(r["four_pi_M_G"]), float(r["n_over_n0"])))
    out = {}
    for key, rows in data.items():
        a = np.array(sorted(rows))
        out[key] = (a[:, 0], a[:, 1] * ERG_PER_MEV_FM3, a[:, 2], a[:, 3])
    return out


def vacuum(b, xi):
    """(4 pi eps_vac, H_vac) in G^2 and G."""
    if not math.isfinite(xi):
        return 0.5 * b * b, b
    return xi * xi * np.log1p(b * b / (2 * xi * xi)), b / (1 + b * b / (2 * xi * xi))


def lower_hull(x, y):
    """Indices of the lower convex hull of (x, y), x increasing (monotone chain)."""
    hull = []
    for i in range(len(x)):
        while len(hull) >= 2:
            i0, i1 = hull[-2], hull[-1]
            # remove i1 if it lies on or above the chord i0 -> i
            if (x[i1] - x[i0]) * (y[i] - y[i0]) - (y[i1] - y[i0]) * (x[i] - x[i0]) <= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    return hull


def coexistence(b, p, four_pi_m, xi, min_points=10):
    """List of dicts (B1, B2, H*, residuals, open) for every bridging hull edge."""
    gv, hv = vacuum(b, xi)
    p = p - p[0]  # constant offset: irrelevant for the tangent, avoids cancellation
    g = gv - 4 * math.pi * p
    h = hv - four_pi_m
    # work relative to the first point so the hull test is well conditioned
    hull = lower_hull(b, g - g[0])
    lb = np.log(b)

    def gi(x):
        return np.interp(np.log(x), lb, g)

    def hi(x):
        return np.interp(np.log(x), lb, h)

    out = []
    for i0, i1 in zip(hull, hull[1:]):
        if i1 - i0 < min_points:
            continue
        is_open = i1 == len(b) - 1
        b1, b2 = b[i0], b[i1]
        if not is_open:
            scale = b1

            def eqs(v):
                x1, x2 = v * scale
                return [(hi(x1) - hi(x2)) / scale,
                        (gi(x2) - gi(x1) - hi(x1) * (x2 - x1)) / scale**2]

            sol, info, ok, _ = fsolve(eqs, [b1 / scale, b2 / scale], full_output=True, xtol=1e-12)
            if ok == 1 and b[0] < sol[0] * scale < sol[1] * scale < b[-1]:
                b1, b2 = sol * scale
        hs = hi(b1)
        r1 = (hi(b1) - hi(b2)) / hs
        r2 = (gi(b2) - gi(b1) - hs * (b2 - b1)) / (hs * (b2 - b1))
        out.append(dict(B1=b1, B2=b2, H=hs, dB_over_B=(b2 - b1) / math.sqrt(b1 * b2),
                        res_H=r1, res_area=r2, open=is_open,
                        converged=is_open or max(abs(r1), abs(r2)) < 1e-6))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hb", type=Path, default=Path("results/study_hb/hb.csv"))
    ap.add_argument("--xi", type=float, default=1e17)
    ap.add_argument("--out", type=Path, default=Path("results/paper"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = load(args.hb)
    cases = [("maxwell", math.inf), (f"log:{args.xi:.2e}", args.xi)]

    rows = []
    for (model, nref), (b, p, fm, n) in sorted(data.items()):
        for label, xi in cases:
            # global construction over the whole scan
            for d in coexistence(b, p, fm, xi):
                rows.append(dict(model=model, n_ref_over_n0=nref, nlem=label, scope="global", **d))
            if math.isfinite(xi):
                # Log: H_vac peaks at sqrt(2) xi; below it the dHvA loops close
                # locally (metastable domains, since the global bridge is open)
                m = b < math.sqrt(2) * xi
                for d in coexistence(b[m], p[m], fm[m], xi):
                    if not d["open"]:
                        rows.append(dict(model=model, n_ref_over_n0=nref, nlem=label, scope="below_sqrt2_xi", **d))
    path = args.out / "hb_domains_phase.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model", "n_ref_over_n0", "nlem", "B1_G", "B2_G", "H_G", "dB_over_B",
                    "residual_H", "residual_area", "open", "converged", "scope"])
        for r in rows:
            w.writerow([r["model"], r["n_ref_over_n0"], r["nlem"], f"{r['B1']:.8e}", f"{r['B2']:.8e}",
                        f"{r['H']:.8e}", f"{r['dB_over_B']:.6e}", f"{r['res_H']:.2e}",
                        f"{r['res_area']:.2e}", str(r["open"]).lower(),
                        str(r["converged"]).lower(), r["scope"]])
    print(path, len(rows), "coexistence pairs")

    # Summary per case
    for (model, nref), (b, *_rest) in sorted(data.items()):
        for label, _ in cases:
            sel = [r for r in rows if r["model"] == model and r["n_ref_over_n0"] == nref and r["nlem"] == label]
            closed = [r for r in sel if not r["open"] and r["converged"]]
            unconverged = sum(1 for r in sel if not r["converged"])
            glob = [r for r in sel if r["scope"] == "global" and r["converged"]]
            span = sum(math.log(r["B2"] / r["B1"]) for r in glob) / math.log(b[-1] / b[0])
            above = [r for r in closed if r["B1"] > 1e17]
            opened = [r for r in sel if r["open"]]
            print(f"{model} n={nref} {label:14s} closed pairs={len(closed):4d} (not converged: {unconverged}) "
                  f"global log-coverage={100 * span:5.1f}% "
                  f"max dB/B(B1>1e17)={max((r['dB_over_B'] for r in above), default=0):.3g} "
                  f"max|res|={max((max(abs(r['res_H']), abs(r['res_area'])) for r in closed), default=0):.1e}"
                  + (f" open bridge from B1={min(r['B1'] for r in opened):.3e}" if opened else ""))

    # Figure: (a) Maxwell, (b) Log (closed pairs below sqrt(2) xi; above it no B2)
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(PAGE_W, 3.0), sharey=True)
    for ax, (label, xi) in zip(axes, cases):
        for (model, nref) in sorted(data):
            sel = [r for r in rows if r["model"] == model and r["n_ref_over_n0"] == nref and r["nlem"] == label
                   and r["converged"] and not r["open"]]
            bm = np.array([math.sqrt(r["B1"] * r["B2"]) for r in sel])
            db = np.array([r["dB_over_B"] for r in sel])
            ax.scatter(bm, db, s=5, color=MODEL_COLOR.get(model, "k"),
                       marker="o" if nref.startswith("4") else "^", alpha=0.6, lw=0)
        if math.isfinite(xi):
            ax.axvspan(math.sqrt(2) * xi, b_max(data), color="0.85", lw=0, zorder=0)
            ax.axvline(math.sqrt(2) * xi, color="0.3", ls="-.", lw=0.7)
            ax.text(0.97, 0.95, "$B>\\sqrt{2}\\,\\xi$:\nno closing $B_2$", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5, color="0.25")
        ax.set(xscale="log", yscale="log", xlabel=r"$\sqrt{B_1B_2}$ [G]", xlim=(8e15, 1.2e19), ylim=(1e-3, 0.15))
    axes[0].set_ylabel(r"domain contrast $(B_2-B_1)/\sqrt{B_1B_2}$")
    axes[0].set_title("(a) Maxwell", fontsize=9)
    axes[1].set_title(rf"(b) Log, $\xi=10^{{{int(round(math.log10(args.xi)))}}}$ G", fontsize=9)
    handles = [Line2D([], [], color=c, marker="o", ls="", ms=4, label=m) for m, c in MODEL_COLOR.items()]
    handles += [Line2D([], [], color="0.3", marker="o", ls="", ms=4, label=r"$n_B\approx4\,n_0$"),
                Line2D([], [], color="0.3", marker="^", ls="", ms=4, label=r"$n_B\approx2\,n_0$")]
    axes[0].legend(handles=handles, fontsize=7, frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(args.out / f"fig_hb_phase.{ext}", dpi=200)
    print(args.out / "fig_hb_phase.png")


def b_max(data):
    return max(v[0][-1] for v in data.values())


if __name__ == "__main__":
    main()
