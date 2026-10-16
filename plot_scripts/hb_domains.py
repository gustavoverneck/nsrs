"""H(B) of magnetized matter at fixed mu_n from `nsrs study hb`, and the
Condon-domain (equal-area) construction on a real T = 0 de Haas-van Alphen loop.

H = H_vac(B) - 4 pi M(B), with H_vac = B (Maxwell) or B/(1 + B^2/2 xi^2) (Log),
M = dP_m/dB at fixed mu_n (column four_pi_M_G = 4 pi M in G).

Panels:
 (a) 4 pi M/B over the full scan (dHvA oscillations of the matter);
 (b) zoom on H(B) with the equal-area construction on one loop: B1, B2;
 (c) local stability s = dH/dB = f_vac - c_m for Maxwell and Log;
 (d) Log: the vacuum H_vac dominates H above sqrt(2) xi.

Usage:
  nsrs study hb --models GM1 --densities 4 --bmin 1e16 --bmax 1e19 --per-decade 3000
  nsrs study hb --models GM1 --densities 4 --bmin 1.385e17 --bmax 1.415e17 --per-decade 1000000 --out results/study_hb_zoom
  python plot_scripts/hb_domains.py --hb results/study_hb/hb.csv --zoom results/study_hb_zoom/hb.csv
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hb_phase import ERG_PER_MEV_FM3, coexistence

STYLE = Path(__file__).with_name("paper.mplstyle")
if STYLE.exists():
    plt.style.use(STYLE)

PAGE_W = 7.0


def load(path: Path, model: str, density: str):
    rows = [r for r in csv.DictReader(open(path, encoding="utf-8"))
            if r["model"] == model and float(r["n_ref_over_n0"]) == float(density)]
    b = np.array([float(r["B_G"]) for r in rows])
    four_pi_m = np.array([float(r["four_pi_M_G"]) for r in rows])
    p_m = np.array([float(r["P_m_MeV_fm3"]) for r in rows]) * ERG_PER_MEV_FM3
    order = np.argsort(b)
    return b[order], four_pi_m[order], float(rows[0]["mu_n"]) if rows else math.nan, p_m[order]


def h_vac(b, xi):
    return b if not math.isfinite(xi) else b / (1.0 + b * b / (2.0 * xi * xi))


def f_vac(b, xi):
    if not math.isfinite(xi):
        return np.ones_like(b)
    x = b * b / (2.0 * xi * xi)
    return (1.0 - x) / (1.0 + x) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hb", type=Path, default=Path("results/study_hb/hb.csv"))
    ap.add_argument("--zoom", type=Path, default=Path("results/study_hb_zoom/hb.csv"))
    ap.add_argument("--model", default="GM1")
    ap.add_argument("--density", default="4")
    ap.add_argument("--xi", type=float, default=1e17)
    ap.add_argument("--out", type=Path, default=Path("results/paper"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    b, fm, mu, _ = load(args.hb, args.model, args.density)
    bz, fmz, _, pz = load(args.zoom, args.model, args.density)
    xi = args.xi

    fig, axes = plt.subplots(2, 2, figsize=(PAGE_W, 5.2))
    (a, bx), (cx, dx) = axes

    # (a) 4 pi M / B
    a.plot(b, fm / b, color="k", lw=0.4)
    a.axhline(0, color="0.5", lw=0.5)
    a.set(xscale="log", xlabel=r"$B$ [G]", ylabel=r"$4\pi\mathcal{M}/B$")
    a.axvspan(bz.min(), bz.max(), color="#e6550d", alpha=0.3, lw=0)
    a.set_title("(a) Matter magnetization (dHvA)", fontsize=9)

    # (b) zoom: H = B - 4 pi M (Maxwell vacuum), equal-area construction
    h = bz - fmz
    unit = 1e17
    bx.plot(bz / unit, h / unit, color="k", lw=1.0)
    # common tangent of g = B^2/2 - 4 pi P_m (2x2 system, see hb_phase.py); widest closed pair
    pairs = [d for d in coexistence(bz, pz, fmz, math.inf) if d["converged"] and not d["open"]]
    eq = max(pairs, key=lambda d: d["B2"] - d["B1"]) if pairs else None
    if eq:
        b1, b2, hs = eq["B1"], eq["B2"], eq["H"]
        m = (bz >= b1) & (bz <= b2)
        bx.fill_between(bz[m] / unit, h[m] / unit, hs / unit, where=h[m] > hs, color="#9ecae1", lw=0)
        bx.fill_between(bz[m] / unit, h[m] / unit, hs / unit, where=h[m] < hs, color="#fdae6b", lw=0)
        bx.plot([b1 / unit, b2 / unit], [hs / unit, hs / unit], color="#08519c", lw=1.4)
        bx.plot([b1 / unit, b2 / unit], [hs / unit, hs / unit], "o", color="#08519c", ms=4)
        bx.annotate(r"$B_1$", (b1 / unit, hs / unit), (-4, -11), textcoords="offset points", color="#08519c", fontsize=8)
        bx.annotate(r"$B_2$", (b2 / unit, hs / unit), (2, 4), textcoords="offset points", color="#08519c", fontsize=8)
        bx.text(0.03, 0.95, rf"$B_1={b1 / unit:.4f}\times10^{{17}}$ G" "\n"
                rf"$B_2={b2 / unit:.4f}\times10^{{17}}$ G" "\n"
                rf"$(B_2-B_1)/B={100 * eq['dB_over_B']:.2f}\%$", transform=bx.transAxes, va="top", fontsize=7)
        print(f"equal-area: B1={b1:.6e} B2={b2:.6e} H*={hs:.6e} dB/B={(b2 - b1) / hs:.3e}")
    bx.set(xlabel=r"$B$ [$10^{17}$ G]", ylabel=r"$H$ [$10^{17}$ G]")
    bx.set_title("(b) Condon domains from the matter (Maxwell)", fontsize=9)

    # (c) s = f_vac - c_m
    c_m = np.gradient(fm, b)
    for lab, x, col, alpha in [("Maxwell", math.inf, "k", 1.0),
                               (rf"Log, $\xi=10^{{{int(round(math.log10(xi)))}}}$ G", xi, "#d95f02", 0.35)]:
        s = f_vac(b, x) - c_m
        cx.plot(b, s, color=col, lw=0.4, alpha=alpha, label=f"{lab}: $s<0$ in {100 * (s < 0).mean():.0f}% of $B$")
    cx.axhline(0, color="0.5", lw=0.5)
    cx.axvline(math.sqrt(2) * xi, color="#d95f02", ls="-.", lw=0.7)
    cx.set(xscale="log", yscale="symlog", xlabel=r"$B$ [G]", ylabel=r"$s=dH/dB=f_{\rm vac}-c_m$")
    leg = cx.legend(fontsize=7, loc="upper right", frameon=True, framealpha=0.9)
    for line in leg.get_lines():
        line.set_alpha(1.0)
        line.set_linewidth(1.2)
    cx.set_title("(c) Local stability of a uniform field", fontsize=9)

    # (d) Log: total H with matter, vacuum dominates
    for lab, x, col, ls in [("Maxwell", math.inf, "k", ":"), ("Log", xi, "#d95f02", "-")]:
        dx.plot(b / xi, (h_vac(b, x) - fm) / xi, color=col, ls=ls, lw=1.0, label=lab + r", $H=H_{\rm vac}-4\pi\mathcal{M}$")
    dx.axvline(math.sqrt(2), color="#d95f02", ls="-.", lw=0.7)
    dx.set(xscale="log", yscale="log", xlabel=r"$B/\xi$", ylabel=r"$H/\xi$", ylim=(1e-2, 1e2))
    dx.legend(fontsize=7, loc="upper left", frameon=False)
    dx.set_title(r"(d) Log vacuum: no closing $B_2$", fontsize=9)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(args.out / f"fig_hb_domains.{ext}", dpi=200)
    print(args.out / "fig_hb_domains.png", f"mu_n={mu}")


if __name__ == "__main__":
    main()
