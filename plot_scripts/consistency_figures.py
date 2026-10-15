#!/usr/bin/env python3
"""Figuras de consistência da Log NLEM: onde o campo uniforme suposto pela EoS é estável.

O vácuo Log tem dH/dB = f_vac = (1 - x)/(1 + x)^2, x = B^2/(2 xi^2), negativo para
B > sqrt(2) xi. Acima disso um campo uniforme é localmente instável (formação de
domínios). Estas figuras mostram onde as estrelas do paper caem em relação a esse
limite. Só o termo do vácuo entra; o termo da matéria (c_m) está em `study profile`.

Campo no centro da estrela de massa máxima, como em `study profile` (campo que a
matéria sente): perfil constante, B_c = B0; perfil BDD, B_c = B(n_c). Linhas com
true_maximum=false (máximo além da EoS) ficam de fora, como em paper_figures.py.

Gera, em --out (padrão results/paper):
- fig_consistency_center: B_c/(sqrt(2) xi) em função de B0, um painel por modelo.
  Acima de 1 o centro da estrela de massa máxima tem f_vac < 0.
- fig_consistency_ximin: xi_min(B0), o menor xi com f_vac >= 0 no centro da estrela
  de massa máxima, resolvido entre os xi calculados (B_c(xi) = sqrt(2) xi);
  tracejado onde é extrapolado a partir do xi mais próximo.
- fig_consistency_radial: f_vac(r) e índices de refração do vácuo n - 1 ao longo do
  raio da estrela de massa máxima (perfil BDD), um painel por B0. n_perp diverge
  onde f_vac = 0 e não se propaga dentro desse raio.

Lê results/study_b/stars.csv e results/study_profile/radial.csv.
Uso: python plot_scripts/consistency_figures.py [--results results] [--out results/paper]
     [--model GM1] [--topology iso]
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from nlem_birefringence import log_indices  # noqa: E402
from paper_figures import MODELS, NLEM_STYLE, PAGE_W, fnum, load_study_b, panel_label  # noqa: E402

COL_W = 3.375
# src/core/magnetic.rs: B_SURFACE_G; src/core/constants.rs: BDD_BETAA, BDD_ALPHAA.
B_SURF, BDD_BETA, BDD_GAMMA = 1e15, 1e-2, 3.0
LOG_NLEMS = [k for k in NLEM_STYLE if k.startswith("log:")]
PROFILE_LS = {"constante": "-", "bdd": "--"}
PROFILE_LABEL = {"constante": "constant", "bdd": "BDD"}
MODEL_COLOR = {"GM1": "#0072B2", "GM3": "#D55E00", "FSU2": "#009E73"}


def xi_of(nlem: str) -> float:
    return float(nlem.split(":")[1])


def bdd_field(b0, n):
    return B_SURF + b0 * (1.0 - np.exp(-BDD_BETA * np.maximum(n, 0.0) ** BDD_GAMMA))


def f_vac(b, xi):
    x = np.asarray(b, dtype=float) ** 2 / (2.0 * xi * xi)
    return (1.0 - x) / (1.0 + x) ** 2


def center_field(row, profile):
    b0, nc = fnum(row["B_G"]), fnum(row["nc_over_n0"])
    return b0 if profile == "constante" else float(bdd_field(b0, nc))


def center_ratios(sb):
    """{(model, profile): {B0: {xi: B_c/(sqrt2 xi)}}} para as linhas válidas."""
    out = defaultdict(lambda: defaultdict(dict))
    for (model, nlem, profile), rows in sb.items():
        if nlem not in LOG_NLEMS:
            continue
        xi = xi_of(nlem)
        for r in rows:
            bc = center_field(r, profile)
            if math.isfinite(bc):
                out[(model, profile)][fnum(r["B_G"])][xi] = bc / (math.sqrt(2.0) * xi)
    return out


def xi_min(ratios: dict):
    """xi com B_c(xi) = sqrt(2) xi, interpolando log r em log xi entre xi calculados.
    Fora do intervalo, xi_min = B_c(xi mais próximo)/sqrt(2). Devolve (xi, extrapolado)."""
    pts = sorted(ratios.items())
    for (x1, r1), (x2, r2) in zip(pts, pts[1:]):
        if (r1 - 1.0) * (r2 - 1.0) <= 0.0:
            t = math.log(r1) / (math.log(r1) - math.log(r2))
            return math.exp(math.log(x1) + t * (math.log(x2) - math.log(x1))), False
    xi_ref, r_ref = pts[-1] if pts[-1][1] > 1.0 else pts[0]
    return r_ref * xi_ref, True


def fig_center(sb, out: Path):
    ratios = center_ratios(sb)
    fig, axes = plt.subplots(1, 3, figsize=(PAGE_W, 2.5), sharey=True)
    for k, (ax, model) in enumerate(zip(axes, MODELS)):
        ax.axhspan(1.0, 1e5, color="0.92", lw=0)
        ax.axhline(1.0, color="0.4", lw=0.7)
        for profile, ls in PROFILE_LS.items():
            by_b = ratios.get((model, profile), {})
            for nlem in LOG_NLEMS:
                xi = xi_of(nlem)
                b = np.array(sorted(b0 for b0, d in by_b.items() if xi in d))
                y = np.array([by_b[b0][xi] for b0 in b])
                ax.plot(b, y, color=NLEM_STYLE[nlem]["color"], ls=ls)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e15, 1e20)
        ax.set_ylim(1e-2, 1e4)
        ax.set_xlabel(r"$B_0$ (G)")
        panel_label(ax, k, model)
    axes[0].set_ylabel(r"$B_c/(\sqrt{2}\,\xi)$")
    axes[0].text(1.5e15, 1.6, r"$f_{\rm vac}<0$ at centre", fontsize=7, color="0.3")
    handles = [Line2D([], [], color=NLEM_STYLE[n]["color"], label=NLEM_STYLE[n]["label"]) for n in LOG_NLEMS]
    handles += [Line2D([], [], color="k", ls=ls, label=PROFILE_LABEL[p]) for p, ls in PROFILE_LS.items()]
    fig.legend(handles=handles, ncol=5, loc="upper center", frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, out, "fig_consistency_center")


def fig_ximin(sb, out: Path):
    ratios = center_ratios(sb)
    fig, ax = plt.subplots(figsize=(COL_W, 2.7))
    for model in MODELS:
        for profile, ls in PROFILE_LS.items():
            by_b = ratios.get((model, profile), {})
            pts = [(b0, *xi_min(d)) for b0, d in sorted(by_b.items()) if b0 > 0 and len(d) >= 1]
            if not pts:
                continue
            b = np.array([p[0] for p in pts])
            xm = np.array([p[1] for p in pts])
            ext = np.array([p[2] for p in pts])
            inner = np.where(ext, np.nan, xm)
            ax.plot(b, xm, color=MODEL_COLOR[model], ls=ls, lw=0.6, alpha=0.5)
            ax.plot(b, inner, color=MODEL_COLOR[model], ls=ls)
    for nlem in LOG_NLEMS:
        ax.axhline(xi_of(nlem), color=NLEM_STYLE[nlem]["color"], lw=0.5, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e15, 1e20)
    ax.set_ylim(1e14, 1e20)
    ax.set_xlabel(r"$B_0$ (G)")
    ax.set_ylabel(r"$\xi_{\min}$ (G)")
    ax.text(3e18, 3e14, "inconsistent\n" + r"($f_{\rm vac}<0$ at centre)", fontsize=7, color="0.3", ha="center")
    handles = [Line2D([], [], color=c, label=m) for m, c in MODEL_COLOR.items()]
    handles += [Line2D([], [], color="k", ls=ls, label=PROFILE_LABEL[p]) for p, ls in PROFILE_LS.items()]
    handles += [Line2D([], [], color="0.5", ls=":", lw=0.5, label=r"$\xi$ values of the scan")]
    ax.legend(handles=handles, loc="upper left", frameon=False)
    fig.tight_layout()
    save(fig, out, "fig_consistency_ximin")


def load_radial(path: Path, model: str):
    """{(nlem, B0): (r_km, n_over_n0, R_km, M)} da estrela de massa máxima, perfil BDD."""
    data = defaultdict(list)
    meta = {}
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if row["model"] != model or row["profile"] != "bdd" or row["star"] != "max":
            continue
        key = (row["nlem"], fnum(row["B0_G"]))
        data[key].append((fnum(row["r_km"]), fnum(row["n_over_n0"])))
        meta[key] = (fnum(row["R_km"]), fnum(row["M_Msun"]))
    out = {}
    for key, pts in data.items():
        pts.sort()
        r = np.array([p[0] for p in pts])
        n = np.nan_to_num(np.array([p[1] for p in pts]), nan=0.0)  # crosta: n = NaN
        out[key] = (r, n, *meta[key])
    return out


def fig_radial(radial_path: Path, out: Path, model: str):
    radial = load_radial(radial_path, model)
    b0s = sorted({b0 for (_, b0) in radial})
    if not b0s:
        print(f"sem perfis radiais BDD de {model} em {radial_path}")
        return
    fig, axes = plt.subplots(2, len(b0s), figsize=(PAGE_W, 4.4), sharex=True, squeeze=False)
    for col, b0 in enumerate(b0s):
        top, bottom = axes[0, col], axes[1, col]
        top.axhline(0.0, color="0.4", lw=0.6)
        for nlem in LOG_NLEMS:
            if (nlem, b0) not in radial:
                continue
            xi, color = xi_of(nlem), NLEM_STYLE[nlem]["color"]
            r, n, radius, mass = radial[(nlem, b0)]
            b = bdd_field(b0, n)
            fv = f_vac(b, xi)
            top.plot(r, fv, color=color, label=NLEM_STYLE[nlem]["label"] + f" ({mass:.2f} $M_\\odot$)")
            dpar, dperp, _ = log_indices(b, xi)
            bottom.plot(r, dpar, color=color, ls="-")
            bottom.plot(r, dperp, color=color, ls="--")
            neg = r[fv < 0]
            if neg.size:  # f_vac < 0 do centro até r0
                r0 = neg.max()
                for ax in (top, bottom):
                    ax.axvline(r0, color=color, lw=0.6, ls=":")
                top.fill_between(r, fv, 0.0, where=fv < 0, color=color, alpha=0.3, lw=0)
        top.set_ylim(-0.15, 1.05)
        top.set_ylabel(r"$f_{\rm vac}=dH/dB$" if col == 0 else "")
        panel_label(top, col, rf"{model}, BDD, $B_0=10^{{{round(math.log10(b0))}}}$ G")
        top.legend(loc="center left", bbox_to_anchor=(0.0, 0.72), frameon=False, fontsize=6.5)
        bottom.set_yscale("log")
        bottom.set_ylim(1e-6, 1e3)
        bottom.set_xlabel(r"$r$ (km)")
        bottom.set_ylabel(r"$n-1$" if col == 0 else "")
        panel_label(bottom, len(b0s) + col)
    axes[1, 0].legend(handles=[Line2D([], [], color="k", ls="-", label=r"$n_\parallel$"),
                               Line2D([], [], color="k", ls="--", label=r"$n_\perp$ (absent where $f_{\rm vac}<0$)")],
                      loc="lower left", frameon=False, fontsize=6.5)
    for ax in axes[1]:
        ax.set_xlim(0, None)
    fig.tight_layout()
    save(fig, out, "fig_consistency_radial")


def save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", dpi=300)
    plt.close(fig)
    print(out / f"{name}.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--out", type=Path, default=Path("results/paper"))
    parser.add_argument("--model", default="GM1", help="modelo do perfil radial")
    parser.add_argument("--topology", default="iso", choices=["iso", "perp"])
    args = parser.parse_args()

    sb = load_study_b(args.results / "study_b" / "stars.csv", args.topology)
    fig_center(sb, args.out)
    fig_ximin(sb, args.out)
    fig_radial(args.results / "study_profile" / "radial.csv", args.out, args.model)


if __name__ == "__main__":
    main()
