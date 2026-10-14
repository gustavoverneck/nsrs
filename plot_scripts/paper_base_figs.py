"""Figuras de base do artigo a partir de `nsrs study profile`:

  eos_B0_1e18.{pdf,png}       P(eps) isotrópica, (P_par + 2 P_perp)/3, por modelo e
                              eletrodinâmica, BDD x constante
  field_profile.{pdf,png}     B(n_B) do perfil BDD e do campo constante

Uso: python plot_scripts/paper_base_figs.py [results/study_profile] [--out results/paper_total]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from magnetic_profiles import MODELS, NLEM_STYLE, PAGE_W, PROFILE_STYLE, load_csv_rows  # noqa: E402


def save(fig, out: Path, stem: str):
    for ext, kw in [("pdf", {}), ("png", {"dpi": 200})]:
        fig.savefig(out / f"{stem}.{ext}", **kw)
    plt.close(fig)
    print(f"figura: {out / stem}.pdf")


def eos_figure(rows, out: Path, b0="1.0000e18"):
    """Em cima P(eps) isotrópica; embaixo o desvio relativo a Maxwell com o mesmo perfil."""
    fig, axes = plt.subplots(2, len(MODELS), figsize=(PAGE_W, 4.2), sharex="col", sharey="row",
                             gridspec_kw=dict(height_ratios=[1.4, 1.0]))
    for col, model in enumerate(MODELS):
        top, bottom = axes[0, col], axes[1, col]
        for profile, pstyle in PROFILE_STYLE.items():
            ref = rows.get((model, b0, profile, "maxwell"))
            for nlem, nstyle in NLEM_STYLE.items():
                r = rows.get((model, b0, profile, nlem))
                if r is None:
                    continue
                top.plot(r["eps"], r["p_iso"], color=nstyle["color"], ls=pstyle["ls"], lw=0.9)
                if ref is not None and nlem != "maxwell":
                    dev = r["p_iso"] / np.interp(r["eps"], ref["eps"], ref["p_iso"]) - 1.0
                    bottom.plot(r["eps"], 100 * dev, color=nstyle["color"], ls=pstyle["ls"], lw=0.9)
        top.set_title(model)
        top.set_xscale("log")
        top.set_yscale("log")
        top.set_ylim(1e-2, None)
        bottom.axhline(0.0, color="0.5", lw=0.5, zorder=0)
        bottom.set_xlabel(r"$\varepsilon$ [MeV fm$^{-3}$]")
    axes[0, 0].set_ylabel(r"$\bar P$ [MeV fm$^{-3}$]")
    axes[1, 0].set_ylabel(r"$\bar P/\bar P^{\rm Maxwell}-1$ [\%]")
    handles = [Line2D([], [], color=s["color"], label=s["label"]) for s in NLEM_STYLE.values()]
    handles += [Line2D([], [], color="0.3", ls=s["ls"], label=s["label"]) for s in PROFILE_STYLE.values()]
    fig.legend(handles=handles, loc="upper center", ncol=6, fontsize=7, bbox_to_anchor=(0.5, 1.0),
               title=r"$B_0=10^{18}$ G", title_fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90), h_pad=0.4)
    save(fig, out, "eos_B0_1e18")


def field_figure(rows, out: Path):
    fig, ax = plt.subplots(figsize=(PAGE_W / 2, 2.6))
    colors = {"1.0000e17": "#0072B2", "1.0000e18": "#D55E00"}
    # O perfil BDD só depende de n_B: os três modelos coincidem.
    model_ls = {"GM1": "-", "GM3": "-", "FSU2": "-"}
    for b0, color in colors.items():
        for model, ls in model_ls.items():
            r = rows.get((model, b0, "bdd", "maxwell"))
            if r is not None:
                ax.plot(r["n"], r["B"], color=color, ls=ls, lw=0.9)
        ax.axhline(float(b0), color=color, lw=0.6, ls="-.", alpha=0.7)
    ax.set_yscale("log")
    ax.set_xlim(0.0, None)
    ax.set_xlabel(r"$n_B/n_0$")
    ax.set_ylabel(r"$B$ [G]")
    handles = [Line2D([], [], color=c, label=rf"$B_0=10^{{{b[-2:]}}}$ G") for b, c in colors.items()]
    handles += [Line2D([], [], color="0.3", label="BDD"), Line2D([], [], color="0.3", ls="-.", lw=0.6, label="constant")]
    ax.legend(handles=handles, fontsize=6.5, loc="lower right")
    fig.tight_layout()
    save(fig, out, "field_profile")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dir", nargs="?", type=Path, default=Path("results/study_profile"))
    parser.add_argument("--out", type=Path, default=Path("results/paper_total"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = load_csv_rows(args.dir / "eos.csv")
    eos_figure(rows, args.out)
    field_figure(rows, args.out)


if __name__ == "__main__":
    main()
