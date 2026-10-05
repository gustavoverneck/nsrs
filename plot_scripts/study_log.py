#!/usr/bin/env python3
"""
Figuras do estudo `nsrs study log` (impacto da eletrodinâmica logarítmica na estrela).

Lê results/study_log/summary.csv e gera, por modelo:
- M_max em função de xi, uma curva por B0, com os limites de Maxwell (xi -> inf,
  pontilhado) e sem tensão do campo (xi -> 0, tracejado); a faixa cinza marca
  os xi para os quais a EoS sai vazia (pressão do campo negativa já na superfície);
- Delta M_max em relação ao caso sem tensão, em função de x_c = B_c^2 / (2 xi^2)
  no centro da estrela de massa máxima, com a linha x = x* (pressão do campo < 0).

Uso: python plot_scripts/study_log.py [summary.csv] [--out results/study_log]
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

TOPOLOGIES = ["anisotropica", "isotropica"]
TITLES = {"anisotropica": r"Anisotrópica ($P_\perp$)", "isotropica": r"Isotrópica ($(P_\parallel+2P_\perp)/3$)"}
COLORS = ["#2a6fdb", "#e0703a", "#2e9e6a", "#9b59b6", "#c0392b"]


def number(text: str) -> float | None:
    return float(text) if text else None


def load(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot_model(rows: list[dict], model: str, out_dir: Path) -> list[Path]:
    groups: dict[tuple[str, float], dict] = defaultdict(lambda: {"log": []})
    for row in rows:
        key = (row["topology"], float(row["b0_G"]))
        if row["case"] == "log":
            groups[key]["log"].append(row)
        else:
            groups[key][row["case"]] = row
    topologies = [t for t in TOPOLOGIES if any(k[0] == t for k in groups)]
    fields = sorted({k[1] for k in groups})
    color = {b0: COLORS[i % len(COLORS)] for i, b0 in enumerate(fields)}
    written = []

    # 1. M_max(xi) com os dois limites.
    fig, axes = plt.subplots(1, len(topologies), figsize=(6.2 * len(topologies), 4.6), squeeze=False)
    for ax, topology in zip(axes[0], topologies):
        empty_xi = []
        for b0 in fields:
            group = groups.get((topology, b0))
            if not group:
                continue
            points = [(float(r["xi_G"]), number(r["m_max_Msun"])) for r in group["log"]]
            empty_xi += [xi for xi, m in points if m is None]
            valid = [(xi, m) for xi, m in points if m is not None]
            label = f"$B_0$ = {b0:.1e} G"
            ax.plot(*zip(*valid), "o-", ms=3.5, lw=1.8, color=color[b0], label=label)
            for case, style in (("maxwell", ":"), ("sem_tensao", "--")):
                ref = number(group.get(case, {}).get("m_max_Msun", ""))
                if ref is not None:
                    ax.axhline(ref, ls=style, lw=1.2, color=color[b0], alpha=0.8)
        if empty_xi:
            ax.axvspan(min(empty_xi) / 1.3, max(empty_xi) * 1.3, color="0.85", zorder=0,
                       label="EoS vazia")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\xi$ [G]")
        ax.set_ylabel(r"$M_{\max}$ [$M_\odot$]")
        ax.set_title(TITLES[topology])
        ax.grid(alpha=0.25, lw=0.6)
    handles, labels = axes[0][0].get_legend_handles_labels()
    handles += [plt.Line2D([], [], ls=":", color="0.3"), plt.Line2D([], [], ls="--", color="0.3")]
    labels += [r"Maxwell ($\xi\to\infty$)", r"sem tensão ($\xi\to0$)"]
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 5), frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{model}: massa máxima com eletrodinâmica logarítmica")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    path = out_dir / f"mmax_vs_xi_{model}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    written.append(path)

    # 2. Delta M_max (vs sem tensão) em função de x_c.
    fig, axes = plt.subplots(1, len(topologies), figsize=(6.2 * len(topologies), 4.6), squeeze=False)
    for ax, topology in zip(axes[0], topologies):
        x_star = None
        for b0 in fields:
            group = groups.get((topology, b0))
            if not group:
                continue
            points = sorted(
                (number(r["x_center"]), 1e3 * number(r["dm_max_vs_no_stress"]))
                for r in group["log"]
                if r["x_center"] and r["dm_max_vs_no_stress"]
            )
            x_star = x_star or next((number(r["x_threshold"]) for r in group["log"]), None)
            if points:
                ax.plot(*zip(*points), "o-", ms=3.5, lw=1.8, color=color[b0],
                        label=f"$B_0$ = {b0:.1e} G")
        if x_star:
            ax.axvline(x_star, color="0.35", lw=1.0, ls="--")
            ax.annotate(f"$x^*$ = {x_star:.2f}", (x_star, 0.02), xycoords=("data", "axes fraction"),
                        xytext=(4, 0), textcoords="offset points", fontsize=9, color="0.3")
        ax.axhline(0.0, color="0.5", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel(r"$x_c = B_c^2/2\xi^2$ (centro da estrela de $M_{\max}$)")
        ax.set_ylabel(r"$M_{\max} - M_{\max}^{\rm sem\,tensão}$ [$10^{-3}\,M_\odot$]")
        ax.set_title(TITLES[topology])
        ax.grid(alpha=0.25, lw=0.6)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"{model}: efeito da tensão do campo em função de $x_c$")
    fig.tight_layout()
    path = out_dir / f"dmmax_vs_xc_{model}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    written.append(path)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("summary", nargs="?", default="results/study_log/summary.csv", type=Path)
    parser.add_argument("--out", default=None, type=Path, help="pasta das figuras (padrão: a do CSV)")
    args = parser.parse_args()
    rows = load(args.summary)
    out_dir = args.out or args.summary.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    for model in sorted({r["model"] for r in rows}):
        for path in plot_model([r for r in rows if r["model"] == model], model, out_dir):
            print(f"figura: {path}")


if __name__ == "__main__":
    main()
