#!/usr/bin/env python3
"""
Figuras da varredura `nsrs study b` (propriedades e estabilidade em função de B).

Lê <pasta>/stars.csv e <pasta>/stability.csv (padrão: results/study_b) e gera, por modelo:
- study_b_stars_<M>.png: M_max, R_1.4, Lambda_1.4 e cobertura da EoS (n_B máximo) em
  função de B, para cada perfil (constante, BDD) e pressão da TOV (P_perp, isotrópica).
  B = 0 aparece como linha horizontal; a linha vertical pontilhada marca o primeiro B com
  término anômalo da EoS de cada perfil.
- study_b_stability_<M>.png: plano (B, n_B/n0) com os pontos magneticamente instáveis
  (s < 0 com os passos delta e delta/2: formação de domínios) e mecanicamente instáveis (dP_perp/dn_B < 0), e as
  trajetórias B(n_B) do perfil BDD para alguns B0.
- study_b_unstable_fraction_<M>.png: fração dos pontos de densidade instáveis por B.

Uso: python plot_scripts/study_b.py [pasta]
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

STYLE = {
    ("constante", "perp"): ("#2a6fdb", "-", r"constante, $P_\perp$"),
    ("constante", "iso"): ("#2a6fdb", "--", "constante, isotrópica"),
    ("bdd", "perp"): ("#e0703a", "-", r"BDD, $P_\perp$"),
    ("bdd", "iso"): ("#e0703a", "--", "BDD, isotrópica"),
}
B_SURFACE = 1e15


def number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return math.nan


def load(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def plot_stars(rows: list[dict], model: str, out: Path) -> Path:
    panels = [
        ("m_max_Msun", r"$M_{\max}$ [$M_\odot$]"),
        ("r14_km", r"$R_{1.4}$ [km]"),
        ("lambda14", r"$\Lambda_{1.4}$"),
        ("n_max_over_n0", r"$n_B$ máximo da EoS [$n_0$]"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True)
    series = defaultdict(list)
    for row in rows:
        series[(row["profile"], row["topology"])].append(row)
    for ax, (column, label) in zip(axes.flat, panels):
        for key, data in sorted(series.items()):
            color, style, name = STYLE.get(key, ("0.4", ":", f"{key[0]}, {key[1]}"))
            if column == "n_max_over_n0" and key[1] == "iso":
                continue  # a cobertura não depende da pressão da TOV
            zero = [number(r[column]) for r in data if number(r["B_G"]) == 0.0]
            points = sorted(
                (number(r["B_G"]), number(r[column])) for r in data if number(r["B_G"]) > 0.0
            )
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            ax.plot(xs, ys, style, color=color, lw=1.8, marker="o", ms=2.5,
                    label=name if column != "n_max_over_n0" else key[0])
            if zero and not math.isnan(zero[0]):
                ax.axhline(zero[0], color="0.55", lw=0.8, ls=":")
            anomalous = [number(r["B_G"]) for r in data if r["anomalous"] == "true" and number(r["B_G"]) > 0]
            if anomalous:
                ax.axvline(min(anomalous), color=color, lw=1.0, ls=":", alpha=0.8)
        ax.set_xscale("log")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25, lw=0.6)
        if column == "lambda14":
            ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("B [G]  (perfil BDD: $B_0$)")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"{model}: propriedades em função de B (linha cinza: B = 0; pontilhado vertical: EoS truncada)")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    path = out / f"study_b_stars_{model}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_stability(rows: list[dict], model: str, out: Path) -> list[Path]:
    b = [number(r["B_G"]) for r in rows]
    n = [number(r["n_over_n0"]) for r in rows]
    # Instável só se s < 0 com os dois passos (delta e delta/2).
    s = [-1.0 if r.get("magnetically_unstable") == "true" else 1.0 for r in rows]
    dp = [number(r["dpperp_dn"]) for r in rows]
    stable = [(x, y) for x, y, si, di in zip(b, n, s, dp) if si >= 0 and not di < 0]
    magnetic = [(x, y) for x, y, si in zip(b, n, s) if si < 0]
    mechanical = [(x, y) for x, y, di in zip(b, n, dp) if di < 0]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    for points, color, size, label in (
        (stable, "0.82", 2, "estável"),
        (mechanical, "#2a6fdb", 5, r"mecânica: $dP_\perp/dn_B < 0$"),
        (magnetic, "#c0392b", 5, r"magnética: $s < 0$ (domínios)"),
    ):
        if points:
            ax.scatter(*zip(*points), s=size, color=color, label=label, linewidths=0)
    n_max = max((y for y in n if not math.isnan(y)), default=8.0)
    for b0, style in ((1e17, ":"), (1e18, "--"), (5e18, "-.")):
        grid = [n_max * k / 200 for k in range(1, 201)]
        field = [B_SURFACE + b0 * (1 - math.exp(-0.01 * x**3)) for x in grid]
        ax.plot(field, grid, style, color="0.2", lw=1.0, label=f"BDD $B_0$ = {b0:.0e} G")
    ax.set_xscale("log")
    ax.set_xlabel("B local [G]")
    ax.set_ylabel(r"$n_B/n_0$")
    ax.set_title(f"{model}: estabilidade local da matéria (campo constante)")
    ax.legend(frameon=False, fontsize=8, markerscale=3)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    path_map = out / f"study_b_stability_{model}.png"
    fig.savefig(path_map, dpi=200, bbox_inches="tight")
    plt.close(fig)

    per_b = defaultdict(lambda: [0, 0, 0])
    for x, si, di in zip(b, s, dp):
        per_b[x][0] += 1
        per_b[x][1] += si < 0
        per_b[x][2] += di < 0
    fields = sorted(per_b)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(fields, [100 * per_b[x][1] / per_b[x][0] for x in fields], "o-", ms=3, color="#c0392b",
            label="magnética ($s < 0$)")
    ax.plot(fields, [100 * per_b[x][2] / per_b[x][0] for x in fields], "o-", ms=3, color="#2a6fdb",
            label=r"mecânica ($dP_\perp/dn_B < 0$)")
    ax.set_xscale("log")
    ax.set_xlabel("B [G]")
    ax.set_ylabel("pontos de densidade instáveis [%]")
    ax.set_title(f"{model}: fração instável da EoS ($n_B \\geq 0.05\\,n_0$)")
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    path_fraction = out / f"study_b_unstable_fraction_{model}.png"
    fig.savefig(path_fraction, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return [path_map, path_fraction]


def main() -> None:
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/study_b")
    stars = load(folder / "stars.csv")
    stability = load(folder / "stability.csv")
    if not stars:
        sys.exit(f"{folder}/stars.csv não encontrado; rode `nsrs study b` antes.")
    for model in sorted({r["model"] for r in stars}):
        print("figura:", plot_stars([r for r in stars if r["model"] == model], model, folder))
        rows = [r for r in stability if r["model"] == model]
        if rows:
            for path in plot_stability(rows, model, folder):
                print("figura:", path)


if __name__ == "__main__":
    main()
