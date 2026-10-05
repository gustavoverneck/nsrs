#!/usr/bin/env python3
"""
Comparação sem x com momentos magnéticos anômalos (AMM) para a eletrodinâmica logarítmica.

Lê as saídas de duas rodadas (sem e com `--amm`) e gera, por modelo:
- amm_study_log_<M>.png: M_max(xi) e R_1.4(xi) para cada B0 (`study log`), com Maxwell e
  sem tensão como linhas horizontais; contínuo = sem AMM, tracejado = com AMM; a última
  linha mostra Delta M_max = M_max(AMM) - M_max(sem AMM).
- amm_study_b_<M>.png: M_max(B), R_1.4(B) e Lambda_1.4(B) (`study b`), uma cor por
  eletrodinâmica (Maxwell, Log com cada xi), perfil e pressão escolhidos por opção;
  contínuo = sem AMM, tracejado = com AMM.
- amm_summary.csv: as diferenças (AMM - sem AMM) ponto a ponto das duas tabelas.

Uso:
  python plot_scripts/amm_compare.py \\
      --study-log results/study_log/summary.csv results/study_log_amm/summary.csv \\
      --study-b results/study_b results/study_b_amm --out results/amm_compare
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
AMM_STYLE = {False: ("-", "sem AMM"), True: ("--", "com AMM")}
TOPOLOGY_TITLE = {
    "anisotropica": r"Anisotrópica ($P_\perp$)",
    "isotropica": r"Isotrópica ($(P_\parallel+2P_\perp)/3$)",
    "perp": r"Anisotrópica ($P_\perp$)",
    "iso": r"Isotrópica ($(P_\parallel+2P_\perp)/3$)",
}


def number(text: str | None) -> float:
    try:
        return float(text) if text else math.nan
    except ValueError:
        return math.nan


def load(path: Path) -> list[dict]:
    if not path.exists():
        print(f"[aviso] {path} não encontrado")
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def amm_legend() -> list[Line2D]:
    return [Line2D([], [], ls=ls, color="0.3", label=label) for ls, label in AMM_STYLE.values()]


def save(fig, path: Path) -> None:
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {path}")


# --------------------------------------------------------------------------- study log
def study_log(off: list[dict], on: list[dict], out: Path, summary: list[list]) -> None:
    runs = {False: off, True: on}
    models = sorted({r["model"] for r in off + on})
    for model in models:
        topologies = [t for t in ("anisotropica", "isotropica") if any(r["topology"] == t for r in off + on)]
        fields = sorted({float(r["b0_G"]) for r in off + on if r["model"] == model})
        color = {b0: COLORS[i % len(COLORS)] for i, b0 in enumerate(fields)}
        fig, axes = plt.subplots(3, len(topologies), figsize=(6.2 * len(topologies), 11), squeeze=False)
        for col, topology in enumerate(topologies):
            ax_m, ax_r, ax_d = axes[0][col], axes[1][col], axes[2][col]
            for b0 in fields:
                curves = {}
                for amm, rows in runs.items():
                    group = [r for r in rows if r["model"] == model and r["topology"] == topology
                             and float(r["b0_G"]) == b0]
                    logs = sorted((number(r["xi_G"]), number(r["m_max_Msun"]), number(r["r14_km"]))
                                  for r in group if r["case"] == "log")
                    if not logs:
                        continue
                    curves[amm] = {xi: (m, r14) for xi, m, r14 in logs}
                    ls = AMM_STYLE[amm][0]
                    xs = [p[0] for p in logs]
                    ax_m.plot(xs, [p[1] for p in logs], ls=ls, marker="o", ms=2.5, lw=1.6, color=color[b0])
                    ax_r.plot(xs, [p[2] for p in logs], ls=ls, marker="o", ms=2.5, lw=1.6, color=color[b0])
                    for case, style in (("maxwell", ":"), ("sem_tensao", (0, (1, 3)))):
                        ref = next((number(r["m_max_Msun"]) for r in group if r["case"] == case), math.nan)
                        if not math.isnan(ref) and not amm:
                            ax_m.axhline(ref, ls=style, lw=0.9, color=color[b0], alpha=0.7)
                if len(curves) == 2:
                    common = sorted(set(curves[False]) & set(curves[True]))
                    dm = [1e3 * (curves[True][xi][0] - curves[False][xi][0]) for xi in common]
                    ax_d.plot(common, dm, "o-", ms=2.5, lw=1.6, color=color[b0], label=f"$B_0$ = {b0:.1e} G")
                    for xi in common:
                        (m0, r0), (m1, r1) = curves[False][xi], curves[True][xi]
                        summary.append(["study_log", model, topology, f"{b0:.4e}", f"log:{xi:.4e}", "",
                                        f"{m0:.5f}", f"{m1:.5f}", f"{m1 - m0:.5f}", f"{r0:.4f}", f"{r1:.4f}",
                                        f"{r1 - r0:.4f}", "", "", ""])
            for ax in (ax_m, ax_r, ax_d):
                ax.set_xscale("log")
                ax.grid(alpha=0.25, lw=0.6)
            ax_m.set_title(TOPOLOGY_TITLE[topology])
            ax_m.set_ylabel(r"$M_{\max}$ [$M_\odot$]")
            ax_r.set_ylabel(r"$R_{1.4}$ [km]")
            ax_d.set_ylabel(r"$\Delta M_{\max}^{\rm AMM}$ [$10^{-3}\,M_\odot$]")
            ax_d.axhline(0.0, color="0.5", lw=0.8)
            ax_d.set_xlabel(r"$\xi$ [G]")
            ax_d.legend(frameon=False, fontsize=8)
        handles = [Line2D([], [], color=color[b0], label=f"$B_0$ = {b0:.1e} G") for b0 in fields]
        handles += amm_legend()
        handles += [Line2D([], [], ls=":", color="0.3", label=r"Maxwell ($\xi\to\infty$)"),
                    Line2D([], [], ls=(0, (1, 3)), color="0.3", label=r"sem tensão ($\xi\to0$)")]
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.03))
        fig.suptitle(f"{model}: Log($\\xi$) sem e com momentos anômalos")
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        save(fig, out / f"amm_study_log_{model}.png")


# --------------------------------------------------------------------------- study b
def study_b(off: list[dict], on: list[dict], out: Path, profile: str, summary: list[list]) -> None:
    runs = {False: off, True: on}
    models = sorted({r["model"] for r in off + on})
    nlems = list(dict.fromkeys(r.get("nlem", "maxwell") for r in off + on))
    color = {n: COLORS[i % len(COLORS)] for i, n in enumerate(nlems)}
    quantities = [("m_max_Msun", r"$M_{\max}$ [$M_\odot$]"), ("r14_km", r"$R_{1.4}$ [km]"),
                  ("lambda14", r"$\Lambda_{1.4}$")]
    for model in models:
        topologies = [t for t in ("perp", "iso") if any(r["topology"] == t for r in off + on)]
        fig, axes = plt.subplots(len(quantities), len(topologies),
                                 figsize=(6.2 * len(topologies), 3.4 * len(quantities)), squeeze=False)
        for col, topology in enumerate(topologies):
            for nlem in nlems:
                series = {}
                for amm, rows in runs.items():
                    pts = sorted((number(r["B_G"]), r) for r in rows
                                 if r["model"] == model and r.get("nlem", "maxwell") == nlem
                                 and r["profile"] == profile and r["topology"] == topology
                                 and number(r["B_G"]) > 0)
                    if not pts:
                        continue
                    series[amm] = {b: r for b, r in pts}
                    for row, (key, _) in enumerate(quantities):
                        axes[row][col].plot([b for b, _ in pts], [number(r[key]) for _, r in pts],
                                            ls=AMM_STYLE[amm][0], lw=1.6, color=color[nlem])
                if len(series) == 2:
                    for b in sorted(set(series[False]) & set(series[True])):
                        a, c = series[False][b], series[True][b]
                        d = lambda k: number(c[k]) - number(a[k])
                        summary.append(["study_b", model, f"{profile}/{topology}", f"{b:.4e}", nlem, "",
                                        a["m_max_Msun"], c["m_max_Msun"], f"{d('m_max_Msun'):.5f}",
                                        a["r14_km"], c["r14_km"], f"{d('r14_km'):.4f}",
                                        a["lambda14"], c["lambda14"], f"{d('lambda14'):.2f}"])
            for row, (_, label) in enumerate(quantities):
                ax = axes[row][col]
                ax.set_xscale("log")
                ax.set_ylabel(label)
                ax.grid(alpha=0.25, lw=0.6)
            axes[0][col].set_title(TOPOLOGY_TITLE[topology])
            axes[-1][col].set_xlabel(r"$B$ [G]")
        handles = [Line2D([], [], color=color[n], label=n) for n in nlems] + amm_legend()
        fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), frameon=False,
                   bbox_to_anchor=(0.5, -0.03))
        fig.suptitle(f"{model}: perfil {profile}, sem e com momentos anômalos")
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        save(fig, out / f"amm_study_b_{model}_{profile}.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study-log", nargs=2, type=Path, metavar=("SEM_AMM", "COM_AMM"),
                        help="summary.csv de `study log` sem e com --amm")
    parser.add_argument("--study-b", nargs=2, type=Path, metavar=("SEM_AMM", "COM_AMM"),
                        help="pastas de `study b` sem e com --amm")
    parser.add_argument("--profiles", nargs="+", default=["constante", "bdd"])
    parser.add_argument("--out", type=Path, default=Path("results/amm_compare"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    summary: list[list] = []
    if args.study_log:
        study_log(load(args.study_log[0]), load(args.study_log[1]), args.out, summary)
    if args.study_b:
        off, on = (load(folder / "stars.csv") for folder in args.study_b)
        for profile in args.profiles:
            study_b(off, on, args.out, profile, summary)

    path = args.out / "amm_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["study", "model", "topology", "b_G", "nlem", "", "m_max_sem", "m_max_amm", "dm_max",
                         "r14_sem", "r14_amm", "dr14", "lambda14_sem", "lambda14_amm", "dlambda14"])
        writer.writerows(summary)
    print(f"tabela: {path}")


if __name__ == "__main__":
    main()
