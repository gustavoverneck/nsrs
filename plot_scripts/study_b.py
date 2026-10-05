#!/usr/bin/env python3
"""
Figuras da varredura `nsrs study b` (propriedades e estabilidade em função de B).

Lê <pasta>/stars.csv e <pasta>/stability.csv (padrão: results/study_b) e gera, por modelo:
- study_b_compare_<M>.png: comparação de todas as configurações. Colunas: pressão da TOV
  (anisotrópica P_perp | isotrópica); linhas: M_max, R_1.4, Lambda_1.4; cor: eletrodinâmica
  (Maxwell, Log com cada xi); traço: perfil do campo (constante contínuo, BDD tracejado).
  B = 0 aparece como linha horizontal cinza.
- study_b_stars_<M>_<nlem>.png: M_max, R_1.4, Lambda_1.4 e cobertura da EoS (n_B máximo) em
  função de B, para cada perfil (constante, BDD) e pressão da TOV (P_perp, isotrópica), uma
  figura por eletrodinâmica. A linha vertical pontilhada marca o primeiro B com término
  anômalo da EoS de cada perfil.
- study_b_stability_<M>_<nlem>.png: plano (B, n_B/n0) com os pontos magneticamente instáveis
  (s = f_vac - c_m < 0 com os passos delta e delta/2: formação de domínios) e mecanicamente
  instáveis (dP_perp/dn_B < 0), e as trajetórias B(n_B) do perfil BDD para alguns B0. No
  caso Log, a linha vertical marca B = sqrt(2) xi (dH/dB do vácuo = 0).
- study_b_unstable_fraction_<M>.png: fração dos pontos de densidade instáveis por B, uma
  curva por eletrodinâmica.

CSVs antigos, sem a coluna `nlem`, são lidos como Maxwell.

Uso: python plot_scripts/study_b.py [pasta]
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STYLE = {
    ("constante", "perp"): ("#2a6fdb", "-", r"constante, $P_\perp$"),
    ("constante", "iso"): ("#2a6fdb", "--", "constante, isotrópica"),
    ("bdd", "perp"): ("#e0703a", "-", r"BDD, $P_\perp$"),
    ("bdd", "iso"): ("#e0703a", "--", "BDD, isotrópica"),
}
# Cores categóricas em ordem fixa (uma por eletrodinâmica, na ordem do CSV).
NLEM_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PROFILE_STYLE = {"constante": ("-", "campo constante"), "bdd": ("--", "perfil BDD")}
TOPOLOGY_TITLE = {"perp": r"Anisotrópica ($P_\perp$)", "iso": r"Isotrópica ($(P_\parallel + 2P_\perp)/3$)"}
B_SURFACE = 1e15


def number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return math.nan


def load(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row.setdefault("nlem", "maxwell")
    return rows


def nlem_name(nlem: str) -> str:
    if nlem == "maxwell":
        return "Maxwell"
    kind, _, value = nlem.partition(":")
    if kind == "log":
        mantissa, _, exponent = f"{number(value):.1e}".partition("e")
        return rf"Log, $\xi = {mantissa}\times10^{{{int(exponent)}}}$ G"
    if kind == "modmax":
        return rf"ModMax, $\gamma = {value}$"
    return nlem


def nlem_tag(nlem: str) -> str:
    return nlem.replace(":", "_").replace("+", "")


def nlem_xi(nlem: str) -> float | None:
    kind, _, value = nlem.partition(":")
    return number(value) if kind == "log" else None


def ordered(values) -> list[str]:
    seen: list[str] = []
    for value in values:
        if value not in seen:
            seen.append(value)
    return seen


def series_xy(data: list[dict], column: str) -> tuple[list[float], list[float]]:
    points = sorted((number(r["B_G"]), number(r[column])) for r in data if number(r["B_G"]) > 0.0)
    return [p[0] for p in points], [p[1] for p in points]


def plot_compare(rows: list[dict], model: str, out: Path) -> Path:
    panels = [
        ("m_max_Msun", r"$M_{\max}$ [$M_\odot$]"),
        ("r14_km", r"$R_{1.4}$ [km]"),
        ("lambda14", r"$\Lambda_{1.4}$"),
    ]
    nlems = ordered(r["nlem"] for r in rows)
    profiles = [p for p in PROFILE_STYLE if any(r["profile"] == p for r in rows)]
    series = defaultdict(list)
    for row in rows:
        series[(row["nlem"], row["profile"], row["topology"])].append(row)

    fig, axes = plt.subplots(len(panels), 2, figsize=(12, 10.5), sharex=True, sharey="row")
    for col, topology in enumerate(("perp", "iso")):
        axes[0][col].set_title(TOPOLOGY_TITLE[topology])
        for row_index, (column, label) in enumerate(panels):
            ax = axes[row_index][col]
            zero = [number(r[column]) for r in rows
                    if number(r["B_G"]) == 0.0 and r["topology"] == topology and r["nlem"] == nlems[0]]
            if zero and not math.isnan(zero[0]):
                ax.axhline(zero[0], color="0.6", lw=0.9, ls=":")
            for k, nlem in enumerate(nlems):
                color = NLEM_COLORS[k % len(NLEM_COLORS)]
                for profile in profiles:
                    xs, ys = series_xy(series.get((nlem, profile, topology), []), column)
                    if xs:
                        ax.plot(xs, ys, PROFILE_STYLE[profile][0], color=color, lw=2.0)
                xi = nlem_xi(nlem)
                if xi is not None and row_index == 0:
                    ax.axvline(math.sqrt(2.0) * xi, color=color, lw=0.8, ls=(0, (1, 3)), alpha=0.9)
            ax.set_xscale("log")
            ax.grid(alpha=0.25, lw=0.6)
            if column == "lambda14":
                ax.set_yscale("log")
            if col == 0:
                ax.set_ylabel(label)
    for ax in axes[-1]:
        ax.set_xlabel("B [G]  (perfil BDD: $B_0$)")
    handles = [Line2D([], [], color=NLEM_COLORS[k % len(NLEM_COLORS)], lw=2.4, label=nlem_name(n))
               for k, n in enumerate(nlems)]
    handles += [Line2D([], [], color="0.25", ls=PROFILE_STYLE[p][0], lw=2.0, label=PROFILE_STYLE[p][1])
                for p in profiles]
    handles.append(Line2D([], [], color="0.6", ls=":", lw=0.9, label="B = 0"))
    if any(nlem_xi(n) is not None for n in nlems):
        handles.append(Line2D([], [], color="0.4", ls=(0, (1, 3)), lw=0.8, label=r"$B = \sqrt{2}\,\xi$ (Log)"))
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 4), frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{model}: propriedades estelares em função de B — eletrodinâmica, perfil e pressão da TOV")
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    path = out / f"study_b_compare_{model}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_stars(rows: list[dict], model: str, nlem: str, out: Path) -> Path:
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
            xs, ys = series_xy(data, column)
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
    fig.suptitle(f"{model}, {nlem_name(nlem)}: propriedades em função de B "
                 "(linha cinza: B = 0; pontilhado vertical: EoS truncada)")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    path = out / f"study_b_stars_{model}_{nlem_tag(nlem)}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def unstable_flags(rows: list[dict]) -> tuple[list[float], list[float], list[bool], list[bool]]:
    b = [number(r["B_G"]) for r in rows]
    n = [number(r["n_over_n0"]) for r in rows]
    # Instável só se s < 0 com os dois passos (delta e delta/2).
    magnetic = [r.get("magnetically_unstable") == "true" for r in rows]
    mechanical = [number(r["dpperp_dn"]) < 0 for r in rows]
    return b, n, magnetic, mechanical


def plot_stability(rows: list[dict], model: str, nlem: str, out: Path) -> Path:
    b, n, magnetic, mechanical = unstable_flags(rows)
    stable = [(x, y) for x, y, m, d in zip(b, n, magnetic, mechanical) if not m and not d]
    mag = [(x, y) for x, y, m in zip(b, n, magnetic) if m]
    mech = [(x, y) for x, y, d in zip(b, n, mechanical) if d]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    for points, color, size, label in (
        (stable, "0.82", 2, "estável"),
        (mech, "#2a6fdb", 5, r"mecânica: $dP_\perp/dn_B < 0$"),
        (mag, "#c0392b", 5, r"magnética: $s = f_{vac} - c_m < 0$ (domínios)"),
    ):
        if points:
            ax.scatter(*zip(*points), s=size, color=color, label=label, linewidths=0)
    n_max = max((y for y in n if not math.isnan(y)), default=8.0)
    for b0, style in ((1e17, ":"), (1e18, "--"), (5e18, "-.")):
        grid = [n_max * k / 200 for k in range(1, 201)]
        field = [B_SURFACE + b0 * (1 - math.exp(-0.01 * x**3)) for x in grid]
        ax.plot(field, grid, style, color="0.2", lw=1.0, label=f"BDD $B_0$ = {b0:.0e} G")
    xi = nlem_xi(nlem)
    if xi is not None:
        ax.axvline(math.sqrt(2.0) * xi, color="#c0392b", lw=1.2, ls="-", alpha=0.6,
                   label=r"$B = \sqrt{2}\,\xi$: $dH_{vac}/dB = 0$")
    ax.set_xscale("log")
    ax.set_xlabel("B local [G]")
    ax.set_ylabel(r"$n_B/n_0$")
    ax.set_title(f"{model}, {nlem_name(nlem)}: estabilidade local (campo constante)")
    ax.legend(frameon=False, fontsize=8, markerscale=3)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    path = out / f"study_b_stability_{model}_{nlem_tag(nlem)}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_fraction(rows: list[dict], model: str, out: Path) -> Path:
    nlems = ordered(r["nlem"] for r in rows)
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for k, nlem in enumerate(nlems):
        color = NLEM_COLORS[k % len(NLEM_COLORS)]
        b, _, magnetic, mechanical = unstable_flags([r for r in rows if r["nlem"] == nlem])
        per_b = defaultdict(lambda: [0, 0, 0])
        for x, m, d in zip(b, magnetic, mechanical):
            per_b[x][0] += 1
            per_b[x][1] += m
            per_b[x][2] += d
        fields = sorted(per_b)
        ax.plot(fields, [100 * per_b[x][1] / per_b[x][0] for x in fields], "-", marker="o", ms=3,
                color=color, lw=2.0, label=f"{nlem_name(nlem)}: magnética")
        if any(per_b[x][2] for x in fields):
            ax.plot(fields, [100 * per_b[x][2] / per_b[x][0] for x in fields], "--", color=color, lw=1.4,
                    label=f"{nlem_name(nlem)}: mecânica")
    ax.set_xscale("log")
    ax.set_xlabel("B [G]")
    ax.set_ylabel("pontos de densidade instáveis [%]")
    ax.set_title(f"{model}: fração instável da EoS ($n_B \\geq 0.05\\,n_0$, campo constante)")
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path = out / f"study_b_unstable_fraction_{model}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/study_b")
    stars = load(folder / "stars.csv")
    stability = load(folder / "stability.csv")
    if not stars:
        sys.exit(f"{folder}/stars.csv não encontrado; rode `nsrs study b` antes.")
    for model in ordered(r["model"] for r in stars):
        of_model = [r for r in stars if r["model"] == model]
        print("figura:", plot_compare(of_model, model, folder))
        for nlem in ordered(r["nlem"] for r in of_model):
            print("figura:", plot_stars([r for r in of_model if r["nlem"] == nlem], model, nlem, folder))
        rows = [r for r in stability if r["model"] == model]
        if rows:
            for nlem in ordered(r["nlem"] for r in rows):
                print("figura:", plot_stability([r for r in rows if r["nlem"] == nlem], model, nlem, folder))
            print("figura:", plot_fraction(rows, model, folder))


if __name__ == "__main__":
    main()
