"""Figuras de publicação (Physical Review D) para a NLEM log, sem e com AMM.

Lê apenas resultados já calculados:
  results/study_log{,_amm}/summary.csv     (study log)
  output/study_log{,_amm}/.../*_stars.txt   (curvas M-R, gravadas com --save-eos)
  results/study_b{,_amm}/stars.csv          (study b)
  output/nlem_log{,_amm}/.../eos.dat        (populações)

Gera PDF (vetorial) e PNG em results/paper/, com larguras de coluna (3.375")
e de página (7.0") do revtex. Os textos das figuras são em inglês.

Uso: python plot_scripts/paper_figures.py [--out results/paper]
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

HERE = Path(__file__).resolve().parent
plt.style.use(HERE / "paper.mplstyle")

PAGE_W = 7.0
MODELS = ["GM1", "GM3", "FSU2"]
B0S = ["1.0000e17", "5.0000e17", "1.0000e18"]
# Cor e marcador fixos por entidade (Okabe-Ito, validado para CVD).
B0_STYLE = {
    "1.0000e17": dict(color="#0072B2", marker="o"),
    "5.0000e17": dict(color="#D55E00", marker="s"),
    "1.0000e18": dict(color="#009E73", marker="^"),
}
NLEM_STYLE = {
    "maxwell": dict(color="#000000", label="Maxwell"),
    "log:1.00e16": dict(color="#0072B2", label=r"Log, $\xi=10^{16}$ G"),
    "log:1.00e17": dict(color="#D55E00", label=r"Log, $\xi=10^{17}$ G"),
    "log:1.00e18": dict(color="#009E73", label=r"Log, $\xi=10^{18}$ G"),
}
TOPO_LABEL = {"anisotropica": r"$P_\perp$", "isotropica": r"$(P_\parallel+2P_\perp)/3$"}


def sci_tex(value: float) -> str:
    exp = math.floor(math.log10(value))
    mant = value / 10**exp
    if abs(mant - 1) < 1e-9:
        return rf"10^{{{exp}}}"
    return rf"{mant:g}\times10^{{{exp}}}"


def b0_label(b0: str) -> str:
    return rf"$B_0={sci_tex(float(b0))}$ G"


def fnum(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def panel_label(ax, k: int, text: str = ""):
    ax.set_title(f"({chr(97 + k)}) {text}".rstrip(), loc="left", pad=3)


def top_legend(fig, handles, ncol: int, top: float, **layout):
    """tight_layout reservando o topo e legenda logo acima dos painéis."""
    fig.tight_layout(rect=(0, 0, 1, top), **layout)
    fig.legend(handles=handles, loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, top))


def save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.pdf")
    fig.savefig(out / f"{name}.png", dpi=300)
    plt.close(fig)
    print(f"figura: {out / name}.pdf")


# ---------------------------------------------------------------- study log


def load_study_log(path: Path):
    """{(model, topology, b0): {"none": row, "maxwell": row, "log": [(xi, row)]}}"""
    data = defaultdict(lambda: {"log": []})
    for row in csv.DictReader(open(path, encoding="utf-8")):
        key = (row["model"], row["topology"], row["b0_G"])
        if row["case"] == "sem_tensao":
            data[key]["none"] = row
        elif row["case"] == "maxwell":
            data[key]["maxwell"] = row
        else:
            data[key]["log"].append((float(row["xi_G"]), row))
    for entry in data.values():
        entry["log"].sort(key=lambda item: item[0])
    return data


def fig_mmax_vs_xi(sl, out: Path):
    """Delta M_max(xi) relativo ao caso sem tensão; Maxwell como assíntota."""
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 4.0), sharex=True)
    for col, model in enumerate(MODELS):
        for row, topo in enumerate(["anisotropica", "isotropica"]):
            ax = axes[row, col]
            for b0 in B0S:
                entry = sl[(model, topo, b0)]
                ref = fnum(entry["none"]["m_max_Msun"])
                xi = np.array([x for x, _ in entry["log"]])
                dm = np.array([fnum(r["m_max_Msun"]) - ref for _, r in entry["log"]]) * 1e3
                st = B0_STYLE[b0]
                ax.plot(xi, dm, color=st["color"], marker=st["marker"], markevery=4,
                        label=b0_label(b0))
                dmax = (fnum(entry["maxwell"]["m_max_Msun"]) - ref) * 1e3
                ax.axhline(dmax, color=st["color"], ls=":", lw=0.8)
            ax.axhline(0.0, color="0.5", lw=0.5, zorder=0)
            ax.set_xscale("log")
            ax.set_xlim(1e15, 1e21)
            if col == 0:
                ax.set_ylabel(r"$\Delta M_{\max}$ [$10^{-3}\,M_\odot$]")
            if row == 1:
                ax.set_xlabel(r"$\xi$ [G]")
            panel_label(ax, row * 3 + col, f"{model}, {TOPO_LABEL[topo]}")
    handles = [Line2D([], [], color=B0_STYLE[b]["color"], marker=B0_STYLE[b]["marker"],
                      label=b0_label(b)) for b in B0S]
    handles.append(Line2D([], [], color="0.3", ls=":", lw=0.8, label=r"Maxwell ($\xi\to\infty$)"))
    top_legend(fig, handles, ncol=4, top=0.93, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig1_dmmax_vs_xi")


def fig_amm_shifts(sl, sl_amm, out: Path):
    """Desvio causado pelos momentos anômalos em M_max e R_1.4 ao longo de xi."""
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 4.0), sharex=True)
    for col, model in enumerate(MODELS):
        for topo, ls in [("anisotropica", "-"), ("isotropica", "--")]:
            for b0 in B0S:
                off, on = sl[(model, topo, b0)], sl_amm[(model, topo, b0)]
                xi = np.array([x for x, _ in off["log"]])
                st = B0_STYLE[b0]
                for row, (field, scale) in enumerate([("m_max_Msun", 1e3), ("r14_km", 1.0)]):
                    d = np.array([fnum(b[field]) - fnum(a[field])
                                  for (_, a), (_, b) in zip(off["log"], on["log"])]) * scale
                    axes[row, col].plot(xi, d, color=st["color"], ls=ls, marker=st["marker"],
                                        markevery=6, mfc="white" if ls == "--" else st["color"])
        for row in range(2):
            ax = axes[row, col]
            ax.axhline(0.0, color="0.5", lw=0.5, zorder=0)
            ax.set_xscale("log")
            ax.set_xlim(1e15, 1e21)
            panel_label(ax, row * 3 + col, model)
        axes[1, col].set_xlabel(r"$\xi$ [G]")
    axes[0, 0].set_ylabel(r"$\Delta M_{\max}^{\rm AMM}$ [$10^{-3}\,M_\odot$]")
    axes[1, 0].set_ylabel(r"$\Delta R_{1.4}^{\rm AMM}$ [km]")
    handles = [Line2D([], [], color=B0_STYLE[b]["color"], marker=B0_STYLE[b]["marker"],
                      label=b0_label(b)) for b in B0S]
    handles += [Line2D([], [], color="0.3", ls="-", label=r"$P_\perp$"),
                Line2D([], [], color="0.3", ls="--", label=r"$(P_\parallel+2P_\perp)/3$")]
    top_legend(fig, handles, ncol=5, top=0.93, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig2_amm_shift_vs_xi")


def read_stars(path: Path):
    """Ramo estável (até M_max) de uma sequência *_stars.txt: (R, M)."""
    rows = np.loadtxt(path, comments="#")
    rows = rows[np.argsort(rows[:, 3])]  # pressão central crescente
    k = int(np.argmax(rows[:, 0]))
    stable = rows[: k + 1]
    stable = stable[stable[:, 1] < 25.0]
    return stable[:, 1], stable[:, 0]


def fig_mass_radius(root: Path, out: Path, b0="1.00e18", topo="anisotropica"):
    cases = [("sem_tensao", "0.45", "No field stress"),
             ("maxwell", "#000000", "Maxwell"),
             ("xi_1.00e17", "#D55E00", r"Log, $\xi=10^{17}$ G")]
    fig, axes = plt.subplots(1, 3, figsize=(PAGE_W, 2.5), sharey=True)
    for col, model in enumerate(MODELS):
        ax = axes[col]
        for amm, ls in [("", "-"), ("_amm", "--")]:
            base = root / f"study_log{amm}" / model / topo / f"B0_{b0}"
            for name, color, _ in cases:
                path = base / f"{name}_stars.txt"
                if path.exists():
                    r, m = read_stars(path)
                    ax.plot(r, m, color=color, ls=ls, lw=1.0)
        ax.set_xlim(10.5, 15.5)
        ax.set_ylim(0.8, 2.15)
        ax.set_xlabel(r"$R$ [km]")
        panel_label(ax, col, model)
    axes[0].set_ylabel(r"$M$ [$M_\odot$]")
    handles = [Line2D([], [], color=c, label=lab) for _, c, lab in cases]
    handles += [Line2D([], [], color="0.3", ls="-", label="without AMM"),
                Line2D([], [], color="0.3", ls="--", label="with AMM")]
    top_legend(fig, handles, ncol=5, top=0.89, w_pad=0.4)
    save(fig, out, "fig3_mass_radius_B0_1e18")


# ------------------------------------------------------------------ study b


def load_study_b(path: Path, topology="perp"):
    """{(model, nlem, profile): [rows ordenadas por B > 0 até a 1a EoS anômala]}"""
    groups = defaultdict(list)
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if row["topology"] == topology and fnum(row["B_G"]) > 0:
            groups[(row["model"], row["nlem"], row["profile"])].append(row)
    valid = {}
    for key, rows in groups.items():
        rows.sort(key=lambda r: fnum(r["B_G"]))
        kept = []
        for r in rows:
            if r["anomalous"] == "true" or not math.isfinite(fnum(r["m_max_Msun"])):
                break
            kept.append(r)
        valid[key] = kept
    return valid


def fig_b_scan(sb, out: Path):
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 4.0), sharex=True)
    for col, model in enumerate(MODELS):
        for nlem, st in NLEM_STYLE.items():
            for profile, ls in [("constante", "-"), ("bdd", "--")]:
                rows = sb.get((model, nlem, profile), [])
                b = np.array([fnum(r["B_G"]) for r in rows])
                for row, field in enumerate(["m_max_Msun", "r14_km"]):
                    y = np.array([fnum(r[field]) for r in rows])
                    axes[row, col].plot(b, y, color=st["color"], ls=ls)
                    if len(b):  # fim da EoS válida
                        axes[row, col].plot(b[-1], y[-1], color=st["color"], marker="x", ms=4)
        for row in range(2):
            ax = axes[row, col]
            ax.set_xscale("log")
            ax.set_xlim(1e16, 3e19)
            panel_label(ax, row * 3 + col, model)
        axes[1, col].set_xlabel(r"$B$ (constant) or $B_0$ (BDD) [G]")
    axes[0, 0].set_ylabel(r"$M_{\max}$ [$M_\odot$]")
    axes[1, 0].set_ylabel(r"$R_{1.4}$ [km]")
    handles = [Line2D([], [], color=s["color"], label=s["label"]) for s in NLEM_STYLE.values()]
    handles += [Line2D([], [], color="0.3", ls="-", label="constant $B$"),
                Line2D([], [], color="0.3", ls="--", label="BDD profile"),
                Line2D([], [], color="0.3", ls="none", marker="x", label="last valid EoS")]
    top_legend(fig, handles, ncol=4, top=0.89, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig4_stars_vs_B")


def fig_b_amm(sb, sb_amm, out: Path):
    """Desvio do AMM em função de B, só onde as duas EoS são válidas."""
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 4.0), sharex=True)
    for col, model in enumerate(MODELS):
        for nlem in ["maxwell", "log:1.00e17"]:
            st = NLEM_STYLE[nlem]
            for profile, ls in [("constante", "-"), ("bdd", "--")]:
                off = {r["B_G"]: r for r in sb.get((model, nlem, profile), [])}
                on = sb_amm.get((model, nlem, profile), [])
                pairs = [(off[r["B_G"]], r) for r in on if r["B_G"] in off]
                b = np.array([fnum(a["B_G"]) for a, _ in pairs])
                for row, (field, scale) in enumerate([("m_max_Msun", 1e3), ("r14_km", 1.0)]):
                    d = np.array([fnum(bb[field]) - fnum(a[field]) for a, bb in pairs]) * scale
                    axes[row, col].plot(b, d, color=st["color"], ls=ls)
                    if len(b):
                        axes[row, col].plot(b[-1], d[-1], color=st["color"], marker="x", ms=4)
        for row in range(2):
            ax = axes[row, col]
            ax.axhline(0.0, color="0.5", lw=0.5, zorder=0)
            ax.set_xscale("log")
            ax.set_xlim(1e16, 3e19)
            panel_label(ax, row * 3 + col, model)
        axes[1, col].set_xlabel(r"$B$ (constant) or $B_0$ (BDD) [G]")
    axes[0, 0].set_ylabel(r"$\Delta M_{\max}^{\rm AMM}$ [$10^{-3}\,M_\odot$]")
    axes[1, 0].set_ylabel(r"$\Delta R_{1.4}^{\rm AMM}$ [km]")
    handles = [Line2D([], [], color=NLEM_STYLE[n]["color"], label=NLEM_STYLE[n]["label"])
               for n in ["maxwell", "log:1.00e17"]]
    handles += [Line2D([], [], color="0.3", ls="-", label="constant $B$"),
                Line2D([], [], color="0.3", ls="--", label="BDD profile"),
                Line2D([], [], color="0.3", ls="none", marker="x", label="last valid EoS pair")]
    top_legend(fig, handles, ncol=5, top=0.93, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig5_amm_shift_vs_B")


# --------------------------------------------------------------- populações

SPECIES = [  # (coluna, rótulo, cor, bárion?)
    (5, "n", "#000000", True),
    (6, "p", "#0072B2", True),
    (3, r"e^-", "#56B4E9", False),
    (4, r"\mu^-", "#009E73", False),
    (7, r"\Lambda", "#D55E00", True),
    (8, r"\Sigma^-", "#E69F00", True),
    (9, r"\Sigma^0", "#8C510A", True),
    (10, r"\Sigma^+", "#999999", True),
    (11, r"\Xi^-", "#CC79A7", True),
    (12, r"\Xi^0", "#7B3294", True),
]


def fig_populations(root: Path, out: Path, b="1.00e18", csi="1.00e17"):
    fig, axes = plt.subplots(1, 3, figsize=(PAGE_W, 2.6), sharey=True)
    present = set()
    for col, model in enumerate(MODELS):
        ax = axes[col]
        for amm, ls in [("", "-"), ("_amm", "--")]:
            path = root / f"nlem_log{amm}" / model / f"B_{b}" / "default" / f"csi_{csi}" / "eos.dat"
            d = np.loadtxt(path, comments="#")
            d = d[d[:, 0] > 0]
            nb = sum(d[:, c] for c, *_, baryon in SPECIES if baryon)
            for c, lab, color, _ in SPECIES:
                y = d[:, c] / nb
                if np.nanmax(y) < 1e-3:
                    continue
                ax.plot(d[:, 0], y, color=color, ls=ls, lw=0.9)
                present.add(c)
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 1.2)
        ax.set_xlim(0, 8)
        ax.set_xlabel(r"$n_B/n_0$")
        panel_label(ax, col, model)
    axes[0].set_ylabel(r"$Y_i = n_i/n_B$")
    handles = [Line2D([], [], color=color, label=f"${lab}$")
               for c, lab, color, _ in SPECIES if c in present]
    handles += [Line2D([], [], color="0.3", ls="-", label="without AMM"),
                Line2D([], [], color="0.3", ls="--", label="with AMM")]
    top_legend(fig, handles, ncol=6, top=0.80, w_pad=0.4)
    save(fig, out, "fig6_populations_B1e18_xi1e17")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path("results/paper"))
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()

    sl = load_study_log(args.results / "study_log" / "summary.csv")
    sl_amm = load_study_log(args.results / "study_log_amm" / "summary.csv")
    sb = load_study_b(args.results / "study_b" / "stars.csv")
    sb_amm = load_study_b(args.results / "study_b_amm" / "stars.csv")

    fig_mmax_vs_xi(sl, args.out)
    fig_amm_shifts(sl, sl_amm, args.out)
    fig_mass_radius(args.output, args.out)
    fig_b_scan(sb, args.out)
    fig_b_amm(sb, sb_amm, args.out)
    fig_populations(args.output, args.out)


if __name__ == "__main__":
    main()
