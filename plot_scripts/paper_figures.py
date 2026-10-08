"""Figuras de publicação (Physical Review D) para a NLEM log, sem e com AMM.

Só as figuras que sustentam a análise do artigo:
  fig1  tensões do campo Log em função de B/xi (analítica, explica as demais)
  fig2  Delta M_max(xi) em relação ao caso sem tensão (study log)
  fig3  M_max, R_1.4 e Lambda_1.4 em função de B, Maxwell x Log (study b)
  fig4  B/xi e P_campo/P_c no centro da estrela de massa máxima (study b):
        liga a fig3 aos limiares da fig1 e mede a validade do TOV esférico
  fig5  curvas M-R em B0 = 1e18 G, sem e com AMM (study log --save-eos)
  fig6  desvio do AMM em função de B (study b, sem e com --amm)
  table_endpoints.csv  por que cada curva da fig3 termina (study b)
O efeito do AMM ao longo de xi é constante e as populações quase não mudam com
o AMM; esses números vão para o texto, não para figuras.

Lê apenas resultados já calculados:
  results/study_log/summary.csv             (study log)
  output/study_log{,_amm}/.../*_stars.txt   (curvas M-R, gravadas com --save-eos)
  results/study_b{,_amm}/stars.csv          (study b)

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
from matplotlib.patches import Patch  # noqa: E402

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
# PSR J0740+6620 (Fonseca et al., ApJL 915, L12 (2021)): 2.08 +- 0.07 M_sun.
J0740 = (2.08, 0.07)
# Figuras de versões anteriores deste script, apagadas de --out para não confundir.
# GW170817 (Abbott et al., PRL 121, 161101 (2018)): Lambda_1.4 = 190 +390 -120 (90%).
GW170817_L14 = (70.0, 580.0)
# Perfil BDD do campo local usado na energia e nas tensões (core::magnetic).
B_SURF, BDD_BETA, BDD_GAMMA = 1e15, 0.01, 3.0
# Figuras de versões anteriores deste script, apagadas de --out para não confundir.
OBSOLETE = ["fig1_dmmax_vs_xi", "fig2_amm_shift_vs_xi", "fig3_mass_radius_B0_1e18",
            "fig4_stars_vs_B", "fig4_mass_radius_B0_1e18", "fig5_amm_shift_vs_B",
            "fig6_populations_B1e18_xi1e17"]
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


def local_field(b0: float, nb_over_n0: float) -> float:
    """Campo local (G) do perfil BDD, também usado pela energia no perfil constante."""
    return B_SURF + b0 * (1.0 - math.exp(-BDD_BETA * max(nb_over_n0, 0.0) ** BDD_GAMMA))


def nlem_xi(nlem: str) -> float:
    """xi (G) de um rótulo 'log:1.00e17'; NaN para Maxwell."""
    return float(nlem.split(":")[1]) if nlem.startswith("log:") else math.nan


def j0740_band(ax):
    m, dm = J0740
    ax.axhspan(m - dm, m + dm, color="0.85", lw=0, zorder=0)


# ------------------------------------------------------------ campo Log


def fig_field_stress(out: Path):
    """eps, P_perp e f_vac da Log normalizados por Maxwell, em função de B/xi.

    eps = xi^2 ln(1+x)/(4 pi), x = B^2/(2 xi^2) (core::magnetic::magnetic_stress);
    P_perp = H B - eps, P_par = -eps; f_vac = (1-x)/(1+x)^2 (NlemModel::vacuum_curvature).
    """
    u = np.logspace(-1, 2, 600)  # B/xi
    x = u**2 / 2
    eps = np.log1p(x) / x
    p_perp = 2 / (1 + x) - eps
    f_vac = (1 - x) / (1 + x) ** 2
    fig, ax = plt.subplots(figsize=(3.375, 2.5))
    ax.axhline(1.0, color="0.5", lw=0.6, ls=":", zorder=0)
    ax.axhline(0.0, color="0.5", lw=0.5, zorder=0)
    ax.plot(u, eps, color="#0072B2", label=r"$\varepsilon_B/\varepsilon_B^{\rm Max}$")
    ax.plot(u, p_perp, color="#D55E00", label=r"$P_\perp/\varepsilon_B^{\rm Max}$")
    ax.plot(u, -eps, color="#009E73", label=r"$P_\parallel/\varepsilon_B^{\rm Max}$")
    ax.plot(u, f_vac, color="#000000", ls="--", label=r"$f_{\rm vac}$")
    # f_vac = 0 em x = 1; P_perp = 0 em x = 3.92 (B/xi = 2.80).
    for xc in (np.sqrt(2.0), 2.8006):
        ax.axvline(xc, color="0.6", lw=0.5, ls="-.", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(u[0], u[-1])
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel(r"$B/\xi$")
    ax.legend(loc="upper right")
    save(fig, out, "fig1_log_field_stress")


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
    save(fig, out, "fig2_dmmax_vs_xi")


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
        j0740_band(ax)
        ax.set_xlim(10.5, 15.5)
        ax.set_ylim(0.8, 2.2)
        ax.set_xlabel(r"$R$ [km]")
        panel_label(ax, col, model)
    axes[0].set_ylabel(r"$M$ [$M_\odot$]")
    handles = [Line2D([], [], color=c, label=lab) for _, c, lab in cases]
    handles += [Line2D([], [], color="0.3", ls="-", label="without AMM"),
                Line2D([], [], color="0.3", ls="--", label="with AMM"),
                Patch(color="0.85", label="PSR J0740+6620")]
    top_legend(fig, handles, ncol=6, top=0.89, w_pad=0.4)
    save(fig, out, "fig5_mass_radius_B0_1e18")


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
    fig, axes = plt.subplots(3, 3, figsize=(PAGE_W, 5.6), sharex=True)
    for col, model in enumerate(MODELS):
        for nlem, st in NLEM_STYLE.items():
            for profile, ls in [("constante", "-"), ("bdd", "--")]:
                rows = sb.get((model, nlem, profile), [])
                b = np.array([fnum(r["B_G"]) for r in rows])
                for row, field in enumerate(["m_max_Msun", "r14_km", "lambda14"]):
                    y = np.array([fnum(r[field]) for r in rows])
                    axes[row, col].plot(b, y, color=st["color"], ls=ls)
                    if len(b):  # fim da EoS válida
                        axes[row, col].plot(b[-1], y[-1], color=st["color"], marker="x", ms=4)
        for row in range(3):
            ax = axes[row, col]
            ax.set_xscale("log")
            ax.set_xlim(1e16, 3e19)
            panel_label(ax, row * 3 + col, model)
        j0740_band(axes[0, col])
        axes[2, col].axhspan(*GW170817_L14, color="#cfe3f2", lw=0, zorder=0)
        axes[2, col].set_yscale("log")
        axes[2, col].set_xlabel(r"$B$ (constant) or $B_0$ (BDD) [G]")
    axes[0, 0].set_ylabel(r"$M_{\max}$ [$M_\odot$]")
    axes[1, 0].set_ylabel(r"$R_{1.4}$ [km]")
    axes[2, 0].set_ylabel(r"$\Lambda_{1.4}$")
    handles = [Line2D([], [], color=s["color"], label=s["label"]) for s in NLEM_STYLE.values()]
    handles += [Line2D([], [], color="0.3", ls="-", label="constant $B$"),
                Line2D([], [], color="0.3", ls="--", label="BDD profile"),
                Line2D([], [], color="0.3", ls="none", marker="x", label="last valid EoS"),
                Patch(color="0.85", label="PSR J0740+6620"),
                Patch(color="#cfe3f2", label="GW170817 (90%)")]
    top_legend(fig, handles, ncol=4, top=0.89, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig3_stars_vs_B")


def fig_center_field(sb, out: Path):
    """No centro da estrela de massa máxima: B_c/xi (onde as estrelas caem na
    fig1) e P_campo/P_c (validade do TOV esférico; ~10% como referência)."""
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 4.0), sharex=True)
    for col, model in enumerate(MODELS):
        for nlem, st in NLEM_STYLE.items():
            xi = nlem_xi(nlem)
            for profile, ls in [("constante", "-"), ("bdd", "--")]:
                rows = sb.get((model, nlem, profile), [])
                b = np.array([fnum(r["B_G"]) for r in rows])
                bc = np.array([local_field(fnum(r["B_G"]), fnum(r["nc_over_n0"])) for r in rows])
                ratio = np.array([fnum(r.get("p_field_c_MeV_fm3")) / fnum(r.get("p_c_MeV_fm3"))
                                  for r in rows])
                if math.isfinite(xi):
                    axes[0, col].plot(b, bc / xi, color=st["color"], ls=ls)
                axes[1, col].plot(b, ratio, color=st["color"], ls=ls)
        axes[0, col].axhline(math.sqrt(2.0), color="0.4", lw=0.6, ls="-.")
        axes[0, col].axhline(2.8006, color="0.4", lw=0.6, ls=":")
        axes[1, col].axhline(0.1, color="0.4", lw=0.6, ls=":")
        axes[1, col].axhline(0.0, color="0.5", lw=0.5, zorder=0)
        axes[0, col].set_yscale("log")
        for row in range(2):
            ax = axes[row, col]
            ax.set_xscale("log")
            ax.set_xlim(1e16, 3e19)
            panel_label(ax, row * 3 + col, model)
        axes[1, col].set_xlabel(r"$B$ (constant) or $B_0$ (BDD) [G]")
    axes[0, 0].set_ylabel(r"$B_c/\xi$")
    axes[1, 0].set_ylabel(r"$P_{\rm field}/P$ at center")
    handles = [Line2D([], [], color=s["color"], label=s["label"]) for s in NLEM_STYLE.values()]
    handles += [Line2D([], [], color="0.3", ls="-", label="constant $B$"),
                Line2D([], [], color="0.3", ls="--", label="BDD profile"),
                Line2D([], [], color="0.4", ls="-.", lw=0.6, label=r"$f_{\rm vac}=0$"),
                Line2D([], [], color="0.4", ls=":", lw=0.6, label=r"$P_\perp=0$; 10% of $P$")]
    top_legend(fig, handles, ncol=4, top=0.88, h_pad=0.4, w_pad=0.6)
    save(fig, out, "fig4_center_field")


def table_endpoints(path: Path, stability: Path, out: Path, topology="perp"):
    """Por que cada curva (modelo, NLEM, perfil) termina: primeiro B inválido,
    seu término, B/xi local na maior densidade alcançada e qual critério físico
    vale ali. 'solver' indica que nenhum critério físico explica o fim. Uma
    curva Log que termina no mesmo B que a Maxwell do mesmo modelo e perfil
    recebe 'as Maxwell': o fim vem da matéria ou do solver, não da Log, e os
    critérios f_vac/P_perp ficam só como sinalização."""
    groups = defaultdict(list)
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if row["topology"] == topology and fnum(row["B_G"]) > 0:
            groups[(row["model"], row["nlem"], row["profile"])].append(row)
    unstable = set()
    if stability.exists():
        for row in csv.DictReader(open(stability, encoding="utf-8")):
            if row["magnetically_unstable"] == "true":
                unstable.add((row["model"], row["nlem"], row["B_G"]))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "table_endpoints.csv", "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow(["model", "nlem", "profile", "B_last_valid_G", "B_first_invalid_G", "termination",
                    "n_max_over_n0", "B_local_over_xi", "f_vac_negative", "pperp_negative_rows",
                    "magnetically_unstable", "ends_with_maxwell", "cause"])

        def first_bad(rows):
            return next((i for i, r in enumerate(rows) if r["anomalous"] == "true"
                         or not math.isfinite(fnum(r["m_max_Msun"]))), None)

        for rows in groups.values():
            rows.sort(key=lambda r: fnum(r["B_G"]))
        end_maxwell = {}
        for (model, nlem, profile), rows in groups.items():
            if nlem == "maxwell":
                i = first_bad(rows)
                end_maxwell[(model, profile)] = rows[i]["B_G"] if i is not None else None
        for (model, nlem, profile), rows in sorted(groups.items()):
            bad = first_bad(rows)
            if bad is None:
                w.writerow([model, nlem, profile, rows[-1]["B_G"], "", "", "", "", "", "", "", "", "none"])
                continue
            r = rows[bad]
            xi = nlem_xi(nlem)
            u = local_field(fnum(r["B_G"]), fnum(r["n_max_over_n0"])) / xi
            fvac_neg = math.isfinite(u) and u > math.sqrt(2.0)
            pneg = int(fnum(r["pperp_negative_rows"]) or 0)
            mag = profile == "constante" and (model, nlem, r["B_G"]) in unstable
            same = nlem != "maxwell" and end_maxwell.get((model, profile)) == r["B_G"]
            cause = ("as Maxwell" if same else "f_vac<0" if fvac_neg else "P_perp<0" if pneg > 0
                     else "magnetic instability" if mag else "solver")
            w.writerow([model, nlem, profile, rows[bad - 1]["B_G"] if bad else "", r["B_G"],
                        r["termination"], r["n_max_over_n0"], f"{u:.3g}" if math.isfinite(u) else "",
                        fvac_neg, pneg, mag, same if nlem != "maxwell" else "", cause])
    print(f"tabela: {out / 'table_endpoints.csv'}")


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
    save(fig, out, "fig6_amm_shift_vs_B")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path("results/paper"))
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()

    for name in OBSOLETE:
        for ext in ("pdf", "png"):
            (args.out / f"{name}.{ext}").unlink(missing_ok=True)

    sl = load_study_log(args.results / "study_log" / "summary.csv")
    sb = load_study_b(args.results / "study_b" / "stars.csv")
    sb_amm = load_study_b(args.results / "study_b_amm" / "stars.csv")

    fig_field_stress(args.out)
    fig_mmax_vs_xi(sl, args.out)
    fig_b_scan(sb, args.out)
    fig_center_field(sb, args.out)
    fig_mass_radius(args.output, args.out)
    fig_b_amm(sb, sb_amm, args.out)
    table_endpoints(args.results / "study_b" / "stars.csv",
                    args.results / "study_b" / "stability.csv", args.out)


if __name__ == "__main__":
    main()
