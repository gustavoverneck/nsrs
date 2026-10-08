"""Permeabilidade, magnetização, razões campo/total e estabilidade magnética
em função da densidade e do raio, perfil constante x BDD (study profile).

Lê <dir>/eos.csv e <dir>/radial.csv de `nsrs study profile` e gera, para cada B0:

  profiles_density_B0_<B0>.{pdf,png}  em função de n_B/n0
  profiles_radius_B0_<B0>.{pdf,png}   ao longo do raio da estrela de massa máxima

Linhas: cor pela eletrodinâmica, contínua = BDD, tracejada = campo constante.
Painéis (de cima para baixo):
  H/B = 1/mu            inverso da permeabilidade, (H/B)_vac - 4 pi M/B
  4 pi M/B              magnetização da matéria (Landau)
  f_vac e c_m           parte do vácuo (cores) e da matéria (cinza) de
                        s = dH/dB = f_vac - c_m; instável (domínios) onde c_m > f_vac
  eps_campo/eps         fração da energia no campo
  P_perp,campo/P_perp   fração da pressão (topologia anisotrópica)

Gera também unstable_fraction.csv: fração do volume e da massa de cada estrela
(M_max e 1.4 M_sun) em que s < 0 com os dois passos da derivada.

Uso: python plot_scripts/magnetic_profiles.py [results/study_profile] [--out results/study_profile]
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
NLEM_STYLE = {
    "maxwell": dict(color="#000000", label="Maxwell"),
    "log:1.00e16": dict(color="#0072B2", label=r"Log, $\xi=10^{16}$ G"),
    "log:1.00e17": dict(color="#D55E00", label=r"Log, $\xi=10^{17}$ G"),
    "log:1.00e18": dict(color="#009E73", label=r"Log, $\xi=10^{18}$ G"),
}
PROFILE_STYLE = {"bdd": dict(ls="-", label="BDD"), "constante": dict(ls="--", label="constante")}
PANELS = [
    ("h_over_b", r"$H/B=\mu^{-1}$", "symlog"),
    ("m_over_b", r"$4\pi\mathcal{M}/B$", "symlog"),
    ("f_vac", r"$f_{\rm vac}$ (cores), $c_m$ (cinza)", "symlog"),
    ("eps_ratio", r"$\varepsilon_{\rm campo}/\varepsilon$", "log"),
    ("p_ratio", r"$P_{\perp,{\rm campo}}/P_\perp$", "symlog"),
]


def load_eos(path: Path):
    """{(model, b0, profile, nlem): dict de arrays ordenados por n}"""
    raw = defaultdict(list)
    for r in csv.DictReader(open(path, encoding="utf-8")):
        raw[(r["model"], r["B0_G"], r["profile"], r["nlem"])].append(r)
    data = {}
    for key, rows in raw.items():
        rows.sort(key=lambda r: float(r["n_over_n0"]))
        f = lambda name: np.array([float(r[name]) for r in rows])  # noqa: E731
        eps, p = f("eps_MeV_fm3"), f("p_perp_MeV_fm3")
        data[key] = dict(
            n=f("n_over_n0"),
            h_over_b=f("h_over_b_vac") - f("four_pi_M_over_B"),
            m_over_b=f("four_pi_M_over_B"),
            s=f("s"),
            f_vac=f("f_vac"),
            c_m=f("c_matter"),
            eps_ratio=f("eps_field_MeV_fm3") / eps,
            p_ratio=f("p_perp_field_MeV_fm3") / p,
            unstable=np.array([r["magnetically_unstable"] == "true" for r in rows]),
        )
    return data


def load_radial(path: Path):
    """{(model, b0, profile, nlem, star): dict(M, R, r, m, n)}"""
    raw = defaultdict(list)
    for r in csv.DictReader(open(path, encoding="utf-8")):
        raw[(r["model"], r["B0_G"], r["profile"], r["nlem"], r["star"])].append(r)
    out = {}
    for key, rows in raw.items():
        f = lambda name: np.array([float(r[name]) for r in rows])  # noqa: E731
        out[key] = dict(M=float(rows[0]["M_Msun"]), R=float(rows[0]["R_km"]), r=f("r_km"), m=f("m_Msun"), n=f("n_over_n0"))
    return out


def along_radius(eos, star):
    """Grandezas da EoS na linha de n_B mais próxima de cada ponto do raio
    (NaN na crosta)."""
    idx = np.searchsorted(eos["n"], star["n"]).clip(1, len(eos["n"]) - 1)
    left = np.abs(star["n"] - eos["n"][idx - 1]) < np.abs(star["n"] - eos["n"][idx])
    idx = np.where(left, idx - 1, idx)
    core = np.isfinite(star["n"])
    out = {}
    for name in [p[0] for p in PANELS] + ["c_m"]:
        out[name] = np.where(core, eos[name][idx], np.nan)
    out["unstable"] = core & eos["unstable"][idx]
    return out


def setup_axis(ax, scale, name):
    if scale == "symlog":
        ax.set_yscale("symlog", linthresh=1e-3, linscale=0.6)
    elif scale == "log":
        ax.set_yscale("log")
    if name in ("f_vac", "h_over_b", "m_over_b", "p_ratio"):
        ax.axhline(0.0, color="0.5", lw=0.5, zorder=0)
    if name in ("f_vac", "h_over_b"):
        ax.axhline(1.0, color="0.5", lw=0.5, ls=":", zorder=0)


def figure(data, radial, b0, out: Path, by_radius: bool):
    nlems = [k for k in NLEM_STYLE if any(key[3] == k for key in data)]
    fig, axes = plt.subplots(len(PANELS), len(MODELS), figsize=(PAGE_W, 8.6), sharex="col", sharey="row")
    for col, model in enumerate(MODELS):
        for row, (name, ylabel, scale) in enumerate(PANELS):
            ax = axes[row, col]
            setup_axis(ax, scale, name)
            for nlem in nlems:
                for profile, pstyle in PROFILE_STYLE.items():
                    eos = data.get((model, b0, profile, nlem))
                    if eos is None:
                        continue
                    if by_radius:
                        star = radial.get((model, b0, profile, nlem, "max"))
                        if star is None:
                            continue
                        x, values = star["r"], along_radius(eos, star)
                    else:
                        x, values = eos["n"], eos
                    # A matéria não depende da NLEM (acoplamento mínimo):
                    # magnetização e c_m uma vez por perfil, em cinza.
                    matter_once = nlem == nlems[0]
                    if name == "m_over_b":
                        if matter_once:
                            ax.plot(x, values[name], color="0.35", ls=pstyle["ls"], lw=0.8)
                        continue
                    if name == "f_vac" and matter_once:
                        ax.plot(x, values["c_m"], color="0.6", ls=pstyle["ls"], lw=0.6, zorder=1)
                    ax.plot(x, values[name], color=NLEM_STYLE[nlem]["color"], ls=pstyle["ls"], lw=0.9, zorder=2)
            if row == 0:
                ax.set_title(model)
            if col == 0:
                ax.set_ylabel(ylabel)
        axes[-1, col].set_xlabel(r"$r$ [km]" if by_radius else r"$n_B/n_0$")
        if not by_radius:
            axes[-1, col].set_xlim(0.0, None)
    handles = [Line2D([], [], color=NLEM_STYLE[k]["color"], label=NLEM_STYLE[k]["label"]) for k in nlems]
    handles += [Line2D([], [], color="0.3", ls=s["ls"], label=s["label"]) for s in PROFILE_STYLE.values()]
    exp = int(round(math.log10(float(b0))))
    where = "estrela de massa máxima" if by_radius else "EoS"
    fig.legend(handles=handles, loc="upper center", ncol=6, fontsize=7, bbox_to_anchor=(0.5, 1.0),
               title=rf"$B_0=10^{{{exp}}}$ G ({where})", title_fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=0.4, w_pad=0.6)
    stem = f"profiles_{'radius' if by_radius else 'density'}_B0_1e{exp}"
    for ext, kw in [("pdf", {}), ("png", {"dpi": 200})]:
        fig.savefig(out / f"{stem}.{ext}", **kw)
    plt.close(fig)
    print(f"figura: {out / stem}.pdf")


def unstable_fractions(data, radial, out: Path):
    """Fração do volume e da massa com s < 0, e faixa de raio instável."""
    rows = []
    for key, star in sorted(radial.items()):
        model, b0, profile, nlem, which = key
        eos = data.get((model, b0, profile, nlem))
        if eos is None:
            continue
        flag = along_radius(eos, star)["unstable"]
        r, m = star["r"], star["m"]
        mid = flag[1:] & flag[:-1]
        dvol = (r[1:] ** 3 - r[:-1] ** 3)
        vol = float(np.sum(dvol[mid]) / star["R"] ** 3)
        mass = float(np.sum(np.diff(m)[mid]) / star["M"])
        r_unst = r[flag]
        rows.append([model, f"{float(b0):.0e}", profile, nlem, which, f"{star['M']:.4f}", f"{star['R']:.3f}",
                     f"{vol:.4f}", f"{mass:.4f}",
                     f"{r_unst.min():.3f}" if r_unst.size else "", f"{r_unst.max():.3f}" if r_unst.size else ""])
    path = out / "unstable_fraction.csv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow(["model", "B0_G", "profile", "nlem", "star", "M_Msun", "R_km", "volume_fraction_unstable",
                    "mass_fraction_unstable", "r_unstable_min_km", "r_unstable_max_km"])
        w.writerows(rows)
    print(f"tabela: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dir", nargs="?", type=Path, default=Path("results/study_profile"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    out = args.out or args.dir
    out.mkdir(parents=True, exist_ok=True)
    data = load_eos(args.dir / "eos.csv")
    radial = load_radial(args.dir / "radial.csv")
    for b0 in sorted({key[1] for key in data}, key=float):
        figure(data, radial, b0, out, by_radius=False)
        figure(data, radial, b0, out, by_radius=True)
    unstable_fractions(data, radial, out)


if __name__ == "__main__":
    main()
