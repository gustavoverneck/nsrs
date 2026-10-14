"""Birrefringência do vácuo na eletrodinâmica logarítmica, comparada com a QED.

Índices dos dois modos de um fóton fraco que se propaga perpendicular a um
campo magnético estático B, com x = B^2/(2 xi^2) (mesma convenção de
NlemModel::photon_indices_perp em src/core/physics.rs):

  Log de Gaete e Helayël-Neto: n_par^2 = 1 + 2x,  n_perp^2 = (1 + x)/(1 - x)
  Log sem o termo G^2 (Soleng): n_par = 1,         n_perp^2 = (1 + x)/(1 - x)

n_par e n_perp: campo elétrico da onda paralelo e perpendicular a B. n_perp
diverge em B = sqrt(2) xi (f_vac = 0) e não existe acima. QED (Euler-Heisenberg),
só as assíntotas: Delta n = (alpha/30 pi)(B/B_c)^2 para B << B_c e
(alpha/6 pi)(B/B_c) para B >> B_c, com B_c = 4.414e13 G.

Uso: python plot_scripts/nlem_birefringence.py [--out results/paper]
Gera fig7_birefringence.{pdf,png} e birefringence.csv.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.style.use(HERE / "paper.mplstyle")

PAGE_W = 7.0
ALPHA = 1.0 / 137.035999
B_CRIT = 4.414e13  # G
XIS = [1e15, 1e16, 1e17, 1e18]
COLORS = ["#CC79A7", "#0072B2", "#D55E00", "#009E73"]
MAGNETAR = (1e14, 1e15)  # campos de superfície de magnetares, G


def log_indices(b, xi, g_term=True):
    """(n_par - 1, n_perp - 1, Delta n) da Log, sem cancelamento numérico em
    campo fraco (n - 1 = (n^2 - 1)/(n + 1)); NaN onde n_perp não se propaga."""
    x = np.asarray(b, dtype=float) ** 2 / (2.0 * xi * xi)
    with np.errstate(divide="ignore", invalid="ignore"):
        good = x < 1.0
        n_perp = np.where(good, np.sqrt(np.abs((1.0 + x) / (1.0 - x))), np.nan)
        dperp = np.where(good, 2.0 * x / (1.0 - x) / (n_perp + 1.0), np.nan)
        if g_term:
            n_par = np.sqrt(1.0 + 2.0 * x)
            dpar = 2.0 * x / (n_par + 1.0)
            dn = np.where(good, 2.0 * x * x / (1.0 - x) / (n_perp + n_par), np.nan)
        else:
            dpar = np.zeros_like(x)
            dn = dperp
    return dpar, dperp, dn


def qed_dn(b):
    """Assíntotas de Euler-Heisenberg para Delta n (campo fraco e forte)."""
    r = b / B_CRIT
    return ALPHA / (30.0 * math.pi) * r**2, ALPHA / (6.0 * math.pi) * r


def fig_birefringence(out: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(PAGE_W, 3.2))

    # (a) universal em B/xi.
    u = np.logspace(-3, math.log10(1.4142), 400)
    gh_par, gh_perp, gh_dn = log_indices(u, 1.0)
    so_dn = log_indices(u, 1.0, g_term=False)[2]
    ax1.loglog(u, gh_perp, color="#0072B2", lw=2.0, label=r"$n_\perp-1$")
    ax1.loglog(u, gh_par, color="#D55E00", label=r"$n_\parallel-1$ (GH)")
    ax1.loglog(u, gh_dn, color="k", label=r"$\Delta n$ (GH)")
    ax1.loglog(u, so_dn, color="k", ls="--", label=r"$\Delta n$ (Soleng)")
    ax1.axvline(math.sqrt(2.0), color="0.4", ls="-.", lw=0.7)
    ax1.text(1.25, 1e-6, r"$f_{\rm vac}=0$", rotation=90, ha="right", va="bottom", fontsize=7)
    ax1.set_xlim(1e-3, 3.0)
    ax1.set_ylim(1e-12, 30.0)
    ax1.set_xlabel(r"$B/\xi$")
    ax1.set_ylabel(r"$n-1$, $\Delta n$")
    ax1.legend(loc="upper left")
    ax1.set_title("(a)", loc="left")

    # (b) Delta n em função de B: Log (GH) para vários xi, contra a QED.
    b = np.logspace(12, 19, 600)
    weak, strong = qed_dn(b)
    ax2.loglog(b[b < B_CRIT], weak[b < B_CRIT], color="0.3", lw=1.6, label="QED, $B\\ll B_c$")
    ax2.loglog(b[b > B_CRIT], strong[b > B_CRIT], color="0.3", lw=1.6, ls=":", label="QED, $B\\gg B_c$")
    for xi, c in zip(XIS, COLORS):
        ax2.loglog(b, log_indices(b, xi)[2], color=c, label=rf"Log, $\xi=10^{{{int(round(math.log10(xi)))}}}$ G")
    ax2.loglog(b, log_indices(b, 1e16, g_term=False)[2], color=COLORS[1], ls="--", lw=0.8, label=r"Soleng, $\xi=10^{16}$ G")
    ax2.axvspan(*MAGNETAR, color="0.88", lw=0, zorder=0)
    ax2.axvline(B_CRIT, color="0.5", lw=0.6, ls="-.")
    ax2.text(B_CRIT * 1.15, 3e-16, r"$B_c$", fontsize=7)
    ax2.set_xlim(1e12, 1e19)
    ax2.set_ylim(1e-16, 10.0)
    ax2.set_xlabel(r"$B$ [G]")
    ax2.set_ylabel(r"$\Delta n=n_\perp-n_\parallel$")
    ax2.set_title("(b)", loc="left")
    fig.legend(*ax2.get_legend_handles_labels(), loc="upper center", ncol=4, fontsize=7,
               bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(w_pad=1.2, rect=(0, 0, 1, 0.88))
    out.mkdir(parents=True, exist_ok=True)
    for ext, kw in [("pdf", {}), ("png", {"dpi": 300})]:
        fig.savefig(out / f"fig7_birefringence.{ext}", **kw)
    plt.close(fig)
    print(f"figura: {out / 'fig7_birefringence'}.pdf")


def write_table(out: Path):
    """Delta n na superfície de magnetares e o ponto em que Log = QED."""
    rows = []
    for xi in XIS:
        for bs in [1e13, 1e14, 1e15]:
            n_par, n_perp, dn = (v[0] for v in log_indices(np.array([bs]), xi))
            so = log_indices(np.array([bs]), xi, g_term=False)[2][0]
            weak, strong = qed_dn(bs)
            qed = weak if bs < B_CRIT else strong
            rows.append([f"{xi:.0e}", f"{bs:.0e}", f"{bs * bs / (2 * xi * xi):.3e}",
                         f"{n_par:.3e}", f"{n_perp:.3e}", f"{dn:.3e}",
                         f"{so:.3e}", f"{qed:.3e}"])
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "birefringence.csv", "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow(["xi_G", "B_G", "x", "n_par_minus_1", "n_perp_minus_1", "dn_GH",
                    "dn_Soleng", "dn_QED_asymptotic"])
        w.writerows(rows)
    print(f"tabela: {out / 'birefringence.csv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path("results/paper"))
    args = parser.parse_args()
    fig_birefringence(args.out)
    write_table(args.out)


if __name__ == "__main__":
    main()
