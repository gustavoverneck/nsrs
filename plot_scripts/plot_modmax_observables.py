#!/usr/bin/env python3
"""
High-quality scientific plots for ModMax exterior observables.

Reads results/modmax/observables_selected_csi.csv and produces publication-ready
figures grouped by (model, B, topology).
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np


@dataclass
class Row:
    model: str
    b_label: str
    b_value: float
    topology: str
    csi: float
    m_star_msun: float
    r_star_km: float
    f_surf: float
    b_pole: float
    h_eff: float
    n_perp: float
    n_par: float
    delta_phi: float
    stokes_i: float
    stokes_q: float
    stokes_u: float
    stokes_v: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ModMax exterior observables at high scientific standard"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/modmax/observables_selected_csi.csv"),
        help="CSV with observables",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/modmax/observables_plots"),
        help="Directory to store plots",
    )
    parser.add_argument("--dpi", type=int, default=600, help="Figure DPI")
    return parser.parse_args()


def load_rows(path: Path) -> List[Row]:
    rows: List[Row] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                Row(
                    model=r["model"],
                    b_label=r["b_label"],
                    b_value=float(r["b_value"]),
                    topology=r["topology"],
                    csi=float(r["csi"]),
                    m_star_msun=float(r["m_star_msun"]),
                    r_star_km=float(r["r_star_km"]),
                    f_surf=float(r["f_surf"]),
                    b_pole=float(r["b_pole"]),
                    h_eff=float(r["h_eff"]),
                    n_perp=float(r["n_perp"]),
                    n_par=float(r["n_par"]),
                    delta_phi=float(r["delta_phi"]),
                    stokes_i=float(r["stokes_i"]),
                    stokes_q=float(r["stokes_q"]),
                    stokes_u=float(r["stokes_u"]),
                    stokes_v=float(r["stokes_v"]),
                )
            )
    return rows


def group_rows(rows: List[Row]) -> Dict[Tuple[str, str, str], List[Row]]:
    grouped: Dict[Tuple[str, str, str], List[Row]] = {}
    for r in rows:
        key = (r.model, r.b_label, r.topology)
        grouped.setdefault(key, []).append(r)
    for k in grouped:
        grouped[k].sort(key=lambda x: x.csi)
    return grouped


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.size": 11,
            "font.family": "serif",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": "--",
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )


CSI_LABEL = r"$\xi$ (ModMax)"


def save_fig(fig: Figure, out_dir: Path, name: str, dpi: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # fig.savefig(out_dir / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_indices(rows: List[Row], out_dir: Path, tag: str, dpi: int) -> None:
    csi = np.array([r.csi for r in rows])
    n_perp = np.array([r.n_perp for r in rows])
    n_par = np.array([r.n_par for r in rows])

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(csi, n_perp, "-", label=r"$n_{\perp}$", linewidth=2)
    ax.plot(
        csi,
        n_par,
        "--",
        label=r"$n_{\parallel}$",
        linewidth=1.5,
    )
    ax.axhline(y=1.0, color="gray", linestyle=":", label="Maxwell ($n=1$)")
    ax.set_xscale("log")
    ax.set_xlabel(CSI_LABEL)
    ax.set_ylabel("Refractive index")
    ax.legend()
    ax.set_title(tag)
    save_fig(fig, out_dir, f"indices_{tag}", dpi)


def plot_phase(rows: List[Row], out_dir: Path, tag: str, dpi: int) -> None:
    csi = np.array([r.csi for r in rows])
    delta_phi = np.array([r.delta_phi for r in rows])

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(csi, delta_phi, "-", color="#2C7BB6")
    ax.axhline(y=0.0, color="gray", linestyle=":", label="Maxwell")
    ax.set_xscale("log")
    ax.set_xlabel(CSI_LABEL)
    ax.set_ylabel(r"$\Delta\phi$ [rad]")
    ax.set_ylim(-0.1, 0.1)
    ax.legend()
    ax.set_title(tag)
    save_fig(fig, out_dir, f"delta_phi_{tag}", dpi)


def plot_h_field(rows: List[Row], out_dir: Path, tag: str, dpi: int) -> None:
    csi = np.array([r.csi for r in rows])
    h = np.array([r.h_eff for r in rows])
    b_pole = np.array([r.b_pole for r in rows])

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(csi, h, "-", label=r"$H_{\mathrm{eff}}$")
    ax.plot(csi, b_pole, "--", label=r"$B_{\mathrm{pole}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(CSI_LABEL)
    ax.set_ylabel(r"Field strength [G]")
    ax.legend()
    ax.set_title(tag)
    save_fig(fig, out_dir, f"field_{tag}", dpi)


def plot_stokes(rows: List[Row], out_dir: Path, tag: str, dpi: int) -> None:
    csi = np.array([r.csi for r in rows])
    i = np.array([r.stokes_i for r in rows])
    q = np.array([r.stokes_q for r in rows])
    u = np.array([r.stokes_u for r in rows])
    v = np.array([r.stokes_v for r in rows])

    p_lin = np.sqrt(q * q + u * u) / np.maximum(i, 1e-30)
    v_frac = v / np.maximum(i, 1e-30)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(csi, p_lin, "-", label=r"$P_{\mathrm{lin}}$")
    ax.plot(csi, v_frac, "--", label=r"$V/I$")
    ax.set_xscale("log")
    ax.set_xlabel(CSI_LABEL)
    ax.set_ylabel("Polarization fraction")
    ax.set_ylim(-0.1, 1.1)
    ax.legend()
    ax.set_title(tag)
    save_fig(fig, out_dir, f"stokes_{tag}", dpi)


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"Arquivo não encontrado: {args.input}")

    setup_style()
    rows = load_rows(args.input)
    grouped = group_rows(rows)

    for (model, b_label, topology), group in grouped.items():
        tag = f"{model}_B_{b_label}_{topology}"
        plot_indices(group, args.output_dir, tag, args.dpi)
        plot_phase(group, args.output_dir, tag, args.dpi)
        plot_h_field(group, args.output_dir, tag, args.dpi)
        plot_stokes(group, args.output_dir, tag, args.dpi)

    print(f"[OK] Plots saved to {args.output_dir}")


if __name__ == "__main__":
    main()
