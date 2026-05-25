#!/usr/bin/env python3
"""
Compute ModMax exterior observables from existing output/modmax datasets.

Reads eos.dat files, extracts M and R at maximum mass, and computes
surface/optical observables using the exterior model (Schwarzschild + ModMax).
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np


@dataclass
class DatasetMeta:
    model: str
    b_label: str
    b_value: float
    topology: str
    csi: float
    eos_path: Path


@dataclass
class Observables:
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
        description="Compute ModMax observables from output/modmax"
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("output/modmax"),
        help="Root with ModMax outputs",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/modmax/observables_from_output.csv"),
        help="CSV output path",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=200,
        help="Steps for phase integral",
    )
    parser.add_argument(
        "--r-factor",
        type=float,
        default=10.0,
        help="Outer radius in units of R",
    )
    parser.add_argument(
        "--omega",
        type=float,
        default=1.0,
        help="Photon angular frequency (arbitrary units)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots after computing observables",
    )
    parser.add_argument(
        "--plot-output-dir",
        type=Path,
        default=Path("results/modmax/observables_plots"),
        help="Directory to store plots",
    )
    parser.add_argument(
        "--plot-dpi",
        type=int,
        default=600,
        help="Plot DPI",
    )
    return parser.parse_args()


def discover_datasets(root: Path) -> List[DatasetMeta]:
    metas: List[DatasetMeta] = []
    for eos in root.rglob("eos.dat"):
        meta = parse_metadata(root, eos)
        if meta:
            metas.append(meta)
    return metas


def parse_metadata(root: Path, eos_path: Path) -> Optional[DatasetMeta]:
    try:
        rel = eos_path.relative_to(root)
    except ValueError:
        return None

    parts = list(rel.parts)
    if len(parts) < 3:
        return None

    model = parts[0]
    b_label_raw = parts[1]
    if not b_label_raw.startswith("B_"):
        return None
    b_label = b_label_raw[2:]

    try:
        b_value = float(b_label)
    except ValueError:
        return None

    if parts[2].startswith("csi_"):
        topology = "isotropic"
        csi_part = parts[2]
    else:
        if len(parts) < 4:
            return None
        topology = parts[2]
        csi_part = parts[3]

    if not csi_part.startswith("csi_"):
        return None

    try:
        csi = float(csi_part[4:])
    except ValueError:
        return None

    return DatasetMeta(
        model=model,
        b_label=b_label,
        b_value=b_value,
        topology=topology,
        csi=csi,
        eos_path=eos_path,
    )


def load_eos(path: Path) -> Optional[np.ndarray]:
    try:
        data = np.genfromtxt(path, comments="#")
    except Exception:
        return None

    if data is None or np.size(data) == 0:
        return None
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def max_mass_pair(data: np.ndarray) -> Optional[Tuple[float, float, int]]:
    if data.shape[1] < 2:
        return None
    mass = data[:, -2]
    radius = data[:, -1]
    mask = (
        np.isfinite(mass)
        & np.isfinite(radius)
        & (mass > 0.0)
        & (radius > 0.0)
    )
    if not np.any(mask):
        return None
    idx = np.argmax(mass[mask])
    true_indices = np.nonzero(mask)[0]
    i_max = int(true_indices[idx])
    return float(mass[i_max]), float(radius[i_max]), i_max


class MagnetarContext:
    def __init__(self, mass: float, radius: float, b_surf: float) -> None:
        self.mass = mass
        self.radius = radius
        self.b_surf = b_surf

    def schwarzschild_f(self, r: float) -> float:
        return 1.0 - 2.0 * self.mass / r

    def dipole_field(self, r: float, theta: float) -> float:
        r_ratio = self.radius / r
        return self.b_surf * (r_ratio ** 3) * math.sqrt(1.0 + 3.0 * math.cos(theta) ** 2)


class ModmaxNled:
    def __init__(self, csi: float) -> None:
        self.csi = csi

    def constitutive_h(self, b: float) -> float:
        return b * math.exp(-self.csi)

    def refractive_indices(self, _b: float) -> Tuple[float, float]:
        n = math.exp(self.csi)
        return n, n


def compute_delta_phi(
    ctx: MagnetarContext, nled: ModmaxNled, steps: int, r_factor: float, omega: float
) -> float:
    if steps <= 0:
        return 0.0
    r0 = ctx.radius
    r1 = ctx.radius * max(1.0, r_factor)
    dr = (r1 - r0) / steps

    n_perp, n_par = nled.refractive_indices(0.0)
    delta_phi = 0.0

    for i in range(steps):
        r = r0 + dr * (i + 0.5)
        g_rr = 1.0 / ctx.schwarzschild_f(r)
        ds = math.sqrt(g_rr) * abs(dr)
        delta_phi += omega * (n_par - n_perp) * ds

    return delta_phi


def compute_observables(
    meta: DatasetMeta, data: np.ndarray, steps: int, r_factor: float, omega: float
) -> Optional[Observables]:
    max_pair = max_mass_pair(data)
    if max_pair is None:
        return None
    m_star, r_star, _ = max_pair

    ctx = MagnetarContext(m_star, r_star, meta.b_value)
    nled = ModmaxNled(meta.csi)

    b_pole = ctx.dipole_field(ctx.radius, 0.0)
    h_eff = nled.constitutive_h(b_pole)
    n_perp, n_par = nled.refractive_indices(b_pole)
    f_surf = ctx.schwarzschild_f(ctx.radius)

    delta_phi = compute_delta_phi(ctx, nled, steps, r_factor, omega)

    stokes_i = 1.0
    stokes_q = 1.0
    stokes_u = 0.0
    stokes_v = 0.0

    if abs(delta_phi) > 0.0:
        cos_phi = math.cos(delta_phi)
        sin_phi = math.sin(delta_phi)
        stokes_u, stokes_v = (
            stokes_u * cos_phi - stokes_v * sin_phi,
            stokes_u * sin_phi + stokes_v * cos_phi,
        )

    return Observables(
        model=meta.model,
        b_label=meta.b_label,
        b_value=meta.b_value,
        topology=meta.topology,
        csi=meta.csi,
        m_star_msun=m_star,
        r_star_km=r_star,
        f_surf=f_surf,
        b_pole=b_pole,
        h_eff=h_eff,
        n_perp=n_perp,
        n_par=n_par,
        delta_phi=delta_phi,
        stokes_i=stokes_i,
        stokes_q=stokes_q,
        stokes_u=stokes_u,
        stokes_v=stokes_v,
    )


def write_csv(path: Path, rows: Iterable[Observables]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "model",
                "b_label",
                "b_value",
                "topology",
                "csi",
                "m_star_msun",
                "r_star_km",
                "f_surf",
                "b_pole",
                "h_eff",
                "n_perp",
                "n_par",
                "delta_phi",
                "stokes_i",
                "stokes_q",
                "stokes_u",
                "stokes_v",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.model,
                    r.b_label,
                    r.b_value,
                    r.topology,
                    r.csi,
                    r.m_star_msun,
                    r.r_star_km,
                    r.f_surf,
                    r.b_pole,
                    r.h_eff,
                    r.n_perp,
                    r.n_par,
                    r.delta_phi,
                    r.stokes_i,
                    r.stokes_q,
                    r.stokes_u,
                    r.stokes_v,
                ]
            )


def main() -> None:
    args = parse_args()
    datasets = discover_datasets(args.input_root)
    if not datasets:
        raise SystemExit(f"Nenhum dataset em {args.input_root}")

    rows: List[Observables] = []
    for meta in datasets:
        data = load_eos(meta.eos_path)
        if data is None:
            continue
        obs = compute_observables(meta, data, args.steps, args.r_factor, args.omega)
        if obs is not None:
            rows.append(obs)

    write_csv(args.output, rows)
    print(f"[OK] Observables saved to {args.output}")

    if args.plot:
        try:
            import plot_modmax_observables as plotter
        except Exception as exc:
            raise SystemExit(f"Falha ao carregar plot_modmax_observables: {exc}")

        plotter.setup_style()
        plot_rows = plotter.load_rows(args.output)
        grouped = plotter.group_rows(plot_rows)
        for (model, b_label, topology), group in grouped.items():
            tag = f"{model}_B_{b_label}_{topology}"
            plotter.plot_indices(group, args.plot_output_dir, tag, args.plot_dpi)
            plotter.plot_phase(group, args.plot_output_dir, tag, args.plot_dpi)
            plotter.plot_h_field(group, args.plot_output_dir, tag, args.plot_dpi)
            plotter.plot_stokes(group, args.plot_output_dir, tag, args.plot_dpi)

        print(f"[OK] Plots saved to {args.plot_output_dir}")


if __name__ == "__main__":
    main()
