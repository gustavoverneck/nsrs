#!/usr/bin/env python3
"""
Publication-quality plots for dark photon EoS and TOV scan outputs.

Expected columns (or close matches):
['epsilon', 'm_x', 'g_d', 'n_chi', 'mun', 'n_n0', 'ener', 'press', 'n_e', 'n_mu',
 'n_n', 'n_p', 'L0', 'Sm', 'S0', 'Sp', 'Xm', 'X0', 'sigma', 'omega', 'rho',
 'm_eff', 'ebsd', 'bdd', 'Radius_km', 'Mass_Msun']
"""

from __future__ import annotations

import argparse
import os
import re
import glob
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 14,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "legend.fontsize": 11,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "figure.dpi": 300,
            "savefig.dpi": 300,
        }
    )
    sns.set_theme(style="whitegrid")


def read_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path, sep=",", engine="python")

    if df.empty:
        raise ValueError("Input file is empty or could not be parsed.")

    return df


def parse_params_from_filename(filename: str) -> dict:
    match = re.search(
        r"eps_(?P<eps>[^_]+)_mx_(?P<mx>[^_]+)_gd_(?P<gd>[^_]+)_nchi_(?P<nchi>[^.]+)",
        filename,
    )
    if not match:
        return {}

    def to_float(value: str) -> float:
        try:
            return float(value.replace("e", "e"))
        except ValueError:
            return float("nan")

    return {
        "epsilon": to_float(match.group("eps")),
        "m_x": to_float(match.group("mx")),
        "g_d": to_float(match.group("gd")),
        "n_chi": to_float(match.group("nchi")),
    }


def load_eos_file(path: str) -> pd.DataFrame:
    columns = [
        "n_n0",
        "ener",
        "press",
        "n_e",
        "n_mu",
        "n_n",
        "n_p",
        "L0",
        "Sm",
        "S0",
        "Sp",
        "Xm",
        "X0",
        "sigma",
        "omega",
        "rho",
        "m_eff",
        "mun",
        "mue",
        "ebsd",
        "bdd",
        "Mass_Msun",
        "Radius_km",
        "Mbaryon_Msun",
    ]

    records: List[List[float]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                values = [float(val) for val in parts]
            except ValueError:
                continue
            records.append(values)

    if not records:
        return pd.DataFrame()

    max_len = max(len(row) for row in records)
    used_cols = columns[:max_len]
    return pd.DataFrame(records, columns=used_cols)


def load_eos_directory(path: str) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    patterns = ("eos_", "eos_hadrons")

    for root, _, files in os.walk(path):
        for name in files:
            if not name.endswith(".dat"):
                continue
            if not name.startswith(patterns):
                continue
            file_path = os.path.join(root, name)
            params = parse_params_from_filename(name)
            df = load_eos_file(file_path)
            if df.empty:
                continue
            for key, value in params.items():
                df[key] = value
            df["eos_file"] = name
            df["eos_dir"] = root
            df["is_baseline"] = name.startswith("eos_hadrons")
            rows.append(df)

    if not rows:
        raise ValueError("No eos_*.dat files found under the directory.")

    return pd.concat(rows, ignore_index=True)


def coerce_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def group_label(row: pd.Series, group_cols: Sequence[str]) -> str:
    parts = []
    for col in group_cols:
        if col not in row:
            continue
        val = row[col]
        if pd.isna(val):
            continue
        if col == "epsilon":
            parts.append(f"eps={val:.1e}")
        else:
            parts.append(f"{col}={val:g}")
    return ", ".join(parts)


def plot_mr_diagram(
    df: pd.DataFrame,
    group_cols: Sequence[str],
    out_path: str,
    n_chi_col: str = "n_chi",
) -> None:
    required = ["Radius_km", "Mass_Msun"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    fig, ax = plt.subplots(figsize=(9.0, 5.2))

    baseline = df[df.get("is_baseline", False)].dropna(subset=required)
    if not baseline.empty:
        base_curve = baseline.sort_values("Radius_km")
        ax.plot(
            base_curve["Radius_km"],
            base_curve["Mass_Msun"],
            color="black",
            linewidth=2.2,
            label="hadrons baseline",
        )

    dark = df[~df.get("is_baseline", False)].dropna(subset=required)
    valid_cols = [col for col in group_cols if col in dark.columns]
    if not valid_cols:
        valid_cols = ["eos_file"]
    grouped = dark.groupby(valid_cols)
    for _, block in grouped:
        block = block.sort_values("Radius_km")
        ax.plot(
            block["Radius_km"],
            block["Mass_Msun"],
            linewidth=1.0,
            alpha=0.2,
            color="#1f77b4",
        )

    ax.set_xlim(5.0, 15.0)
    ax.set_ylim(0.0, 2.5)
    ax.set_xlabel("Radius (km)")
    ax.set_ylabel("Mass (M_sun)")

    # Observational constraints
    ax.axhspan(2.01 - 0.04, 2.01 + 0.04, color="#1f77b4", alpha=0.2, label="PSR J0348+0432")
    ax.axhspan(2.08 - 0.07, 2.08 + 0.07, color="#ff7f0e", alpha=0.2, label="PSR J0740+6620")

    ax.legend(loc="lower left", ncol=1, frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def plot_eos(
    df: pd.DataFrame,
    group_col: str,
    out_path: str,
    max_eps: float = 1500.0,
) -> None:
    required = ["ener", "press", group_col]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    data = df.dropna(subset=required).copy()
    data = data[(data["ener"] > 0.0) & (data["press"] > 0.0)]
    if data.empty:
        raise ValueError("No positive EoS data available for plotting.")

    baseline = data[data.get("is_baseline", False)]
    if not baseline.empty:
        base_curve = baseline.sort_values("ener")
        ax.plot(
            base_curve["ener"],
            base_curve["press"],
            color="black",
            linewidth=2.2,
            label="hadrons baseline",
        )

    dark = data[~data.get("is_baseline", False)]
    grouped = dark.groupby(group_col)
    values = sorted([value for value, _ in grouped])
    if values:
        vmin = min(values)
        vmax = max(values)
    else:
        vmin = 0.0
        vmax = 1.0
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")
    for value, block in grouped:
        block = block.sort_values("ener")
        ax.plot(
            block["ener"],
            block["press"],
            linewidth=1.2,
            color=cmap(norm(value)),
            alpha=0.85,
        )

    ax.set_xlabel("Energy density (MeV/fm^3)")
    ax.set_ylabel("Pressure (MeV/fm^3)")
    min_eps = max(1e-3, float(data["ener"].min()))
    max_eps_val = float(data["ener"].max())
    ax.set_xlim(min_eps, min(max_eps_val, max_eps))
    ax.set_yscale("log")
    ax.set_xscale("log")
    if values:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label(group_col)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def plot_particle_fractions(
    df: pd.DataFrame,
    out_path: str,
    selector: Optional[Tuple[str, float]] = None,
) -> None:
    required = ["n_n0"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    particle_cols = [
        "n_n",
        "n_p",
        "n_e",
        "n_mu",
        "L0",
        "Sm",
        "S0",
        "Sp",
        "Xm",
        "X0",
    ]
    available = [c for c in particle_cols if c in df.columns]
    if not available:
        raise ValueError("No particle density columns found.")

    data = df.copy()
    if selector is not None:
        key, value = selector
        if key in data.columns:
            data = data.loc[np.isclose(data[key], value)]

    data = data.dropna(subset=["n_n0"] + available)
    data = data.sort_values("n_n0")

    total = data[available].sum(axis=1)
    total = total.replace(0.0, np.nan)

    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    baryons = {"n_n", "n_p", "L0", "Sm", "S0", "Sp", "Xm", "X0"}
    for col in available:
        y = data[col] / total
        if col in baryons:
            ax.plot(data["n_n0"], y, linewidth=1.4, label=col)
        else:
            ax.plot(data["n_n0"], y, linewidth=1.4, linestyle="--", label=col)

    ax.set_xlabel("n/n0")
    ax.set_ylabel("Particle fraction")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.0)
    ax.legend(loc="upper right", ncol=2, frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def plot_sound_speed(
    df: pd.DataFrame,
    group_col: str,
    out_path: str,
) -> None:
    required = ["n_n0", "ener", "press", group_col]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    grouped = df.dropna(subset=required).groupby(group_col)
    values = sorted([value for value, _ in grouped])
    if values:
        vmin = min(values)
        vmax = max(values)
    else:
        vmin = 0.0
        vmax = 1.0
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")
    for value, block in grouped:
        block = block.sort_values("n_n0").dropna(subset=["ener", "press", "n_n0"])
        ener = block["ener"].to_numpy()
        press = block["press"].to_numpy()
        n_n0 = block["n_n0"].to_numpy()
        if len(ener) < 3:
            continue
        de = np.diff(ener)
        dp = np.diff(press)
        valid = de != 0.0
        if not np.any(valid):
            continue
        vs2 = np.full_like(dp, np.nan, dtype=float)
        vs2[valid] = dp[valid] / de[valid]
        vs2 = np.clip(vs2, -10.0, 10.0)
        n_mid = 0.5 * (n_n0[1:] + n_n0[:-1])
        ax.plot(n_mid, vs2, linewidth=1.2, color=cmap(norm(value)), alpha=0.85)

    ax.axhline(1.0, color="red", linestyle="--", linewidth=1.0, label="causality")
    ax.axhline(1.0 / 3.0, color="#555555", linestyle=":", linewidth=1.0, label="conformal")

    ax.set_xlabel("n/n0")
    ax.set_ylabel("v_s^2 / c^2")
    ax.set_ylim(-0.05, 1.2)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(group_col)
    fig.tight_layout()
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Publication plots for dark photon EoS/TOV scan")
    parser.add_argument("--input", default="output/darkphotons_scan")
    parser.add_argument("--out-dir", default="plots")
    parser.add_argument("--group-cols", nargs="*", default=["n_chi", "epsilon"])
    parser.add_argument("--eos-group", default="n_chi")
    parser.add_argument("--sound-group", default="n_chi")
    parser.add_argument("--fchi-slice", type=float, default=None)
    args = parser.parse_args()

    configure_matplotlib()

    input_path = args.input
    if not os.path.exists(input_path):
        candidates = glob.glob(os.path.join("output", "darkphotons_scan", "*"))
        candidates = [c for c in candidates if os.path.isdir(c)]
        if len(candidates) == 1:
            input_path = candidates[0]
        else:
            raise FileNotFoundError(
                "Input not found. Provide --input as a directory with eos_*.dat files or a single file."
            )

    if os.path.isdir(input_path):
        df = load_eos_directory(input_path)
    else:
        df = read_table(input_path)
    df = coerce_numeric(
        df,
        [
            "epsilon",
            "m_x",
            "g_d",
            "n_chi",
            "n_chi",
            "mun",
            "n_n0",
            "ener",
            "press",
            "n_e",
            "n_mu",
            "n_n",
            "n_p",
            "L0",
            "Sm",
            "S0",
            "Sp",
            "Xm",
            "X0",
            "sigma",
            "omega",
            "rho",
            "m_eff",
            "ebsd",
            "bdd",
            "Radius_km",
            "Mass_Msun",
        ],
    )

    os.makedirs(args.out_dir, exist_ok=True)

    plot_mr_diagram(
        df,
        group_cols=args.group_cols,
        out_path=os.path.join(args.out_dir, "mr_diagram.pdf"),
    )

    plot_eos(
        df,
        group_col=args.eos_group,
        out_path=os.path.join(args.out_dir, "eos_pressure.pdf"),
    )

    selector = None
    if args.fchi_slice is not None:
        selector = ("n_chi", args.fchi_slice)

    plot_particle_fractions(
        df,
        out_path=os.path.join(args.out_dir, "particle_fractions.pdf"),
        selector=selector,
    )

    plot_sound_speed(
        df,
        group_col=args.sound_group,
        out_path=os.path.join(args.out_dir, "sound_speed.pdf"),
    )


if __name__ == "__main__":
    main()
