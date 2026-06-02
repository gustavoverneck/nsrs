#!/usr/bin/env python3
"""
Analyze dark photon scan outputs, compare to hadrons baseline, and rank parameters.

Usage examples:
  python darkphotons.py --base-dir output/darkphotons_scan/GM1
  python darkphotons.py --base-dir output/darkphotons_scan/GM1 --top 8 --no-plots
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


@dataclass
class Metrics:
    epsilon: float
    m_x: float
    g_d: float
    n_chi: float
    eos_file: str
    rms_rel_p: float
    mean_gain: float
    max_cs2: float
    frac_invalid: float
    n_overlap: int
    score: float


def read_eos(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read EOS file with at least three columns: nB, eps, P."""
    rho_vals: List[float] = []
    eps_vals: List[float] = []
    p_vals: List[float] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rho = float(parts[0])
                eps = float(parts[1])
                pres = float(parts[2])
            except ValueError:
                continue
            if eps > 0.0 and pres >= 0.0:
                rho_vals.append(rho)
                eps_vals.append(eps)
                p_vals.append(pres)

    if not eps_vals:
        raise ValueError(f"No valid EOS data found in {path}")

    rho_arr = np.asarray(rho_vals, dtype=float)
    eps_arr = np.asarray(eps_vals, dtype=float)
    p_arr = np.asarray(p_vals, dtype=float)

    # Sort by pressure for stable interpolation
    order = np.argsort(p_arr)
    return rho_arr[order], eps_arr[order], p_arr[order]


def read_mr_from_eos(path: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Read MR data from eos_with_mr files (expects >= 24 columns)."""
    masses: List[float] = []
    radii: List[float] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 24:
                continue
            try:
                mass = float(parts[-3])
                radius = float(parts[-2])
            except ValueError:
                continue
            if math.isfinite(mass) and math.isfinite(radius):
                masses.append(mass)
                radii.append(radius)

    if not masses:
        return None

    return np.asarray(masses, dtype=float), np.asarray(radii, dtype=float)


def compute_metrics(
    eps_h: np.ndarray,
    p_h: np.ndarray,
    eps_d: np.ndarray,
    p_d: np.ndarray,
    min_overlap: int,
) -> Tuple[float, float, float, float, int]:
    """Compare dark EOS to hadrons EOS using overlapping eps range."""
    eps_min = max(eps_h.min(), eps_d.min())
    eps_max = min(eps_h.max(), eps_d.max())

    mask_d = (eps_d >= eps_min) & (eps_d <= eps_max)
    eps_d = eps_d[mask_d]
    p_d = p_d[mask_d]

    if eps_d.size < min_overlap:
        return float("inf"), 0.0, float("inf"), 1.0, eps_d.size

    # Interpolate hadrons pressure at dark eps values
    p_h_interp = np.interp(eps_d, eps_h, p_h)

    # Relative pressure difference, avoid near-zero baseline
    denom = np.where(p_h_interp > 1e-8, p_h_interp, 1e-8)
    rel = (p_d - p_h_interp) / denom
    rms_rel = float(np.sqrt(np.mean(rel * rel)))
    mean_gain = float(np.mean(rel))

    # Causality/stability proxy: cs^2 = dP/dE
    dp = np.gradient(p_d)
    de = np.gradient(eps_d)
    cs2 = np.where(de != 0.0, dp / de, np.nan)
    invalid = (cs2 < 0.0) | (cs2 > 1.0) | ~np.isfinite(cs2)
    frac_invalid = float(np.mean(invalid))
    max_cs2 = float(np.nanmax(cs2)) if np.isfinite(cs2).any() else float("inf")

    return rms_rel, mean_gain, max_cs2, frac_invalid, eps_d.size


def score_metrics(rms_rel: float, mean_gain: float, frac_invalid: float) -> float:
    """Lower is better. Emphasize physical consistency then proximity to baseline."""
    return rms_rel + 3.0 * frac_invalid + 0.2 * abs(mean_gain)


def parse_summary(path: str) -> Tuple[str, List[Dict[str, str]]]:
    with open(path, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    hadrons = [r for r in rows if r.get("label") == "hadrons"]
    if not hadrons:
        raise ValueError("summary.csv missing hadrons row")
    hadrons_file = hadrons[0].get("eos_file")
    if not hadrons_file:
        raise ValueError("summary.csv hadrons row missing eos_file")

    dark_rows = [r for r in rows if r.get("label") == "darkphotons"]
    return hadrons_file, dark_rows


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def analyze(base_dir: str, min_overlap: int) -> Tuple[List[Metrics], str]:
    summary_path = os.path.join(base_dir, "summary.csv")
    hadrons_file, rows = parse_summary(summary_path)

    hadrons_path = os.path.join(base_dir, hadrons_file)
    _, eps_h, p_h = read_eos(hadrons_path)

    metrics_list: List[Metrics] = []

    for row in rows:
        epsilon = safe_float(row.get("epsilon", ""))
        m_x = safe_float(row.get("m_x", ""))
        g_d = safe_float(row.get("g_d", ""))
        n_chi = safe_float(row.get("n_chi", ""))
        eos_file = row.get("eos_file", "")
        if not eos_file:
            continue

        eos_path = os.path.join(base_dir, eos_file)
        try:
            _, eps_d, p_d = read_eos(eos_path)
        except ValueError:
            continue

        rms_rel, mean_gain, max_cs2, frac_invalid, n_overlap = compute_metrics(
            eps_h, p_h, eps_d, p_d, min_overlap
        )
        score = score_metrics(rms_rel, mean_gain, frac_invalid)

        metrics_list.append(
            Metrics(
                epsilon=epsilon,
                m_x=m_x,
                g_d=g_d,
                n_chi=n_chi,
                eos_file=eos_file,
                rms_rel_p=rms_rel,
                mean_gain=mean_gain,
                max_cs2=max_cs2,
                frac_invalid=frac_invalid,
                n_overlap=n_overlap,
                score=score,
            )
        )

    return metrics_list, hadrons_path


def write_analysis_csv(out_path: str, rows: Iterable[Metrics]) -> None:
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "epsilon",
                "m_x",
                "g_d",
                "n_chi",
                "eos_file",
                "rms_rel_p",
                "mean_gain",
                "max_cs2",
                "frac_invalid",
                "n_overlap",
                "score",
            ]
        )
        for m in rows:
            writer.writerow(
                [
                    f"{m.epsilon:.6e}",
                    f"{m.m_x:.6e}",
                    f"{m.g_d:.6e}",
                    f"{m.n_chi:.6e}",
                    m.eos_file,
                    f"{m.rms_rel_p:.6e}",
                    f"{m.mean_gain:.6e}",
                    f"{m.max_cs2:.6e}",
                    f"{m.frac_invalid:.6e}",
                    str(m.n_overlap),
                    f"{m.score:.6e}",
                ]
            )


def plot_eos(
    out_path: str,
    hadrons_path: str,
    base_dir: str,
    best: List[Metrics],
    log_scale: bool,
) -> None:
    _, eps_h, p_h = read_eos(hadrons_path)
    plt.figure(figsize=(7.0, 5.0))
    plt.plot(eps_h, p_h, color="black", linewidth=2.0, label="hadrons")

    for m in best:
        eos_path = os.path.join(base_dir, m.eos_file)
        _, eps_d, p_d = read_eos(eos_path)
        label = f"eps={m.epsilon:.1e}, mx={m.m_x:.1e}, gd={m.g_d:.2f}, nchi={m.n_chi:.3f}"
        plt.plot(eps_d, p_d, linewidth=1.0, alpha=0.8, label=label)

    plt.xlabel("Energy density (MeV/fm^3)")
    plt.ylabel("Pressure (MeV/fm^3)")
    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
    plt.legend(fontsize=7, ncol=1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_score_scatter(out_path: str, rows: List[Metrics]) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.0))
    params = [
        ("epsilon", [m.epsilon for m in rows]),
        ("m_x", [m.m_x for m in rows]),
        ("g_d", [m.g_d for m in rows]),
        ("n_chi", [m.n_chi for m in rows]),
    ]
    scores = [m.score for m in rows]

    for ax, (name, values) in zip(axs.ravel(), params):
        ax.scatter(values, scores, s=12, alpha=0.7)
        ax.set_xlabel(name)
        ax.set_ylabel("score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_mr_all(out_path: str, hadrons_path: str, base_dir: str, rows: List[Metrics]) -> bool:
    hadrons_mr = read_mr_from_eos(hadrons_path)
    if hadrons_mr is None:
        return False

    plt.figure(figsize=(7.0, 5.0))
    m_h, r_h = hadrons_mr
    plt.plot(r_h, m_h, color="black", linewidth=2.0, label="hadrons")

    for m in rows:
        eos_path = os.path.join(base_dir, m.eos_file)
        mr = read_mr_from_eos(eos_path)
        if mr is None:
            continue
        mass, radius = mr
        plt.plot(radius, mass, color="#1f77b4", alpha=0.25, linewidth=0.8)

    plt.xlabel("Radius (km)")
    plt.ylabel("Mass (M_sun)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze dark photon scan outputs")
    parser.add_argument("--base-dir", default="output/darkphotons_scan/GM1")
    parser.add_argument("--top", type=int, default=8)
    parser.add_argument("--min-overlap", type=int, default=30)
    parser.add_argument("--max-rms", type=float, default=0.5)
    parser.add_argument("--max-invalid", type=float, default=0.02)
    parser.add_argument("--log-scale", action="store_true")
    parser.add_argument("--linear-scale", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    base_dir = args.base_dir
    os.makedirs(base_dir, exist_ok=True)
    out_dir = os.path.join(base_dir, "analysis")
    os.makedirs(out_dir, exist_ok=True)

    rows, hadrons_path = analyze(base_dir, args.min_overlap)
    if not rows:
        raise SystemExit("No dark photons rows found in summary.csv")

    rows_sorted = sorted(rows, key=lambda m: m.score)
    probable = [
        m
        for m in rows_sorted
        if m.rms_rel_p <= args.max_rms
        and m.frac_invalid <= args.max_invalid
        and m.n_overlap >= args.min_overlap
    ]

    best = probable[0] if probable else rows_sorted[0]

    analysis_csv = os.path.join(out_dir, "analysis.csv")
    write_analysis_csv(analysis_csv, rows_sorted)

    print("Best parameters:")
    print(
        "  epsilon={:.3e} m_x={:.3e} g_d={:.3e} n_chi={:.3e} score={:.3e}".format(
            best.epsilon, best.m_x, best.g_d, best.n_chi, best.score
        )
    )
    print("  eos_file={}".format(best.eos_file))

    top_n = rows_sorted[: max(1, args.top)]
    print("Top candidates:")
    for m in top_n:
        print(
            "  score={:.3e} eps={:.1e} mx={:.1e} gd={:.2f} nchi={:.3f}"
            .format(m.score, m.epsilon, m.m_x, m.g_d, m.n_chi)
        )

    if not args.no_plots:
        use_log = args.log_scale or not args.linear_scale
        plot_eos(
            os.path.join(out_dir, "eos_compare_top.png"),
            hadrons_path,
            base_dir,
            top_n,
            use_log,
        )
        plot_score_scatter(os.path.join(out_dir, "score_scatter.png"), rows_sorted)
        mr_written = plot_mr_all(
            os.path.join(out_dir, "mr_all.png"),
            hadrons_path,
            base_dir,
            rows_sorted,
        )
        if not mr_written:
            print("MR plot skipped: eos files do not include MR columns.")


if __name__ == "__main__":
    main()
