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


EOS_RESULTS_SIZE = 34
EOS_WITH_MR_COLUMNS = EOS_RESULTS_SIZE + 3


@dataclass
class Metrics:
    epsilon: float
    m_x_over_mn: float
    g_d: float
    y_chi: float
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
    """Read the final three M-R columns appended to an EOS row."""
    masses: List[float] = []
    radii: List[float] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < EOS_WITH_MR_COLUMNS:
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

    mass_arr = np.asarray(masses, dtype=float)
    radius_arr = np.asarray(radii, dtype=float)
    order = np.argsort(radius_arr)
    return mass_arr[order], radius_arr[order]


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
        m_x_over_mn = safe_float(row.get("m_x_over_mN", ""))
        g_d = safe_float(row.get("g_d", ""))
        y_chi = safe_float(row.get("y_chi", ""))
        eos_file = row.get("eos_file", "")
        if not eos_file:
            continue

        eos_path = os.path.join(base_dir, eos_file)
        try:
            _, eps_d, p_d = read_eos(eos_path)
        except (OSError, ValueError):
            continue

        rms_rel, mean_gain, max_cs2, frac_invalid, n_overlap = compute_metrics(
            eps_h, p_h, eps_d, p_d, min_overlap
        )
        score = score_metrics(rms_rel, mean_gain, frac_invalid)

        metrics_list.append(
            Metrics(
                epsilon=epsilon,
                m_x_over_mn=m_x_over_mn,
                g_d=g_d,
                y_chi=y_chi,
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
                "m_x_over_mN",
                "g_d",
                "y_chi",
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
                    f"{m.m_x_over_mn:.6e}",
                    f"{m.g_d:.6e}",
                    f"{m.y_chi:.6e}",
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
        label = f"eps={m.epsilon:.1e}, mx/MN={m.m_x_over_mn:.1e}, gd={m.g_d:.2f}, Ychi={m.y_chi:.3f}"
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
        ("m_x_over_mN", [m.m_x_over_mn for m in rows]),
        ("g_d", [m.g_d for m in rows]),
        ("y_chi", [m.y_chi for m in rows]),
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
    plt.xlim(8.0, 16.0)
    plt.ylim(0.0, 2.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return True


def plot_pairwise_scores(out_path: str, rows: List[Metrics], top_frac: float) -> None:
    eps = np.array([m.epsilon for m in rows])
    mx = np.array([m.m_x_over_mn for m in rows])
    gd = np.array([m.g_d for m in rows])
    ychi = np.array([m.y_chi for m in rows])
    score = np.array([m.score for m in rows])

    n_top = max(1, int(len(rows) * top_frac))
    top_mask = np.zeros(len(rows), dtype=bool)
    top_mask[:n_top] = True

    params = {
        "epsilon": eps,
        "m_x_over_mN": mx,
        "g_d": gd,
        "y_chi": ychi,
    }
    pairs = [
        ("epsilon", "m_x_over_mN"),
        ("epsilon", "g_d"),
        ("epsilon", "y_chi"),
        ("m_x_over_mN", "g_d"),
        ("m_x_over_mN", "y_chi"),
        ("g_d", "y_chi"),
    ]

    fig, axs = plt.subplots(2, 3, figsize=(11.0, 6.0))
    mappable = None

    for ax, (x_name, y_name) in zip(axs.ravel(), pairs):
        x = params[x_name]
        y = params[y_name]
        mappable = ax.scatter(x, y, c=score, cmap="viridis", s=12, alpha=0.6)
        ax.scatter(
            x[top_mask],
            y[top_mask],
            facecolors="none",
            edgecolors="#d62728",
            s=36,
            linewidths=0.8,
        )
        if x_name == "epsilon":
            ax.set_xscale("log")
        if y_name == "epsilon":
            ax.set_yscale("log")
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)

    if mappable is not None:
        fig.colorbar(mappable, ax=axs.ravel().tolist(), label="score", shrink=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def write_optimal_ranges(out_path: str, rows: List[Metrics], top_frac: float) -> None:
    n_top = max(1, int(len(rows) * top_frac))
    top = rows[:n_top]

    eps = np.array([m.epsilon for m in top])
    mx = np.array([m.m_x_over_mn for m in top])
    gd = np.array([m.g_d for m in top])
    ychi = np.array([m.y_chi for m in top])
    score = np.array([m.score for m in top])

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Optimal ranges from top {:.1f}% (by score)\n".format(top_frac * 100.0))
        f.write("count={}\n\n".format(len(top)))
        f.write("epsilon: min={:.3e} max={:.3e} mean={:.3e}\n".format(eps.min(), eps.max(), eps.mean()))
        f.write("m_x/MN:  min={:.3e} max={:.3e} mean={:.3e}\n".format(mx.min(), mx.max(), mx.mean()))
        f.write("g_d:     min={:.3e} max={:.3e} mean={:.3e}\n".format(gd.min(), gd.max(), gd.mean()))
        f.write("y_chi:   min={:.3e} max={:.3e} mean={:.3e}\n".format(ychi.min(), ychi.max(), ychi.mean()))
        f.write("\nscore:  min={:.3e} max={:.3e} mean={:.3e}\n".format(score.min(), score.max(), score.mean()))


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze dark photon scan outputs")
    parser.add_argument("--base-dir", default="output/darkphotons_scan/GM1")
    parser.add_argument("--top", type=int, default=8)
    parser.add_argument("--min-overlap", type=int, default=30)
    parser.add_argument("--max-rms", type=float, default=0.5)
    parser.add_argument("--max-invalid", type=float, default=0.02)
    parser.add_argument("--top-frac", type=float, default=0.1)
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
    write_optimal_ranges(os.path.join(out_dir, "optimal_ranges.txt"), rows_sorted, args.top_frac)

    print("Best parameters:")
    print(
        "  epsilon={:.3e} m_x/MN={:.3e} g_d={:.3e} y_chi={:.3e} score={:.3e}".format(
            best.epsilon, best.m_x_over_mn, best.g_d, best.y_chi, best.score
        )
    )
    print("  eos_file={}".format(best.eos_file))

    top_n = rows_sorted[: max(1, args.top)]
    print("Top candidates:")
    for m in top_n:
        print(
            "  score={:.3e} eps={:.1e} mx/MN={:.1e} gd={:.2f} Ychi={:.3f}"
            .format(m.score, m.epsilon, m.m_x_over_mn, m.g_d, m.y_chi)
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
        plot_pairwise_scores(
            os.path.join(out_dir, "pairwise_score.png"),
            rows_sorted,
            args.top_frac,
        )
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
