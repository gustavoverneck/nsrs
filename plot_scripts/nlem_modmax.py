#!/usr/bin/env python3
"""
Análise científica completa dos dados NLEM gerados em output/modmax.

Objetivos:
- Consolidar EoS, M-R e populações para GM1/GM3 e topologias isotrópica/anisotrópica.
- Quantificar o efeito de log10(csi) na estrutura estelar.
- Avaliar os limites de causalidade e estabilidade através da velocidade do som (c_s^2).
- Mapear os limiares de surgimento (onset thresholds) das partículas.
- Registrar como indisponível o diagnóstico de Landau enquanto B(n) não for exportado.
- Gerar tabelas e figuras prontas para publicação/inspeção posterior em Python.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


# ----------------------------
# Índices de colunas no eos.dat
# ----------------------------
COL_NB = 0
COL_EPS = 1
COL_P = 2

# Populações (conforme usado em src/bin/magtop.rs)
COL_NE = 3
COL_NMU = 4
COL_NN = 5
COL_NP = 6
COL_NL0 = 7
COL_NSM = 8
COL_NS0 = 9
COL_NSP = 10
COL_NXM = 11
COL_NX0 = 12

COL_MEFF = 16
COL_MUN = 17
COL_MUE = 18
COL_EMAG = 19
COL_MU_TOTAL = 20
EOS_RESULTS_SIZE = 34
COL_MR_MASS = EOS_RESULTS_SIZE
COL_MR_RADIUS = EOS_RESULTS_SIZE + 1

# Fator padrão. Se necessário, o código também detecta automaticamente
# se o mu_e parece normalizado e converte por 939.0 para MeV.
MUE_TO_MEV_FACTOR = 1.0

LOG_CSI_LABEL = r"$\log_{10}(\gamma)$"
CSI_LABEL = r"$\gamma$"
HBARC_MEV_FM = 197.3269804
MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0
COMPARISON_MR_TARGETS = [1e-22, 1e-20, 1e-15, 1e-10, 1e-5, 1e-1]


@dataclass
class Dataset:
    model: str
    b_label: str
    b_value: float
    topology: str
    csi: float
    log_csi: float
    file_path: Path
    data: np.ndarray


@dataclass
class SummaryRow:
    model: str
    b_label: str
    b_value: float
    topology: str
    csi: float
    log_csi: float
    n_eos_points: int
    n_mr_points: int
    max_mass_msun: float
    radius_at_max_km: float
    central_nb_n0: float
    central_eps_mevfm3: float
    central_p_mevfm3: float
    central_meff: float
    central_emag_mevfm3: float
    cs2_min: float
    cs2_max: float
    cs2_p95: float
    cs2_frac_negative: float
    cs2_frac_superluminal: float
    landau_nu_max: float
    landau_nu_p95: float
    eos_path: str

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Análise NLEM completa para dados em output/modmax"
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("output/modmax"),
        help="Pasta raiz dos dados gerados",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/modmax"),
        help="Pasta de saída para tabelas/figuras/relatório",
    )
    parser.add_argument(
        "--baseline-root",
        type=Path,
        default=Path("output/b"),
        help="Raiz dos dados baseline sem csi (B=0.0)",
    )
    parser.add_argument(
        "--max-curves-per-family",
        type=int,
        default=1000,
        help="Máximo de curvas por figura de família (evita poluição visual)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Resolução dos gráficos",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _safe_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except ValueError:
        return None


def parse_metadata_from_path(path: Path) -> Optional[Tuple[str, str, float, str, float, float]]:
    rx = re.compile(
        r".*/output/modmax/(?P<model>GM\d+)/B_(?P<b>[^/]+)/(?:"
        r"(?P<topology>isotropic|anisotropic)/)?csi_(?P<csi>[^/]+)/eos\.dat$"
    )
    m = rx.match(path.resolve().as_posix())
    if not m:
        return None

    model = m.group("model")
    b_label = m.group("b")
    topology = m.group("topology") or "isotropic"

    b_value = _safe_float(b_label)
    csi = _safe_float(m.group("csi"))
    if b_value is None or csi is None:
        return None

    log_csi = math.log10(csi) if csi > 0.0 else math.nan
    return model, b_label, b_value, topology, csi, log_csi


def load_eos(path: Path) -> Optional[np.ndarray]:
    try:
        arr = np.genfromtxt(path, comments="#")
    except Exception:
        return None

    if arr is None or np.size(arr) == 0:
        return None

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    if arr.shape[1] <= COL_MR_RADIUS:
        return None

    return arr


def discover_datasets(input_root: Path) -> List[Dataset]:
    datasets: List[Dataset] = []
    for eos_file in input_root.rglob("eos.dat"):
        meta = parse_metadata_from_path(eos_file)
        if meta is None:
            continue
        loaded = load_eos(eos_file)
        if loaded is None:
            continue

        model, b_label, b_value, topology, csi, log_csi = meta
        datasets.append(
            Dataset(
                model=model,
                b_label=b_label,
                b_value=b_value,
                topology=topology,
                csi=csi,
                log_csi=log_csi,
                file_path=eos_file,
                data=loaded,
            )
        )

    datasets.sort(key=lambda d: (d.model, d.topology, d.b_value, d.csi))
    return datasets


def compute_cs2_bounds(eps: np.ndarray, p: np.ndarray) -> Tuple[float, float]:
    mask = np.isfinite(eps) & np.isfinite(p)
    eps = eps[mask]
    p = p[mask]

    if eps.size < 3:
        return math.nan, math.nan

    idx = np.argsort(eps)
    eps = eps[idx]
    p = p[idx]

    de = np.diff(eps)
    dp = np.diff(p)
    good = de > 1e-14
    if not np.any(good):
        return math.nan, math.nan

    cs2 = dp[good] / de[good]
    cs2 = cs2[np.isfinite(cs2)]
    if cs2.size == 0:
        return math.nan, math.nan

    return float(np.min(cs2)), float(np.max(cs2))


def compute_cs2_stats(eps: np.ndarray, p: np.ndarray) -> Tuple[float, float, float, float, float]:
    """
    Retorna estatísticas robustas de c_s^2:
    (min, max, p95, fração<0, fração>1)
    """
    mask = np.isfinite(eps) & np.isfinite(p)
    eps = eps[mask]
    p = p[mask]
    if eps.size < 3:
        return math.nan, math.nan, math.nan, math.nan, math.nan

    idx = np.argsort(eps)
    eps = eps[idx]
    p = p[idx]

    de = np.diff(eps)
    dp = np.diff(p)
    good = de > 1e-14
    if not np.any(good):
        return math.nan, math.nan, math.nan, math.nan, math.nan

    cs2 = dp[good] / de[good]
    cs2 = cs2[np.isfinite(cs2)]
    if cs2.size == 0:
        return math.nan, math.nan, math.nan, math.nan, math.nan

    cmin = float(np.min(cs2))
    cmax = float(np.max(cs2))
    cp95 = float(np.percentile(cs2, 95.0))
    frac_neg = float(np.mean(cs2 < 0.0))
    frac_sup = float(np.mean(cs2 > 1.0))
    return cmin, cmax, cp95, frac_neg, frac_sup


def valid_mr_mask(mass_col: np.ndarray, radius_col: np.ndarray) -> np.ndarray:
    return (
        np.isfinite(mass_col)
        & np.isfinite(radius_col)
        & (mass_col > 0.0)
        & (radius_col > 0.0)
        & (mass_col <= MAX_VALID_MASS_MSUN)
        & (radius_col <= MAX_VALID_RADIUS_KM)
    )


def _load_mr_from_eos(path: Path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    try:
        arr = np.genfromtxt(path, comments="#")
    except Exception:
        return None

    if arr is None or np.size(arr) == 0:
        return None

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    if arr.shape[1] < 5:
        return None

    if arr.shape[1] <= COL_MR_RADIUS:
        return None
    mass = arr[:, COL_MR_MASS].astype(float)
    radius = arr[:, COL_MR_RADIUS].astype(float)
    mask = valid_mr_mask(mass, radius)
    if not np.any(mask):
        return None

    return radius[mask], mass[mask]


def _load_baseline_b0_curve(baseline_root: Path, model: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    eos_path = baseline_root / model / "B_0.00e0" / "default" / "eos.dat"
    if not eos_path.exists():
        return None
    return _load_mr_from_eos(eos_path)


def summarize_dataset(ds: Dataset) -> SummaryRow:
    arr = ds.data

    mass_col = arr[:, COL_MR_MASS]
    radius_col = arr[:, COL_MR_RADIUS]
    mr_mask = valid_mr_mask(mass_col, radius_col)

    n_mr_points = int(np.count_nonzero(mr_mask))
    n_eos_points = int(arr.shape[0])

    if n_mr_points > 0:
        mr_indices = np.nonzero(mr_mask)[0]
        local = np.argmax(mass_col[mr_mask])
        i_max = int(mr_indices[local])

        max_mass = float(arr[i_max, COL_MR_MASS])
        r_at_max = float(arr[i_max, COL_MR_RADIUS])
        central_nb = float(arr[i_max, COL_NB])
        central_eps = float(arr[i_max, COL_EPS])
        central_p = float(arr[i_max, COL_P])
        central_meff = float(arr[i_max, COL_MEFF]) if arr.shape[1] > COL_MEFF else math.nan
        central_emag = float(arr[i_max, COL_EMAG]) if arr.shape[1] > COL_EMAG else math.nan
    else:
        max_mass = math.nan
        r_at_max = math.nan
        central_nb = math.nan
        central_eps = math.nan
        central_p = math.nan
        central_meff = math.nan
        central_emag = math.nan

    cs2_min, cs2_max, cs2_p95, cs2_frac_negative, cs2_frac_superluminal = compute_cs2_stats(
        arr[:, COL_EPS], arr[:, COL_P]
    )

    # The EOS schema does not export local B(n); column 20 is mu_total/M_N.
    landau_nu_max = math.nan
    landau_nu_p95 = math.nan

    return SummaryRow(
        model=ds.model,
        b_label=ds.b_label,
        b_value=ds.b_value,
        topology=ds.topology,
        csi=ds.csi,
        log_csi=ds.log_csi,
        n_eos_points=n_eos_points,
        n_mr_points=n_mr_points,
        max_mass_msun=max_mass,
        radius_at_max_km=r_at_max,
        central_nb_n0=central_nb,
        central_eps_mevfm3=central_eps,
        central_p_mevfm3=central_p,
        central_meff=central_meff,
        central_emag_mevfm3=central_emag,
        cs2_min=cs2_min,
        cs2_max=cs2_max,
        cs2_p95=cs2_p95,
        cs2_frac_negative=cs2_frac_negative,
        cs2_frac_superluminal=cs2_frac_superluminal,
        landau_nu_max=landau_nu_max,
        landau_nu_p95=landau_nu_p95,
        eos_path=ds.file_path.as_posix(),
    )


def write_summary_csv(rows: Sequence[SummaryRow], out_csv: Path) -> None:
    ensure_dir(out_csv.parent)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "b_label", "b_value", "topology", "csi", "log_csi",
            "n_eos_points", "n_mr_points", "max_mass_msun", "radius_at_max_km",
            "central_nb_n0", "central_eps_mevfm3", "central_p_mevfm3", "central_meff",
            "central_emag_mevfm3", "cs2_min", "cs2_max", "cs2_p95",
            "cs2_frac_negative", "cs2_frac_superluminal",
            "landau_nu_max", "landau_nu_p95", "eos_path",
        ])
        for r in rows:
            w.writerow([
                r.model, r.b_label, r.b_value, r.topology, r.csi, r.log_csi,
                r.n_eos_points, r.n_mr_points, r.max_mass_msun, r.radius_at_max_km,
                r.central_nb_n0, r.central_eps_mevfm3, r.central_p_mevfm3,
                r.central_meff, r.central_emag_mevfm3, r.cs2_min, r.cs2_max, r.cs2_p95,
                r.cs2_frac_negative, r.cs2_frac_superluminal,
                r.landau_nu_max, r.landau_nu_p95, r.eos_path,
            ])


def _linear_slope_r2(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return math.nan, math.nan
    x = x[mask]
    y = y[mask]
    a, b = np.polyfit(x, y, 1)
    y_hat = a * x + b
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0.0 else math.nan
    return float(a), float(r2)


def write_group_stats_csv(rows: Sequence[SummaryRow], out_csv: Path) -> None:
    """
    Estatísticas agregadas por (model, topology, B), incluindo sensibilidade
    a log10(csi) e qualidade de ajuste linear (R²).
    """
    ensure_dir(out_csv.parent)

    grouped: Dict[Tuple[str, str, str], List[SummaryRow]] = {}
    for r in rows:
        grouped.setdefault((r.model, r.topology, r.b_label), []).append(r)

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "topology", "b_label", "n_samples",
            "mass_slope", "mass_r2", "radius_slope", "radius_r2",
            "cs2max_slope", "cs2max_r2", "landau_max_slope", "landau_max_r2",
            "mean_cs2_frac_negative", "mean_cs2_frac_superluminal",
            "median_landau_nu_max", "median_landau_nu_p95",
        ])

        for (model, topo, b_label), sub in sorted(
            grouped.items(),
            key=lambda kv: (
                kv[0][0],
                kv[0][1],
                _safe_float(kv[0][2]) if _safe_float(kv[0][2]) is not None else math.inf,
            ),
        ):
            sub = sorted(sub, key=lambda x: x.log_csi)
            x = np.array([r.log_csi for r in sub], dtype=float)
            y_mass = np.array([r.max_mass_msun for r in sub], dtype=float)
            y_rad = np.array([r.radius_at_max_km for r in sub], dtype=float)
            y_cs2max = np.array([r.cs2_max for r in sub], dtype=float)
            y_landau = np.array([r.landau_nu_max for r in sub], dtype=float)

            mass_slope, mass_r2 = _linear_slope_r2(x, y_mass)
            radius_slope, radius_r2 = _linear_slope_r2(x, y_rad)
            cs2max_slope, cs2max_r2 = _linear_slope_r2(x, y_cs2max)
            landau_slope, landau_r2 = _linear_slope_r2(x, y_landau)

            frac_neg = np.array([r.cs2_frac_negative for r in sub], dtype=float)
            frac_sup = np.array([r.cs2_frac_superluminal for r in sub], dtype=float)
            landau_max = np.array([r.landau_nu_max for r in sub], dtype=float)
            landau_p95 = np.array([r.landau_nu_p95 for r in sub], dtype=float)

            w.writerow([
                model,
                topo,
                b_label,
                len(sub),
                mass_slope,
                mass_r2,
                radius_slope,
                radius_r2,
                cs2max_slope,
                cs2max_r2,
                landau_slope,
                landau_r2,
                float(np.nanmean(frac_neg)) if np.any(np.isfinite(frac_neg)) else math.nan,
                float(np.nanmean(frac_sup)) if np.any(np.isfinite(frac_sup)) else math.nan,
                float(np.nanmedian(landau_max)) if np.any(np.isfinite(landau_max)) else math.nan,
                float(np.nanmedian(landau_p95)) if np.any(np.isfinite(landau_p95)) else math.nan,
            ])


def group_key(ds: Dataset) -> Tuple[str, str, str]:
    return ds.model, ds.topology, ds.b_label


def _subsample_sorted(items: Sequence[Dataset], max_items: int) -> List[Dataset]:
    if len(items) <= max_items:
        return list(items)
    idx = np.linspace(0, len(items) - 1, max_items).round().astype(int)
    idx = np.unique(idx)
    return [items[i] for i in idx]


def _compute_electron_landau_nu_max(mu_e: np.ndarray, b_tesla: np.ndarray) -> np.ndarray:
    m_e_mev = 0.511
    b_crit_tesla = 4.414e9
    nu_max = np.zeros_like(mu_e, dtype=float)

    valid = (
        np.isfinite(mu_e)
        & np.isfinite(b_tesla)
        & (b_tesla > 0.0)
        & (mu_e * mu_e > m_e_mev * m_e_mev)
    )
    if not np.any(valid):
        return nu_max

    e_b_mev2 = (b_tesla[valid] / b_crit_tesla) * (m_e_mev * m_e_mev)
    denom_landau = 2.0 * e_b_mev2
    good = denom_landau > 0.0
    if not np.any(good):
        return nu_max

    valid_idx = np.nonzero(valid)[0]
    target_idx = valid_idx[good]
    numer = (mu_e[target_idx] * mu_e[target_idx]) - (m_e_mev * m_e_mev)
    nu_vals = np.floor(numer / denom_landau[good])
    nu_max[target_idx] = np.maximum(0.0, nu_vals)
    return nu_max


def _infer_mue_to_mev_factor(mu_e_raw: np.ndarray) -> float:
    """
    Heurística: se o valor máximo típico de mu_e for muito pequeno (ordem < 5),
    assume-se que está normalizado por m_n e converte para MeV.
    """
    finite = mu_e_raw[np.isfinite(mu_e_raw)]
    if finite.size == 0:
        return MUE_TO_MEV_FACTOR
    max_val = float(np.nanmax(finite))
    return 939.0 if max_val < 5.0 else 1.0


def _plot_single_landau_curve(
    ds: Dataset,
    cmap: mcolors.Colormap,
    cmin: float,
    denom: float,
) -> bool:
    # No local B(n) column is available in the current EOS schema.
    return False


def plot_landau_levels(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Calcula e plota o número máximo de níveis de Landau ocupados por elétrons
    em função da densidade bariônica para cada família (model, topology, B).
    """
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, topology, b_label = key
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        plot_set = arr_sorted
        if not plot_set:
            continue

        cvals = np.array([d.csi for d in plot_set], dtype=float)
        finite_cvals = cvals[np.isfinite(cvals)]
        if finite_cvals.size == 0:
            continue
        cmin, cmax = float(np.min(finite_cvals)), float(np.max(finite_cvals))
        denom = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0
        cmap = plt.get_cmap("plasma")

        plt.figure(figsize=(8, 6))
        plotted = False

        for d in plot_set:
            plotted = _plot_single_landau_curve(d, cmap=cmap, cmin=cmin, denom=denom) or plotted

        if plotted:
            plt.xlabel(r"Densidade bariônica $n_B/n_0$")
            plt.ylabel(r"Nível máximo de Landau eletrônico $\nu_{\max}^e$")
            plt.title(f"Quantização de Landau | {model} | {topology} | B={b_label} G")
            plt.grid(alpha=0.3)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
            cbar = plt.colorbar(sm, ax=plt.gca())
            cbar.set_label(CSI_LABEL)
            plt.tight_layout()
            plt.savefig(out_dir / f"landau_levels_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()


def plot_family_eos_mr(
    datasets: Sequence[Dataset],
    out_dir: Path,
    max_curves: int,
    dpi: int,
) -> None:
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, topology, b_label = key
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        plot_set = _subsample_sorted(arr_sorted, max_curves)
        cmap_eos = plt.get_cmap("viridis")
        cmap_mr = plt.get_cmap("plasma")
        cmap_cs2 = plt.get_cmap("cool")

        cvals = np.array([d.csi for d in plot_set], dtype=float)
        finite_cvals = cvals[np.isfinite(cvals) & (cvals > 0.0)]
        if finite_cvals.size == 0:
            continue
        cmin, cmax = float(np.min(finite_cvals)), float(np.max(finite_cvals))
        norm = mcolors.LogNorm(vmin=max(cmin, 1e-300), vmax=cmax)

        def _color_value(value: float):
            return norm(value if value > 0.0 else cmin)

        # --- EoS family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            color = cmap_eos(_color_value(d.csi))
            eps = d.data[:, COL_EPS]
            p = d.data[:, COL_P]
            mask = np.isfinite(eps) & np.isfinite(p) & (eps > 0.0) & (p > 0.0)
            if np.count_nonzero(mask) < 2:
                continue
            eps_fm4 = eps[mask] / HBARC_MEV_FM
            p_fm4 = p[mask] / HBARC_MEV_FM
            plt.plot(eps_fm4, p_fm4, color=color, alpha=0.8, lw=1.0)

        plt.xlabel(r"$\epsilon$ [fm$^{-4}$]")
        plt.ylabel(r"$P$ [fm$^{-4}$]")
        plt.title(f"EoS family | {model} | {topology} | B={b_label} G")
        plt.xlim(0, 7)
        plt.ylim(0, 4)
        plt.grid(alpha=0.25)

        sm = plt.cm.ScalarMappable(cmap=cmap_eos, norm=norm)
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(LOG_CSI_LABEL)
        plt.tight_layout()
        plt.savefig(out_dir / f"eos_family_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()

        # --- M-R family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            m = d.data[:, COL_MR_MASS]
            r = d.data[:, COL_MR_RADIUS]
            p_central = d.data[:, COL_P]

            mask = valid_mr_mask(m, r)
            if np.count_nonzero(mask) < 3:
                continue

            m_valid = m[mask]
            r_valid = r[mask]
            p_valid = p_central[mask]

            sort_idx = np.argsort(p_valid)
            color = cmap_mr(_color_value(d.csi))
            plt.plot(r_valid[sort_idx], m_valid[sort_idx], color=color, alpha=0.8, lw=1.0)

        plt.xlabel("Radius [km]")
        plt.ylabel(r"Mass [$M_\odot$]")
        plt.title(f"M-R family | {model} | {topology} | B={b_label} G")
        plt.xlim(9, 14.5)
        plt.ylim(1.0, 2.25)
        plt.grid(alpha=0.25)

        sm_mr = plt.cm.ScalarMappable(cmap=cmap_mr, norm=norm)
        cbar_mr = plt.colorbar(sm_mr, ax=plt.gca())
        cbar_mr.set_label(LOG_CSI_LABEL)
        plt.tight_layout()
        plt.savefig(out_dir / f"mr_family_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()

        # --- Speed of Sound (c_s^2) family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            eps = d.data[:, COL_EPS]
            p = d.data[:, COL_P]
            mask = np.isfinite(eps) & np.isfinite(p)
            eps_val = eps[mask]
            p_val = p[mask]

            if eps_val.size > 2:
                idx = np.argsort(eps_val)
                eps_sorted = eps_val[idx]
                p_sorted = p_val[idx]

                de = np.diff(eps_sorted)
                dp = np.diff(p_sorted)

                good = de > 1e-14
                if np.any(good):
                    cs2 = np.zeros_like(de)
                    cs2[good] = dp[good] / de[good]

                    color = cmap_cs2(_color_value(d.csi))
                    plt.plot(eps_sorted[:-1][good], cs2[good], color=color, alpha=0.8, lw=1.0)

        plt.axhline(1.0, color='red', linestyle='--', alpha=0.7, label='Causalidade ($c_s^2 = 1$)')
        plt.axhline(0.0, color='black', linestyle='--', alpha=0.7, label='Estabilidade ($c_s^2 = 0$)')
        plt.axhline(1/3, color='gray', linestyle=':', alpha=0.7, label='Limite Conforme ($c_s^2 = 1/3$)')

        plt.xlabel(r"$\epsilon$ [MeV/fm$^3$]")
        plt.ylabel(r"$c_s^2$ (Unidades de $c=1$)")
        plt.title(f"Speed of Sound | {model} | {topology} | B={b_label} G")
        plt.grid(alpha=0.25)
        plt.legend(loc='upper right', fontsize=9)

        sm_cs2 = plt.cm.ScalarMappable(cmap=cmap_cs2, norm=norm)
        cbar_cs2 = plt.colorbar(sm_cs2, ax=plt.gca())
        cbar_cs2.set_label(LOG_CSI_LABEL)

        plt.tight_layout()
        plt.savefig(out_dir / f"cs2_family_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()


def _select_nearest_mr_targets(arr_sorted: Sequence[Dataset], targets: Sequence[float]) -> List[Tuple[Dataset, float]]:
    selected: List[Tuple[Dataset, float]] = []
    used_ids: set[int] = set()

    for target in targets:
        chosen = min(arr_sorted, key=lambda d, target=target: abs(d.csi - target))
        chosen_id = id(chosen)
        if chosen_id in used_ids:
            continue
        used_ids.add(chosen_id)
        selected.append((chosen, target))

    return selected


def plot_mr_specific_targets(
    datasets: Sequence[Dataset],
    out_dir: Path,
    baseline_root: Path,
    targets: Sequence[float],
    dpi: int,
) -> None:
    ensure_dir(out_dir)
    comparison_dir = out_dir / "mr_specific"
    ensure_dir(comparison_dir)

    by_combo: Dict[Tuple[str, str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    cmap = plt.get_cmap("tab10")

    for key, arr in by_combo.items():
        model, topology, b_label = key
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        if not arr_sorted:
            continue

        baseline_mr = _load_baseline_b0_curve(baseline_root, model)
        # initialize zoom center variables (will be set if baseline exists)
        _zoom_center_r = None
        _zoom_center_m = None

        selected = _select_nearest_mr_targets(arr_sorted, targets)
        if not selected:
            selected = []

        plt.figure(figsize=(8, 6))
        plotted_curves: List[Tuple[np.ndarray, np.ndarray, str, Tuple[float, float, float, float]]] = []

        # prepare colormap for gamma-based coloring (like family plots)
        csi_values = [d.csi for d, _ in selected]
        if csi_values:
            cmin = min(csi_values)
            cmax = max(csi_values)
        else:
            cmin = 1e-25
            cmax = 1e-1

        if cmin > 0:
            norm = mcolors.LogNorm(vmin=max(cmin, 1e-300), vmax=cmax)
        else:
            norm = mcolors.Normalize(vmin=0, vmax=1)
        cmap_colors = plt.get_cmap("viridis")

        if baseline_mr is not None:
            baseline_radius, baseline_mass = baseline_mr
            baseline_label = r"$B=0.0$ G"
            baseline_color = (0.0, 0.0, 0.0, 0.95)
            plt.plot(
                baseline_radius,
                baseline_mass,
                color=baseline_color,
                alpha=0.95,
                lw=1.8,
                label=baseline_label,
            )
            plotted_curves.append((baseline_radius, baseline_mass, baseline_label, baseline_color))
            # compute center of zoom at the maximum-mass point of the B=0 curve
            try:
                if getattr(baseline_mass, "size", 0) > 0:
                    max_idx = int(np.nanargmax(baseline_mass))
                    _zoom_center_r = float(baseline_radius[max_idx])
                    _zoom_center_m = float(baseline_mass[max_idx])
                else:
                    _zoom_center_r = None
                    _zoom_center_m = None
            except Exception:
                _zoom_center_r = None
                _zoom_center_m = None

        for idx, (d, target) in enumerate(selected):
            m = d.data[:, COL_MR_MASS]
            r = d.data[:, COL_MR_RADIUS]
            p_central = d.data[:, COL_P]

            mask = valid_mr_mask(m, r)
            if np.count_nonzero(mask) < 3:
                continue

            m_valid = m[mask]
            r_valid = r[mask]
            p_valid = p_central[mask]
            sort_idx = np.argsort(p_valid)
            color = cmap_colors(norm(d.csi))
            curve_label = rf"$\gamma={d.csi:.2e}$"
            plt.plot(
                r_valid[sort_idx],
                m_valid[sort_idx],
                color=color,
                alpha=0.9,
                lw=1.5,
                linestyle="--",
                label=curve_label,
            )
            plotted_curves.append((r_valid[sort_idx], m_valid[sort_idx], curve_label, color))

        plt.xlabel("Radius [km]")
        plt.ylabel(r"Mass [$M_\odot$]")
        plt.title(f"MR comparison | {model} | {topology} | B={b_label} G")
        plt.xlim(8, 14.5)
        plt.ylim(1.0, 2.25)
        plt.grid(alpha=0.25)
        plt.axhline(1.6, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)
        plt.legend(fontsize=8)

        ax = plt.gca()
        axins = inset_axes(ax, width="42%", height="42%", loc="lower left", borderpad=3.0)
        for x_curve, y_curve, _label, color in plotted_curves:
            axins.scatter(x_curve, y_curve, color=color, alpha=0.8, s=20, edgecolors='none')
        # center inset on baseline B=0 maximum if available, otherwise keep default window
        # fetch computed zoom center if available
        try:
            if '_zoom_center_r' in locals() or '_zoom_center_r' in globals():
                center_r = _zoom_center_r
                center_m = _zoom_center_m
            else:
                center_r = None
                center_m = None
        except Exception:
            center_r = None
            center_m = None

        if center_r is not None and center_m is not None:
            x_span = ax.get_xlim()[1] - ax.get_xlim()[0]
            y_span = ax.get_ylim()[1] - ax.get_ylim()[0]
            # high-magnification zoom: 3% of full axes span
            dx = 0.03 * x_span
            dy = 0.03 * y_span
            # shift center up so the max appears near the bottom
            y_offset = +0.5 * dy
            axins.set_xlim(center_r - dx, center_r + dx)
            axins.set_ylim(center_m - dy + y_offset, center_m + dy + y_offset)
        else:
            axins.set_xlim(11.2, 14.5)
            axins.set_ylim(1.75, 2.25)
        axins.grid(alpha=0.2)
        axins.tick_params(axis='both', which='major', labelsize=8)
        # draw connection lines from zoom window to inset
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="gray", lw=1.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(comparison_dir / f"mr_specific_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()


def plot_population_snapshots(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    ensure_dir(out_dir)
    by_combo: Dict[Tuple[str, str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    species = [
        (COL_NE, "e-"),
        (COL_NMU, "mu-"),
        (COL_NN, "n"),
        (COL_NP, "p"),
        (COL_NL0, "Lambda0"),
        (COL_NSM, "Sigma-"),
        (COL_NS0, "Sigma0"),
        (COL_NSP, "Sigma+"),
        (COL_NXM, "Xi-"),
        (COL_NX0, "Xi0"),
    ]

    for key, arr in by_combo.items():
        model, topology, b_label = key
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        if not arr_sorted:
            continue

        picks = np.unique(np.linspace(0, len(arr_sorted) - 1, min(3, len(arr_sorted))).round().astype(int)).tolist()

        ncols = len(picks)
        fig, axs = plt.subplots(1, ncols, figsize=(6 * ncols, 5), squeeze=False)

        for j, idx in enumerate(picks):
            ds = arr_sorted[idx]
            ax = axs[0, j]
            x = ds.data[:, COL_NB]

            for col, label in species:
                if ds.data.shape[1] > col:
                    y = ds.data[:, col]
                    ax.plot(x, y, lw=1.1, label=label)

            ax.set_title(f"csi={ds.csi:.3f}")
            ax.set_xlabel(r"$n_B/n_0$")
            ax.set_ylabel(r"$n_i$ [fm$^{-3}$]")
            ax.grid(alpha=0.25)

        handles, labels = axs[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=5)

        fig.suptitle(f"Population snapshots | {model} | {topology} | B={b_label} G", y=1.03)
        fig.tight_layout()
        fig.savefig(out_dir / f"population_snapshots_{model}_{topology}_B_{b_label}.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)


def find_onset_threshold(nb_array: np.ndarray, pop_array: np.ndarray, threshold_val: float = 1e-6) -> float:
    """Retorna a densidade nb_array onde pop_array ultrapassa o limite (threshold_val)."""
    mask = pop_array > threshold_val
    if not np.any(mask):
        return math.nan
    idx = np.argmax(mask)
    return float(nb_array[idx])


def plot_population_thresholds(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plota a densidade bariônica de surgimento (onset) das partículas em função de csi.
    Foco nas partículas que 'nascem' em altas densidades (Múons e Hiperons).
    """
    ensure_dir(out_dir)
    by_combo: Dict[Tuple[str, str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    # Excluímos n, p, e- porque já existem desde o início (crosta/baixa densidade)
    species_to_track = [
        (COL_NMU, "mu-"),
        (COL_NL0, "Lambda0"),
        (COL_NSM, "Sigma-"),
        (COL_NS0, "Sigma0"),
        (COL_NSP, "Sigma+"),
        (COL_NXM, "Xi-"),
        (COL_NX0, "Xi0"),
    ]

    for key, arr in by_combo.items():
        model, topology, b_label = key
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        if not arr_sorted:
            continue

        x_vals = np.array([d.csi for d in arr_sorted], dtype=float)

        plt.figure(figsize=(8, 6))
        plotted_any = False

        for col, label in species_to_track:
            y_vals = []
            for ds in arr_sorted:
                if ds.data.shape[1] > col:
                    # Avalia o threshold para cada valor de csi
                    onset = find_onset_threshold(ds.data[:, COL_NB], ds.data[:, col])
                    y_vals.append(onset)
                else:
                    y_vals.append(math.nan)

            y_arr = np.array(y_vals)
            mask = np.isfinite(y_arr) & (x_vals > 0.0)

            if np.count_nonzero(mask) > 0:
                plt.plot(x_vals[mask], y_arr[mask], marker='o', ms=4, lw=1.5, label=label)
                plotted_any = True

        if plotted_any:
            plt.xlabel(CSI_LABEL)
            plt.ylabel(r"Onset Density ($n_B$) [fm$^{-3}$]")
            plt.title(f"Particle Onset Thresholds | {model} | {topology} | B={b_label} G")
            plt.grid(alpha=0.3, linestyle='--')
            plt.xscale("log")
            plt.legend(fontsize=9, loc='best')
            plt.tight_layout()
            plt.savefig(out_dir / f"onset_thresholds_{model}_{topology}_B_{b_label}.png", dpi=dpi)
        plt.close()


def _plot_metric_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    topo: str,
    metric_name: str,
    y_label: str,
    use_log_csi: bool,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model and r.topology == topo]
    if not subset:
        return

    bvals = sorted({r.b_value for r in subset})
    plt.figure(figsize=(8, 6))
    x_label = LOG_CSI_LABEL if use_log_csi else CSI_LABEL

    for b in bvals:
        part = sorted([r for r in subset if r.b_value == b], key=lambda x: x.csi)
        if use_log_csi:
            x = np.array([r.log_csi for r in part], dtype=float)
        else:
            x = np.array([r.csi for r in part], dtype=float)
        y = np.array([getattr(r, metric_name) for r in part], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not use_log_csi:
            mask &= x > 0.0

        if np.count_nonzero(mask) < 1:
            continue
        plt.plot(x[mask], y[mask], marker="o", ms=3, lw=1.2, label=f"B={part[0].b_label} G")

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if use_log_csi:
        plt.title(f"{y_label} vs log10(csi) | {model} | {topo}")
    else:
        plt.title(f"{y_label} vs csi | {model} | {topo}")
        plt.xscale("log")

    # Para max_mass e radius_at_max, evita eixo com offset/diferença
    # e mantém escala padrão absoluta.
    if metric_name in {"max_mass_msun", "radius_at_max_km"}:
        ax = plt.gca()
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def plot_trends(rows: Sequence[SummaryRow], out_dir: Path, dpi: int) -> None:
    ensure_dir(out_dir)

    models = sorted({r.model for r in rows})
    topologies = sorted({r.topology for r in rows})

    metrics = [
        ("max_mass_msun", r"$M_{\max}$ [$M_\odot$]", "max_mass_vs_csi", False),
        ("radius_at_max_km", r"$R(M_{\max})$ [km]", "radius_at_max_vs_csi", False),
        ("central_nb_n0", r"$n_c/n_0$", "central_density_vs_logcsi", True),
        ("central_emag_mevfm3", r"$\epsilon_{mag,c}$ [MeV/fm$^3$]", "central_emag_vs_logcsi", True),
        ("cs2_max", r"$\max(c_s^2)$", "cs2max_vs_logcsi", True),
        ("cs2_min", r"$\min(c_s^2)$", "cs2min_vs_logcsi", True),
        ("cs2_p95", r"$P95(c_s^2)$", "cs2p95_vs_logcsi", True),
        ("cs2_frac_superluminal", r"$f(c_s^2>1)$", "cs2frac_superluminal_vs_logcsi", True),
        ("cs2_frac_negative", r"$f(c_s^2<0)$", "cs2frac_negative_vs_logcsi", True),
        ("landau_nu_max", r"$\nu_{max}^{e}$", "landau_numax_vs_logcsi", True),
        ("landau_nu_p95", r"$P95(\nu_{max}^{e})$", "landau_nup95_vs_logcsi", True),
    ]

    for metric_name, y_label, prefix, use_log_csi in metrics:
        for model in models:
            for topo in topologies:
                _plot_metric_for_subset(
                    rows=rows,
                    model=model,
                    topo=topo,
                    metric_name=metric_name,
                    y_label=y_label,
                    use_log_csi=use_log_csi,
                    out_path=out_dir / f"{prefix}_{model}_{topo}.png",
                    dpi=dpi,
                )


def _slope_for_rows(sub: Sequence[SummaryRow], metric: str) -> float:
    x = np.array([r.log_csi for r in sub], dtype=float)
    y = np.array([getattr(r, metric) for r in sub], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return math.nan
    a, _b = np.polyfit(x[mask], y[mask], 1)
    return float(a)


def _group_rows_by_model_topology_b(rows: Sequence[SummaryRow]) -> Dict[Tuple[str, str, str], List[SummaryRow]]:
    grouped: Dict[Tuple[str, str, str], List[SummaryRow]] = {}
    for r in rows:
        grouped.setdefault((r.model, r.topology, r.b_label), []).append(r)
    return grouped


def _append_scope_section(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    models = sorted({r.model for r in rows})
    tops = sorted({r.topology for r in rows})
    lines.append("# NLEM neutron-star analysis report")
    lines.append("")
    lines.append("## Scope")
    lines.append("This report quantifies how $\\log_{10}(\\xi)$ modifies stellar structure using generated EoS/M-R data.")
    lines.append("")
    lines.append(f"- Total datasets analyzed: **{len(rows)}**")
    lines.append(f"- Models: **{', '.join(models)}**")
    lines.append(f"- Topologies: **{', '.join(tops)}**")
    lines.append("")


def _append_cs2_section(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    lines.append("## Causality and Stability Diagnostics")
    lines.append("Analyzing the speed of sound $c_s^2 = dP/d\\epsilon$ to verify physical validity constraints.")
    lines.append("")

    causality_violations = sum(1 for r in rows if r.cs2_max > 1.0)
    stability_violations = sum(1 for r in rows if r.cs2_min < 0.0)
    lines.append(f"- **Causality Violations ($c_s^2 > 1$):** {causality_violations} datasets")
    lines.append(f"- **Stability Violations ($c_s^2 < 0$):** {stability_violations} datasets")
    lines.append("")

    if causality_violations > 0 or stability_violations > 0:
        lines.append("**Warning:** Some datasets exhibit non-physical behavior. High $\\xi$ combined with extreme magnetic fields may lead to $c_s^2 > 1.0$ (superluminal sound speed) or $c_s^2 < 0$ (mechanical instability). Inspect the `cs2_max` and `cs2_min` trend graphs to determine the $\\xi$ cutoff.")
    else:
        lines.append("All analyzed datasets respect the physical limits of causality and thermodynamic stability.")
    lines.append("")

    cs2_frac_neg = np.array([r.cs2_frac_negative for r in rows], dtype=float)
    cs2_frac_sup = np.array([r.cs2_frac_superluminal for r in rows], dtype=float)
    if np.any(np.isfinite(cs2_frac_neg)) and np.any(np.isfinite(cs2_frac_sup)):
        lines.append("Fractional diagnostics over each EoS curve:")
        lines.append(f"- Mean $f(c_s^2<0)$: {float(np.nanmean(cs2_frac_neg)):.6f}")
        lines.append(f"- Mean $f(c_s^2>1)$: {float(np.nanmean(cs2_frac_sup)):.6f}")
        lines.append("")


def _append_landau_section(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    landau_global = np.array([r.landau_nu_max for r in rows], dtype=float)
    if not np.any(np.isfinite(landau_global)):
        return
    lines.append("## Landau Quantization Diagnostics")
    lines.append("Global statistics for electron Landau occupancy:")
    lines.append(f"- Median $\\nu_{{max}}^e$: {float(np.nanmedian(landau_global)):.3f}")
    lines.append(f"- P95 $\\nu_{{max}}^e$: {float(np.nanpercentile(landau_global, 95.0)):.3f}")
    lines.append("")


def _append_group_trend_section(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    lines.append("## Trend diagnostics by $(model, topology, B)$")
    lines.append("")
    lines.append("Linear slopes (first-order sensitivity):")
    lines.append("- $dM_{max}/d\\log_{10}(\\xi)$")
    lines.append("- $dR(M_{max})/d\\log_{10}(\\xi)$")
    lines.append("")

    grouped = _group_rows_by_model_topology_b(rows)
    for (model, topo, b_label), sub in sorted(
        grouped.items(),
        key=lambda kv: (
            kv[0][0],
            kv[0][1],
            _safe_float(kv[0][2]) if _safe_float(kv[0][2]) is not None else math.inf,
        ),
    ):
        sub = sorted(sub, key=lambda x: x.log_csi)
        s_mass = _slope_for_rows(sub, "max_mass_msun")
        s_rad = _slope_for_rows(sub, "radius_at_max_km")

        mvals = np.array([r.max_mass_msun for r in sub], dtype=float)
        rvals = np.array([r.radius_at_max_km for r in sub], dtype=float)
        cmask = np.isfinite(mvals) & np.isfinite(rvals)

        m_min = float(np.nanmin(mvals)) if np.any(np.isfinite(mvals)) else math.nan
        m_max = float(np.nanmax(mvals)) if np.any(np.isfinite(mvals)) else math.nan
        r_min = float(np.nanmin(rvals)) if np.any(np.isfinite(rvals)) else math.nan
        r_max = float(np.nanmax(rvals)) if np.any(np.isfinite(rvals)) else math.nan

        lines.append(f"### {model} | {topo} | B={b_label} G")
        lines.append(f"- Samples: {len(sub)}")
        lines.append(f"- $M_{{max}}$ range: {m_min:.4f} to {m_max:.4f} $M_\\odot$")
        lines.append(f"- $R(M_{{max}})$ range: {r_min:.4f} to {r_max:.4f} km")
        lines.append(f"- Slope $dM_{{max}}/d\\log_{{10}}(\\xi)$: {s_mass:.6f}")
        lines.append(f"- Slope $dR(M_{{max}})/d\\log_{{10}}(\\xi)$: {s_rad:.6f}")
        if np.any(cmask):
            lines.append(
                f"- Median pair $(M_{{max}}, R)$: ({float(np.nanmedian(mvals)):.4f}, {float(np.nanmedian(rvals)):.4f})"
            )
        lines.append("")


def write_scientific_report(rows: Sequence[SummaryRow], out_file: Path) -> None:
    ensure_dir(out_file.parent)
    lines: List[str] = []
    _append_scope_section(lines, rows)
    _append_cs2_section(lines, rows)
    _append_landau_section(lines, rows)
    _append_group_trend_section(lines, rows)

    lines.append("## Interpretation guide")
    lines.append("- Positive $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ tends to stiffen the effective sequence in your setup.")
    lines.append("- Negative $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ softens the sequence.")
    lines.append("- Compare isotropic vs anisotropic at fixed $(model, B)$ to isolate topology effects.")
    lines.append("- Use population thresholds to map how $\\xi$ delays or anticipates hyperon onset.")

    out_file.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_root)

    datasets = discover_datasets(args.input_root)
    if not datasets:
        print(f"[ERRO] Nenhum eos.dat válido encontrado em: {args.input_root}")
        return

    print(f"[INFO] Datasets encontrados: {len(datasets)}")

    summary_rows = [summarize_dataset(ds) for ds in datasets]

    summary_csv = args.output_root / "summary_all_datasets.csv"
    write_summary_csv(summary_rows, summary_csv)
    print(f"[OK] CSV consolidado: {summary_csv}")

    group_stats_csv = args.output_root / "summary_group_stats.csv"
    write_group_stats_csv(summary_rows, group_stats_csv)
    print(f"[OK] Estatísticas agregadas por grupo: {group_stats_csv}")

    trends_dir = args.output_root / "trends"
    plot_trends(summary_rows, trends_dir, dpi=args.dpi)
    print(f"[OK] Tendências globais: {trends_dir}")

    families_dir = args.output_root / "families"
    plot_family_eos_mr(
        datasets,
        families_dir,
        max_curves=max(1, args.max_curves_per_family),
        dpi=args.dpi,
    )
    print(f"[OK] Famílias EoS/MR e Velocidade do Som (cs2): {families_dir}")

    comparison_dir = args.output_root / "comparison"
    plot_mr_specific_targets(
        datasets,
        comparison_dir,
        baseline_root=args.baseline_root,
        targets=COMPARISON_MR_TARGETS,
        dpi=args.dpi,
    )
    print(f"[OK] Comparação M-R para valores específicos: {comparison_dir / 'mr_specific'}")

    pops_dir = args.output_root / "populations"
    plot_population_snapshots(datasets, pops_dir, dpi=args.dpi)
    print(f"[OK] Snapshots de populações: {pops_dir}")

    # --- NOVO: Geração dos gráficos de limiar de surgimento (Thresholds) ---
    thresholds_dir = args.output_root / "thresholds"
    plot_population_thresholds(datasets, thresholds_dir, dpi=args.dpi)
    print(f"[OK] Gráficos de Thresholds de Partículas: {thresholds_dir}")

    print(
        "[SKIP] Quantização de Landau indisponível: "
        "o esquema EOS atual não exporta B(n)."
    )

    report_md = args.output_root / "scientific_report.md"
    write_scientific_report(summary_rows, report_md)
    print(f"[OK] Relatório científico: {report_md}")

    print("[DONE] Análise completa finalizada.")


if __name__ == "__main__":
    main()
