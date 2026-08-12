#!/usr/bin/env python3
"""
Análise científica completa dos dados NLEM gerados em output/modmax.

Objetivos:
- Consolidar EoS, M-R e populações para GM1/GM3.
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
import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

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

# Populações
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

# M-R columns appended by write_eos_with_mr
EOS_RESULTS_SIZE = 34
COL_MR_MASS = EOS_RESULTS_SIZE
COL_MR_RADIUS = EOS_RESULTS_SIZE + 1
COL_MR_BARYONIC = EOS_RESULTS_SIZE + 2

# Fator padrão. Se necessário, o código também detecta automaticamente
# se o mu_e parece normalizado e converte por 939.0 para MeV.
MUE_TO_MEV_FACTOR = 1.0

LOG_CSI_LABEL = r"$\log_{10}(\gamma)$"
CSI_LABEL = r"$\gamma$"
HBARC_MEV_FM = 197.3269804
MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0
COMPARISON_MR_TARGETS = [1e-22, 1e-20, 1e-15, 1e-10, 1e-5, 1e-1]
MR_XLIM_BY_MODEL = {
    "GM3": (9.0, 15.0),
}
DEFAULT_MR_XLIM = (10.0, 16.5)
PARTICLE_SPECIES = [
    (COL_NE, r"$e^-$"),
    (COL_NMU, r"$\mu^-$"),
    (COL_NN, r"$n$"),
    (COL_NP, r"$p$"),
    (COL_NL0, r"$\Lambda^0$"),
    (COL_NSM, r"$\Sigma^-$"),
    (COL_NS0, r"$\Sigma^0$"),
    (COL_NSP, r"$\Sigma^+$"),
    (COL_NXM, r"$\Xi^-$"),
    (COL_NX0, r"$\Xi^0$"),
]
THRESHOLD_SPECIES = [
    (COL_NMU, r"$\mu^-$"),
    (COL_NL0, r"$\Lambda^0$"),
    (COL_NSM, r"$\Sigma^-$"),
    (COL_NS0, r"$\Sigma^0$"),
    (COL_NSP, r"$\Sigma^+$"),
    (COL_NXM, r"$\Xi^-$"),
    (COL_NX0, r"$\Xi^0$"),
]


@dataclass
class Dataset:
    model: str
    b_label: str
    b_value: float
    csi: float
    log_csi: float
    is_baseline: bool
    file_path: Path
    data: np.ndarray


@dataclass
class SummaryRow:
    model: str
    b_label: str
    b_value: float
    csi: float
    log_csi: float
    is_baseline: bool
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
        default=Path("output/modmax"),
        help="Raiz dos dados baseline eos_baseline.dat",
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
    except (TypeError, ValueError):
        return None


def _latex_sci(value: float, precision: int = 2) -> str:
    """Format numbers for plot labels as mantissa times power of ten."""
    if not np.isfinite(value):
        return r"\mathrm{nan}"
    if value == 0.0:
        return "0"

    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10.0 ** exponent)
    return rf"{mantissa:.{precision}f}\times10^{{{exponent}}}"


def _b_plot_label(b_value: float, *, precision: int = 2) -> str:
    return rf"$B={_latex_sci(b_value, precision)}\,\mathrm{{G}}$"


def _gamma_plot_label(csi: float, *, precision: int = 2) -> str:
    return rf"$\gamma={_latex_sci(csi, precision)}$"


def _has_columns(arr: np.ndarray, last_col: int) -> bool:
    return arr.ndim == 2 and arr.shape[1] > last_col


def parse_metadata_from_path(path: Path) -> Optional[Tuple[str, str, float, float, float, bool]]:
    """
    Extract metadata from the actual output/modmax layout:

        output/modmax/GM1/B_1.00e16/eos_csi_1.00e-3.dat
        output/modmax/GM1/B_1.00e16/eos_baseline.dat

    The parser intentionally uses only Path components and the filename. It
    does not inspect or regex-match the absolute path.
    """
    if path.suffix != ".dat":
        return None
    if path.parent.parent == path.parent:
        return None

    model = path.parent.parent.name
    b_dir = path.parent.name
    if not model.startswith("GM") or not b_dir.startswith("B_"):
        return None

    b_label = b_dir.removeprefix("B_")
    b_value = _safe_float(b_label)
    if b_value is None:
        return None

    is_baseline = path.name == "eos_baseline.dat"
    if is_baseline:
        csi = math.nan
        log_csi = math.nan
        return model, b_label, b_value, csi, log_csi, True

    prefix = "eos_csi_"
    if not (path.name.startswith(prefix) and path.suffix == ".dat"):
        return None

    gamma_text = path.stem.removeprefix(prefix)
    csi = _safe_float(gamma_text)
    if csi is None:
        return None

    log_csi = math.log10(csi) if csi > 0.0 else math.nan
    return model, b_label, b_value, csi, log_csi, False


def load_eos(path: Path) -> Optional[np.ndarray]:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            arr = np.genfromtxt(path, comments="#", invalid_raise=False)
    except (OSError, ValueError):
        return None

    if arr is None or np.size(arr) == 0:
        return None

    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 0:
        return None

    if arr.ndim == 1:
        if arr.size == 0:
            return None
        arr = arr.reshape(1, -1)

    if arr.ndim != 2 or arr.shape[1] <= COL_MU_TOTAL:
        return None

    return arr


def discover_datasets(input_root: Path) -> List[Dataset]:
    datasets: List[Dataset] = []

    for data_file in input_root.rglob("*.dat"):
        if not data_file.is_file():
            continue

        meta = parse_metadata_from_path(data_file)
        if meta is None:
            continue

        loaded = load_eos(data_file)
        if loaded is None:
            continue

        model, b_label, b_value, csi, log_csi, is_baseline = meta
        datasets.append(
            Dataset(
                model=model,
                b_label=b_label,
                b_value=b_value,
                csi=csi,
                log_csi=log_csi,
                is_baseline=is_baseline,
                file_path=data_file,
                data=loaded,
            )
        )

    datasets.sort(key=lambda d: (d.model, d.b_value, d.is_baseline, d.csi if np.isfinite(d.csi) else math.inf))
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


def mr_xlim_for_model(model: str) -> Tuple[float, float]:
    return MR_XLIM_BY_MODEL.get(model, DEFAULT_MR_XLIM)


def _load_mr_from_eos(path: Path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    arr = load_eos(path)
    if arr is None or not _has_columns(arr, COL_MR_RADIUS):
        return None

    mass = arr[:, COL_MR_MASS].astype(float)
    radius = arr[:, COL_MR_RADIUS].astype(float)
    mask = valid_mr_mask(mass, radius)
    if not np.any(mask):
        return None

    return radius[mask], mass[mask]

def _load_baseline_curve(
    baseline_root: Path,
    model: str,
    b_label: str,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Carrega a curva eos_baseline.dat correspondente ao par (modelo, B).
    """
    model_dir = baseline_root / model
    if not model_dir.exists() or not model_dir.is_dir():
        return None

    b_dir_name = f"B_{b_label}"
    exact_path = model_dir / b_dir_name / "eos_baseline.dat"
    if exact_path.exists():
        curve = _load_mr_from_eos(exact_path)
        if curve is not None:
            return curve

    for eos_path in sorted(model_dir.rglob("eos_baseline.dat")):
        if eos_path.parent.name != b_dir_name:
            continue
        curve = _load_mr_from_eos(eos_path)
        if curve is not None:
            return curve

    return None


def summarize_dataset(ds: Dataset) -> SummaryRow:
    arr = ds.data

    if not _has_columns(arr, COL_MR_RADIUS):
        return SummaryRow(
            model=ds.model,
            b_label=ds.b_label,
            b_value=ds.b_value,
            csi=ds.csi,
            log_csi=ds.log_csi,
            is_baseline=ds.is_baseline,
            n_eos_points=int(arr.shape[0]) if arr.ndim == 2 else 0,
            n_mr_points=0,
            max_mass_msun=math.nan,
            radius_at_max_km=math.nan,
            central_nb_n0=math.nan,
            central_eps_mevfm3=math.nan,
            central_p_mevfm3=math.nan,
            central_meff=math.nan,
            central_emag_mevfm3=math.nan,
            cs2_min=math.nan,
            cs2_max=math.nan,
            cs2_p95=math.nan,
            cs2_frac_negative=math.nan,
            cs2_frac_superluminal=math.nan,
            landau_nu_max=math.nan,
            landau_nu_p95=math.nan,
            eos_path=ds.file_path.as_posix(),
        )

    mass_col = arr[:, COL_MR_MASS]
    radius_col = arr[:, COL_MR_RADIUS]
    mr_mask = valid_mr_mask(mass_col, radius_col)

    n_mr_points = int(np.count_nonzero(mr_mask))
    n_eos_points = arr.shape[0]

    if n_mr_points > 0:
        mr_indices = np.nonzero(mr_mask)[0]
        local = np.argmax(mass_col[mr_mask])
        i_max = int(mr_indices[local])

        max_mass = float(arr[i_max, COL_MR_MASS])
        r_at_max = float(arr[i_max, COL_MR_RADIUS])
        central_nb = float(arr[i_max, COL_NB])
        central_eps = float(arr[i_max, COL_EPS])
        central_p = float(arr[i_max, COL_P])
        central_meff = float(arr[i_max, COL_MEFF]) if _has_columns(arr, COL_MEFF) else math.nan
        central_emag = float(arr[i_max, COL_EMAG]) if _has_columns(arr, COL_EMAG) else math.nan
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

    # The current EOS schema does not export local B(n). Column 20 is
    # mu_total/M_N, so Landau diagnostics cannot be reconstructed safely.
    landau_nu_max = math.nan
    landau_nu_p95 = math.nan

    return SummaryRow(
        model=ds.model,
        b_label=ds.b_label,
        b_value=ds.b_value,
        csi=ds.csi,
        log_csi=ds.log_csi,
        is_baseline=ds.is_baseline,
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
            "model", "b_label", "b_value", "csi", "log_csi", "is_baseline",
            "n_eos_points", "n_mr_points", "max_mass_msun", "radius_at_max_km",
            "central_nb_n0", "central_eps_mevfm3", "central_p_mevfm3", "central_meff",
            "central_emag_mevfm3", "cs2_min", "cs2_max", "cs2_p95",
            "cs2_frac_negative", "cs2_frac_superluminal",
            "landau_nu_max", "landau_nu_p95", "eos_path",
        ])
        for r in rows:
            w.writerow([
                r.model, r.b_label, r.b_value, r.csi, r.log_csi, r.is_baseline,
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
    Estatísticas agregadas por (model, B), incluindo sensibilidade
    a log10(csi) e qualidade de ajuste linear (R²).
    """
    ensure_dir(out_csv.parent)

    grouped: Dict[Tuple[str, str], List[SummaryRow]] = {}
    for r in rows:
        if r.is_baseline or not np.isfinite(r.log_csi):
            continue
        grouped.setdefault((r.model, r.b_label), []).append(r)

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "b_label", "n_samples",
            "mass_slope", "mass_r2", "radius_slope", "radius_r2",
            "cs2max_slope", "cs2max_r2", "landau_max_slope", "landau_max_r2",
            "mean_cs2_frac_negative", "mean_cs2_frac_superluminal",
            "median_landau_nu_max", "median_landau_nu_p95",
        ])

        for (model, b_label), sub in sorted(
            grouped.items(),
            key=lambda kv: (
                kv[0][0],
                _safe_float(kv[0][1]) if _safe_float(kv[0][1]) is not None else math.inf,
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


def group_key(ds: Dataset) -> Tuple[str, str]:
    return ds.model, ds.b_label


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
    em função da densidade bariônica para cada família (model, B).
    """
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, b_label = key
        b_value = _safe_float(b_label)
        b_text = _b_plot_label(b_value) if b_value is not None else rf"$B={b_label}\,\mathrm{{G}}$"
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        plot_set = [d for d in arr_sorted if not d.is_baseline and np.isfinite(d.csi)]
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
            plt.title(f"Quantização de Landau | {model} | {b_text}")
            plt.grid(alpha=0.3)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
            cbar = plt.colorbar(sm, ax=plt.gca())
            cbar.set_label(CSI_LABEL)
            plt.tight_layout()
            plt.savefig(out_dir / f"landau_levels_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()


def plot_family_eos_mr(
    datasets: Sequence[Dataset],
    out_dir: Path,
    max_curves: int,
    dpi: int,
) -> None:
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, b_label = key
        b_value = _safe_float(b_label)
        b_text = _b_plot_label(b_value) if b_value is not None else rf"$B={b_label}\,\mathrm{{G}}$"
        arr_sorted = sorted(arr, key=lambda d: d.csi)
        baseline_set = [d for d in arr_sorted if d.is_baseline]
        csi_set = [d for d in arr_sorted if not d.is_baseline and np.isfinite(d.csi)]
        plot_set = baseline_set[:1] + _subsample_sorted(csi_set, max_curves)
        cmap_eos = plt.get_cmap("viridis")
        cmap_mr = plt.get_cmap("plasma")
        cmap_cs2 = plt.get_cmap("cool")

        cvals = np.array([d.csi for d in csi_set], dtype=float)
        finite_cvals = cvals[np.isfinite(cvals) & (cvals > 0.0)]
        if finite_cvals.size == 0:
            continue
        cmin, cmax = float(np.min(finite_cvals)), float(np.max(finite_cvals))
        vmax = cmax if cmax > cmin else cmin * 10.0
        norm = mcolors.LogNorm(vmin=max(cmin, 1e-300), vmax=vmax)

        def _color_value(value: float):
            return norm(value if value > 0.0 else cmin)

        # --- EoS family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            color = "black" if d.is_baseline else cmap_eos(_color_value(d.csi))
            eps = d.data[:, COL_EPS]
            p = d.data[:, COL_P]
            mask = np.isfinite(eps) & np.isfinite(p) & (eps > 0.0) & (p > 0.0)
            if np.count_nonzero(mask) < 2:
                continue
            eps_fm4 = eps[mask] / HBARC_MEV_FM
            p_fm4 = p[mask] / HBARC_MEV_FM
            plt.plot(eps_fm4, p_fm4, color=color, alpha=0.55, lw=2.2)

        plt.xlabel(r"$\epsilon$ [fm$^{-4}$]")
        plt.ylabel(r"$P$ [fm$^{-4}$]")
        plt.title(f"EoS family | {model} | {b_text}")
        plt.xlim(0, 7)
        plt.ylim(0, 4)
        plt.grid(alpha=0.25)

        sm = plt.cm.ScalarMappable(cmap=cmap_eos, norm=norm)
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(LOG_CSI_LABEL)
        plt.tight_layout()
        plt.savefig(out_dir / f"eos_family_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()

        # --- M-R family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            if not _has_columns(d.data, COL_MR_RADIUS):
                continue

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
            color = "black" if d.is_baseline else cmap_mr(_color_value(d.csi))
            plt.plot(r_valid[sort_idx], m_valid[sort_idx], color=color, alpha=0.55, lw=2.2)

        plt.xlabel("Radius [km]")
        plt.ylabel(r"Mass [$M_\odot$]")
        plt.title(f"M-R family | {model} | {b_text}")
        plt.xlim(*mr_xlim_for_model(model))
        plt.ylim(1.0, 2.25)
        plt.grid(alpha=0.25)

        sm_mr = plt.cm.ScalarMappable(cmap=cmap_mr, norm=norm)
        cbar_mr = plt.colorbar(sm_mr, ax=plt.gca())
        cbar_mr.set_label(LOG_CSI_LABEL)
        plt.tight_layout()
        plt.savefig(out_dir / f"mr_family_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()

        # --- Speed of Sound (c_s^2) family ---
        plt.figure(figsize=(8, 6))
        for d in plot_set:
            color = "black" if d.is_baseline else cmap_cs2(_color_value(d.csi))
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

                    plt.plot(eps_sorted[:-1][good], cs2[good], color=color, alpha=0.55, lw=2.0)

        plt.axhline(1.0, color='red', linestyle='--', alpha=0.7, label='Causalidade ($c_s^2 = 1$)')
        plt.axhline(0.0, color='black', linestyle='--', alpha=0.7, label='Estabilidade ($c_s^2 = 0$)')
        plt.axhline(1/3, color='gray', linestyle=':', alpha=0.7, label='Limite Conforme ($c_s^2 = 1/3$)')

        plt.xlabel(r"$\epsilon$ [MeV/fm$^3$]")
        plt.ylabel(r"$c_s^2$ (Unidades de $c=1$)")
        plt.title(f"Speed of Sound | {model} | {b_text}")
        plt.grid(alpha=0.25)
        plt.legend(loc='upper right', fontsize=9)

        sm_cs2 = plt.cm.ScalarMappable(cmap=cmap_cs2, norm=norm)
        cbar_cs2 = plt.colorbar(sm_cs2, ax=plt.gca())
        cbar_cs2.set_label(LOG_CSI_LABEL)

        plt.tight_layout()
        plt.savefig(out_dir / f"cs2_family_{model}_B_{b_label}.png", dpi=dpi)
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

    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, b_label = key
        b_value = _safe_float(b_label)
        b_text = _b_plot_label(b_value) if b_value is not None else rf"$B={b_label}\,\mathrm{{G}}$"
        arr_sorted = sorted(
            [d for d in arr if not d.is_baseline and np.isfinite(d.csi)],
            key=lambda d: d.csi,
        )
        if not arr_sorted:
            continue

        baseline_mr = _load_baseline_curve(baseline_root, model, b_label)
        # initialize zoom center variables (will be set if baseline exists)
        _zoom_center_r = None
        _zoom_center_m = None

        selected = _select_nearest_mr_targets(arr_sorted, targets)
        if not selected:
            selected = []

        plt.figure(figsize=(8, 6))
        plotted_curves: List[Tuple[np.ndarray, np.ndarray, str, Tuple[float, float, float, float]]] = []

        # prepare colormap for gamma-based coloring (like family plots)
        csi_values = [d.csi for d, _ in selected if np.isfinite(d.csi) and d.csi > 0.0]
        if csi_values:
            cmin = min(csi_values)
            cmax = max(csi_values)
        else:
            cmin = 1e-25
            cmax = 1e-1

        if cmin > 0:
            vmax = cmax if cmax > cmin else cmin * 10.0
            norm = mcolors.LogNorm(vmin=max(cmin, 1e-300), vmax=vmax)
        else:
            norm = mcolors.Normalize(vmin=0, vmax=1)
        cmap_colors = plt.get_cmap("viridis")

        if baseline_mr is not None:
            baseline_radius, baseline_mass = baseline_mr
            baseline_label = r"$B=0\,\mathrm{G}$"
            baseline_color = (0.0, 0.0, 0.0, 0.95)
            plt.plot(
                baseline_radius,
                baseline_mass,
                color=baseline_color,
                alpha=0.6,
                lw=2.4,
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
            if not _has_columns(d.data, COL_MR_RADIUS):
                continue

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
            curve_label = _gamma_plot_label(d.csi)
            plt.plot(
                r_valid[sort_idx],
                m_valid[sort_idx],
                color=color,
                alpha=0.6,
                lw=2.2,
                linestyle="--",
                label=curve_label,
            )
            plotted_curves.append((r_valid[sort_idx], m_valid[sort_idx], curve_label, color))

        plt.xlabel("Radius [km]")
        plt.ylabel(r"Mass [$M_\odot$]")
        plt.title(f"MR comparison | {model} | {b_text}")
        plt.xlim(*mr_xlim_for_model(model))
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
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x_span = xlim[1] - xlim[0]
            y_span = ylim[1] - ylim[0]
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
        plt.subplots_adjust(left=0.12, right=0.96, bottom=0.12, top=0.92)
        plt.savefig(comparison_dir / f"mr_specific_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()


def plot_population_snapshots(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    ensure_dir(out_dir)
    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        model, b_label = key
        b_value = _safe_float(b_label)
        b_text = _b_plot_label(b_value) if b_value is not None else rf"$B={b_label}\,\mathrm{{G}}$"
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

            for col, label in PARTICLE_SPECIES:
                if _has_columns(ds.data, col):
                    y = ds.data[:, col]
                    ax.plot(x, y, lw=2.0, alpha=0.6, label=label)

            ax.set_title("baseline" if ds.is_baseline else _gamma_plot_label(ds.csi))
            ax.set_xlabel(r"$n_B/n_0$")
            ax.set_ylabel(r"$n_i$ [fm$^{-3}$]")
            ax.grid(alpha=0.25)

        handles, labels = axs[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=5)

        fig.suptitle(f"Population snapshots | {model} | {b_text}", y=1.03)
        fig.tight_layout()
        fig.savefig(out_dir / f"population_snapshots_{model}_B_{b_label}.png", dpi=dpi, bbox_inches="tight")
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
    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    # Excluímos n, p, e- porque já existem desde o início (crosta/baixa densidade)
    for key, arr in by_combo.items():
        model, b_label = key
        b_value = _safe_float(b_label)
        b_text = _b_plot_label(b_value) if b_value is not None else rf"$B={b_label}\,\mathrm{{G}}$"
        arr_sorted = sorted(
            [d for d in arr if not d.is_baseline and np.isfinite(d.csi)],
            key=lambda d: d.csi,
        )
        if not arr_sorted:
            continue

        x_vals = np.array([d.csi for d in arr_sorted], dtype=float)

        plt.figure(figsize=(8, 6))
        plotted_any = False

        for col, label in THRESHOLD_SPECIES:
            y_vals = []
            for ds in arr_sorted:
                if _has_columns(ds.data, col):
                    # Avalia o threshold para cada valor de csi
                    onset = find_onset_threshold(ds.data[:, COL_NB], ds.data[:, col])
                    y_vals.append(onset)
                else:
                    y_vals.append(math.nan)

            y_arr = np.array(y_vals)
            mask = np.isfinite(y_arr) & (x_vals > 0.0)

            if np.count_nonzero(mask) > 0:
                plt.plot(x_vals[mask], y_arr[mask], marker='o', ms=4, lw=2.2, alpha=0.65, label=label)
                plotted_any = True

        if plotted_any:
            plt.xlabel(CSI_LABEL)
            plt.ylabel(r"Onset Density ($n_B$) [fm$^{-3}$]")
            plt.title(f"Particle Onset Thresholds | {model} | {b_text}")
            plt.grid(alpha=0.3, linestyle='--')
            plt.xscale("log")
            plt.legend(fontsize=9, loc='best')
            plt.tight_layout()
            plt.savefig(out_dir / f"onset_thresholds_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()


def _plot_metric_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    metric_name: str,
    y_label: str,
    use_log_csi: bool,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model and not r.is_baseline]
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
        plt.plot(
            x[mask],
            y[mask],
            marker="o",
            ms=3,
            lw=2.0,
            alpha=0.65,
            label=_b_plot_label(part[0].b_value),
        )

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if use_log_csi:
        plt.title(f"{y_label} vs log10(csi) | {model}")
    else:
        plt.title(f"{y_label} vs csi | {model}")
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


def _plot_delta_from_first_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    metric_name: str,
    y_label: str,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model and not r.is_baseline]
    if not subset:
        return

    bvals = sorted({r.b_value for r in subset})
    plt.figure(figsize=(8, 6))

    for b in bvals:
        part = sorted([r for r in subset if r.b_value == b], key=lambda x: x.csi)
        x = np.array([r.csi for r in part], dtype=float)
        y = np.array([getattr(r, metric_name) for r in part], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0)
        if np.count_nonzero(mask) < 1:
            continue

        x_valid = x[mask]
        y_valid = y[mask]
        y0 = y_valid[0]
        y_delta = y_valid - y0
        plt.plot(
            x_valid,
            y_delta,
            marker="o",
            ms=3,
            lw=2.0,
            alpha=0.65,
            label=_b_plot_label(part[0].b_value),
        )

    plt.axhline(0.0, color="black", lw=1.0, alpha=0.5)
    plt.xscale("log")
    plt.xlabel(CSI_LABEL)
    plt.ylabel(rf"$\Delta$({y_label})")
    plt.title(f"Delta from first point vs csi | {model}")
    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def _plot_derivative_vs_log_csi_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    metric_name: str,
    y_label: str,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model and not r.is_baseline]
    if not subset:
        return

    bvals = sorted({r.b_value for r in subset})
    plt.figure(figsize=(8, 6))

    for b in bvals:
        part = sorted([r for r in subset if r.b_value == b], key=lambda x: x.log_csi)
        x = np.array([r.log_csi for r in part], dtype=float)
        y = np.array([getattr(r, metric_name) for r in part], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if np.count_nonzero(mask) < 2:
            continue

        x_valid = x[mask]
        y_valid = y[mask]
        sort_idx = np.argsort(x_valid)
        x_valid = x_valid[sort_idx]
        y_valid = y_valid[sort_idx]

        unique_mask = np.concatenate(([True], np.diff(x_valid) > 0.0))
        x_valid = x_valid[unique_mask]
        y_valid = y_valid[unique_mask]

        if x_valid.size < 2:
            continue

        dy_dlogcsi = np.gradient(y_valid, x_valid)
        plt.plot(
            x_valid,
            dy_dlogcsi,
            marker="o",
            ms=3,
            lw=2.0,
            alpha=0.65,
            label=_b_plot_label(b),
        )

    plt.axhline(0.0, color="black", lw=1.0, alpha=0.5)
    plt.xlabel(LOG_CSI_LABEL)
    plt.ylabel(y_label)
    plt.title(f"{y_label} vs log10(gamma) | {model}")
    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def _plot_slope_vs_b_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    metric_name: str,
    y_label: str,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model and not r.is_baseline]
    if not subset:
        return

    bvals = sorted({r.b_value for r in subset if np.isfinite(r.b_value) and r.b_value > 0.0})
    x_vals: List[float] = []
    y_vals: List[float] = []

    for b in bvals:
        sub = sorted([r for r in subset if r.b_value == b], key=lambda x: x.log_csi)
        slope = _slope_for_rows(sub, metric_name)
        if np.isfinite(slope):
            x_vals.append(b)
            y_vals.append(slope)

    if not x_vals:
        return

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, y_vals, marker="o", ms=4, lw=2.2, alpha=0.65, color="tab:blue")
    plt.axhline(0.0, color="black", lw=1.0, alpha=0.5)
    plt.xscale("log")
    plt.xlabel(r"$B$ [G]")
    plt.ylabel(y_label)
    plt.title(f"{y_label} vs B | {model}")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def plot_trends(rows: Sequence[SummaryRow], out_dir: Path, dpi: int) -> None:
    ensure_dir(out_dir)

    models = sorted({r.model for r in rows})

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
            _plot_metric_for_subset(
                rows=rows,
                model=model,
                metric_name=metric_name,
                y_label=y_label,
                use_log_csi=use_log_csi,
                out_path=out_dir / f"{prefix}_{model}.png",
                dpi=dpi,
            )

            if metric_name in {"max_mass_msun", "radius_at_max_km"} and not use_log_csi:
                _plot_delta_from_first_for_subset(
                    rows=rows,
                    model=model,
                    metric_name=metric_name,
                    y_label=y_label,
                    out_path=out_dir / f"{prefix}_delta_from_first_{model}.png",
                    dpi=dpi,
                )

    derivative_metrics = [
        ("max_mass_msun", r"$dM_{max}/d\log_{10}(\gamma)$ [$M_\odot$]", "dmmax_dloggamma_vs_loggamma"),
        ("radius_at_max_km", r"$dR(M_{max})/d\log_{10}(\gamma)$ [km]", "drmax_dloggamma_vs_loggamma"),
    ]

    for metric_name, y_label, prefix in derivative_metrics:
        for model in models:
            _plot_derivative_vs_log_csi_for_subset(
                rows=rows,
                model=model,
                metric_name=metric_name,
                y_label=y_label,
                out_path=out_dir / f"{prefix}_{model}.png",
                dpi=dpi,
            )


def _slope_for_rows(sub: Sequence[SummaryRow], metric: str) -> float:
    usable = [r for r in sub if not r.is_baseline]
    x = np.array([r.log_csi for r in usable], dtype=float)
    y = np.array([getattr(r, metric) for r in usable], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return math.nan
    a, _b = np.polyfit(x[mask], y[mask], 1)
    return float(a)


def _group_rows_by_model_b(rows: Sequence[SummaryRow]) -> Dict[Tuple[str, str], List[SummaryRow]]:
    grouped: Dict[Tuple[str, str], List[SummaryRow]] = {}
    for r in rows:
        if r.is_baseline:
            continue
        grouped.setdefault((r.model, r.b_label), []).append(r)
    return grouped


def _append_scope_section(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    models = sorted({r.model for r in rows})
    lines.append("# NLEM neutron-star analysis report")
    lines.append("")
    lines.append("## Scope")
    lines.append("This report quantifies how $\\log_{10}(\\xi)$ modifies stellar structure using generated EoS/M-R data.")
    lines.append("")
    lines.append(f"- Total datasets analyzed: **{len(rows)}**")
    lines.append(f"- Models: **{', '.join(models)}**")
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
    lines.append("## Trend diagnostics by $(model, B)$")
    lines.append("")
    lines.append("First derivative diagnostics plotted against $\\log_{10}(\\gamma)$:")
    lines.append("- $dM_{max}/d\\log_{10}(\\gamma)$")
    lines.append("- $dR(M_{max})/d\\log_{10}(\\gamma)$")
    lines.append("")

    grouped = _group_rows_by_model_b(rows)
    for (model, b_label), sub in sorted(
        grouped.items(),
        key=lambda kv: (
            kv[0][0],
            _safe_float(kv[0][1]) if _safe_float(kv[0][1]) is not None else math.inf,
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

        lines.append(f"### {model} | B={b_label} G")
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
    lines.append("- Use population thresholds to map how $\\xi$ delays or anticipates hyperon onset.")

    out_file.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_root)

    datasets = discover_datasets(args.input_root)
    if not datasets:
        print(f"[ERRO] Nenhum eos_csi_*.dat ou eos_baseline.dat válido encontrado em: {args.input_root}")
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
