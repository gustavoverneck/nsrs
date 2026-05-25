#!/usr/bin/env python3
"""
Análise científica completa dos dados NLEM gerados em output/nlem_log.

Objetivos:
- Consolidar EoS, M-R e populações para GM1/GM3.
- Quantificar o efeito de log10(csi) na estrutura estelar.
- Quantificar como o campo magnético efetivo H varia com csi no modelo Log.
- Avaliar os limites de causalidade e estabilidade através da velocidade do som (c_s^2).
- Mapear os limiares de surgimento (onset thresholds) das partículas.
- Mapear a quantização de Landau para elétrons via mu_e e B(n).
- Analisar a variação dos potenciais químicos com a densidade bariônica.
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
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


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
COL_BDD = 20

# Fator padrão. Se necessário, o código também detecta automaticamente
# se o mu_e parece normalizado e converte por 939.0 para MeV.
MUE_TO_MEV_FACTOR = 1.0

LOG_CSI_LABEL = r"$\log_{10}(\xi)$"
NB_DENSITY_LABEL = r"$n_B/n_0$"
POPULATION_DENSITY_LABEL = r"$Y_i$ [%]"
MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0
MR_PLOT_MIN_RADIUS_KM = 8.0
MR_PLOT_MAX_RADIUS_KM = 15.0
MR_PLOT_MIN_MASS_MSUN = 1.0
MR_PLOT_MAX_MASS_MSUN = 2.5
H_COLOR = "tab:blue"
HB_RATIO_COLOR = "tab:orange"
LEGEND_LOC_UPPER_RIGHT = "upper right"
LEGEND_LOC_UPPER_LEFT = "upper left"
G_C2_KM_PER_MSUN = 1.4766
TARGET_CANONICAL_MASS_MSUN = 1.4
HYPERON_COLUMNS = (COL_NL0, COL_NSM, COL_NS0, COL_NSP, COL_NXM, COL_NX0)
OBSERVATIONAL_CONSTRAINTS_PATH = Path("input/observations/neutron_star_constraints.csv")


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
    radius_at_1p4_msun_km: float
    compactness_at_max_mass: float
    central_nb_n0: float
    central_eps_mevfm3: float
    central_p_mevfm3: float
    central_meff: float
    central_emag_mevfm3: float
    hyperon_onset_n0: float
    cs2_min: float
    cs2_max: float
    eos_path: str


@dataclass
class ObservationConstraint:
    kind: str
    label: str
    source: str
    mass_msun: float
    mass_err_minus: float
    mass_err_plus: float
    radius_km: float
    radius_err_minus: float
    radius_err_plus: float
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Análise NLEM completa para dados em output/nlem_log"
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("output/nlem_log"),
        help="Pasta raiz dos dados gerados",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/nlem_log"),
        help="Pasta de saída para tabelas/figuras/relatório",
    )
    parser.add_argument(
        "--max-curves-per-family",
        type=int,
        default=40,
        help="Máximo de curvas por figura de família (evita poluição visual)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolução dos gráficos",
    )
    parser.add_argument(
        "--comparison-log-csi",
        type=float,
        nargs="*",
        default=[24.0, 17.0, 15.0],
        help="Pontos de log10(csi) usados na pasta comparison",
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
        r".*/output/(?:limits|nlem_log)/(?P<model>GM\d+)/B_(?P<b>[^/]+)/(?:"
        r"(?P<topology>default|isotropic|anisotropic)/)?csi_(?P<csi>[^/]+)/eos\.dat$"
    )
    m = rx.match(path.resolve().as_posix())
    if not m:
        return None

    model = m.group("model")
    b_label = m.group("b")
    topology = m.group("topology") or "default"

    b_value = _safe_float(b_label)
    csi = _safe_float(m.group("csi"))
    if b_value is None or csi is None or csi <= 0:
        return None

    return model, b_label, b_value, topology, csi, math.log10(csi)


def load_eos(path: Path) -> Optional[np.ndarray]:
    try:
        arr = np.genfromtxt(path, comments="#")
    except Exception:
        return None

    if arr is None or np.size(arr) == 0:
        return None

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    if arr.shape[1] < 22:
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

    datasets.sort(key=lambda d: (d.model, d.topology, d.b_value, d.log_csi))
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


def valid_mr_mask(mass_col: np.ndarray, radius_col: np.ndarray) -> np.ndarray:
    return (
        np.isfinite(mass_col)
        & np.isfinite(radius_col)
        & (mass_col > 0.0)
        & (radius_col > 0.0)
        & (mass_col <= MAX_VALID_MASS_MSUN)
        & (radius_col <= MAX_VALID_RADIUS_KM)
    )


def summarize_dataset(ds: Dataset) -> SummaryRow:
    arr = ds.data

    mass_col = arr[:, -2]
    radius_col = arr[:, -1]
    mr_mask = valid_mr_mask(mass_col, radius_col)

    n_mr_points = int(np.count_nonzero(mr_mask))
    n_eos_points = int(arr.shape[0])

    if n_mr_points > 0:
        mr_indices = np.nonzero(mr_mask)[0]
        local = np.argmax(mass_col[mr_mask])
        i_max = int(mr_indices[local])

        max_mass = float(arr[i_max, -2])
        r_at_max = float(arr[i_max, -1])
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

    radius_at_14 = _radius_at_canonical_mass(mass_col, radius_col)
    compactness = _compactness_from_mass_radius(max_mass, r_at_max)
    hyperon_onset = _hyperon_threshold_density(arr)

    cs2_min, cs2_max = compute_cs2_bounds(arr[:, COL_EPS], arr[:, COL_P])

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
        radius_at_1p4_msun_km=radius_at_14,
        compactness_at_max_mass=compactness,
        central_nb_n0=central_nb,
        central_eps_mevfm3=central_eps,
        central_p_mevfm3=central_p,
        central_meff=central_meff,
        central_emag_mevfm3=central_emag,
        hyperon_onset_n0=hyperon_onset,
        cs2_min=cs2_min,
        cs2_max=cs2_max,
        eos_path=ds.file_path.as_posix(),
    )


def _stable_branch_indices(mass_col: np.ndarray, radius_col: np.ndarray) -> np.ndarray:
    mr_mask = valid_mr_mask(mass_col, radius_col)
    if not np.any(mr_mask):
        return np.array([], dtype=int)

    mr_indices = np.nonzero(mr_mask)[0]
    local_max = int(np.argmax(mass_col[mr_mask]))
    return mr_indices[: local_max + 1]


def _radius_at_canonical_mass(mass_col: np.ndarray, radius_col: np.ndarray, target_mass: float = TARGET_CANONICAL_MASS_MSUN) -> float:
    stable_idx = _stable_branch_indices(mass_col, radius_col)
    if stable_idx.size < 2:
        return math.nan

    masses = np.asarray(mass_col[stable_idx], dtype=float)
    radii = np.asarray(radius_col[stable_idx], dtype=float)
    mask = np.isfinite(masses) & np.isfinite(radii)
    masses = masses[mask]
    radii = radii[mask]
    if masses.size < 2:
        return math.nan

    order = np.argsort(masses)
    masses = masses[order]
    radii = radii[order]
    masses, unique_idx = np.unique(masses, return_index=True)
    radii = radii[unique_idx]

    if masses.size < 2 or target_mass < masses[0] or target_mass > masses[-1]:
        return math.nan

    return float(np.interp(target_mass, masses, radii))


def _compactness_from_mass_radius(mass_msun: float, radius_km: float) -> float:
    if not np.isfinite(mass_msun) or not np.isfinite(radius_km) or radius_km <= 0.0:
        return math.nan
    return float(G_C2_KM_PER_MSUN * mass_msun / radius_km)


def _hyperon_threshold_density(arr: np.ndarray) -> float:
    if arr.shape[1] <= max(HYPERON_COLUMNS):
        return math.nan

    thresholds = [
        find_onset_threshold(arr[:, COL_NB], arr[:, col])
        for col in HYPERON_COLUMNS
        if arr.shape[1] > col
    ]
    finite = [value for value in thresholds if np.isfinite(value)]
    if not finite:
        return math.nan
    return float(np.min(finite))


def write_summary_csv(rows: Sequence[SummaryRow], out_csv: Path) -> None:
    ensure_dir(out_csv.parent)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "b_label", "b_value", "csi", "log_csi",
            "n_eos_points", "n_mr_points", "max_mass_msun", "radius_at_max_km",
            "radius_at_1p4_msun_km", "compactness_at_max_mass",
            "central_nb_n0", "central_eps_mevfm3", "central_p_mevfm3", "central_meff",
            "central_emag_mevfm3", "hyperon_onset_n0", "cs2_min", "cs2_max", "eos_path",
        ])
        for r in rows:
            w.writerow([
                r.model, r.b_label, r.b_value, r.csi, r.log_csi,
                r.n_eos_points, r.n_mr_points, r.max_mass_msun, r.radius_at_max_km,
                r.radius_at_1p4_msun_km, r.compactness_at_max_mass,
                r.central_nb_n0, r.central_eps_mevfm3, r.central_p_mevfm3,
                r.central_meff, r.central_emag_mevfm3, r.hyperon_onset_n0,
                r.cs2_min, r.cs2_max, r.eos_path,
            ])


def group_key(ds: Dataset) -> Tuple[str, str]:
    return ds.model, ds.b_label


def _subsample_sorted(items: Sequence[Dataset], max_items: int) -> List[Dataset]:
    if len(items) <= max_items:
        return list(items)
    idx = np.linspace(0, len(items) - 1, max_items).round().astype(int)
    idx = np.unique(idx)
    return [items[i] for i in idx]


def _population_percent_of_total(pop: np.ndarray, total: np.ndarray) -> np.ndarray:
    pop = np.asarray(pop, dtype=float)
    total = np.asarray(total, dtype=float)
    percent = np.full_like(pop, np.nan, dtype=float)
    valid = np.isfinite(pop) & np.isfinite(total) & (total > 0.0)
    percent[valid] = 100.0 * pop[valid] / total[valid]
    return percent


def _population_total_for_dataset(ds: Dataset, species: Sequence[Tuple[int, str]]) -> np.ndarray:
    total = np.zeros(ds.data.shape[0], dtype=float)
    for col, _label in species:
        if ds.data.shape[1] <= col:
            continue
        values = np.asarray(ds.data[:, col], dtype=float)
        total += np.where(np.isfinite(values) & (values > 0.0), values, 0.0)
    total[~np.isfinite(total)] = np.nan
    return total


def _threshold_percent(nb_value: float) -> float:
    if not np.isfinite(nb_value) or nb_value <= 0.0:
        return math.nan
    return float(nb_value)


def load_observational_constraints(csv_path: Path = OBSERVATIONAL_CONSTRAINTS_PATH) -> List[ObservationConstraint]:
    if not csv_path.exists():
        return []

    constraints: List[ObservationConstraint] = []
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                constraints.append(
                    ObservationConstraint(
                        kind=row.get("kind", "").strip(),
                        label=row.get("label", "").strip(),
                        source=row.get("source", "").strip(),
                        mass_msun=float(row.get("mass_msun", "nan") or "nan"),
                        mass_err_minus=float(row.get("mass_err_minus", "nan") or "nan"),
                        mass_err_plus=float(row.get("mass_err_plus", "nan") or "nan"),
                        radius_km=float(row.get("radius_km", "nan") or "nan"),
                        radius_err_minus=float(row.get("radius_err_minus", "nan") or "nan"),
                        radius_err_plus=float(row.get("radius_err_plus", "nan") or "nan"),
                        notes=row.get("notes", "").strip(),
                    )
                )
            except ValueError:
                continue
    return constraints


def _observation_style(kind: str) -> Tuple[str, str, str]:
    if kind == "radius_constraint":
        return "tab:green", "D", "R1.4 constraint"
    if kind == "nicer_point":
        return "tab:purple", "o", "NICER"
    return "tab:gray", "s", "Obs."


def overlay_observational_constraints(ax: Axes, constraints: Sequence[ObservationConstraint]) -> List[Line2D]:
    handles: List[Line2D] = []
    for obs in constraints:
        color, marker, _short = _observation_style(obs.kind)
        has_mass = np.isfinite(obs.mass_msun)
        has_radius = np.isfinite(obs.radius_km)
        if not (has_mass and has_radius):
            continue

        xerr = None
        if np.isfinite(obs.mass_err_minus) or np.isfinite(obs.mass_err_plus):
            xerr = np.array([[0.0 if not np.isfinite(obs.mass_err_minus) else obs.mass_err_minus], [0.0 if not np.isfinite(obs.mass_err_plus) else obs.mass_err_plus]])
        yerr = None
        if np.isfinite(obs.radius_err_minus) or np.isfinite(obs.radius_err_plus):
            yerr = np.array([[0.0 if not np.isfinite(obs.radius_err_minus) else obs.radius_err_minus], [0.0 if not np.isfinite(obs.radius_err_plus) else obs.radius_err_plus]])

        ax.errorbar(
            obs.radius_km,
            obs.mass_msun,
            xerr=yerr,
            yerr=xerr,
            fmt=marker,
            ms=7,
            mfc="white",
            mec=color,
            ecolor=color,
            elinewidth=1.2,
            capsize=3,
            linestyle="none",
            zorder=5,
        )
        handles.append(
            Line2D([0], [0], color=color, marker=marker, linestyle="none", markersize=7, markerfacecolor="white", markeredgecolor=color, label=f"{obs.label} ({obs.source})")
        )

    return handles


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


def _effective_h_log(bg: np.ndarray, csi: np.ndarray) -> np.ndarray:
    """
    Campo efetivo para o modelo Log (Soleng):
    H = bg / (1 + bg^2 / (2*csi^2))
    """
    bg = np.asarray(bg, dtype=float)
    csi = np.asarray(csi, dtype=float)
    csi_safe = np.maximum(csi, 1e-300)
    denom = 1.0 + (bg * bg) / (2.0 * csi_safe * csi_safe)
    return bg / denom


def plot_effective_h_vs_csi(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plota H(csi) para o modelo Log usando a fórmula de Soleng,
    agrupando por (model, B).
    """
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault((ds.model, ds.b_label), []).append(ds)

    for (model, b_label), arr in by_combo.items():
        if not arr:
            continue

        arr_sorted = sorted(arr, key=lambda d: d.log_csi)
        csi = np.array([d.csi for d in arr_sorted], dtype=float)
        log_csi = np.array([d.log_csi for d in arr_sorted], dtype=float)
        bg = np.array([d.b_value for d in arr_sorted], dtype=float)

        h_eff = _effective_h_log(bg, csi)
        ratio = np.divide(h_eff, bg, out=np.zeros_like(h_eff), where=np.abs(bg) > 1e-300)

        finite = np.isfinite(log_csi) & np.isfinite(h_eff) & np.isfinite(ratio)
        if np.count_nonzero(finite) < 2:
            continue

        x = log_csi[finite]
        y_h = h_eff[finite]
        y_r = ratio[finite]

        fig, ax1 = plt.subplots(figsize=(8, 6))
        ax1.plot(x, y_h, color=H_COLOR, marker="o", ms=4, lw=1.4, label="H")
        ax1.set_xlabel(LOG_CSI_LABEL)
        ax1.set_ylabel("H efetivo [G]", color=H_COLOR)
        ax1.set_yscale("log")
        ax1.tick_params(axis="y", labelcolor=H_COLOR)
        ax1.grid(alpha=0.25)

        ax2 = ax1.twinx()
        ax2.plot(x, y_r, color=HB_RATIO_COLOR, linestyle="--", lw=1.4, label="H/B")
        ax2.set_ylabel("Razão H/B", color=HB_RATIO_COLOR)
        ax2.tick_params(axis="y", labelcolor=HB_RATIO_COLOR)

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="best", fontsize=9)

        plt.title(f"Campo efetivo H(csi) | {model} | B={b_label} G")
        fig.tight_layout()
        fig.savefig(out_dir / f"h_vs_csi_{model}_B_{b_label}.png", dpi=dpi)
        plt.close(fig)


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


def _infer_mu_to_mev_factor(mu_raw: np.ndarray) -> float:
    """
    Heurística genérica para potenciais químicos: se o valor típico for muito
    pequeno, assume-se normalização por m_n e converte para MeV.
    """
    finite = mu_raw[np.isfinite(mu_raw)]
    if finite.size == 0:
        return 1.0
    max_val = float(np.nanmax(finite))
    return 939.0 if max_val < 5.0 else 1.0


def _linestyle_for_index(index: int) -> str:
    styles = ["-", "--", "-.", ":"]
    return styles[index % len(styles)]


def _linestyle_for_index_distinct(index: int):
    styles = [
        "-",
        "--",
        ":",
        "-.",
        (0, (1, 1)),
        (0, (5, 1)),
        (0, (3, 1, 1, 1)),
        (0, (8, 2, 1, 2)),
    ]
    return styles[index % len(styles)]


def _species_color_map(species: Sequence[Tuple[int, str]]) -> Dict[str, Tuple[float, float, float, float]]:
    return {label: plt.get_cmap("tab10")(idx % 10) for idx, (_, label) in enumerate(species)}


def _build_csi_handles(arr_sorted: Sequence[Dataset], distinct_styles: bool = False) -> List[Line2D]:
    handles: List[Line2D] = []
    style_fn = _linestyle_for_index_distinct if distinct_styles else _linestyle_for_index
    for idx, ds in enumerate(arr_sorted):
        linestyle = style_fn(idx)
        handles.append(
            Line2D([0], [0], color="black", linestyle=linestyle, lw=1.4, label=f"log10(csi)={ds.log_csi:.3f}")
        )
    return handles


def _style_axis(ax: Axes, xlabel: str, ylabel: str, title: str, grid_alpha: float = 0.25, grid_ls: str = "-") -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=grid_alpha, linestyle=grid_ls)


def _fixed_color_lookup(datasets: Sequence[Dataset]) -> Dict[float, Tuple[float, float, float, float]]:
    logs = sorted({d.log_csi for d in datasets})
    cmap = plt.get_cmap("tab10")
    return {log: cmap(i % 10) for i, log in enumerate(logs)}


def _fixed_color_handles(color_lookup: Dict[float, Tuple[float, float, float, float]]) -> List[Line2D]:
    return [
        Line2D([0], [0], color=color, lw=1.6, label=f"log10(csi)={log:.3f}")
        for log, color in sorted(color_lookup.items(), key=lambda kv: kv[0])
    ]


def _plot_single_landau_curve(
    ds: Dataset,
    cmap: mcolors.Colormap,
    cmin: float,
    denom: float,
    fixed_color: Optional[Tuple[float, float, float, float]] = None,
) -> bool:
    if ds.data.shape[1] <= COL_BDD:
        return False

    nb = ds.data[:, COL_NB]
    mu_e_raw = ds.data[:, COL_MUE]
    mu_factor = _infer_mue_to_mev_factor(mu_e_raw)
    mu_e = mu_e_raw * mu_factor
    b_tesla = ds.data[:, COL_BDD]
    nu_max = _compute_electron_landau_nu_max(mu_e, b_tesla)
    if not np.any(np.isfinite(nu_max)):
        return False

    color = fixed_color if fixed_color is not None else cmap((ds.log_csi - cmin) / denom)
    plt.step(
        nb,
        nu_max,
        where="post",
        color=color,
        alpha=0.85,
        lw=1.3,
    )
    return True


def plot_landau_levels(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
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
        arr_sorted = sorted(arr, key=lambda d: d.log_csi)
        plot_set = arr_sorted
        if not plot_set:
            continue

        cvals = np.array([d.log_csi for d in plot_set], dtype=float)
        cmin, cmax = float(np.min(cvals)), float(np.max(cvals))
        denom = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0
        cmap = plt.get_cmap("plasma")
        color_lookup = _fixed_color_lookup(plot_set) if fixed_colors else None

        plt.figure(figsize=(8, 6))
        plotted = False

        for d in plot_set:
            fixed_color = color_lookup.get(d.log_csi) if color_lookup is not None else None
            plotted = _plot_single_landau_curve(d, cmap=cmap, cmin=cmin, denom=denom, fixed_color=fixed_color) or plotted

        if plotted:
            plt.xlabel(r"Densidade bariônica $n_B/n_0$")
            plt.ylabel(r"Nível máximo de Landau eletrônico $\nu_{\max}^e$")
            plt.title(f"Quantização de Landau | {model} | B={b_label} G")
            plt.grid(alpha=0.3)
            if fixed_colors and color_lookup is not None:
                plt.legend(handles=_fixed_color_handles(color_lookup), fontsize=8, title="csi", loc="best")
            elif show_colorbar:
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
                cbar = plt.colorbar(sm, ax=plt.gca())
                cbar.set_label(LOG_CSI_LABEL)
            plt.tight_layout()
            plt.savefig(out_dir / f"landau_levels_{model}_B_{b_label}.png", dpi=dpi)
        plt.close()


def plot_chemical_potential_profiles(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
) -> None:
    """
    Plota os potenciais químicos nucleônico e eletrônico em função da densidade
    bariônica, colorindo cada curva por log10(csi).
    """
    ensure_dir(out_dir)

    by_combo: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        by_combo.setdefault(group_key(ds), []).append(ds)

    for key, arr in by_combo.items():
        _plot_chemical_potential_family(
            key,
            arr,
            out_dir,
            dpi,
            fixed_colors=fixed_colors,
            show_colorbar=show_colorbar,
        )


def _plot_single_chemical_potential_curve(
    ds: Dataset,
    ax_mu_n: Axes,
    ax_mu_e: Axes,
    color: Tuple[float, float, float, float],
) -> bool:
    plotted = False

    nb = ds.data[:, COL_NB]

    if ds.data.shape[1] > COL_MUN:
        mun_raw = ds.data[:, COL_MUN]
        mun = mun_raw * _infer_mu_to_mev_factor(mun_raw)
        mask_mun = np.isfinite(nb) & np.isfinite(mun)
        if np.count_nonzero(mask_mun) > 1:
            ax_mu_n.plot(nb[mask_mun], mun[mask_mun], color=color, alpha=0.85, lw=1.0)
            plotted = True

    if ds.data.shape[1] > COL_MUE:
        mue_raw = ds.data[:, COL_MUE]
        mue = mue_raw * _infer_mue_to_mev_factor(mue_raw)
        mask_mue = np.isfinite(nb) & np.isfinite(mue)
        if np.count_nonzero(mask_mue) > 1:
            ax_mu_e.plot(nb[mask_mue], mue[mask_mue], color=color, alpha=0.85, lw=1.0)
            plotted = True

    return plotted


def _plot_chemical_potential_family(
    key: Tuple[str, str],
    arr: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
) -> None:
    model, b_label = key
    arr_sorted = sorted(arr, key=lambda d: d.log_csi)
    if not arr_sorted:
        return

    cvals = np.array([d.log_csi for d in arr_sorted], dtype=float)
    cmin, cmax = float(np.min(cvals)), float(np.max(cvals))
    denom = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0
    cmap = plt.get_cmap("viridis")
    color_lookup = _fixed_color_lookup(arr_sorted) if fixed_colors else None

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    plotted_any = False

    for ds in arr_sorted:
        color = color_lookup[ds.log_csi] if color_lookup is not None else cmap((ds.log_csi - cmin) / denom)
        plotted_any = _plot_single_chemical_potential_curve(ds, axs[0], axs[1], color) or plotted_any

    if not plotted_any:
        plt.close(fig)
        return

    axs[0].set_xlabel(NB_DENSITY_LABEL)
    axs[0].set_ylabel(r"$\mu_n$ [MeV]")
    axs[0].set_title(r"Potencial químico do nêutron")
    axs[0].grid(alpha=0.25)

    axs[1].set_xlabel(NB_DENSITY_LABEL)
    axs[1].set_ylabel(r"$\mu_e$ [MeV]")
    axs[1].set_title(r"Potencial químico do elétron")
    axs[1].grid(alpha=0.25)

    if fixed_colors and color_lookup is not None:
        axs[0].legend(handles=_fixed_color_handles(color_lookup), fontsize=8, title="csi", loc="best")
    elif show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axs, pad=0.02)
        cbar.set_label(LOG_CSI_LABEL)

    fig.suptitle(f"Potenciais químicos | {model} | B={b_label} G", y=1.03)
    fig.tight_layout()
    fig.savefig(out_dir / f"chemical_potentials_{model}_B_{b_label}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _log_csi_scale(datasets: Sequence[Dataset]) -> Tuple[float, float, float]:
    cvals = np.array([d.log_csi for d in datasets], dtype=float)
    cmin = float(np.min(cvals))
    cmax = float(np.max(cvals))
    span = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0
    return cmin, cmax, span


def _log_csi_color(log_csi: float, cmap: mcolors.Colormap, cmin: float, span: float):
    return cmap((log_csi - cmin) / span)


def _select_comparison_datasets(
    datasets: Sequence[Dataset],
    targets: Sequence[float],
) -> List[Tuple[Dataset, float]]:
    selected: List[Tuple[Dataset, float]] = []
    by_combo = _group_datasets_by_combo(datasets)

    for arr in by_combo.values():
        arr_sorted = sorted(arr, key=lambda d: d.log_csi)
        for target in targets:
            if not arr_sorted:
                continue
            chosen = min(arr_sorted, key=lambda d, target=target: abs(d.log_csi - target))
            pair = (chosen, target)
            if pair not in selected:
                selected.append(pair)

    selected.sort(key=lambda item: (item[0].model, item[0].b_label, item[0].log_csi, item[1]))
    return selected


def _group_datasets_by_combo(datasets: Sequence[Dataset]) -> Dict[Tuple[str, str], List[Dataset]]:
    grouped: Dict[Tuple[str, str], List[Dataset]] = {}
    for ds in datasets:
        grouped.setdefault(group_key(ds), []).append(ds)
    return grouped


def _sorted_plot_set(arr: Sequence[Dataset], max_curves: int) -> List[Dataset]:
    arr_sorted = sorted(arr, key=lambda d: d.log_csi)
    return _subsample_sorted(arr_sorted, max_curves)


def _plot_eos_family_panel(
    plot_set: Sequence[Dataset],
    model: str,
    b_label: str,
    cmin: float,
    cmax: float,
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
) -> None:
    cmap_eos = plt.get_cmap("viridis")
    color_lookup = _fixed_color_lookup(plot_set) if fixed_colors else None

    plt.figure(figsize=(8, 6))
    span = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0
    for d in plot_set:
        if color_lookup is not None:
            color = color_lookup[d.log_csi]
        else:
            color = _log_csi_color(d.log_csi, cmap_eos, cmin, span)
        plt.plot(d.data[:, COL_EPS], d.data[:, COL_P], color=color, alpha=0.8, lw=1.0)

    plt.xlabel(r"$\epsilon$ [MeV/fm$^3$]")
    plt.ylabel(r"$P$ [MeV/fm$^3$]")
    plt.title(f"EoS family | {model} | B={b_label} G")
    plt.grid(alpha=0.25)

    if fixed_colors and color_lookup is not None:
        plt.legend(handles=_fixed_color_handles(color_lookup), fontsize=8, title="csi", loc="best")
    elif show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap_eos, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(LOG_CSI_LABEL)
    plt.tight_layout()
    plt.savefig(out_dir / f"eos_family_{model}_B_{b_label}.png", dpi=dpi)
    plt.close()


def _plot_mr_family_panel(
    plot_set: Sequence[Dataset],
    model: str,
    b_label: str,
    cmin: float,
    cmax: float,
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
    truncate_at_max_mass: bool = False,
) -> None:
    cmap_mr = plt.get_cmap("plasma")
    color_lookup = _fixed_color_lookup(plot_set) if fixed_colors else None
    span = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0

    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    for d in plot_set:
        m = d.data[:, -2]
        r = d.data[:, -1]
        p_central = d.data[:, COL_P]

        mask = valid_mr_mask(m, r)
        if np.count_nonzero(mask) < 3:
            continue

        m_valid = m[mask]
        r_valid = r[mask]
        p_valid = p_central[mask]
        sort_idx = np.argsort(p_valid)
        m_valid = m_valid[sort_idx]
        r_valid = r_valid[sort_idx]
        if truncate_at_max_mass and m_valid.size > 0:
            max_idx = int(np.nanargmax(m_valid))
            m_valid = m_valid[: max_idx + 1]
            r_valid = r_valid[: max_idx + 1]
        if color_lookup is not None:
            color = color_lookup[d.log_csi]
        else:
            color = _log_csi_color(d.log_csi, cmap_mr, cmin, span)
        ax.plot(r_valid, m_valid, color=color, alpha=0.8, lw=1.0)

    obs_handles = overlay_observational_constraints(ax, load_observational_constraints())

    ax.set_xlabel("Radius [km]")
    ax.set_ylabel(r"Mass [$M_\odot$]")
    ax.set_title(f"M-R family | {model} | B={b_label} G")
    ax.set_xlim(MR_PLOT_MIN_RADIUS_KM, MR_PLOT_MAX_RADIUS_KM)
    ax.set_ylim(MR_PLOT_MIN_MASS_MSUN, MR_PLOT_MAX_MASS_MSUN)
    ax.grid(alpha=0.25)

    if fixed_colors and color_lookup is not None:
        model_legend = ax.legend(handles=_fixed_color_handles(color_lookup), fontsize=8, title="csi", loc=LEGEND_LOC_UPPER_LEFT)
        ax.add_artist(model_legend)
        if obs_handles:
            ax.legend(handles=obs_handles, fontsize=7, loc=LEGEND_LOC_UPPER_LEFT, title="Observações")
    elif show_colorbar:
        sm_mr = plt.cm.ScalarMappable(cmap=cmap_mr, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
        cbar_mr = plt.colorbar(sm_mr, ax=ax)
        cbar_mr.set_label(LOG_CSI_LABEL)
        if obs_handles:
            ax.legend(handles=obs_handles, fontsize=7, loc=LEGEND_LOC_UPPER_LEFT, title="Observações")
    plt.tight_layout()
    plt.savefig(out_dir / f"mr_family_{model}_B_{b_label}.png", dpi=dpi)
    plt.close()


def _plot_cs2_family_panel(
    plot_set: Sequence[Dataset],
    model: str,
    b_label: str,
    cmin: float,
    cmax: float,
    out_dir: Path,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
) -> None:
    cmap_cs2 = plt.get_cmap("cool")
    color_lookup = _fixed_color_lookup(plot_set) if fixed_colors else None
    span = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0

    plt.figure(figsize=(8, 6))
    for d in plot_set:
        eps = d.data[:, COL_EPS]
        p = d.data[:, COL_P]
        mask = np.isfinite(eps) & np.isfinite(p)
        eps_val = eps[mask]
        p_val = p[mask]

        if eps_val.size <= 2:
            continue

        idx = np.argsort(eps_val)
        eps_sorted = eps_val[idx]
        p_sorted = p_val[idx]
        de = np.diff(eps_sorted)
        dp = np.diff(p_sorted)

        good = de > 1e-14
        if not np.any(good):
            continue

        cs2 = np.zeros_like(de)
        cs2[good] = dp[good] / de[good]
        if color_lookup is not None:
            color = color_lookup[d.log_csi]
        else:
            color = _log_csi_color(d.log_csi, cmap_cs2, cmin, span)
        plt.plot(eps_sorted[:-1][good], cs2[good], color=color, alpha=0.8, lw=1.0)

    plt.axhline(1.0, color="red", linestyle="--", alpha=0.7, label="Causalidade ($c_s^2 = 1$)")
    plt.axhline(0.0, color="black", linestyle="--", alpha=0.7, label="Estabilidade ($c_s^2 = 0$)")
    plt.axhline(1 / 3, color="gray", linestyle=":", alpha=0.7, label="Limite Conforme ($c_s^2 = 1/3$)")

    plt.xlabel(r"$\epsilon$ [MeV/fm$^3$]")
    plt.ylabel(r"$c_s^2$ (Unidades de $c=1$)")
    plt.title(f"Speed of Sound | {model} | B={b_label} G")
    plt.grid(alpha=0.25)
    plt.legend(loc=LEGEND_LOC_UPPER_RIGHT, fontsize=9)

    if fixed_colors and color_lookup is not None:
        plt.legend(handles=_fixed_color_handles(color_lookup), loc=LEGEND_LOC_UPPER_RIGHT, fontsize=8, title="csi")
    elif show_colorbar:
        sm_cs2 = plt.cm.ScalarMappable(cmap=cmap_cs2, norm=mcolors.Normalize(vmin=cmin, vmax=cmax))
        cbar_cs2 = plt.colorbar(sm_cs2, ax=plt.gca())
        cbar_cs2.set_label(LOG_CSI_LABEL)

    plt.tight_layout()
    plt.savefig(out_dir / f"cs2_family_{model}_B_{b_label}.png", dpi=dpi)
    plt.close()


def write_comparison_selection_csv(
    selected: Sequence[Tuple[Dataset, float]],
    out_file: Path,
) -> None:
    ensure_dir(out_file.parent)

    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "model",
            "b_label",
            "b_value",
            "log_b_value",
            "selected_log_csi",
            "target_log_csi",
            "max_mass_msun",
            "radius_at_max_km",
            "radius_at_1p4_msun_km",
            "compactness_at_max_mass",
            "central_nb_n0",
            "central_eps_mevfm3",
            "central_p_mevfm3",
            "hyperon_onset_n0",
            "eos_path",
        ])

        for ds, target in selected:
            arr = ds.data
            mass_col = arr[:, -2]
            radius_col = arr[:, -1]
            stable_idx = _stable_branch_indices(mass_col, radius_col)

            if stable_idx.size > 0:
                local_max = int(np.argmax(mass_col[stable_idx]))
                i_max = int(stable_idx[local_max])
                max_mass = float(mass_col[i_max])
                radius_at_max = float(radius_col[i_max])
                central_nb = float(arr[i_max, COL_NB])
                central_eps = float(arr[i_max, COL_EPS])
                central_p = float(arr[i_max, COL_P])
            else:
                max_mass = math.nan
                radius_at_max = math.nan
                central_nb = math.nan
                central_eps = math.nan
                central_p = math.nan

            radius_at_14 = _radius_at_canonical_mass(mass_col, radius_col)
            compactness = _compactness_from_mass_radius(max_mass, radius_at_max)
            n_hyp = _hyperon_threshold_density(arr)

            writer.writerow([
                ds.model,
                ds.b_label,
                ds.b_value,
                math.log10(ds.b_value) if ds.b_value > 0 else math.nan,
                ds.log_csi,
                target,
                max_mass,
                radius_at_max,
                radius_at_14,
                compactness,
                central_nb,
                central_eps,
                central_p,
                n_hyp,
                ds.file_path.as_posix(),
            ])


def plot_comparison_h_vs_csi(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Versão de comparação para H(csi): destaca os pontos selecionados e usa
    Maxwell como referência explícita com H/B = 1.
    """
    ensure_dir(out_dir)

    by_combo = _group_datasets_by_combo(datasets)

    for (model, b_label), arr in by_combo.items():
        if not arr:
            continue

        arr_sorted = sorted(arr, key=lambda d: d.log_csi)
        csi = np.array([d.csi for d in arr_sorted], dtype=float)
        log_csi = np.array([d.log_csi for d in arr_sorted], dtype=float)
        bg = np.array([d.b_value for d in arr_sorted], dtype=float)

        h_eff = _effective_h_log(bg, csi)
        ratio = np.divide(h_eff, bg, out=np.zeros_like(h_eff), where=np.abs(bg) > 1e-300)

        finite = np.isfinite(log_csi) & np.isfinite(h_eff) & np.isfinite(ratio)
        if np.count_nonzero(finite) < 2:
            continue

        x = log_csi[finite]
        y_h = h_eff[finite]
        y_r = ratio[finite]

        fig, ax1 = plt.subplots(figsize=(8, 6))
        ax1.plot(x, y_h, color=H_COLOR, marker="o", ms=4, lw=1.4, label="NLEM: H")
        ax1.set_xlabel(LOG_CSI_LABEL)
        ax1.set_ylabel("H efetivo [G]", color=H_COLOR)
        ax1.set_yscale("log")
        ax1.tick_params(axis="y", labelcolor=H_COLOR)
        ax1.grid(alpha=0.25)

        ax2 = ax1.twinx()
        ax2.plot(x, y_r, color=HB_RATIO_COLOR, linestyle="--", lw=1.4, label="NLEM: H/B")
        ax2.axhline(1.0, color="black", linestyle=":", alpha=0.75, label="Maxwell: H/B = 1")
        ax2.set_ylabel("Razão H/B", color=HB_RATIO_COLOR)
        ax2.tick_params(axis="y", labelcolor=HB_RATIO_COLOR)

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="best", fontsize=9)

        plt.title(f"Comparison H(csi) | {model} | B={b_label} G")
        fig.tight_layout()
        fig.savefig(out_dir / f"comparison_h_vs_csi_{model}_B_{b_label}.png", dpi=dpi)
        plt.close(fig)


def plot_comparison_folder(
    datasets: Sequence[Dataset],
    out_dir: Path,
    comparison_targets: Sequence[float],
    max_curves_per_family: int,
    dpi: int,
) -> None:
    ensure_dir(out_dir)
    selected = _select_comparison_datasets(datasets, comparison_targets)
    if not selected:
        return

    write_comparison_selection_csv(selected, out_dir / "comparison_selection.csv")

    selected_datasets = [ds for ds, _target in selected]

    plot_family_eos_mr(
        selected_datasets,
        out_dir / "families",
        max_curves=max(1, max_curves_per_family),
        dpi=dpi,
        fixed_colors=True,
        show_colorbar=False,
        truncate_at_max_mass=True,
    )
    plot_population_snapshots(selected_datasets, out_dir / "populations", dpi=dpi, distinct_styles=True)
    plot_population_thresholds(selected_datasets, out_dir / "thresholds", dpi=dpi)
    plot_landau_levels(selected_datasets, out_dir / "landau_levels", dpi=dpi, fixed_colors=True, show_colorbar=False)
    plot_chemical_potential_profiles(
        selected_datasets,
        out_dir / "chemical_potentials",
        dpi=dpi,
        fixed_colors=True,
        show_colorbar=False,
    )


def plot_family_eos_mr(
    datasets: Sequence[Dataset],
    out_dir: Path,
    max_curves: int,
    dpi: int,
    fixed_colors: bool = False,
    show_colorbar: bool = True,
    truncate_at_max_mass: bool = False,
) -> None:
    ensure_dir(out_dir)

    by_combo = _group_datasets_by_combo(datasets)

    for key, arr in by_combo.items():
        model, b_label = key
        plot_set = _sorted_plot_set(arr, max_curves)
        if not plot_set:
            continue

        cmin, cmax, _span = _log_csi_scale(plot_set)
        _plot_eos_family_panel(plot_set, model, b_label, cmin, cmax, out_dir, dpi, fixed_colors=fixed_colors, show_colorbar=show_colorbar)
        _plot_mr_family_panel(
            plot_set,
            model,
            b_label,
            cmin,
            cmax,
            out_dir,
            dpi,
            fixed_colors=fixed_colors,
            show_colorbar=show_colorbar,
            truncate_at_max_mass=truncate_at_max_mass,
        )
        _plot_cs2_family_panel(plot_set, model, b_label, cmin, cmax, out_dir, dpi, fixed_colors=fixed_colors, show_colorbar=show_colorbar)


def plot_population_snapshots(
    datasets: Sequence[Dataset],
    out_dir: Path,
    dpi: int,
    distinct_styles: bool = False,
    show_colorbar: bool = True,
) -> None:
    ensure_dir(out_dir)
    by_combo = _group_datasets_by_combo(datasets)

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
        model, b_label = key
        if distinct_styles:
            _plot_population_overlay_family(
                arr,
                species,
                model,
                b_label,
                out_dir,
                dpi,
                distinct_styles=True,
                show_colorbar=show_colorbar,
            )
        else:
            _plot_population_snapshot_family(arr, species, model, b_label, out_dir, dpi)


def _plot_population_snapshot_family(
    arr: Sequence[Dataset],
    species: Sequence[Tuple[int, str]],
    model: str,
    b_label: str,
    out_dir: Path,
    dpi: int,
) -> None:
    arr_sorted = sorted(arr, key=lambda d: d.log_csi)
    if not arr_sorted:
        return

    picks = np.unique(np.linspace(0, len(arr_sorted) - 1, min(3, len(arr_sorted))).round().astype(int)).tolist()

    ncols = len(picks)
    fig, axs = plt.subplots(1, ncols, figsize=(6 * ncols, 5), squeeze=False)

    for j, idx in enumerate(picks):
        ds = arr_sorted[idx]
        ax = axs[0, j]
        x = ds.data[:, COL_NB]
        total = _population_total_for_dataset(ds, species)

        for col, label in species:
            if ds.data.shape[1] > col:
                y = _population_percent_of_total(ds.data[:, col], total)
                ax.plot(x, y, lw=1.1, label=label)

        ax.set_title(f"log10(csi)={ds.log_csi:.3f}")
        ax.set_xlabel(r"$n_B/n_0$")
        ax.set_ylabel(POPULATION_DENSITY_LABEL)
        ax.grid(alpha=0.25)

    handles, labels = axs[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=5)

    fig.suptitle(f"Population snapshots | {model} | B={b_label} G", y=1.03)
    fig.tight_layout()
    fig.savefig(out_dir / f"population_snapshots_{model}_B_{b_label}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_population_overlay_family(
    arr: Sequence[Dataset],
    species: Sequence[Tuple[int, str]],
    model: str,
    b_label: str,
    out_dir: Path,
    dpi: int,
    distinct_styles: bool = False,
    show_colorbar: bool = True,
) -> None:
    arr_sorted = sorted(arr, key=lambda d: d.log_csi)
    if not arr_sorted:
        return

    fig, ax = plt.subplots(figsize=(9, 6))
    plotted_species = set()

    if distinct_styles:
        species_colors = _species_color_map(species)
        csi_handles = _build_csi_handles(arr_sorted, distinct_styles=True)
        style_fn = _linestyle_for_index_distinct

        for idx, ds in enumerate(arr_sorted):
            linestyle = style_fn(idx)
            x = ds.data[:, COL_NB]
            total = _population_total_for_dataset(ds, species)

            for col, label in species:
                if ds.data.shape[1] <= col:
                    continue
                y = _population_percent_of_total(ds.data[:, col], total)
                color = species_colors[label]
                plotted_label = label if label not in plotted_species else None
                ax.plot(x, y, color=color, linestyle=linestyle, lw=1.1, alpha=0.9, label=plotted_label)
                plotted_species.add(label)

        _style_axis(ax, NB_DENSITY_LABEL, POPULATION_DENSITY_LABEL, f"Population overlay | {model} | B={b_label} G")

        species_handles = [Line2D([0], [0], color=species_colors[label], lw=1.6, label=label) for _, label in species]
        legend_species = ax.legend(handles=species_handles, loc=LEGEND_LOC_UPPER_LEFT, fontsize=8, ncol=2, title="Espécies")
        ax.add_artist(legend_species)
        ax.legend(handles=csi_handles, loc="upper right", fontsize=8, title="csi", ncol=1)
    else:
        cvals = np.array([d.log_csi for d in arr_sorted], dtype=float)
        cvals = cvals[np.isfinite(cvals)]
        if cvals.size == 0:
            plt.close(fig)
            return
        cmin = float(np.min(cvals))
        cmax = float(np.max(cvals))
        norm = mcolors.Normalize(vmin=cmin, vmax=cmax)
        cmap = plt.get_cmap("plasma")
        species_styles = {label: _linestyle_for_index(idx) for idx, (_, label) in enumerate(species)}

        for ds in arr_sorted:
            color = cmap(norm(ds.log_csi))
            x = ds.data[:, COL_NB]
            total = _population_total_for_dataset(ds, species)

            for col, label in species:
                if ds.data.shape[1] <= col:
                    continue
                y = _population_percent_of_total(ds.data[:, col], total)
                linestyle = species_styles[label]
                plotted_label = label if label not in plotted_species else None
                ax.plot(x, y, color=color, linestyle=linestyle, lw=1.1, alpha=0.9, label=plotted_label)
                plotted_species.add(label)

        _style_axis(ax, NB_DENSITY_LABEL, POPULATION_DENSITY_LABEL, f"Population overlay | {model} | B={b_label} G")

        species_handles = [Line2D([0], [0], color="black", linestyle=species_styles[label], lw=1.6, label=label) for _, label in species]
        legend_species = ax.legend(handles=species_handles, loc=LEGEND_LOC_UPPER_LEFT, fontsize=8, ncol=2, title="Espécies")
        ax.add_artist(legend_species)
        if show_colorbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.02)
            cbar.set_label(LOG_CSI_LABEL)

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
    Plota a densidade bariônica de surgimento (onset) das partículas em função de log10(csi).
    Foco nas partículas que 'nascem' em altas densidades (Múons e Hiperons).
    """
    ensure_dir(out_dir)
    by_combo = _group_datasets_by_combo(datasets)

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
        model, b_label = key
        _plot_population_threshold_family(arr, species_to_track, model, b_label, out_dir, dpi)


def _plot_population_threshold_family(
    arr: Sequence[Dataset],
    species_to_track: Sequence[Tuple[int, str]],
    model: str,
    b_label: str,
    out_dir: Path,
    dpi: int,
) -> None:
    arr_sorted = sorted(arr, key=lambda d: d.log_csi)
    if not arr_sorted:
        return

    x_vals = np.array([d.log_csi for d in arr_sorted], dtype=float)
    fig, ax = plt.subplots(figsize=(8, 6))
    plotted_any = False
    raw_series: List[Tuple[str, np.ndarray]] = []

    for col, label in species_to_track:
        y_arr = np.array([
            find_onset_threshold(ds.data[:, COL_NB], ds.data[:, col]) if ds.data.shape[1] > col else math.nan
            for ds in arr_sorted
        ], dtype=float)
        raw_series.append((label, y_arr))

    stacked = np.vstack([series for _label, series in raw_series]) if raw_series else np.empty((0, len(arr_sorted)))
    totals = np.nansum(np.where(np.isfinite(stacked), stacked, 0.0), axis=0) if raw_series else np.array([])

    for label, y_arr in raw_series:
        y_pct = np.full_like(y_arr, np.nan, dtype=float)
        mask = np.isfinite(y_arr) & np.isfinite(totals) & (totals > 0.0)
        y_pct[mask] = 100.0 * y_arr[mask] / totals[mask]
        mask = np.isfinite(y_pct)
        if np.count_nonzero(mask) > 0:
            ax.scatter(x_vals[mask], y_pct[mask], s=22, label=label)
            plotted_any = True

    if plotted_any:
        _style_axis(
            ax,
            LOG_CSI_LABEL,
            r"Onset Density [% of total]",
            f"Particle Onset Thresholds | {model} | B={b_label} G",
            grid_alpha=0.25,
            grid_ls="--",
        )
        ax.legend(fontsize=8, loc='best')
        fig.tight_layout()
        fig.savefig(out_dir / f"onset_thresholds_{model}_B_{b_label}.png", dpi=dpi)
    plt.close(fig)


def _plot_metric_for_subset(
    rows: Sequence[SummaryRow],
    model: str,
    metric_name: str,
    y_label: str,
    out_path: Path,
    dpi: int,
) -> None:
    subset = [r for r in rows if r.model == model]
    if not subset:
        return

    bvals = sorted({r.b_value for r in subset})
    plt.figure(figsize=(8, 6))

    for b in bvals:
        part = sorted([r for r in subset if r.b_value == b], key=lambda x: x.log_csi)
        x = np.array([r.log_csi for r in part], dtype=float)
        y = np.array([getattr(r, metric_name) for r in part], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        
        if np.count_nonzero(mask) < 1:
            continue
        plt.plot(x[mask], y[mask], marker="o", ms=3, lw=1.2, label=f"B={part[0].b_label} G")

    plt.xlabel(LOG_CSI_LABEL)
    plt.ylabel(y_label)
    plt.title(f"{y_label} vs log10(csi) | {model}")
    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def plot_trends(rows: Sequence[SummaryRow], out_dir: Path, dpi: int) -> None:
    ensure_dir(out_dir)

    models = sorted({r.model for r in rows})
    
    metrics = [
        ("max_mass_msun", r"$M_{\max}$ [$M_\odot$]", "max_mass_vs_logcsi"),
        ("radius_at_max_km", r"$R(M_{\max})$ [km]", "radius_at_max_vs_logcsi"),
        ("radius_at_1p4_msun_km", r"$R_{1.4}$ [km]", "radius_at_1p4_vs_logcsi"),
        ("compactness_at_max_mass", r"$C(M_{\max})$", "compactness_vs_logcsi"),
        ("central_nb_n0", r"$n_c/n_0$", "central_density_vs_logcsi"),
        ("hyperon_onset_n0", r"$n_{hyp}/n_0$", "hyperon_threshold_vs_logcsi"),
        ("central_emag_mevfm3", r"$\epsilon_{mag,c}$ [MeV/fm$^3$]", "central_emag_vs_logcsi"),
        ("cs2_max", r"$\max(c_s^2)$", "cs2max_vs_logcsi"),
        ("cs2_min", r"$\min(c_s^2)$", "cs2min_vs_logcsi"),
    ]

    for metric_name, y_label, prefix in metrics:
        for model in models:
            _plot_metric_for_subset(
                rows=rows,
                model=model,
                metric_name=metric_name,
                y_label=y_label,
                out_path=out_dir / f"{prefix}_{model}.png",
                dpi=dpi,
            )


def write_scientific_report(rows: Sequence[SummaryRow], out_file: Path) -> None:
    ensure_dir(out_file.parent)
    models = sorted({r.model for r in rows})
    lines: List[str] = []
    _append_report_scope(lines, rows, models)
    _append_report_diagnostics(lines, rows)
    _append_report_trends(lines, rows)
    _append_report_interpretation(lines)
    out_file.write_text("\n".join(lines), encoding="utf-8")


def _append_report_scope(lines: List[str], rows: Sequence[SummaryRow], models: Sequence[str]) -> None:
    lines.append("# NLEM neutron-star analysis report")
    lines.append("")
    lines.append("## Scope")
    lines.append("This report quantifies how $\\log_{10}(\\xi)$ modifies stellar structure using generated EoS/M-R data.")
    lines.append("")
    lines.append(f"- Total datasets analyzed: **{len(rows)}**")
    lines.append(f"- Models: **{', '.join(models)}**")
    lines.append("- A pasta `comparison` reúne cortes selecionados de log10(csi) para comparação direta entre curvas de referência e curvas NLEM.")
    lines.append("")


def _append_report_diagnostics(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    causality_violations = sum(1 for r in rows if r.cs2_max > 1.0)
    stability_violations = sum(1 for r in rows if r.cs2_min < 0.0)

    lines.append("## Causality and Stability Diagnostics")
    lines.append("Analyzing the speed of sound $c_s^2 = dP/d\\epsilon$ to verify physical validity constraints.")
    lines.append("")
    lines.append(f"- **Causality Violations ($c_s^2 > 1$):** {causality_violations} datasets")
    lines.append(f"- **Stability Violations ($c_s^2 < 0$):** {stability_violations} datasets")
    lines.append("")
    if causality_violations > 0 or stability_violations > 0:
        lines.append("**Warning:** Some datasets exhibit non-physical behavior. High $\\xi$ combined with extreme magnetic fields may lead to $c_s^2 > 1.0$ (superluminal sound speed) or $c_s^2 < 0$ (mechanical instability). Inspect the `cs2_max` and `cs2_min` trend graphs to determine the $\\xi$ cutoff.")
    else:
        lines.append("All analyzed datasets respect the physical limits of causality and thermodynamic stability.")
    lines.append("")


def _report_slope_for(sub: Sequence[SummaryRow], metric: str) -> float:
    x = np.array([r.log_csi for r in sub], dtype=float)
    y = np.array([getattr(r, metric) for r in sub], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return math.nan
    a, _b = np.polyfit(x[mask], y[mask], deg=1, full=False)
    return float(a)


def _append_report_trends(lines: List[str], rows: Sequence[SummaryRow]) -> None:
    grouped = _group_rows_by_combo(rows)

    lines.append("## Trend diagnostics by $(model, B)$")
    lines.append("")
    lines.append("Linear slopes (first-order sensitivity):")
    lines.append("- $dM_{max}/d\\log_{10}(\\xi)$")
    lines.append("- $dR(M_{max})/d\\log_{10}(\\xi)$")
    lines.append("")

    for (model, b_label), sub in sorted(
        grouped.items(),
        key=lambda kv: (
            kv[0][0],
            _safe_float(kv[0][1]) if _safe_float(kv[0][1]) is not None else math.inf,
        ),
    ):
        sub_sorted = sorted(sub, key=lambda x: x.log_csi)
        s_mass = _report_slope_for(sub_sorted, "max_mass_msun")
        s_rad = _report_slope_for(sub_sorted, "radius_at_max_km")

        mvals = np.array([r.max_mass_msun for r in sub_sorted], dtype=float)
        rvals = np.array([r.radius_at_max_km for r in sub_sorted], dtype=float)
        cmask = np.isfinite(mvals) & np.isfinite(rvals)

        m_min = float(np.nanmin(mvals)) if np.any(np.isfinite(mvals)) else math.nan
        m_max = float(np.nanmax(mvals)) if np.any(np.isfinite(mvals)) else math.nan
        r_min = float(np.nanmin(rvals)) if np.any(np.isfinite(rvals)) else math.nan
        r_max = float(np.nanmax(rvals)) if np.any(np.isfinite(rvals)) else math.nan

        lines.append(f"### {model} | B={b_label} G")
        lines.append(f"- Samples: {len(sub_sorted)}")
        lines.append(f"- $M_{{max}}$ range: {m_min:.4f} to {m_max:.4f} $M_\\odot$")
        lines.append(f"- $R(M_{{max}})$ range: {r_min:.4f} to {r_max:.4f} km")
        lines.append(f"- Slope $dM_{{max}}/d\\log_{{10}}(\\xi)$: {s_mass:.6f}")
        lines.append(f"- Slope $dR(M_{{max}})/d\\log_{{10}}(\\xi)$: {s_rad:.6f}")
        if np.any(cmask):
            lines.append(f"- Median pair $(M_{{max}}, R)$: ({float(np.nanmedian(mvals)):.4f}, {float(np.nanmedian(rvals)):.4f})")
        lines.append("")


def _group_rows_by_combo(rows: Sequence[SummaryRow]) -> Dict[Tuple[str, str], List[SummaryRow]]:
    grouped: Dict[Tuple[str, str], List[SummaryRow]] = {}
    for row in rows:
        grouped.setdefault((row.model, row.b_label), []).append(row)
    return grouped


def _append_report_interpretation(lines: List[str]) -> None:
    lines.append("## Interpretation guide")
    lines.append("- Positive $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ tends to stiffen the effective sequence in your setup.")
    lines.append("- Negative $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ softens the sequence.")
    lines.append("- Use population thresholds to map how $\\xi$ delays or anticipates hyperon onset.")
    lines.append("- Inspect the chemical-potential figures to compare how $\\mu_n$ and $\\mu_e$ evolve with density for each $(model, B)$ family.")


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

    pops_dir = args.output_root / "populations"
    plot_population_snapshots(datasets, pops_dir, dpi=args.dpi)
    print(f"[OK] Snapshots de populações: {pops_dir}")

    # --- NOVO: Geração dos gráficos de limiar de surgimento (Thresholds) ---
    thresholds_dir = args.output_root / "thresholds"
    plot_population_thresholds(datasets, thresholds_dir, dpi=args.dpi)
    print(f"[OK] Gráficos de Thresholds de Partículas: {thresholds_dir}")

    mu_dir = args.output_root / "chemical_potentials"
    plot_chemical_potential_profiles(datasets, mu_dir, dpi=args.dpi)
    print(f"[OK] Potenciais químicos vs densidade: {mu_dir}")

    comparison_dir = args.output_root / "comparison"
    plot_comparison_folder(
        datasets,
        comparison_dir,
        comparison_targets=args.comparison_log_csi,
        max_curves_per_family=max(1, min(2, args.max_curves_per_family)),
        dpi=args.dpi,
    )
    print(f"[OK] Comparação direta Maxwell/NLEM (pontos selecionados): {comparison_dir}")

    landau_dir = args.output_root / "landau_levels"
    plot_landau_levels(datasets, landau_dir, dpi=args.dpi)
    print(f"[OK] Níveis de Quantização de Landau: {landau_dir}")

    h_dir = args.output_root / "h_vs_csi"
    plot_effective_h_vs_csi(datasets, h_dir, dpi=args.dpi)
    print(f"[OK] Campo efetivo H(csi): {h_dir}")

    report_md = args.output_root / "scientific_report.md"
    write_scientific_report(summary_rows, report_md)
    print(f"[OK] Relatório científico: {report_md}")

    print("[DONE] Análise completa finalizada.")


if __name__ == "__main__":
    main()