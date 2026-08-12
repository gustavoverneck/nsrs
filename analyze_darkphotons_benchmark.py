# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "jinja2>=3.1,<4",
#   "matplotlib>=3.8,<4",
#   "numpy>=1.26,<3",
#   "pandas>=2.1,<4",
# ]
# ///

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import matplotlib
    import numpy as np
    import pandas as pd

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "Dependências ausentes. Execute o script com:\n"
        "  uv run --script analyze_darkphotons_benchmark.py\n"
        f"\nErro original: {exc}"
    )


# =============================================================================
# Configuração científica e gráfica
# =============================================================================

MODELS = ("GM1", "GM3")
SCENARIOS = ("H0", "S1", "S2", "S3")

SCENARIO_LABELS = {
    "H0": "H0 — hadrônica",
    "S1": "S1 — mediador pesado",
    "S2": "S2 — DM pesada",
    "S3": "S3 — mediador leve",
}

# Não especificamos cores manualmente: o matplotlib usa o ciclo padrão.
LINESTYLES = {
    "H0": "-",
    "S1": "--",
    "S2": "-.",
    "S3": ":",
}

N0_FM3 = 0.153

# Mapeamento padrão da saída darkphotons.
DEFAULT_COLUMNS = {
    "nb_over_n0": 0,
    "energy": 1,
    "pressure": 2,
    "n_e": 3,
    "n_mu": 4,
    "n_n": 5,
    "n_p": 6,
    "n_lambda": 7,
    "n_sigma_minus": 8,
    "n_sigma_zero": 9,
    "n_sigma_plus": 10,
    "n_xi_minus": 11,
    "n_xi_zero": 12,
    "sigma": 13,
    "omega": 14,
    "rho": 15,
    "mstar_over_mn": 16,
    "mu_n": 17,
    "mu_e": 18,
    "energy_b": 19,
    "thermo_diag": 20,
    "n_chi": 21,
    "y_chi": 22,
    "m_chi": 23,
    "m_x": 24,
    "epsilon": 25,
    "g_d": 26,
    "x0": 27,
    "kf_chi": 28,
    "mu_chi": 29,
    "energy_chi_kin": 30,
    "pressure_chi_kin": 31,
    "energy_x": 32,
    "pressure_x": 33,
    "mass": 34,
    "radius": 35,
    "baryonic_mass": 36,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analisa os quatro benchmarks H0/S1/S2/S3 para GM1 e GM3."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output/darkphotons_benchmarks"),
        help="Diretório com as saídas do bin Rust.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/darkphotons_benchmark"),
        help="Diretório onde figuras, tabelas e relatório serão gravados.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI dos PNGs.",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Não salvar versões PDF das figuras.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Falhar se faltar um dos cenários esperados.",
    )
    return parser.parse_args()


# =============================================================================
# I/O
# =============================================================================

def scenario_from_text(text: str) -> Optional[str]:
    upper = text.upper()
    for scenario in SCENARIOS:
        if re.search(rf"(^|[^A-Z0-9]){scenario}([^A-Z0-9]|$)", upper):
            return scenario

    if "HADRON" in upper or "PURE_HADRONIC" in upper:
        return "H0"

    return None


def read_summary(model_dir: Path) -> pd.DataFrame:
    path = model_dir / "summary.csv"
    if not path.exists():
        raise FileNotFoundError(f"summary.csv não encontrado em {model_dir}")

    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]

    if "scenario" not in df.columns:
        if "label" in df.columns:
            df["scenario"] = df["label"].astype(str).map(scenario_from_text)
        else:
            df["scenario"] = None

    if "eos_file" not in df.columns:
        raise ValueError(f"{path}: coluna 'eos_file' ausente.")

    df["scenario"] = df.apply(
        lambda row: row["scenario"]
        if isinstance(row["scenario"], str) and row["scenario"] in SCENARIOS
        else scenario_from_text(str(row["eos_file"])),
        axis=1,
    )

    return df


def find_scenario_files(model_dir: Path, summary: pd.DataFrame) -> Dict[str, Path]:
    out: Dict[str, Path] = {}

    for _, row in summary.iterrows():
        scenario = row.get("scenario")
        if scenario not in SCENARIOS:
            continue

        filename = str(row["eos_file"]).strip()
        candidate = model_dir / filename
        if candidate.exists():
            out[scenario] = candidate

    # Fallback por nomes dos arquivos.
    for path in sorted(model_dir.glob("*.dat")):
        scenario = scenario_from_text(path.name)
        if scenario and scenario not in out:
            out[scenario] = path

    return out


def _first_data_line(path: Path) -> Tuple[Optional[str], Optional[str]]:
    """
    Retorna (possible_header, first_data_line).
    """
    possible_header = None
    first_data = None

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("#"):
                content = line.lstrip("#").strip()
                if any(ch.isalpha() for ch in content):
                    possible_header = content
                continue

            # Não procure letras para distinguir cabeçalho de dados: o ``e``
            # da notação científica (por exemplo, 1.0e-3) também é uma letra.
            tokens = [t for t in re.split(r"[\s,]+", line) if t]
            try:
                for token in tokens:
                    float(token)
            except ValueError:
                possible_header = line
                continue

            if tokens:
                first_data = line
                break

    return possible_header, first_data


def read_eos(path: Path, scenario: str) -> pd.DataFrame:
    header, first_data = _first_data_line(path)
    if first_data is None:
        raise ValueError(f"{path}: arquivo sem dados numéricos.")

    # Separador: qualquer whitespace ou vírgula.
    df = pd.read_csv(
        path,
        sep=r"[\s,]+",
        engine="python",
        comment="#",
        header=None,
    )

    # Remove linhas totalmente vazias ou não numéricas.
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if df.empty:
        raise ValueError(f"{path}: nenhuma linha numérica válida.")

    ncol = df.shape[1]

    names = [f"col_{i}" for i in range(ncol)]

    # Se houver cabeçalho compatível, use-o como informação auxiliar apenas.
    if header:
        tokens = [t for t in re.split(r"[\s,]+", header.strip()) if t]
        if len(tokens) == ncol:
            names = [sanitize_column_name(t, i) for i, t in enumerate(tokens)]

    df.columns = names

    # Cria aliases canônicos pelas posições conhecidas.
    for canonical, idx in DEFAULT_COLUMNS.items():
        if idx < ncol:
            df[canonical] = pd.to_numeric(df.iloc[:, idx], errors="coerce")

    df["scenario"] = scenario
    df["source_file"] = str(path)

    return df


def sanitize_column_name(token: str, idx: int) -> str:
    token = token.strip()
    token = token.replace("[", "_").replace("]", "")
    token = token.replace("(", "_").replace(")", "")
    token = token.replace("/", "_per_")
    token = re.sub(r"[^0-9A-Za-z_]+", "_", token)
    token = token.strip("_").lower()
    return token or f"col_{idx}"


# =============================================================================
# Validação/QC
# =============================================================================

def finite_series(df: pd.DataFrame, key: str) -> np.ndarray:
    if key not in df.columns:
        return np.array([], dtype=float)
    arr = pd.to_numeric(df[key], errors="coerce").to_numpy(dtype=float)
    return arr[np.isfinite(arr)]


def qc_eos(model: str, scenario: str, df: pd.DataFrame) -> Dict[str, object]:
    result: Dict[str, object] = {
        "model": model,
        "scenario": scenario,
        "rows": len(df),
        "finite_energy_pressure": False,
        "energy_monotonic": False,
        "pressure_monotonic": False,
        "min_cs2_fd": np.nan,
        "max_cs2_fd": np.nan,
        "causal_fd": False,
        "has_mr": False,
        "status": "FAIL",
    }

    if "energy" not in df.columns or "pressure" not in df.columns:
        result["status"] = "MISSING_EOS_COLUMNS"
        return result

    energy = df["energy"].to_numpy(dtype=float)
    pressure = df["pressure"].to_numpy(dtype=float)

    mask = np.isfinite(energy) & np.isfinite(pressure)
    e = energy[mask]
    p = pressure[mask]

    if len(e) < 3:
        result["status"] = "TOO_FEW_POINTS"
        return result

    result["finite_energy_pressure"] = True

    de = np.diff(e)
    dp = np.diff(p)

    scale_e = max(1.0, float(np.nanmax(np.abs(e))))
    scale_p = max(1.0, float(np.nanmax(np.abs(p))))
    tol_e = 1e-10 * scale_e
    tol_p = 1e-10 * scale_p

    result["energy_monotonic"] = bool(np.all(de >= -tol_e))
    result["pressure_monotonic"] = bool(np.all(dp >= -tol_p))

    good = de > max(tol_e, 1e-14)
    if np.any(good):
        cs2 = dp[good] / de[good]
        cs2 = cs2[np.isfinite(cs2)]
        if len(cs2):
            result["min_cs2_fd"] = float(np.min(cs2))
            result["max_cs2_fd"] = float(np.max(cs2))
            # Pequena tolerância numérica.
            result["causal_fd"] = bool(np.min(cs2) >= -1e-2 and np.max(cs2) <= 1.01)

    if "mass" in df.columns and "radius" in df.columns:
        m = df["mass"].to_numpy(dtype=float)
        r = df["radius"].to_numpy(dtype=float)
        result["has_mr"] = bool(np.any(np.isfinite(m) & np.isfinite(r) & (m > 0) & (r > 0)))

    if (
        result["finite_energy_pressure"]
        and result["energy_monotonic"]
        and result["pressure_monotonic"]
        and result["causal_fd"]
    ):
        result["status"] = "PASS"
    else:
        result["status"] = "CHECK"

    return result


# =============================================================================
# Métricas
# =============================================================================

def stable_mr_branch(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Retorna (M, R, indices_originais) no ramo estável até o primeiro máximo de M.

    Se M/R não estiverem disponíveis, retorna arrays vazios.
    """
    if "mass" not in df.columns or "radius" not in df.columns:
        return np.array([]), np.array([]), np.array([], dtype=int)

    m = df["mass"].to_numpy(dtype=float)
    r = df["radius"].to_numpy(dtype=float)
    idx = np.arange(len(df))

    mask = np.isfinite(m) & np.isfinite(r) & (m > 0) & (r > 0)
    m = m[mask]
    r = r[mask]
    idx = idx[mask]

    if len(m) < 3:
        return np.array([]), np.array([]), np.array([], dtype=int)

    imax = int(np.argmax(m))
    return m[: imax + 1], r[: imax + 1], idx[: imax + 1]


def interpolate_radius_at_mass(df: pd.DataFrame, target_mass: float) -> float:
    m, r, _ = stable_mr_branch(df)
    if len(m) < 2:
        return np.nan

    # Remove duplicatas e ordena em massa.
    order = np.argsort(m)
    m = m[order]
    r = r[order]

    m_unique, unique_idx = np.unique(m, return_index=True)
    r_unique = r[unique_idx]

    if target_mass < m_unique.min() or target_mass > m_unique.max():
        return np.nan

    return float(np.interp(target_mass, m_unique, r_unique))


def max_mass_metrics(df: pd.DataFrame) -> Tuple[float, float]:
    m, r, _ = stable_mr_branch(df)
    if len(m) == 0:
        return np.nan, np.nan

    i = int(np.argmax(m))
    return float(m[i]), float(r[i])


def dark_energy_fraction(df: pd.DataFrame) -> np.ndarray:
    if "energy" not in df.columns:
        return np.full(len(df), np.nan)

    if "energy_chi_kin" not in df.columns or "energy_x" not in df.columns:
        return np.full(len(df), np.nan)

    total = df["energy"].to_numpy(dtype=float)
    dark = (
        df["energy_chi_kin"].to_numpy(dtype=float)
        + df["energy_x"].to_numpy(dtype=float)
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        frac = dark / total

    frac[~np.isfinite(frac)] = np.nan
    return frac


def interpolate_pressure_at_energy(
    reference: pd.DataFrame,
    other: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Retorna energy_grid e P_other/P_reference no domínio comum.
    """
    e0 = reference["energy"].to_numpy(dtype=float)
    p0 = reference["pressure"].to_numpy(dtype=float)
    e1 = other["energy"].to_numpy(dtype=float)
    p1 = other["pressure"].to_numpy(dtype=float)

    mask0 = np.isfinite(e0) & np.isfinite(p0) & (p0 > 0)
    mask1 = np.isfinite(e1) & np.isfinite(p1)

    e0, p0 = e0[mask0], p0[mask0]
    e1, p1 = e1[mask1], p1[mask1]

    if len(e0) < 3 or len(e1) < 3:
        return np.array([]), np.array([])

    o0 = np.argsort(e0)
    o1 = np.argsort(e1)
    e0, p0 = e0[o0], p0[o0]
    e1, p1 = e1[o1], p1[o1]

    emin = max(float(e0.min()), float(e1.min()))
    emax = min(float(e0.max()), float(e1.max()))

    if emax <= emin:
        return np.array([]), np.array([])

    grid = np.linspace(emin, emax, 400)
    pref = np.interp(grid, e0, p0)
    pother = np.interp(grid, e1, p1)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = pother / pref

    mask = np.isfinite(ratio) & (pref > 0)
    return grid[mask], ratio[mask]


def collect_metrics(
    all_data: Dict[str, Dict[str, pd.DataFrame]]
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    for model, scenarios in all_data.items():
        for scenario, df in scenarios.items():
            mmax, rmax = max_mass_metrics(df)
            r14 = interpolate_radius_at_mass(df, 1.4)
            r20 = interpolate_radius_at_mass(df, 2.0)

            energy = finite_series(df, "energy")
            pressure = finite_series(df, "pressure")
            nb = finite_series(df, "nb_over_n0")
            frac = dark_energy_fraction(df)
            frac_finite = frac[np.isfinite(frac)]

            row = {
                "model": model,
                "scenario": scenario,
                "label": SCENARIO_LABELS.get(scenario, scenario),
                "rows": len(df),
                "energy_max_MeV_fm3": float(np.max(energy)) if len(energy) else np.nan,
                "pressure_max_MeV_fm3": float(np.max(pressure)) if len(pressure) else np.nan,
                "nb_over_n0_max": float(np.max(nb)) if len(nb) else np.nan,
                "Mmax_Msun": mmax,
                "R_at_Mmax_km": rmax,
                "R_1p4_km": r14,
                "R_2p0_km": r20,
                "dark_energy_fraction_max": (
                    float(np.max(frac_finite)) if len(frac_finite) else np.nan
                ),
                "dark_energy_fraction_last": (
                    float(frac_finite[-1]) if len(frac_finite) else np.nan
                ),
            }

            for key in ("epsilon", "m_chi", "m_x", "g_d", "y_chi"):
                arr = finite_series(df, key)
                row[key] = float(np.nanmedian(arr)) if len(arr) else np.nan

            rows.append(row)

    columns = [
        "model",
        "scenario",
        "label",
        "rows",
        "energy_max_MeV_fm3",
        "pressure_max_MeV_fm3",
        "nb_over_n0_max",
        "Mmax_Msun",
        "R_at_Mmax_km",
        "R_1p4_km",
        "R_2p0_km",
        "dark_energy_fraction_max",
        "dark_energy_fraction_last",
        "epsilon",
        "m_chi",
        "m_x",
        "g_d",
        "y_chi",
    ]
    return pd.DataFrame(rows, columns=columns)


# =============================================================================
# Plotting
# =============================================================================

def save_figure(fig, basename: Path, dpi: int, save_pdf: bool) -> None:
    basename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(basename.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(basename.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_eos(
    model: str,
    scenarios: Dict[str, pd.DataFrame],
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    any_curve = False

    for scenario in SCENARIOS:
        df = scenarios.get(scenario)
        if df is None or "energy" not in df or "pressure" not in df:
            continue

        energy = pd.to_numeric(df["energy"], errors="coerce").to_numpy(dtype=float)
        pressure = pd.to_numeric(df["pressure"], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(energy) & np.isfinite(pressure)
        if not np.any(finite):
            continue

        any_curve = True
        ax.plot(
            energy[finite],
            pressure[finite],
            linestyle=LINESTYLES[scenario],
            linewidth=1.8,
            label=SCENARIO_LABELS[scenario],
        )

    if not any_curve:
        plt.close(fig)
        return

    ax.set_xlabel(r"$\varepsilon\;[\mathrm{MeV/fm^3}]$")
    ax.set_ylabel(r"$P\;[\mathrm{MeV/fm^3}]$")
    ax.set_title(f"Equação de estado — {model}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    save_figure(fig, out_dir / f"eos_{model}", dpi, save_pdf)


def plot_mr(
    model: str,
    scenarios: Dict[str, pd.DataFrame],
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> bool:
    any_curve = False
    fig, ax = plt.subplots(figsize=(6.4, 4.8))

    for scenario in SCENARIOS:
        df = scenarios.get(scenario)
        if df is None:
            continue

        m, r, _ = stable_mr_branch(df)
        if len(m) < 2:
            continue

        any_curve = True
        ax.plot(
            r,
            m,
            linestyle=LINESTYLES[scenario],
            linewidth=1.8,
            label=SCENARIO_LABELS[scenario],
        )

    if not any_curve:
        plt.close(fig)
        return False

    ax.set_xlabel(r"$R\;[\mathrm{km}]$")
    ax.set_ylabel(r"$M\;[M_\odot]$")
    ax.set_title(f"Relação massa–raio — {model}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    save_figure(fig, out_dir / f"mr_{model}", dpi, save_pdf)
    return True


def plot_dark_energy_fraction(
    model: str,
    scenarios: Dict[str, pd.DataFrame],
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    plotted = False

    for scenario in ("S1", "S2", "S3"):
        df = scenarios.get(scenario)
        if df is None or "nb_over_n0" not in df:
            continue

        frac = dark_energy_fraction(df)
        if not np.any(np.isfinite(frac)):
            continue

        plotted = True
        ax.plot(
            df["nb_over_n0"],
            100.0 * frac,
            linestyle=LINESTYLES[scenario],
            linewidth=1.8,
            label=SCENARIO_LABELS[scenario],
        )

    if not plotted:
        plt.close(fig)
        return

    ax.set_xlabel(r"$n_B/n_0$")
    ax.set_ylabel(r"$100\,\varepsilon_{\rm dark}/\varepsilon_{\rm tot}\;[\%]$")
    ax.set_title(f"Fração local de energia escura — {model}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    save_figure(
        fig,
        out_dir / f"dark_energy_fraction_{model}",
        dpi,
        save_pdf,
    )


def plot_dark_quantity(
    model: str,
    scenarios: Dict[str, pd.DataFrame],
    key: str,
    ylabel: str,
    title: str,
    filename: str,
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    plotted = False

    for scenario in ("S1", "S2", "S3"):
        df = scenarios.get(scenario)
        if df is None or "nb_over_n0" not in df or key not in df:
            continue

        x = df["nb_over_n0"].to_numpy(dtype=float)
        y = df[key].to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)

        if np.count_nonzero(mask) < 2:
            continue

        plotted = True
        ax.plot(
            x[mask],
            y[mask],
            linestyle=LINESTYLES[scenario],
            linewidth=1.8,
            label=SCENARIO_LABELS[scenario],
        )

    if not plotted:
        plt.close(fig)
        return

    ax.set_xlabel(r"$n_B/n_0$")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} — {model}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    save_figure(fig, out_dir / f"{filename}_{model}", dpi, save_pdf)


def plot_pressure_ratio(
    model: str,
    scenarios: Dict[str, pd.DataFrame],
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> None:
    reference = scenarios.get("H0")
    if reference is None:
        return

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    plotted = False

    for scenario in ("S1", "S2", "S3"):
        df = scenarios.get(scenario)
        if df is None:
            continue

        energy, ratio = interpolate_pressure_at_energy(reference, df)
        if len(energy) < 2:
            continue

        plotted = True
        ax.plot(
            energy,
            ratio,
            linestyle=LINESTYLES[scenario],
            linewidth=1.8,
            label=SCENARIO_LABELS[scenario],
        )

    if not plotted:
        plt.close(fig)
        return

    ax.axhline(1.0, linewidth=1.0)
    ax.set_xlabel(r"$\varepsilon\;[\mathrm{MeV/fm^3}]$")
    ax.set_ylabel(r"$P_{\rm cenário}/P_{\rm H0}$")
    ax.set_title(f"Modificação relativa da pressão — {model}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    save_figure(fig, out_dir / f"pressure_ratio_{model}", dpi, save_pdf)


def plot_metric_comparison(
    metrics: pd.DataFrame,
    value_col: str,
    ylabel: str,
    title: str,
    filename: str,
    out_dir: Path,
    dpi: int,
    save_pdf: bool,
) -> None:
    data = metrics[np.isfinite(pd.to_numeric(metrics[value_col], errors="coerce"))].copy()
    if data.empty:
        return

    # Cada combinação modelo/cenário é uma categoria.
    data["category"] = data["model"] + " " + data["scenario"]
    x = np.arange(len(data))

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(x, data[value_col].to_numpy(dtype=float))
    ax.set_xticks(x)
    ax.set_xticklabels(data["category"], rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)

    save_figure(fig, out_dir / filename, dpi, save_pdf)


# =============================================================================
# Tabelas e relatório
# =============================================================================

def build_input_scenarios(summary_by_model: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []

    for model, df in summary_by_model.items():
        for _, row in df.iterrows():
            scenario = row.get("scenario")
            if scenario not in SCENARIOS:
                continue

            record = {
                "model": model,
                "scenario": scenario,
                "eos_file": row.get("eos_file", ""),
            }

            for col in (
                "description",
                "epsilon",
                "m_chi_MeV",
                "m_x_MeV",
                "g_d",
                "y_chi",
                "kf_chi_ref_MeV",
                "nB_ref_fm-3",
                "b_field_G",
                "status",
            ):
                if col in row.index:
                    record[col] = row[col]

            rows.append(record)

    return pd.DataFrame(rows)


def write_latex_table(metrics: pd.DataFrame, path: Path) -> None:
    cols = [
        "model",
        "scenario",
        "Mmax_Msun",
        "R_at_Mmax_km",
        "R_1p4_km",
        "dark_energy_fraction_max",
    ]

    available = [c for c in cols if c in metrics.columns]
    table = metrics[available].copy()

    rename = {
        "model": "Modelo",
        "scenario": "Cenário",
        "Mmax_Msun": r"$M_{\max}/M_\odot$",
        "R_at_Mmax_km": r"$R(M_{\max})$ [km]",
        "R_1p4_km": r"$R_{1.4}$ [km]",
        "dark_energy_fraction_max": r"$f_{\rm dark}^{\max}$",
    }
    table = table.rename(columns=rename)

    latex = table.to_latex(
        index=False,
        float_format=lambda x: f"{x:.5g}",
        escape=False,
        na_rep="--",
        caption=(
            "Propriedades dos quatro cenários de benchmark para GM1 e GM3. "
            "As quantidades estelares são apresentadas apenas quando as colunas "
            "massa--raio estão disponíveis nos arquivos de saída."
        ),
        label="tab:darkphotons_benchmark_metrics",
    )
    path.write_text(latex, encoding="utf-8")


def write_report(
    path: Path,
    metrics: pd.DataFrame,
    qc: pd.DataFrame,
    missing: List[str],
    mr_available: Dict[str, bool],
) -> None:
    lines: List[str] = []

    lines.append("ANÁLISE DOS BENCHMARKS DE FÓTON ESCURO")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Cenários: H0, S1, S2, S3")
    lines.append("Modelos hadrônicos: GM1, GM3")
    lines.append("")

    if missing:
        lines.append("ARQUIVOS/CENÁRIOS AUSENTES")
        lines.append("-" * 72)
        lines.extend(f"- {item}" for item in missing)
        lines.append("")

    lines.append("CONTROLE DE QUALIDADE")
    lines.append("-" * 72)
    if qc.empty:
        lines.append("Nenhum arquivo pôde ser validado.")
    else:
        for _, row in qc.iterrows():
            lines.append(
                f"{row['model']} {row['scenario']}: "
                f"{row['status']}; rows={row['rows']}; "
                f"c_s^2=[{fmt(row['min_cs2_fd'])}, {fmt(row['max_cs2_fd'])}]"
            )
    lines.append("")

    lines.append("MÉTRICAS PRINCIPAIS")
    lines.append("-" * 72)
    if not metrics.empty:
        for _, row in metrics.iterrows():
            lines.append(
                f"{row['model']} {row['scenario']}: "
                f"Mmax={fmt(row['Mmax_Msun'])} Msun; "
                f"R(Mmax)={fmt(row['R_at_Mmax_km'])} km; "
                f"R1.4={fmt(row['R_1p4_km'])} km; "
                f"f_dark,max={fmt(row['dark_energy_fraction_max'])}"
            )
    lines.append("")

    lines.append("DISPONIBILIDADE MASSA–RAIO")
    lines.append("-" * 72)
    for model in MODELS:
        lines.append(
            f"{model}: {'sim' if mr_available.get(model, False) else 'não'}"
        )
    lines.append("")

    lines.append("INTERPRETAÇÃO RECOMENDADA")
    lines.append("-" * 72)
    lines.append(
        "1. Compare H0 com S1/S2 para isolar o regime de mediador pesado, "
        "no qual a contribuição vetorial tende a ser suprimida."
    )
    lines.append(
        "2. Compare H0 com S3 para quantificar o impacto do mediador leve "
        "na rigidez da EOS."
    )
    lines.append(
        "3. O gráfico P_cenário/P_H0 mostra diretamente enrijecimento (>1) "
        "ou amolecimento (<1) a energia fixa."
    )
    lines.append(
        "4. A fração de energia escura deve ser interpretada como "
        "(epsilon_chi^kin + epsilon_X)/epsilon_total; ela não é igual a Y_chi."
    )
    lines.append(
        "5. Não interpretar estes benchmarks como reprodução numérica direta "
        "do artigo de portal Z': o modelo atual usa kinetic mixing com a "
        "corrente eletromagnética, enquanto o artigo de referência acopla o "
        "mediador diretamente a quarks/núcleons."
    )
    lines.append(
        "6. Os testes de causalidade aqui são diagnósticos por diferenças "
        "finitas da tabela P(epsilon); para publicação, confirme com o "
        "validador Rust dedicado."
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fmt(value: object) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(x):
        return "--"
    return f"{x:.6g}"


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    args = parse_args()

    input_root = args.input
    output_root = args.output
    figures_dir = output_root / "figures"
    tables_dir = output_root / "tables"

    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    all_data: Dict[str, Dict[str, pd.DataFrame]] = {}
    summary_by_model: Dict[str, pd.DataFrame] = {}
    qc_rows: List[Dict[str, object]] = []
    missing: List[str] = []

    for model in MODELS:
        model_dir = input_root / model

        if not model_dir.exists():
            message = f"{model}: diretório ausente ({model_dir})"
            if args.strict:
                raise FileNotFoundError(message)
            missing.append(message)
            continue

        try:
            summary = read_summary(model_dir)
        except Exception as exc:
            if args.strict:
                raise
            missing.append(f"{model}: falha ao ler summary.csv: {exc}")
            continue

        summary_by_model[model] = summary
        files = find_scenario_files(model_dir, summary)

        model_data: Dict[str, pd.DataFrame] = {}

        for scenario in SCENARIOS:
            path = files.get(scenario)

            if path is None:
                message = f"{model} {scenario}: arquivo EOS não encontrado"
                if args.strict:
                    raise FileNotFoundError(message)
                missing.append(message)
                continue

            try:
                df = read_eos(path, scenario)
            except Exception as exc:
                if args.strict:
                    raise
                missing.append(f"{model} {scenario}: erro de leitura: {exc}")
                continue

            model_data[scenario] = df
            qc_rows.append(qc_eos(model, scenario, df))

        if model_data:
            all_data[model] = model_data

    if not all_data:
        raise SystemExit(
            f"Nenhuma EOS com dados numéricos foi encontrada em {input_root}. "
            "Verifique --input e regenere arquivos que contenham apenas o cabeçalho."
        )

    # -------------------------------------------------------------------------
    # Figuras por modelo
    # -------------------------------------------------------------------------

    mr_available: Dict[str, bool] = {}

    for model, scenarios in all_data.items():
        plot_eos(
            model,
            scenarios,
            figures_dir,
            args.dpi,
            not args.no_pdf,
        )

        mr_available[model] = plot_mr(
            model,
            scenarios,
            figures_dir,
            args.dpi,
            not args.no_pdf,
        )

        plot_dark_energy_fraction(
            model,
            scenarios,
            figures_dir,
            args.dpi,
            not args.no_pdf,
        )

        plot_dark_quantity(
            model,
            scenarios,
            key="kf_chi",
            ylabel=r"$k_{F\chi}\;[\mathrm{MeV}]$",
            title=r"Momento de Fermi da matéria escura",
            filename="kf_chi",
            out_dir=figures_dir,
            dpi=args.dpi,
            save_pdf=not args.no_pdf,
        )

        plot_dark_quantity(
            model,
            scenarios,
            key="x0",
            ylabel=r"$X_0\;[\mathrm{MeV}]$",
            title=r"Campo médio do fóton escuro",
            filename="x0",
            out_dir=figures_dir,
            dpi=args.dpi,
            save_pdf=not args.no_pdf,
        )

        plot_pressure_ratio(
            model,
            scenarios,
            figures_dir,
            args.dpi,
            not args.no_pdf,
        )

    # -------------------------------------------------------------------------
    # Métricas e comparações globais
    # -------------------------------------------------------------------------

    metrics = collect_metrics(all_data)
    if metrics.empty:
        raise SystemExit(
            "Nenhuma métrica pôde ser calculada a partir das EOS carregadas."
        )
    metrics = metrics.sort_values(["model", "scenario"]).reset_index(drop=True)

    qc = pd.DataFrame(qc_rows)
    if not qc.empty:
        qc = qc.sort_values(["model", "scenario"]).reset_index(drop=True)

    metrics.to_csv(
        tables_dir / "benchmark_metrics.csv",
        index=False,
        float_format="%.10e",
    )

    qc.to_csv(
        tables_dir / "qc_report.csv",
        index=False,
        float_format="%.10e",
    )

    input_scenarios = build_input_scenarios(summary_by_model)
    input_scenarios.to_csv(
        tables_dir / "input_scenarios.csv",
        index=False,
    )

    write_latex_table(
        metrics,
        tables_dir / "benchmark_metrics.tex",
    )

    plot_metric_comparison(
        metrics,
        value_col="Mmax_Msun",
        ylabel=r"$M_{\max}\;[M_\odot]$",
        title="Massa máxima dos benchmarks",
        filename="max_mass_comparison",
        out_dir=figures_dir,
        dpi=args.dpi,
        save_pdf=not args.no_pdf,
    )

    plot_metric_comparison(
        metrics,
        value_col="R_1p4_km",
        ylabel=r"$R_{1.4}\;[\mathrm{km}]$",
        title=r"Raio da estrela de $1.4\,M_\odot$",
        filename="radius_1p4_comparison",
        out_dir=figures_dir,
        dpi=args.dpi,
        save_pdf=not args.no_pdf,
    )

    write_report(
        output_root / "analysis_report.txt",
        metrics,
        qc,
        missing,
        mr_available,
    )

    # -------------------------------------------------------------------------
    # Resumo no terminal
    # -------------------------------------------------------------------------

    print("\nAnálise concluída.")
    print(f"Entrada : {input_root}")
    print(f"Saída   : {output_root}")
    print(f"Figuras : {figures_dir}")
    print(f"Tabelas : {tables_dir}")

    if missing:
        print("\nAvisos:")
        for item in missing:
            print(f"  - {item}")

    if not qc.empty:
        print("\nQC:")
        print(qc[["model", "scenario", "status", "rows", "min_cs2_fd", "max_cs2_fd"]].to_string(index=False))

    if not metrics.empty:
        display_cols = [
            "model",
            "scenario",
            "Mmax_Msun",
            "R_at_Mmax_km",
            "R_1p4_km",
            "dark_energy_fraction_max",
        ]
        print("\nMétricas:")
        print(metrics[display_cols].to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
