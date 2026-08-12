#!/usr/bin/env python3
"""
Análise detalhada do início dos efeitos de CSI no dataset nlem_log.

Este script identifica onde os efeitos do CSI (campo magnético effectivo relativo)
começam a impactar significativamente cada grandeza física na equação de estado,
em todo o dataset nlem_log para diferentes modelos (GM1, GM3) e campos magnéticos.

Analisa:
- Populações de partículas (onset thresholds)
- Potenciais químicos
- Campos magnéticos efetivos
- Energia de magnetização
- Velocidade do som
- Massa e raio estelares
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


# Índices de colunas conforme nlem_log.py
COL_NB = 0
COL_EPS = 1
COL_P = 2
COL_NE = 3
COL_NMU = 4
COL_NN = 5
COL_NP = 6
COL_NL0 = 7      # Lambda
COL_NSM = 8      # Sigma minus
COL_NS0 = 9      # Sigma zero
COL_NSP = 10     # Sigma plus
COL_NXM = 11     # Xi minus
COL_NX0 = 12     # Xi zero

COL_MEFF = 16    # Massa efetiva
COL_MUN = 17     # Potencial químico do nêutron
COL_MUE = 18     # Potencial químico do elétron
COL_EMAG = 19    # Energia de magnetização
COL_MU_TOTAL = 20  # Potencial químico fermiônico total por bárion / M_N


# Mapeamento de coluna para nome legível
COLUMN_NAMES = {
    COL_NB: "n_B",
    COL_EPS: "ε",
    COL_P: "P",
    COL_NE: "n_e",
    COL_NMU: "n_μ",
    COL_NN: "n_n",
    COL_NP: "n_p",
    COL_NL0: "n_Λ",
    COL_NSM: "n_Σ⁻",
    COL_NS0: "n_Σ⁰",
    COL_NSP: "n_Σ⁺",
    COL_NXM: "n_Ξ⁻",
    COL_NX0: "n_Ξ⁰",
    COL_MEFF: "M_eff",
    COL_MUN: "μ_n",
    COL_MUE: "μ_e",
    COL_EMAG: "E_mag",
    COL_MU_TOTAL: "mu_total_over_mN",
}

# Populações de hiperanos
HYPERON_COLUMNS = (COL_NL0, COL_NSM, COL_NS0, COL_NSP, COL_NXM, COL_NX0)


@dataclass
class EoSDataset:
    """Representa um arquivo eos.dat carregado"""
    model: str
    b_label: str
    b_value: float
    topology: str
    csi: float
    log_csi: float
    data: np.ndarray
    file_path: Path


@dataclass
class ColumnOnsetAnalysis:
    """Análise de onset para uma coluna específica"""
    column_idx: int
    column_name: str
    log_csi_values: List[float] = field(default_factory=list)
    onset_log_csi: Optional[float] = None
    onset_csi: Optional[float] = None
    onset_density_nb: Optional[float] = None
    onset_relative_change: Optional[float] = None
    values_at_onset: List[float] = field(default_factory=list)
    onset_location: str = ""  # "superfície", "centro", "intermediário"


@dataclass
class DatasetOnsetReport:
    """Relatório completo de onset para um dataset"""
    model: str
    b_label: str
    b_value: float
    topology: str
    log_csi: float
    csi: float
    column_analyses: Dict[int, ColumnOnsetAnalysis] = field(default_factory=dict)
    hyperon_onset_log_csi: Optional[float] = None
    first_hyperon_appearance: Optional[str] = None


def parse_metadata_from_path(path: Path) -> Optional[Tuple[str, str, float, str, float, float]]:
    """Extrai metadados do caminho do arquivo"""
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


def _safe_float(text: str) -> Optional[float]:
    """Conversão segura para float"""
    try:
        return float(text)
    except (ValueError, TypeError):
        return None


def load_eos(path: Path) -> Optional[np.ndarray]:
    """Carrega um arquivo eos.dat"""
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


def discover_datasets(input_root: Path) -> List[EoSDataset]:
    """Descobre todos os datasets em um diretório"""
    datasets: List[EoSDataset] = []

    for eos_file in sorted(input_root.rglob("eos.dat")):
        meta = parse_metadata_from_path(eos_file)
        if meta is None:
            continue

        data = load_eos(eos_file)
        if data is None:
            continue

        model, b_label, b_value, topology, csi, log_csi = meta
        datasets.append(
            EoSDataset(
                model=model,
                b_label=b_label,
                b_value=b_value,
                topology=topology,
                csi=csi,
                log_csi=log_csi,
                data=data,
                file_path=eos_file,
            )
        )

    return datasets


def detect_onset_for_column(
    log_csi_sorted: List[float],
    column_data: Dict[float, np.ndarray],
    column_idx: int,
    threshold_change: float = 0.05,
) -> Optional[ColumnOnsetAnalysis]:
    """
    Detecta o onset de efeitos para uma coluna específica.

    Analisa mudanças entre valores consecutivos de log_csi.
    O onset é quando uma mudança relativa significativa (threshold_change) ocorre.
    """
    if len(log_csi_sorted) < 2:
        return None

    analysis = ColumnOnsetAnalysis(
        column_idx=column_idx,
        column_name=COLUMN_NAMES.get(column_idx, f"col_{column_idx}"),
    )

    log_csi_sorted = sorted(log_csi_sorted)
    analysis.log_csi_values = log_csi_sorted

    # Calcular valores médios em cada ponto de CSI
    means_by_log_csi = {}
    for log_csi in log_csi_sorted:
        col_data = column_data[log_csi]
        col_data = col_data[np.isfinite(col_data)]
        if col_data.size > 0:
            means_by_log_csi[log_csi] = np.mean(col_data)
        else:
            means_by_log_csi[log_csi] = np.nan

    # Detectar mudanças relativas
    for i in range(1, len(log_csi_sorted)):
        prev_log_csi = log_csi_sorted[i-1]
        curr_log_csi = log_csi_sorted[i]

        prev_val = means_by_log_csi[prev_log_csi]
        curr_val = means_by_log_csi[curr_log_csi]

        if not np.isfinite(prev_val) or not np.isfinite(curr_val):
            continue

        if prev_val == 0:
            continue

        relative_change = abs((curr_val - prev_val) / prev_val)

        # Se primeira mudança significativa ou mudança em população (non-zero onset)
        if relative_change > threshold_change:
            analysis.onset_log_csi = curr_log_csi
            analysis.onset_csi = 10 ** curr_log_csi
            analysis.onset_relative_change = relative_change
            analysis.values_at_onset = [prev_val, curr_val]
            break

    return analysis


def analyze_hyperon_onset(datasets_by_b_model: Dict[str, List[EoSDataset]]) -> Dict[str, Optional[float]]:
    """Analisa o onset de hiperanos (primeira aparição em densidade)."""
    hyperon_onsets = {}

    for key, datasets in datasets_by_b_model.items():
        sorted_datasets = sorted(datasets, key=lambda d: d.log_csi)
        onset_log_csi = _find_first_hyperon_onset(sorted_datasets)
        hyperon_onsets[key] = onset_log_csi

    return hyperon_onsets


def _find_first_hyperon_onset(sorted_datasets: List[EoSDataset]) -> Optional[float]:
    """Encontra o primeiro onset de hiperão em uma lista ordenada."""
    for ds in sorted_datasets:
        for col_idx in HYPERON_COLUMNS:
            hyperon_pop = ds.data[:, col_idx]
            if np.any(hyperon_pop > 1e-8):
                return ds.log_csi
    return None


def analyze_dataset_group(
    datasets_by_b_model: Dict[str, List[EoSDataset]],
) -> Dict[str, List[DatasetOnsetReport]]:
    """Analisa um grupo de datasets."""
    reports = defaultdict(list)

    for key, datasets in datasets_by_b_model.items():
        sorted_datasets = sorted(datasets, key=lambda d: d.log_csi)
        data_by_column = _build_column_data_map(sorted_datasets)
        log_csi_values = sorted([ds.log_csi for ds in sorted_datasets])

        for ds in sorted_datasets:
            report = _build_dataset_report(ds, log_csi_values, data_by_column)
            reports[key].append(report)

    return dict(reports)


def _build_column_data_map(
    sorted_datasets: List[EoSDataset],
) -> Dict[int, Dict[float, np.ndarray]]:
    """Constrói mapa de dados por coluna e log_csi."""
    data_by_column: Dict[int, Dict[float, np.ndarray]] = defaultdict(dict)
    for ds in sorted_datasets:
        for col_idx in range(ds.data.shape[1]):
            data_by_column[col_idx][ds.log_csi] = ds.data[:, col_idx]
    return dict(data_by_column)


def _build_dataset_report(
    ds: EoSDataset,
    log_csi_values: List[float],
    data_by_column: Dict[int, Dict[float, np.ndarray]],
) -> DatasetOnsetReport:
    """Constrói relatório de onset para um dataset."""
    report = DatasetOnsetReport(
        model=ds.model,
        b_label=ds.b_label,
        b_value=ds.b_value,
        topology=ds.topology,
        log_csi=ds.log_csi,
        csi=ds.csi,
    )

    important_columns = [
        COL_NE, COL_NMU, COL_NN, COL_NP,
        COL_NL0, COL_NSM, COL_NS0, COL_NSP, COL_NXM, COL_NX0,
        COL_MEFF, COL_MUN, COL_MUE, COL_EMAG, COL_MU_TOTAL
    ]

    for col_idx in important_columns:
        if col_idx < ds.data.shape[1]:
            analysis = detect_onset_for_column(
                log_csi_values,
                data_by_column[col_idx],
                col_idx,
            )
            if analysis:
                report.column_analyses[col_idx] = analysis

    _detect_hyperon_onsets_in_report(ds, report)
    return report


def _detect_hyperon_onsets_in_report(
    ds: EoSDataset,
    report: DatasetOnsetReport,
) -> None:
    """Detecta onsets de hiperões no relatório."""
    for col_idx in HYPERON_COLUMNS:
        hyperon_data = ds.data[:, col_idx]
        if np.any(hyperon_data > 1e-8) and report.hyperon_onset_log_csi is None:
            report.hyperon_onset_log_csi = ds.log_csi
            report.first_hyperon_appearance = COLUMN_NAMES.get(col_idx, f"col_{col_idx}")


def write_detailed_report(
    reports: Dict[str, List[DatasetOnsetReport]],
    output_file: Path,
) -> None:
    """Escreve relatório detalhado em arquivo CSV"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # Cabeçalho
        writer.writerow([
            "Modelo",
            "B (Gauss)",
            "B_label",
            "Topologia",
            "log10(ξ)",
            "ξ",
            "Coluna",
            "Nome Coluna",
            "Onset log10(ξ)",
            "Onset ξ",
            "Mudança Relativa (%)",
            "Caminho Arquivo",
        ])

        for key in sorted(reports.keys()):
            for report in sorted(reports[key], key=lambda r: r.log_csi):
                for col_idx, analysis in sorted(report.column_analyses.items()):
                    if analysis.onset_log_csi is not None:
                        writer.writerow([
                            report.model,
                            f"{report.b_value:.3e}",
                            report.b_label,
                            report.topology,
                            f"{report.log_csi:.2f}",
                            f"{report.csi:.3e}",
                            col_idx,
                            analysis.column_name,
                            f"{analysis.onset_log_csi:.2f}",
                            f"{analysis.onset_csi:.3e}",
                            f"{analysis.onset_relative_change*100:.2f}" if analysis.onset_relative_change else "N/A",
                            str(report.model),
                        ])


def write_summary_report(
    reports: Dict[str, List[DatasetOnsetReport]],
    output_file: Path,
) -> None:
    """Escreve sumário de onsets por modelo e campo magnético"""
    summary = defaultdict(lambda: defaultdict(list))

    # Agregar dados
    for key, report_list in reports.items():
        for report in report_list:
            model_b_key = f"{report.model}_{report.b_label}"
            for col_idx, analysis in report.column_analyses.items():
                if analysis.onset_log_csi is not None:
                    summary[model_b_key][col_idx].append(analysis.onset_log_csi)

    # Escrever sumário
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        writer.writerow([
            "Modelo_B",
            "Coluna",
            "Nome Coluna",
            "Onset log10(ξ) Mínimo",
            "Onset log10(ξ) Máximo",
            "Onset log10(ξ) Médio",
            "Contagem",
        ])

        for model_b in sorted(summary.keys()):
            for col_idx in sorted(summary[model_b].keys()):
                onsets = summary[model_b][col_idx]
                col_name = COLUMN_NAMES.get(col_idx, f"col_{col_idx}")
                writer.writerow([
                    model_b,
                    col_idx,
                    col_name,
                    f"{min(onsets):.2f}",
                    f"{max(onsets):.2f}",
                    f"{np.mean(onsets):.2f}",
                    len(onsets),
                ])


def plot_onset_analysis(
    reports: Dict[str, List[DatasetOnsetReport]],
    output_dir: Path,
) -> None:
    """Gera gráficos da análise de onset"""
    output_dir.mkdir(parents=True, exist_ok=True)

    all_onsets = _collect_onset_data(reports)
    labels = sorted(all_onsets.keys())
    data = [all_onsets[label] for label in labels]

    _, ax = plt.subplots(figsize=(14, 8))

    bp = ax.boxplot(data, patch_artist=True)
    for idx, patch in enumerate(bp['boxes']):
        patch.set_facecolor('lightblue')
        ax.set_xticklabels(labels, rotation=45, ha='right')

    ax.set_xlabel("Modelo e Campo Magnético", fontsize=12)
    ax.set_ylabel(r"$\log_{10}(\xi)$ no Onset", fontsize=12)
    ax.set_title("Distribuição de Onsets de CSI por Grandeza Física", fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "onset_distribution.png", dpi=150)
    plt.close()

    print(f"Gráfico salvo em: {output_dir / 'onset_distribution.png'}")


def _collect_onset_data(
    reports: Dict[str, List[DatasetOnsetReport]],
) -> Dict[str, List[float]]:
    """Coleta dados de onset por modelo e coluna."""
    all_onsets: Dict[str, List[float]] = defaultdict(list)

    for key, report_list in reports.items():
        model_b = f"{report_list[0].model}_{report_list[0].b_label}"

        for report in report_list:
            for col_idx, analysis in report.column_analyses.items():
                if analysis.onset_log_csi is not None:
                    label = f"{model_b}_{COLUMN_NAMES.get(col_idx)}"
                    all_onsets[label].append(analysis.onset_log_csi)

    return all_onsets


def main():
    input_root = Path("output/nlem_log")
    output_root = Path("results/nlem_log/csi_onset_analysis")
    output_root.mkdir(parents=True, exist_ok=True)

    print("📊 Análise de Onset de CSI no Dataset NLEM_LOG")
    print("=" * 60)

    # 1. Descobrir datasets
    print("\n🔍 Descobrindo datasets...")
    all_datasets = discover_datasets(input_root)
    print(f"   ✓ {len(all_datasets)} datasets encontrados")

    # 2. Agrupar por modelo e campo magnético
    datasets_by_model_b = defaultdict(list)
    for ds in all_datasets:
        key = f"{ds.model}_{ds.b_label}"
        datasets_by_model_b[key].append(ds)

    print(f"   ✓ Agrupados em {len(datasets_by_model_b)} grupos")

    # 3. Analisar onset para cada grupo
    print("\n📈 Analisando onsets de CSI...")
    reports = analyze_dataset_group(datasets_by_model_b)

    total_analyses = sum(
        len(report.column_analyses)
        for report_list in reports.values()
        for report in report_list
    )
    print(f"   ✓ {total_analyses} análises de onset concluídas")

    # 4. Gerar relatórios
    print("\n📄 Gerando relatórios...")

    detailed_report_path = output_root / "detailed_onset_analysis.csv"
    write_detailed_report(reports, detailed_report_path)
    print(f"   ✓ Relatório detalhado: {detailed_report_path}")

    summary_report_path = output_root / "onset_summary.csv"
    write_summary_report(reports, summary_report_path)
    print(f"   ✓ Sumário: {summary_report_path}")

    # 5. Gerar gráficos
    print("\n🎨 Gerando visualizações...")
    plot_onset_analysis(reports, output_root / "plots")
    print(f"   ✓ Gráficos salvos em: {output_root / 'plots'}")

    # 6. Análise de hiperanos
    print("\n🔬 Analisando onset de hiperanos...")
    hyperon_onsets = analyze_hyperon_onset(datasets_by_model_b)

    with open(output_root / "hyperon_onset.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Modelo_B", "log10(ξ) Onset Hiperão"])
        for key, onset_log_csi in sorted(hyperon_onsets.items()):
            if onset_log_csi is not None:
                writer.writerow([key, f"{onset_log_csi:.2f}"])
            else:
                writer.writerow([key, "N/A"])

    print(f"   ✓ Onsets de hiperão: {output_root / 'hyperon_onset.csv'}")

    print("\n" + "=" * 60)
    print("✅ Análise concluída com sucesso!")
    print(f"📁 Resultados em: {output_root}")


if __name__ == "__main__":
    main()
