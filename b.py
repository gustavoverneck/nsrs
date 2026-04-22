#!/usr/bin/env python3
"""
Análise dos dados sem csi (output/b) com comparação direta contra NLEM.

Entregas principais:
- consolidação de métricas (Mmax e R em Mmax) para output/b;
- gráficos vs log10(B) para o caso sem csi;
- gráficos comparativos maxM_vs_logcsi e maxR_vs_logcsi para vários B (NLEM);
- comparação com/sem csi via curvas vs log10(B) e tabelas de relação direta;
- extração de extremos (mín./máx.) do campo magnético para cada métrica.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


COL_NB = 0
COL_EPS = 1
COL_P = 2

MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0

LABEL_NO_CSI = "sem csi"
X_LABEL_LOG_B = r"$\log_{10}(B\,[G])$"
X_LABEL_LOG_CSI = r"$\log_{10}(\xi)$"
LABEL_LOG_CSI = X_LABEL_LOG_CSI
Y_LABEL_MAXM = r"$M_{max}\,[M_\odot]$"
Y_LABEL_MAXR = r"$R(M_{max})\,[km]$"
CSI_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]


@dataclass
class BDataset:
	model: str
	b_label: str
	b_value: float
	file_path: Path
	data: np.ndarray


@dataclass
class BSummaryRow:
	model: str
	b_label: str
	b_value: float
	log_b_value: float
	n_eos_points: int
	n_mr_points: int
	max_mass_msun: float
	radius_at_max_km: float
	central_nb_n0: float
	central_eps_mevfm3: float
	central_p_mevfm3: float
	eos_path: str


@dataclass
class NlemSummaryRow:
	model: str
	b_label: str
	b_value: float
	csi: float
	log_csi: float
	max_mass_msun: float
	radius_at_max_km: float


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Análise de output/b com comparação NLEM")
	parser.add_argument("--b-input-root", type=Path, default=Path("output/b"))
	parser.add_argument(
		"--nlem-summary-csv",
		type=Path,
		default=Path("results/nlem_log/summary_all_datasets.csv"),
		help="CSV consolidado do nlem_log.py",
	)
	parser.add_argument("--output-root", type=Path, default=Path("results/b"))
	parser.add_argument("--dpi", type=int, default=150)
	parser.add_argument(
		"--max-b-curves",
		type=int,
		default=10,
		help="Máximo de curvas de B em maxM/maxR vs logcsi",
	)
	return parser.parse_args()


def ensure_dir(path: Path) -> None:
	path.mkdir(parents=True, exist_ok=True)


def _safe_float(text: str) -> Optional[float]:
	try:
		return float(text)
	except Exception:
		return None


def load_eos(path: Path) -> Optional[np.ndarray]:
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

	return arr


def parse_b_metadata(path: Path) -> Optional[Tuple[str, str, float]]:
	rx = re.compile(r".*/output/b/(?P<model>GM\d+)/B_(?P<b>[^/]+)/default/eos\.dat$")
	m = rx.match(path.resolve().as_posix())
	if not m:
		return None

	model = m.group("model")
	b_label = m.group("b")
	b_value = _safe_float(b_label)
	if b_value is None:
		return None
	return model, b_label, b_value


def discover_b_datasets(input_root: Path) -> List[BDataset]:
	out: List[BDataset] = []
	for eos_file in input_root.rglob("eos.dat"):
		meta = parse_b_metadata(eos_file)
		if meta is None:
			continue
		arr = load_eos(eos_file)
		if arr is None:
			continue
		model, b_label, b_value = meta
		out.append(BDataset(model=model, b_label=b_label, b_value=b_value, file_path=eos_file, data=arr))

	out.sort(key=lambda d: (d.model, d.b_value))
	return out


def valid_mr_mask(mass_col: np.ndarray, radius_col: np.ndarray) -> np.ndarray:
	return (
		np.isfinite(mass_col)
		& np.isfinite(radius_col)
		& (mass_col > 0.0)
		& (radius_col > 0.0)
		& (mass_col <= MAX_VALID_MASS_MSUN)
		& (radius_col <= MAX_VALID_RADIUS_KM)
	)


def summarize_b_dataset(ds: BDataset) -> BSummaryRow:
	arr = ds.data
	mass_col = arr[:, -2]
	radius_col = arr[:, -1]
	mr_mask = valid_mr_mask(mass_col, radius_col)

	n_mr_points = int(np.count_nonzero(mr_mask))
	n_eos_points = int(arr.shape[0])

	if n_mr_points > 0:
		valid_idx = np.nonzero(mr_mask)[0]
		local = int(np.argmax(mass_col[mr_mask]))
		i_max = int(valid_idx[local])

		max_mass = float(arr[i_max, -2])
		radius_at_max = float(arr[i_max, -1])
		central_nb = float(arr[i_max, COL_NB]) if arr.shape[1] > COL_NB else math.nan
		central_eps = float(arr[i_max, COL_EPS]) if arr.shape[1] > COL_EPS else math.nan
		central_p = float(arr[i_max, COL_P]) if arr.shape[1] > COL_P else math.nan
	else:
		max_mass = math.nan
		radius_at_max = math.nan
		central_nb = math.nan
		central_eps = math.nan
		central_p = math.nan

	return BSummaryRow(
		model=ds.model,
		b_label=ds.b_label,
		b_value=ds.b_value,
		log_b_value=(math.log10(ds.b_value) if ds.b_value > 0 else math.nan),
		n_eos_points=n_eos_points,
		n_mr_points=n_mr_points,
		max_mass_msun=max_mass,
		radius_at_max_km=radius_at_max,
		central_nb_n0=central_nb,
		central_eps_mevfm3=central_eps,
		central_p_mevfm3=central_p,
		eos_path=ds.file_path.as_posix(),
	)


def write_b_summary_csv(rows: Sequence[BSummaryRow], out_csv: Path) -> None:
	ensure_dir(out_csv.parent)
	with out_csv.open("w", newline="", encoding="utf-8") as f:
		w = csv.writer(f)
		w.writerow(
			[
				"model",
				"b_label",
				"b_value",
				"log_b_value",
				"n_eos_points",
				"n_mr_points",
				"max_mass_msun",
				"radius_at_max_km",
				"central_nb_n0",
				"central_eps_mevfm3",
				"central_p_mevfm3",
				"eos_path",
			]
		)
		for r in rows:
			w.writerow(
				[
					r.model,
					r.b_label,
					r.b_value,
					r.log_b_value,
					r.n_eos_points,
					r.n_mr_points,
					r.max_mass_msun,
					r.radius_at_max_km,
					r.central_nb_n0,
					r.central_eps_mevfm3,
					r.central_p_mevfm3,
					r.eos_path,
				]
			)


def load_nlem_summary_rows(path: Path) -> List[NlemSummaryRow]:
	if not path.exists():
		return []

	out: List[NlemSummaryRow] = []
	with path.open("r", newline="", encoding="utf-8") as f:
		reader = csv.DictReader(f)
		for row in reader:
			model = row.get("model", "").strip()
			b_label = row.get("b_label", "").strip()
			b_value = _safe_float(row.get("b_value", ""))
			csi = _safe_float(row.get("csi", ""))
			log_csi = _safe_float(row.get("log_csi", ""))
			max_mass = _safe_float(row.get("max_mass_msun", ""))
			radius_at_max = _safe_float(row.get("radius_at_max_km", ""))

			if (
				not model
				or not b_label
				or b_value is None
				or csi is None
				or log_csi is None
				or max_mass is None
				or radius_at_max is None
			):
				continue

			out.append(
				NlemSummaryRow(
					model=model,
					b_label=b_label,
					b_value=b_value,
					csi=csi,
					log_csi=log_csi,
					max_mass_msun=max_mass,
					radius_at_max_km=radius_at_max,
				)
			)

	out.sort(key=lambda r: (r.model, r.b_value, r.log_csi))
	return out


def _subsample_sorted_by_value(values: Sequence[float], max_items: int) -> List[float]:
	vals = sorted(set(values))
	if len(vals) <= max_items:
		return vals
	idx = np.linspace(0, len(vals) - 1, max_items).round().astype(int)
	idx = np.unique(idx)
	return [vals[i] for i in idx]


def plot_b_trends(rows: Sequence[BSummaryRow], out_dir: Path, dpi: int) -> None:
	ensure_dir(out_dir)
	models = sorted({r.model for r in rows})

	for model in models:
		sub = [r for r in rows if r.model == model and np.isfinite(r.log_b_value)]
		if not sub:
			continue

		sub = sorted(sub, key=lambda x: x.log_b_value)
		x = np.array([r.log_b_value for r in sub], dtype=float)
		y_m = np.array([r.max_mass_msun for r in sub], dtype=float)
		y_r = np.array([r.radius_at_max_km for r in sub], dtype=float)

		plt.figure(figsize=(8, 6))
		plt.plot(x, y_m, marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(Y_LABEL_MAXM)
		plt.title(f"maxM vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"maxM_vs_logB_{model}.png", dpi=dpi)
		plt.close()

		plt.figure(figsize=(8, 6))
		plt.plot(x, y_r, marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(Y_LABEL_MAXR)
		plt.title(f"maxR vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"maxR_vs_logB_{model}.png", dpi=dpi)
		plt.close()


def _plot_compare_metric_vs_logb(
	model: str,
	x_no: np.ndarray,
	y_no: np.ndarray,
	chosen_rows: Sequence[NlemSummaryRow],
	csi_targets: Sequence[float],
	y_metric: str,
	y_label: str,
	title: str,
	out_file: Path,
	dpi: int,
) -> None:
	plt.figure(figsize=(9, 6))
	plt.plot(x_no, y_no, color="black", lw=2.0, label=LABEL_NO_CSI)

	for idx, csi_t in enumerate(csi_targets):
		arr = sorted([r for r in chosen_rows if abs(r.log_csi - csi_t) < 1e-9], key=lambda x: x.b_value)
		if not arr:
			continue
		x = np.log10(np.array([r.b_value for r in arr], dtype=float))
		y = np.array([getattr(r, y_metric) for r in arr], dtype=float)
		marker = CSI_MARKERS[idx % len(CSI_MARKERS)]
		plt.plot(
			x,
			y,
			marker=marker,
			ms=6,
			lw=1.2,
			mec="black",
			mew=0.6,
			label=fr"csi={arr[0].log_csi:.2f} ({marker})",
		)

	plt.xlabel(X_LABEL_LOG_B)
	plt.ylabel(y_label)
	plt.title(f"Comparação direta com/sem csi | {title} | {model}")
	plt.grid(alpha=0.25)
	plt.legend(fontsize=8)
	plt.tight_layout()
	plt.savefig(out_file, dpi=dpi)
	plt.close()


def _safe_delta(a: float, b: float) -> float:
	if np.isfinite(a) and np.isfinite(b):
		return a - b
	return math.nan


def _safe_ratio(a: float, b: float) -> float:
	if np.isfinite(a) and np.isfinite(b) and abs(b) > 1e-14:
		return a / b
	return math.nan


def plot_nlem_metric_vs_logcsi(
	nlem_rows: Sequence[NlemSummaryRow],
	metric: str,
	y_label: str,
	title_prefix: str,
	out_path: Path,
	model: str,
	max_b_curves: int,
	dpi: int,
) -> None:
	sub = [r for r in nlem_rows if r.model == model]
	if not sub:
		return

	b_values = _subsample_sorted_by_value([r.b_value for r in sub], max_b_curves)
	if not b_values:
		return

	plt.figure(figsize=(9, 6))
	plotted = False
	for b in b_values:
		arr = sorted([r for r in sub if r.b_value == b], key=lambda x: x.log_csi)
		if not arr:
			continue
		x = np.array([r.log_csi for r in arr], dtype=float)
		y = np.array([getattr(r, metric) for r in arr], dtype=float)
		mask = np.isfinite(x) & np.isfinite(y)
		if np.count_nonzero(mask) < 2:
			continue
		plt.plot(x[mask], y[mask], marker="o", ms=2.5, lw=1.0, label=f"B={arr[0].b_label} G")
		plotted = True

	if not plotted:
		plt.close()
		return

	plt.xlabel(X_LABEL_LOG_CSI)
	plt.ylabel(y_label)
	plt.title(f"{title_prefix} | {model}")
	plt.grid(alpha=0.25)
	plt.legend(fontsize=7, ncol=2)
	plt.tight_layout()
	ensure_dir(out_path.parent)
	plt.savefig(out_path, dpi=dpi)
	plt.close()


def _nearest_rows_for_targets(rows: Sequence[NlemSummaryRow], targets: Sequence[float]) -> List[NlemSummaryRow]:
	if not rows:
		return []
	rows_sorted = sorted(rows, key=lambda r: r.log_csi)
	out: List[NlemSummaryRow] = []
	for t in targets:
		pick = min(rows_sorted, key=lambda r, tt=t: abs(r.log_csi - tt))
		if pick not in out:
			out.append(pick)
	out.sort(key=lambda r: r.log_csi)
	return out


def _select_csi_targets(nlem_rows: Sequence[NlemSummaryRow], model: str) -> List[float]:
	vals = sorted({r.log_csi for r in nlem_rows if r.model == model and np.isfinite(r.log_csi)})
	if not vals:
		return []
	return [vals[0], vals[len(vals) // 2], vals[-1]]


def _interp_metric_no_csi(
	b_rows: Sequence[BSummaryRow],
	model: str,
	b_values_target: np.ndarray,
	metric: str,
) -> np.ndarray:
	sub = [r for r in b_rows if r.model == model and r.b_value > 0 and np.isfinite(getattr(r, metric))]
	if len(sub) < 2:
		return np.full_like(b_values_target, np.nan, dtype=float)

	sub = sorted(sub, key=lambda r: r.b_value)
	x = np.log10(np.array([r.b_value for r in sub], dtype=float))
	y = np.array([getattr(r, metric) for r in sub], dtype=float)
	x_t = np.log10(np.asarray(b_values_target, dtype=float))

	valid = np.isfinite(x) & np.isfinite(y)
	if np.count_nonzero(valid) < 2:
		return np.full_like(x_t, np.nan, dtype=float)
	return np.interp(x_t, x[valid], y[valid], left=np.nan, right=np.nan)


def _build_series_by_csi_target(
	groups_by_b: Dict[float, List[NlemSummaryRow]],
	csi_targets: Sequence[float],
) -> Dict[float, List[NlemSummaryRow]]:
	series_by_target: Dict[float, List[NlemSummaryRow]] = {target: [] for target in csi_targets}

	for b_value, rows in groups_by_b.items():
		selected = _nearest_rows_for_targets(rows, csi_targets)
		for row in selected:
			series_by_target.setdefault(row.log_csi, []).append(row)

	for target in list(series_by_target.keys()):
		series_by_target[target] = sorted(series_by_target[target], key=lambda r: r.b_value)

	return series_by_target


def _plot_delta_and_ratio_for_metric(
	model: str,
	b_rows: Sequence[BSummaryRow],
	series_by_target: Dict[float, List[NlemSummaryRow]],
	metric: str,
	delta_label: str,
	ratio_label: str,
	output_dir: Path,
	dpi: int,
) -> None:
	ensure_dir(output_dir)
	metric_name = metric.replace("_", "")

	for target_log_csi, rows in sorted(series_by_target.items(), key=lambda kv: kv[0]):
		if not rows:
			continue

		rows = sorted(rows, key=lambda r: r.b_value)
		b_values = np.array([r.b_value for r in rows], dtype=float)
		x = np.log10(b_values)
		y_nlem = np.array([getattr(r, metric) for r in rows], dtype=float)
		y_no = _interp_metric_no_csi(b_rows, model, b_values, metric)
		mask = np.isfinite(x) & np.isfinite(y_nlem) & np.isfinite(y_no)
		if np.count_nonzero(mask) < 2:
			continue

		delta = y_nlem[mask] - y_no[mask]
		ratio = np.divide(y_nlem[mask], y_no[mask], out=np.full_like(y_nlem[mask], np.nan), where=np.abs(y_no[mask]) > 1e-14)
		xm = x[mask]

		fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
		marker = CSI_MARKERS[int(abs(target_log_csi * 10)) % len(CSI_MARKERS)]

		ax1.axhline(0.0, color="black", linestyle="--", alpha=0.7)
		ax1.plot(xm, delta, marker=marker, lw=1.4, ms=5, label=fr"log10(csi)={target_log_csi:.2f}")
		ax1.set_ylabel(delta_label)
		ax1.set_title(f"Delta comparison | {model} | {metric_name}")
		ax1.grid(alpha=0.25)
		ax1.legend(fontsize=8)

		ax2.axhline(1.0, color="black", linestyle="--", alpha=0.7)
		ax2.plot(xm, ratio, marker=marker, lw=1.4, ms=5, label=fr"log10(csi)={target_log_csi:.2f}")
		ax2.set_xlabel(X_LABEL_LOG_B)
		ax2.set_ylabel(ratio_label)
		ax2.grid(alpha=0.25)
		ax2.legend(fontsize=8)

		fig.tight_layout()
		fig.savefig(output_dir / f"delta_ratio_{metric_name}_{model}_csi_{target_log_csi:.2f}.png", dpi=dpi)
		plt.close(fig)


def plot_direct_comparison_diagnostics(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_dir: Path,
	dpi: int,
) -> None:
	ensure_dir(out_dir)
	models = sorted({r.model for r in nlem_rows} | {r.model for r in b_rows})

	for model in models:
		no_csi = sorted(
			[r for r in b_rows if r.model == model and r.b_value > 0 and np.isfinite(r.log_b_value)],
			key=lambda x: x.log_b_value,
		)
		nlem_model = [r for r in nlem_rows if r.model == model and r.b_value > 0]
		if not nlem_model or not no_csi:
			continue

		csi_targets = _select_csi_targets(nlem_model, model)
		if not csi_targets:
			continue

		groups_by_b: Dict[float, List[NlemSummaryRow]] = {}
		for r in nlem_model:
			groups_by_b.setdefault(r.b_value, []).append(r)

		series_by_target = _build_series_by_csi_target(groups_by_b, csi_targets)
		model_dir = out_dir / model

		_plot_delta_and_ratio_for_metric(
			model=model,
			b_rows=b_rows,
			series_by_target=series_by_target,
			metric="max_mass_msun",
			delta_label=r"$\Delta M_{max}$ [M$_\odot$]",
			ratio_label=r"$M_{max}^{csi} / M_{max}^{sem}$",
			output_dir=model_dir,
			dpi=dpi,
		)
		_plot_delta_and_ratio_for_metric(
			model=model,
			b_rows=b_rows,
			series_by_target=series_by_target,
			metric="radius_at_max_km",
			delta_label=r"$\Delta R(M_{max})$ [km]",
			ratio_label=r"$R^{csi}(M_{max}) / R^{sem}(M_{max})$",
			output_dir=model_dir,
			dpi=dpi,
		)


def plot_compare_with_without_csi(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_dir: Path,
	dpi: int,
) -> None:
	ensure_dir(out_dir)
	models = sorted({r.model for r in nlem_rows} | {r.model for r in b_rows})

	for model in models:
		no_csi = sorted(
			[r for r in b_rows if r.model == model and r.b_value > 0 and np.isfinite(r.log_b_value)],
			key=lambda x: x.log_b_value,
		)
		nlem_model = [r for r in nlem_rows if r.model == model and r.b_value > 0]
		if not nlem_model or not no_csi:
			continue

		csi_targets = _select_csi_targets(nlem_model, model)
		if not csi_targets:
			continue

		groups_by_b: Dict[float, List[NlemSummaryRow]] = {}
		for r in nlem_model:
			groups_by_b.setdefault(r.b_value, []).append(r)

		chosen_rows: List[NlemSummaryRow] = []
		for arr in groups_by_b.values():
			chosen_rows.extend(_nearest_rows_for_targets(arr, csi_targets))

		chosen_rows.sort(key=lambda r: (r.log_csi, r.b_value))

		x_no = np.array([r.log_b_value for r in no_csi], dtype=float)
		y_no_m = np.array([r.max_mass_msun for r in no_csi], dtype=float)
		y_no_r = np.array([r.radius_at_max_km for r in no_csi], dtype=float)

		_plot_compare_metric_vs_logb(
			model=model,
			x_no=x_no,
			y_no=y_no_m,
			chosen_rows=chosen_rows,
			csi_targets=csi_targets,
			y_metric="max_mass_msun",
			y_label=Y_LABEL_MAXM,
			title="maxM",
			out_file=out_dir / f"compare_with_without_csi_maxM_vs_logB_{model}.png",
			dpi=dpi,
		)

		_plot_compare_metric_vs_logb(
			model=model,
			x_no=x_no,
			y_no=y_no_r,
			chosen_rows=chosen_rows,
			csi_targets=csi_targets,
			y_metric="radius_at_max_km",
			y_label=Y_LABEL_MAXR,
			title="maxR",
			out_file=out_dir / f"compare_with_without_csi_maxR_vs_logB_{model}.png",
			dpi=dpi,
		)


def _relation_rows_for_model(
	model: str,
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
) -> List[List[float | str]]:
	sub_n = [r for r in nlem_rows if r.model == model and r.b_value > 0]
	if not sub_n:
		return []

	b_targets = np.array([r.b_value for r in sub_n], dtype=float)
	m_no = _interp_metric_no_csi(b_rows, model, b_targets, "max_mass_msun")
	r_no = _interp_metric_no_csi(b_rows, model, b_targets, "radius_at_max_km")

	out: List[List[float | str]] = []
	for i, rr in enumerate(sub_n):
		m_n = rr.max_mass_msun
		r_n = rr.radius_at_max_km
		m0 = float(m_no[i])
		r0 = float(r_no[i])

		dm = _safe_delta(m_n, m0)
		dr = _safe_delta(r_n, r0)
		rm = _safe_ratio(m_n, m0)
		rr_ = _safe_ratio(r_n, r0)

		out.append(
			[
				rr.model,
				rr.b_value,
				(math.log10(rr.b_value) if rr.b_value > 0 else math.nan),
				rr.csi,
				rr.log_csi,
				m_n,
				m0,
				dm,
				rm,
				r_n,
				r0,
				dr,
				rr_,
			]
		)

	return out


def write_direct_relation_csv(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_csv: Path,
) -> None:
	ensure_dir(out_csv.parent)
	models = sorted({r.model for r in nlem_rows} & {r.model for r in b_rows})

	with out_csv.open("w", newline="", encoding="utf-8") as f:
		w = csv.writer(f)
		w.writerow(
			[
				"model",
				"b_value",
				"log_b_value",
				"csi",
				"log_csi",
				"max_mass_nlem",
				"max_mass_no_csi_interp",
				"delta_max_mass",
				"ratio_max_mass",
				"radius_nlem",
				"radius_no_csi_interp",
				"delta_radius",
				"ratio_radius",
			]
		)

		for model in models:
			for row in _relation_rows_for_model(model, b_rows, nlem_rows):
				w.writerow(row)


def _extrema_row(items: Sequence[Tuple[float, float]]) -> Tuple[Tuple[float, float], Tuple[float, float]]:
	valid = [(x, y) for x, y in items if np.isfinite(x) and np.isfinite(y)]
	if not valid:
		return (math.nan, math.nan), (math.nan, math.nan)
	return min(valid, key=lambda p: p[1]), max(valid, key=lambda p: p[1])


def write_extrema_report(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_md: Path,
) -> None:
	ensure_dir(out_md.parent)
	models = sorted({r.model for r in b_rows} | {r.model for r in nlem_rows})

	lines: List[str] = []
	lines.append("# Extremos de campo magnético (com e sem csi)")
	lines.append("")
	lines.append("Relações extraídas para $M_{max}$ e $R(M_{max})$.")
	lines.append("")

	for model in models:
		lines.append(f"## {model}")

		no_csi = [r for r in b_rows if r.model == model and r.b_value > 0]
		with_csi = [r for r in nlem_rows if r.model == model and r.b_value > 0]

		mn_no, mx_no = _extrema_row([(r.b_value, r.max_mass_msun) for r in no_csi])
		rn_no, rx_no = _extrema_row([(r.b_value, r.radius_at_max_km) for r in no_csi])

		mn_csi, mx_csi = _extrema_row([(r.b_value, r.max_mass_msun) for r in with_csi])
		rn_csi, rx_csi = _extrema_row([(r.b_value, r.radius_at_max_km) for r in with_csi])

		lines.append("### Sem csi")
		lines.append(f"- min $M_{{max}}$: B={mn_no[0]:.4e} G, $M_{{max}}$={mn_no[1]:.6f}")
		lines.append(f"- max $M_{{max}}$: B={mx_no[0]:.4e} G, $M_{{max}}$={mx_no[1]:.6f}")
		lines.append(f"- min $R(M_{{max}})$: B={rn_no[0]:.4e} G, R={rn_no[1]:.6f} km")
		lines.append(f"- max $R(M_{{max}})$: B={rx_no[0]:.4e} G, R={rx_no[1]:.6f} km")
		lines.append("")

		lines.append("### Com csi (NLEM)")
		lines.append(f"- min $M_{{max}}$: B={mn_csi[0]:.4e} G, $M_{{max}}$={mn_csi[1]:.6f}")
		lines.append(f"- max $M_{{max}}$: B={mx_csi[0]:.4e} G, $M_{{max}}$={mx_csi[1]:.6f}")
		lines.append(f"- min $R(M_{{max}})$: B={rn_csi[0]:.4e} G, R={rn_csi[1]:.6f} km")
		lines.append(f"- max $R(M_{{max}})$: B={rx_csi[0]:.4e} G, R={rx_csi[1]:.6f} km")
		lines.append("")

		# Relação direta aproximada via interpolação (média dos deltas por csi)
		if no_csi and with_csi:
			b_targets = np.array([r.b_value for r in with_csi], dtype=float)
			m_no_interp = _interp_metric_no_csi(b_rows, model, b_targets, "max_mass_msun")
			r_no_interp = _interp_metric_no_csi(b_rows, model, b_targets, "radius_at_max_km")
			dmass = np.array([r.max_mass_msun for r in with_csi], dtype=float) - m_no_interp
			drad = np.array([r.radius_at_max_km for r in with_csi], dtype=float) - r_no_interp

			m_mask = np.isfinite(dmass)
			r_mask = np.isfinite(drad)
			dmass_mean = float(np.nanmean(dmass[m_mask])) if np.any(m_mask) else math.nan
			drad_mean = float(np.nanmean(drad[r_mask])) if np.any(r_mask) else math.nan

			lines.append("### Relação direta média (com csi - sem csi)")
			lines.append(f"- $\\Delta M_{{max}}$ médio: {dmass_mean:.6e}")
			lines.append(f"- $\\Delta R(M_{{max}})$ médio: {drad_mean:.6e} km")
			lines.append("")

	out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
	args = parse_args()
	ensure_dir(args.output_root)

	b_datasets = discover_b_datasets(args.b_input_root)
	if not b_datasets:
		print(f"[ERRO] Nenhum eos.dat válido em {args.b_input_root}")
		return
	print(f"[INFO] datasets sem csi: {len(b_datasets)}")

	b_rows = [summarize_b_dataset(ds) for ds in b_datasets]
	b_summary_csv = args.output_root / "summary_b_no_csi.csv"
	write_b_summary_csv(b_rows, b_summary_csv)
	print(f"[OK] resumo sem csi: {b_summary_csv}")

	plot_b_trends(b_rows, args.output_root / "trends", dpi=args.dpi)
	print(f"[OK] tendências sem csi: {args.output_root / 'trends'}")

	nlem_rows = load_nlem_summary_rows(args.nlem_summary_csv)
	if not nlem_rows:
		print(f"[AVISO] CSV NLEM não encontrado/válido: {args.nlem_summary_csv}")
		return
	print(f"[INFO] linhas NLEM carregadas: {len(nlem_rows)}")

	comp_dir = args.output_root / "comparison"
	ensure_dir(comp_dir)

	models = sorted({r.model for r in nlem_rows})
	for model in models:
		plot_nlem_metric_vs_logcsi(
			nlem_rows=nlem_rows,
			metric="max_mass_msun",
			y_label=r"$M_{max}\,[M_\odot]$",
			title_prefix="maxM vs log10(csi) para vários B",
			out_path=comp_dir / f"maxM_vs_logcsi_{model}.png",
			model=model,
			max_b_curves=max(1, args.max_b_curves),
			dpi=args.dpi,
		)
		plot_nlem_metric_vs_logcsi(
			nlem_rows=nlem_rows,
			metric="radius_at_max_km",
			y_label=r"$R(M_{max})\,[km]$",
			title_prefix="maxR vs log10(csi) para vários B",
			out_path=comp_dir / f"maxR_vs_logcsi_{model}.png",
			model=model,
			max_b_curves=max(1, args.max_b_curves),
			dpi=args.dpi,
		)

	plot_compare_with_without_csi(
		b_rows=b_rows,
		nlem_rows=nlem_rows,
		out_dir=comp_dir,
		dpi=args.dpi,
	)
	print(f"[OK] plots comparativos com/sem csi: {comp_dir}")

	diagnostics_dir = comp_dir / "diagnostics"
	plot_direct_comparison_diagnostics(
		b_rows=b_rows,
		nlem_rows=nlem_rows,
		out_dir=diagnostics_dir,
		dpi=args.dpi,
	)
	print(f"[OK] gráficos diagnósticos com/sem csi: {diagnostics_dir}")

	relation_csv = comp_dir / "relation_with_without_csi.csv"
	write_direct_relation_csv(b_rows, nlem_rows, relation_csv)
	print(f"[OK] relação direta em CSV: {relation_csv}")

	extrema_md = args.output_root / "extrema_report.md"
	write_extrema_report(b_rows, nlem_rows, extrema_md)
	print(f"[OK] relatório de extremos: {extrema_md}")

	print("[DONE] Análise de B finalizada.")


if __name__ == "__main__":
	main()
