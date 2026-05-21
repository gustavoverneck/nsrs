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
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

try:
	import plotly.graph_objects as go  # type: ignore[import-not-found]
	PLOTLY_AVAILABLE = True
except ImportError:
	go = None  # type: ignore[assignment]
	PLOTLY_AVAILABLE = False


COL_NB = 0
COL_EPS = 1
COL_P = 2
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
COL_VSIGMA = 13
COL_VOMEGA = 14
COL_VRHO = 15
COL_MEFF = 16
COL_MUN = 17
COL_MUE = 18
COL_EMAG = 19
COL_BDD = 20

MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0

LABEL_NO_CSI = "sem csi"
X_LABEL_LOG_B = r"$\log_{10}(B\,[G])$"
X_LABEL_LOG_CSI = r"$\log_{10}(\xi)$"
LABEL_LOG_CSI = X_LABEL_LOG_CSI
Y_LABEL_MAXM = r"$M_{max}\,[M_\odot]$"
Y_LABEL_MAXR = r"$R(M_{max})\,[km]$"
Y_LABEL_MAXM_SHORT = r"$M_{max}$"
Z_LABEL_MAXR_SHORT = r"$R(M_{max})$"
TARGET_CANONICAL_MASS_MSUN = 1.4
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



def plot_all_variables_vs_logb(rows: Sequence[BSummaryRow], out_dir: Path, dpi: int) -> None:
	"""Plot all summary variables vs log10(B)."""
	ensure_dir(out_dir)
	models = sorted({r.model for r in rows})
	metrics = [
		("max_mass_msun", Y_LABEL_MAXM, "max_mass_msun"),
		("radius_at_max_km", Y_LABEL_MAXR, "radius_at_max_km"),
		("central_nb_n0", r"$n_{B,\,c}$ [n0]", "central_nb_n0"),
		("central_eps_mevfm3", r"$\epsilon_c$ [MeV/fm$^3$]", "central_eps_mevfm3"),
		("central_p_mevfm3", r"$P_c$ [MeV/fm$^3$]", "central_p_mevfm3"),
		("n_eos_points", r"$N_{\mathrm{EOS}}$", "n_eos_points"),
		("n_mr_points", r"$N_{\mathrm{MR}}$", "n_mr_points"),
	]

	for model in models:
		sub = [r for r in rows if r.model == model and np.isfinite(r.log_b_value)]
		if not sub:
			continue
		sub = sorted(sub, key=lambda x: x.log_b_value)
		x = np.array([r.log_b_value for r in sub], dtype=float)

		for attr, y_label, suffix in metrics:
			y = np.array([getattr(r, attr) for r in sub], dtype=float)
			mask = np.isfinite(x) & np.isfinite(y)
			if np.count_nonzero(mask) < 2:
				continue

			plt.figure(figsize=(8, 6))
			plt.plot(x[mask], y[mask], marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
			plt.xlabel(X_LABEL_LOG_B)
			plt.ylabel(y_label)
			plt.title(f"{suffix} vs log10(B) | {model} | sem csi")
			plt.grid(alpha=0.25)
			plt.legend(fontsize=9)
			plt.tight_layout()
			plt.savefig(out_dir / f"{suffix}_vs_logB_{model}.png", dpi=dpi)
			plt.close()


def plot_all_core_variables_vs_logb(
	b_datasets: Sequence[BDataset],
	out_dir: Path,
	dpi: int,
) -> None:  # noqa: C901
	"""Plot all core variables (columns 0..20) at max-mass point vs log10(B)."""
	ensure_dir(out_dir)
	models = sorted({d.model for d in b_datasets})
	variables = [
		(COL_NB, r"$n/n_0$", "nb_over_n0"),
		(COL_EPS, r"$\epsilon_{tot}$ [MeV/fm$^3$]", "eps_total"),
		(COL_P, r"$P_{tot}$ [MeV/fm$^3$]", "p_total"),
		(COL_NE, r"$n_e$ [fm$^{-3}$]", "ne"),
		(COL_NMU, r"$n_\mu$ [fm$^{-3}$]", "nmu"),
		(COL_NN, r"$n_n$ [fm$^{-3}$]", "nn"),
		(COL_NP, r"$n_p$ [fm$^{-3}$]", "np"),
		(COL_NL0, r"$n_{\Lambda^0}$ [fm$^{-3}$]", "nL0"),
		(COL_NSM, r"$n_{\Sigma^-}$ [fm$^{-3}$]", "nSm"),
		(COL_NS0, r"$n_{\Sigma^0}$ [fm$^{-3}$]", "nS0"),
		(COL_NSP, r"$n_{\Sigma^+}$ [fm$^{-3}$]", "nSp"),
		(COL_NXM, r"$n_{\Xi^-}$ [fm$^{-3}$]", "nXm"),
		(COL_NX0, r"$n_{\Xi^0}$ [fm$^{-3}$]", "nX0"),
		(COL_VSIGMA, r"$\sigma$ [MeV]", "sigma"),
		(COL_VOMEGA, r"$\omega$ [MeV]", "omega"),
		(COL_VRHO, r"$\rho$ [MeV]", "rho"),
		(COL_MEFF, r"$m^*/m_N$", "mstar_over_mn"),
		(COL_MUN, r"$\mu_n$", "mu_n"),
		(COL_MUE, r"$\mu_e$", "mu_e"),
		(COL_EMAG, r"$\epsilon_B$ [MeV/fm$^3$]", "emag"),
		(COL_BDD, r"$B(n)$ [T]", "b_local"),
	]

	for model in models:
		sub = [d for d in b_datasets if d.model == model and d.b_value > 0]
		if not sub:
			continue
		sub = sorted(sub, key=lambda d: d.b_value)
		x = np.array([math.log10(d.b_value) for d in sub], dtype=float)

		for col, y_label, suffix in variables:
			y_vals: List[float] = []
			for d in sub:
				arr = d.data
				if arr.shape[1] <= col:
					y_vals.append(math.nan)
					continue
				mass_col = arr[:, -2]
				radius_col = arr[:, -1]
				mr_mask = valid_mr_mask(mass_col, radius_col)
				if not np.any(mr_mask):
					y_vals.append(math.nan)
					continue
				valid_idx = np.nonzero(mr_mask)[0]
				local = int(np.argmax(mass_col[mr_mask]))
				i_max = int(valid_idx[local])
				y_vals.append(float(arr[i_max, col]))

			y = np.array(y_vals, dtype=float)
			mask = np.isfinite(x) & np.isfinite(y)
			if np.count_nonzero(mask) < 2:
				continue

			plt.figure(figsize=(8, 6))
			plt.plot(x[mask], y[mask], marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
			plt.xlabel(X_LABEL_LOG_B)
			plt.ylabel(y_label)
			plt.title(f"{suffix} vs log10(B) | {model} | sem csi")
			plt.grid(alpha=0.25)
			plt.legend(fontsize=9)
			plt.tight_layout()
			plt.savefig(out_dir / f"{suffix}_vs_logB_{model}.png", dpi=dpi)
			plt.close()


def plot_b_trends(rows: Sequence[BSummaryRow], b_datasets: Sequence[BDataset], out_dir: Path, dpi: int) -> None:  # noqa: C901
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

		# Plotar M_max vs log_B
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

		# Plotar R_max vs log_B
		plt.figure(figsize=(8, 6))
		plt.plot(x, y_r, marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(Y_LABEL_MAXR)
		plt.title(f"R(M_max) vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"R_max_vs_logB_{model}.png", dpi=dpi)
		plt.close()
		
		# Calcular e plotar taxas de crescimento de R_max
		_plot_growth_rates(x, y_m, y_r, model, out_dir, dpi)
		
		# Extrair raio em massa constante (1.4 M_solar)
		target_mass = TARGET_CANONICAL_MASS_MSUN
		y_r_at_const_m = _extract_radius_at_fixed_mass(model, sub, b_datasets, target_mass)
		
		if len(y_r_at_const_m) > 0:
			valid_r_mask = np.isfinite(y_r_at_const_m)
			if np.count_nonzero(valid_r_mask) > 1:
				# Plotar R(M=const) vs log_B
				plt.figure(figsize=(8, 6))
				plt.plot(x[valid_r_mask], np.array(y_r_at_const_m)[valid_r_mask], marker="o", ms=3, lw=1.2, label=LABEL_NO_CSI)
				plt.xlabel(X_LABEL_LOG_B)
				plt.ylabel(Y_LABEL_MAXR)
				plt.title(f"R(M={target_mass:.1f} M☉) vs log10(B) | {model} | sem csi")
				plt.grid(alpha=0.25)
				plt.legend(fontsize=9)
				plt.tight_layout()
				plt.savefig(out_dir / f"R_at_const_M_vs_logB_{model}.png", dpi=dpi)
				plt.close()
				
				# Calcular e plotar taxas de crescimento de R(M=const)
				_plot_growth_rates_const_mass(x, np.array(y_r_at_const_m), model, target_mass, out_dir, dpi)


def _extract_radius_at_fixed_mass(
	model: str, summary_rows: Sequence[BSummaryRow], b_datasets: Sequence[BDataset], target_mass: float
) -> List[float]:
	"""Extrai o raio para uma massa fixa interpolando os dados M-R de cada dataset."""
	radii = []
	
	for s_row in summary_rows:
		# Encontrar dataset correspondente
		matching_ds = [
			d for d in b_datasets 
			if d.model == model and abs(d.b_value - s_row.b_value) < 1e-6
		]
		
		if not matching_ds:
			radii.append(np.nan)
			continue
		
		ds = matching_ds[0]
		mass_col = ds.data[:, -2]
		radius_col = ds.data[:, -1]
		
		# Filtrar dados válidos
		valid_mask = valid_mr_mask(mass_col, radius_col)
		if np.count_nonzero(valid_mask) < 2:
			radii.append(np.nan)
			continue
		
		mass_valid = mass_col[valid_mask]
		radius_valid = radius_col[valid_mask]
		
		# Ordenar por massa
		sorted_idx = np.argsort(mass_valid)
		mass_sorted = mass_valid[sorted_idx]
		radius_sorted = radius_valid[sorted_idx]
		
		# Interpolar para encontrar raio em massa alvo
		if mass_sorted.min() <= target_mass <= mass_sorted.max():
			r_at_target = float(np.interp(target_mass, mass_sorted, radius_sorted))
			radii.append(r_at_target)
		else:
			radii.append(np.nan)
	
	return radii


def _plot_growth_rates_const_mass(
	x: np.ndarray, y_r: np.ndarray, model: str, target_mass: float, out_dir: Path, dpi: int
) -> None:
	"""Plota a taxa de crescimento de R em massa constante."""
	valid_r = np.isfinite(y_r)
	
	if np.count_nonzero(valid_r) > 2:
		x_valid = x[valid_r]
		y_r_valid = y_r[valid_r]
		# Calcular derivada numérica usando diferenças finitas
		dx = np.diff(x_valid)
		dy_r = np.diff(y_r_valid)
		# Usar ponto médio para x
		x_mid_r = x_valid[:-1] + dx / 2.0
		growth_r = dy_r / np.where(np.abs(dx) > 1e-10, dx, 1.0)
		
		# Plotar taxa de crescimento de R(M=const)
		plt.figure(figsize=(8, 6))
		plt.plot(x_mid_r, growth_r, marker="^", ms=4, lw=1.2, color="darkred", label="dR/d(log B)")
		plt.axhline(0.0, color="black", linestyle="--", alpha=0.5)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(r"$\frac{dR(M=const)}{d\log_{10}(B)}$ [km]", fontsize=11)
		plt.title(f"Taxa de crescimento de R(M={target_mass:.1f} M☉) vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"growth_rate_R_const_M_vs_logB_{model}.png", dpi=dpi)
		plt.close()



def _plot_growth_rates(
	x: np.ndarray, y_m: np.ndarray, y_r: np.ndarray, model: str, out_dir: Path, dpi: int
) -> None:
	"""Plota as taxas de crescimento (derivadas) de M_max e R_max em relação a log_B."""
	# Máscaras de valores finitos
	valid_m = np.isfinite(y_m)
	valid_r = np.isfinite(y_r)
	
	if np.count_nonzero(valid_m) > 2:
		x_valid = x[valid_m]
		y_m_valid = y_m[valid_m]
		# Calcular derivada numérica usando diferenças finitas
		dx = np.diff(x_valid)
		dy_m = np.diff(y_m_valid)
		# Usar ponto médio para x
		x_mid_m = x_valid[:-1] + dx / 2.0
		growth_m = dy_m / np.where(np.abs(dx) > 1e-10, dx, 1.0)
		
		# Plotar taxa de crescimento de M_max
		plt.figure(figsize=(8, 6))
		plt.plot(x_mid_m, growth_m, marker="s", ms=4, lw=1.2, color="darkblue", label="dM/d(log B)")
		plt.axhline(0.0, color="black", linestyle="--", alpha=0.5)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(r"$\frac{dM_{max}}{d\log_{10}(B)}$ [$M_\odot$]", fontsize=11)
		plt.title(f"Taxa de crescimento de maxM vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"growth_rate_maxM_vs_logB_{model}.png", dpi=dpi)
		plt.close()
	
	if np.count_nonzero(valid_r) > 2:
		x_valid = x[valid_r]
		y_r_valid = y_r[valid_r]
		# Calcular derivada numérica usando diferenças finitas
		dx = np.diff(x_valid)
		dy_r = np.diff(y_r_valid)
		# Usar ponto médio para x
		x_mid_r = x_valid[:-1] + dx / 2.0
		growth_r = dy_r / np.where(np.abs(dx) > 1e-10, dx, 1.0)
		
		# Plotar taxa de crescimento de R_max
		plt.figure(figsize=(8, 6))
		plt.plot(x_mid_r, growth_r, marker="^", ms=4, lw=1.2, color="darkred", label="dR/d(log B)")
		plt.axhline(0.0, color="black", linestyle="--", alpha=0.5)
		plt.xlabel(X_LABEL_LOG_B)
		plt.ylabel(r"$\frac{dR(M_{max})}{d\log_{10}(B)}$ [km]", fontsize=11)
		plt.title(f"Taxa de crescimento de R(M_max) vs log10(B) | {model} | sem csi")
		plt.grid(alpha=0.25)
		plt.legend(fontsize=9)
		plt.tight_layout()
		plt.savefig(out_dir / f"growth_rate_R_max_vs_logB_{model}.png", dpi=dpi)
		plt.close()
		
	# Plotar ambas as taxas de crescimento juntas (eixos normalizados)
	if np.count_nonzero(valid_m) > 2 and np.count_nonzero(valid_r) > 2:
		_plot_combined_growth_rates(x, y_m, y_r, model, out_dir, dpi)



def _plot_combined_growth_rates(
	x: np.ndarray, y_m: np.ndarray, y_r: np.ndarray, model: str, out_dir: Path, dpi: int
) -> None:
	"""Plota taxas de crescimento de M_max e R_max combinadas."""
	valid_m = np.isfinite(y_m)
	valid_r = np.isfinite(y_r)
	
	x_valid_m = x[valid_m]
	y_m_valid = y_m[valid_m]
	dx_m = np.diff(x_valid_m)
	dy_m = np.diff(y_m_valid)
	x_mid_m = x_valid_m[:-1] + dx_m / 2.0
	growth_m = dy_m / np.where(np.abs(dx_m) > 1e-10, dx_m, 1.0)
	
	x_valid_r = x[valid_r]
	y_r_valid = y_r[valid_r]
	dx_r = np.diff(x_valid_r)
	dy_r = np.diff(y_r_valid)
	x_mid_r = x_valid_r[:-1] + dx_r / 2.0
	growth_r = dy_r / np.where(np.abs(dx_r) > 1e-10, dx_r, 1.0)
	
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=False)
	
	ax1.plot(x_mid_m, growth_m, marker="s", ms=4, lw=1.2, color="darkblue")
	ax1.axhline(0.0, color="black", linestyle="--", alpha=0.5)
	ax1.set_ylabel(r"$\frac{dM_{max}}{d\log_{10}(B)}$ [$M_\odot$]", fontsize=10)
	ax1.set_title(f"Taxas de crescimento vs log10(B) | {model} | sem csi", fontsize=11)
	ax1.grid(alpha=0.25)
	
	ax2.plot(x_mid_r, growth_r, marker="^", ms=4, lw=1.2, color="darkred")
	ax2.axhline(0.0, color="black", linestyle="--", alpha=0.5)
	ax2.set_xlabel(X_LABEL_LOG_B)
	ax2.set_ylabel(r"$\frac{dR(M_{max})}{d\log_{10}(B)}$ [km]", fontsize=10)
	ax2.grid(alpha=0.25)
	
	fig.tight_layout()
	fig.savefig(out_dir / f"growth_rates_combined_max_{model}.png", dpi=dpi)
	plt.close(fig)




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


def _report_gap_between(
	items: Sequence[float],
	start: float,
	end: float,
) -> Tuple[List[float], bool]:
	"""Return (values_in_range, has_gap) for (start, end)."""
	vals = sorted(v for v in items if np.isfinite(v))
	in_range = [v for v in vals if start < v < end]
	has_gap = len(in_range) == 0
	return in_range, has_gap


def _append_gap_section(
	lines: List[str],
	no_csi: Sequence[BSummaryRow],
	with_csi: Sequence[NlemSummaryRow],
) -> None:
	b_vals_no = [r.b_value for r in no_csi if r.b_value > 0]
	b_vals_nlem = [r.b_value for r in with_csi if r.b_value > 0]
	in_range_no, gap_no = _report_gap_between(b_vals_no, 1e16, 1e17)
	in_range_nlem, gap_nlem = _report_gap_between(b_vals_nlem, 1e16, 1e17)

	lines.append("### Cobertura de B entre 1e16 e 1e17")
	lines.append(f"- Sem csi: {len(in_range_no)} pontos no intervalo (1e16, 1e17)")
	if in_range_no:
		lines.append(f"  - valores: {', '.join(f'{v:.2e}' for v in in_range_no)}")
	if gap_no:
		lines.append("  - gap identificado (nenhum valor intermediário)")
	lines.append(f"- Com csi: {len(in_range_nlem)} pontos no intervalo (1e16, 1e17)")
	if in_range_nlem:
		lines.append(f"  - valores: {', '.join(f'{v:.2e}' for v in in_range_nlem)}")
	if gap_nlem:
		lines.append("  - gap identificado (nenhum valor intermediário)")
	lines.append("")


def plot_3d_surface_m_r_vs_b(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_dir: Path,
	dpi: int,
) -> None:
	"""Gera plots 3D de superfície: Massa vs Raio vs Campo Magnético"""
	ensure_dir(out_dir)
	models = sorted({r.model for r in nlem_rows} | {r.model for r in b_rows})

	for model in models:
		# ===== Sem CSI =====
		no_csi = [r for r in b_rows if r.model == model and r.b_value > 0 and 
		         np.isfinite(r.log_b_value) and np.isfinite(r.max_mass_msun) and 
		         np.isfinite(r.radius_at_max_km)]
		
		if no_csi:
			no_csi = sorted(no_csi, key=lambda x: x.b_value)
			b_vals_no = np.array([r.b_value for r in no_csi], dtype=float)
			m_vals_no = np.array([r.max_mass_msun for r in no_csi], dtype=float)
			r_vals_no = np.array([r.radius_at_max_km for r in no_csi], dtype=float)
			
			fig = plt.figure(figsize=(12, 9))
			ax = fig.add_subplot(111, projection='3d')
			
			scatter = ax.scatter(xs=b_vals_no, ys=m_vals_no, zs=r_vals_no, c=np.log10(b_vals_no),  # type: ignore
			           cmap='viridis', s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
			
			ax.set_xlabel(r'$B\,[G]$', fontsize=11)
			ax.set_ylabel(Y_LABEL_MAXM, fontsize=11)
			ax.set_zlabel(Y_LABEL_MAXR, fontsize=11)
			ax.set_title(f'{Y_LABEL_MAXM_SHORT} vs {Z_LABEL_MAXR_SHORT} vs B | {model} | sem csi', fontsize=12)
			
			cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, pad=0.1)
			cbar.set_label(r'$\log_{10}(B)$ [log10(G)]', fontsize=10)
			
			fig.tight_layout()
			fig.savefig(out_dir / f"surface_3d_no_csi_{model}.png", dpi=dpi)
			plt.close(fig)
		
		# ===== Com CSI (NLEM) =====
		with_csi = [r for r in nlem_rows if r.model == model and r.b_value > 0 and
		           np.isfinite(r.max_mass_msun) and np.isfinite(r.radius_at_max_km)]
		
		if with_csi:
			with_csi = sorted(with_csi, key=lambda x: (x.b_value, x.log_csi))
			b_vals = np.array([r.b_value for r in with_csi], dtype=float)
			m_vals = np.array([r.max_mass_msun for r in with_csi], dtype=float)
			r_vals = np.array([r.radius_at_max_km for r in with_csi], dtype=float)
			
			fig = plt.figure(figsize=(12, 9))
			ax = fig.add_subplot(111, projection='3d')
			
			scatter = ax.scatter(xs=b_vals, ys=m_vals, zs=r_vals, c=np.array([r.log_csi for r in with_csi]),  # type: ignore
			           cmap='plasma', s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
			
			ax.set_xlabel(r'$B\,[G]$', fontsize=11)
			ax.set_ylabel(Y_LABEL_MAXM, fontsize=11)
			ax.set_zlabel(Y_LABEL_MAXR, fontsize=11)
			ax.set_title(f'{Y_LABEL_MAXM_SHORT} vs {Z_LABEL_MAXR_SHORT} vs B | {model} | com csi', fontsize=12)
			
			cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, pad=0.1)
			cbar.set_label(r'$\log_{10}(\xi)$ [log10]', fontsize=10)
			
			fig.tight_layout()
			fig.savefig(out_dir / f"surface_3d_with_csi_{model}.png", dpi=dpi)
			plt.close(fig)


def plot_3d_surface_interactive(
	b_rows: Sequence[BSummaryRow],
	nlem_rows: Sequence[NlemSummaryRow],
	out_dir: Path,
) -> None:
	"""Gera plots 3D interativos (Plotly) de Massa vs Raio vs Campo Magnético"""
	if not PLOTLY_AVAILABLE:
		print("Plotly não está instalado. Pulando gráficos 3D interativos.")
		return
	
	ensure_dir(out_dir)
	models = sorted({r.model for r in nlem_rows} | {r.model for r in b_rows})

	for model in models:
		# ===== Sem CSI =====
		no_csi = [r for r in b_rows if r.model == model and r.b_value > 0 and 
		         np.isfinite(r.log_b_value) and np.isfinite(r.max_mass_msun) and 
		         np.isfinite(r.radius_at_max_km)]
		
		if no_csi:
			no_csi = sorted(no_csi, key=lambda x: x.b_value)
			b_vals_no = np.array([r.b_value for r in no_csi], dtype=float)
			m_vals_no = np.array([r.max_mass_msun for r in no_csi], dtype=float)
			r_vals_no = np.array([r.radius_at_max_km for r in no_csi], dtype=float)
			log_b_vals = np.log10(b_vals_no)
			
			fig = go.Figure(data=[go.Scatter3d(  # type: ignore[union-attr,attr-defined]
				x=b_vals_no,
				y=m_vals_no,
				z=r_vals_no,
				mode='markers',
				marker={
					'size': 6,
					'color': log_b_vals,
					'colorscale': 'Viridis',
					'showscale': True,
					'colorbar': {'title': r'log₁₀(B)<br>[log10(G)]'},
					'line': {'width': 0.5, 'color': 'black'},
					'opacity': 0.8,
				},
				text=[f"B={b:.2e} G<br>M={m:.3f} M☉<br>R={r:.2f} km<br>log₁₀(B)={lb:.2f}" 
				      for b, m, r, lb in zip(b_vals_no, m_vals_no, r_vals_no, log_b_vals)],
				hoverinfo='text'
			)])
			
			fig.update_layout(  # type: ignore[union-attr]
				title=f'Massa vs Raio vs Campo Magnético | {model} | sem CSI',
				scene={
					'xaxis_title': 'B [G]',
					'yaxis_title': 'M_max [M☉]',
					'zaxis_title': 'R(M_max) [km]',
					'xaxis': {'type': 'log'},
				},
				width=1000,
				height=800,
				hovermode='closest'
			)
			
			fig.write_html(out_dir / f"surface_3d_interactive_no_csi_{model}.html")  # type: ignore[union-attr]
		
		# ===== Com CSI =====
		with_csi = [r for r in nlem_rows if r.model == model and r.b_value > 0 and
		           np.isfinite(r.max_mass_msun) and np.isfinite(r.radius_at_max_km)]
		
		if with_csi:
			with_csi = sorted(with_csi, key=lambda x: (x.b_value, x.log_csi))
			b_vals = np.array([r.b_value for r in with_csi], dtype=float)
			m_vals = np.array([r.max_mass_msun for r in with_csi], dtype=float)
			r_vals = np.array([r.radius_at_max_km for r in with_csi], dtype=float)
			log_csi_vals = np.array([r.log_csi for r in with_csi])
			
			fig = go.Figure(data=[go.Scatter3d(  # type: ignore[union-attr,attr-defined]
				x=b_vals,
				y=m_vals,
				z=r_vals,
				mode='markers',
				marker={
					'size': 6,
					'color': log_csi_vals,
					'colorscale': 'Plasma',
					'showscale': True,
					'colorbar': {'title': r'log₁₀(ξ)<br>[log10]'},
					'line': {'width': 0.5, 'color': 'black'},
					'opacity': 0.8,
				},
				text=[f"B={b:.2e} G<br>ξ={csi:.2e}<br>M={m:.3f} M☉<br>R={r:.2f} km<br>log₁₀(ξ)={lc:.2f}" 
				      for b, csi, m, r, lc in zip(b_vals, 10**log_csi_vals, m_vals, r_vals, log_csi_vals)],
				hoverinfo='text'
			)])
			
			fig.update_layout(  # type: ignore[union-attr]
				title=f'Massa vs Raio vs Campo Magnético | {model} | com CSI',
				scene={
					'xaxis_title': 'B [G]',
					'yaxis_title': 'M_max [M☉]',
					'zaxis_title': 'R(M_max) [km]',
					'xaxis': {'type': 'log'},
				},
				width=1000,
				height=800,
				hovermode='closest'
			)
			
			fig.write_html(out_dir / f"surface_3d_interactive_with_csi_{model}.html")  # type: ignore[union-attr]


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

		_append_gap_section(lines, no_csi, with_csi)

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

	plot_b_trends(b_rows, b_datasets, args.output_root / "trends", dpi=args.dpi)
	print(f"[OK] tendências sem csi: {args.output_root / 'trends'}")

	plot_all_variables_vs_logb(b_rows, args.output_root / "trends" / "all_variables", dpi=args.dpi)
	print(f"[OK] variáveis vs logB: {args.output_root / 'trends' / 'all_variables'}")

	plot_all_core_variables_vs_logb(
		b_datasets,
		args.output_root / "trends" / "all_core_variables",
		dpi=args.dpi,
	)
	print(f"[OK] variáveis do core vs logB: {args.output_root / 'trends' / 'all_core_variables'}")

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
	
	# Gráficos 3D de superfície
	surface_dir = comp_dir / "surfaces"
	plot_3d_surface_m_r_vs_b(
		b_rows=b_rows,
		nlem_rows=nlem_rows,
		out_dir=surface_dir,
		dpi=args.dpi,
	)
	print(f"[OK] plots 3D de superfícies: {surface_dir}")
	
	# Gráficos 3D interativos (Plotly)
	plot_3d_surface_interactive(
		b_rows=b_rows,
		nlem_rows=nlem_rows,
		out_dir=surface_dir,
	)
	if PLOTLY_AVAILABLE:
		print(f"[OK] plots 3D interativos (Plotly): {surface_dir}")

	relation_csv = comp_dir / "relation_with_without_csi.csv"
	write_direct_relation_csv(b_rows, nlem_rows, relation_csv)
	print(f"[OK] relação direta em CSV: {relation_csv}")

	extrema_md = args.output_root / "extrema_report.md"
	write_extrema_report(b_rows, nlem_rows, extrema_md)
	print(f"[OK] relatório de extremos: {extrema_md}")

	print("[DONE] Análise de B finalizada.")


if __name__ == "__main__":
	main()
