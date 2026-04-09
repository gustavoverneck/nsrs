#!/usr/bin/env python3
"""
Análise científica completa dos dados NLEM gerados em output/limits.

Objetivos:
- Consolidar EoS, M-R e populações para GM1/GM3 e topologias isotrópica/anisotrópica.
- Quantificar o efeito de log10(csi) na estrutura estelar.
- Gerar tabelas e figuras prontas para publicação/inspeção posterior em Python.

Estrutura esperada dos dados (nova):
	output/limits/<MODEL>/B_<BVAL>/<isotropic|anisotropic>/csi_<CSI>/eos.dat

Também aceita estrutura legada sem pasta de topologia:
	output/limits/<MODEL>/B_<BVAL>/csi_<CSI>/eos.dat
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
COL_EMAG = 19


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
	eos_path: str


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Análise NLEM completa para dados em output/limits"
	)
	parser.add_argument(
		"--input-root",
		type=Path,
		default=Path("output/limits"),
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
	return parser.parse_args()


def ensure_dir(path: Path) -> None:
	path.mkdir(parents=True, exist_ok=True)


def _safe_float(text: str) -> Optional[float]:
	try:
		return float(text)
	except ValueError:
		return None


def parse_metadata_from_path(path: Path) -> Optional[Tuple[str, str, float, str, float, float]]:
	"""
	Retorna: (model, b_label, b_value, topology, csi, log_csi)
	"""
	rx = re.compile(
		r".*/output/limits/(?P<model>GM[13])/B_(?P<b>[^/]+)/(?:"
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

	# Esperamos pelo menos as 20 colunas físicas + 2 colunas MR
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


def summarize_dataset(ds: Dataset) -> SummaryRow:
	arr = ds.data

	mass_col = arr[:, -2]
	radius_col = arr[:, -1]
	mr_mask = np.isfinite(mass_col) & np.isfinite(radius_col) & (mass_col > 0.0) & (radius_col > 0.0)

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
		central_nb_n0=central_nb,
		central_eps_mevfm3=central_eps,
		central_p_mevfm3=central_p,
		central_meff=central_meff,
		central_emag_mevfm3=central_emag,
		cs2_min=cs2_min,
		cs2_max=cs2_max,
		eos_path=ds.file_path.as_posix(),
	)


def write_summary_csv(rows: Sequence[SummaryRow], out_csv: Path) -> None:
	ensure_dir(out_csv.parent)
	with out_csv.open("w", newline="", encoding="utf-8") as f:
		w = csv.writer(f)
		w.writerow([
			"model",
			"b_label",
			"b_value",
			"topology",
			"csi",
			"log_csi",
			"n_eos_points",
			"n_mr_points",
			"max_mass_msun",
			"radius_at_max_km",
			"central_nb_n0",
			"central_eps_mevfm3",
			"central_p_mevfm3",
			"central_meff",
			"central_emag_mevfm3",
			"cs2_min",
			"cs2_max",
			"eos_path",
		])
		for r in rows:
			w.writerow([
				r.model,
				r.b_label,
				r.b_value,
				r.topology,
				r.csi,
				r.log_csi,
				r.n_eos_points,
				r.n_mr_points,
				r.max_mass_msun,
				r.radius_at_max_km,
				r.central_nb_n0,
				r.central_eps_mevfm3,
				r.central_p_mevfm3,
				r.central_meff,
				r.central_emag_mevfm3,
				r.cs2_min,
				r.cs2_max,
				r.eos_path,
			])


def group_key(ds: Dataset) -> Tuple[str, str, str]:
	return ds.model, ds.topology, ds.b_label


def _subsample_sorted(items: Sequence[Dataset], max_items: int) -> List[Dataset]:
	if len(items) <= max_items:
		return list(items)
	idx = np.linspace(0, len(items) - 1, max_items).round().astype(int)
	idx = np.unique(idx)
	return [items[i] for i in idx]


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
		arr_sorted = sorted(arr, key=lambda d: d.log_csi)
		plot_set = _subsample_sorted(arr_sorted, max_curves)
		cmap_eos = plt.get_cmap("viridis")
		cmap_mr = plt.get_cmap("plasma")

		cvals = np.array([d.log_csi for d in plot_set], dtype=float)
		cmin, cmax = float(np.min(cvals)), float(np.max(cvals))
		denom = (cmax - cmin) if (cmax - cmin) > 1e-12 else 1.0

		# EoS family
		plt.figure(figsize=(8, 6))
		for d in plot_set:
			color = cmap_eos((d.log_csi - cmin) / denom)
			plt.plot(d.data[:, COL_EPS], d.data[:, COL_P], color=color, alpha=0.8, lw=1.0)
		plt.xlabel(r"$\epsilon$ [MeV/fm$^3$]")
		plt.ylabel(r"$P$ [MeV/fm$^3$]")
		plt.title(f"EoS family | {model} | {topology} | B={b_label} G")
		plt.grid(alpha=0.25)
		plt.tight_layout()
		plt.savefig(out_dir / f"eos_family_{model}_{topology}_B_{b_label}.png", dpi=dpi)
		plt.close()

		# M-R family
		plt.figure(figsize=(8, 6))
		for d in plot_set:
			m = d.data[:, -2]
			r = d.data[:, -1]
			mask = np.isfinite(m) & np.isfinite(r) & (m > 0.0) & (r > 0.0)
			if np.count_nonzero(mask) < 3:
				continue
			color = cmap_mr((d.log_csi - cmin) / denom)
			plt.plot(r[mask], m[mask], color=color, alpha=0.8, lw=1.0)
		plt.xlabel("Radius [km]")
		plt.ylabel(r"Mass [$M_\odot$]")
		plt.title(f"M-R family | {model} | {topology} | B={b_label} G")
		plt.grid(alpha=0.25)
		plt.tight_layout()
		plt.savefig(out_dir / f"mr_family_{model}_{topology}_B_{b_label}.png", dpi=dpi)
		plt.close()


def plot_population_snapshots(
	datasets: Sequence[Dataset],
	out_dir: Path,
	dpi: int,
) -> None:
	"""
	Para cada combinação (modelo, topologia, B), salva 3 snapshots de populações
	(csi baixo, intermediário e alto), com todas as espécies principais.
	"""
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
		arr_sorted = sorted(arr, key=lambda d: d.log_csi)
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

			ax.set_title(f"log10(csi)={ds.log_csi:.3f}")
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


def _plot_metric_for_subset(
	rows: Sequence[SummaryRow],
	model: str,
	topo: str,
	metric_name: str,
	y_label: str,
	out_path: Path,
	dpi: int,
) -> None:
	subset = [r for r in rows if r.model == model and r.topology == topo]
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
		plt.plot(x[mask], y[mask], marker="o", ms=3, lw=1.2, label=f"B={part[0].b_label} G")

	plt.xlabel(r"$\log_{10}(\xi)$")
	plt.ylabel(y_label)
	plt.title(f"{y_label} vs log10(csi) | {model} | {topo}")
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
		("max_mass_msun", r"$M_{\max}$ [$M_\odot$]", "max_mass_vs_logcsi"),
		("radius_at_max_km", r"$R(M_{\max})$ [km]", "radius_at_max_vs_logcsi"),
		("central_nb_n0", r"$n_c/n_0$", "central_density_vs_logcsi"),
		("central_emag_mevfm3", r"$\epsilon_{mag,c}$ [MeV/fm$^3$]", "central_emag_vs_logcsi"),
	]

	for metric_name, y_label, prefix in metrics:
		for model in models:
			for topo in topologies:
				_plot_metric_for_subset(
					rows=rows,
					model=model,
					topo=topo,
					metric_name=metric_name,
					y_label=y_label,
					out_path=out_dir / f"{prefix}_{model}_{topo}.png",
					dpi=dpi,
				)


def write_scientific_report(rows: Sequence[SummaryRow], out_file: Path) -> None:
	ensure_dir(out_file.parent)

	models = sorted({r.model for r in rows})
	tops = sorted({r.topology for r in rows})

	lines: List[str] = []
	lines.append("# NLEM neutron-star analysis report")
	lines.append("")
	lines.append("## Scope")
	lines.append("This report quantifies how $\\log_{10}(\\xi)$ modifies stellar structure using generated EoS/M-R data.")
	lines.append("")
	lines.append(f"- Total datasets analyzed: **{len(rows)}**")
	lines.append(f"- Models: **{', '.join(models)}**")
	lines.append(f"- Topologies: **{', '.join(tops)}**")
	lines.append("")

	def slope_for(sub: List[SummaryRow], metric: str) -> float:
		x = np.array([r.log_csi for r in sub], dtype=float)
		y = np.array([getattr(r, metric) for r in sub], dtype=float)
		mask = np.isfinite(x) & np.isfinite(y)
		if np.count_nonzero(mask) < 2:
			return math.nan
		a, _b = np.polyfit(x[mask], y[mask], 1)
		return float(a)

	lines.append("## Trend diagnostics by $(model, topology, B)$")
	lines.append("")
	lines.append("Linear slopes (first-order sensitivity):")
	lines.append("- $dM_{max}/d\\log_{10}(\\xi)$")
	lines.append("- $dR(M_{max})/d\\log_{10}(\\xi)$")
	lines.append("")

	grouped: Dict[Tuple[str, str, str], List[SummaryRow]] = {}
	for r in rows:
		grouped.setdefault((r.model, r.topology, r.b_label), []).append(r)

	for (model, topo, b_label), sub in sorted(grouped.items(), key=lambda kv: (kv[0][0], kv[0][1], float(kv[0][2]))):
		sub = sorted(sub, key=lambda x: x.log_csi)
		s_mass = slope_for(sub, "max_mass_msun")
		s_rad = slope_for(sub, "radius_at_max_km")

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

	lines.append("## Interpretation guide")
	lines.append("- Positive $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ tends to stiffen the effective sequence in your setup.")
	lines.append("- Negative $dM_{max}/d\\log_{10}(\\xi)$: larger $\\xi$ softens the sequence.")
	lines.append("- Compare isotropic vs anisotropic at fixed $(model, B)$ to isolate topology effects.")
	lines.append("- Use population snapshots to identify threshold shifts for hyperon onset with changing $\\xi$.")

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
	print(f"[OK] Famílias EoS/MR: {families_dir}")

	pops_dir = args.output_root / "populations"
	plot_population_snapshots(datasets, pops_dir, dpi=args.dpi)
	print(f"[OK] Snapshots de populações: {pops_dir}")

	report_md = args.output_root / "scientific_report.md"
	write_scientific_report(summary_rows, report_md)
	print(f"[OK] Relatório científico: {report_md}")

	print("[DONE] Análise completa finalizada.")


if __name__ == "__main__":
	main()

