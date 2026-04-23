#!/usr/bin/env python3
"""
Gera gráficos massa-raio (M-R) separados por modelo (GM1 e GM3), com:

i)   sem csi em B = 0 (dados de output/b)
ii)  MODMAX em B = 1e17, 1e18 com gamma = 0
iii) MODMAX em B = 1e17, 1e18 com gamma = 0.5

Entradas:
- sem csi: output/b
- modmax:  output/modmax
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


MAX_VALID_MASS_MSUN = 3.0
MAX_VALID_RADIUS_KM = 20.0
EPS = 1e-12
EOS_FILENAME = "eos.dat"


@dataclass
class Curve:
	label: str
	source: str
	b_target: float
	b_used: float
	gamma: Optional[float]
	path: Path
	radius_km: np.ndarray
	mass_msun: np.ndarray


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Geração de gráficos M-R para GM1/GM3")
	parser.add_argument("--b-root", type=Path, default=Path("output/b"), help="Raiz dos dados sem csi")
	parser.add_argument("--modmax-root", type=Path, default=Path("output/modmax"), help="Raiz dos dados MODMAX")
	parser.add_argument("--output-dir", type=Path, default=Path("results/graficos_mario"))
	parser.add_argument("--dpi", type=int, default=180)
	return parser.parse_args()


def _safe_float(text: str) -> Optional[float]:
	try:
		return float(text)
	except Exception:
		return None


def _sci_folder_label(value: float, zero_label: str) -> str:
	if abs(value) < EPS:
		return zero_label
	label = f"{value:.2e}"
	label = label.replace("e+", "e")
	label = label.replace("e-0", "e-")
	label = label.replace("e0", "e0")
	return label


def _latex_sci(value: float, digits: int = 1) -> str:
	if abs(value) < EPS:
		return "0"
	s = f"{value:.{digits}e}"
	mant_str, exp_str = s.split("e")
	mant = float(mant_str)
	exp = int(exp_str)
	if abs(mant - 1.0) < EPS:
		return rf"10^{{{exp}}}"
	if abs(mant + 1.0) < EPS:
		return rf"-10^{{{exp}}}"
	return rf"{mant:g}\times10^{{{exp}}}"


def _valid_mr_mask(mass_col: np.ndarray, radius_col: np.ndarray) -> np.ndarray:
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

	# Convenção do projeto: últimas 2 colunas = massa e raio
	if arr.shape[1] < 5:
		return None

	mass = arr[:, -2].astype(float)
	radius = arr[:, -1].astype(float)
	mask = _valid_mr_mask(mass, radius)
	if not np.any(mask):
		return None

	return radius[mask], mass[mask]


def _parse_b_label_to_value(b_label: str) -> Optional[float]:
	return _safe_float(b_label)


def _discover_b_dirs_for_model(base_dir: Path, model: str) -> List[Tuple[float, Path]]:
	model_dir = base_dir / model
	if not model_dir.exists():
		return []

	out: List[Tuple[float, Path]] = []
	for child in model_dir.iterdir():
		if not child.is_dir() or not child.name.startswith("B_"):
			continue
		b_label = child.name[2:]
		b_value = _parse_b_label_to_value(b_label)
		if b_value is None:
			continue
		out.append((b_value, child))

	out.sort(key=lambda x: x[0])
	return out


def _closest_b_dir(base_dir: Path, model: str, b_target: float) -> Optional[Tuple[float, Path]]:
	candidates = _discover_b_dirs_for_model(base_dir, model)
	if not candidates:
		return None

	if abs(b_target) < EPS:
		for b_val, b_dir in candidates:
			if abs(b_val) < EPS:
				return b_val, b_dir

	return min(candidates, key=lambda x: abs(x[0] - b_target))


def _build_curve_sem_csi(b_root: Path, model: str, b_target: float, label: str) -> Optional[Curve]:
	hit = _closest_b_dir(b_root, model, b_target)
	if hit is None:
		return None

	b_used, b_dir = hit
	eos_path = b_dir / "default" / EOS_FILENAME
	if not eos_path.exists():
		return None

	mr = _load_mr_from_eos(eos_path)
	if mr is None:
		return None

	radius, mass = mr
	return Curve(
		label=label,
		source="sem_csi",
		b_target=b_target,
		b_used=b_used,
		gamma=None,
		path=eos_path,
		radius_km=radius,
		mass_msun=mass,
	)


def _resolve_modmax_eos_path(modmax_root: Path, model: str, b_value: float, gamma: float) -> Optional[Path]:
	b_label = _sci_folder_label(b_value, zero_label="0.00e0")
	b_dir = modmax_root / model / f"B_{b_label}"
	csi_label = _sci_folder_label(gamma, zero_label="0.00e0")
	direct = b_dir / "anisotropic" / f"csi_{csi_label}" / EOS_FILENAME
	if direct.exists():
		return direct

	anisotropic = b_dir / "anisotropic"
	if not anisotropic.exists():
		return None

	for child in anisotropic.iterdir():
		if not child.is_dir() or not child.name.startswith("csi_"):
			continue
		val = _safe_float(child.name[4:])
		if val is None:
			continue
		if abs(val - gamma) < EPS:
			candidate = child / EOS_FILENAME
			if candidate.exists():
				return candidate

	return None


def _build_curve_modmax(modmax_root: Path, model: str, b_value: float, gamma: float, label: str) -> Optional[Curve]:
	eos_path = _resolve_modmax_eos_path(modmax_root=modmax_root, model=model, b_value=b_value, gamma=gamma)
	if eos_path is None:
		return None

	mr = _load_mr_from_eos(eos_path)
	if mr is None:
		return None

	radius, mass = mr
	return Curve(
		label=label,
		source="modmax",
		b_target=b_value,
		b_used=b_value,
		gamma=gamma,
		path=eos_path,
		radius_km=radius,
		mass_msun=mass,
	)


def _ensure_dir(path: Path) -> None:
	path.mkdir(parents=True, exist_ok=True)


def _collect_model_curves(model: str, b_root: Path, modmax_root: Path, b_values: Sequence[float]) -> List[Curve]:
	curves: List[Curve] = []

	curve_b0 = _build_curve_sem_csi(
		b_root=b_root,
		model=model,
		b_target=0.0,
		label=r"$B=0\,\mathrm{G}$",
	)
	if curve_b0 is not None:
		curves.append(curve_b0)

	for b in b_values:
		c0 = _build_curve_modmax(
			modmax_root=modmax_root,
			model=model,
			b_value=b,
			gamma=0.0,
			label=rf"MODMAX | $B={_latex_sci(b)}\,\mathrm{{G}}$ | $\gamma=0.0$",
		)
		if c0 is not None:
			curves.append(c0)

		c05 = _build_curve_modmax(
			modmax_root=modmax_root,
			model=model,
			b_value=b,
			gamma=0.5,
			label=rf"MODMAX | $B={_latex_sci(b)}\,\mathrm{{G}}$ | $\gamma=0.5$",
		)
		if c05 is not None:
			curves.append(c05)

	return curves


def _plot_model_for_field(model: str, b_value: float, curves: Sequence[Curve], out_dir: Path, dpi: int) -> Optional[Path]:
	if not curves:
		return None

	plt.figure(figsize=(9, 7))
	for curve in curves:
		if curve.source == "sem_csi":
			color = "black"
		elif abs(curve.gamma or 0.0) < EPS:
			color = "tab:green"
		else:
			color = "tab:red"
		if curve.gamma is None:
			style = "-"
		elif abs(curve.gamma) < EPS:
			style = "--"
		else:
			style = "-."
		plt.plot(curve.radius_km, curve.mass_msun, style, lw=2.0, color=color, label=curve.label)

	plt.xlabel("R [km]")
	plt.ylabel(r"M [$M_\odot$]")
	plt.xlim(9, 14.5)
	plt.ylim(1.0, 2.25)
	plt.title(rf"Mass-Radius Diagram | {model} | $B={_latex_sci(b_value)}\,\mathrm{{G}}$")
	plt.grid(alpha=0.25)
	plt.legend(fontsize=9)
	plt.tight_layout()

	b_tag = f"{b_value:.0e}".replace("+", "")
	out_file = out_dir / f"mr_{model}_B_{b_tag}.png"
	plt.savefig(out_file, dpi=dpi)
	plt.close()
	return out_file


def _curves_for_field(curves: Sequence[Curve], b_value: float) -> List[Curve]:
	selected: List[Curve] = []
	for c in curves:
		if c.source == "sem_csi":
			selected.append(c)
			continue
		if abs(c.b_target - b_value) < EPS * max(1.0, abs(b_value)):
			selected.append(c)
	return selected


def _curve_stats(curve: Curve) -> Tuple[float, float, float, float, float, float]:
	mass_min = float(np.min(curve.mass_msun))
	mass_max = float(np.max(curve.mass_msun))
	radius_min = float(np.min(curve.radius_km))
	radius_max = float(np.max(curve.radius_km))
	idx_max_mass = int(np.argmax(curve.mass_msun))
	radius_at_max_mass = float(curve.radius_km[idx_max_mass])
	return mass_min, mass_max, radius_min, radius_max, radius_at_max_mass, float(curve.mass_msun[idx_max_mass])



def _write_values_table_csv(output_dir: Path, curves_by_model: Mapping[str, Sequence[Curve]]) -> Path:
	out_file = output_dir / "tabela_valores_utilizados.csv"
	with out_file.open("w", newline="", encoding="utf-8") as f:
		writer = csv.writer(f)
		writer.writerow(
			[
				"modelo",
				"fonte",
				"b_alvo_g",
				"b_utilizado_g",
				"gamma",
				"massa_min_msun",
				"massa_max_msun",
				"raio_min_km",
				"raio_max_km",
				"raio_em_massa_max_km",
				"arquivo",
			]
		)

		for model in sorted(curves_by_model.keys()):
			curves = curves_by_model[model]
			ordered = sorted(
				curves,
				key=lambda c: (
					0 if c.source == "sem_csi" else 1,
					c.b_target,
					-1.0 if c.gamma is None else c.gamma,
				),
			)
			for c in ordered:
				mass_min, mass_max, radius_min, radius_max, radius_at_max_mass, _ = _curve_stats(c)
				writer.writerow(
					[
						model,
						c.source,
						f"{c.b_target:.3e}",
						f"{c.b_used:.3e}",
						"" if c.gamma is None else f"{c.gamma:.1f}",
						f"{mass_min:.6f}",
						f"{mass_max:.6f}",
						f"{radius_min:.6f}",
						f"{radius_max:.6f}",
						f"{radius_at_max_mass:.6f}",
						c.path.as_posix(),
					]
				)

	return out_file


def main() -> None:
	args = parse_args()
	_ensure_dir(args.output_dir)

	models = ["GM1", "GM3"]
	b_values = [1e17, 1e18, 3e18]
	curves_by_model: Dict[str, List[Curve]] = {}

	created_files: List[Path] = []
	for model in models:
		curves = _collect_model_curves(
			model=model,
			b_root=args.b_root,
			modmax_root=args.modmax_root,
			b_values=b_values,
		)
		curves_by_model[model] = curves

		for b in b_values:
			curves_b = _curves_for_field(curves, b)
			out = _plot_model_for_field(model=model, b_value=b, curves=curves_b, out_dir=args.output_dir, dpi=args.dpi)
			if out is not None:
				created_files.append(out)

	if created_files:
		print("Arquivos gerados:")
		for p in created_files:
			print(f" - {p.as_posix()}")
	else:
		print("Nenhum gráfico foi gerado. Verifique os diretórios de entrada.")

	csv_file = _write_values_table_csv(args.output_dir, curves_by_model)
	print(f"Tabela CSV: {csv_file.as_posix()}")


if __name__ == "__main__":
	main()
