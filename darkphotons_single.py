from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext, ScalarFormatter
import numpy as np

plt.rcParams["axes.formatter.use_mathtext"] = True
plt.rcParams["axes.formatter.useoffset"] = False

N0 = 0.153
M_NUCLEON = 939.56534623

INPUT_ROOT = Path("output/darkphotons")
EOS_FILES = [
    ("darkphoton.dat", "Com Foton Escuro (DM)", "purple", "-"),
    ("eos_hadrons.dat", "Hadronica Pura", "black", "--"),
]

BARYON_COLUMNS = {
    5: 0.0,   # n
    6: 1.0,   # p
    7: 0.0,   # L0
    8: -1.0,  # S-
    9: 0.0,   # S0
    10: 1.0,  # S+
    11: -1.0, # X-
    12: 0.0,  # X0
}

PARTICLE_COLUMNS = {
    "e": 3,
    "mu": 4,
    "n": 5,
    "p": 6,
    "Lambda0": 7,
    "Sigma-": 8,
    "Sigma0": 9,
    "Sigma+": 10,
    "Xi-": 11,
    "Xi0": 12,
}

PARTICLE_STYLES = {
    "e": "#1f77b4",
    "mu": "#17becf",
    "n": "#111111",
    "p": "#d62728",
    "Lambda0": "#2ca02c",
    "Sigma-": "#9467bd",
    "Sigma0": "#8c564b",
    "Sigma+": "#e377c2",
    "Xi-": "#ff7f0e",
    "Xi0": "#7f7f7f",
}


def math_scalar_formatter() -> ScalarFormatter:
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))
    formatter.set_useOffset(False)
    return formatter


def apply_scientific_ticks(*, log_x: bool = False, log_y: bool = False) -> None:
    ax = plt.gca()
    if log_x:
        ax.xaxis.set_major_formatter(LogFormatterMathtext())
    else:
        ax.xaxis.set_major_formatter(math_scalar_formatter())
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2), useMathText=True, useOffset=False)

    if log_y:
        ax.yaxis.set_major_formatter(LogFormatterMathtext())
    else:
        ax.yaxis.set_major_formatter(math_scalar_formatter())
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True, useOffset=False)


def parse_header(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    return line
                if line.strip():
                    break
    except OSError:
        pass
    return ""


def reconstruct_mu_total_over_mn(values: list[float]) -> float:
    nb_total = sum(values[col] for col in BARYON_COLUMNS)
    if nb_total <= 0.0:
        return np.nan

    mu_n = values[17]
    mu_e = values[18]
    baryon_mu_density = 0.0

    for col, charge in BARYON_COLUMNS.items():
        mu_i = mu_n - charge * mu_e
        baryon_mu_density += mu_i * values[col]

    lepton_mu_density = mu_e * (values[3] + values[4])
    return (baryon_mu_density + lepton_mu_density) / nb_total


def read_curve(path: Path) -> dict[str, np.ndarray]:
    header = parse_header(path)
    col20_is_mu_total = "20:mu_total" in header

    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) not in (24, 37):
                continue

            try:
                values = [float(part) for part in parts]
            except ValueError:
                continue

            mass = values[-3]
            radius = values[-2]
            if not (np.isfinite(mass) and np.isfinite(radius)):
                continue
            if mass <= 0.0 or mass >= 3.0 or radius <= 0.0 or radius >= 30.0:
                continue

            if col20_is_mu_total:
                mu_total = values[20]
            else:
                mu_total = reconstruct_mu_total_over_mn(values)

            rows.append(
                [
                    values[0],   # nB / n0
                    values[1],   # eps
                    values[2],   # P
                    values[17],  # mu_n / mN
                    values[18],  # mu_e / mN
                    mu_total,    # mu_total / mN
                    mass,
                    radius,
                    *[values[col] for col in PARTICLE_COLUMNS.values()],
                ]
            )

    if not rows:
        raise ValueError(f"No usable data in {path}")

    arr = np.asarray(rows, dtype=float)
    return {
        "dens": arr[:, 0],
        "energ": arr[:, 1],
        "press": arr[:, 2],
        "mu_b": arr[:, 3],
        "mu_e": arr[:, 4],
        "mu_total": arr[:, 5],
        "masses": arr[:, 6],
        "radii": arr[:, 7],
        "particles": {
            name: arr[:, 8 + idx]
            for idx, name in enumerate(PARTICLE_COLUMNS)
        },
    }


def save_current(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def add_particle_legends(active_particles: list[str], source_styles: list[tuple[str, str]]) -> None:
    ax = plt.gca()
    particle_handles = [
        Line2D([0], [0], color=PARTICLE_STYLES[particle], lw=2.4, label=particle)
        for particle in active_particles
    ]
    source_handles = [
        Line2D([0], [0], color="black", lw=2.4, linestyle=linestyle, label=label)
        for label, linestyle in source_styles
    ]

    particle_legend = ax.legend(
        handles=particle_handles,
        title="Particulas",
        fontsize=8,
        title_fontsize=9,
        ncol=2,
        loc="lower left",
    )
    ax.add_artist(particle_legend)
    ax.legend(
        handles=source_handles,
        title="Fonte",
        fontsize=8,
        title_fontsize=9,
        loc="upper right",
    )


def plot_model(model_dir: Path) -> None:
    curves = []
    for filename, label, color, linestyle in EOS_FILES:
        path = model_dir / filename
        if not path.exists():
            print(f"Skipping missing file: {path}")
            continue

        try:
            curve = read_curve(path)
        except ValueError as exc:
            print(exc)
            continue

        curves.append((curve, label, color, linestyle))

    if not curves:
        return

    model = model_dir.name

    plt.figure(figsize=(8, 6))
    for curve, label, color, linestyle in curves:
        plt.plot(
            curve["radii"],
            curve["masses"],
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=2.5,
        )

        max_mass_idx = int(np.argmax(curve["masses"]))
        max_m = curve["masses"][max_mass_idx]
        corresponding_r = curve["radii"][max_mass_idx]
        plt.scatter(corresponding_r, max_m, color=color, s=50, zorder=5)
        plt.annotate(
            f"$M_{{max}} = {max_m:.5f} M_\\odot$\n$R = {corresponding_r:.5f}$ km",
            xy=(corresponding_r, max_m),
            xytext=(-25 if color == "purple" else 50, 35),
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            color=color,
            arrowprops=dict(
                arrowstyle="->",
                color=color,
                lw=1.5,
                connectionstyle="arc3,rad=.2",
            ),
        )

    plt.axhline(y=2.08, color="gray", linestyle=":", label="PSR J0740+6620")
    plt.axhspan(2.01, 2.15, color="gray", alpha=0.2)
    plt.xlabel(r"Raio $R$ (km)", fontsize=14)
    plt.ylabel(r"Massa $M$ ($M_\odot$)", fontsize=14)
    plt.title(f"Relacao Massa-Raio: Portal Vetorial - {model}", fontsize=16)
    all_radii = np.concatenate([curve["radii"] for curve, *_ in curves])
    all_masses = np.concatenate([curve["masses"] for curve, *_ in curves])
    radius_pad = max(0.5, 0.05 * (np.nanmax(all_radii) - np.nanmin(all_radii)))
    mass_pad = max(0.1, 0.05 * (np.nanmax(all_masses) - np.nanmin(all_masses)))
    plt.xlim(max(0.0, np.nanmin(all_radii) - radius_pad), np.nanmax(all_radii) + radius_pad)
    plt.ylim(max(0.0, np.nanmin(all_masses) - mass_pad), max(2.3, np.nanmax(all_masses) + mass_pad))
    plt.tick_params(axis="both", which="major", labelsize=12)
    plt.legend(loc="lower left", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.5)
    apply_scientific_ticks()
    save_current(model_dir / "mr_comparativo.png")

    plt.figure(figsize=(8, 6))
    for curve, label, color, linestyle in curves:
        plt.plot(curve["energ"], curve["press"], label=label, color=color, linestyle=linestyle, linewidth=2.5)
    plt.xlabel(r"energia [MeV/fm$^3$]", fontsize=14)
    plt.ylabel(r"pressao [MeV/fm$^3$]", fontsize=14)
    plt.title(f"EoS - {model}")
    plt.legend()
    apply_scientific_ticks()
    save_current(model_dir / "eos.png")

    plt.figure(figsize=(8, 6))
    for curve, label, color, linestyle in curves:
        positive = curve["press"] > 0.0
        plt.plot(
            curve["dens"][positive],
            np.log10(curve["press"][positive]),
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=2.5,
        )
    plt.ylabel(r"$\log_{10}(P)$", fontsize=14)
    plt.xlabel(r"$n_B / n_0$", fontsize=14)
    plt.title(f"Flow experiment - {model}")
    plt.legend()
    apply_scientific_ticks()
    save_current(model_dir / "pxrho.png")

    plt.figure(figsize=(8, 6))
    for curve, label, color, linestyle in curves:
        n_b_fm3 = curve["dens"] * N0
        mu_euler = (curve["press"] + curve["energ"]) / n_b_fm3
        plt.plot(curve["dens"], mu_euler, label=label, color=color, linestyle=linestyle, linewidth=2.5)
    plt.ylabel(r"$(P + E) / n_B$ [MeV]", fontsize=14)
    plt.xlabel(r"$n_B / n_0$", fontsize=14)
    plt.title(f"Euler chemical potential - {model}")
    plt.legend()
    apply_scientific_ticks()
    save_current(model_dir / "perho.png")

    plt.figure(figsize=(8, 6))
    for curve, label, color, linestyle in curves:
        n_b_fm3 = curve["dens"] * N0
        mu_euler = (curve["press"] + curve["energ"]) / n_b_fm3
        mu_total_mev = curve["mu_total"] * M_NUCLEON
        valid = np.isfinite(mu_total_mev) & (mu_total_mev > 0.0)
        plt.plot(
            curve["dens"][valid],
            mu_euler[valid] / mu_total_mev[valid],
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=2.5,
        )
    plt.axhline(1.0, color="gray", linestyle=":", linewidth=1.5)
    plt.ylabel(r"$[(P + E)/n_B] / \mu_{\rm total}$", fontsize=14)
    plt.xlabel(r"$n_B / n_0$", fontsize=14)
    plt.title(f"Thermodynamic consistency - {model}")
    plt.legend()
    apply_scientific_ticks()
    save_current(model_dir / "mu_txrho.png")

    source_styles = [(label, linestyle) for _, label, _, linestyle in curves]

    plt.figure(figsize=(9, 6))
    active_particles: list[str] = []
    for curve, source_label, _, source_linestyle in curves:
        for particle, density in curve["particles"].items():
            if np.nanmax(density) <= 0.0:
                continue
            if particle not in active_particles:
                active_particles.append(particle)
            plt.plot(
                curve["dens"],
                density,
                color=PARTICLE_STYLES[particle],
                linestyle=source_linestyle,
                linewidth=1.8,
                alpha=0.9,
            )
    plt.yscale("log")
    plt.ylabel(r"$n_i$ [fm$^{-3}$]", fontsize=14)
    plt.xlabel(r"$n_B / n_0$", fontsize=14)
    plt.title(f"Particle densities - {model}")
    add_particle_legends(active_particles, source_styles)
    apply_scientific_ticks(log_y=True)
    save_current(model_dir / "particle_densities.png")

    plt.figure(figsize=(9, 6))
    active_particles = []
    for curve, source_label, _, source_linestyle in curves:
        n_b_fm3 = curve["dens"] * N0
        valid_nb = n_b_fm3 > 0.0
        for particle, density in curve["particles"].items():
            fraction_percent = np.full_like(density, np.nan)
            fraction_percent[valid_nb] = 100.0 * density[valid_nb] / n_b_fm3[valid_nb]
            if np.nanmax(fraction_percent) <= 1.0e-6:
                continue
            if particle not in active_particles:
                active_particles.append(particle)
            plt.plot(
                curve["dens"],
                fraction_percent,
                color=PARTICLE_STYLES[particle],
                linestyle=source_linestyle,
                linewidth=1.8,
                alpha=0.9,
            )
    plt.yscale("log")
    plt.ylabel(r"$100\,n_i / n_B$ [%]", fontsize=14)
    plt.xlabel(r"$n_B / n_0$", fontsize=14)
    plt.title(f"Particle fractions [%] - {model}")
    add_particle_legends(active_particles, source_styles)
    apply_scientific_ticks(log_y=True)
    save_current(model_dir / "particle_fractions.png")


def main() -> None:
    if not INPUT_ROOT.exists():
        raise SystemExit(f"Missing input directory: {INPUT_ROOT}")

    model_dirs = sorted(path for path in INPUT_ROOT.iterdir() if path.is_dir())
    if not model_dirs:
        raise SystemExit(f"No model directories found in {INPUT_ROOT}")

    for model_dir in model_dirs:
        print(f"Analyzing {model_dir}")
        plot_model(model_dir)


if __name__ == "__main__":
    main()
