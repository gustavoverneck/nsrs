import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib.colors as mcolors

# Configurações de estilo científico
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 10,
    'axes.labelsize': 11, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'figure.dpi': 200
})

def calculate_gamma(n_n0, press):
    """Calcula Gamma = (n/P) * (dP/dn)"""
    mask = (press > 1e-4) & (n_n0 > 0.05)
    n, p = n_n0[mask], press[mask]
    if len(n) < 10: return None, None
    dp_dn = np.gradient(p, n)
    gamma = (n / p) * dp_dn
    return n, gamma

def find_threshold(n_n0, nb_tot, species_density, threshold=1e-5):
    """Identifica a densidade crítica n/n0 onde a fração excede o limiar."""
    yi = species_density / nb_tot
    idx = np.where(yi > threshold)[0]
    return n_n0[idx[0]] if len(idx) > 0 else np.nan

def main():
    hyp_map = {'$\Lambda^0$': 7, '$\Sigma^-$': 8, '$\Sigma^0$': 9, '$\Sigma^+$': 10, '$\Xi^-$': 11, '$\Xi^0$': 12}
    
    input_root = Path('output/lsv')
    output_dir = Path('results/lsv/insights')
    output_dir.mkdir(parents=True, exist_ok=True)

    for model_path in input_root.iterdir():
        if not model_path.is_dir(): continue
        for b_path in model_path.iterdir():
            if not b_path.is_dir(): continue
            
            b_str = b_path.name
            entries = []
            for csi_path in b_path.glob('csi_*'):
                eos_file = csi_path / 'eos.dat'
                if eos_file.exists():
                    try:
                        csi_val = float(csi_path.name.split('_')[1])
                        entries.append((csi_val, eos_file))
                    except: continue
            
            entries.sort(key=lambda x: x[0])
            if not entries: continue

            # --- PLOT 1: Limiares de Hiperonização NORMALIZADOS ---
            plt.figure(figsize=(8, 5)) # Aumentado para acomodar a legenda lateral
            thresholds_results = defaultdict(list)
            log_csis = [np.log10(e[0]) if e[0] > 0 else -25 for e in entries]

            for csi, path in entries:
                df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
                nb_tot = df.iloc[:, 5:13].sum(axis=1).values
                for name, col in hyp_map.items():
                    thresholds_results[name].append(find_threshold(df[0].values, nb_tot, df[col].values))

            for name, vals in thresholds_results.items():
                v = np.array(vals)
                mask = ~np.isnan(v)
                if np.any(mask):
                    # --- NORMALIZAÇÃO PERCENTUAL ---
                    # Define o primeiro valor válido (caso Maxwell ou xi min) como 100%
                    baseline = v[mask][0]
                    v_percent = (v[mask] / baseline) * 100
                    plt.plot(np.array(log_csis)[mask], v_percent, marker='o', markersize=4, label=name)

            plt.xlabel(r'$\log_{10}(\xi)$ [MeV$^{-1}$]')
            plt.ylabel(r'Threshold Shift [$n_{crit}(\xi) / n_{crit,0}$] (%)')
            plt.title(f'Normalized Hyperonization Thresholds - {b_str}')
            
            # --- LEGENDA FORA DO GRÁFICO ---
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=9)
            
            plt.grid(True, alpha=0.3, ls=':')
            plt.tight_layout() # Ajusta para a legenda não ser cortada
            plt.savefig(output_dir / f"thresholds_normalized_{model_path.name}_{b_str}.pdf")
            plt.close()

            # --- PLOT 2: Índice Adiabático (Geminado) ---
            fig, ax = plt.subplots(figsize=(7, 5))
            norm = mcolors.Normalize(vmin=min(log_csis), vmax=max(log_csis))
            cmap = plt.get_cmap('plasma')

            sample_indices = np.linspace(0, len(entries)-1, 6, dtype=int)
            for idx in sample_indices:
                csi, path = entries[idx]
                df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
                n_sub, gamma = calculate_gamma(df[0].values, df[2].values)
                if gamma is not None:
                    ax.plot(n_sub, gamma, color=cmap(norm(np.log10(csi) if csi > 0 else -25)), 
                            label=rf'$\log\xi={np.log10(csi):.1f}$' if csi > 0 else "Maxwell")

            ax.axhline(4/3, color='red', ls='--', alpha=0.6, label=r'$\Gamma = 4/3$')
            ax.set_xlabel(r'Density $n_B/n_0$')
            ax.set_ylabel(r'Adiabatic Index $\Gamma$')
            ax.set_ylim(0, 4)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            ax.grid(True, alpha=0.3, ls=':')
            plt.tight_layout()
            fig.savefig(output_dir / f"adiabatic_index_insights_{model_path.name}_{b_str}.pdf")
            plt.close()

    print(f"✅ Insights normalizados gerados em: {output_dir}")

if __name__ == "__main__":
    main()