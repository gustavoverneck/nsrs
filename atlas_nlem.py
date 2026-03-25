import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import math
from pathlib import Path
from collections import defaultdict
import sys

# ==========================================
# 1. Configurações de Estilo
# ==========================================
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 10,
    'axes.labelsize': 11, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'figure.dpi': 200, 'lines.linewidth': 1.1
})

M_NUC = 939.565346 
HYP_MAP = {r'$\Lambda^0$': 7, r'$\Sigma^-$': 8, r'$\Sigma^0$': 9, r'$\Sigma^+$': 10, r'$\Xi^-$': 11, r'$\Xi^0$': 12}

LINE_STYLES = {
    'n': '-', 'p': '--',          
    'e-': ':', 'mu-': '-.',       
    'L0': '-', 'S-': '--', 'S0': ':', 'S+': '-.', 'X-': '-', 'X0': '--' 
}

def calc_nu_max(mue, B_gauss):
    me, Bc = 0.51099895, 4.414e13
    eB_MeV2 = (B_gauss / Bc) * (me**2)
    return np.floor(np.maximum((mue**2 - me**2) / (2 * eB_MeV2), 0))

def find_threshold(n_n0, nb_tot, species_density, threshold=1e-5):
    with np.errstate(divide='ignore', invalid='ignore'):
        yi = np.true_divide(species_density, nb_tot)
        yi[~np.isfinite(yi)] = 0
    idx = np.where(yi > threshold)[0]
    return n_n0[idx[0]] if len(idx) > 0 else np.nan

def find_dir(base_path, prefix, target_val):
    if not base_path.exists(): return None
    for d in base_path.iterdir():
        if d.is_dir() and d.name.startswith(prefix):
            try:
                val = float(d.name.replace(prefix, ""))
                if abs(val - target_val) / (target_val + 1e-20) < 1e-3: return d
            except: continue
    return None

def main():
    if len(sys.argv) < 5:
        print("Uso: uv run atlas_nlem.py <MODELO> <EXP_MIN> <EXP_MAX> <B1> <B2> ...")
        return

    model_name = sys.argv[1]
    exp_min, exp_max = int(sys.argv[2]), int(sys.argv[3])
    b_fields = [float(b) for b in sys.argv[4:]]

    input_root = Path(f'output/limits/{model_name}')
    output_root = Path('results/analise_final')
    output_root.mkdir(parents=True, exist_ok=True)

    csi_vals = []
    for exp in range(exp_min, exp_max + 1):
        csi_vals.append(10.0**exp)
        if exp < exp_max: csi_vals.append(5.0 * 10.0**exp)
    
    indices_plot = np.linspace(0, len(csi_vals) - 1, 6, dtype=int)

    for B_val in b_fields:
        b_str = f"{B_val:.2e}"
        b_dir = find_dir(input_root, "B_", B_val)
        if not b_dir: continue

        print(f"📊 Processando B = {b_str} G")
        norm = mcolors.Normalize(vmin=exp_min, vmax=exp_max)
        cmap = plt.get_cmap('magma')

        fig_triple, (ax_eps, ax_p, ax_yi_t) = plt.subplots(3, 1, figsize=(7, 10), sharex=True, gridspec_kw={'hspace': 0.05})
        fig_pop_cs, (ax_pop_c, ax_cs2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True, gridspec_kw={'hspace': 0.05})
        fig_lan, ax_lan = plt.subplots(figsize=(7, 5))
        
        thresholds_data = defaultdict(list)

        for i, csi in enumerate(csi_vals):
            csi_dir = find_dir(b_dir, "csi_", csi)
            if not csi_dir: continue
            
            eos_file = csi_dir / "eos.dat"
            if not eos_file.exists(): continue

            df = pd.read_csv(eos_file, sep=r'\s+', header=None, comment='#')
            l_csi = math.log10(float(csi))
            color = cmap(norm(l_csi))
            n_n0, eps, press = df[0].values, df[1].values, df[2].values
            nb_tot = df.iloc[:, 5:13].sum(axis=1).values

            if i in indices_plot:
                mask = press > 1e-3
                if len(n_n0[mask]) > 0:
                    ax_eps.plot(n_n0[mask], eps[mask], color=color, alpha=0.8)
                    ax_p.plot(n_n0[mask], press[mask], color=color, alpha=0.8)
                ax_lan.step(n_n0, calc_nu_max(df[18].values*M_NUC, B_val), where='post', color=color, label=rf"$10^{{{l_csi:.1f}}}$")

            with np.errstate(divide='ignore', invalid='ignore'):
                parts_to_plot = {'n':5, 'p':6, 'e-':3, 'mu-':4, 'L0':7, 'S-':8, 'S0':9, 'S+':10, 'X-':11, 'X0':12}
                for p_name, col in parts_to_plot.items():
                    yi = np.true_divide(df[col].values, nb_tot)
                    yi[~np.isfinite(yi)] = 1e-6
                    yi = np.maximum(yi, 1e-6)
                    if i in indices_plot:
                        ls = LINE_STYLES.get(p_name, '-')
                        ax_yi_t.plot(n_n0, yi, color=color, alpha=0.5, lw=1.0, linestyle=ls)
                        ax_pop_c.plot(n_n0, yi, color=color, alpha=0.5, lw=1.0, linestyle=ls)

            mask_cs = (n_n0 > 0.1) & (eps > 0) & (press > 0)
            if len(n_n0[mask_cs]) > 2:
                try:
                    cs2 = np.gradient(press[mask_cs], eps[mask_cs])
                    ax_cs2.plot(n_n0[mask_cs], cs2, color=color, alpha=0.4)
                except ValueError: pass

            for name, col in HYP_MAP.items():
                thresholds_data[name].append(find_threshold(n_n0, nb_tot, df[col].values))

        # --- FINALIZAÇÃO TRIPLE ---
        ax_eps.set_title(f"Model: {model_name} | Field: {b_str} G", fontsize=12)
        ax_eps.set_yscale('log'); ax_p.set_yscale('log'); ax_yi_t.set_yscale('log')
        ax_eps.set_xlim(0, 8.0); ax_cs2.set_ylabel(r'$c_s^2$')
        ax_eps.set_ylabel(r'$\varepsilon$ [MeV/fm$^3$]'); ax_p.set_ylabel(r'$P$ [MeV/fm$^3$]'); ax_yi_t.set_ylabel(r'$Y_i$')
        fig_triple.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=[ax_eps, ax_p, ax_yi_t], label=r'$\log_{10}(\xi)$')
        fig_triple.savefig(output_root / f"atlas_macro_{b_str}.pdf", bbox_inches='tight')
        
        # --- FINALIZAÇÃO POP/CS2 ---
        ax_pop_c.set_title(f"Populations & $c_s^2$ - {model_name} ({b_str} G)")
        ax_pop_c.set_yscale('log'); ax_pop_c.set_ylim(1e-5, 1.2)
        ax_cs2.set_ylim(0, 1.0); ax_cs2.set_ylabel(r'$c_s^2$')
        ax_cs2.set_xlim(0, 8.0); ax_cs2.set_ylabel(r'$c_s^2$')
        ax_cs2.axhline(1.0, color='r', ls='--', alpha=0.5)
        ax_cs2.axhline(1/3, color='gray', ls=':', alpha=0.5)
        fig_pop_cs.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=[ax_pop_c, ax_cs2], label=r'$\log_{10}(\xi)$')
        fig_pop_cs.savefig(output_root / f"pop_cs2_{b_str}.pdf", bbox_inches='tight')

        # --- FINALIZAÇÃO LANDAU ---
        ax_lan.set_title(f"Landau Levels - {model_name} ({b_str} G)")
        ax_lan.set_ylabel(r'$\nu_{max}^{(e)}$'); ax_lan.set_xlabel(r'$n_B/n_0$')
        ax_lan.legend(ncol=2, fontsize=8)
        fig_lan.savefig(output_root / f"landau_sweep_{b_str}.pdf")

        # --- FINALIZAÇÃO THRESHOLDS (CORREÇÃO DO GRID AQUI) ---
        fig_t, ax_t = plt.subplots(figsize=(8, 5))
        ax_t.set_title(f"Hyperon Threshold Shift - {model_name} ({b_str} G)")
        log_axis = [math.log10(float(c)) for c in csi_vals]
        for name, vals in thresholds_data.items():
            v = np.array(vals)
            mask_v = ~np.isnan(v)
            if np.any(mask_v):
                ax_t.plot(np.array(log_axis)[mask_v], (v[mask_v] / v[mask_v][0])*100, marker='o', ms=4, label=name)
        ax_t.set_ylabel("Threshold Shift (%)"); ax_t.set_xlabel(r"$\log_{10}(\xi)$")
        ax_t.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        ax_t.grid(True, alpha=0.2) # CORRIGIDO: ax_t em vez de fig_t
        fig_t.savefig(output_root / f"thresholds_norm_{b_str}.pdf", bbox_inches='tight')
        
        plt.close('all')

    print(f"✅ Análise completa concluída em: {output_root}")

if __name__ == "__main__":
    main()