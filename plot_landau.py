import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
import os

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11,
    'axes.labelsize': 12, 'xtick.direction': 'in', 'ytick.direction': 'in'
})

def calc_nu_max(mue, B_gauss):
    me = 0.51099895 
    Bc = 4.414e13   
    eB_MeV2 = (B_gauss / Bc) * (me**2)
    nu_raw = (mue**2 - me**2) / (2 * eB_MeV2)
    return np.floor(np.maximum(nu_raw, 0))

# --- BUSCA INTELIGENTE DE PASTAS ---
def find_b_dir(base_dir, target_b):
    """ Encontra a pasta B_... que corresponde numericamente ao campo B desejado. """
    for d in base_dir.iterdir():
        if d.is_dir() and d.name.startswith("B_"):
            try:
                val = float(d.name[2:])
                if abs(val - target_b) / target_b < 1e-3: # Tolerância de 0.1%
                    return d
            except ValueError: continue
    return None

def find_csi_dir(b_dir, target_csi):
    """ Encontra a pasta csi_... que corresponde numericamente ao csi desejado. """
    for d in b_dir.iterdir():
        if d.is_dir() and d.name.startswith("csi_"):
            try:
                val = float(d.name[4:])
                if target_csi == 0.0 and val == 0.0: return d
                if target_csi != 0.0 and val != 0.0:
                    if abs(val - target_csi) / target_csi < 1e-3: # Tolerância de 0.1%
                        return d
            except ValueError: continue
    return None
# -----------------------------------

def main():
    if len(sys.argv) < 6:
        print("Uso: python lsv_landau_plots.py <GM1|GM3> <exp_min> <exp_max> <num_pontos> <B1> <B2> ...")
        print("Exemplo: python lsv_landau_plots.py GM1 -5.0 -1.5 80 1e17 5e17 1e18")
        return

    model_name = sys.argv[1]
    exp_min = float(sys.argv[2])
    exp_max = float(sys.argv[3])
    num_points = int(sys.argv[4])
    b_fields = [float(b) for b in sys.argv[5:]]

    base_dir = Path(f"output/lsv/{model_name}")
    if not base_dir.exists():
        print(f"Pasta {base_dir} não encontrada. Rode o script Rust primeiro!")
        return
    
    # Reconstrói os valores de csi exatamente como o Rust fez
    csi_vals = []
    if num_points == 1:
        csi_vals.append(10**exp_min)
    else:
        step = (exp_max - exp_min) / (num_points - 1)
        for i in range(num_points):
            csi_vals.append(10**(exp_min + i * step))
            
    # Filtra para desenhar apenas ~6 curvas, para não poluir o gráfico
    num_curves_to_plot = min(6, num_points)
    indices_to_plot = np.linspace(0, num_points - 1, num_curves_to_plot, dtype=int)
    
    fig, axes = plt.subplots(len(b_fields), 1, figsize=(8, 4 * len(b_fields)), sharex=True)
    if len(b_fields) == 1: axes = [axes]
    fig.subplots_adjust(hspace=0.1)
    
    cmap = plt.get_cmap('plasma')
    colors = [cmap(i) for i in np.linspace(0, 0.9, num_curves_to_plot)]
    M_NUC = 939.565346 # Constante para desnormalizar o potencial químico
    
    for ax, B_val in zip(axes, b_fields):
        b_dir = find_b_dir(base_dir, B_val)
        
        if b_dir is None:
            print(f"Aviso: Pasta para B = {B_val} G não encontrada.")
            continue
            
        for color_idx, i in enumerate(indices_to_plot):
            csi = csi_vals[i]
            csi_dir = find_csi_dir(b_dir, csi)
            
            if csi_dir is not None:
                path = csi_dir / "eos.dat"
                if path.exists():
                    df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
                    n_std, mue_norm = df[0].values, df[18].values
                    mue_MeV = mue_norm * M_NUC
                    
                    nu_max = calc_nu_max(mue_MeV, B_val)
                    
                    c = colors[color_idx]
                    ls = '-' if color_idx == 0 else '--'
                    lw = 2.0 if color_idx == 0 else 1.5
                    label = rf'$\xi = 10^{{{np.log10(csi):.2f}}}$'
                    
                    ax.step(n_std, nu_max, where='post', color=c, linestyle=ls, lw=lw, label=label)
            else:
                print(f"Aviso: Dados faltando para csi = {csi} no campo B = {B_val}")
        
        ax.set_ylabel(r'$\nu_{max}^{(e)}$')
        # Formata o texto do painel
        b_text = f"10^{{{np.log10(B_val):.2f}}}" if np.log10(B_val).is_integer() else f"{B_val:.1e}"
        ax.text(0.05, 0.85, f'$B = {b_text}$ G', transform=ax.transAxes, 
                fontsize=11, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
        if B_val >= 1e18: ax.set_ylim(-0.5, 12)
        elif B_val >= 1e17: ax.set_ylim(-1, 50)
        
        ax.legend(loc='lower right', fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3, ls=':')

    axes[-1].set_xlabel(r'Densidade Bariônica $n_B/n_0$')
    axes[-1].set_xlim(0.0, 9.0)
    
    out_dir = Path("results/lsv")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"landau_shift_sweep_{model_name}.pdf"
    
    plt.savefig(out_file, bbox_inches='tight')
    print(f"✅ Gráfico de quantização de Landau (Sweep) gerado: {out_file}")

if __name__ == "__main__":
    main()