import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path
from collections import defaultdict

# ==========================================
# 1. Configurações de Estilo Científico
# ==========================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'figure.dpi': 200,
    'lines.linewidth': 1.0
})

def clamp(val):
    return np.maximum(val, 1e-7)

# Paletas por partícula
particle_cmaps = {
    'n': 'Greys', 'p': 'Blues', 'e-': 'Reds', 'mu-': 'Greens',
    'L0': 'Purples', 'S-': 'Oranges', 'S0': 'YlGn', 'S+': 'PuRd',
    'X-': 'YlOrBr', 'X0': 'BuPu'
}

# Rótulos em LaTeX
latex_labels = {
    'n': r'$n$', 'p': r'$p$', 'e-': r'$e^-$', 'mu-': r'$\mu^-$',
    'L0': r'$\Lambda^0$', 'S-': r'$\Sigma^-$', 'S0': r'$\Sigma^0$',
    'S+': r'$\Sigma^+$', 'X-': r'$\Xi^-$', 'X0': r'$\Xi^0$'
}

# POSIÇÕES CALIBRADAS PELA IMAGEM DE REFERÊNCIA
fixed_label_positions = {
    'n':    (1.2, 0.85),
    'p':    (2.1, 0.45),
    'e-':   (3.3, 0.03),
    'mu-':  (4.3, 0.02),
    'L0':   (5.4, 0.25),
    'S+':   (5.1, 0.001),
    'S-':   (6.2, 0.007),
    'S0':   (6.8, 0.004),
    'X-':   (7.4, 0.002),
    'X0':   (7.8, 0.0008)
}

# ==========================================
# 2. Processamento Principal
# ==========================================
def main():
    input_root = Path('output/lsv')
    output_root = Path('results/lsv/python_plots')
    output_root.mkdir(parents=True, exist_ok=True)

    simulations = defaultdict(lambda: defaultdict(list))
    for eos_path in input_root.rglob('eos.dat'):
        parts = eos_path.parts
        try:
            model, b_str = parts[2], parts[3]
            csi_val = float(parts[4].replace('csi_', ''))
            simulations[model][b_str].append((np.log10(csi_val), eos_path))
        except: continue

    if not simulations:
        print("❌ Nenhum dado encontrado em output/lsv/")
        return

    for model, b_fields in simulations.items():
        for b_str, entries in b_fields.items():
            entries.sort(key=lambda x: x[0])
            print(f"📊 Processando {model} | {b_str} ({len(entries)} arquivos)")

            log_csi_list = [e[0] for e in entries]
            norm = mcolors.Normalize(vmin=min(log_csi_list), vmax=max(log_csi_list))
            global_cmap = plt.get_cmap('magma')

            # --- Criação das Figuras ---
            
            # 1. Antigo: Populações e Velocidade do Som (2 linhas)
            fig_combined, (ax_pop, ax_cs2) = plt.subplots(
                nrows=2, ncols=1, 
                figsize=(7, 8), 
                sharex=True,
                gridspec_kw={'height_ratios': [1.5, 1], 'hspace': 0.05} 
            )

            # 2. NOVO: Densidade de Energia, Pressão e Populações (3 linhas)
            fig_triple, (ax_eps, ax_press, ax_pop_triple) = plt.subplots(
                nrows=3, ncols=1, 
                figsize=(7, 10), 
                sharex=True, 
                gridspec_kw={'height_ratios': [1, 1, 1.5], 'hspace': 0.05} 
            )
            
            # Figuras avulsas restantes
            figs = {
                'eos': plt.subplots(figsize=(6, 5)),
                'meff': plt.subplots(figsize=(6, 5)),
                'mun': plt.subplots(figsize=(6, 5))
            }

            meff_data_all, mun_data_all = [], []

            for i, (log_csi, path) in enumerate(entries):
                try:
                    df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
                    num_cols = df.shape[1]
                    
                    n_n0, eps, press = df[0].values, df[1].values, df[2].values
                    nb_tot = df.iloc[:, 5:13].sum(axis=1).values
                    color = global_cmap(norm(log_csi))

                    # 1. Populações (Aplicado em ambos os gráficos)
                    mask_pop = (nb_tot > 1e-8)
                    parts_data = {
                        'n': df[5], 'p': df[6], 'e-': df[3], 'mu-': df[4],
                        'L0': df[7], 'S-': df[8], 'S0': df[9], 'S+': df[10],
                        'X-': df[11], 'X0': df[12]
                    }
                    for name, dens in parts_data.items():
                        yi = clamp(dens.values[mask_pop] / nb_tot[mask_pop])
                        if np.any(yi > 5e-4):
                            p_cmap = plt.get_cmap(particle_cmaps.get(name, 'viridis'))
                            p_color = p_cmap(norm(log_csi))
                            
                            # Plota no gráfico de 2 painéis
                            ax_pop.plot(n_n0[mask_pop], yi, color=p_color, alpha=0.3, lw=0.7)
                            # Plota no NOVO gráfico de 3 painéis
                            ax_pop_triple.plot(n_n0[mask_pop], yi, color=p_color, alpha=0.3, lw=0.7)

                    # 2. Velocidade do Som (Gráfico de baixo 2 painéis)
                    _, idx_unq = np.unique(eps, return_index=True)
                    idx_unq = np.sort(idx_unq)
                    if len(idx_unq) > 5:
                        cs2 = np.gradient(press[idx_unq], eps[idx_unq])
                        cs2 = np.where((cs2 > 0) & (cs2 < 1.2), cs2, np.nan)
                        ax_cs2.plot(n_n0[idx_unq], cs2, color=color, alpha=0.3, lw=0.8)

                    # 3. Macro EoS no Novo Gráfico Triplo (Eps e Pressão vs Densidade Bariônica)
                    # Usamos valores estritamente positivos para não quebrar a escala logarítmica
                    mask_macro = (eps > 1e-5) & (press > 1e-5)
                    ax_eps.plot(n_n0[mask_macro], eps[mask_macro], color=color, alpha=0.4, lw=0.8)
                    ax_press.plot(n_n0[mask_macro], press[mask_macro], color=color, alpha=0.4, lw=0.8)

                    # 4. EoS log-log
                    mask_eos = (eps > 0.1) & (press > 0.1)
                    figs['eos'][1].plot(eps[mask_eos], press[mask_eos], color=color, alpha=0.4, lw=0.8)

                    # 5. m* e mun
                    if num_cols >= 18:
                        m_eff_vals = df[16].values
                        mun_vals = df[17].values
                        meff_data_all.extend(m_eff_vals)
                        mun_data_all.extend(mun_vals)
                        figs['meff'][1].plot(n_n0, m_eff_vals, color=color, alpha=0.4)
                        figs['mun'][1].plot(n_n0, mun_vals, color=color, alpha=0.4)

                except Exception as e: 
                    continue

            # ========== FINALIZAÇÃO DOS GRÁFICOS COMBINADOS ==========
            
            # --- Formatação do Gráfico de 2 Painéis (Original) ---
            # for name, (px, py) in fixed_label_positions.items():
            #     ax_pop.text(px, py, latex_labels.get(name, name), fontsize=10, fontweight='bold', 
            #             ha='center', va='center', zorder=40,
            #             bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5))
            
            ax_pop.set_yscale('log')
            ax_pop.set_ylim(1e-5, 1.3)
            ax_pop.set_ylabel(r'Fraction $Y_i$')
            ax_pop.tick_params(labelbottom=False) 
            fig_combined.colorbar(plt.cm.ScalarMappable(cmap='Greys', norm=norm), ax=ax_pop, label=r'$\log_{10}(\xi)$', pad=0.02)

            ax_cs2.axhline(1.0, color='red', ls='--', lw=1, alpha=0.6, label='Causality')
            ax_cs2.axhline(1/3, color='gray', ls=':', lw=1, alpha=0.6, label='Conformal')
            ax_cs2.set_ylim(0, 1.1)
            ax_cs2.set_xlim(0, 9.25)
            ax_cs2.set_ylabel(r'$c_s^2$')
            ax_cs2.set_xlabel(r'Density $n_B/n_0$')
            ax_cs2.legend(loc='upper left', fontsize=8)
            fig_combined.colorbar(plt.cm.ScalarMappable(cmap=global_cmap, norm=norm), ax=ax_cs2, label=r'$\log_{10}(\xi)$', pad=0.02)

            fig_combined.savefig(output_root / f"pop_and_cs2_aligned_{model}_{b_str}.pdf", bbox_inches='tight')


            # --- Formatação do NOVO Gráfico de 3 Painéis ---
            ax_eps.set_ylabel(r'log $\varepsilon$ [MeV/fm$^3$]')
            ax_press.set_ylabel(r'$log P$ [MeV/fm$^3$]')
            
            # Aplicação da escala logarítmica (Descomentado!)
            ax_eps.set_yscale('log')
            ax_press.set_yscale('log')
            
            # Oculta os rótulos do eixo X nos dois gráficos de cima
            ax_eps.tick_params(labelbottom=False)
            ax_press.tick_params(labelbottom=False)

            # for name, (px, py) in fixed_label_positions.items():
            #     ax_pop_triple.text(px, py, latex_labels.get(name, name), fontsize=10, fontweight='bold', 
            #             ha='center', va='center', zorder=40,
            #             bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5))

            ax_pop_triple.set_yscale('log')
            ax_pop_triple.set_ylim(1e-5, 1.3)
            ax_pop_triple.set_ylabel(r'Fraction $Y_i$')
            ax_pop_triple.set_xlabel(r'Density $n_B/n_0$')
            ax_pop_triple.set_xlim(0, 9.25) # Limite padrão para alinhar com as populações

            # Adiciona uma única barra de cores lateral para a figura inteira
            fig_triple.colorbar(plt.cm.ScalarMappable(cmap=global_cmap, norm=norm), ax=[ax_eps, ax_press, ax_pop_triple], label=r'$\log_{10}(\xi)$', pad=0.02)

            fig_triple.savefig(output_root / f"macro_and_pop_aligned_{model}_{b_str}.pdf", bbox_inches='tight')


            # ========== FINALIZAÇÃO DOS OUTROS GRÁFICOS ==========
            
            # EoS
            ax = figs['eos'][1]
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlabel(r'$\varepsilon$ [MeV/fm$^3$]'); ax.set_ylabel(r'$P$ [MeV/fm$^3$]')
            figs['eos'][0].colorbar(plt.cm.ScalarMappable(cmap=global_cmap, norm=norm), ax=ax, label=r'$\log_{10}(\xi)$')
            figs['eos'][0].savefig(output_root / f"eos_loglog_{model}_{b_str}.pdf", bbox_inches='tight')

            # Massa Efetiva
            ax = figs['meff'][1]
            if meff_data_all:
                ax.set_ylim(np.min(meff_data_all)*0.95, np.max(meff_data_all)*1.05)
            ax.set_xlabel(r'$n_B/n_0$'); ax.set_ylabel(r'$m^*/m_N$')
            figs['meff'][0].colorbar(plt.cm.ScalarMappable(cmap=global_cmap, norm=norm), ax=ax)
            figs['meff'][0].savefig(output_root / f"meff_gradient_{model}_{b_str}.pdf")

            # Potencial Químico
            ax = figs['mun'][1]
            if mun_data_all:
                ax.set_ylim(np.min(mun_data_all)*0.95, np.max(mun_data_all)*1.05)
            ax.set_xlabel(r'$n_B/n_0$'); ax.set_ylabel(r'$\mu_n$ [MeV]')
            figs['mun'][0].colorbar(plt.cm.ScalarMappable(cmap=global_cmap, norm=norm), ax=ax)
            figs['mun'][0].savefig(output_root / f"mun_gradient_{model}_{b_str}.pdf")

            plt.close('all')

    print(f"\n✅ Atlas físico atualizado (incluindo painel triplo log) salvo em: {output_root}")

if __name__ == "__main__":
    main()