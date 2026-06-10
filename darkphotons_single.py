import matplotlib.pyplot as plt
import numpy as np

# Listas para armazenar Raio e Massa de cada modelo
radii = [[], []]
masses = [[], []]

base_dir = "output/darkphotons/GM1"
eos_filenames = ["darkphoton.dat", "eos_hadrons.dat"]
labels = ["Com Fóton Escuro (DM)", "Hadrônica Pura"]
colors = ["purple", "black"]
linestyles = ["-", "--"]

# Processamento e leitura dos arquivos
for i, n in enumerate(eos_filenames):
    with open(f"{base_dir}/{n}", "r") as f:
        for line in f:
            # Pula cabeçalhos ou linhas vazias
            if line.startswith("#") or not line.strip():
                continue
            
            data = line.split()
            
            try:
                r_val = float(data[22]) 
                m_val = float(data[21]) 
                if m_val < 2.1 and r_val < 15:
                    radii[i].append(r_val)
                    masses[i].append(m_val)
            except ValueError:
                # Ignora linhas que não puderam ser convertidas para float
                continue

# Configuração do Gráfico
plt.figure(figsize=(8, 6))

# Plotando as EoS e adicionando marcadores nos pontos máximos
for i in range(len(eos_filenames)):
    # Plota a curva normal
    plt.plot(radii[i], masses[i], label=labels[i], 
             color=colors[i], linestyle=linestyles[i], linewidth=2.5)
    
    # Se a lista não estiver vazia, encontra o ponto de massa máxima
    if masses[i]:
        max_mass_idx = np.argmax(masses[i]) # Índice da maior massa
        max_m = masses[i][max_mass_idx]
        corresponding_r = radii[i][max_mass_idx]
        
        # Plota um ponto (bolinha) exatamente no topo da curva
        plt.scatter(corresponding_r, max_m, color=colors[i], s=50, zorder=5)
        
        # Define deslocamentos para o texto não ficar em cima da linha
        # (Ajuste o xytext se o texto cobrir as curvas)
        if i == 0:
            offset_x, offset_y = -25, 35
        else:
            offset_x, offset_y = 50, 35
        
        # Adiciona a anotação com seta
        plt.annotate(
            f"$M_{{max}} = {max_m:.4f} M_\\odot$\n$R = {corresponding_r:.4f}$ km",
            xy=(corresponding_r, max_m),
            xytext=(offset_x, offset_y),
            textcoords="offset points",
            fontsize=10,
            fontweight='bold',
            color=colors[i],
            arrowprops=dict(
                arrowstyle="->",
                color=colors[i],
                lw=1.5,
                connectionstyle="arc3,rad=.2"  # Seta levemente curvada
            )
        )

# Adicionando o limite observacional do PSR J0740+6620 (~2.08 M_sol)
plt.axhline(y=2.08, color='gray', linestyle=':', label='PSR J0740+6620')
plt.axhspan(2.01, 2.15, color='gray', alpha=0.2)

# Embelezamento do gráfico
plt.xlabel(r'Raio $R$ (km)', fontsize=14)
plt.ylabel(r'Massa $M$ ($M_\odot$)', fontsize=14)
plt.title('Relação Massa-Raio: Portal Vetorial (Matéria Escura)', fontsize=16)

plt.xlim(9, 15) 
plt.ylim(0.5, 2.6)

plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(loc='lower left', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# Salva e exibe o gráfico
plt.savefig(f"{base_dir}/mr_comparativo.png", dpi=300)
# plt.show()


