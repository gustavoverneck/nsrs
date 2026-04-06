import os
import pickle
import numpy as np
import pandas as pd
import emcee
import corner
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import LinearNDInterpolator

# ==============================================================================
# 1. CONFIGURAÇÕES GERAIS E DADOS OBSERVACIONAIS
# ==============================================================================
INPUT_ROOT = 'output/limits'      # Raiz onde o seu código Rust salva os dados
MODELOS = ['GM1', 'GM3']          # Modelos a serem comparados

# Dados de várias estrelas (Massa, ErroM, Raio, ErroR)
OBSERVACOES = [
    (2.08, 0.07, 12.35, 0.75),  # PSR J0740+6620 (Estrela Pesada)
    (1.44, 0.15, 13.02, 1.24)   # PSR J0030+0451 (Estrela Canónica)
]

# Cores para o gráfico (padrão de publicações científicas)
CORES_MODELOS = {'GM1': 'darkblue', 'GM3': 'darkred'}

# ==============================================================================
# 2. CARREGAMENTO DE DADOS (CRAWLER DOS DIRETÓRIOS)
# ==============================================================================
def build_dataframe_from_dirs(input_root_str, target_model):
    print(f"Buscando dados na árvore de diretórios: {input_root_str}/{target_model}/...")
    input_root = Path(input_root_str)
    data_rows = []
    
    # Varre todos os eos.dat
    for eos_path in input_root.rglob('eos.dat'):
        parts = eos_path.parts
        try:
            # Estrutura baseada no seu script: output/limits/MODELO/CAMPO/CSI/eos.dat
            model = parts
            if model != target_model:
                continue
            
            # Limpeza das strings para extrair os floats
            b_val = float(parts.replace('B_', ''))
            csi_val = float(parts.replace('csi_', ''))
            
            # Buscar mr.dat na mesma pasta
            mr_path = eos_path.parent / 'mr.dat'
            
            if mr_path.exists():
                mr_data = np.loadtxt(mr_path) # Assumindo 0: R(km), 1: M(M_sun)
                if len(mr_data) > 0:
                    radii = mr_data[:, 0]
                    masses = mr_data[:, 1]
                    
                    # Estrela de Massa Máxima
                    max_idx = np.argmax(masses)
                    data_rows.append({
                        'b_field': b_val,
                        'csi': csi_val,
                        'max_mass': masses[max_idx],
                        'radius_at_max': radii[max_idx]
                    })
        except Exception as e:
            continue
            
    if not data_rows:
        raise ValueError(f"❌ Nenhum dado válido encontrado para {target_model}.")
        
    df = pd.DataFrame(data_rows)
    print(f"[{target_model}] Sucesso! {len(df)} simulações carregadas.")
    return df

# ==============================================================================
# 3. TREINAMENTO DOS MODELOS SURROGATOS (INTERPOLADORES)
# ==============================================================================
def train_surrogate_model(input_root, target_model):
    model_save_path = f"surrogate_{target_model.lower()}.pkl"
    
    if os.path.exists(model_save_path):
        print(f"[{target_model}] Lendo modelo surrogato compilado ({model_save_path})...")
        with open(model_save_path, 'rb') as f:
            return pickle.load(f)

    df = build_dataframe_from_dirs(input_root, target_model)
    print(f"[{target_model}] Treinando a malha matemática de interpolação...")
    
    df['log10_B'] = np.log10(df['b_field'])
    df['log10_csi'] = np.log10(df['csi'])
    df = df.drop_duplicates(subset=['log10_B', 'log10_csi'])
    
    X = df[['log10_B', 'log10_csi']].values
    interp_M = LinearNDInterpolator(X, df['max_mass'].values)
    interp_R = LinearNDInterpolator(X, df['radius_at_max'].values)
    
    bounds = {
        'min_B': df['log10_B'].min(), 'max_B': df['log10_B'].max(),
        'min_csi': df['log10_csi'].min(), 'max_csi': df['log10_csi'].max()
    }
    
    model = {'interp_M': interp_M, 'interp_R': interp_R, 'bounds': bounds, 'name': target_model}
    
    with open(model_save_path, 'wb') as f:
        pickle.dump(model, f)
    
    return model

# ==============================================================================
# 4. FUNÇÕES DA ESTATÍSTICA BAYESIANA
# ==============================================================================
def log_prior(theta, bounds):
    log10_B, log10_csi = theta
    if (bounds['min_B'] <= log10_B <= bounds['max_B'] and 
        bounds['min_csi'] <= log10_csi <= bounds['max_csi']):
        return 0.0
    return -np.inf

def log_likelihood(theta, surrogate):
    log10_B, log10_csi = theta
    m_pred = surrogate['interp_M'](log10_B, log10_csi)
    r_pred = surrogate['interp_R'](log10_B, log10_csi)
    
    if np.isnan(m_pred) or np.isnan(r_pred):
        return -np.inf
    
    log_l_total = 0.0
    
    # Soma o ajustamento estatístico para todas as estrelas medidas
    for obs_m, err_m, obs_r, err_r in OBSERVACOES:
        log_l_m = -0.5 * ((m_pred - obs_m) / err_m)**2
        log_l_r = -0.5 * ((r_pred - obs_r) / err_r)**2
        log_l_total += (log_l_m + log_l_r)
        
    return log_l_total

def log_probability(theta, surrogate):
    lp = log_prior(theta, surrogate['bounds'])
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, surrogate)

# ==============================================================================
# 5. MCMC RUNNER
# ==============================================================================
def run_mcmc(surrogate):
    print(f"\nIniciando MCMC para o modelo {surrogate['name']}...")
    bounds = surrogate['bounds']
    
    nwalkers, ndim, nsteps = 32, 2, 5000 
    
    initial_B = np.random.uniform(bounds['min_B'], bounds['max_B'], nwalkers)
    initial_csi = np.random.uniform(bounds['min_csi'], bounds['max_csi'], nwalkers)
    pos = np.column_stack((initial_B, initial_csi))
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(surrogate,))
    sampler.run_mcmc(pos, nsteps, progress=True)
    
    return sampler.get_chain(discard=1000, thin=15, flat=True)

# ==============================================================================
# 6. PLOTAGEM GERAL COMPARATIVA (CORNER PLOT)
# ==============================================================================
def plot_combined_posteriors(todas_samples):
    print("\nGerando Corner Plot partilhado...")
    labels = [r"$\log_{10}(B)$ [G]", r"$\log_{10}(\xi)$"]
    fig = None

    import matplotlib.lines as mlines
    legend_handles = []

    for modelo, samples in todas_samples.items():
        cor = CORES_MODELOS[modelo]
        
        # O primeiro modelo cria a figura, os seguintes desenham por cima
        fig = corner.corner(
            samples, 
            labels=labels, 
            quantiles=[0.16, 0.5, 0.84], 
            show_titles=False, # Desativado aqui para não sobrepor texto
            color=cor,
            fig=fig,
            hist_kwargs={'density': True, 'linewidth': 2},
            contour_kwargs={'linewidths': 1.5}
        )
        # Adiciona à legenda
        legend_handles.append(mlines.Line2D([], [], color=cor, label=modelo, linewidth=3))

        # Calcula os resultados para imprimir no terminal
        mcmc_B = np.percentile(samples[:, 0],)
        mcmc_csi = np.percentile(samples[:, 1],)
        print(f"--- Inferência para {modelo} ---")
        print(f"log10(B)   = {mcmc_B:.2f} (+{mcmc_B-mcmc_B:.2f} / -{mcmc_B-mcmc_B:.2f})")
        print(f"log10(csi) = {mcmc_csi:.2f} (+{mcmc_csi-mcmc_csi:.2f} / -{mcmc_csi-mcmc_csi:.2f})\n")

    # Título Geral e Legenda
    fig.suptitle(f"Seleção de Modelos (Obs: M={OBS_MASS}±{ERR_MASS}, R={OBS_RADIUS}±{ERR_RADIUS})", fontsize=15, y=1.02)
    fig.legend(handles=legend_handles, loc='upper right', fontsize=14, bbox_to_anchor=(0.95, 0.95))

    plt.savefig("corner_plot_comparativo.png", dpi=300, bbox_inches='tight')
    print("Gráfico salvo como 'corner_plot_comparativo.png'")

# ==============================================================================
# MAIN 
# ==============================================================================
if __name__ == "__main__":
    todas_samples = {}
    
    # 1. Corre o treino e a inferência para cada modelo separadamente
    for modelo in MODELOS:
        print(f"\n{'='*40}\nProcessando modelo: {modelo}\n{'='*40}")
        surrogate = train_surrogate_model(INPUT_ROOT, modelo)
        samples = run_mcmc(surrogate)
        todas_samples[modelo] = samples
        
    # 2. Desenha tudo no mesmo painel
    plot_combined_posteriors(todas_samples)