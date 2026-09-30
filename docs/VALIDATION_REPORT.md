# Relatório de validação do NSRS

> Gerado automaticamente por `cargo run --release --bin nsrs -- report validation` em 2026-09-30 (revisão 86558d8), em 25 s. Não edite à mão: rode o comando de novo.

Legenda: ✅ dentro do critério · ⚠️ tensão ou desvio conhecido · ❌ fora do critério · ℹ️ informativo (sem critério).

Critério nas comparações com incerteza publicada: |Δ| ≤ 1σ ✅, ≤ 2σ ⚠️, > 2σ ❌.

## Resumo

| Seção | ✅ | ⚠️ | ❌ | ℹ️ |
|---|:-:|:-:|:-:|:-:|
| 1. Verificação numérica (soluções exatas e identidades) | 21 | 0 | 0 | 1 |
| 2. Parametrizações contra os artigos originais | 25 | 1 | 0 | 0 |
| 3. Relações universais e propriedades derivadas | 12 | 0 | 0 | 3 |
| 4. Vínculos observacionais e empíricos | 37 | 26 | 15 | 0 |
| 5. Campo magnético: perfis e acoplamento | 3 | 0 | 0 | 1 |
| 6. Eletrodinâmica não linear (NLEM) | 10 | 1 | 0 | 2 |
| **Total** | **108** | **28** | **15** | **7** |

Falhas de implementação (Seções 1–3, 5 e 6): **0**. Os ⚠️ e ❌ da Seção 4 (26 e 15) medem os **modelos** contra observações, não o código.

## 1. Verificação numérica (soluções exatas e identidades)

Estas verificações não dependem de dados empíricos: uma falha indica erro de implementação, não de modelo físico.


### 1.1 Integrador TOV

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| Raio R, densidade uniforme (12 estrelas, P_c/ε de 0.01 a 2) | erro rel. máx. 4.0e-11 | solução exata de Schwarzschild [1] | < 1.0e-8 | ✅ |
| Massa M, densidade uniforme (12 estrelas, P_c/ε de 0.01 a 2) | erro rel. máx. 1.2e-10 | solução exata de Schwarzschild [1] | < 1.0e-8 | ✅ |
| Massa própria, densidade uniforme (12 estrelas, P_c/ε de 0.01 a 2) | erro rel. máx. 1.2e-10 | solução exata de Schwarzschild [1] | < 1.0e-8 | ✅ |


### 1.2 Número de Love e momento de inércia no limite newtoniano

k₂ pela equação de Hinderer [2], com a correção de descontinuidade de densidade na superfície de Damour & Nagar [3]; I pela aproximação de rotação lenta de Hartle [4]. Em compacidade C → 0 ambos devem recuperar os valores newtonianos analíticos.

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| k₂, polítropo n = 1 (C ≈ 2×10⁻⁴) | 0.25981 | 0.25991 (15 − π²)/2π² (newtoniano) | erro rel. < 5.0e-3 | ✅ |
| I/MR², polítropo n = 1 | 0.26144 | 0.26138 (2/3)(1 − 6/π²) (newtoniano) | erro rel. < 5.0e-3 | ✅ |
| k₂, densidade uniforme (C ≈ 2×10⁻⁵) | 0.75000 | 0.75000 3/4 (newtoniano, fluido incompressível) | erro rel. < 1.0e-3 | ✅ |
| I/MR², densidade uniforme | 0.40001 | 0.40000 2/5 (newtoniano) | erro rel. < 1.0e-3 | ✅ |


### 1.3 Consistência termodinâmica e neutralidade de carga

Gibbs–Duhem a T = 0: dP/dμ_n = n_B, com P = P∥ (a pressão termodinâmica, −Ω), derivada por diferença central na malha (passo ≈ 1.4 MeV; erro de truncamento esperado O(10⁻⁴), maior nos limiares de partículas e níveis de Landau). Neutralidade: |Σ q_i n_i| < 10⁻⁷ + 10⁻⁶ n_B fm⁻³.

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| Gibbs–Duhem, GM1, B = 0 (531 pontos com n_B ≥ 0.5 n₀) | p95 3.7e-5, máx. 9.0e-4 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, GM1, B = 0 | máx. \|q\|/tolerância = 0.00 | 0 | < 1 | ✅ |
| Gibbs–Duhem, GM1, B = 10¹⁷ G (533 pontos com n_B ≥ 0.5 n₀) | p95 1.4e-4, máx. 1.5e-3 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, GM1, B = 10¹⁷ G | máx. \|q\|/tolerância = 0.01 | 0 | < 1 | ✅ |
| Gibbs–Duhem, GM3, B = 0 (531 pontos com n_B ≥ 0.5 n₀) | p95 2.8e-5, máx. 7.0e-4 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, GM3, B = 0 | máx. \|q\|/tolerância = 0.00 | 0 | < 1 | ✅ |
| Gibbs–Duhem, GM3, B = 10¹⁷ G (532 pontos com n_B ≥ 0.5 n₀) | p95 1.3e-4, máx. 1.1e-3 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, GM3, B = 10¹⁷ G | máx. \|q\|/tolerância = 0.01 | 0 | < 1 | ✅ |
| Gibbs–Duhem, FSU2, B = 0 (301 pontos com n_B ≥ 0.5 n₀) | p95 4.7e-5, máx. 3.5e-4 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, FSU2, B = 0 | máx. \|q\|/tolerância = 0.00 | 0 | < 1 | ✅ |
| Gibbs–Duhem, FSU2, B = 10¹⁷ G (301 pontos com n_B ≥ 0.5 n₀) | p95 2.3e-4, máx. 1.0e-3 | identidade exata | p95 < 5e-4, máx. < 5e-3 | ✅ |
| Neutralidade de carga, FSU2, B = 10¹⁷ G | máx. \|q\|/tolerância = 0.01 | 0 | < 1 | ✅ |


### 1.4 Níveis de Landau e magnetização

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| Limite B → 0 da soma sobre níveis de Landau (GM1, B = 10¹⁵ G, 532 pontos) | dif. rel. máx. em P e ε: 9.8e-7 | gás de Fermi isotrópico | < 1e-5 | ✅ |
| Magnetização 𝓜B = B ∂P∥/∂B exportada (GM1, B = 3×10¹⁷ G, μ_n = 1.25 M_N) | dif. rel. 7.3e-7 | diferença finita com passo 10⁻⁴ (interno: 10⁻⁵); P⊥ = P∥ − 𝓜B [5] [6] | < 1e-5 | ✅ |
| Oscilações de de Haas–van Alphen: 𝓜B com passo 10⁻³ vs exportado (μ_n = 1.40 M_N, n_B ≈ 4 n₀) | dif. rel. 3.4e-2 | P(B) oscila quando níveis de Landau cruzam a superfície de Fermi | informativo | ℹ️ |


## 2. Parametrizações contra os artigos originais

Propriedades da matéria nuclear simétrica na saturação, calculadas pelo mesmo motor usado nas estrelas, e massas máximas de estrelas só com núcleons (npeμ, B = 0, crosta BPS). Para GM1 e GM3 os artigos dão valores arredondados, sem incerteza; o critério reflete esse arredondamento. Para FSU2 o critério é a incerteza (1σ) publicada.


### 2.1 Saturação

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| GM1: n₀ [fm⁻³] | 0.1530 | 0.153 [7] | \|Δ\| = 0.0000 ≤ 0.002 | ✅ |
| GM1: E/A [MeV] | -16.31 | -16.3 [7] | \|Δ\| = 0.01 ≤ 0.1 | ✅ |
| GM1: K [MeV] | 299.7 | 300 [7] | \|Δ\| = 0.3 ≤ 3 | ✅ |
| GM1: M*/M | 0.700 | 0.7 [7] | \|Δ\| = 0.000 ≤ 0.005 | ✅ |
| GM1: J [MeV] | 32.48 | 32.5 [7] | \|Δ\| = 0.02 ≤ 0.3 | ✅ |
| GM1: J [MeV] | 32.48 | 32.52 [8] | \|Δ\| = 0.04 ≤ 0.1 | ✅ |
| GM1: L [MeV] | 93.88 | 94.04 [8] | \|Δ\| = 0.16 ≤ 0.5 | ✅ |
| GM1: K [MeV] | 299.7 | 300.5 [8] | \|Δ\| = 0.8 ≤ 1.5 | ✅ |
| GM3: n₀ [fm⁻³] | 0.1529 | 0.153 [7] | \|Δ\| = 0.0001 ≤ 0.002 | ✅ |
| GM3: E/A [MeV] | -16.30 | -16.3 [7] | \|Δ\| = 0.00 ≤ 0.1 | ✅ |
| GM3: K [MeV] | 239.8 | 240 [7] | \|Δ\| = 0.2 ≤ 2.4 | ✅ |
| GM3: M*/M | 0.780 | 0.78 [7] | \|Δ\| = 0.000 ≤ 0.005 | ✅ |
| GM3: J [MeV] | 32.47 | 32.5 [7] | \|Δ\| = 0.03 ≤ 0.3 | ✅ |
| GM3: J [MeV] | 32.47 | 32.51 [8] | \|Δ\| = 0.04 ≤ 0.1 | ✅ |
| GM3: L [MeV] | 89.63 | 89.75 [8] | \|Δ\| = 0.12 ≤ 0.5 | ✅ |
| GM3: K [MeV] | 239.8 | 240.04 [8] | \|Δ\| = 0.3 ≤ 1.5 | ✅ |
| FSU2: n₀ [fm⁻³] | 0.1503 | 0.1505 ± 0.0007 [9] | \|Δ\| = 0.22σ | ✅ |
| FSU2: E/A [MeV] | -16.26 | -16.28 ± 0.02 [9] | \|Δ\| = 0.88σ | ✅ |
| FSU2: M*/M | 0.593 | 0.593 ± 0.004 [9] | \|Δ\| = 0.04σ | ✅ |
| FSU2: K [MeV] | 237.5 | 238 ± 2.8 [9] | \|Δ\| = 0.17σ | ✅ |
| FSU2: J [MeV] | 37.56 | 37.62 ± 1.11 [9] | \|Δ\| = 0.05σ | ✅ |
| FSU2: L [MeV] | 112.6 | 112.8 ± 16.1 [9] | \|Δ\| = 0.01σ | ✅ |


### 2.2 Estrelas só com núcleons

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| GM1: M_max [M☉] | 2.359 | 2.363 [8] | \|Δ\| < 0.01 M☉ (M_N e crosta) | ✅ |
| GM3: M_max [M☉] | 2.015 | 2.018 [8] | \|Δ\| < 0.01 M☉ (M_N e crosta) | ✅ |
| FSU2: M_max [M☉] | 2.071 | 2.07 ± 0.02 [9] | \|Δ\| = 0.05σ | ✅ |
| FSU2: R₁.₄ [km] | 13.95 | 14.42 ± 0.26 [9] | \|Δ\| = 1.8σ | ⚠️ |

O raio de FSU2 difere do artigo porque o NSRS junta a tabela BPS [10] diretamente ao núcleo, enquanto Chen & Piekarewicz [9] interpolam a crosta interna com um polítropo. A massa máxima, dominada pelo núcleo, não é afetada. Desvio conhecido; uma crosta unificada resolveria.


## 3. Relações universais e propriedades derivadas

Λ e k₂ vêm da equação de Hinderer [2] [12]; Ī = I/M³ da aproximação de rotação lenta de Hartle [4]. A relação universal I-Love [11] (Tabela I: ln Ī = 1.47 + 0.0817x + 0.0149x² + 2.87×10⁻⁴x³ − 3.64×10⁻⁵x⁴, x = ln Λ) tem precisão declarada < 1% e é independente da EoS; ela testa, em conjunto, o integrador de maré e o de inércia.

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| I-Love, GM1 com hyperons (379 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.84% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |
| I-Love, GM3 com hyperons (332 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.74% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |
| I-Love, FSU2 com hyperons (205 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.78% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |
| I-Love, GM1 só núcleons (584 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.84% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |
| I-Love, GM3 só núcleons (515 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.74% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |
| I-Love, FSU2 só núcleons (384 estrelas, 1 M☉ ≤ M ≤ M_max) | desvio máx. 0.78% | ajuste universal [11] | < 1% (declarado); ⚠️ até 1.5% | ✅ |


### 3.1 URCA direto, causalidade e propriedades de 1.4 M☉

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| GM1: Y_p no limiar do URCA direto | 13.25% (n_B = 1.80 n₀) | 11.1% (npe) a 14.8% (npeμ) [13] | dentro da faixa ± 0.5% | ✅ |
| GM1: c_s² máx. até o centro de M_max | 0.440 | causalidade | c_s² ≤ 1 e Γ > 0 | ✅ |
| GM1: estrela de 1.4 M☉ (com hyperons) | R = 13.71 km, k₂ = 0.1040, Λ = 890, I = 1.827×10⁴⁵ g cm², z = 0.197 | ver Seção 4 (observações) | informativo | ℹ️ |
| GM3: Y_p no limiar do URCA direto | 13.37% (n_B = 1.91 n₀) | 11.1% (npe) a 14.8% (npeμ) [13] | dentro da faixa ± 0.5% | ✅ |
| GM3: c_s² máx. até o centro de M_max | 0.400 | causalidade | c_s² ≤ 1 e Γ > 0 | ✅ |
| GM3: estrela de 1.4 M☉ (com hyperons) | R = 13.03 km, k₂ = 0.0860, Λ = 571, I = 1.616×10⁴⁵ g cm², z = 0.210 | ver Seção 4 (observações) | informativo | ℹ️ |
| FSU2: Y_p no limiar do URCA direto | 13.00% (n_B = 1.40 n₀) | 11.1% (npe) a 14.8% (npeμ) [13] | dentro da faixa ± 0.5% | ✅ |
| FSU2: c_s² máx. até o centro de M_max | 0.230 | causalidade | c_s² ≤ 1 e Γ > 0 | ✅ |
| FSU2: estrela de 1.4 M☉ (com hyperons) | R = 13.70 km, k₂ = 0.0865, Λ = 738, I = 1.731×10⁴⁵ g cm², z = 0.197 | ver Seção 4 (observações) | informativo | ℹ️ |


## 4. Vínculos observacionais e empíricos

Vínculos de `input/observations/constraints.csv`. Distância d em desvios-padrão, com barras assimétricas (intervalos de 90% convertidos para 1σ dividindo por 1.645): d ≤ 1 compatível ✅, 1 < d ≤ 2 tensão ⚠️, d > 2 excluído ❌. Massa máxima: só conta se o modelo não a alcança. Pontos M-R: menor distância ao ramo estável. Estrelas com B = 0 e crosta BPS; "H" = com hyperons, "N" = só núcleons.

| Vínculo | Observado | GM1 H | GM1 N | GM3 H | GM3 N | FSU2 H | FSU2 N |
|---|---|:-:|:-:|:-:|:-:|:-:|:-:|
| PSR J0348+0432 [14] | 2.01 (−0.04, +0.04), 68% | ✅ 1.994 (d = 0.4) | ✅ 2.359 (d = 0.0) | ❌ 1.700 (d = 7.8) | ✅ 2.015 (d = 0.0) | ❌ 1.598 (d = 10.3) | ✅ 2.071 (d = 0.0) |
| PSR J0740+6620 (timing) [15] | 2.08 (−0.07, +0.07), 68.3% | ⚠️ 1.994 (d = 1.2) | ✅ 2.359 (d = 0.0) | ❌ 1.700 (d = 5.4) | ✅ 2.015 (d = 0.9) | ❌ 1.598 (d = 6.9) | ✅ 2.071 (d = 0.1) |
| PSR J0030+0451 (Riley+2019) [16] | 12.71 (−1.19, +1.14), 68% | ✅ 13.718 (d = 0.9) | ✅ 13.721 (d = 0.9) | ✅ 13.117 (d = 0.4) | ✅ 13.173 (d = 0.4) | ✅ 13.734 (d = 1.0) | ⚠️ 13.981 (d = 1.1) |
| PSR J0030+0451 (Miller+2019) [17] | 13.02 (−1.06, +1.24), 68% | ✅ 13.696 (d = 0.5) | ✅ 13.717 (d = 0.6) | ✅ 12.952 (d = 0.1) | ✅ 13.073 (d = 0.0) | ✅ 13.459 (d = 0.4) | ✅ 13.905 (d = 0.7) |
| PSR J0740+6620 (Riley+2021) [18] | 12.39 (−0.98, +1.3), 68% | ⚠️ 11.942 (d = 1.3) | ✅ 13.270 (d = 0.7) | ❌ 11.085 (d = 5.8) | ⚠️ 11.325 (d = 1.5) | ❌ 12.141 (d = 7.2) | ✅ 12.321 (d = 0.1) |
| PSR J0740+6620 (Miller+2021) [19] | 13.7 (−1.5, +2.6), 68% | ⚠️ 12.032 (d = 1.7) | ✅ 13.267 (d = 0.3) | ❌ 11.078 (d = 5.7) | ⚠️ 11.317 (d = 1.9) | ❌ 12.172 (d = 7.0) | ✅ 12.514 (d = 0.9) |
| R(1.4 Msun) combined [19] | 12.45 (−0.65, +0.65), 68% (range of frameworks) | ⚠️ 13.709 (d = 1.9) | ⚠️ 13.720 (d = 2.0) | ✅ 13.030 (d = 0.9) | ⚠️ 13.121 (d = 1.0) | ⚠️ 13.700 (d = 1.9) | ❌ 13.948 (d = 2.3) |
| GW170817 Lambda(1.4 Msun) [20] | 190 (−120, +390), 90% | ❌ 889.833 (d = 3.0) | ❌ 896.595 (d = 3.0) | ⚠️ 570.768 (d = 1.6) | ⚠️ 608.038 (d = 1.8) | ❌ 737.536 (d = 2.3) | ❌ 866.075 (d = 2.9) |
| saturation density [fm^-3] [21] | 0.155 (−0.005, +0.005), 1 sigma | ✅ 0.153 (d = 0.4) | ✅ 0.153 (d = 0.4) | ✅ 0.153 (d = 0.4) | ✅ 0.153 (d = 0.4) | ✅ 0.150 (d = 0.9) | ✅ 0.150 (d = 0.9) |
| binding energy E/A [MeV] [21] | -15.8 (−0.3, +0.3), 1 sigma | ⚠️ -16.310 (d = 1.7) | ⚠️ -16.310 (d = 1.7) | ⚠️ -16.304 (d = 1.7) | ⚠️ -16.304 (d = 1.7) | ⚠️ -16.262 (d = 1.5) | ⚠️ -16.262 (d = 1.5) |
| incompressibility K [MeV] [21] | 230 (−20, +20), 1 sigma | ❌ 299.725 (d = 3.5) | ❌ 299.725 (d = 3.5) | ✅ 239.762 (d = 0.5) | ✅ 239.762 (d = 0.5) | ✅ 237.536 (d = 0.4) | ✅ 237.536 (d = 0.4) |
| symmetry energy J [MeV] [22] | 31.7 (−3.2, +3.2), 1 sigma | ✅ 32.481 (d = 0.2) | ✅ 32.481 (d = 0.2) | ✅ 32.473 (d = 0.2) | ✅ 32.473 (d = 0.2) | ⚠️ 37.561 (d = 1.8) | ⚠️ 37.561 (d = 1.8) |
| symmetry slope L [MeV] [22] | 58.7 (−28.1, +28.1), 1 sigma | ⚠️ 93.876 (d = 1.3) | ⚠️ 93.876 (d = 1.3) | ⚠️ 89.626 (d = 1.1) | ⚠️ 89.626 (d = 1.1) | ⚠️ 112.639 (d = 1.9) | ⚠️ 112.639 (d = 1.9) |

¹ Não comparável: a massa pedida está acima da massa máxima do modelo (a exclusão já aparece no vínculo de massa máxima).

Este confronto avalia os **modelos**, não o código: tensões e exclusões são resultados físicos (p.ex. GM1 e FSU2 rígidos demais para GW170817; hyperons reduzindo M_max abaixo de 2 M☉).


## 5. Campo magnético: perfis e acoplamento

Perfis do campo local: BDD [23], B = B_surf + B₀[1 − exp(−β(n_B/n₀)^γ)] com β = 0.01, γ = 3, B_surf = 10¹⁵ G; e o ajuste polar de Dexheimer et al. [24], Eq. (1), B = (a + bμ_B + cμ_B²)μ/B_c, com os coeficientes da Tabela 2.

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| Dexheimer, Eq. (1): B(μ_B = 1000 MeV), μ = 3×10³² A m², M_B = 2.2 M☉ | 5.7771e17 G | 5.777×10¹⁷ G (coeficientes da Tabela 2) [24] | dif. rel. < 1e-3 | ✅ |
| Dexheimer: dP/dμ = n_B + 𝓜 dB/dμ (GM1, n_B = 1, 2, 4, 6 n₀) | dif. rel. máx. 1.6e-7 | identidade termodinâmica com B = B(μ_B) | < 1e-6 | ✅ |
| BDD: campo usado = B(n_B da própria solução) (GM1, B₀ = 10¹⁸ G) | dif. rel. máx. 2.7e-10 | perfil [23] | < 1e-8 | ✅ |
| Dexheimer: efeito do campo na matéria sobre M_max (GM1, sem tensões do campo) | 1.9908 vs 1.9940 M☉ (-0.16%) | campo tratado na estrutura por Einstein–Maxwell | informativo | ℹ️ |

Com o perfil de Dexheimer a energia e as tensões do campo ficam fora da EoS da TOV por padrão: o ajuste vem de soluções de Einstein–Maxwell, nas quais o campo já está na estrutura, e é destinado à EoS microscópica.


## 6. Eletrodinâmica não linear (NLEM)

A NLEM altera só a energia e as tensões do próprio campo: P∥ = −ε_B e P⊥ = HB − ε_B, com H = dε_B/dB [25]. No modelo logarítmico ε_B = ξ² ln(1 + x), x = B²/2ξ². Não há, até onde sabemos, resultados publicados de estrelas de nêutrons com a NLEM logarítmica para comparar; as verificações abaixo são de consistência interna, de limites analíticos e de sistemáticos.


### 6.1 Tensor de tensões e acoplamento

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| P⊥ = HB − ε_B, Log(ξ = 2×10¹⁷ G), B de 10¹⁶ a 3×10¹⁸ G | dif. rel. máx. 3.1e-10 | H = dε_B/dB numérico [25] | < 1e-7 | ✅ |
| Log com ξ ≫ B recupera Maxwell (B = 10¹⁵ G, ξ = 10²⁵ G) | dif. rel. 0.0e0 | limite x → 0 | < 1e-14 | ✅ |
| Matéria idêntica a Maxwell no mesmo B (GM1, 3×10¹⁷ G; ModMax(1), Log(10¹⁶), Log(10¹⁸)) | dif. rel. máx. em n_B e ε_matéria: 1.1e-16 | acoplamento mínimo: Landau usa B | < 1e-12 | ✅ |


### 6.2 Log(ξ): limites, sinal da pressão do campo e sistemáticos

| Grandeza | NSRS | Referência | Critério | |
|---|---|---|---|:-:|
| x* (P⊥ do campo < 0 para x > x*), anisotrópica | 3.92155 | raiz de 2x = (1+x) ln(1+x): 3.92155 | dif. rel. < 1e-6 | ✅ |
| x* ((P∥+2P⊥)/3 do campo < 0), isotrópica | 0.83283 | raiz de 4x = 3(1+x) ln(1+x): 0.83283 | dif. rel. < 1e-6 | ✅ |
| EoS vazia para ξ < B_surf/√(2x*) = 3.571e14 G (anisotrópica) | 0.9ξ_min: vazia; 1.1ξ_min: 1.9903 | pressão do campo já negativa na superfície (B_surf = 10¹⁵ G) | vazia abaixo, válida acima | ✅ |
| Limite ξ → ∞ (ξ = 10²⁰ G) = Maxwell, anisotrópica | 2.0156 vs 2.0156 M☉ | limite analítico | \|ΔM_max\| < 2×10⁻⁴ M☉ | ✅ |
| Limite ξ → 0 (ξ = 10¹⁵ G) = sem tensão do campo, anisotrópica | 1.9903 vs 1.9903 M☉ | ε_B = ξ² ln(1+x) → 0 | \|ΔM_max\| < 2×10⁻⁴ M☉ | ✅ |
| M_max mínimo na varredura ξ = 10¹⁶–10¹⁹ G, anisotrópica | 1.9861 M☉ em ξ = 1.78e17 G | mínimo entre Maxwell e sem tensão: 1.9903 | resultado (janela com pressão do campo < 0) | ℹ️ |
| Limite ξ → ∞ (ξ = 10²⁰ G) = Maxwell, isotrópica | 1.9907 vs 1.9907 M☉ | limite analítico | \|ΔM_max\| < 2×10⁻⁴ M☉ | ✅ |
| Limite ξ → 0 (ξ = 10¹⁵ G) = sem tensão do campo, isotrópica | 1.9924 vs 1.9925 M☉ | ε_B = ξ² ln(1+x) → 0 | \|ΔM_max\| < 2×10⁻⁴ M☉ | ✅ |
| M_max mínimo na varredura ξ = 10¹⁶–10¹⁹ G, isotrópica | 1.9811 M☉ em ξ = 3.16e17 G | mínimo entre Maxwell e sem tensão: 1.9907 | resultado (janela com pressão do campo < 0) | ℹ️ |
| Sistemático de geometria: M_max(P⊥) − M_max(média isotrópica), Maxwell, B₀ = 10¹⁸ G | +0.0250 M☉ | P⊥ numa TOV esférica é inconsistente [26] | incerteza sistemática | ⚠️ |

GM1 com hyperons, perfil `Constant` (campo da energia pelo perfil BDD com B_surf = 10¹⁵ G). A média isotrópica (P∥ + 2P⊥)/3 é a média angular das tensões de um campo de direção fixa e é consistente com a simetria esférica; usar P⊥ em todas as direções não corresponde a nenhuma geometria esférica consistente [26]. Os valores absolutos de ΔM na topologia anisotrópica devem ser lidos como cota heurística.


## 7. Limitações conhecidas

- **Crosta.** A tabela BPS [10] é unida diretamente ao núcleo; o raio de FSU2 fica ~0.5 km abaixo do de Chen & Piekarewicz [9]. Uma crosta unificada resolveria.
- **Malha em μ_n.** O padrão (≤ 1.8 M_N) é curto para EoS nucleônicas rígidas; este relatório usa μ_n ≤ 3 M_N nas sequências estelares.
- **Campo constante forte.** Com o perfil `Constant` e B ≳ 10¹⁸ G há uma transição de primeira ordem na entrada da matéria; a continuação em μ_n pode atravessá-la ou parar, conforme a plataforma. Uma construção de Maxwell tornaria o resultado único.
- **Tensões do campo na TOV.** A TOV é esférica; a topologia anisotrópica (P⊥ em todas as direções) é heurística [26]. Resultados quantitativos em B ≳ 10¹⁸ G exigem Einstein–Maxwell axissimétrico.
- **Maré e inércia com campo.** Λ e I usam as equações de perturbação isotrópicas mesmo quando a EoS inclui tensões anisotrópicas.
- **Perfil BDD.** β e γ são parâmetros fenomenológicos; o perfil não resolve as equações de Maxwell.
- **NLEM logarítmica.** Sem resultados publicados de estrelas para comparar; ver Seção 6.


## Referências

1. K. Schwarzschild, *Sitzungsber. Preuss. Akad. Wiss. Berlin*, 424 (1916) — solução interior de densidade uniforme.
2. T. Hinderer, *ApJ* **677**, 1216 (2008). arXiv:0711.2420
3. T. Damour & A. Nagar, *Phys. Rev. D* **80**, 084035 (2009). arXiv:0906.0096
4. J. B. Hartle, *ApJ* **150**, 1005 (1967).
5. E. J. Ferrer *et al.*, *Phys. Rev. C* **82**, 065802 (2010).
6. M. Strickland, V. Dexheimer & D. P. Menezes, *Phys. Rev. D* **86**, 125032 (2012).
7. N. K. Glendenning & S. A. Moszkowski, *Phys. Rev. Lett.* **67**, 2414 (1991). doi:10.1103/PhysRevLett.67.2414
8. G. Nam, Y. Lim & J. W. Holt, *Universal Relation for the Neutron Star Maximum Mass within Relativistic Mean-Field Theories*, arXiv:2510.15356 (2025), Tabela III.
9. W.-C. Chen & J. Piekarewicz, *Phys. Rev. C* **90**, 044305 (2014). doi:10.1103/PhysRevC.90.044305, arXiv:1408.4159
10. G. Baym, C. Pethick & P. Sutherland, *ApJ* **170**, 299 (1971) — crosta BPS.
11. K. Yagi & N. Yunes, *Science* **341**, 365 (2013). doi:10.1126/science.1236462, arXiv:1302.4499
12. S. Postnikov, M. Prakash & J. M. Lattimer, *Phys. Rev. D* **82**, 024016 (2010). arXiv:1004.5098
13. J. M. Lattimer, C. J. Pethick, M. Prakash & P. Haensel, *Phys. Rev. Lett.* **66**, 2701 (1991). doi:10.1103/PhysRevLett.66.2701
14. Antoniadis et al. 2013 Science 340 1233232. doi:10.1126/science.1233232, arXiv:1304.6875
15. Fonseca et al. 2021 ApJL. doi:10.3847/2041-8213/ac03b8, arXiv:2104.00880
16. Riley et al. 2019 ApJL 887 L21. doi:10.3847/2041-8213/ab481c, arXiv:1912.05702
17. Miller et al. 2019 ApJL 887 L24. doi:10.3847/2041-8213/ab50c5, arXiv:1912.05705
18. Riley et al. 2021 ApJL 918. doi:10.3847/2041-8213/ac0a81, arXiv:2105.06980
19. Miller et al. 2021 ApJL. doi:10.3847/2041-8213/ac089b, arXiv:2105.06979
20. Abbott et al. (LVC) 2018 PRL 121 161101. doi:10.1103/PhysRevLett.121.161101, arXiv:1805.11581
21. Margueron Hoffmann Casali & Gulminelli 2018 PRC 97 025805. doi:10.1103/PhysRevC.97.025805, arXiv:1708.06894
22. Oertel et al. 2017 RMP 89 015007. doi:10.1103/RevModPhys.89.015007, arXiv:1610.03361
23. D. Bandyopadhyay, S. Chakrabarty & S. Pal, *Phys. Rev. Lett.* **79**, 2176 (1997). arXiv:astro-ph/9703066
24. V. Dexheimer, B. Franzon, R. O. Gomes, R. L. S. Farias, S. S. Avancini & S. Schramm, *Phys. Lett. B* **773**, 487 (2017). arXiv:1612.05795
25. H. H. Soleng, *Phys. Rev. D* **52**, 6178 (1995). arXiv:hep-th/9509033
26. D. Chatterjee, T. Elghozi, J. Novak & M. Oertel, *MNRAS* **447**, 3785 (2015).
