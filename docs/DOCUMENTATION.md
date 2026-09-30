# NSRS Documentation (Foco Hadrônico)

> Documentação técnica com foco no pipeline hadrônico: EoS RMF, campo magnético, topologia magnética, NLEM e integração TOV.

---

## Índice

- [NSRS Documentation (Foco Hadrônico)](#nsrs-documentation-foco-hadrônico)
  - [Índice](#índice)
  - [Escopo atual](#escopo-atual)
  - [Visão geral do pipeline hadrônico](#visão-geral-do-pipeline-hadrônico)
  - [Arquitetura do projeto](#arquitetura-do-projeto)
    - [`src/core/solver.rs`](#srccoresolverrs)
    - [`src/core/physics.rs`](#srccorephysicsrs)
    - [`src/core/tov_solver.rs`](#srccoretov_solverrs)
    - [`src/core/plotting.rs`](#srccoreplottingrs)
    - [`src/core/io_utils.rs`](#srccoreio_utilsrs)
  - [Fluxo numérico principal](#fluxo-numérico-principal)
  - [Formato dos dados da EoS hadrônica](#formato-dos-dados-da-eos-hadrônica)
  - [Executável `nsrs` (`src/bin/nsrs`)](#executável-nsrs-srcbinnsrs)
  - [Módulos principais (`src/core`)](#módulos-principais-srccore)
  - [Quarks e híbridas (em desenvolvimento)](#quarks-e-híbridas-em-desenvolvimento)
  - [Diretórios de saída](#diretórios-de-saída)
  - [Como executar (foco hadrônico)](#como-executar-foco-hadrônico)
  - [Dicas de estabilidade numérica](#dicas-de-estabilidade-numérica)

---

## Escopo atual

Este documento prioriza **estrelas hadrônicas**. O pipeline operacional principal hoje está centrado em:

- `HadronsMatter` (modelo RMF),
- campo magnético dependente da densidade,
- topologia magnética isotrópica/anisotrópica,
- NLEM (`Maxwell`, `Modmax`, `Log`),
- resolução da EoS e integração TOV para curvas M-R.

---

## Visão geral do pipeline hadrônico

Fluxo alto nível:

1. Configurar um motor hadrônico (`HadronsMatter`).
2. Resolver a EoS ponto a ponto em $\mu_n$ com `Solver`.
3. Extrair $\epsilon(P)$ e integrar TOV.
4. Gerar curvas M-R e gráficos auxiliares.
5. Exportar dados em `.dat`, `.csv`, `.svg`.

---

## Arquitetura do projeto

### `src/core/solver.rs`

- `EngineMode` define o modo físico (hadrons/quarks/hybrid).
- Para o escopo atual, o caminho principal é `EngineMode::Hadrons(HadronsMatter)`.
- `Solver::solve()` implementa varredura adaptativa em $\mu_n$.
- `Solver::solve_parallel()` acelera varreduras grandes em paralelo.
- `Solver::write_eos()` salva tabelas `eos.dat`.

### `src/core/physics.rs`

Centro do setor hadrônico:

- `HadronsMatter`:
  - resolve equilíbrio beta e neutralidade elétrica,
  - calcula densidades de partículas e campos médios,
  - monta energia e pressão totais.
- `NlemModel`:
  - `Maxwell`
  - `Modmax(csi)`
  - `Log(csi)`
- `MagneticTopology`:
  - `Isotropic`: $P_{mag} = \epsilon_{mag}/3$ (compatível com TOV 1D)
  - `Anisotropic`: $P_{mag} = \epsilon_{mag}$

### `src/core/tov_solver.rs`

- `generate_mr_curve(eps, p, with_crust)` gera sequência massa-raio.
- `integrate_star(pc, eps, p)` integra uma estrela para pressão central fixa.
- `unify_with_crust()` permite costura com crosta BPS.

### `src/core/plotting.rs`

- `Artist`: infraestrutura de plot multi-curva.
- suporte a escala log, autoscale e múltiplas séries por cenário.

### `src/core/io_utils.rs`

- `read_eos_file()` lê EoS externa e organiza em ordem crescente de pressão.

---

## Fluxo numérico principal

No caminho hadrônico (`HadronsMatter`), cada ponto em $\mu_n$ passa por:

1. solução numérica de variáveis de campo e potencial químico eletrônico,
2. definição do campo local $B$ pelo perfil (os níveis de Landau usam $B$; a NLEM entra só nas tensões do campo),
3. cálculo de densidades bariônicas e leptônicas,
4. cálculo de $\epsilon$ e $P$ totais,
5. armazenamento em linha de saída (`RESULTS_SIZE = 34`).

Após a EoS pronta:

1. limpeza/ordenação de dados para interpolação estável,
2. integração TOV para várias pressões centrais,
3. construção da curva M-R,
4. extração de observáveis (massa máxima, raio, etc.).

---

## Formato dos dados da EoS hadrônica

`RESULTS_SIZE = 34`

Cada linha de saída do solver contém:

| Índice | Significado (foco hadrônico) |
| ---: | --- |
| 0 | $n/n_0$ |
| 1 | Energia total $\epsilon$ [MeV/fm³] |
| 2 | Pressão total $P$ [MeV/fm³] |
| 3 | $n_{e^-}$ [fm⁻³] |
| 4 | $n_{\mu^-}$ [fm⁻³] |
| 5 | $n_n$ [fm⁻³] |
| 6 | $n_p$ [fm⁻³] |
| 7 | $n_{\Lambda^0}$ [fm⁻³] |
| 8 | $n_{\Sigma^-}$ [fm⁻³] |
| 9 | $n_{\Sigma^0}$ [fm⁻³] |
| 10 | $n_{\Sigma^+}$ [fm⁻³] |
| 11 | $n_{\Xi^-}$ [fm⁻³] |
| 12 | $n_{\Xi^0}$ [fm⁻³] |
| 13 | Potencial escalar $g_\sigma\sigma$ [MeV] |
| 14 | Potencial vetorial $g_\omega\omega_0$ [MeV] |
| 15 | Potencial isovetorial $g_\rho\rho_{03}$ [MeV] |
| 16 | $m^*/m_N$ |
| 17 | $\mu_n/M_N$ |
| 18 | $\mu_e/M_N$ |
| 19 | Energia magnética $\epsilon_{mag}$ [MeV/fm³] |
| 20 | Potencial químico fermiônico total por bárion dividido por $M_N$ |
| 21 | $n_\chi$ [fm⁻³] |
| 22 | $Y_\chi=n_\chi/n_B$ |
| 23 | $m_\chi$ [MeV] |
| 24 | $m_X$ [MeV] |
| 25 | $\epsilon$ |
| 26 | $g_D$ |
| 27 | $X_0$ [MeV] |
| 28 | $k_{F\chi}$ [MeV] |
| 29 | $\mu_\chi$ [MeV] |
| 30 | $\epsilon_\chi^{\rm kin}$ [MeV/fm³] |
| 31 | $P_\chi^{\rm kin}$ [MeV/fm³] |
| 32 | $\epsilon_X$ [MeV/fm³] |
| 33 | $P_X$ [MeV/fm³] |

As colunas 21–33 valem zero para motores sem setor escuro. Internamente,
massas, momentos e potenciais são normalizados por `M_NUCLEON`; densidades por
`M_NUCLEON³`; e energias/pressões por `M_NUCLEON⁴`. A saída aplica os fatores
de `HBAR_C` necessários para as unidades indicadas acima.

Quando `write_eos_with_mr` é usado, ele acrescenta `34:mr_mass_msun`,
`35:mr_radius_km` e `36:mr_baryonic_mass_msun`, totalizando 37 colunas. Arquivos
antigos com 21+3 colunas devem ser regenerados; apenas
`darkphotons_single.py` mantém leitura legada explícita desse formato.

### API de parâmetros escuros

- `with_m_x(m_x)` e `with_m_chi(m_chi)`: massas físicas divididas por
  `M_NUCLEON`.
- `with_m_x_mev(m_x_mev)` e `with_m_chi_mev(m_chi_mev)`: massas em MeV,
  convertidas uma única vez para a normalização interna.
- `with_g_d(g_d)`, `with_epsilon(epsilon)` e `with_y_chi(y_chi)`: parâmetros
  adimensionais independentes.
- Os antigos `with_n_chi` e `with_n_chi_natural` foram removidos; `n_chi` é
  agora estado calculado por `n_chi = y_chi * n_B`, não parâmetro constante.

`nsrs validate` aceita somente linhas uniformes com 34 colunas EOS ou 37 colunas
EOS+M-R. Para linhas escuras, ele também verifica a relação de fração, a
conversão entre `n_chi` e `kF_chi`, a definição de `mu_chi`, a solução neutra de
Proca e `eps_X = P_X`; isso evita validar silenciosamente arquivos truncados ou
com diagnósticos escuros corrompidos.

---

## Executável `nsrs` (`src/bin/nsrs`)

Todas as campanhas e ferramentas ficam num único executável com subcomandos:

```
cargo run --release --bin nsrs -- <comando> [argumentos]
cargo run --release --bin nsrs -- help
```

| Comando | Faz | Saída | Substitui |
|---|---|---|---|
| `scan b [--models GM1,GM3,FSU2] [--points 100] [--bmin Bc] [--bmax 3e18]` | EoS em função do campo constante (0 e malha log de `bmin` a `bmax`) | `output/b/` | `b` |
| `scan log <exp_min> <exp_max> <pontos_por_década> <B...>` | NLEM logarítmica, $\xi=10^{exp}$ | `output/nlem_log/` | `nlem_log` |
| `scan modmax <B...>` | ModMax, $\gamma = \{1..9\}\times10^{-10..-1}$, com `summary.csv` | `output/modmax/` | `nlem_modmax` |
| `scan topology <modelo> <B...> [--prefix TAG] [--plot-only]` | isotrópica x anisotrópica: EoS, M-R e populações | `output/magtop/`, `results/magtop/` | `magtop` |
| `study log <exp_min> <exp_max> <por_década> <B0...> [--topology aniso\|iso\|ambas] [--save-eos]` | impacto de Log($\xi$) na estrela: $M_{max}$, $R_{1.4}$, $\Lambda_{1.4}$, $B_c$, $x_c=B_c^2/2\xi^2$, sinal da pressão do campo, comparados com Maxwell e sem tensão; figuras com `plot_scripts/study_log.py` | `results/study_log/summary.csv` | — |
| `dark scan [--models GM1,GM3] [--b 1e17]` | grade $10^4$ em $(\epsilon, m_X, g_D, Y_\chi)$ | `output/darkphotons_scan/` | `darkphotons` |
| `dark benchmarks [--models GM1,GM3]` | cenários H0, S1-S3 (Kumar et al.), `summary.csv` transacional | `output/darkphotons_benchmarks/` | `darkphotons_base` |
| `dark single` | GM1, $B=10^{17}$ G, um ponto do setor escuro | `output/darkphotons/GM1/` | `single_darkphotons` |
| `quarks bag` | estrelas de quarks (MIT bag), varrendo $B_{bag}$ e $g_v$ | `results/bag_var_*.svg`, `results/gv_var_*.svg` | `bag_model` |
| `quarks hybrid` | hádrons x quarks x híbrida (Maxwell) | `results/comparison_*.svg` | `hybrid` |
| `report properties` | saturação, $M_{max}$, $R_{1.4}$, $\Lambda_{1.4}$, $I_{1.4}$, URCA, hyperons, $c_s^2$ | terminal | `properties` |
| `report validation [--out docs/VALIDATION_REPORT.md] [--constraints csv]` | relatório completo da validação em Markdown: soluções exatas, identidades termodinâmicas, parametrizações contra os artigos originais, relações universais, observações, perfis de campo e NLEM, com critérios e referências numeradas | `docs/VALIDATION_REPORT.md` | — |
| `report observations [constraints.csv]` | confronto com vínculos observacionais | `results/observations_report.csv` | `observations` |
| `validate <eos.dat\|pasta>... [opções]` | verificações de arquivos de EoS (`validate --help`) | terminal, `--csv` | `validate_eos` |
| `tov <eos.dat>` | curva M-R de uma EoS em arquivo | `results/mr_<nome>.svg` | `tov` |
| `plot <eos.dat> <saida.png> <col_x> <col_y>` | uma coluna contra outra (índices a partir de 0) | PNG | `bdd` |

Opção comum às varreduras: `--threads N` (padrão: todos os núcleos). Os
caminhos de saída e as malhas são os mesmos dos executáveis antigos. O
antigo `nlem_log_limits` foi removido: ele escolhia $\xi$ a partir de um campo
efetivo $B_{\rm ef}=B/(1+B^2/2\xi^2)$ nos níveis de Landau, premissa que não
vale mais (os níveis de Landau usam $B$ e a NLEM entra só na tensão do campo;
ver PHYSICS.md). `scan log` cobre a varredura em $\xi$.

Organização do código: `main.rs` (despacho e ajuda), `cli.rs` (argumentos),
`scan.rs`, `dark.rs`, `quarks.rs`, `report.rs`, `tools.rs` (`tov`, `plot`) e
`validate.rs`.

### `scan topology`

- compara topologias magnéticas (`Isotropic` vs `Anisotropic`) em vários campos $B$;
- `--plot-only` reaproveita EoS já calculadas; `--prefix TAG` organiza campanhas;
- salva `eos.dat` por topologia, extrai massa máxima e raio correspondente e
  plota EoS, M-R e população de partículas.

---

## Módulos principais (`src/core`)

- `constants.rs`: constantes físicas e tamanhos de vetor.
- `model.rs`: parametrizações hadrônicas (`GM1`, `GM3`, `FSU2`).
- `particles.rs`: densidades e estrutura de níveis de Landau.
- `eos.rs`: composição da EoS hadrônica.
- `physics.rs`: motor físico hadrônico `HadronsMatter` (NLEM, topologia, perfil de campo e setor escuro opcional).
- `magnetic.rs`: perfis do campo local (`Constant`, `Bdd`, `Dexheimer2017`) e tensões do campo (Maxwell, ModMax, Log).
- `darkphotons.rs`: gás de Dirac escuro; `DarkPhotonsMatter` é um apelido de `HadronsMatter`.
- `nuclear.rs`: propriedades de saturação (n0, E/A, K, J, L, M*/M).
- `tov_solver.rs` também integra maré (k2, Λ) e rotação lenta (I); `generate_star_sequence` devolve `StarProperties`.
- Saídas com `with_eos_output("x.dat")`: `x.dat` (EoS, 34 + 3 colunas M-R sem crosta), `x_stars.txt` (estrelas com crosta: M, R, M_B, P_c, C, z, k2, Λ, I, Ī) e `x_diag.txt` (M·B, c_s², Γ, frações, URCA direto).
- `nsrs report properties`: relatório de saturação e propriedades estelares dos modelos.
- `observations.rs` e `nsrs report observations`: confronto com vínculos observacionais e empíricos (`input/observations/constraints.csv`), relatório em `results/observations_report.csv`.
- `solver.rs`: varredura em $\mu_n$ e controle adaptativo.
- `tov_solver.rs`: integração de TOV e curva M-R.
- `plotting.rs`: infraestrutura de gráficos.
- `io_utils.rs`: leitura de EoS externa.

---

## Quarks e híbridas (em desenvolvimento)

As rotas de **quarks** e **híbridas** existem no código, mas neste momento são tratadas como trilhas em desenvolvimento nesta documentação:

- `quarks` (`src/core/quarks.rs`, `nsrs quarks bag`):
  - implementação baseada em MIT Bag com acoplamento vetorial,
  - foco em matéria de quarks pura e varreduras de parâmetros (`bag_constant`, `gv`).

- `hybrid` (`src/core/hybrid.rs`, `nsrs quarks hybrid`):
  - combina fase hadrônica e fase de quarks,
  - usa construção de Maxwell para decidir a fase estável ao longo de $\mu_n$.

---

## Diretórios de saída

- `output/`: dados intermediários (`eos.dat`) por campanha.
- `results/`: gráficos finais e tabelas-resumo (`.svg`, `.csv`, `.png`).

Regra prática:

- `output` para reuso em `--plot-only`,
- `results` para análise final.

---

## Como executar (foco hadrônico)

Exemplos:

- `cargo run --release --bin nsrs -- scan topology GM1 1e16 5e17 1e18`
- `cargo run --release --bin nsrs -- scan topology GM1 1e17 1e18 --plot-only --prefix novo`
- `cargo run --release --bin nsrs -- scan log 15 20 2 1e17 1e18 --models GM1`
- `cargo run --release --bin nsrs -- tov output/magtop/GM1/B_1.00e17/isotropic/eos.dat`
- `cargo run --release --bin nsrs -- validate output/magtop/GM1 --csv results/validacao.csv`
- `cargo run --release --bin nsrs -- report observations`

---

## Dicas de estabilidade numérica

1. Em falhas de convergência hadrônica:
   - aumente `n_points`,
   - reduza o intervalo de `mu_n` em `with_limits`.

2. Se M-R vier vazia:
   - valide monotonicidade de $P(\epsilon)$,
   - confirme pressão positiva no domínio relevante.

3. Se TOV falhar por interpolação:
   - use dados estritamente crescentes em pressão,
   - reaproveite a limpeza já implementada no `tov_solver`.

4. Em varreduras grandes:
   - execute primeiro sem `--plot-only`,
   - depois itere gráficos reaproveitando `output/`.
