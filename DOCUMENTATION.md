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
  - [Binaries principais (`src/bin`)](#binaries-principais-srcbin)
  - [`magtop`](#magtop)
  - [`nlem_limits`](#nlem_limits)
  - [`tov`](#tov)
  - [`bdd`](#bdd)
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
- `Solver::solve_batch()` acelera varreduras grandes em paralelo.
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
2. atualização do campo magnético efetivo (incluindo NLEM),
3. cálculo de densidades bariônicas e leptônicas,
4. cálculo de $\epsilon$ e $P$ totais,
5. armazenamento em linha de saída (`RESULTS_SIZE = 21`).

Após a EoS pronta:

1. limpeza/ordenação de dados para interpolação estável,
2. integração TOV para várias pressões centrais,
3. construção da curva M-R,
4. extração de observáveis (massa máxima, raio, etc.).

---

## Formato dos dados da EoS hadrônica

`RESULTS_SIZE = 21`

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
| 13 | Campo $\sigma$ |
| 14 | Campo $\omega$ |
| 15 | Campo $\rho$ |
| 16 | $m^*/m_N$ |
| 17 | $\mu_n$ |
| 18 | $\mu_e$ |
| 19 | Energia magnética $\epsilon_{mag}$ [MeV/fm³] |
| 20 | Campo magnético local $B(n)$ |

---

## Binaries principais (`src/bin`)

## `magtop`

Arquivo: `src/bin/magtop.rs`

Objetivo:

- comparar topologias magnéticas (`Isotropic` vs `Anisotropic`) em vários campos $B$,
- produzir EoS, curvas M-R e população de partículas por topologia.

CLI:

- `--plot-only`: reaproveita EoS já calculadas.
- `--prefix TAG`: organiza campanhas com prefixo de nome.

Fluxo:

1. cria dois motores hadrônicos por valor de $B$ (iso/aniso),
2. resolve em paralelo com `solve_batch`,
3. salva `eos.dat` por topologia,
4. reconstrói vetores para TOV,
5. extrai massa máxima e raio correspondente,
6. plota EoS, M-R e população de partículas.

---

## `nlem_limits`

Arquivo: `src/bin/nlem_limits.rs`

Objetivo:

- varrer o parâmetro NLEM logarítmico $\xi$ para diferentes campos magnéticos,
- estudar impacto em observáveis estelares.

CLI:

- `--plot-only`
- `<GM1|GM3> <exp_min> <exp_max> <B1> <B2> ...`

Malha de $\xi$:

- usa pontos intercalados $1.0\times 10^k$ e $5.0\times 10^k$.

Observáveis extraídos:

- massa máxima,
- raio no ponto de massa máxima,
- densidade central,
- massa efetiva central,
- energia magnética central.

---

## `tov`

Arquivo: `src/bin/tov.rs`

Objetivo:

- ler um `eos.dat` e gerar curva M-R diretamente.

Fluxo:

1. leitura com `read_eos_file()`,
2. integração TOV com `generate_mr_curve`,
3. saída em `results/mr_<nome>.svg`.

---

## `bdd`

Arquivo: `src/bin/bdd.rs`

Objetivo:

- utilitário de plot coluna vs coluna (PNG), útil para inspeção rápida de tabelas numéricas.

---

## Módulos principais (`src/core`)

- `constants.rs`: constantes físicas e tamanhos de vetor.
- `model.rs`: parametrizações hadrônicas (`GM1`, `GM3`).
- `particles.rs`: densidades e estrutura de níveis de Landau.
- `eos.rs`: composição da EoS hadrônica.
- `physics.rs`: motor físico hadrônico (inclui NLEM/topologia).
- `solver.rs`: varredura em $\mu_n$ e controle adaptativo.
- `tov_solver.rs`: integração de TOV e curva M-R.
- `plotting.rs`: infraestrutura de gráficos.
- `io_utils.rs`: leitura de EoS externa.

---

## Quarks e híbridas (em desenvolvimento)

As rotas de **quarks** e **híbridas** existem no código, mas neste momento são tratadas como trilhas em desenvolvimento nesta documentação:

- `quarks` (`src/core/quarks.rs`, `src/bin/bag_model.rs`):
  - implementação baseada em MIT Bag com acoplamento vetorial,
  - foco em matéria de quarks pura e varreduras de parâmetros (`bag_constant`, `gv`).

- `hybrid` (`src/core/hybrid.rs`, `src/bin/hybrid.rs`):
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

- `cargo run --bin magtop -- GM1 1e16 5e17 1e18`
- `cargo run --bin magtop -- --plot-only --prefix novo GM1 1e17 1e18`
- `cargo run --bin nlem_limits -- GM1 8 28 1e15 5e16 1e17`
- `cargo run --bin nlem_limits -- --plot-only GM1 8 28 1e15`
- `cargo run --bin tov output/magtop/GM1/B_1.00e17/isotropic/eos.dat`

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
