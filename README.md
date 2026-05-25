# NSRS

NSRS é um solver em Rust para estudar estrelas de nêutrons, com foco atual em **estrelas hadrônicas magnetizadas**.

O projeto resolve a **equação de estado (EoS)** para um modelo hadrônico, integra as equações de **Tolman–Oppenheimer–Volkoff (TOV)** e gera curvas **massa-raio (M-R)**, além de gráficos e tabelas para análise.

## 📖 Documentação

Para informações detalhadas, consulte os arquivos na pasta [docs/](docs/):

- **[Documentação Técnica](docs/DOCUMENTATION.md)**: Detalhes sobre a arquitetura do código, módulos, fluxo numérico e como executar o projeto.
- **[Fundamentos Físicos](docs/PHYSICS.md)**: Detalhes sobre os modelos físicos e equações utilizadas.

## O que o projeto faz (visão superficial)

1. Configura um modelo físico hadrônico (GM1/GM3, campo magnético, topologia magnética e NLEM).
2. Resolve a EoS ponto a ponto em potencial químico bariônico.
3. Usa a EoS para integrar TOV e construir a sequência estelar.
4. Exporta resultados em `output/` (dados intermediários) e `results/` (gráficos/CSV finais).

## Estrutura em alto nível

- [docs/](docs/): documentação detalhada técnica e física.
- `src/core/`: núcleo físico e numérico (EoS, solver, TOV, plot, I/O)
- `src/bin/`: executáveis para análises e varreduras
