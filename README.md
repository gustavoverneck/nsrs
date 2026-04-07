# NSRS

NSRS é um solver em Rust para estudar estrelas de nêutrons, com foco atual em **estrelas hadrônicas magnetizadas**.

O projeto resolve a **equação de estado (EoS)** para um modelo hadrônico, integra as equações de **Tolman–Oppenheimer–Volkoff (TOV)** e gera curvas **massa-raio (M-R)**, além de gráficos e tabelas para análise.

## O que o projeto faz (visão superficial)

1. Configura um modelo físico hadrônico (GM1/GM3, campo magnético, topologia magnética e NLEM).
2. Resolve a EoS ponto a ponto em potencial químico bariônico.
3. Usa a EoS para integrar TOV e construir a sequência estelar.
4. Exporta resultados em `output/` (dados intermediários) e `results/` (gráficos/CSV finais).

## Estrutura em alto nível

- `src/core/`: núcleo físico e numérico (EoS, solver, TOV, plot, I/O)
- `src/bin/`: executáveis para análises e varreduras

Executáveis mais usados no fluxo hadrônico:

- `magtop`: compara topologia isotrópica vs anisotrópica
- `nlem_limits`: varre parâmetros NLEM e extrai limites observáveis
- `tov`: gera curva M-R a partir de um `eos.dat`

## Estado atual

- **Hadrônico**: principal trilha de uso e documentação.
- **Quarks e híbridas**: presentes no código como desenvolvimento em andamento.

## Documentação

- Documentação técnica do pipeline (foco hadrônico): [DOCUMENTATION.md](DOCUMENTATION.md)
- Modelos e equações físicas: [PHYSICS.md](PHYSICS.md)

## Execução rápida

Exemplos:

- `cargo run --bin magtop -- GM1 1e16 5e17 1e18`
- `cargo run --bin nlem_limits -- GM1 8 28 1e15 5e16 1e17`
- `cargo run --bin tov output/magtop/GM1/B_1.00e17/isotropic/eos.dat`
