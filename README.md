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

## Dependências e Pré-requisitos

Para compilar o núcleo em Rust e executar os scripts de análise avançada, você precisará das seguintes ferramentas instaladas no seu sistema:

### 1. Dependências do Sistema (C/Rust)
O NSRS utiliza a biblioteca `rgsl` (bindings de Rust para a GNU Scientific Library) para integrações numéricas (Runge-Kutta) e interpolações. É necessário ter a biblioteca original do GSL instalada no sistema operacional:

* **Ubuntu / Debian / Linux Mint:**
    ```bash
    sudo apt-get install libgsl-dev
    ```
* **Arch Linux / Manjaro:**
    ```bash
    sudo pacman -S gsl
    ```
* **macOS (via Homebrew):**
    ```bash
    brew install gsl
    ```
* **Fedora / CentOS:**
    ```bash
    sudo dnf install gsl-devel
    ```