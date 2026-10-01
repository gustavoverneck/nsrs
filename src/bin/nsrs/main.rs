// src/bin/nsrs/main.rs
//
// Ponto de entrada único do NSRS. Uso:
//
//   cargo run --release --bin nsrs -- <comando> [argumentos]
//
// Execute sem argumentos (ou com `help`) para a lista de comandos.

mod cli;
mod dark;
mod quarks;
mod report;
mod scan;
mod study;
mod study_b;
mod tools;
mod validate;
mod validation;

const USAGE: &str = "\
NSRS - Neutron Stars Rust Solver

Uso: cargo run --release --bin nsrs -- <comando> [argumentos]

Varreduras (saída em output/):
  scan b        [--models GM1,GM3,FSU2] [--points 100] [--bmin Bc] [--bmax 3e18]
                Campo magnético de 0 e Bc até bmax (log).           -> output/b/
  scan log      <exp_min> <exp_max> <pontos_por_década> <B1> [B2 ...] [--models ...]
                NLEM logarítmica, xi = 10^exp_min .. 10^exp_max.     -> output/nlem_log/
  scan modmax   <B1> [B2 ...] [--models ...]
                ModMax, gamma = 1e-10 .. 9e-1.                      -> output/modmax/
  scan topology <GM1|GM3|FSU2> <B1> [B2 ...] [--prefix TAG] [--plot-only]
                Topologia isotrópica x anisotrópica, com gráficos.   -> output/magtop/, results/magtop/

Estudos (tabela para análise e publicação):
  study log     <exp_min> <exp_max> <pontos_por_década> <B0_1> [B0_2 ...]
                [--models GM1] [--topology aniso|iso|ambas] [--points 1500]
                [--no-hyperons] [--save-eos] [--out results/study_log/summary.csv]
                Impacto de Log(xi) na estrela, com Maxwell e sem tensão
                do campo como referências.                         -> results/study_log/

  study b       [--models GM1] [--bmin 1e14] [--bmax 1e20] [--per-decade 8]
                [--profiles constante,bdd] [--points 1500] [--mu-max 3.0] [--no-hyperons]
                [--no-stability] [--delta 1e-4] [--landau-max 20000] [--out results/study_b]
                Varredura completa em B (de 0 até quebrar): estrelas, cobertura da EoS
                e estabilidade mecânica e magnética local.            -> results/study_b/

Setor escuro:
  dark scan        [--models GM1,GM3] [--b 1e17]  Grade 10^4 (epsilon, m_X, g_D, Y_chi). -> output/darkphotons_scan/
  dark benchmarks  [--models GM1,GM3]             Cenários H0, S1-S3 (Kumar et al.).     -> output/darkphotons_benchmarks/
  dark single                                     GM1, B = 1e17 G, um ponto do setor.    -> output/darkphotons/GM1/

Quarks (em desenvolvimento):
  quarks bag       Estrelas de quarks (MIT bag): varre B_bag e g_v.  -> results/*.svg
  quarks hybrid    Hádrons (GM1) x quarks x híbrida (Maxwell).       -> results/comparison_*.svg

Relatórios de validação:
  report properties    Saturação, estrelas (M_max, R, Lambda, I), URCA, hyperons.
  report observations  [constraints.csv]  Confronto com vínculos observacionais. -> results/observations_report.csv
  report validation    [--out docs/VALIDATION_REPORT.md] [--constraints csv]
                       Relatório completo da validação, com literatura e citações. -> docs/VALIDATION_REPORT.md

Ferramentas:
  validate <eos.dat|pasta>... [opções]   Verifica arquivos de EoS (use 'validate --help').
  tov <eos.dat>                          Curva M-R de uma EoS externa.   -> results/mr_<nome>.svg
  plot <eos.dat> <saida.png> <col_x> <col_y>  Coluna contra coluna.

Opção comum às varreduras: --threads N (padrão: todos os núcleos).
";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (command, rest) = match args.split_first() {
        Some((c, rest)) => (c.as_str(), rest),
        None => ("help", &args[..0]),
    };
    let sub = rest.first().map(String::as_str).unwrap_or("");
    let tail = if rest.is_empty() { rest } else { &rest[1..] };

    let result = match (command, sub) {
        ("scan", "b") => scan::b(tail),
        ("scan", "log") => scan::log(tail),
        ("scan", "modmax") => scan::modmax(tail),
        ("scan", "topology") => scan::topology(tail),
        ("study", "log") => study::log(tail),
        ("study", "b") => study_b::run(tail),
        ("dark", "scan") => dark::scan(tail),
        ("dark", "benchmarks") => dark::benchmarks(tail),
        ("dark", "single") => dark::single(tail),
        ("quarks", "bag") => quarks::bag(tail),
        ("quarks", "hybrid") => quarks::hybrid(tail),
        ("report", "properties") => report::properties(tail),
        ("report", "observations") => report::observations(tail),
        ("report", "validation") => validation::run(tail),
        ("validate", _) => std::process::exit(validate::run(rest)),
        ("tov", _) => tools::tov(rest),
        ("plot", _) => tools::plot(rest),
        ("help" | "-h" | "--help", _) => {
            print!("{USAGE}");
            return;
        }
        _ => Err(format!("comando desconhecido: {}", args.join(" "))),
    };

    if let Err(message) = result {
        eprintln!("erro: {message}\n");
        eprint!("{USAGE}");
        std::process::exit(2);
    }
}
