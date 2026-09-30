use crate::core::constants::{M_NUCLEON, N0, RESULTS_SIZE};
// src/core/solver.rs
use crate::core::darkphotons::DarkPhotonsMatter;
use crate::core::hybrid::HybridMatter;
use crate::core::io_utils::write_eos_with_mr;
use crate::core::physics::HadronsMatter;
use crate::core::quarks::QuarksMatter;
use crate::core::tov_solver::generate_mr_curve;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

// 1. Enum que define o modelo físico a ser resolvido
pub enum EngineMode {
    Hadrons(HadronsMatter),
    Quarks(QuarksMatter),
    Hybrid(HybridMatter),
    DarkPhotons(DarkPhotonsMatter),
}

/// Motivo pelo qual a varredura em mu_n terminou. `mu_n_mev` é o potencial
/// químico do último ponto aceito e `nb_over_n0` a densidade correspondente.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EosTermination {
    /// A malha chegou a `mun_sup`.
    ReachedUpperLimit,
    /// Corte intencional em densidade (n_B > 11 n0).
    DensityCap,
    /// Corte intencional por causalidade (c_s^2 = dP/deps > 1).
    Acausal { mu_n_mev: f64, nb_over_n0: f64 },
    /// Corte intencional: a massa efetiva do nêutron (coluna 16) chegou a
    /// M*/M <= 0, fora da validade física do modelo de campo médio.
    NonPositiveEffectiveMass { mu_n_mev: f64, nb_over_n0: f64 },
    /// Queda resolvida de eps ou P entre pontos consecutivos.
    NonMonotonic { mu_n_mev: f64, nb_over_n0: f64 },
    /// O Newton não convergiu nem com o passo mínimo; a EoS está truncada.
    ConvergenceFailure { mu_n_mev: f64, nb_over_n0: f64 },
}

impl EosTermination {
    /// Terminações que truncam a EoS sem que isso tenha sido pedido.
    pub fn is_anomalous(&self) -> bool {
        matches!(
            self,
            EosTermination::NonMonotonic { .. } | EosTermination::ConvergenceFailure { .. }
        )
    }
}

/// Grandezas por ponto da EoS que não cabem no formato de 34 colunas.
/// Exportadas em `<saída>_diag.dat` (ver `io_utils::write_diagnostics`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointDiagnostics {
    /// M B = B dP/dB a mu fixo (MeV/fm^3); já descontado da pressão
    /// perpendicular exportada na coluna 2.
    pub magnetization_b: f64,
}

pub struct Solver {
    engine: EngineMode,
    termination: Option<EosTermination>,
    diagnostics: Vec<PointDiagnostics>,
}

impl Solver {
    pub fn new(engine: EngineMode) -> Self {
        Solver {
            engine,
            termination: None,
            diagnostics: Vec::new(),
        }
    }

    /// Motivo do término da última chamada a `solve()`.
    pub fn termination(&self) -> Option<EosTermination> {
        self.termination
    }

    /// Diagnósticos por linha da última chamada a `solve()` (mesma ordem das
    /// linhas da EoS).
    pub fn diagnostics(&self) -> &[PointDiagnostics] {
        &self.diagnostics
    }

    pub fn solve(&mut self) -> Vec<[f64; RESULTS_SIZE]> {
        // Obtém limites e metadados conforme o modo da engine
        let (mun_inf, mun_sup, n, bg_val) = match &self.engine {
            EngineMode::Hadrons(h) => (h.mun_inf, h.mun_sup, h.n_points, h.bg),
            EngineMode::Quarks(q) => (q.mun_inf, q.mun_sup, q.n_points, q.bg),
            EngineMode::Hybrid(hyb) => (
                hyb.hadrons.mun_inf,
                hyb.hadrons.mun_sup,
                hyb.hadrons.n_points,
                hyb.hadrons.bg,
            ),
            EngineMode::DarkPhotons(d) => (d.mun_inf, d.mun_sup, d.n_points, 0.0),
        };

        let initial_dmub = (mun_sup - mun_inf) / (n - 1) as f64;
        let mut dmub = initial_dmub;
        let min_dmub = 1e-6; // passo mínimo aceitável

        let mut results: Vec<[f64; RESULTS_SIZE]> = Vec::with_capacity(n);
        let mut diagnostics: Vec<PointDiagnostics> = Vec::with_capacity(n);
        let mut last_visible_x = [0.0; 5];
        let mut last_dark_x = [0.0; 5];
        let mut last_mun = mun_inf; // último mun que convergiu
        let mut mun = mun_inf;
        let mut termination = EosTermination::ReachedUpperLimit;
        // Apenas estas engines exportam M*/M do nêutron na coluna 16.
        let exports_effective_mass = matches!(
            self.engine,
            EngineMode::Hadrons(_) | EngineMode::DarkPhotons(_)
        );
        // Último ponto aceito, em (MeV, n_B/n0), para os diagnósticos.
        let last_point = |results: &[[f64; RESULTS_SIZE]]| {
            results
                .last()
                .map_or((mun_inf * M_NUCLEON, 0.0), |r| (r[17] * M_NUCLEON, r[0]))
        };

        while mun <= mun_sup + 1e-9 {
            // Tenta resolver o ponto com o mun atual
            let point_data =
                match &mut self.engine {
                    EngineMode::Hadrons(h_engine) => h_engine
                        .solve_point(mun, &last_visible_x)
                        .map(|(x, result)| {
                            last_visible_x = x;
                            result
                        }),
                    EngineMode::Quarks(q_engine) => q_engine.solve_point(mun),
                    EngineMode::Hybrid(hyb_engine) => hyb_engine
                        .solve_point(mun, &last_visible_x)
                        .map(|(x, result)| {
                            last_visible_x = x;
                            result
                        }),
                    EngineMode::DarkPhotons(d_engine) => {
                        d_engine.solve_point(mun, &last_dark_x).map(|(x, result)| {
                            last_dark_x = x;
                            result
                        })
                    }
                };

            if let Some(point_result) = point_data {
                last_mun = mun;
                let point_diagnostics = PointDiagnostics {
                    magnetization_b: match &self.engine {
                        EngineMode::Hadrons(h) | EngineMode::DarkPhotons(h) => h.magnetization_b,
                        _ => 0.0,
                    },
                };

                // Para a integração se a densidade bariônica ultrapassar 11 N0.
                if point_result[0] * N0 > 11.0 * N0 {
                    termination = EosTermination::DensityCap;
                    break;
                }

                if exports_effective_mass && point_result[16] <= 0.0 {
                    let (mu_n_mev, nb_over_n0) = last_point(&results);
                    termination = EosTermination::NonPositiveEffectiveMass {
                        mu_n_mev,
                        nb_over_n0,
                    };
                    break;
                }

                // Verificação de estabilidade (dP/dE). Perto do limiar
                // vácuo-matéria, diferenças no nível do arredondamento não
                // representam uma instabilidade física e não devem encerrar a
                // continuação da EOS.
                if !results.is_empty() {
                    let prev = results.last().unwrap();
                    let de = point_result[1] - prev[1];
                    let dp = point_result[2] - prev[2];
                    let resolved_matter = prev[0] > 1e-6 && point_result[0] > 1e-6;

                    let de_tol = 1e-10 * point_result[1].abs().max(prev[1].abs()).max(1.0);
                    let dp_tol = 1e-10 * point_result[2].abs().max(prev[2].abs()).max(1.0);

                    if resolved_matter && de > de_tol && dp > dp_tol {
                        let cs2 = dp / de;
                        if cs2 > 1.0 {
                            let (mu_n_mev, nb_over_n0) = last_point(&results);
                            termination = EosTermination::Acausal {
                                mu_n_mev,
                                nb_over_n0,
                            };
                            break;
                        }
                    } else if resolved_matter && (de < -de_tol || dp < -dp_tol) {
                        // Queda resolvida de energia ou pressão: esta sim é
                        // incompatível com a ramificação monótona usada pela TOV.
                        let (mu_n_mev, nb_over_n0) = last_point(&results);
                        termination = EosTermination::NonMonotonic {
                            mu_n_mev,
                            nb_over_n0,
                        };
                        break;
                    }
                }

                results.push(point_result);
                diagnostics.push(point_diagnostics);

                // Avança para o próximo mun
                mun += dmub;

                // Se o passo estava reduzido, tenta aumentá-lo gradualmente
                if dmub < initial_dmub {
                    dmub = (dmub * 1.5).min(initial_dmub);
                }
            } else {
                // Falha na convergência: reduz o passo e tenta novamente a partir do último sucesso
                dmub *= 0.5;
                if dmub < min_dmub {
                    let (mu_n_mev, nb_over_n0) = last_point(&results);
                    termination = EosTermination::ConvergenceFailure {
                        mu_n_mev,
                        nb_over_n0,
                    };
                    break;
                }
                mun = if results.is_empty() {
                    mun_inf
                } else {
                    last_mun + dmub
                };
            }
        }

        let output_path = match &self.engine {
            EngineMode::Hadrons(h) => h.eos_output.clone(),
            EngineMode::Quarks(q) => q.eos_output.clone(),
            EngineMode::Hybrid(h) => h.eos_output.clone(),
            EngineMode::DarkPhotons(d) => d.eos_output.clone(),
        };

        self.termination = Some(termination);
        if termination.is_anomalous() {
            eprintln!(
                "Aviso: EoS truncada antes de mun_sup = {:.2} MeV (B = {:.3e} G, saída = {}): {:?}",
                mun_sup * M_NUCLEON,
                bg_val,
                output_path.as_deref().unwrap_or("-"),
                termination
            );
        }

        if let Some(path) = &output_path {
            if let Err(error) =
                crate::core::io_utils::write_diagnostics(&results, &diagnostics, path)
            {
                eprintln!("failed to write diagnostics for '{}': {error}", path);
            }
        }
        self.diagnostics = diagnostics;

        if let Some(path) = output_path {
            let eps_arr: Vec<f64> = results.iter().map(|r| r[1]).collect();
            let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();
            let rho_arr: Vec<f64> = results.iter().map(|r| r[0]).collect();
            let (masses, radii, b_masses, pc_list) =
                generate_mr_curve(&eps_arr, &p_arr, &rho_arr, false);
            if let Err(error) =
                write_eos_with_mr(&results, &masses, &radii, &b_masses, &pc_list, &path)
            {
                eprintln!("failed to write EOS output '{}': {error}", path);
            }
            // Parallel production scans intentionally do not retain every EOS
            // in memory after it has been written.
            return Vec::new();
        }

        results
    }

    /// Resolve múltiplas EoS de forma paralela usando Rayon.
    pub fn solve_parallel(
        engines: Vec<EngineMode>,
        num_threads: usize,
    ) -> Vec<Vec<[f64; RESULTS_SIZE]>> {
        let pb = ProgressBar::new(engines.len() as u64);
        let style = ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar());
        pb.set_style(style);
        pb.set_message("NSRS");

        // 1. Criamos um construtor de pool de threads personalizado
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("Falha ao criar o ThreadPool do Rayon");

        // 2. Executamos o processamento paralelo dentro deste pool específico
        let results = pool.install(|| {
            engines
                .into_par_iter()
                .map(|engine| {
                    let mut solver = Solver::new(engine);
                    let result = solver.solve();
                    pb.inc(1);
                    result
                })
                .collect()
        });

        pb.finish_and_clear();
        results
    }

    /// Exporta os resultados da EoS para um arquivo formatado.
    pub fn write_eos(results: &[[f64; RESULTS_SIZE]], filename: &str) -> std::io::Result<()> {
        use std::io::Write; // Garante que o trait Write está no escopo
        let mut file = std::fs::File::create(filename)?;

        for data in results.iter() {
            // Transforma todas as colunas EOS em uma String separada por espaços.
            let line = data
                .iter()
                .map(|val| format!("{:12.5e}", val))
                .collect::<Vec<String>>()
                .join(" ");

            // 2. Escreve a linha inteira no arquivo
            writeln!(file, "{}", line)?;
        }
        Ok(())
    }
}
