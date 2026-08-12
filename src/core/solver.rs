use crate::core::constants::{N0, RESULTS_SIZE};
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

pub struct Solver {
    engine: EngineMode,
}

impl Solver {
    pub fn new(engine: EngineMode) -> Self {
        Solver { engine }
    }

    pub fn solve(&mut self) -> Vec<[f64; RESULTS_SIZE]> {
        // Obtém limites e metadados conforme o modo da engine
        let (mun_inf, mun_sup, n, _bg_val) = match &self.engine {
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
        let mut last_visible_x = [0.0; 4];
        let mut last_dark_x = [0.0; 5];
        let mut last_mun = mun_inf; // último mun que convergiu
        let mut mun = mun_inf;

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

                // Para a integração se a densidade bariônica ultrapassar 15 N0.
                if point_result[0] * N0 > 11.0 * N0 {
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
                        if cs2 > 1.1 {
                            // println!(
                            //     "Aviso: EoS não-causal (dp/de = {:.2e}) em mun = {:.4}",
                            //     cs2, mun
                            // );
                            break;
                        }
                    } else if resolved_matter && (de < -de_tol || dp < -dp_tol) {
                        // Queda resolvida de energia ou pressão: esta sim é
                        // incompatível com a ramificação monótona usada pela TOV.
                        break;
                    }
                }

                results.push(point_result);

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
