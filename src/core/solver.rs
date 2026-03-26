// src/core/solver.rs
use crate::core::physics::HadronsMatter;
use crate::core::quarks::QuarksMatter;
use crate::core::hybrid::HybridMatter;
use rayon::prelude::*;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use console::style;
use std::io::Write;

// 1. Enum que define o modelo físico a ser resolvido
pub enum EngineMode {
    Hadrons(HadronsMatter),
    Quarks(QuarksMatter),
    Hybrid(HybridMatter),
}

pub struct Solver {
    engine: EngineMode,
}

impl Solver {
    pub fn new(engine: EngineMode) -> Self {
        Solver { engine }
    }

    pub fn solve(&mut self) -> Vec<[f64; 20]> {
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
        };

        let initial_dmub = (mun_sup - mun_inf) / (n - 1) as f64;
        let mut dmub = initial_dmub;
        let min_dmub = 1e-6;   // passo mínimo aceitável

        let mut results: Vec<[f64; 20]> = Vec::with_capacity(n);
        let mut last_x = [0.0, 0.0, 0.0, 0.0];
        let mut last_mun = mun_inf;   // último mun que convergiu
        let mut mun = mun_inf;

        while mun <= mun_sup + 1e-9 {
            // Tenta resolver o ponto com o mun atual
            // Agora todos os braços retornam Option<([f64;4], [f64;20])>
            let point_data = match &mut self.engine {
                EngineMode::Hadrons(h_engine) => h_engine.solve_point(mun, &last_x),
                // Para Quarks, transformamos o resultado em tupla com estado dummy
                EngineMode::Quarks(q_engine) => q_engine.solve_point(mun).map(|result| ([0.0; 4], result)),
                EngineMode::Hybrid(hyb_engine) => hyb_engine.solve_point(mun, &last_x),
            };

            if let Some((x_converged, point_result)) = point_data {
                // Sucesso: atualiza chute, faz verificações e armazena
                last_x = x_converged;
                last_mun = mun;

                // Verificação de estabilidade (dp/de)
                if !results.is_empty() {
                    let prev = results.last().unwrap();
                    let de = point_result[1] - prev[1];
                    let dp = point_result[2] - prev[2];

                    if de > 0.0 {
                        let cs2 = dp / de;
                        if cs2 < 1e-6 {
                            println!(
                                "EoS instável (dp/de = {:.2e}) em mun = {:.4}. Encerrando.",
                                cs2, mun
                            );
                            break;
                        }
                        if cs2 > 1.1 {
                            println!(
                                "Aviso: EoS não-causal (dp/de = {:.2e}) em mun = {:.4}",
                                cs2, mun
                            );
                        }
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
                    println!(
                        "Limite atingido: não foi possível avançar após mun = {:.4} (passo mínimo).",
                        last_mun
                    );
                    break;
                }
                println!(
                    "Dificuldade de convergência em mun = {:.4}. Reduzindo passo para {:.1e}",
                    mun, dmub
                );
                // Retrocede para o último ponto bem-sucedido e tenta com passo menor
                mun = if results.is_empty() {
                    mun_inf
                } else {
                    last_mun + dmub
                };
            }
        }

        results
    }

    /// Resolve múltiplas EoS de forma paralela usando Rayon.
    pub fn solve_batch(engines: Vec<EngineMode>) -> Vec<Vec<[f64; 20]>> {
        engines
            .into_par_iter()
            .map(|engine| {
                let mut solver = Solver::new(engine);
                solver.solve()
            })
            .collect()
    }

    /// Exporta os resultados da EoS para um arquivo formatado.
    pub fn write_eos(results: &[[f64; 20]], filename: &str) -> std::io::Result<()> {
        use std::io::Write; // Garante que o trait Write está no escopo
        let mut file = std::fs::File::create(filename)?;

        for data in results.iter() {
            // 1. Transforma os 20 valores em uma única String separada por espaços
            let line = data.iter()
                .map(|val| format!("{:12.5e}", val))
                .collect::<Vec<String>>()
                .join(" ");

            // 2. Escreve a linha inteira no arquivo
            writeln!(file, "{}", line)?;
        }
        Ok(())
    }
}