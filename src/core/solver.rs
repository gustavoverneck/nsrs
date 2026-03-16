// src/core/solver.rs
use crate::core::physics::HadronsMatter;
use crate::core::quarks::QuarksMatter;
use crate::core::hybrid::HybridMatter;
use rayon::prelude::*;
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
        // Puxa os limites e metadados dependendo do modelo ativo
        let (mun_inf, mun_sup, n, bg_val) = match &self.engine {
            EngineMode::Hadrons(h) => (h.mun_inf, h.mun_sup, h.n_points, h.bg),
            EngineMode::Quarks(q) => (q.mun_inf, q.mun_sup, q.n_points, q.bg),
            // No modo híbrido, usamos os limites definidos no motor de hádrons
            EngineMode::Hybrid(hyb) => (
                hyb.hadrons.mun_inf, 
                hyb.hadrons.mun_sup, 
                hyb.hadrons.n_points, 
                hyb.hadrons.bg
            ),
        };

        let dmub = (mun_sup - mun_inf) / (n - 1) as f64;
        let mut results: Vec<[f64; 20]> = Vec::with_capacity(n);
        let mut last_x = [0.0, 0.0, 0.0, 0.0];

        for i in 0..n {
            let mun = mun_inf + i as f64 * dmub;
            
            // 2. Escolha do método de resolução baseado no modo da Engine
            let point_data = match &mut self.engine {
                EngineMode::Hadrons(h_engine) => {
                    if let Some((x_converged, result)) = h_engine.solve_point(mun, &last_x) {
                        last_x = x_converged; // Atualiza chute para o próximo mun
                        Some(result)
                    } else {
                        None
                    }
                },
                EngineMode::Quarks(q_engine) => {
                    // Quarks são resolvidos analiticamente (Beta-Equilibrium)
                    q_engine.solve_point(mun)
                },
                EngineMode::Hybrid(hyb_engine) => {
                    // Híbrido compara Maxwell e também gerencia o estado x dos hádrons
                    if let Some((x_converged, result)) = hyb_engine.solve_point(mun, &last_x) {
                        last_x = x_converged;
                        Some(result)
                    } else {
                        None
                    }
                }
            };

            // 3. Processamento de estabilidade e armazenamento
            if let Some(point_result) = point_data {
                
                // --- VERIFICAÇÃO DE ESTABILIDADE (dp/de) ---
                if !results.is_empty() {
                    let prev = results.last().unwrap();
                    let de = point_result[1] - prev[1]; // Delta Energia
                    let dp = point_result[2] - prev[2]; // Delta Pressão

                    if de > 0.0 {
                        let cs2 = dp / de; // Velocidade do som ao quadrado
                        
                        if cs2 < 1e-6 {
                            println!(
                                "Abortando: EoS instável (dp/de = {:.2e}) | mun = {:.4} | B = {:.2e} G",
                                cs2, mun, bg_val
                            );
                            break;
                        }

                        if cs2 > 1.1 {
                            println!("Aviso: EoS não-causal (dp/de = {:.2e}) em mun = {:.4}", cs2, mun);
                        }
                    }
                }

                results.push(point_result); 
            } else {
                println!(
                    "Abortando: Falha na convergência da EoS | mun = {:.4} | B = {:.2e} G",
                    mun, bg_val
                );
                break;
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