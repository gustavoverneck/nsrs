// src/solver/solver.rs

use crate::core::physics::PhysicsEngine;
use rayon::prelude::*;
use std::io::Write;

pub struct Solver {
    engine: PhysicsEngine,
}

impl Solver {
    pub fn new(engine: PhysicsEngine) -> Self {
        Solver { engine }
    }

    /// Resolve a EoS sequencialmente
    pub fn solve(&mut self) -> Vec<[f64; 13]> {
        let n = self.engine.n_points;
        let dmub = (self.engine.mun_sup - self.engine.mun_inf) / (n - 1) as f64;
        let mut results: Vec<[f64; 13]> = Vec::with_capacity(n);
        
        let mut last_x = [0.0, 0.0, 0.0, 0.0];

        for i in 0..n {
            let mun = self.engine.mun_inf + i as f64 * dmub;
            
            if let Some((x, point_result)) = self.engine.solve_point(mun, &last_x) {
                
                // --- VERIFICAÇÃO DE ESTABILIDADE (dp/de) ---
                if !results.is_empty() {
                    let prev = results.last().unwrap();
                    let de = point_result[1] - prev[1]; // Variação da Densidade de Energia
                    let dp = point_result[2] - prev[2]; // Variação da Pressão

                    if de > 0.0 {
                        let cs2 = dp / de; // Velocidade do som ao quadrado c_s^2
                        
                        // Se dp/de for quase zero ou negativo, a estrela é instável.
                        // Abortamos para evitar que o TOV gere resultados errôneos.
                        if cs2 < 1e-6 {
                            println!(
                                "Abortando: EoS instável (dp/de = {:.2e}) | mun = {:.4} | B = {:.2e} G",
                                cs2, mun, self.engine.bg
                            );
                            break;
                        }

                        // Verificação opcional de causalidade (cs2 <= 1)
                        if cs2 > 1.1 {
                            println!("Aviso: EoS não-causal (dp/de = {:.2e}) em mun = {:.4}", cs2, mun);
                        }
                    }
                }

                results.push(point_result); 
                last_x = x;
            } else {
                println!(
                    "Abortando: mun = {:.4} | Modelo (gs) = {:.2} | B = {:.2e} G",
                    mun, self.engine.model.gs, self.engine.bg
                );
                break;
            }
        }
        results
    }

    /// Resolve múltiplas EoS de forma paralela usando Rayon.
    pub fn solve_batch(engines: Vec<PhysicsEngine>) -> Vec<Vec<[f64; 13]>> {
        engines
            .into_par_iter()
            .map(|engine| {
                let mut solver = Solver::new(engine);
                solver.solve()
            })
            .collect()
    }

    /// Escreve os resultados em um arquivo.
    pub fn write_eos(results: &[[f64; 13]], filename: &str) -> std::io::Result<()> {
        let mut file = std::fs::File::create(filename)?;
        for data in results.iter() {
            writeln!(
                file,
                "{:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}",
                data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12]
            )?;
        }
        Ok(())
    }
}