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

    /// Resolve uma única EoS sequencialmente (necessário por causa do last_x)
    pub fn solve(&mut self) -> Vec<[f64; 13]> {
        let n = self.engine.n_points;
        let dmub = (self.engine.mun_sup - self.engine.mun_inf) / (n - 1) as f64;
        let mut results = Vec::with_capacity(n);
        
        // Chute inicial ideal para o vácuo (densidade inicial próxima de zero)
        let mut last_x = [0.0, 0.0, 0.0, 0.0]; 

        for i in 0..n {
            // Agora vamos do MENOR para o MAIOR (Crescente)
            let mun = self.engine.mun_inf + i as f64 * dmub;
            
            if let Some((x, point_result)) = self.engine.solve_point(mun, &last_x) {
                results.push(point_result);
                last_x = x; // A solução atual guia suavemente a próxima
            } else {
                println!(
                    "Finalizado: mun = {:.4} |  B = {:.2e} G",
                    mun, 
                    self.engine.bg
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