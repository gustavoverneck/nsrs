// solver/solver.rs

use crate::solver::physics::PhysicsEngine;
use std::fs::File;
use std::io::Write;

pub struct Solver {
    engine: PhysicsEngine,
}

impl Solver {
    pub fn new(engine: PhysicsEngine) -> Self {
        Solver { engine }
    }

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
                println!("Solver parou de convergir no núcleo em mun = {:.4}.", mun);
                break;
            }
        }
        results
    }

    pub fn write_eos(&self, results: &[[f64; 13]], filename: &str) -> std::io::Result<()> {
        let mut file = std::fs::File::create(filename)?;
        // Removemos o .rev() para manter a ordem natural gravada no ficheiro
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