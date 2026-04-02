// src/core/hybrid.rs

use crate::core::constants::RESULTS_SIZE;
use crate::core::physics::HadronsMatter;
use crate::core::quarks::QuarksMatter;

pub struct HybridMatter {
    pub hadrons: HadronsMatter,
    pub quarks: QuarksMatter,
}

impl HybridMatter {
    pub fn new(hadrons: HadronsMatter, quarks: QuarksMatter) -> Self {
        Self { hadrons, quarks }
    }

    /// Resolve o ponto híbrido comparando as pressões (Construção de Maxwell).
    /// Retorna Option<([Chute_Newton_Raphson], [Dados_Termodinamicos])>
    pub fn solve_point(&mut self, mun: f64, last_x_hadrons: &[f64; 4]) -> Option<([f64; 4], [f64; RESULTS_SIZE])> {
        // 1. Calcula a termodinâmica para os Hádrons (retorna tupla com estado numérico)
        let had_res = self.hadrons.solve_point(mun, last_x_hadrons);
        
        // 2. Calcula a termodinâmica para os Quarks (retorna apenas dados físicos)
        let qrk_res = self.quarks.solve_point(mun);
        
        // 3. A CONSTRUÇÃO DE MAXWELL: Quem ganha a batalha da Pressão?
        match (had_res, qrk_res) {
            (Some((x_converged, h)), Some(q)) => {
                let press_had = h[2]; // Índice 2 é a pressão total (termo + magnética)
                let press_qrk = q[2];
                
                if press_had > press_qrk {
                    // Fase Hadrônica é mais estável
                    Some((x_converged, h))
                } else {
                    // Fase de Quarks venceu. 
                    // Mantemos o last_x_hadrons para não "perder o fio da meada" numérico
                    // caso a densidade mude novamente.
                    Some((*last_x_hadrons, q))
                }
            },
            
            // Se apenas uma fase convergir, usamos a que estiver disponível
            (Some((x, h)), None) => Some((x, h)),
            (None, Some(q)) => Some((*last_x_hadrons, q)),
            
            _ => None
        }
    }
}