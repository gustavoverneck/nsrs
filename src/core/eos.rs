// src/solver/eos.rs

use crate::core::physics::PhysicsEngine;
use crate::core::constants::PI2;

pub fn compute(engine: &PhysicsEngine, mue: f64, vsigma: f64, vomega: f64, vrho: f64) -> (f64, f64) {
    // 1. Energia dos mésons (Potenciais de campo)
    // Inclui termos de massa e auto-interações (rb, rc para sigma e rxi para omega)
    let enerf = (vsigma / engine.model.gs).powi(2) / 2.0
        + (vomega / engine.model.gv).powi(2) / 2.0
        + (vrho / engine.model.gr).powi(2) / 2.0
        + engine.model.rb * vsigma.powi(3) / 3.0
        + engine.model.rc * vsigma.powi(4) / 4.0
        + engine.model.rxi * vomega.powi(4) / 4.0;

    let mut enerbar = 0.0;

    // --- LOOP DE BARIÕES (0:n, 1:p, 2:L0, 3:S-, 4:S0, 5:S+, 6:X-, 7:X0) ---
    for i in 0..8 {
        let ef = engine.ef_b[i];
        if ef <= 0.0 { continue; }

        // Se a partícula for neutra OU B=0, usa a fórmula contínua!
        if engine.charges_b[i] == 0.0 || engine.b == 0.0 {
            // --- Partículas Neutras (n, L0, S0, X0) ---
            // O AMM desdobra a partícula em 2 estados de spin (Up e Down)
            for &kf in &[engine.kf_b_up[i][0], engine.kf_b_down[i][0]] {
                if kf > 0.0 {
                    // Massa efetiva de spin derivada da cinemática: m_s^2 = Ef^2 - kf^2
                    let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                    
                    // Fórmula para um único estado de spin (g=1), fator 1/4pi^2
                    enerbar += (1.0 / (4.0 * PI2)) * (
                        ef.powi(3) * kf / 2.0
                        - (m_spin / 4.0) * (m_spin * kf * ef + m_spin.powi(3) * ((kf + ef) / m_spin.abs()).ln())
                    );
                }
            }
        } else {
            // --- Partículas Carregadas (p, S-, S+, X-) ---
            // Soma sobre os níveis de Landau (nu) para ambos os spins
            let qb = engine.charges_b[i].abs() * engine.qe * engine.b;
            let factor = qb / (4.0 * PI2);

            // Contribuição Spin Up
            for nu in 0..engine.n_b_up[i] {
                let kf = engine.kf_b_up[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerbar += factor * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }

            // Contribuição Spin Down
            for nu in 0..engine.n_b_down[i] {
                let kf = engine.kf_b_down[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerbar += factor * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }
        }
    }

        // --- LOOP DE LÉPTONS (0:e-, 1:mu-) ---
    let mut enerlep = 0.0;
    for i in 0..2 {
        let ef = engine.ef_l[i];
        if ef <= 0.0 { continue; }

        if engine.b == 0.0 {
            // Fórmula isotrópica para léptons se B=0 (com fator 2 para spin-up e spin-down)
            let kf = engine.f_l[i][0];
            if kf > 0.0 {
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerlep += 2.0 * (1.0 / (4.0 * PI2)) * (
                    ef.powi(3) * kf / 2.0
                    - (m_spin / 4.0) * (m_spin * kf * ef + m_spin.powi(3) * ((kf + ef) / m_spin.abs()).ln())
                );
            }
        } else {
            // Léptons sob Efeito de Landau (B > 0)
            let qb = engine.qe * engine.b;
            
            for nu in 0..engine.n_l[i] {
                let kf = engine.f_l[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                let g = if nu == 0 { 1.0 } else { 2.0 }; // Degenerescência de Landau para Dirac
                
                enerlep += (g * qb / (4.0 * PI2)) * (
                    ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln()
                );
            }
        }
    }

    // Energia Total (Mésons + Bariões + Léptons)
    let ener = enerf + enerbar + enerlep;

    // Pressão via relação termodinâmica: P = sum(mu_i * n_i) - epsilon
    let mut press_sum = 0.0;
    for i in 0..8 {
        press_sum += engine.mu_b[i] * engine.nb[i];
    }
    for i in 0..2 {
        press_sum += mue * engine.nl[i]; // mu_e = mu_mu = mue
    }
    
    let press = press_sum - ener;

    (ener, press)
}