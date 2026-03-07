// src/solver/particles.rs

use crate::core::physics::PhysicsEngine;
use crate::core::constants::PI2;

pub fn calculate_all_densities(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) {
    for i in 0..8 {
        let (rs, rb) = if engine.charges_b[i] == 0.0 {
            density_baryon_neutral(engine, i, vomega, vrho)
        } else {
            density_baryon_charged(engine, i, vomega, vrho)
        };
        engine.rhos_b[i] = rs;
        engine.nb[i] = rb;
    }

    for i in 0..2 {
        let (rl, nl) = density_lepton(engine, i);
        engine.rhos_l[i] = rl;
        engine.nl[i] = nl;
    }

    let rhos_n_p = engine.rhos_b[0] + engine.rhos_b[1];
    let rhos_h = engine.rhos_b[2..8].iter().sum::<f64>();
    
    // CORREÇÃO: O / PI2 foi removido. As densidades já vieram divididas das funções abaixo!
    engine.rhosb = rhos_n_p + engine.xs * rhos_h; 
    engine.nbt = engine.nb.iter().sum();
}

pub fn density_baryon_neutral(
    engine: &mut PhysicsEngine, 
    idx: usize, 
    vomega: f64, 
    vrho: f64
) -> (f64, f64) {
    let ef = engine.mu_b[idx] 
           - (engine.xv_v[idx] * vomega) 
           - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx]);
    
    engine.ef_b[idx] = ef;
    
    // BLINDAGEM 1: Se a energia de Fermi efetiva é negativa, não existe matéria!
    if ef <= 0.0 { return (0.0, 0.0); }

    let m_star = engine.m_eff[idx];    
    let amm = engine.amm_b[idx];      
    let b = engine.b;                 

    let m_up = m_star - amm * b;
    let m_down = m_star + amm * b;

    let mut rhos_total = 0.0;
    let mut dens_total = 0.0;

    for &m_spin in &[m_up, m_down] {
        let kf2 = ef.powi(2) - m_spin.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            // BLINDAGEM 2: Previne divisão por zero ou logaritmo infinito se a massa colapsar
            let m_safe = m_spin.abs().max(1e-15); 
            
            rhos_total += (m_spin / (4.0 * PI2)) * (
                ef * kf - m_spin.powi(2) * ((kf + ef) / m_safe).ln()
            );
            
            dens_total += kf.powi(3) / (6.0 * PI2);
            
            if m_spin == m_up { engine.kf_b_up[idx][0] = kf; } 
            else { engine.kf_b_down[idx][0] = kf; }
        }
    }
    (rhos_total, dens_total)
}

fn density_baryon_charged(engine: &mut PhysicsEngine, idx: usize, vomega: f64, vrho: f64) -> (f64, f64) {
    let q = engine.charges_b[idx].abs() * engine.qe;
    let b = engine.b;
    let m = engine.m_eff[idx]; // m*
    let amm = engine.amm_b[idx];
    
    let ef = engine.mu_b[idx] 
             - (engine.xv_v[idx] * vomega) 
             - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx]); // Correção do sinal do isospin (Item 6B)
    
    engine.ef_b[idx] = ef;
    
    if ef <= 0.0 { return (0.0, 0.0); }

    // TRATAMENTO PARA O CASO ISOTRÓPICO (B=0)
    if b == 0.0 {
        let kf2 = ef.powi(2) - m.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let dens = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos = (m / (2.0 * PI2)) * (ef * kf - m.powi(2) * ((kf + ef) / m_safe).ln());
            
            // Salvamos no índice 0 como se fosse um único "nível macro" para a EoS usar
            engine.kf_b_up[idx][0] = kf;
            engine.kf_b_down[idx][0] = kf;
            engine.n_b_up[idx] = 1;
            engine.n_b_down[idx] = 1;
            
            return (rhos, dens);
        }
        return (0.0, 0.0);
    }

    let nu_max_approx_up = ((ef + amm * b).powi(2) - m.powi(2)) / (2.0 * q * b);
    let nu_max_approx_down = ((ef - amm * b).powi(2) - m.powi(2)) / (2.0 * q * b);
    
    let nu_max = if nu_max_approx_up > 0.0 || nu_max_approx_down > 0.0 { 
        let max_nu = nu_max_approx_up.max(nu_max_approx_down);
        (max_nu.floor() as usize + 1).min(engine.max_landau_limit) 
    } else { 
        0 
    };

    let q_sign = engine.charges_b[idx].signum();
    let (nu_start_up, nu_start_down) = if q_sign > 0.0 {
        (0, 1) // Cargas positivas: Spin UP tem nu=0
    } else {
        (1, 0) // Cargas negativas: Spin DOWN tem nu=0
    };

    let (mut rhos, mut dens) = (0.0, 0.0);
    let mut n_up = 0;
    let mut n_down = 0;

    // --- Spin UP ---
    for nu in nu_start_up..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau - amm * b;
        
        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 { break; } 
        
        let kf = kf2.sqrt();
        // CORREÇÃO: Usa n_up como índice para gravar os dados sequencialmente
        engine.kf_b_up[idx][n_up] = kf;
        n_up += 1;
        
        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15); 
        
        rhos += (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    // --- Spin DOWN ---
    for nu in nu_start_down..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau + amm * b;
        
        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 { break; }
        
        let kf = kf2.sqrt();
        // CORREÇÃO: Usa n_down como índice. Isso evita o underflow de (nu - 1)
        engine.kf_b_down[idx][n_down] = kf;
        n_down += 1;
        
        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15);
        
        rhos += (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    engine.n_b_up[idx] = n_up;
    engine.n_b_down[idx] = n_down;

    (rhos, dens)
}

pub fn density_lepton(engine: &mut PhysicsEngine, idx: usize) -> (f64, f64) {
    let mue = engine.mue;      
    
    // BLINDAGEM: Lêptons com energia negativa não existem
    if mue <= 0.0 { 
        engine.n_l[idx] = 0;
        engine.ef_l[idx] = 0.0;
        return (0.0, 0.0); 
    }

    let b = engine.b;
    let q = engine.qe; 
    let m = engine.ml[idx];

    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut n_occupied = 0; 

    if b == 0.0 {
        let kf2 = mue.powi(2) - m.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let dens = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos = (m / (2.0 * PI2)) * (mue * kf - m.powi(2) * ((kf + mue) / m_safe).ln());
            
            engine.f_l[idx][0] = kf;
            engine.n_l[idx] = 1;
            engine.ef_l[idx] = mue;
            
            return (rhos, dens);
        }
        return (0.0, 0.0);
    }


    let nu_max_approx = (mue.powi(2) - m.powi(2)) / (2.0 * q * b);
    let nu_max = if nu_max_approx > 0.0 { 
        (nu_max_approx.floor() as usize + 1).min(engine.max_landau_limit) 
    } else { 0 };

    for nu in 0..nu_max {
        let m_landau_2 = m.powi(2) + 2.0 * q * b * nu as f64;
        let kf2 = mue.powi(2) - m_landau_2;
        
        if kf2 <= 0.0 { break; } 
        
        let kf = kf2.sqrt();
        let g = if nu == 0 { 1.0 } else { 2.0 }; 

        engine.f_l[idx][nu] = kf;

        let factor = (g * q * b) / (2.0 * PI2);
        let m_safe = m_landau_2.sqrt().max(1e-15);
        
        rhos += factor * m * ((kf + mue) / m_safe).ln();
        dens += factor * kf;
        
        n_occupied += 1;
    }

    engine.n_l[idx] = n_occupied; 
    engine.ef_l[idx] = mue;   

    (rhos, dens)
}