// src/solver/eos.rs

use crate::core::constants::PI2;
use crate::core::particles::neutral_amm_spin;
use crate::core::physics::HadronsMatter;

pub fn compute(
    engine: &HadronsMatter,
    mue: f64,
    vsigma: f64,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    // 1. Energia dos mésons (Potenciais de campo)
    // Inclui termos de massa e auto-interações (rb, rc para sigma e rxi para omega).
    // Para os campos vetoriais, eps = g_v*omega*n_B + g_rho*rho*n_3 - L_mesons, logo
    // zeta/24 (g_v omega)^4 -> zeta/8 = 3 rxi/4 e Lambda_v(...)(...) -> 3 Lambda_v.
    let enerf = (vsigma / engine.model.gs).powi(2) / 2.0
        + (vomega / engine.model.gv).powi(2) / 2.0
        + (vrho / engine.model.gr).powi(2) / 2.0
        + engine.model.rb * vsigma.powi(3) / 3.0
        + engine.model.rc * vsigma.powi(4) / 4.0
        + 3.0 * engine.model.rxi * vomega.powi(4) / 4.0
        + 3.0 * engine.model.lambda_v * vomega.powi(2) * vrho.powi(2);

    let mut enerbar = 0.0;

    // --- LOOP DE BARIÕES (0:n, 1:p, 2:L0, 3:S-, 4:S0, 5:S+, 6:X-, 7:X0) ---
    for i in 0..8 {
        let ef = engine.ef_b[i];
        if ef <= 0.0 {
            continue;
        }

        // Se a partícula for neutra OU B=0, usa a fórmula contínua!
        if engine.charges_b[i] == 0.0 || engine.b == 0.0 {
            // --- Partículas Neutras (n, L0, S0, X0) ---
            // Com AMM os dois estados de spin têm a = +/- kappa mu_N B; com
            // B = 0 ou kappa = 0 ambos coincidem com o gás isotrópico.
            let a = engine.amm_b[i] * engine.b;
            for a_s in [a, -a] {
                if let Some(state) = neutral_amm_spin(engine.m_eff[i], ef, a_s) {
                    enerbar += state.energy;
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
        if ef <= 0.0 {
            continue;
        }

        if engine.b == 0.0 {
            // Fórmula isotrópica para léptons se B=0 (com fator 2 para spin-up e spin-down)
            let kf = engine.f_l[i].first().copied().unwrap_or(0.0);
            if kf > 0.0 {
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerlep += 2.0
                    * (1.0 / (4.0 * PI2))
                    * (ef.powi(3) * kf / 2.0
                        - (m_spin / 4.0)
                            * (m_spin * kf * ef
                                + m_spin.powi(3) * ((kf + ef) / m_spin.abs()).ln()));
            }
        } else {
            // Léptons sob Efeito de Landau (B > 0)
            let qb = engine.qe * engine.b;

            for nu in 0..engine.n_l[i] {
                let kf = engine.f_l[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                let g = if nu == 0 { 1.0 } else { 2.0 }; // Degenerescência de Landau para Dirac

                enerlep += (g * qb / (4.0 * PI2))
                    * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }
        }
    }

    // Energia Total (Mésons + Bariões + Léptons + setor escuro, nulo por padrão)
    let ener = enerf + enerbar + enerlep + engine.ener_chi_kin + engine.dark_vector_energy_density();

    // Pressão via relação termodinâmica: P = sum(mu_i * n_i) - epsilon
    let mut press_sum = 0.0;
    for i in 0..8 {
        press_sum += engine.mu_b[i] * engine.nb[i];
    }
    for i in 0..2 {
        press_sum += mue * engine.nl[i]; // mu_e = mu_mu = mue
    }
    press_sum += engine.mu_chi * engine.n_chi;

    let press = press_sum - ener;

    (ener, press)
}
