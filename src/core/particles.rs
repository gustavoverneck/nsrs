// src/core/particles.rs

use crate::core::constants::PI2;
use crate::core::physics::HadronsMatter;

pub fn calculate_all_densities(engine: &mut HadronsMatter, vomega: f64, vrho: f64) {
    for i in 0..8 {
        if i >= 2 && !engine.include_hyperons {
            // Hyperon excluído: sem estados ocupados (ef = 0 zera também a
            // contribuição à energia em eos::compute).
            engine.ef_b[i] = 0.0;
            engine.kf_b_up[i].clear();
            engine.kf_b_down[i].clear();
            engine.n_b_up[i] = 0;
            engine.n_b_down[i] = 0;
            engine.rhos_b[i] = 0.0;
            engine.nb[i] = 0.0;
            continue;
        }
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

    engine.rhosb = rhos_n_p + engine.xs * rhos_h;
    engine.nbt = engine.nb.iter().sum();
}

pub fn density_baryon_neutral(
    engine: &mut HadronsMatter,
    idx: usize,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    let ef = engine.mu_b[idx]
        - (engine.xv_v[idx] * vomega)
        - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx]);

    engine.ef_b[idx] = ef;

    // ZERA OS MOMENTOS PARA EVITAR "FANTASMAS" DO NEWTON-RAPHSON
    engine.kf_b_up[idx].clear();
    engine.kf_b_down[idx].clear();

    if ef <= 0.0 {
        return (0.0, 0.0);
    }

    let m_star = engine.m_eff[idx];
    let a = engine.amm_b[idx] * engine.b;

    let mut rhos_total = 0.0;
    let mut dens_total = 0.0;

    // Spin up: a = +kappa mu_N B; spin down: a = -kappa mu_N B (mesma convenção
    // dos bárions carregados, m_eff_spin = m_landau -/+ amm b).
    for (spin_idx, a_s) in [a, -a].into_iter().enumerate() {
        if let Some(state) = neutral_amm_spin(m_star, ef, a_s) {
            rhos_total += state.scalar_density;
            dens_total += state.density;
            if spin_idx == 0 {
                engine.kf_b_up[idx].push(state.kf);
            } else {
                engine.kf_b_down[idx].push(state.kf);
            }
        }
    }
    (rhos_total, dens_total)
}

/// Um estado de spin de um bárion neutro com momento magnético anômalo.
pub(crate) struct NeutralSpinState {
    /// k_F ao longo de B (p_perp = 0): sqrt(E_F^2 - m_bar^2), m_bar = m* - a.
    pub kf: f64,
    pub density: f64,
    pub scalar_density: f64,
    pub energy: f64,
}

/// Bárion neutro, um estado de spin, com acoplamento de Pauli a = s kappa mu_N B
/// (unidades de M_N). Espectro E = sqrt(p_z^2 + (sqrt(m*^2 + p_perp^2) - a)^2).
/// Com m_bar = m* - a, k_F = sqrt(E_F^2 - m_bar^2), A = asin(m_bar/E_F) - pi/2 e
/// L = ln((E_F + k_F)/|m_bar|):
///   n   = [k_F^3/3 - (a/2)(m_bar k_F + E_F^2 A)] / 2pi^2
///   eps = [E_F^3 k_F/2 - (m_bar/4)(m_bar k_F E_F + m_bar^3 L)
///          - (a/3)(E_F m_bar k_F + m_bar^3 L) - (2/3) a E_F^3 A] / 4pi^2
///   n_s = m* [E_F k_F - m_bar^2 L] / 4pi^2
/// (cf. Broderick, Prakash & Lattimer, ApJ 537, 351 (2000)). Formas fechadas
/// derivadas e conferidas contra quadratura numérica e contra dP/dE_F = n e
/// d(eps - E_F n)/dm* = n_s; com a = 0 reduzem-se ao gás isotrópico com g = 1.
/// Retorna None se o estado está vazio (E_F <= |m_bar|).
pub(crate) fn neutral_amm_spin(m_star: f64, ef: f64, a: f64) -> Option<NeutralSpinState> {
    let m_bar = m_star - a;
    let kf2 = ef * ef - m_bar * m_bar;
    if ef <= 0.0 || kf2 <= 0.0 {
        return None;
    }
    let kf = kf2.sqrt();
    let m2 = m_bar * m_bar;
    let log = if m_bar == 0.0 { 0.0 } else { ((ef + kf) / m_bar.abs()).ln() };
    let angle = (m_bar / ef).clamp(-1.0, 1.0).asin() - 0.5 * std::f64::consts::PI;

    let density = (kf * kf2 / 3.0 - 0.5 * a * (m_bar * kf + ef * ef * angle)) / (2.0 * PI2);
    let energy = (0.5 * ef.powi(3) * kf
        - 0.25 * m_bar * (m_bar * kf * ef + m2 * m_bar * log)
        - (a / 3.0) * (ef * m_bar * kf + m2 * m_bar * log)
        - (2.0 / 3.0) * a * ef.powi(3) * angle)
        / (4.0 * PI2);
    let scalar_density = m_star * (ef * kf - m2 * log) / (4.0 * PI2);
    Some(NeutralSpinState { kf, density, scalar_density, energy })
}

fn density_baryon_charged(
    engine: &mut HadronsMatter,
    idx: usize,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    let q = engine.charges_b[idx].abs() * engine.qe;
    let b = engine.b;
    let m = engine.m_eff[idx];
    let amm = engine.amm_b[idx];

    // Deslocamento pelo fóton escuro (nulo sem mistura cinética).
    let dark_shift = engine.dark_shift_for_charge(engine.charges_b[idx]);
    let ef = engine.mu_b[idx]
        - (engine.xv_v[idx] * vomega)
        - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx])
        - dark_shift;

    engine.ef_b[idx] = ef;

    // ZERA OS MOMENTOS PARA EVITAR FANTASMAS
    engine.n_b_up[idx] = 0;
    engine.n_b_down[idx] = 0;
    engine.kf_b_up[idx].clear();
    engine.kf_b_down[idx].clear();

    if ef <= 0.0 {
        return (0.0, 0.0);
    }

    // TRATAMENTO PARA O CASO ISOTRÓPICO (B=0)
    if b == 0.0 {
        let kf2 = ef.powi(2) - m.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let dens = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos = (m / (2.0 * PI2)) * (ef * kf - m.powi(2) * ((kf + ef) / m_safe).ln());

            engine.kf_b_up[idx].push(kf);
            engine.kf_b_down[idx].push(kf);
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

    engine.kf_b_up[idx].reserve(nu_max);
    engine.kf_b_down[idx].reserve(nu_max);

    // --- Spin UP ---
    for nu in nu_start_up..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau - amm * b;

        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        engine.kf_b_up[idx].push(kf);
        n_up += 1;

        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15);

        rhos +=
            (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    // --- Spin DOWN ---
    for nu in nu_start_down..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau + amm * b;

        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        engine.kf_b_down[idx].push(kf);
        n_down += 1;

        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15);

        rhos +=
            (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    engine.n_b_up[idx] = n_up;
    engine.n_b_down[idx] = n_down;

    (rhos, dens)
}

pub fn density_lepton(engine: &mut HadronsMatter, idx: usize) -> (f64, f64) {
    // Energia de Fermi efetiva: mu_e menos o deslocamento do fóton escuro
    // para carga -1 (nulo sem mistura cinética).
    let mue = engine.mue - engine.dark_shift_for_charge(-1.0);

    // ZERA OS ESTADOS PARA EVITAR FANTASMAS
    engine.n_l[idx] = 0;
    engine.f_l[idx].clear();
    engine.ef_l[idx] = mue;

    if mue <= 0.0 {
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
            let dens_val = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos_val = (m / (2.0 * PI2)) * (mue * kf - m.powi(2) * ((kf + mue) / m_safe).ln());

            engine.f_l[idx].push(kf);
            engine.n_l[idx] = 1;

            return (rhos_val, dens_val);
        }
        return (0.0, 0.0);
    }

    let nu_max_approx = (mue.powi(2) - m.powi(2)) / (2.0 * q * b);
    let nu_max = if nu_max_approx > 0.0 {
        (nu_max_approx.floor() as usize + 1).min(engine.max_landau_limit)
    } else {
        0
    };

    engine.f_l[idx].reserve(nu_max);

    for nu in 0..nu_max {
        let m_landau_2 = m.powi(2) + 2.0 * q * b * nu as f64;
        let kf2 = mue.powi(2) - m_landau_2;

        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        let g = if nu == 0 { 1.0 } else { 2.0 };

        engine.f_l[idx].push(kf);

        let factor = (g * q * b) / (2.0 * PI2);
        let m_safe = m_landau_2.sqrt().max(1e-15);

        rhos += factor * m * ((kf + mue) / m_safe).ln();
        dens += factor * kf;

        n_occupied += 1;
    }

    engine.n_l[idx] = n_occupied;

    (rhos, dens)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// n, eps, n_s de um estado de spin por quadratura. Com u = sqrt(m^2 + p_perp^2) - a,
    /// p_perp dp_perp = (u + a) du e a integral em p_z é analítica; u = E - (E - m_bar) s^2
    /// remove a raiz em u = E. Requer m_bar > 0.
    fn quadrature(m: f64, ef: f64, a: f64) -> (f64, f64, f64) {
        let m_bar = m - a;
        let steps = 20_000;
        let h = 1.0 / steps as f64;
        let (mut n, mut eps, mut ns) = (0.0, 0.0, 0.0);
        for k in 0..=steps {
            let s = k as f64 * h;
            let w = if k == 0 || k == steps { 1.0 } else if k % 2 == 1 { 4.0 } else { 2.0 };
            let u = ef - (ef - m_bar) * s * s;
            let jac = 2.0 * (ef - m_bar) * s;
            let kz = (ef * ef - u * u).max(0.0).sqrt();
            let log = ((kz + ef) / u).ln();
            n += w * jac * (u + a) * 2.0 * kz;
            eps += w * jac * (u + a) * (kz * ef + u * u * log);
            ns += w * jac * 2.0 * m * u * log;
        }
        let norm = h / 3.0 / (4.0 * PI2);
        (n * norm, eps * norm, ns * norm)
    }

    const CASES: [(f64, f64, f64); 5] =
        [(0.65, 0.9, 0.01), (0.65, 0.9, -0.01), (0.6, 1.1, 0.05), (0.7, 0.75, -0.04), (0.3, 0.5, 0.08)];

    #[test]
    fn neutral_amm_closed_forms_match_quadrature() {
        for (m, ef, a) in CASES {
            let state = neutral_amm_spin(m, ef, a).unwrap();
            let (n, eps, ns) = quadrature(m, ef, a);
            for (name, closed, num) in
                [("n", state.density, n), ("eps", state.energy, eps), ("n_s", state.scalar_density, ns)]
            {
                assert!(((closed - num) / num).abs() < 1e-9, "{name}({m}, {ef}, {a}): {closed} vs {num}");
            }
        }
    }

    #[test]
    fn neutral_amm_without_field_is_isotropic_gas() {
        let (m, ef) = (0.65, 0.9);
        let state = neutral_amm_spin(m, ef, 0.0).unwrap();
        let kf = (ef * ef - m * m).sqrt();
        let log = ((ef + kf) / m).ln();
        let n = kf.powi(3) / (6.0 * PI2);
        let eps = (ef.powi(3) * kf / 2.0 - (m / 4.0) * (m * kf * ef + m.powi(3) * log)) / (4.0 * PI2);
        let ns = m * (ef * kf - m * m * log) / (4.0 * PI2);
        assert!((state.density / n - 1.0).abs() < 1e-14);
        assert!((state.energy / eps - 1.0).abs() < 1e-14);
        assert!((state.scalar_density / ns - 1.0).abs() < 1e-14);
    }

    #[test]
    fn neutral_amm_is_thermodynamically_consistent() {
        // P = E_F n - eps: dP/dE_F = n e d(eps - E_F n)/dm* = n_s.
        let h = 1e-6;
        for (m, ef, a) in CASES {
            let at = |m: f64, ef: f64| neutral_amm_spin(m, ef, a).unwrap();
            let p = |ef: f64| ef * at(m, ef).density - at(m, ef).energy;
            let omega = |m: f64| at(m, ef).energy - ef * at(m, ef).density;
            let dp = (p(ef + h) - p(ef - h)) / (2.0 * h);
            let domega = (omega(m + h) - omega(m - h)) / (2.0 * h);
            let state = at(m, ef);
            assert!((dp / state.density - 1.0).abs() < 1e-7, "dP/dE_F = {dp}, n = {}", state.density);
            assert!((domega / state.scalar_density - 1.0).abs() < 1e-7, "{domega} vs {}", state.scalar_density);
        }
    }
}
