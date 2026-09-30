// src/core/nuclear.rs
//
// Propriedades da matéria nuclear simétrica na saturação, extraídas
// diretamente do motor `HadronsMatter` (mesmas equações de campo e energia
// usadas nas estrelas).

use crate::core::constants::{HBAR_C, M_NUCLEON, MB, PI2};
use crate::core::model::ModelParams;
use crate::core::physics::HadronsMatter;
use nalgebra::{Matrix3, Vector3};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturationProperties {
    /// Densidade de saturação (fm^-3).
    pub n0: f64,
    /// Energia por nucleon, medida da massa nucleônica média (MeV).
    pub energy_per_nucleon: f64,
    /// Incompressibilidade K = 9 n0^2 d^2(E/A)/dn^2 (MeV).
    pub incompressibility: f64,
    /// Energia de simetria J (MeV).
    pub symmetry_energy: f64,
    /// Inclinação L = 3 n0 dE_sym/dn (MeV).
    pub symmetry_slope: f64,
    /// Massa efetiva de Dirac M*/M.
    pub effective_mass: f64,
}

/// Estado de matéria simétrica (mu_e = 0, sem neutralidade) a mu_n fixo.
struct SymmetricPoint {
    fields: [f64; 3],
    /// n_B em unidades de M_N^3.
    n_b: f64,
    /// Pressão em M_N^4.
    pressure: f64,
}

fn solve_symmetric(engine: &mut HadronsMatter, mun: f64, guess: [f64; 3]) -> Option<SymmetricPoint> {
    engine.mun = mun;
    let residual = |engine: &mut HadronsMatter, y: &Vector3<f64>| {
        let f = engine.funcv(&[0.0, y[0], y[1], y[2]]);
        Vector3::new(f[0], f[1], f[2])
    };
    let mut y = Vector3::from(guess);
    for _ in 0..200 {
        let f = residual(engine, &y);
        if f.norm() < 1e-13 {
            let _ = engine.funcv(&[0.0, y[0], y[1], y[2]]);
            let (_, pressure) = crate::core::eos::compute(engine, 0.0, y[0], y[1], y[2]);
            return Some(SymmetricPoint {
                fields: [y[0], y[1], y[2]],
                n_b: engine.nbt,
                pressure,
            });
        }
        let mut jac = Matrix3::zeros();
        for i in 0..3 {
            let h = 1e-7 * (y[i].abs() + 1e-3);
            let mut yh = y;
            yh[i] += h;
            jac.set_column(i, &((residual(engine, &yh) - f) / h));
        }
        let step = jac.lu().solve(&(-f))?;
        let mut alpha = 1.0;
        while alpha > 1e-6 && residual(engine, &(y + alpha * step)).norm() >= f.norm() {
            alpha *= 0.5;
        }
        y += alpha * step;
    }
    None
}

/// Energia de simetria em RMF (sem méson delta), em unidades de M_N:
/// E_sym = k_F^2 / (6 E_F*) + C_rho,ef^2 n / 8, com o acoplamento omega-rho
/// incluído na massa efetiva do rho: C_rho,ef^2 = C_rho^2 / (1 + 2 Lambda_v
/// C_rho^2 v_omega^2).
fn symmetry_energy(model: &ModelParams, point: &SymmetricPoint) -> f64 {
    let kf = (1.5 * PI2 * point.n_b).cbrt();
    let m_star = 0.5 * (MB[0] + MB[1]) - point.fields[0];
    let ef = (kf * kf + m_star * m_star).sqrt();
    let c_rho2 = model.gr.powi(2);
    let c_eff2 = c_rho2 / (1.0 + 2.0 * model.lambda_v * c_rho2 * point.fields[1].powi(2));
    kf * kf / (6.0 * ef) + c_eff2 * point.n_b / 8.0
}

/// Propriedades de saturação da matéria simétrica: P(mu) = 0 no ramo denso.
pub fn saturation_properties(model: ModelParams) -> Option<SaturationProperties> {
    let mut engine = HadronsMatter::new(model, 0.0);

    // Continuação descendo pelo ramo denso até a pressão trocar de sinal.
    let (mut mu, mut guess) = (1.2, [0.4, 0.4, 0.0]);
    let (mut mu_hi, mut g_hi) = (mu, guess);
    let mut mu_lo = loop {
        let point = solve_symmetric(&mut engine, mu, guess)?;
        if point.n_b <= 1e-5 {
            return None;
        }
        guess = point.fields;
        if point.pressure < 0.0 {
            break mu;
        }
        (mu_hi, g_hi) = (mu, point.fields);
        mu -= if mu > 1.0 { 5e-3 } else { 5e-4 };
        if mu < 0.9 {
            return None;
        }
    };
    for _ in 0..50 {
        let mid = 0.5 * (mu_hi + mu_lo);
        let point = solve_symmetric(&mut engine, mid, g_hi)?;
        if point.pressure > 0.0 {
            (mu_hi, g_hi) = (mid, point.fields);
        } else {
            mu_lo = mid;
        }
    }

    let mu_sat = mu_hi;
    let sat = solve_symmetric(&mut engine, mu_sat, g_hi)?;
    let h = 1e-5;
    let plus = solve_symmetric(&mut engine, mu_sat + h, g_hi)?;
    let minus = solve_symmetric(&mut engine, mu_sat - h, g_hi)?;

    let to_fm3 = (M_NUCLEON / HBAR_C).powi(3);
    let n0 = sat.n_b * to_fm3;
    // P = 0: E/A = mu_n; K = 9 dP/dn = 9 n / (dn/dmu) (P = n^2 d(E/A)/dn).
    let dn_dmu = (plus.n_b - minus.n_b) / (2.0 * h);
    let incompressibility = 9.0 * sat.n_b / dn_dmu * M_NUCLEON;
    let symmetry_energy_sat = symmetry_energy(&model, &sat) * M_NUCLEON;
    let d_esym_dn = (symmetry_energy(&model, &plus) - symmetry_energy(&model, &minus))
        / (plus.n_b - minus.n_b);
    let symmetry_slope = 3.0 * sat.n_b * d_esym_dn * M_NUCLEON;

    Some(SaturationProperties {
        n0,
        energy_per_nucleon: (mu_sat - 0.5 * (MB[0] + MB[1])) * M_NUCLEON,
        incompressibility,
        symmetry_energy: symmetry_energy_sat,
        symmetry_slope,
        effective_mass: 0.5 * (MB[0] + MB[1]) - sat.fields[0],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::model::{FSU2, GM1, GM3};

    /// Chen & Piekarewicz, PRC 90, 044305 (2014), Tabela: n0 = 0.1505(7),
    /// E/A = -16.28(2), M*/M = 0.593(4), K = 238.0(2.8), J = 37.62(1.11),
    /// L = 112.8(16.1).
    #[test]
    fn fsu2_bulk_properties_match_chen_piekarewicz() {
        let s = saturation_properties(FSU2).unwrap();
        assert!((s.n0 - 0.1505).abs() < 0.0015, "n0 = {}", s.n0);
        assert!((s.energy_per_nucleon + 16.28).abs() < 0.1, "E/A = {}", s.energy_per_nucleon);
        assert!((s.effective_mass - 0.593).abs() < 0.005, "M* = {}", s.effective_mass);
        assert!((s.incompressibility - 238.0).abs() < 2.8, "K = {}", s.incompressibility);
        assert!((s.symmetry_energy - 37.62).abs() < 1.11, "J = {}", s.symmetry_energy);
        assert!((s.symmetry_slope - 112.8).abs() < 2.0, "L = {}", s.symmetry_slope);
    }

    /// Glendenning & Moszkowski, PRL 67, 2414 (1991): n0 = 0.153, E/A = -16.3,
    /// K = 300 (GM1) e 240 (GM3), M*/M = 0.70 e 0.78, a_sym = 32.5 MeV.
    #[test]
    fn gm_bulk_properties_match_glendenning_moszkowski() {
        for (model, k, m_star) in [(GM1, 300.0, 0.70), (GM3, 240.0, 0.78)] {
            let s = saturation_properties(model).unwrap();
            assert!((s.n0 - 0.153).abs() < 0.002, "n0 = {}", s.n0);
            assert!((s.energy_per_nucleon + 16.3).abs() < 0.1, "E/A = {}", s.energy_per_nucleon);
            assert!((s.incompressibility - k).abs() < 0.01 * k, "K = {}", s.incompressibility);
            assert!((s.effective_mass - m_star).abs() < 0.005, "M* = {}", s.effective_mass);
            assert!((s.symmetry_energy - 32.5).abs() < 0.3, "J = {}", s.symmetry_energy);
        }
    }
}
