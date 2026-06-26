// solver/model.rs

use crate::core::constants::{HBAR_C, M_NUCLEON};

#[derive(Debug, Clone, Copy)]
pub struct ModelParams {
    pub gs: f64,
    pub gv: f64,
    pub gr: f64,
    pub rb: f64,
    pub rc: f64,
    pub rxi: f64,
    pub lambda_v: f64,
}

pub const GM1: ModelParams = ModelParams {
    gs: (3.43292877875 / HBAR_C) * M_NUCLEON, // sqrt(11.785)
    gv: (2.67357438647 / HBAR_C) * M_NUCLEON, // sqrt(7.148)
    gr: (2.1 / HBAR_C) * M_NUCLEON,           // sqrt(4.41)
    rb: 0.002948,
    rc: -0.001071,
    rxi: 0.0,
    lambda_v: 0.0,
};

pub const GM3: ModelParams = ModelParams {
    gs: (3.15071420475 / HBAR_C) * M_NUCLEON, // sqrt(9.927)
    gv: (2.19544984001 / HBAR_C) * M_NUCLEON, // sqrt(4.820)
    gr: (2.18883530673 / HBAR_C) * M_NUCLEON, // sqrt(4.791)
    rb: 0.008659,
    rc: -0.002421,
    rxi: 0.0,
    lambda_v: 0.0,
};

/// FSUGold2 / FSU2 from Chen & Piekarewicz, Phys. Rev. C 90, 044305 (2014).
///
/// The table gives m_s, m_v, m_rho, g_s^2, g_v^2, g_rho^2, kappa, lambda,
/// zeta, and Lambda_v.  The solver uses the same normalized field variables as
/// GM1/GM3, so kappa and lambda are converted into Boguta-Bodmer b and c by
/// b = kappa / (2 M_N) and c = lambda / 6.  The omega self-coupling is stored
/// as zeta / 6 so that the omega field equation contains -rxi * omega^3.
pub const FSU2: ModelParams = ModelParams {
    gs: (10.39684086634012 / 497.479) * M_NUCLEON,
    gv: (13.556891236563049 / 782.500) * M_NUCLEON,
    gr: (8.970261980566677 / 763.000) * M_NUCLEON,
    rb: 3.0029 / (2.0 * M_NUCLEON),
    rc: -0.000533 / 6.0,
    rxi: 0.0256 / 6.0,
    lambda_v: 0.000823,
};
