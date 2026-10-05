#![allow(dead_code)]

pub const PI: f64 = std::f64::consts::PI;
pub const PI2: f64 = std::f64::consts::PI * std::f64::consts::PI;

pub const HBAR_C: f64 = 197.3269804; // MeV * fm
pub const M_NUCLEON: f64 = 939.56534623; //938.9187; // Massa média do nucleon (MeV)
pub const N0: f64 = 0.153; // Densidade de saturação nuclear (fm^-3)
pub const MS_TOV: f64 = 5660.57; // fator para converter MeV/fm³ -> M_sun/km³
pub const GMS_TOV: f64 = 1.47556; // GM_sun/c²  [km]
pub const QE: f64 = 0.302822120868846; // (4.0 * PI / 137.0).sqrt()

// massas das partículas
pub const MN: f64 = 939.565; // Massa do nêutron (MeV)
pub const MP: f64 = 938.272; // Massa do próton (MeV)
pub const ME: f64 = 0.510998; // Massa do elétron (MeV)

pub const MEV_FM3_TO_MSUN_KM3: f64 = 8.9653e-7;
pub const G_C2: f64 = 1.4766; // km / M_sol

// Baryon masses
pub const MB: [f64; 8] = [
    939.56534623 / M_NUCLEON,
    938.272081323 / M_NUCLEON,
    1116.0 / M_NUCLEON,
    1193.0 / M_NUCLEON,
    1193.0 / M_NUCLEON,
    1193.0 / M_NUCLEON,
    1318.0 / M_NUCLEON,
    1318.0 / M_NUCLEON,
];

pub const ML: [f64; 2] = [0.511 / M_NUCLEON, 105.66 / M_NUCLEON];

// Meson Massses
pub const MS: f64 = 400.0 / M_NUCLEON; // Scalar meson (sigma)
pub const MV: f64 = 783.0 / M_NUCLEON; // Vector meson (Omega)
pub const MR: f64 = 770.0 / M_NUCLEON; // Isovector meson (Rho)

// Magneton nuclear mu_N = e/(2 m_p) em unidades de M_N: a energia de Pauli é
// kappa * RNCM * b, com b o campo nas unidades do código (b = B/B_ce * BCE).
pub const RNCM: f64 = QE * M_NUCLEON / (2.0 * MP);

// Momentos magnéticos anômalos kappa_b = mu_b/mu_N - q_b m_p/m_b (PDG), na ordem
// [n, p, Lambda, Sigma-, Sigma0, Sigma+, Xi-, Xi0]. O valor de Sigma0 é a média de
// Sigma+ e Sigma- (não medido). Os elétrons e múons ficam sem AMM. Só entram com
// `HadronsMatter::with_anomalous_moments(true)`; o padrão é kappa = 0.
pub const KAPPA_B: [f64; 8] = [
    -1.913,
    2.793 - MP / 938.272081323,
    -0.613,
    -1.160 + MP / 1193.0,
    0.649,
    2.458 - MP / 1193.0,
    -0.650 + MP / 1318.0,
    -1.250,
];

// BDD Constants
pub const BCE: f64 = ML[0] * ML[0] / QE;
pub const BCE_G: f64 = 4.41e13;
pub const BDD_BETAA: f64 = 1e-2;
pub const BDD_ALPHAA: f64 = 3.0;

/// Number of EOS diagnostic columns. Columns 21..=33 describe the fermionic
/// dark sector; non-dark engines leave them at zero.
pub const RESULTS_SIZE: usize = 34;
pub const DATA_SIZE: usize = RESULTS_SIZE + 3;

pub const MAX_LANDAU_LIMIT: usize = 20_000;
