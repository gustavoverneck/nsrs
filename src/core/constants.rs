#![allow(dead_code)]

pub const PI: f64 = std::f64::consts::PI;
pub const PI2: f64 = std::f64::consts::PI * std::f64::consts::PI;

pub const HBAR_C: f64 = 197.3269804; // MeV * fm
pub const M_NUCLEON: f64 = 938.9187; // Massa média do nucleon (MeV)
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

pub const ML: [f64; 2] = [
            0.511 / M_NUCLEON,
            105.66 / M_NUCLEON,
];

// Meson Massses
pub const MS: f64 = 400.0 / M_NUCLEON; // Scalar meson (sigma)
pub const MV: f64 = 783.0 / M_NUCLEON; // Vector meson (Omega)
pub const MR: f64 = 770.0 / M_NUCLEON; // Isovector meson (Rho)

pub const BCE: f64 = ML[0] * ML[0] / QE;

// Ative o magneton nuclear removendo o * 0.0
pub const RNCM: f64 = QE / 2.0; 

// Valores baseados no Particle Data Group (PDG) para kappa = mu - q
pub const AMMN: f64  = RNCM * -1.913; // Neutrão
pub const AMMP: f64  = RNCM * 1.793; // Protão
pub const AMML0: f64 = RNCM * -0.613; // Lambda0
pub const AMMSM: f64 = RNCM * -0.160; // Sigma- (-1.16 - (-1))
pub const AMMS0: f64 = RNCM * 0.649; // Sigma0 (Teórico)
pub const AMMSP: f64 = RNCM * 1.458; // Sigma+ (2.458 - 1)
pub const AMMXM: f64 = RNCM * 0.349; // Xi-    (-0.65 - (-1))
pub const AMMX0: f64 = RNCM * -1.250; // Xi0
