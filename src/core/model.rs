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
}

pub const GM1: ModelParams = ModelParams {
    gs: (3.43292877875 / HBAR_C) * M_NUCLEON, // sqrt(11.785)
    gv: (2.67357438647 / HBAR_C) * M_NUCLEON, // sqrt(7.148)
    gr: (2.1 / HBAR_C) * M_NUCLEON,           // sqrt(4.41)
    rb: 0.002948,
    rc: -0.001071,
    rxi: 0.0,
};

pub const GM3: ModelParams = ModelParams {
    gs: (3.15071420475 / HBAR_C) * M_NUCLEON, // sqrt(9.927)
    gv: (2.19544984001 / HBAR_C) * M_NUCLEON, // sqrt(4.820)
    gr: (2.18883530673 / HBAR_C) * M_NUCLEON, // sqrt(4.791)
    rb: 0.008659,
    rc: -0.002421,
    rxi: 0.0,
};