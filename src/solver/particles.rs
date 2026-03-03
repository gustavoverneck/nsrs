// solver/particles.rs

use crate::solver::physics::PhysicsEngine;
use crate::solver::constants::PI2;

// Fator de degenerescência para níveis de Landau (g=1 para n=0, g=2 para n>0)
pub fn dg(nu: usize) -> f64 {
    if nu == 0 { 1.0 } else { 2.0 }
}

pub fn calculate_all_densities(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) {
    let (rhosn, densn) = density_neutron(engine, vomega, vrho);
    let (rhosp, densp) = density_proton(engine, vomega, vrho);
    let (rhosl0, densl0) = density_lambda0(engine, vomega);
    let (rhossm, denssm) = density_sigma_minus(engine, vomega, vrho);
    let (rhoss0, denss0) = density_sigma0(engine, vomega);
    let (rhossp, denssp) = density_sigma_plus(engine, vomega, vrho);
    let (rhosxm, densxm) = density_xi_minus(engine, vomega, vrho);
    let (rhosx0, densx0) = density_xi0(engine, vomega, vrho);
    let (rhose, dense) = density_electron(engine);
    let (rhosmu, densmu) = density_muon(engine);

    // total escalar dos híperons
    let rhosh = rhosl0 + rhossm + rhoss0 + rhossp + rhosxm + rhosx0;
    let cahs = engine.xs * rhosh;
    engine.rhosb = rhosn + rhosp + cahs;

    engine.nb = [densn, densp, densl0, denssm, denss0, denssp, densxm, densx0];
    engine.nl = [dense, densmu];
    engine.nbt = engine.nb.iter().sum();
}

// ---------- Nêutron ----------
fn density_neutron(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let efn = engine.mun - vomega + vrho / 2.0;
    engine.efn = efn;
    let kf2 = efn.powi(2) - engine.mns.powi(2);
    if kf2 <= 0.0 {
        return (0.0, 0.0);
    }
    let kf = kf2.sqrt();
    let rhos = (engine.mns / (2.0 * PI2)) * (efn * kf - engine.mns.powi(2) * ((kf + efn) / engine.mns).ln());
    let dens = kf.powi(3) / (3.0 * PI2);
    (rhos, dens)
}

// ---------- Próton (carregado) ----------
pub fn density_proton(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let ef = engine.mup - vomega - vrho / 2.0;
    engine.efp = ef;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.mns;
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut npu = 0;
    let mut npd = 0;

    // spin up (ν = 0,1,2,...)
    for nu in 0..engine.fpu.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = ef.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fpu[nu] = kf;
        
        rhos += (q * b * m / (2.0 * PI2)) * ((kf + ef) / m_eff).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
        npu += 1;
    }

    // spin down (ν = 1,2,...)
    for nu in 1..engine.fpd.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = ef.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fpd[nu-1] = kf;
        
        rhos += (q * b * m / (2.0 * PI2)) * ((kf + ef) / m_eff).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
        npd += 1;
    }

    engine.npu = npu;
    engine.npd = npd;
    (rhos, dens)
}


// ---------- Elétron ----------
pub fn density_electron(engine: &mut PhysicsEngine) -> (f64, f64) {
    let mue = engine.mue;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.ml[0];
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut ne = 0;

    for nu in 0..engine.fe.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = mue.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fe[nu] = kf;
        let g = dg(nu);
        rhos += (g * q * b * m / (2.0 * PI2)) * ((kf + mue) / m_eff).ln();
        dens += (g * q * b / (2.0 * PI2)) * kf;
        ne += 1;
    }
    engine.ne = ne;
    if ne > 0 {
        engine.efe = (engine.ml[0].powi(2) + engine.fe[0].powi(2)).sqrt();
    } else {
        engine.efe = 0.0;
    }
    (rhos, dens)
}


// ---------- Muon ----------
pub fn density_muon(engine: &mut PhysicsEngine) -> (f64, f64) {
    let mue = engine.mue;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.ml[1];
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut nu = 0;

    for nu_idx in 0..engine.fmu.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu_idx as f64).sqrt();
        let kf2 = mue.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fmu[nu_idx] = kf;
        let g = dg(nu_idx);
        rhos += (g * q * b * m / (2.0 * PI2)) * ((kf + mue) / m_eff).ln();
        dens += (g * q * b / (2.0 * PI2)) * kf;
        nu += 1;
    }
    engine.nu = nu;
    if nu > 0 {
        engine.efmu = (engine.ml[1].powi(2) + engine.fmu[0].powi(2)).sqrt();
    } else {
        engine.efmu = 0.0;
    }
    (rhos, dens)
}

// ---------- Lambda0 (neutro) ----------
pub fn density_lambda0(engine: &mut PhysicsEngine, vomega: f64) -> (f64, f64) {
    let ef = engine.mul0 - engine.xv * vomega;
    let kf2 = ef.powi(2) - engine.ml0s.powi(2);
    if kf2 <= 0.0 {
        engine.rkfl0 = 0.0;
        engine.efl0 = 0.0;
        return (0.0, 0.0);
    }
    let kf = kf2.sqrt();
    engine.rkfl0 = kf;
    engine.efl0 = (engine.ml0s.powi(2) + kf.powi(2)).sqrt();
    let rhos = (engine.ml0s / (2.0 * PI2)) * (ef * kf - engine.ml0s.powi(2) * ((kf + ef) / engine.ml0s).ln());
    let dens = kf.powi(3) / (3.0 * PI2);
    (rhos, dens)
}

// ---------- Sigma0 (neutro) ----------
pub fn density_sigma0(engine: &mut PhysicsEngine, vomega: f64) -> (f64, f64) {
    let ef = engine.mus0 - engine.xv * vomega;
    let kf2 = ef.powi(2) - engine.ms0s.powi(2);
    if kf2 <= 0.0 {
        engine.rkfs0 = 0.0;
        engine.efs0 = 0.0;
        return (0.0, 0.0);
    }
    let kf = kf2.sqrt();
    engine.rkfs0 = kf;
    engine.efs0 = (engine.ms0s.powi(2) + kf.powi(2)).sqrt();
    let rhos = (engine.ms0s / (2.0 * PI2)) * (ef * kf - engine.ms0s.powi(2) * ((kf + ef) / engine.ms0s).ln());
    let dens = kf.powi(3) / (3.0 * PI2);
    (rhos, dens)
}

// ---------- Xi0 (neutro) ----------
pub fn density_xi0(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let ef = engine.mux0 - engine.xv * vomega - engine.xv * vrho / 2.0;
    let kf2 = ef.powi(2) - engine.mx0s.powi(2);
    if kf2 <= 0.0 {
        engine.rkfx0 = 0.0;
        engine.efx0 = 0.0;
        return (0.0, 0.0);
    }
    let kf = kf2.sqrt();
    engine.rkfx0 = kf;
    engine.efx0 = (engine.mx0s.powi(2) + kf.powi(2)).sqrt();
    let rhos = (engine.mx0s / (2.0 * PI2)) * (ef * kf - engine.mx0s.powi(2) * ((kf + ef) / engine.mx0s).ln());
    let dens = kf.powi(3) / (3.0 * PI2);
    (rhos, dens)
}

// ---------- Sigma- (carregado) ----------
pub fn density_sigma_minus(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let ef = engine.musm - engine.xv * vomega + engine.xv * vrho;
    engine.efsm = ef;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.msms;
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut nsm = 0;

    for nu in 0..engine.fsm.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = ef.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fsm[nu] = kf;
        let g = dg(nu);
        rhos += (g * q * b * m / (2.0 * PI2)) * ((kf + ef) / m_eff).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
        nsm += 1;
    }
    engine.nsm = nsm;
    (rhos, dens)
}

// ---------- Sigma+ (carregado) ----------
pub fn density_sigma_plus(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let ef = engine.musp - engine.xv * vomega - engine.xv * vrho;
    engine.efsp = ef;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.msps;
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut nsp = 0;

    for nu in 0..engine.fsp.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = ef.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fsp[nu] = kf;
        let g = dg(nu);
        rhos += (g * q * b * m / (2.0 * PI2)) * ((kf + ef) / m_eff).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
        nsp += 1;
    }
    engine.nsp = nsp;
    (rhos, dens)
}

// ---------- Xi- (carregado) ----------
pub fn density_xi_minus(engine: &mut PhysicsEngine, vomega: f64, vrho: f64) -> (f64, f64) {
    let ef = engine.muxm - engine.xv * vomega + engine.xv * vrho / 2.0;
    engine.efxm = ef;
    let b = engine.b;
    let q = engine.qe;
    let m = engine.mxms;
    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut nxm = 0;

    for nu in 0..engine.fxm.len() {
        let m_eff = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let kf2 = ef.powi(2) - m_eff.powi(2);
        if kf2 <= 0.0 { break; }
        let kf = kf2.sqrt();
        engine.fxm[nu] = kf;
        let g = dg(nu);
        rhos += (g * q * b * m / (2.0 * PI2)) * ((kf + ef) / m_eff).ln();
        dens += (g * q * b / (2.0 * PI2)) * kf;
        nxm += 1;
    }
    engine.nxm = nxm;
    (rhos, dens)
}