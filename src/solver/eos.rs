// src/solver/eos.rs

use crate::solver::physics::PhysicsEngine;
use crate::solver::constants::PI2;
use crate::solver::particles::dg;

pub fn compute(engine: &PhysicsEngine, mue: f64, vsigma: f64, vomega: f64, vrho: f64) -> (f64, f64) {
    // energia dos mésons
    let enerf = (vsigma / engine.model.gs).powi(2) / 2.0
        + (vomega / engine.model.gv).powi(2) / 2.0
        + (vrho / engine.model.gr).powi(2) / 2.0
        + engine.model.rb * vsigma.powi(3) / 3.0
        + engine.model.rc * vsigma.powi(4) / 4.0
        + engine.model.rxi * vomega.powi(4) / 4.0;

    let mut enerbar = 0.0;

    // ----- NÊUTRON -----
    let kfn = (engine.efn.powi(2) - engine.mns.powi(2)).sqrt().max(0.0);
    if kfn > 0.0 {
        enerbar += (1.0 / (2.0 * PI2)) * (
            engine.efn.powi(3) * kfn / 2.0
            - (engine.mns / 4.0) * (engine.mns * kfn * engine.efn + engine.mns.powi(3) * ((kfn + engine.efn) / engine.mns).ln())
        );
    }

    // ----- PRÓTON -----
    for nu in 0..engine.npu {
        let kf = engine.fpu[nu];
        let m_eff = (engine.mns.powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        enerbar += (engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efp * kf + m_eff.powi(2) * ((kf + engine.efp) / m_eff).ln());
    }
    for nu in 0..engine.npd {
        let kf = engine.fpd[nu];
        let m_eff = (engine.mns.powi(2) + 2.0 * engine.qe * engine.b * (nu+1) as f64).sqrt();
        enerbar += (engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efp * kf + m_eff.powi(2) * ((kf + engine.efp) / m_eff).ln());
    }

    // ----- LAMBDA 0 -----
    let kfl0 = engine.rkfl0;
    if kfl0 > 0.0 {
        enerbar += (1.0 / (2.0 * PI2)) * (
            engine.efl0.powi(3) * kfl0 / 2.0
            - (engine.ml0s / 4.0) * (engine.ml0s * kfl0 * engine.efl0 + engine.ml0s.powi(3) * ((kfl0 + engine.efl0) / engine.ml0s).ln())
        );
    }

    // ----- SIGMA 0 -----
    let kfs0 = engine.rkfs0;
    if kfs0 > 0.0 {
        enerbar += (1.0 / (2.0 * PI2)) * (
            engine.efs0.powi(3) * kfs0 / 2.0
            - (engine.ms0s / 4.0) * (engine.ms0s * kfs0 * engine.efs0 + engine.ms0s.powi(3) * ((kfs0 + engine.efs0) / engine.ms0s).ln())
        );
    }

    // ----- XI 0 -----
    let kfx0 = engine.rkfx0;
    if kfx0 > 0.0 {
        enerbar += (1.0 / (2.0 * PI2)) * (
            engine.efx0.powi(3) * kfx0 / 2.0
            - (engine.mx0s / 4.0) * (engine.mx0s * kfx0 * engine.efx0 + engine.mx0s.powi(3) * ((kfx0 + engine.efx0) / engine.mx0s).ln())
        );
    }

    // ----- SIGMA- (carregado) -----
    for nu in 0..engine.nsm {
        let kf = engine.fsm[nu];
        let m_eff = (engine.msms.powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        let g = dg(nu);
        enerbar += (g * engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efsm * kf + m_eff.powi(2) * ((kf + engine.efsm) / m_eff).ln());
    }

    // ----- SIGMA+ (carregado) -----
    for nu in 0..engine.nsp {
        let kf = engine.fsp[nu];
        let m_eff = (engine.msps.powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        let g = dg(nu);
        enerbar += (g * engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efsp * kf + m_eff.powi(2) * ((kf + engine.efsp) / m_eff).ln());
    }

    // ----- XI- (carregado) -----
    for nu in 0..engine.nxm {
        let kf = engine.fxm[nu];
        let m_eff = (engine.mxms.powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        let g = dg(nu);
        enerbar += (g * engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efxm * kf + m_eff.powi(2) * ((kf + engine.efxm) / m_eff).ln());
    }

    // ----- LÉPTONS -----
    let mut enerlep = 0.0;
    for nu in 0..engine.ne {
        let kf = engine.fe[nu];
        let m_eff = (engine.ml[0].powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        let g = dg(nu);
        enerlep += (g * engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efe * kf + m_eff.powi(2) * ((kf + engine.efe) / m_eff).ln());
    }
    for nu in 0..engine.nu {
        let kf = engine.fmu[nu];
        let m_eff = (engine.ml[1].powi(2) + 2.0 * engine.qe * engine.b * nu as f64).sqrt();
        let g = dg(nu);
        enerlep += (g * engine.qe * engine.b / (4.0 * PI2)) *
            (engine.efmu * kf + m_eff.powi(2) * ((kf + engine.efmu) / m_eff).ln());
    }

    let ener = enerf + enerbar + enerlep;

    // pressão (a sua fórmula estava matematicamente idêntica à do Fortran, então podemos mantê-la)
    let press = engine.mun * (engine.nb[0] + engine.nb[2] + engine.nb[4] + engine.nb[7])
        + (engine.mun - mue) * (engine.nb[1] + engine.nb[5])
        + mue * (engine.nl[0] + engine.nl[1])
        + (engine.mun + mue) * (engine.nb[3] + engine.nb[6])
        - ener;

    (ener, press)
}