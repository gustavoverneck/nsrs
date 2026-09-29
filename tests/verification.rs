//! Nível 1 — verificação numérica contra soluções exatas e identidades
//! termodinâmicas. Estes testes não dependem de valores da literatura: uma
//! falha aqui indica erro de implementação, não de modelo físico.

use nsrs::constants::{G_C2, M_NUCLEON, MEV_FM3_TO_MSUN_KM3, N0, RESULTS_SIZE};
use nsrs::core::model::ModelParams;
use nsrs::core::tov_solver::integrate_star;
use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver};
use std::f64::consts::PI;

type Row = [f64; RESULTS_SIZE];

const COL_NB: usize = 0;
const COL_EPS: usize = 1;
const COL_P: usize = 2;
const COL_NE: usize = 3;
const COL_NMU: usize = 4;
const COL_NP: usize = 6;
const COL_NSM: usize = 8;
const COL_NSP: usize = 10;
const COL_NXM: usize = 11;
const COL_MU_N: usize = 17;
const COL_EPS_MAG: usize = 19;

fn solve_eos(model: ModelParams, b_gauss: f64) -> Vec<Row> {
    Solver::new(EngineMode::Hadrons(HadronsMatter::new(model, b_gauss))).solve()
}

/// Pressão da matéria sem a contribuição macroscópica do campo. Com a
/// topologia padrão (anisotrópica), P_mag = eps_mag (coluna 19).
fn matter_pressure(row: &Row) -> f64 {
    row[COL_P] - row[COL_EPS_MAG]
}

/// Estrela de densidade uniforme: solução exata de Schwarzschild para o
/// interior. Para eps constante e razão p = P_c/eps,
///   sqrt(1 - 2GM/R) = (1 + p)/(1 + 3p),   M = 4 pi R^3 eps / 3,
/// e a massa própria é eps * V_própria com
///   V = 2 pi a^3 [asin(x) - x sqrt(1 - x^2)],  a^2 = 3/(8 pi G eps), x = R/a.
#[test]
fn tov_reproduces_uniform_density_schwarzschild_interior() {
    for eps_mev in [150.0, 500.0, 1000.0] {
        let eps = eps_mev * MEV_FM3_TO_MSUN_KM3;
        for ratio in [0.01, 0.1, 0.5, 2.0] {
            let pc = ratio * eps;
            let p_min = pc * 1e-12;
            let n = 400;
            let pressures: Vec<f64> = (0..n)
                .map(|i| p_min * (10.0 * pc / p_min).powf(i as f64 / (n - 1) as f64))
                .collect();
            let energies = vec![eps; n];

            let (mass, radius, baryonic_mass, _) =
                integrate_star(pc, p_min, &pressures, &energies, &energies)
                    .expect("uniform-density star must reach its surface");

            let s = (1.0 + ratio) / (1.0 + 3.0 * ratio);
            let radius_exact = ((1.0 - s * s) / (8.0 * PI * G_C2 * eps / 3.0)).sqrt();
            let mass_exact = 4.0 / 3.0 * PI * radius_exact.powi(3) * eps;
            let a = (3.0 / (8.0 * PI * G_C2 * eps)).sqrt();
            let x = radius_exact / a;
            let proper_mass_exact = eps * 2.0 * PI * a.powi(3) * (x.asin() - x * (1.0 - x * x).sqrt());

            let tag = format!("eps={eps_mev} MeV/fm3, Pc/eps={ratio}");
            assert!((radius / radius_exact - 1.0).abs() < 1e-6, "R: {tag}");
            assert!((mass / mass_exact - 1.0).abs() < 1e-6, "M: {tag}");
            assert!(
                (baryonic_mass / proper_mass_exact - 1.0).abs() < 1e-6,
                "M_proper: {tag}"
            );
        }
    }
}

/// Gibbs–Duhem a T = 0 com neutralidade de carga: dP = n_B dmu_n. A derivada
/// é tomada por diferença central na malha em mu_n (passo ~1.4 MeV); o erro
/// de truncamento esperado é O(1e-4) e maior perto dos limiares de novas
/// partículas e níveis de Landau, onde n_B(mu) tem quinas.
fn assert_gibbs_duhem(rows: &[Row], label: &str) {
    let mut errors = Vec::new();
    for w in rows.windows(3) {
        let (a, b, c) = (&w[0], &w[1], &w[2]);
        if b[COL_NB] < 0.5 {
            continue;
        }
        let dp_dmu =
            (matter_pressure(c) - matter_pressure(a)) / ((c[COL_MU_N] - a[COL_MU_N]) * M_NUCLEON);
        let n_b = b[COL_NB] * N0;
        errors.push(((dp_dmu - n_b) / n_b).abs());
    }
    assert!(errors.len() > 100, "{label}: EoS too short ({} dense rows)", errors.len());
    errors.sort_by(f64::total_cmp);
    let p95 = errors[errors.len() * 95 / 100];
    let max = *errors.last().unwrap();
    assert!(p95 < 5e-4, "{label}: 95th percentile |dP/dmu - n_B|/n_B = {p95:e}");
    assert!(max < 5e-3, "{label}: max |dP/dmu - n_B|/n_B = {max:e}");
}

/// Neutralidade de carga nas colunas exportadas. O Newton do solver usa
/// tolerância absoluta 1e-10 em unidades de M_N^3 (~1.1e-8 fm^-3).
fn assert_charge_neutral(rows: &[Row], label: &str) {
    for row in rows {
        let charge = row[COL_NP] + row[COL_NSP]
            - row[COL_NSM]
            - row[COL_NXM]
            - row[COL_NE]
            - row[COL_NMU];
        let tol = 1e-7 + 1e-6 * row[COL_NB] * N0;
        assert!(
            charge.abs() < tol,
            "{label}: net charge {charge:e} fm^-3 at nB/n0 = {}",
            row[COL_NB]
        );
    }
}

#[test]
fn eos_satisfies_gibbs_duhem_and_neutrality_without_field() {
    for (name, model) in [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)] {
        let rows = solve_eos(model, 0.0);
        assert_gibbs_duhem(&rows, name);
        assert_charge_neutral(&rows, name);
    }
}

#[test]
fn eos_satisfies_gibbs_duhem_and_neutrality_with_landau_levels() {
    for (name, model) in [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)] {
        let rows = solve_eos(model, 1e17);
        let label = format!("{name}, B = 1e17 G");
        assert_gibbs_duhem(&rows, &label);
        assert_charge_neutral(&rows, &label);
    }
}

/// Limite B -> 0: a soma sobre níveis de Landau deve recuperar o gás de
/// Fermi isotrópico. Compara P_matéria(mu_n) e eps_matéria(mu_n) nos mesmos
/// pontos da malha. Com B = 1e15 G o desvio medido é ~4e-7.
#[test]
fn landau_quantization_recovers_isotropic_limit() {
    let reference = solve_eos(GM1, 0.0);
    let weak = solve_eos(GM1, 1e15);
    let mut compared = 0;
    for row in &weak {
        if row[COL_NB] < 0.5 {
            continue;
        }
        if let Some(r0) = reference
            .iter()
            .find(|r| (r[COL_MU_N] - row[COL_MU_N]).abs() < 1e-12)
        {
            compared += 1;
            let rel = (matter_pressure(row) / matter_pressure(r0) - 1.0).abs();
            assert!(rel < LANDAU_LIMIT_TOL, "nB/n0 = {}: rel. diff {rel:e}", row[COL_NB]);
            let rel_e = ((row[COL_EPS] - row[COL_EPS_MAG]) / r0[COL_EPS] - 1.0).abs();
            assert!(rel_e < LANDAU_LIMIT_TOL, "nB/n0 = {}: rel. diff eps {rel_e:e}", row[COL_NB]);
        }
    }
    assert!(compared > 100, "only {compared} common grid points");
}

const LANDAU_LIMIT_TOL: f64 = 1e-5;
