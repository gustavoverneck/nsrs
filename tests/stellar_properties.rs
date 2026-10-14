//! Nível 1/2 — maré (k2, Lambda), momento de inércia e redshift.

use nsrs::constants::MEV_FM3_TO_MSUN_KM3;
use nsrs::core::model::ModelParams;
use nsrs::core::tov_solver::{StarProperties, generate_star_sequence, integrate_star_properties};
use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver};

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// 1 M_sun km^2 em 10^45 g cm^2.
fn m_r2(star: &StarProperties) -> f64 {
    star.mass * star.radius * star.radius * 1.98847e33 * 1e10 / 1e45
}

/// Polítropo n = 1 (P = eps^2/100 em MeV/fm^3) com compacidade ~2e-4.
/// Limites newtonianos analíticos: k2 = (15 - pi^2)/(2 pi^2) e
/// I = (2/3)(1 - 6/pi^2) M R^2.
#[test]
fn n1_polytrope_matches_newtonian_love_number_and_inertia() {
    let pc = 1e-6;
    let eps_c = (100.0 * pc as f64).sqrt();
    let (n, eps_min) = (600, eps_c * 1e-8);
    let eps: Vec<f64> = (0..n)
        .map(|i| eps_min * (3.0 * eps_c / eps_min).powf(i as f64 / (n - 1) as f64))
        .collect();
    let e_t: Vec<f64> = eps.iter().map(|e| e * MEV_FM3_TO_MSUN_KM3).collect();
    let p_t: Vec<f64> = eps.iter().map(|e| e * e / 100.0 * MEV_FM3_TO_MSUN_KM3).collect();
    let star = integrate_star_properties(pc * MEV_FM3_TO_MSUN_KM3, p_t[0], &p_t, &e_t, &e_t)
        .expect("polytrope must reach its surface");
    let pi2 = std::f64::consts::PI.powi(2);
    assert!(star.compactness < 1e-3);
    assert!(rel(star.love_k2, (15.0 - pi2) / (2.0 * pi2)) < 5e-3, "k2 = {}", star.love_k2);
    let inertia = 2.0 / 3.0 * (1.0 - 6.0 / pi2);
    assert!(rel(star.moment_of_inertia / m_r2(&star), inertia) < 5e-3);
}

/// Densidade uniforme, compacidade ~2e-5: k2 = 3/4 (com a correção de
/// descontinuidade na superfície) e I = (2/5) M R^2.
#[test]
fn uniform_density_matches_newtonian_love_number_and_inertia() {
    let eps = 500.0 * MEV_FM3_TO_MSUN_KM3;
    let pc = 1e-5 * eps;
    let (n, p_min) = (400, pc * 1e-12);
    let p: Vec<f64> = (0..n)
        .map(|i| p_min * (10.0 * pc / p_min).powf(i as f64 / (n - 1) as f64))
        .collect();
    let e = vec![eps; n];
    let star = integrate_star_properties(pc, p_min, &p, &e, &e).unwrap();
    assert!(rel(star.love_k2, 0.75) < 1e-3, "k2 = {}", star.love_k2);
    assert!(rel(star.moment_of_inertia / m_r2(&star), 0.4) < 1e-3);
}

fn sequence(model: ModelParams) -> Vec<StarProperties> {
    let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(model, 0.0))).solve();
    let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
    let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    generate_star_sequence(&e, &p, &n, true)
}

/// Relação universal I-Love de Yagi & Yunes, Science 341, 365 (2013),
/// Tabela I: ln I = 1.47 + 0.0817 x + 0.0149 x^2 + 2.87e-4 x^3 - 3.64e-5 x^4,
/// x = ln Lambda, precisão declarada < 1%. Vale para estrelas de 1 M_sun até
/// a massa máxima dos três modelos.
#[test]
fn stars_follow_the_universal_i_love_relation() {
    for model in [GM1, GM3, FSU2] {
        let stars = sequence(model);
        let i_max = (0..stars.len())
            .max_by(|&a, &b| stars[a].mass.total_cmp(&stars[b].mass))
            .unwrap();
        let mut checked = 0;
        for star in stars[..=i_max].iter().filter(|s| s.mass >= 1.0) {
            let x = star.tidal_deformability.ln();
            let fit = (1.47 + 0.0817 * x + 0.0149 * x * x + 2.87e-4 * x.powi(3)
                - 3.64e-5 * x.powi(4))
            .exp();
            assert!(
                rel(star.moment_of_inertia_bar, fit) < 1.5e-2,
                "M = {}: I = {} vs fit {fit}",
                star.mass,
                star.moment_of_inertia_bar
            );
            let z = 1.0 / (1.0 - 2.0 * star.compactness).sqrt() - 1.0;
            assert!(rel(star.redshift, z) < 1e-12);
            checked += 1;
        }
        assert!(checked > 20, "{checked} stars");
    }
}

/// URCA direto nucleônico: em matéria npe o limiar é Y_p = 1/9; com múons
/// sobe até ~14.8% (Lattimer, Pethick, Prakash & Haensel, PRL 66, 2701
/// (1991)). A EoS exportada é causal (c_s^2 <= 1) e tem Gamma > 0.
#[test]
fn direct_urca_threshold_and_causality() {
    use nsrs::core::io_utils::derived_diagnostics;
    for model in [GM1, GM3, FSU2] {
        let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(model, 0.0))).solve();
        let derived = derived_diagnostics(&rows);
        let onset = (0..rows.len())
            .find(|&i| rows[i][0] > 1e-3 && derived[i].direct_urca_electron)
            .expect("direct Urca must open in these models");
        let y_p = derived[onset].proton_fraction;
        assert!(y_p > 1.0 / 9.0 - 0.005 && y_p < 0.148 + 0.005, "Y_p at DU onset = {y_p}");
        for (row, d) in rows.iter().zip(&derived).filter(|(r, _)| r[0] > 0.5) {
            assert!(d.sound_speed_squared <= 1.0 + 1e-9, "c_s^2 = {} at nB/n0 = {}", d.sound_speed_squared, row[0]);
            assert!(d.adiabatic_index > 0.0);
        }
    }
}

/// Passo adaptativo em mu_n: com tol = 1e-3 as estrelas concordam com uma
/// malha adaptativa 4x mais fina (M_max a 1e-4, R_1.4 e Lambda_1.4 a 2e-3).
/// A malha uniforme de 1500 pontos é grossa no limiar e erra R_1.4 em ~0.15%
/// já com B = 0.
#[test]
fn adaptive_grid_converges_star_properties() {
    use nsrs::core::observations::{interpolate_at_mass, stable_branch};
    let summary = |tol: f64| {
        let engine = HadronsMatter::new(GM1, 0.0).with_limits(0.02, 3.0).with_points(1500).with_step_tolerance(tol);
        let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
        let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
        let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
        let stars = generate_star_sequence(&e, &p, &n, true);
        let branch = stable_branch(&stars);
        (
            branch.last().unwrap().mass,
            interpolate_at_mass(branch, 1.4, |s| s.radius).unwrap(),
            interpolate_at_mass(branch, 1.4, |s| s.tidal_deformability).unwrap(),
        )
    };
    let (coarse, fine) = (summary(1e-3), summary(2.5e-4));
    assert!(rel(coarse.0, fine.0) < 1e-4, "M_max {coarse:?} vs {fine:?}");
    assert!(rel(coarse.1, fine.1) < 2e-3, "R_1.4 {coarse:?} vs {fine:?}");
    assert!(rel(coarse.2, fine.2) < 2e-3, "Lambda_1.4 {coarse:?} vs {fine:?}");
}
