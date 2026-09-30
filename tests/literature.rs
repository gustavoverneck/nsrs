//! Nível 2 — reprodução de resultados publicados para as parametrizações
//! usadas no NSRS. Cada valor de referência foi conferido no artigo citado.
//! Estrelas só com núcleons (`with_hyperons(false)`), como nos artigos.

use nsrs::core::model::ModelParams;
use nsrs::core::nuclear::saturation_properties;
use nsrs::core::tov_solver::{StarProperties, generate_star_sequence};
use nsrs::{EngineMode, EosTermination, FSU2, GM1, GM3, HadronsMatter, Solver};

/// Sequência estelar só com núcleons, com crosta BPS. A malha em mu vai até
/// 3 M_N: com o limite padrão (1.8 M_N) a EoS rígida do GM1 termina antes do
/// centro da estrela de massa máxima.
fn nucleonic_stars(model: ModelParams) -> Vec<StarProperties> {
    let engine = HadronsMatter::new(model, 0.0)
        .with_hyperons(false)
        .with_limits(0.02, 3.0)
        .with_points(2000);
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();
    assert!(!solver.termination().unwrap().is_anomalous());
    assert!(rows.iter().all(|r| r[7..13].iter().all(|&n| n == 0.0)));
    let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
    let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    generate_star_sequence(&e, &p, &n, true)
}

/// Massa máxima, exigindo que seja um máximo de fato (há estrelas depois).
fn maximum_mass(stars: &[StarProperties]) -> f64 {
    let i = (0..stars.len())
        .max_by(|&a, &b| stars[a].mass.total_cmp(&stars[b].mass))
        .unwrap();
    assert!(
        i + 10 < stars.len(),
        "maximum at the end of the sequence: EoS too short"
    );
    stars[i].mass
}

/// Nam & Lim, arXiv:2510.15356 (2025), Tabela III (RMF só com núcleons,
/// crosta BPS): GM1 M_max = 2.363, GM3 M_max = 2.018 M_sun.
/// Diferenças esperadas de ~0.2%: M_N e tratamento da crosta.
#[test]
fn gm_nucleonic_maximum_masses_match_nam_lim() {
    for (model, expected) in [(GM1, 2.363), (GM3, 2.018)] {
        let m_max = maximum_mass(&nucleonic_stars(model));
        assert!(
            (m_max - expected).abs() < 0.01,
            "M_max = {m_max} vs {expected}"
        );
    }
}

/// Chen & Piekarewicz, PRC 90, 044305 (2014): FSU2, M_max = 2.07 +- 0.02 M_sun.
/// (O mesmo artigo dá R_1.4 = 14.42 +- 0.26 km com uma interpolação
/// politrópica entre a crosta externa BPS e o núcleo; o NSRS junta a tabela
/// BPS diretamente ao núcleo e obtém ~13.95 km. Não testado.)
#[test]
fn fsu2_nucleonic_maximum_mass_matches_chen_piekarewicz() {
    let m_max = maximum_mass(&nucleonic_stars(FSU2));
    assert!((m_max - 2.07).abs() < 0.02, "M_max = {m_max}");
}

/// Nam & Lim (2025), Tabela III: GM1 (J, L, K) = (32.52, 94.04, 300.50) MeV;
/// GM3 (32.51, 89.75, 240.04) MeV.
#[test]
fn gm_symmetry_slope_and_incompressibility_match_nam_lim() {
    for (model, j, l, k) in [(GM1, 32.52, 94.04, 300.50), (GM3, 32.51, 89.75, 240.04)] {
        let s = saturation_properties(model).unwrap();
        assert!(
            (s.symmetry_energy - j).abs() < 0.1,
            "J = {}",
            s.symmetry_energy
        );
        assert!(
            (s.symmetry_slope - l).abs() < 0.5,
            "L = {}",
            s.symmetry_slope
        );
        assert!(
            (s.incompressibility - k).abs() < 1.5,
            "K = {}",
            s.incompressibility
        );
    }
}

/// Sem hyperons a EoS é mais rígida: M_max cresce (GM1: 1.99 -> 2.36).
#[test]
fn hyperons_soften_the_equation_of_state() {
    let with_hyperons = {
        let engine = HadronsMatter::new(GM1, 0.0)
            .with_limits(0.02, 3.0)
            .with_points(2000);
        let mut solver = Solver::new(EngineMode::Hadrons(engine));
        let rows = solver.solve();
        assert!(matches!(
            solver.termination(),
            Some(
                EosTermination::NonPositiveEffectiveMass { .. } | EosTermination::ReachedUpperLimit
            )
        ));
        let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
        let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
        maximum_mass(&generate_star_sequence(&e, &p, &n, true))
    };
    let nucleonic = maximum_mass(&nucleonic_stars(GM1));
    assert!(
        nucleonic > with_hyperons + 0.3,
        "{nucleonic} vs {with_hyperons}"
    );
}
