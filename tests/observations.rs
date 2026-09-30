//! Nível 3 — arquivo de vínculos e critérios de avaliação.

use nsrs::core::observations::{Constraint, ConstraintKind, Status, assess, load_constraints};
use nsrs::core::tov_solver::StarProperties;

fn star(mass: f64, radius: f64, lambda: f64) -> StarProperties {
    StarProperties {
        mass,
        radius,
        baryonic_mass: mass,
        central_pressure: 0.0,
        compactness: 1.4766 * mass / radius,
        redshift: 0.0,
        love_k2: 0.0,
        tidal_deformability: lambda,
        moment_of_inertia: 0.0,
        moment_of_inertia_bar: 0.0,
    }
}

fn constraint(kind: ConstraintKind, value: f64, em: f64, ep: f64, credibility: &str) -> Constraint {
    Constraint {
        kind,
        label: "test".into(),
        value,
        err_minus: em,
        err_plus: ep,
        at_mass: Some(1.4),
        mass: Some((1.4, 0.1, 0.1)),
        credibility: credibility.into(),
        reference: String::new(),
        doi: String::new(),
        arxiv: String::new(),
    }
}

/// Todas as linhas têm referência, DOI e arXiv, erros positivos e tipos
/// conhecidos.
#[test]
fn constraints_file_is_complete() {
    let constraints = load_constraints("input/observations/constraints.csv").unwrap();
    assert!(constraints.len() >= 13);
    for c in &constraints {
        assert!(c.err_minus > 0.0 && c.err_plus > 0.0, "{}", c.label);
        assert!(!c.reference.is_empty() && c.doi.starts_with("10.") && !c.arxiv.is_empty(), "{}", c.label);
        if c.kind == ConstraintKind::MassRadius {
            assert!(c.mass.is_some(), "{}", c.label);
        }
        if matches!(c.kind, ConstraintKind::RadiusAtMass | ConstraintKind::TidalAtMass) {
            assert!(c.at_mass.is_some(), "{}", c.label);
        }
    }
}

#[test]
fn scoring_uses_asymmetric_errors_and_credibility() {
    // Ramo estável até 2.0 Msun e um ponto instável depois.
    let stars: Vec<StarProperties> = (0..=20)
        .map(|k| {
            let m = 0.1 * k as f64;
            star(m, 13.0 - m, 2000.0 - 1000.0 * m)
        })
        .chain([star(1.9, 10.0, 10.0)])
        .collect();

    // Massa máxima 2.0: 2.01 +- 0.04 -> d = 0.25; 2.2 +- 0.05 -> d = 4.
    let m = constraint(ConstraintKind::MaxMass, 2.01, 0.04, 0.04, "68%");
    assert!((assess(&m, &stars, None).unwrap().distance - 0.25).abs() < 1e-12);
    let m = constraint(ConstraintKind::MaxMass, 2.2, 0.05, 0.05, "68%");
    assert_eq!(assess(&m, &stars, None).unwrap().status, Status::Excluded);

    // R(1.4) = 11.6: 12.0 (-0.2, +1.0) -> d = 2 (lado de baixo).
    let r = constraint(ConstraintKind::RadiusAtMass, 12.0, 0.2, 1.0, "68%");
    let a = assess(&r, &stars, None).unwrap();
    assert!((a.model_value - 11.6).abs() < 1e-12 && (a.distance - 2.0).abs() < 1e-9);

    // Lambda(1.4) = 600 com 190 (-120, +390) a 90%: d = 410/(390/1.645).
    let t = constraint(ConstraintKind::TidalAtMass, 190.0, 120.0, 390.0, "90%");
    let a = assess(&t, &stars, None).unwrap();
    assert!((a.distance - 410.0 * 1.645 / 390.0).abs() < 1e-9);
    assert_eq!(a.status, Status::Tension);

    // Ponto M-R sobre a curva: d = 0.
    let mut p = constraint(ConstraintKind::MassRadius, 11.6, 1.0, 1.0, "68%");
    p.mass = Some((1.4, 0.1, 0.1));
    assert!(assess(&p, &stars, None).unwrap().distance < 1e-12);
}
