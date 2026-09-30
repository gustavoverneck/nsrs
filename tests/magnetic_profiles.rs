//! Perfis do campo magnético local (`FieldProfile`) e tensões da NLEM.

use nsrs::constants::{M_NUCLEON, N0, RESULTS_SIZE};
use nsrs::core::magnetic::{B_SURFACE_G, bdd_field_g};
use nsrs::core::model::ModelParams;
use nsrs::{
    DexheimerFit, EngineMode, EosTermination, FSU2, FieldProfile, GM1, GM3, HadronsMatter,
    MagneticTopology, NlemModel, Solver,
};

type Row = [f64; RESULTS_SIZE];

/// Resolve a EoS; a coluna 2 é devolvida como a pressão exportada
/// (perpendicular) e o vetor extra traz M B por linha.
fn solve_with_magnetization(engine: HadronsMatter) -> (Vec<Row>, Vec<f64>, Option<EosTermination>) {
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();
    let mb = solver.diagnostics().iter().map(|d| d.magnetization_b).collect();
    (rows, mb, solver.termination())
}

fn solve(engine: HadronsMatter) -> (Vec<Row>, Option<EosTermination>) {
    let (rows, _, termination) = solve_with_magnetization(engine);
    (rows, termination)
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs())
}

/// O perfil padrão (`Constant`) deve reproduzir o código anterior à
/// introdução dos perfis. Referências: linha 900 (nB/n0 ~ 3.58) das EoS
/// geradas antes da mudança, para Maxwell (anisotrópico e isotrópico). O
/// caso ModMax (1e18 G) saiu: com campo constante de 1e18 G a pressão
/// perpendicular deixa de ser monótona em baixa densidade.
#[test]
fn constant_profile_reproduces_previous_results() {
    // A pressão exportada passou a ser P_perp = P_par - w M B (w = 1 anisotrópico,
    // 2/3 isotrópico); P_par + campo coincide com a referência anterior.
    let cases: [(HadronsMatter, f64, [f64; 3]); 2] = [
        (
            HadronsMatter::new(GM1, 1e17),
            1.0,
            [3.58072868340806139e0, 5.89643597285503688e2, 1.07139674412909841e2],
        ),
        (
            HadronsMatter::new(GM1, 1e17).with_topology(MagneticTopology::Isotropic),
            2.0 / 3.0,
            [3.58072868340806139e0, 5.89643597285503688e2, 1.07115999268299191e2],
        ),
    ];
    for (engine, weight, expected) in cases {
        // A EoS num dado mu não depende da malha: a varredura termina
        // exatamente no mu de referência (mu_n/M_N = 1.3535166...).
        let (rows, mb, _) =
            solve_with_magnetization(engine.with_limits(0.9, 1.35351666666663761).with_points(301));
        let k = rows.len() - 1;
        assert!((rows[k][17] - 1.35351666666663761).abs() < 1e-12);
        let row = rows[k];
        assert!(mb[k] > 0.0);
        let values = [row[0], row[1], row[2] + weight * mb[k]];
        for (col, value) in expected.iter().enumerate() {
            assert!(rel(values[col], *value) < 1e-12, "col {col}: {} vs {value}", values[col]);
        }
    }
}

/// Com o campo local B(n_B) nos níveis de Landau, a matéria de baixa
/// densidade vê ~B_surf e não há a transição de primeira ordem que trunca a
/// EoS no modo `Constant` (GM3 a 1e18 G; GM1 e FSU2 a 2e18 G).
#[test]
fn bdd_profile_covers_the_grid_where_constant_field_truncates() {
    let cases: [(ModelParams, f64); 3] = [(GM3, 1e18), (GM1, 2e18), (FSU2, 2e18)];
    for (model, b0) in cases {
        let (rows, termination) =
            solve(HadronsMatter::new(model, b0).with_field_profile(FieldProfile::bdd(b0)));
        let termination = termination.expect("solver must report termination");
        assert!(!termination.is_anomalous(), "B0 = {b0:e}: {termination:?}");
        let nb_max = rows.iter().map(|r| r[0]).fold(0.0, f64::max);
        assert!(nb_max > 6.0, "B0 = {b0:e}: EoS stops at nB/n0 = {nb_max}");
    }
}

/// Autoconsistência do perfil BDD: o campo usado em cada linha é
/// B(n_B da própria solução). Com Maxwell, eps_mag = B^2/(8 pi).
#[test]
fn bdd_profile_is_self_consistent() {
    let b0 = 1e18;
    let (rows, _) = solve(HadronsMatter::new(GM1, b0).with_field_profile(FieldProfile::bdd(b0)));
    let mut checked = 0;
    for row in rows.iter().filter(|r| r[0] > 1e-3) {
        let b_used = (row[19] * 8.0 * std::f64::consts::PI * 1.602176634e33).sqrt();
        let b_profile = bdd_field_g(B_SURFACE_G, b0, 0.01, 3.0, row[0]);
        assert!(rel(b_used, b_profile) < 1e-8, "nB/n0 = {}", row[0]);
        checked += 1;
    }
    assert!(checked > 500);
}

fn dexheimer(dipole_am2: f64) -> FieldProfile {
    FieldProfile::Dexheimer2017 {
        dipole_am2,
        fit: DexheimerFit::BaryonMass22,
    }
}

fn matter_pressure(row: &Row) -> f64 {
    // Topologia anisotrópica com Maxwell: P_mag = eps_mag (coluna 19).
    row[2] - row[19]
}

/// Com B = B(mu_B), a derivada total da pressão da matéria é
/// dP/dmu = n_B + M dB/dmu, com M = dP/dB a mu fixo. Verifica que cada ponto
/// é resolvido com o campo declarado. Para mu = 3e32 A m^2 o termo de
/// magnetização vale ~1e-3 n_B, muito acima da tolerância.
#[test]
fn dexheimer_profile_satisfies_the_magnetized_gibbs_duhem_relation() {
    let dipole = 3e32;
    let (rows, _) = solve(HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(dipole)));
    // Linha com a coluna 2 = P_par = P_perp + M B.
    let solve_at = |dipole: f64, mu: f64, x: [f64; 4]| {
        let mut engine = HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(dipole));
        let mut row = engine.solve_point(mu, &x).expect("point must converge").1;
        row[2] += engine.magnetization_b;
        row
    };
    let field = |mu: f64| dexheimer(dipole).local_field_g(0.0, mu * M_NUCLEON).unwrap();
    for target in [1.0, 2.0, 4.0, 6.0] {
        let row = rows.iter().find(|r| r[0] >= target).unwrap();
        let mu = row[17];
        let x = [row[18], row[13] / M_NUCLEON, row[14] / M_NUCLEON, row[15] / M_NUCLEON];
        let (h, e) = (1e-5, 1e-5);
        let dp_dmu = (matter_pressure(&solve_at(dipole, mu + h, x))
            - matter_pressure(&solve_at(dipole, mu - h, x)))
            / (2.0 * h * M_NUCLEON);
        // B dP/dB por variação do dipolo (B é linear no dipolo).
        let b_dp_db = (matter_pressure(&solve_at(dipole * (1.0 + e), mu, x))
            - matter_pressure(&solve_at(dipole * (1.0 - e), mu, x)))
            / (2.0 * e);
        let dlnb_dmu = (field(mu + h).ln() - field(mu - h).ln()) / (2.0 * h * M_NUCLEON);
        let n_b = row[0] * N0;
        let magnetization = b_dp_db * dlnb_dmu;
        assert!(magnetization.abs() > 1e-4 * n_b, "magnetization term too small to test");
        assert!(
            rel(dp_dmu, n_b + magnetization) < 1e-6,
            "nB/n0 = {target}: dP/dmu = {dp_dmu}, n_B + M dB/dmu = {}",
            n_b + magnetization
        );
    }
}

/// Dipolo fraco (B <~ 6e15 G): o perfil de Dexheimer recupera a EoS sem campo.
#[test]
fn dexheimer_profile_weak_dipole_recovers_field_free_eos() {
    let (reference, _) = solve(HadronsMatter::new(GM1, 0.0));
    let (mut weak, mb, _) =
        solve_with_magnetization(HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(1e30)));
    for (row, m) in weak.iter_mut().zip(&mb) {
        row[2] += m;
    }
    let mut compared = 0;
    for row in weak.iter().filter(|r| r[0] > 0.5) {
        if let Some(r0) = reference.iter().find(|r| (r[17] - row[17]).abs() < 1e-12) {
            assert!(rel(matter_pressure(row), r0[2]) < 1e-5, "nB/n0 = {}", row[0]);
            compared += 1;
        }
    }
    assert!(compared > 50);
}

/// A NLEM não altera o acoplamento das partículas ao campo: no mesmo B, a
/// matéria (densidades, energia sem o campo) é idêntica à de Maxwell; só a
/// energia e as tensões do campo mudam. Com xi pequeno, a EoS continua
/// cobrindo a malha (antes os níveis de Landau recebiam B(1 + B^2/2xi^2)).
#[test]
fn nlem_changes_only_the_field_stress() {
    // 3e17 G: com campo constante de 1e18 G a pressão perpendicular da matéria
    // (P - M B) deixa de ser monótona em baixa densidade e a EoS é truncada.
    let (maxwell, _) = solve(HadronsMatter::new(GM1, 3e17));
    for nlem in [NlemModel::Modmax(1.0), NlemModel::Log(1e16), NlemModel::Log(1e18)] {
        let (rows, termination) = solve(HadronsMatter::new(GM1, 3e17).with_nlem(nlem));
        assert_eq!(termination, Some(EosTermination::ReachedUpperLimit), "{nlem:?}");
        assert_eq!(rows.len(), maxwell.len(), "{nlem:?}");
        for (row, reference) in rows.iter().zip(&maxwell) {
            assert_eq!(row[0], reference[0], "{nlem:?}: n_B");
            let matter = |r: &Row| r[1] - r[19];
            let (a, b) = (matter(row), matter(reference));
            assert!((a - b).abs() <= 1e-12 * a.abs().max(b.abs()) + 1e-12, "{nlem:?}: eps_matter {a} vs {b}");
        }
    }
}

/// M B exportado nos diagnósticos é B dP_par/dB a mu fixo: confere com uma
/// diferença finita independente entre dois campos constantes.
#[test]
fn magnetization_is_the_field_derivative_of_the_parallel_pressure() {
    let (bg, e) = (3e17, 1e-5);
    let mu = 1.25;
    let rows = Solver::new(EngineMode::Hadrons(
        HadronsMatter::new(GM1, bg).with_limits(0.9, mu).with_points(301),
    ))
    .solve();
    let last = rows.last().expect("continuation");
    assert!((last[17] - mu).abs() < 1e-9);
    let x = [last[18], last[13] / M_NUCLEON, last[14] / M_NUCLEON, last[15] / M_NUCLEON, 0.0];
    let mut reference = HadronsMatter::new(GM1, bg);
    let row = reference.solve_point(mu, &x).unwrap().1;
    let mb = reference.magnetization_b;
    let p_par = |scale: f64| {
        let mut engine = HadronsMatter::new(GM1, bg * scale);
        let r = engine.solve_point(mu, &x).unwrap().1;
        r[2] - r[19] + engine.magnetization_b
    };
    let independent = (p_par(1.0 + e) - p_par(1.0 - e)) / (2.0 * e);
    // O campo da energia (perfil legado) também muda com bg, mas eps_mag é
    // descontado; a diferença restante é a da matéria.
    assert!(rel(mb, independent) < 1e-4, "M B = {mb}, B dP/dB = {independent}");
    assert!(mb.abs() > 1e-6 * (row[2] - row[19] + mb));
}

/// O perfil de Dexheimer é entrada apenas da EoS microscópica: por padrão a
/// energia e as tensões do campo não entram na EoS da TOV (coluna 19 nula e
/// vácuo com P = 0). Com dipolo de 3e32 A m^2 o efeito do campo na matéria
/// muda a massa máxima em menos de 1%.
#[test]
fn dexheimer_profile_excludes_field_stress_by_default() {
    use nsrs::core::tov_solver::generate_mr_curve;
    let max_mass = |rows: &[Row]| {
        let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
        let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
        generate_mr_curve(&e, &p, &n, true).0.into_iter().fold(0.0, f64::max)
    };
    let (field_free, _) = solve(HadronsMatter::new(GM1, 0.0));
    let (rows, termination) =
        solve(HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(3e32)));
    assert_eq!(termination, Some(EosTermination::ReachedUpperLimit));
    assert!(rows.iter().all(|r| r[19] == 0.0));
    assert!(rows.iter().filter(|r| r[0] == 0.0).all(|r| r[2].abs() < 1e-12));
    let (m0, m) = (max_mass(&field_free), max_mass(&rows));
    assert!(rel(m, m0) < 1e-2, "M_max {m} vs {m0}");

    // Com as tensões do campo explicitamente incluídas, B(m_N) ~ 4e17 G
    // cria um envelope sem matéria (vácuo com P > 0).
    let (with_stress, _) = solve(
        HadronsMatter::new(GM1, 0.0)
            .with_field_profile(dexheimer(3e32))
            .with_field_stress(true),
    );
    assert!(with_stress.iter().filter(|r| r[0] == 0.0).all(|r| r[2] > 1.0));
}
