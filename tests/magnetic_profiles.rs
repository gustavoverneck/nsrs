//! Perfis do campo magnético local (`FieldProfile`) e tensões da NLEM.

use nsrs::constants::{M_NUCLEON, N0, RESULTS_SIZE};
use nsrs::core::magnetic::{B_SURFACE_G, bdd_field_g};
use nsrs::core::model::ModelParams;
use nsrs::{
    DexheimerFit, EngineMode, EosTermination, FSU2, FieldProfile, GM1, GM3, HadronsMatter,
    MagneticTopology, NlemModel, Solver,
};

type Row = [f64; RESULTS_SIZE];

fn solve(engine: HadronsMatter) -> (Vec<Row>, Option<EosTermination>) {
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();
    (rows, solver.termination())
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs())
}

/// O perfil padrão (`Constant`) deve reproduzir o código anterior à
/// introdução dos perfis. Referências: linha 900 (nB/n0 ~ 3.58) das EoS
/// geradas antes da mudança, para Maxwell (anisotrópico e isotrópico). A
/// referência ModMax foi regenerada quando os níveis de Landau passaram a
/// usar B em vez de e^{-gamma} B (acoplamento mínimo).
#[test]
fn constant_profile_reproduces_previous_results() {
    let cases: [(HadronsMatter, [f64; 3]); 3] = [
        (
            HadronsMatter::new(GM1, 1e17),
            [3.58072868340806139e0, 5.89643597285503688e2, 1.07139674412909841e2],
        ),
        (
            HadronsMatter::new(GM1, 1e17).with_topology(MagneticTopology::Isotropic),
            [3.58072868340806139e0, 5.89643597285503688e2, 1.07115999268299191e2],
        ),
        (
            HadronsMatter::new(GM1, 1e18).with_nlem(NlemModel::Modmax(1.0)),
            [3.53452528898327190e0, 5.81608724074900010e2, 1.05372571197183376e2],
        ),
    ];
    for (engine, expected) in cases {
        let (rows, _) = solve(engine);
        for (col, value) in expected.iter().enumerate() {
            assert!(rel(rows[899][col], *value) < 1e-12, "col {col}: {} vs {value}", rows[899][col]);
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
    let solve_at = |dipole: f64, mu: f64, x: [f64; 4]| {
        let mut engine = HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(dipole));
        engine.solve_point(mu, &x).expect("point must converge").1
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
    let (weak, _) = solve(HadronsMatter::new(GM1, 0.0).with_field_profile(dexheimer(1e30)));
    let mut compared = 0;
    for row in weak.iter().filter(|r| r[0] > 0.5) {
        if let Some(r0) = reference.iter().find(|r| (r[17] - row[17]).abs() < 1e-12) {
            assert!(rel(matter_pressure(row), r0[2]) < 1e-5, "nB/n0 = {}", row[0]);
            compared += 1;
        }
    }
    assert!(compared > 100);
}

/// A NLEM não altera o acoplamento das partículas ao campo: no mesmo B, a
/// matéria (densidades, energia sem o campo) é idêntica à de Maxwell; só a
/// energia e as tensões do campo mudam. Com xi pequeno, a EoS continua
/// cobrindo a malha (antes os níveis de Landau recebiam B(1 + B^2/2xi^2)).
#[test]
fn nlem_changes_only_the_field_stress() {
    let (maxwell, _) = solve(HadronsMatter::new(GM1, 1e18));
    for nlem in [NlemModel::Modmax(1.0), NlemModel::Log(1e16), NlemModel::Log(1e18)] {
        let (rows, termination) = solve(HadronsMatter::new(GM1, 1e18).with_nlem(nlem));
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
