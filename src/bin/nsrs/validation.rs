// Relatório de validação em Markdown (`nsrs report validation`).
//
// Roda as verificações numéricas, as comparações com os artigos originais das
// parametrizações, as relações universais, o confronto com observações e os
// testes do setor magnético, e grava docs/VALIDATION_REPORT.md com os valores
// do NSRS, os valores de referência, o critério usado e as citações.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt::Write as _;

use rayon::prelude::*;

use nsrs::constants::{G_C2, M_NUCLEON, MEV_FM3_TO_MSUN_KM3, N0, RESULTS_SIZE};
use nsrs::core::io_utils::derived_diagnostics;
use nsrs::core::magnetic::{B_SURFACE_G, bdd_field_g, magnetic_stress};
use nsrs::core::model::ModelParams;
use nsrs::core::nuclear::{SaturationProperties, saturation_properties};
use nsrs::core::observations::{Status as ObsStatus, assess, interpolate_at_mass, load_constraints, stable_branch};
use nsrs::core::tov_solver::{StarProperties, generate_star_sequence, integrate_star, integrate_star_properties};
use nsrs::{
    DexheimerFit, EngineMode, FSU2, FieldProfile, GM1, GM3, HadronsMatter, MagneticTopology, NlemModel,
    Solver,
};

use crate::cli::{Args, create_dir};
use crate::study::negative_pressure_threshold;

type Row = [f64; RESULTS_SIZE];

const MODELS: [(&str, ModelParams); 3] = [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)];

// ---------------------------------------------------------------------------
// Referências (apenas trabalhos conferidos; as de observações vêm do CSV).
// ---------------------------------------------------------------------------

const REFERENCES: &[(&str, &str)] = &[
    ("schwarzschild1916", "K. Schwarzschild, *Sitzungsber. Preuss. Akad. Wiss. Berlin*, 424 (1916) — solução interior de densidade uniforme."),
    ("hinderer2008", "T. Hinderer, *ApJ* **677**, 1216 (2008). arXiv:0711.2420"),
    ("damour_nagar2009", "T. Damour & A. Nagar, *Phys. Rev. D* **80**, 084035 (2009). arXiv:0906.0096"),
    ("postnikov2010", "S. Postnikov, M. Prakash & J. M. Lattimer, *Phys. Rev. D* **82**, 024016 (2010). arXiv:1004.5098"),
    ("hartle1967", "J. B. Hartle, *ApJ* **150**, 1005 (1967)."),
    ("yagi_yunes2013", "K. Yagi & N. Yunes, *Science* **341**, 365 (2013). doi:10.1126/science.1236462, arXiv:1302.4499"),
    ("lpph1991", "J. M. Lattimer, C. J. Pethick, M. Prakash & P. Haensel, *Phys. Rev. Lett.* **66**, 2701 (1991). doi:10.1103/PhysRevLett.66.2701"),
    ("gm1991", "N. K. Glendenning & S. A. Moszkowski, *Phys. Rev. Lett.* **67**, 2414 (1991). doi:10.1103/PhysRevLett.67.2414"),
    ("cp2014", "W.-C. Chen & J. Piekarewicz, *Phys. Rev. C* **90**, 044305 (2014). doi:10.1103/PhysRevC.90.044305, arXiv:1408.4159"),
    ("nlh2025", "G. Nam, Y. Lim & J. W. Holt, *Universal Relation for the Neutron Star Maximum Mass within Relativistic Mean-Field Theories*, arXiv:2510.15356 (2025), Tabela III."),
    ("bps1971", "G. Baym, C. Pethick & P. Sutherland, *ApJ* **170**, 299 (1971) — crosta BPS."),
    ("bcp1997", "D. Bandyopadhyay, S. Chakrabarty & S. Pal, *Phys. Rev. Lett.* **79**, 2176 (1997). arXiv:astro-ph/9703066"),
    ("dex2017", "V. Dexheimer, B. Franzon, R. O. Gomes, R. L. S. Farias, S. S. Avancini & S. Schramm, *Phys. Lett. B* **773**, 487 (2017). arXiv:1612.05795"),
    ("soleng1995", "H. H. Soleng, *Phys. Rev. D* **52**, 6178 (1995). arXiv:hep-th/9509033"),
    ("ferrer2010", "E. J. Ferrer *et al.*, *Phys. Rev. C* **82**, 065802 (2010)."),
    ("strickland2012", "M. Strickland, V. Dexheimer & D. P. Menezes, *Phys. Rev. D* **86**, 125032 (2012)."),
    ("chatterjee2015", "D. Chatterjee, T. Elghozi, J. Novak & M. Oertel, *MNRAS* **447**, 3785 (2015)."),
];

#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Pass,
    Warn,
    Fail,
    Info,
}

impl Status {
    fn icon(self) -> &'static str {
        match self {
            Status::Pass => "✅",
            Status::Warn => "⚠️",
            Status::Fail => "❌",
            Status::Info => "ℹ️",
        }
    }

    fn check(ok: bool) -> Self {
        if ok { Status::Pass } else { Status::Fail }
    }

    /// Distância em desvios-padrão: <= 1 ✅, <= 2 ⚠️, > 2 ❌.
    fn sigma(distance: f64) -> Self {
        if distance <= 1.0 {
            Status::Pass
        } else if distance <= 2.0 {
            Status::Warn
        } else {
            Status::Fail
        }
    }
}

struct Check {
    quantity: String,
    nsrs: String,
    reference: String,
    criterion: String,
    status: Status,
}

fn check(
    quantity: impl Into<String>,
    nsrs: impl Into<String>,
    reference: impl Into<String>,
    criterion: impl Into<String>,
    status: Status,
) -> Check {
    Check {
        quantity: quantity.into(),
        nsrs: nsrs.into(),
        reference: reference.into(),
        criterion: criterion.into(),
        status,
    }
}

/// Documento em construção: texto, contagem por status e citações numeradas
/// na ordem de primeira aparição.
struct Report {
    body: String,
    cited: Vec<String>,
    extra_refs: HashMap<String, String>,
    counts: Vec<(String, [usize; 4])>,
}

impl Report {
    fn new() -> Self {
        Report { body: String::new(), cited: Vec::new(), extra_refs: HashMap::new(), counts: Vec::new() }
    }

    /// Número da citação `key` ("[n]").
    fn cite(&mut self, key: &str) -> String {
        let known = REFERENCES.iter().any(|(k, _)| *k == key) || self.extra_refs.contains_key(key);
        assert!(known, "referência desconhecida: {key}");
        let n = match self.cited.iter().position(|k| k == key) {
            Some(i) => i + 1,
            None => {
                self.cited.push(key.to_string());
                self.cited.len()
            }
        };
        format!("[{n}]")
    }

    fn section(&mut self, title: &str) {
        let _ = write!(self.body, "\n## {title}\n\n");
        self.counts.push((title.to_string(), [0; 4]));
    }

    fn subsection(&mut self, title: &str) {
        let _ = write!(self.body, "\n### {title}\n\n");
    }

    fn text(&mut self, text: &str) {
        let _ = write!(self.body, "{text}\n\n");
    }

    fn table(&mut self, checks: &[Check]) {
        self.body.push_str("| Grandeza | NSRS | Referência | Critério | |\n|---|---|---|---|:-:|\n");
        let cell = |text: &str| text.replace('|', "\\|");
        for c in checks {
            let _ = writeln!(
                self.body,
                "| {} | {} | {} | {} | {} |",
                cell(&c.quantity),
                cell(&c.nsrs),
                cell(&c.reference),
                cell(&c.criterion),
                c.status.icon()
            );
            if let Some((_, counts)) = self.counts.last_mut() {
                let i = match c.status {
                    Status::Pass => 0,
                    Status::Warn => 1,
                    Status::Fail => 2,
                    Status::Info => 3,
                };
                counts[i] += 1;
            }
        }
        self.body.push('\n');
    }
}

// ---------------------------------------------------------------------------
// Auxiliares numéricos
// ---------------------------------------------------------------------------

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

fn sci(x: f64) -> String {
    format!("{x:.1e}")
}

/// EoS com a coluna 2 convertida para P_par = P_perp + M B (a que obedece
/// Gibbs-Duhem).
fn solve_parallel_pressure(engine: HadronsMatter) -> Vec<Row> {
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let mut rows = solver.solve();
    for (row, diag) in rows.iter_mut().zip(solver.diagnostics()) {
        row[2] += diag.magnetization_b;
    }
    rows
}

/// Pressão da matéria (anisotrópica, Maxwell): P_mag = eps_mag.
fn matter_pressure(row: &Row) -> f64 {
    row[2] - row[19]
}

fn columns(rows: &[Row]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    (
        rows.iter().map(|r| r[1]).collect(),
        rows.iter().map(|r| r[2]).collect(),
        rows.iter().map(|r| r[0]).collect(),
    )
}

/// Sequência estelar com crosta BPS, B = 0, malha em mu até 3 M_N.
struct Sequence {
    rows: Vec<Row>,
    stars: Vec<StarProperties>,
}

impl Sequence {
    fn new(model: ModelParams, hyperons: bool) -> Self {
        let engine = HadronsMatter::new(model, 0.0)
            .with_hyperons(hyperons)
            .with_limits(0.02, 3.0)
            .with_points(2000);
        let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
        let (e, p, n) = columns(&rows);
        let stars = generate_star_sequence(&e, &p, &n, true);
        Sequence { rows, stars }
    }

    fn branch(&self) -> &[StarProperties] {
        stable_branch(&self.stars)
    }

    fn m_max(&self) -> f64 {
        self.branch().last().map_or(f64::NAN, |s| s.mass)
    }

    /// Máximo de fato: há estrelas instáveis depois dele.
    fn has_true_maximum(&self) -> bool {
        self.branch().len() + 10 <= self.stars.len()
    }

    fn at_mass(&self, mass: f64, value: impl Fn(&StarProperties) -> f64) -> Option<f64> {
        interpolate_at_mass(self.branch(), mass, value)
    }
}

// ---------------------------------------------------------------------------
// Seção 1: verificação numérica
// ---------------------------------------------------------------------------

fn numerical(report: &mut Report) {
    report.section("1. Verificação numérica (soluções exatas e identidades)");
    report.text(
        "Estas verificações não dependem de dados empíricos: uma falha indica erro de implementação, \
         não de modelo físico.",
    );

    // 1.1 TOV x Schwarzschild.
    report.subsection("1.1 Integrador TOV");
    let mut worst = [0.0f64; 3];
    for eps_mev in [150.0, 500.0, 1000.0] {
        let eps = eps_mev * MEV_FM3_TO_MSUN_KM3;
        for ratio in [0.01, 0.1, 0.5, 2.0] {
            let pc = ratio * eps;
            let p_min = pc * 1e-12;
            let n = 400;
            let pressures: Vec<f64> =
                (0..n).map(|i| p_min * (10.0 * pc / p_min).powf(i as f64 / (n - 1) as f64)).collect();
            let energies = vec![eps; n];
            let Some((mass, radius, proper, _)) = integrate_star(pc, p_min, &pressures, &energies, &energies)
            else {
                worst = [f64::INFINITY; 3];
                continue;
            };
            let s = (1.0 + ratio) / (1.0 + 3.0 * ratio);
            let r_exact = ((1.0 - s * s) / (8.0 * PI * G_C2 * eps / 3.0)).sqrt();
            let m_exact = 4.0 / 3.0 * PI * r_exact.powi(3) * eps;
            let a = (3.0 / (8.0 * PI * G_C2 * eps)).sqrt();
            let x = r_exact / a;
            let proper_exact = eps * 2.0 * PI * a.powi(3) * (x.asin() - x * (1.0 - x * x).sqrt());
            worst[0] = worst[0].max(rel(radius, r_exact));
            worst[1] = worst[1].max(rel(mass, m_exact));
            worst[2] = worst[2].max(rel(proper, proper_exact));
        }
    }
    let schw = report.cite("schwarzschild1916");
    let tol = 1e-8;
    let checks: Vec<Check> = ["Raio R", "Massa M", "Massa própria"]
        .iter()
        .zip(worst)
        .map(|(q, w)| {
            check(
                format!("{q}, densidade uniforme (12 estrelas, P_c/ε de 0.01 a 2)"),
                format!("erro rel. máx. {}", sci(w)),
                format!("solução exata de Schwarzschild {schw}"),
                format!("< {}", sci(tol)),
                Status::check(w < tol),
            )
        })
        .collect();
    report.table(&checks);

    // 1.2 Limites newtonianos de k2 e I.
    report.subsection("1.2 Número de Love e momento de inércia no limite newtoniano");
    let pi2 = PI * PI;
    let m_r2 = |s: &StarProperties| s.mass * s.radius * s.radius * 1.98847e33 * 1e10 / 1e45;
    let polytrope = {
        let pc = 1e-6;
        let eps_c = (100.0f64 * pc).sqrt();
        let (n, eps_min) = (600, eps_c * 1e-8);
        let eps: Vec<f64> =
            (0..n).map(|i| eps_min * (3.0 * eps_c / eps_min).powf(i as f64 / (n - 1) as f64)).collect();
        let e_t: Vec<f64> = eps.iter().map(|e| e * MEV_FM3_TO_MSUN_KM3).collect();
        let p_t: Vec<f64> = eps.iter().map(|e| e * e / 100.0 * MEV_FM3_TO_MSUN_KM3).collect();
        integrate_star_properties(pc * MEV_FM3_TO_MSUN_KM3, p_t[0], &p_t, &e_t, &e_t)
    };
    let uniform = {
        let eps = 500.0 * MEV_FM3_TO_MSUN_KM3;
        let pc = 1e-5 * eps;
        let (n, p_min) = (400, pc * 1e-12);
        let p: Vec<f64> = (0..n).map(|i| p_min * (10.0 * pc / p_min).powf(i as f64 / (n - 1) as f64)).collect();
        let e = vec![eps; n];
        integrate_star_properties(pc, p_min, &p, &e, &e)
    };
    let (hind, dn, hartle) = (report.cite("hinderer2008"), report.cite("damour_nagar2009"), report.cite("hartle1967"));
    report.text(&format!(
        "k₂ pela equação de Hinderer {hind}, com a correção de descontinuidade de densidade na \
         superfície de Damour & Nagar {dn}; I pela aproximação de rotação lenta de Hartle {hartle}. \
         Em compacidade C → 0 ambos devem recuperar os valores newtonianos analíticos."
    ));
    let mut checks = Vec::new();
    let mut push = |q: &str, value: Option<f64>, exact: f64, reference: String, tol: f64| {
        let r = value.map(|v| rel(v, exact));
        checks.push(check(
            q,
            value.map_or("falhou".into(), |v| format!("{v:.5}")),
            format!("{exact:.5} {reference}"),
            format!("erro rel. < {}", sci(tol)),
            Status::check(r.is_some_and(|r| r < tol)),
        ));
    };
    push(
        "k₂, polítropo n = 1 (C ≈ 2×10⁻⁴)",
        polytrope.map(|s| s.love_k2),
        (15.0 - pi2) / (2.0 * pi2),
        "(15 − π²)/2π² (newtoniano)".into(),
        5e-3,
    );
    push(
        "I/MR², polítropo n = 1",
        polytrope.map(|s| s.moment_of_inertia / m_r2(&s)),
        2.0 / 3.0 * (1.0 - 6.0 / pi2),
        "(2/3)(1 − 6/π²) (newtoniano)".into(),
        5e-3,
    );
    push(
        "k₂, densidade uniforme (C ≈ 2×10⁻⁵)",
        uniform.map(|s| s.love_k2),
        0.75,
        "3/4 (newtoniano, fluido incompressível)".into(),
        1e-3,
    );
    push("I/MR², densidade uniforme", uniform.map(|s| s.moment_of_inertia / m_r2(&s)), 0.4, "2/5 (newtoniano)".into(), 1e-3);
    report.table(&checks);

    // 1.3 Termodinâmica e neutralidade.
    report.subsection("1.3 Consistência termodinâmica e neutralidade de carga");
    report.text(
        "Gibbs–Duhem a T = 0: dP/dμ_n = n_B, com P = P∥ (a pressão termodinâmica, −Ω), derivada por \
         diferença central na malha (passo ≈ 1.4 MeV; erro de truncamento esperado O(10⁻⁴), maior \
         nos limiares de partículas e níveis de Landau). Neutralidade: |Σ q_i n_i| < 10⁻⁷ + 10⁻⁶ n_B fm⁻³.",
    );
    let cases: Vec<(String, ModelParams, f64)> = MODELS
        .iter()
        .flat_map(|&(name, model)| {
            [(format!("{name}, B = 0"), model, 0.0), (format!("{name}, B = 10¹⁷ G"), model, 1e17)]
        })
        .collect();
    let results: Vec<(f64, f64, usize, f64)> = cases
        .par_iter()
        .map(|(_, model, b)| {
            let rows = solve_parallel_pressure(HadronsMatter::new(*model, *b));
            let mut errors: Vec<f64> = rows
                .windows(3)
                .filter(|w| w[1][0] >= 0.5)
                .map(|w| {
                    let dp = (matter_pressure(&w[2]) - matter_pressure(&w[0])) / ((w[2][17] - w[0][17]) * M_NUCLEON);
                    let n_b = w[1][0] * N0;
                    ((dp - n_b) / n_b).abs()
                })
                .collect();
            errors.sort_by(f64::total_cmp);
            let (p95, max) = if errors.is_empty() {
                (f64::INFINITY, f64::INFINITY)
            } else {
                (errors[errors.len() * 95 / 100], *errors.last().unwrap())
            };
            let charge = rows
                .iter()
                .map(|r| (r[6] + r[10] - r[8] - r[11] - r[3] - r[4]).abs() / (1e-7 + 1e-6 * r[0] * N0))
                .fold(0.0, f64::max);
            (p95, max, errors.len(), charge)
        })
        .collect();
    let mut checks = Vec::new();
    for ((label, _, _), (p95, max, n, charge)) in cases.iter().zip(results) {
        checks.push(check(
            format!("Gibbs–Duhem, {label} ({n} pontos com n_B ≥ 0.5 n₀)"),
            format!("p95 {}, máx. {}", sci(p95), sci(max)),
            "identidade exata",
            "p95 < 5e-4, máx. < 5e-3",
            Status::check(n > 100 && p95 < 5e-4 && max < 5e-3),
        ));
        checks.push(check(
            format!("Neutralidade de carga, {label}"),
            format!("máx. |q|/tolerância = {charge:.2}"),
            "0",
            "< 1",
            Status::check(charge < 1.0),
        ));
    }
    report.table(&checks);

    // 1.4 Quantização de Landau e magnetização.
    report.subsection("1.4 Níveis de Landau e magnetização");
    let (reference, weak) = rayon::join(
        || solve_parallel_pressure(HadronsMatter::new(GM1, 0.0)),
        || solve_parallel_pressure(HadronsMatter::new(GM1, 1e15)),
    );
    let (mut landau, mut compared) = (0.0f64, 0);
    for row in weak.iter().filter(|r| r[0] >= 0.5) {
        if let Some(r0) = reference.iter().find(|r| (r[17] - row[17]).abs() < 1e-12) {
            landau = landau
                .max(rel(matter_pressure(row), matter_pressure(r0)))
                .max(rel(row[1] - row[19], r0[1]));
            compared += 1;
        }
    }
    let magnetization = magnetization_check(1.25, 1e-4);
    let oscillation = magnetization_check(1.40, 1e-3);
    let (ferrer, strickland) = (report.cite("ferrer2010"), report.cite("strickland2012"));
    report.table(&[
        check(
            format!("Limite B → 0 da soma sobre níveis de Landau (GM1, B = 10¹⁵ G, {compared} pontos)"),
            format!("dif. rel. máx. em P e ε: {}", sci(landau)),
            "gás de Fermi isotrópico",
            "< 1e-5",
            Status::check(compared > 100 && landau < 1e-5),
        ),
        check(
            "Magnetização 𝓜B = B ∂P∥/∂B exportada (GM1, B = 3×10¹⁷ G, μ_n = 1.25 M_N)",
            magnetization.map_or("falhou".into(), |r| format!("dif. rel. {}", sci(r))),
            format!("diferença finita com passo 10⁻⁴ (interno: 10⁻⁵); P⊥ = P∥ − 𝓜B {ferrer} {strickland}"),
            "< 1e-5",
            Status::check(magnetization.is_some_and(|r| r < 1e-5)),
        ),
        check(
            "Oscilações de de Haas–van Alphen: 𝓜B com passo 10⁻³ vs exportado (μ_n = 1.40 M_N, n_B ≈ 4 n₀)",
            oscillation.map_or("falhou".into(), |r| format!("dif. rel. {r:.1e}")),
            "P(B) oscila quando níveis de Landau cruzam a superfície de Fermi",
            "informativo",
            Status::Info,
        ),
    ]);
}

/// |𝓜B exportado − B ΔP∥/ΔB| / 𝓜B, com diferença central de passo relativo
/// `e` em B, a mu_n fixo (GM1, B = 3e17 G).
fn magnetization_check(mu: f64, e: f64) -> Option<f64> {
    let bg = 3e17;
    let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(GM1, bg).with_limits(0.9, mu).with_points(301)))
        .solve();
    let last = rows.last()?;
    let x = [last[18], last[13] / M_NUCLEON, last[14] / M_NUCLEON, last[15] / M_NUCLEON, 0.0];
    let mut reference = HadronsMatter::new(GM1, bg);
    reference.solve_point(mu, &x)?;
    let mb = reference.magnetization_b;
    let p_par = |scale: f64| {
        let mut engine = HadronsMatter::new(GM1, bg * scale);
        let r = engine.solve_point(mu, &x)?.1;
        Some(r[2] - r[19] + engine.magnetization_b)
    };
    let independent = (p_par(1.0 + e)? - p_par(1.0 - e)?) / (2.0 * e);
    Some(rel(mb, independent))
}

// ---------------------------------------------------------------------------
// Seção 2: parametrizações contra os artigos originais
// ---------------------------------------------------------------------------

fn parametrizations(report: &mut Report, saturation: &[Option<SaturationProperties>], nucleonic: &[Sequence]) {
    report.section("2. Parametrizações contra os artigos originais");
    report.text(
        "Propriedades da matéria nuclear simétrica na saturação, calculadas pelo mesmo motor usado nas \
         estrelas, e massas máximas de estrelas só com núcleons (npeμ, B = 0, crosta BPS). Para GM1 e GM3 \
         os artigos dão valores arredondados, sem incerteza; o critério reflete esse arredondamento. \
         Para FSU2 o critério é a incerteza (1σ) publicada.",
    );

    report.subsection("2.1 Saturação");
    let gm = report.cite("gm1991");
    let nlh = report.cite("nlh2025");
    let cp = report.cite("cp2014");
    let mut checks = Vec::new();
    for (i, (name, _)) in MODELS.iter().enumerate() {
        let Some(s) = saturation[i] else {
            checks.push(check(format!("{name}: saturação"), "não encontrada", "—", "—", Status::Fail));
            continue;
        };
        if *name == "FSU2" {
            // Chen & Piekarewicz (2014): valor ± 1 sigma.
            for (q, v, r, e, d) in [
                ("n₀ [fm⁻³]", s.n0, 0.1505, 0.0007, 4),
                ("E/A [MeV]", s.energy_per_nucleon, -16.28, 0.02, 2),
                ("M*/M", s.effective_mass, 0.593, 0.004, 3),
                ("K [MeV]", s.incompressibility, 238.0, 2.8, 1),
                ("J [MeV]", s.symmetry_energy, 37.62, 1.11, 2),
                ("L [MeV]", s.symmetry_slope, 112.8, 16.1, 1),
            ] {
                let dist = (v - r).abs() / e;
                checks.push(check(
                    format!("FSU2: {q}"),
                    format!("{v:.d$}"),
                    format!("{r} ± {e} {cp}"),
                    format!("|Δ| = {dist:.2}σ"),
                    Status::sigma(dist),
                ));
            }
            continue;
        }
        let (k_gm, m_gm, j_nlh, l_nlh, k_nlh) = if *name == "GM1" {
            (300.0, 0.70, 32.52, 94.04, 300.50)
        } else {
            (240.0, 0.78, 32.51, 89.75, 240.04)
        };
        for (q, v, r, tol, d, src) in [
            ("n₀ [fm⁻³]", s.n0, 0.153, 0.002, 4, &gm),
            ("E/A [MeV]", s.energy_per_nucleon, -16.3, 0.1, 2, &gm),
            ("K [MeV]", s.incompressibility, k_gm, 0.01 * k_gm, 1, &gm),
            ("M*/M", s.effective_mass, m_gm, 0.005, 3, &gm),
            ("J [MeV]", s.symmetry_energy, 32.5, 0.3, 2, &gm),
            ("J [MeV]", s.symmetry_energy, j_nlh, 0.1, 2, &nlh),
            ("L [MeV]", s.symmetry_slope, l_nlh, 0.5, 2, &nlh),
            ("K [MeV]", s.incompressibility, k_nlh, 1.5, 1, &nlh),
        ] {
            let delta = v - r;
            checks.push(check(
                format!("{name}: {q}"),
                format!("{v:.d$}"),
                format!("{r} {src}"),
                format!("|Δ| = {:.d$} ≤ {tol}", delta.abs()),
                Status::check(delta.abs() <= tol),
            ));
        }
    }
    report.table(&checks);

    report.subsection("2.2 Estrelas só com núcleons");
    let bps = report.cite("bps1971");
    let mut checks = Vec::new();
    for (i, (name, _)) in MODELS.iter().enumerate() {
        let seq = &nucleonic[i];
        let m = seq.m_max();
        let truncated = if seq.has_true_maximum() { "" } else { " (máximo no fim da EoS)" };
        let (reference, criterion, status) = match *name {
            "GM1" => (format!("2.363 {nlh}"), "|Δ| < 0.01 M☉ (M_N e crosta)".to_string(), Status::check((m - 2.363).abs() < 0.01)),
            "GM3" => (format!("2.018 {nlh}"), "|Δ| < 0.01 M☉ (M_N e crosta)".to_string(), Status::check((m - 2.018).abs() < 0.01)),
            _ => {
                let d = (m - 2.07).abs() / 0.02;
                (format!("2.07 ± 0.02 {cp}"), format!("|Δ| = {d:.2}σ"), Status::sigma(d))
            }
        };
        let status = if seq.has_true_maximum() { status } else { Status::Fail };
        checks.push(check(format!("{name}: M_max [M☉]"), format!("{m:.3}{truncated}"), reference, criterion, status));
    }
    let r14 = nucleonic[2].at_mass(1.4, |s| s.radius);
    let d = r14.map(|r| (r - 14.42).abs() / 0.26);
    checks.push(check(
        "FSU2: R₁.₄ [km]",
        r14.map_or("—".into(), |r| format!("{r:.2}")),
        format!("14.42 ± 0.26 {cp}"),
        d.map_or("—".into(), |d| format!("|Δ| = {d:.1}σ")),
        d.map_or(Status::Fail, Status::sigma),
    ));
    report.table(&checks);
    report.text(&format!(
        "O raio de FSU2 difere do artigo porque o NSRS junta a tabela BPS {bps} diretamente ao núcleo, \
         enquanto Chen & Piekarewicz {cp} interpolam a crosta interna com um polítropo. A massa máxima, \
         dominada pelo núcleo, não é afetada. Desvio conhecido; uma crosta unificada resolveria."
    ));
}

// ---------------------------------------------------------------------------
// Seção 3: relações universais e propriedades derivadas
// ---------------------------------------------------------------------------

fn derived(report: &mut Report, with_hyperons: &[Sequence], nucleonic: &[Sequence]) {
    report.section("3. Relações universais e propriedades derivadas");
    let yy = report.cite("yagi_yunes2013");
    let (hind, post, hartle) = (report.cite("hinderer2008"), report.cite("postnikov2010"), report.cite("hartle1967"));
    report.text(&format!(
        "Λ e k₂ vêm da equação de Hinderer {hind} {post}; Ī = I/M³ da aproximação de rotação lenta \
         de Hartle {hartle}. A relação universal I-Love {yy} (Tabela I: ln Ī = 1.47 + 0.0817x + \
         0.0149x² + 2.87×10⁻⁴x³ − 3.64×10⁻⁵x⁴, x = ln Λ) tem precisão declarada < 1% e é \
         independente da EoS; ela testa, em conjunto, o integrador de maré e o de inércia."
    ));
    let lpph = report.cite("lpph1991");
    let mut checks = Vec::new();
    for (label, set) in [("com hyperons", with_hyperons), ("só núcleons", nucleonic)] {
        for (i, (name, _)) in MODELS.iter().enumerate() {
            let branch = set[i].branch();
            let deviations: Vec<f64> = branch
                .iter()
                .filter(|s| s.mass >= 1.0)
                .map(|s| {
                    let x = s.tidal_deformability.ln();
                    let fit = (1.47 + 0.0817 * x + 0.0149 * x * x + 2.87e-4 * x.powi(3) - 3.64e-5 * x.powi(4)).exp();
                    rel(s.moment_of_inertia_bar, fit)
                })
                .collect();
            let max = deviations.iter().copied().fold(0.0, f64::max);
            let status = if deviations.len() < 20 {
                Status::Fail
            } else if max < 0.01 {
                Status::Pass
            } else if max < 0.015 {
                Status::Warn
            } else {
                Status::Fail
            };
            checks.push(check(
                format!("I-Love, {name} {label} ({} estrelas, 1 M☉ ≤ M ≤ M_max)", deviations.len()),
                format!("desvio máx. {:.2}%", 100.0 * max),
                format!("ajuste universal {yy}"),
                "< 1% (declarado); ⚠️ até 1.5%",
                status,
            ));
        }
    }
    report.table(&checks);

    report.subsection("3.1 URCA direto, causalidade e propriedades de 1.4 M☉");
    let mut checks = Vec::new();
    for (i, (name, _)) in MODELS.iter().enumerate() {
        let seq = &with_hyperons[i];
        let diag = derived_diagnostics(&seq.rows);
        let onset = (0..seq.rows.len()).find(|&k| seq.rows[k][0] > 1e-3 && diag[k].direct_urca_electron);
        let y_p = onset.map(|k| diag[k].proton_fraction);
        checks.push(check(
            format!("{name}: Y_p no limiar do URCA direto"),
            match (y_p, onset) {
                (Some(y), Some(k)) => format!("{:.2}% (n_B = {:.2} n₀)", 100.0 * y, seq.rows[k][0]),
                _ => "não ocorre".into(),
            },
            format!("11.1% (npe) a 14.8% (npeμ) {lpph}"),
            "dentro da faixa ± 0.5%",
            Status::check(y_p.is_some_and(|y| y > 1.0 / 9.0 - 0.005 && y < 0.148 + 0.005)),
        ));
        let p_c = seq.branch().last().map_or(0.0, |s| s.central_pressure);
        let (cs2, gamma_ok) = seq
            .rows
            .iter()
            .zip(&diag)
            .filter(|(r, _)| r[0] > 0.5 && r[2] <= p_c)
            .fold((0.0f64, true), |(c, g), (_, d)| (c.max(d.sound_speed_squared), g && d.adiabatic_index > 0.0));
        checks.push(check(
            format!("{name}: c_s² máx. até o centro de M_max"),
            format!("{cs2:.3}"),
            "causalidade",
            "c_s² ≤ 1 e Γ > 0",
            Status::check(cs2 <= 1.0 + 1e-9 && gamma_ok),
        ));
        let at14 = |f: fn(&StarProperties) -> f64| seq.at_mass(1.4, f);
        let (r, k2, lam, inertia, z) = (
            at14(|s| s.radius),
            at14(|s| s.love_k2),
            at14(|s| s.tidal_deformability),
            at14(|s| s.moment_of_inertia),
            at14(|s| s.redshift),
        );
        checks.push(check(
            format!("{name}: estrela de 1.4 M☉ (com hyperons)"),
            match (r, k2, lam, inertia, z) {
                (Some(r), Some(k2), Some(l), Some(i), Some(z)) => {
                    format!("R = {r:.2} km, k₂ = {k2:.4}, Λ = {l:.0}, I = {i:.3}×10⁴⁵ g cm², z = {z:.3}")
                }
                _ => "—".into(),
            },
            "ver Seção 4 (observações)",
            "informativo",
            Status::Info,
        ));
    }
    report.table(&checks);
}

// ---------------------------------------------------------------------------
// Seção 4: observações
// ---------------------------------------------------------------------------

fn observations(
    report: &mut Report,
    constraints_path: &str,
    saturation: &[Option<SaturationProperties>],
    with_hyperons: &[Sequence],
    nucleonic: &[Sequence],
) -> Result<(), String> {
    report.section("4. Vínculos observacionais e empíricos");
    let constraints = load_constraints(constraints_path).map_err(|e| format!("{constraints_path}: {e}"))?;
    report.text(&format!(
        "Vínculos de `{constraints_path}`. Distância d em desvios-padrão, com barras assimétricas \
         (intervalos de 90% convertidos para 1σ dividindo por 1.645): d ≤ 1 compatível ✅, 1 < d ≤ 2 \
         tensão ⚠️, d > 2 excluído ❌. Massa máxima: só conta se o modelo não a alcança. Pontos M-R: menor \
         distância ao ramo estável. Estrelas com B = 0 e crosta BPS; \"H\" = com hyperons, \"N\" = só núcleons."
    ));
    let mut header = String::from("| Vínculo | Observado |");
    let mut rule = String::from("|---|---|");
    for (name, _) in MODELS {
        let _ = write!(header, " {name} H | {name} N |");
        rule.push_str(":-:|:-:|");
    }
    let _ = writeln!(report.body, "{header}\n{rule}");
    for c in &constraints {
        let key = format!("csv:{}", c.reference);
        let text = format!(
            "{}. {}{}",
            c.reference,
            if c.doi.is_empty() { String::new() } else { format!("doi:{}", c.doi) },
            if c.arxiv.is_empty() { String::new() } else { format!(", arXiv:{}", c.arxiv) }
        );
        report.extra_refs.insert(key.clone(), text);
        let cite = report.cite(&key);
        let mut line = format!(
            "| {} {cite} | {} (−{}, +{}), {} |",
            c.label, c.value, c.err_minus, c.err_plus, c.credibility
        );
        for (i, _) in MODELS.iter().enumerate() {
            for set in [with_hyperons, nucleonic] {
                let cell = match assess(c, &set[i].stars, saturation[i].as_ref()) {
                    Some(a) => {
                        let status = match a.status {
                            ObsStatus::Compatible => Status::Pass,
                            ObsStatus::Tension => Status::Warn,
                            ObsStatus::Excluded => Status::Fail,
                        };
                        if let Some((_, counts)) = report.counts.last_mut() {
                            counts[match status {
                                Status::Pass => 0,
                                Status::Warn => 1,
                                _ => 2,
                            }] += 1;
                        }
                        format!("{} {:.3} (d = {:.1})", status.icon(), a.model_value, a.distance)
                    }
                    None => "— ¹".into(),
                };
                let _ = write!(line, " {cell} |");
            }
        }
        let _ = writeln!(report.body, "{line}");
    }
    report.text(
        "\n¹ Não comparável: a massa pedida está acima da massa máxima do modelo (a exclusão já aparece \
         no vínculo de massa máxima).\n\nEste confronto avalia os **modelos**, não o código: tensões e \
         exclusões são resultados físicos (p.ex. GM1 e FSU2 rígidos demais para GW170817; hyperons \
         reduzindo M_max abaixo de 2 M☉).",
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Seção 5: campo magnético
// ---------------------------------------------------------------------------

fn magnetic(report: &mut Report) {
    report.section("5. Campo magnético: perfis e acoplamento");
    let (bcp, dex) = (report.cite("bcp1997"), report.cite("dex2017"));
    report.text(&format!(
        "Perfis do campo local: BDD {bcp}, B = B_surf + B₀[1 − exp(−β(n_B/n₀)^γ)] com β = 0.01, \
         γ = 3, B_surf = 10¹⁵ G; e o ajuste polar de Dexheimer et al. {dex}, Eq. (1), \
         B = (a + bμ_B + cμ_B²)μ/B_c, com os coeficientes da Tabela 2."
    ));

    // Dexheimer Eq. (1) e Tabela 2.
    let profile = FieldProfile::Dexheimer2017 { dipole_am2: 3e32, fit: DexheimerFit::BaryonMass22 };
    let b1000 = profile.local_field_g(0.0, 1000.0).unwrap_or(f64::NAN);
    // BDD: autoconsistência (campo usado = campo do perfil na densidade da própria solução).
    let b0 = 1e18;
    let rows = Solver::new(EngineMode::Hadrons(
        HadronsMatter::new(GM1, b0).with_field_profile(FieldProfile::bdd(b0)),
    ))
    .solve();
    let bdd = rows
        .iter()
        .filter(|r| r[0] > 1e-3)
        .map(|r| {
            let used = (r[19] * 8.0 * PI * 1.602176634e33).sqrt();
            rel(used, bdd_field_g(B_SURFACE_G, b0, 0.01, 3.0, r[0]))
        })
        .fold(0.0, f64::max);
    // Dexheimer: Gibbs-Duhem magnetizado dP/dmu = n_B + M dB/dmu.
    let gibbs = dexheimer_gibbs_duhem();
    // Efeito do campo (só na matéria) sobre M_max, dipolo 3e32 A m^2.
    let mass = |engine: HadronsMatter| {
        let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
        let (e, p, n) = columns(&rows);
        let stars = generate_star_sequence(&e, &p, &n, true);
        stable_branch(&stars).last().map_or(f64::NAN, |s| s.mass)
    };
    let (m0, m_dex) = rayon::join(
        || mass(HadronsMatter::new(GM1, 0.0)),
        || mass(HadronsMatter::new(GM1, 0.0).with_field_profile(profile)),
    );
    report.table(&[
        check(
            "Dexheimer, Eq. (1): B(μ_B = 1000 MeV), μ = 3×10³² A m², M_B = 2.2 M☉",
            format!("{b1000:.4e} G"),
            format!("5.777×10¹⁷ G (coeficientes da Tabela 2) {dex}"),
            "dif. rel. < 1e-3",
            Status::check(rel(b1000, 5.777e17) < 1e-3),
        ),
        check(
            "Dexheimer: dP/dμ = n_B + 𝓜 dB/dμ (GM1, n_B = 1, 2, 4, 6 n₀)",
            gibbs.map_or("falhou".into(), |g| format!("dif. rel. máx. {}", sci(g))),
            "identidade termodinâmica com B = B(μ_B)",
            "< 1e-6",
            Status::check(gibbs.is_some_and(|g| g < 1e-6)),
        ),
        check(
            "BDD: campo usado = B(n_B da própria solução) (GM1, B₀ = 10¹⁸ G)",
            format!("dif. rel. máx. {}", sci(bdd)),
            format!("perfil {bcp}"),
            "< 1e-8",
            Status::check(bdd < 1e-8),
        ),
        check(
            "Dexheimer: efeito do campo na matéria sobre M_max (GM1, sem tensões do campo)",
            format!("{m_dex:.4} vs {m0:.4} M☉ ({:+.2}%)", 100.0 * (m_dex / m0 - 1.0)),
            "campo tratado na estrutura por Einstein–Maxwell",
            "informativo",
            Status::Info,
        ),
    ]);
    report.text(
        "Com o perfil de Dexheimer a energia e as tensões do campo ficam fora da EoS da TOV por padrão: o \
         ajuste vem de soluções de Einstein–Maxwell, nas quais o campo já está na estrutura, e é destinado \
         à EoS microscópica.",
    );
}

fn dexheimer_gibbs_duhem() -> Option<f64> {
    let dipole = 3e32;
    let profile = |d: f64| FieldProfile::Dexheimer2017 { dipole_am2: d, fit: DexheimerFit::BaryonMass22 };
    let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(GM1, 0.0).with_field_profile(profile(dipole))))
        .solve();
    let solve_at = |d: f64, mu: f64, x: [f64; 4]| -> Option<Row> {
        let mut engine = HadronsMatter::new(GM1, 0.0).with_field_profile(profile(d));
        let mut row = engine.solve_point(mu, &x)?.1;
        row[2] += engine.magnetization_b;
        Some(row)
    };
    let field = |mu: f64| profile(dipole).local_field_g(0.0, mu * M_NUCLEON);
    let mut worst = 0.0f64;
    for target in [1.0, 2.0, 4.0, 6.0] {
        let row = rows.iter().find(|r| r[0] >= target)?;
        let mu = row[17];
        let x = [row[18], row[13] / M_NUCLEON, row[14] / M_NUCLEON, row[15] / M_NUCLEON];
        let (h, e) = (1e-5, 1e-5);
        let dp = (matter_pressure(&solve_at(dipole, mu + h, x)?) - matter_pressure(&solve_at(dipole, mu - h, x)?))
            / (2.0 * h * M_NUCLEON);
        let b_dp_db = (matter_pressure(&solve_at(dipole * (1.0 + e), mu, x)?)
            - matter_pressure(&solve_at(dipole * (1.0 - e), mu, x)?))
            / (2.0 * e);
        let dlnb = (field(mu + h)?.ln() - field(mu - h)?.ln()) / (2.0 * h * M_NUCLEON);
        worst = worst.max(rel(dp, row[0] * N0 + b_dp_db * dlnb));
    }
    Some(worst)
}

// ---------------------------------------------------------------------------
// Seção 6: NLEM
// ---------------------------------------------------------------------------

fn m_max_of(engine: HadronsMatter) -> Option<f64> {
    let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
    if rows.len() < 5 {
        return None;
    }
    let (e, p, n) = columns(&rows);
    let stars = generate_star_sequence(&e, &p, &n, true);
    stable_branch(&stars).last().map(|s| s.mass)
}

fn nlem(report: &mut Report) {
    report.section("6. Eletrodinâmica não linear (NLEM)");
    let soleng = report.cite("soleng1995");
    report.text(&format!(
        "A NLEM altera só a energia e as tensões do próprio campo: P∥ = −ε_B e P⊥ = HB − ε_B, com \
         H = dε_B/dB {soleng}. No modelo logarítmico ε_B = ξ² ln(1 + x), x = B²/2ξ². Não há, até onde \
         sabemos, resultados publicados de estrelas de nêutrons com a NLEM logarítmica para comparar; as \
         verificações abaixo são de consistência interna, de limites analíticos e de sistemáticos."
    ));

    report.subsection("6.1 Tensor de tensões e acoplamento");
    let xi = 2e17;
    let eps = |b: f64| magnetic_stress(NlemModel::Log(xi), b).energy;
    let identity = [1e16, 1e17, 2e17, 5e17, 3e18]
        .iter()
        .map(|&b| {
            let s = magnetic_stress(NlemModel::Log(xi), b);
            let h = 1e-6 * b;
            let hb = b * (eps(b + h) - eps(b - h)) / (2.0 * h);
            rel(s.p_perpendicular, hb - s.energy)
        })
        .fold(0.0, f64::max);
    let weak = rel(
        magnetic_stress(NlemModel::Log(1e25), 1e15).p_perpendicular,
        magnetic_stress(NlemModel::Maxwell, 1e15).p_perpendicular,
    );
    // A matéria não muda com a NLEM (mesmo B).
    let maxwell_rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(GM1, 3e17))).solve();
    let coupling = [NlemModel::Modmax(1.0), NlemModel::Log(1e16), NlemModel::Log(1e18)]
        .par_iter()
        .map(|&n| {
            let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(GM1, 3e17).with_nlem(n))).solve();
            if rows.len() != maxwell_rows.len() {
                return f64::INFINITY;
            }
            rows.iter()
                .zip(&maxwell_rows)
                .filter(|(a, _)| a[0] > 1e-3)
                .map(|(a, b)| rel(a[1] - a[19], b[1] - b[19]).max(rel(a[0], b[0])))
                .fold(0.0, f64::max)
        })
        .reduce(|| 0.0, f64::max);
    report.table(&[
        check(
            "P⊥ = HB − ε_B, Log(ξ = 2×10¹⁷ G), B de 10¹⁶ a 3×10¹⁸ G",
            format!("dif. rel. máx. {}", sci(identity)),
            format!("H = dε_B/dB numérico {soleng}"),
            "< 1e-7",
            Status::check(identity < 1e-7),
        ),
        check(
            "Log com ξ ≫ B recupera Maxwell (B = 10¹⁵ G, ξ = 10²⁵ G)",
            format!("dif. rel. {}", sci(weak)),
            "limite x → 0",
            "< 1e-14",
            Status::check(weak < 1e-14),
        ),
        check(
            "Matéria idêntica a Maxwell no mesmo B (GM1, 3×10¹⁷ G; ModMax(1), Log(10¹⁶), Log(10¹⁸))",
            format!("dif. rel. máx. em n_B e ε_matéria: {}", sci(coupling)),
            "acoplamento mínimo: Landau usa B",
            "< 1e-12",
            Status::check(coupling < 1e-12),
        ),
    ]);

    report.subsection("6.2 Log(ξ): limites, sinal da pressão do campo e sistemáticos");
    let root = |f: &dyn Fn(f64) -> f64| {
        let (mut lo, mut hi) = (1e-3_f64, 1e3_f64);
        for _ in 0..200 {
            let mid = (lo * hi).sqrt();
            if f(mid) > 0.0 { lo = mid } else { hi = mid }
        }
        (lo * hi).sqrt()
    };
    let analytic_aniso = root(&|x: f64| 2.0 * x - (1.0 + x) * x.ln_1p());
    let analytic_iso = root(&|x: f64| 4.0 * x - 3.0 * (1.0 + x) * x.ln_1p());
    let numeric_aniso = negative_pressure_threshold(MagneticTopology::Anisotropic);
    let numeric_iso = negative_pressure_threshold(MagneticTopology::Isotropic);
    let xi_min = B_SURFACE_G / (2.0 * numeric_aniso).sqrt();

    // Varredura para B0 = 1e18 G, GM1 com hyperons (as duas topologias).
    let b0 = 1e18;
    let xis: Vec<f64> = (0..=12).map(|k| 10f64.powf(16.0 + 0.25 * k as f64)).collect();
    let engine = |topology: MagneticTopology, nlem: Option<NlemModel>| {
        let e = HadronsMatter::new(GM1, b0).with_topology(topology).with_limits(0.02, 3.0).with_points(1500);
        match nlem {
            Some(n) => e.with_nlem(n),
            None => e.with_field_stress(false),
        }
    };
    let mut jobs: Vec<(MagneticTopology, Option<NlemModel>)> = Vec::new();
    for topology in [MagneticTopology::Anisotropic, MagneticTopology::Isotropic] {
        jobs.push((topology, None));
        jobs.push((topology, Some(NlemModel::Maxwell)));
        jobs.push((topology, Some(NlemModel::Log(1e20))));
        jobs.push((topology, Some(NlemModel::Log(1e15))));
        jobs.extend(xis.iter().map(|&xi| (topology, Some(NlemModel::Log(xi)))));
    }
    let masses: Vec<Option<f64>> = jobs.par_iter().map(|&(t, n)| m_max_of(engine(t, n))).collect();
    let (below, above) = rayon::join(
        || m_max_of(engine(MagneticTopology::Anisotropic, Some(NlemModel::Log(0.9 * xi_min)))),
        || m_max_of(engine(MagneticTopology::Anisotropic, Some(NlemModel::Log(1.1 * xi_min)))),
    );
    let get = |t: MagneticTopology, n: Option<NlemModel>| {
        jobs.iter().position(|&j| j == (t, n)).and_then(|i| masses[i])
    };
    let fmt = |m: Option<f64>| m.map_or("—".into(), |m| format!("{m:.4}"));

    let mut checks = vec![
        check(
            "x* (P⊥ do campo < 0 para x > x*), anisotrópica",
            format!("{numeric_aniso:.5}"),
            format!("raiz de 2x = (1+x) ln(1+x): {analytic_aniso:.5}"),
            "dif. rel. < 1e-6",
            Status::check(rel(numeric_aniso, analytic_aniso) < 1e-6),
        ),
        check(
            "x* ((P∥+2P⊥)/3 do campo < 0), isotrópica",
            format!("{numeric_iso:.5}"),
            format!("raiz de 4x = 3(1+x) ln(1+x): {analytic_iso:.5}"),
            "dif. rel. < 1e-6",
            Status::check(rel(numeric_iso, analytic_iso) < 1e-6),
        ),
        check(
            format!("EoS vazia para ξ < B_surf/√(2x*) = {xi_min:.3e} G (anisotrópica)"),
            format!("0.9ξ_min: {}; 1.1ξ_min: {}", if below.is_none() { "vazia" } else { "não vazia" }, fmt(above)),
            "pressão do campo já negativa na superfície (B_surf = 10¹⁵ G)",
            "vazia abaixo, válida acima",
            Status::check(below.is_none() && above.is_some()),
        ),
    ];
    for (topology, label) in [(MagneticTopology::Anisotropic, "anisotrópica"), (MagneticTopology::Isotropic, "isotrópica")] {
        let maxwell = get(topology, Some(NlemModel::Maxwell));
        let no_stress = get(topology, None);
        let large = get(topology, Some(NlemModel::Log(1e20)));
        let small = get(topology, Some(NlemModel::Log(1e15)));
        let diff = |a: Option<f64>, b: Option<f64>| a.zip(b).map(|(a, b)| (a - b).abs());
        let d_large = diff(large, maxwell);
        let d_small = diff(small, no_stress);
        checks.push(check(
            format!("Limite ξ → ∞ (ξ = 10²⁰ G) = Maxwell, {label}"),
            format!("{} vs {} M☉", fmt(large), fmt(maxwell)),
            "limite analítico",
            "|ΔM_max| < 2×10⁻⁴ M☉",
            Status::check(d_large.is_some_and(|d| d < 2e-4)),
        ));
        checks.push(check(
            format!("Limite ξ → 0 (ξ = 10¹⁵ G) = sem tensão do campo, {label}"),
            format!("{} vs {} M☉", fmt(small), fmt(no_stress)),
            "ε_B = ξ² ln(1+x) → 0",
            "|ΔM_max| < 2×10⁻⁴ M☉",
            Status::check(d_small.is_some_and(|d| d < 2e-4)),
        ));
        let (i_min, m_min) = xis
            .iter()
            .enumerate()
            .filter_map(|(i, &xi)| get(topology, Some(NlemModel::Log(xi))).map(|m| (i, m)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap_or((0, f64::NAN));
        let floor = maxwell.zip(no_stress).map(|(a, b)| a.min(b));
        checks.push(check(
            format!("M_max mínimo na varredura ξ = 10¹⁶–10¹⁹ G, {label}"),
            format!("{m_min:.4} M☉ em ξ = {:.2e} G", xis[i_min]),
            format!("mínimo entre Maxwell e sem tensão: {}", fmt(floor)),
            "resultado (janela com pressão do campo < 0)",
            Status::Info,
        ));
    }
    let chat = report.cite("chatterjee2015");
    let aniso = get(MagneticTopology::Anisotropic, Some(NlemModel::Maxwell));
    let iso = get(MagneticTopology::Isotropic, Some(NlemModel::Maxwell));
    checks.push(check(
        "Sistemático de geometria: M_max(P⊥) − M_max(média isotrópica), Maxwell, B₀ = 10¹⁸ G",
        aniso.zip(iso).map_or("—".into(), |(a, i)| format!("{:+.4} M☉", a - i)),
        format!("P⊥ numa TOV esférica é inconsistente {chat}"),
        "incerteza sistemática",
        Status::Warn,
    ));
    report.table(&checks);
    report.text(&format!(
        "GM1 com hyperons, perfil `Constant` (campo da energia pelo perfil BDD com B_surf = 10¹⁵ G). \
         A média isotrópica (P∥ + 2P⊥)/3 é a média angular das tensões de um campo de direção fixa e é \
         consistente com a simetria esférica; usar P⊥ em todas as direções não corresponde a nenhuma \
         geometria esférica consistente {chat}. Os valores absolutos de ΔM na topologia anisotrópica \
         devem ser lidos como cota heurística."
    ));
}

// ---------------------------------------------------------------------------
// Cabeçalho, limitações e montagem
// ---------------------------------------------------------------------------

fn git_revision() -> String {
    let run = |args: &[&str]| {
        std::process::Command::new("git")
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
    };
    match run(&["rev-parse", "--short", "HEAD"]) {
        Some(hash) => {
            let dirty = run(&["status", "--porcelain", "--untracked-files=no"]).is_some_and(|s| !s.is_empty());
            if dirty { format!("{hash} (com alterações locais)") } else { hash }
        }
        None => "desconhecida".into(),
    }
}

/// Data UTC (AAAA-MM-DD) a partir do relógio do sistema.
fn utc_date() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let z = secs.div_euclid(86_400) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

fn limitations(report: &mut Report) {
    report.section("7. Limitações conhecidas");
    let (bps, cp, chat) = (report.cite("bps1971"), report.cite("cp2014"), report.cite("chatterjee2015"));
    let items = [
        format!("**Crosta.** A tabela BPS {bps} é unida diretamente ao núcleo; o raio de FSU2 fica ~0.5 km abaixo do de Chen & Piekarewicz {cp}. Uma crosta unificada resolveria."),
        "**Malha em μ_n.** O padrão (≤ 1.8 M_N) é curto para EoS nucleônicas rígidas; este relatório usa μ_n ≤ 3 M_N nas sequências estelares.".to_string(),
        "**Campo constante forte.** Com o perfil `Constant` e B ≳ 10¹⁸ G há uma transição de primeira ordem na entrada da matéria; a continuação em μ_n pode atravessá-la ou parar, conforme a plataforma. Uma construção de Maxwell tornaria o resultado único.".to_string(),
        format!("**Tensões do campo na TOV.** A TOV é esférica; a topologia anisotrópica (P⊥ em todas as direções) é heurística {chat}. Resultados quantitativos em B ≳ 10¹⁸ G exigem Einstein–Maxwell axissimétrico."),
        "**Maré e inércia com campo.** Λ e I usam as equações de perturbação isotrópicas mesmo quando a EoS inclui tensões anisotrópicas.".to_string(),
        "**Perfil BDD.** β e γ são parâmetros fenomenológicos; o perfil não resolve as equações de Maxwell.".to_string(),
        "**NLEM logarítmica.** Sem resultados publicados de estrelas para comparar; ver Seção 6.".to_string(),
    ];
    for item in items {
        let _ = writeln!(report.body, "- {item}");
    }
    report.body.push('\n');
}

/// `report validation [--out docs/VALIDATION_REPORT.md] [--constraints csv]`
pub fn run(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let out = args.value("out").unwrap_or("docs/VALIDATION_REPORT.md").to_string();
    let constraints = args.value("constraints").unwrap_or("input/observations/constraints.csv").to_string();
    if let Some(n) = args.value("threads") {
        let n: usize = n.parse().map_err(|_| format!("--threads: inteiro inválido '{n}'"))?;
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
    }
    let started = std::time::Instant::now();

    println!("Calculando saturação e sequências estelares...");
    let saturation: Vec<Option<SaturationProperties>> =
        MODELS.par_iter().map(|&(_, m)| saturation_properties(m)).collect();
    let sequences: Vec<Sequence> = MODELS
        .iter()
        .flat_map(|&(_, m)| [(m, true), (m, false)])
        .collect::<Vec<_>>()
        .par_iter()
        .map(|&(m, h)| Sequence::new(m, h))
        .collect();
    let (with_hyperons, nucleonic): (Vec<Sequence>, Vec<Sequence>) = {
        let mut hyp = Vec::new();
        let mut nuc = Vec::new();
        for (i, s) in sequences.into_iter().enumerate() {
            if i % 2 == 0 { hyp.push(s) } else { nuc.push(s) }
        }
        (hyp, nuc)
    };

    let mut report = Report::new();
    println!("1/6 verificação numérica...");
    numerical(&mut report);
    println!("2/6 parametrizações...");
    parametrizations(&mut report, &saturation, &nucleonic);
    println!("3/6 relações universais...");
    derived(&mut report, &with_hyperons, &nucleonic);
    println!("4/6 observações...");
    observations(&mut report, &constraints, &saturation, &with_hyperons, &nucleonic)?;
    println!("5/6 campo magnético...");
    magnetic(&mut report);
    println!("6/6 NLEM...");
    nlem(&mut report);
    limitations(&mut report);

    // Resumo e referências.
    let mut summary = String::from("| Seção | ✅ | ⚠️ | ❌ | ℹ️ |\n|---|:-:|:-:|:-:|:-:|\n");
    let mut total = [0usize; 4];
    for (title, c) in &report.counts {
        if c.iter().sum::<usize>() == 0 {
            continue;
        }
        let _ = writeln!(summary, "| {title} | {} | {} | {} | {} |", c[0], c[1], c[2], c[3]);
        for i in 0..4 {
            total[i] += c[i];
        }
    }
    let _ = writeln!(summary, "| **Total** | **{}** | **{}** | **{}** | **{}** |", total[0], total[1], total[2], total[3]);
    let model_section = report.counts.iter().find(|(t, _)| t.starts_with("4."));
    let code_failures: usize = report
        .counts
        .iter()
        .filter(|(t, _)| !t.starts_with("4."))
        .map(|(_, c)| c[2])
        .sum();
    let _ = write!(
        summary,
        "\nFalhas de implementação (Seções 1–3, 5 e 6): **{code_failures}**. Os ⚠️ e ❌ da Seção 4 \
         ({} e {}) medem os **modelos** contra observações, não o código.\n",
        model_section.map_or(0, |(_, c)| c[1]),
        model_section.map_or(0, |(_, c)| c[2]),
    );

    let mut references = String::from("\n## Referências\n\n");
    for (i, key) in report.cited.iter().enumerate() {
        let text = REFERENCES
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, t)| t.to_string())
            .or_else(|| report.extra_refs.get(key).cloned())
            .unwrap_or_default();
        let _ = writeln!(references, "{}. {text}", i + 1);
    }

    let header = format!(
        "# Relatório de validação do NSRS\n\n\
         > Gerado automaticamente por `cargo run --release --bin nsrs -- report validation` \
         em {} (revisão {}), em {:.0} s. Não edite à mão: rode o comando de novo.\n\n\
         Legenda: ✅ dentro do critério · ⚠️ tensão ou desvio conhecido · ❌ fora do critério · \
         ℹ️ informativo (sem critério).\n\n\
         Critério nas comparações com incerteza publicada: |Δ| ≤ 1σ ✅, ≤ 2σ ⚠️, > 2σ ❌.\n\n\
         ## Resumo\n\n{summary}",
        utc_date(),
        git_revision(),
        started.elapsed().as_secs_f64()
    );
    let document = format!("{header}{}{references}", report.body);

    if let Some(parent) = std::path::Path::new(&out).parent().filter(|p| !p.as_os_str().is_empty()) {
        create_dir(&parent.to_string_lossy())?;
    }
    std::fs::write(&out, document).map_err(|e| format!("{out}: {e}"))?;
    println!(
        "\nRelatório gravado em {out}: {} ✅, {} ⚠️, {} ❌, {} ℹ️ ({:.0} s)",
        total[0],
        total[1],
        total[2],
        total[3],
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
