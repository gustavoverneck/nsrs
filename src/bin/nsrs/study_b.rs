// `study b`: varredura completa em B, de B = 0 até o código deixar de produzir
// estrelas, com propriedades estelares, cobertura da EoS e estabilidade local.
//
// Para cada B (e cada perfil de campo):
// - estrelas com crosta BPS nas duas pressões da TOV: P_perp (topologia
//   anisotrópica) e a média isotrópica (P_par + 2 P_perp)/3;
// - cobertura da EoS (n_B máximo, motivo do término);
// - estabilidade mecânica: P_perp deve crescer com n_B;
// - estabilidade magnética (só campo constante, física local): convexidade da
//   energia total em B a mu_n fixo,
//       s = 4 pi d^2 eps_campo/dB^2 - 4 pi d^2 P_m/dB^2 = 1 - 4 pi d(MB/B)/dB
//   (Maxwell, unidades gaussianas). s < 0: instável à formação de domínios
//   magnéticos. É o complemento de Schur da hessiana de eps(n_B, B), logo o
//   critério a mu fixo já inclui a estabilidade em n_B. A derivada segunda é
//   tomada por diferença central com passo relativo `delta` em B e repetida
//   com delta/2; o ponto só conta como instável se s < 0 nos dois (o passo
//   precisa resolver as oscilações de de Haas-van Alphen e o sinal não pode
//   depender dele).

use std::fs;
use std::io::Write;

use rayon::prelude::*;

use nsrs::constants::{M_NUCLEON, RESULTS_SIZE};
use nsrs::core::model::ModelParams;
use nsrs::core::observations::{interpolate_at_mass, stable_branch};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, FieldProfile, HadronsMatter, Solver};

use crate::cli::{Args, create_dir};

type Row = [f64; RESULTS_SIZE];

/// MeV/fm^3 -> erg/cm^3 (= G^2).
const ERG_PER_MEV_FM3: f64 = 1.602176634e33;

#[derive(Clone, Copy, PartialEq)]
enum Profile {
    /// Níveis de Landau com o B central em todas as densidades; energia do
    /// campo pelo perfil BDD (comportamento padrão do código).
    Constant,
    /// B(n_B) de Bandyopadhyay, Chakrabarty & Pal (1997) nos níveis de Landau
    /// e na energia do campo, com B0 = B.
    Bdd,
}

impl Profile {
    fn label(self) -> &'static str {
        match self {
            Profile::Constant => "constante",
            Profile::Bdd => "bdd",
        }
    }
}

struct Config {
    points: usize,
    mu_max: f64,
    hyperons: bool,
}

/// EoS resolvida, com 𝓜B e a pressão sem magnetização por linha.
struct Eos {
    rows: Vec<Row>,
    magnetization_b: Vec<f64>,
    stability_pressure: Vec<f64>,
    termination: String,
    anomalous: bool,
}

fn solve(model: ModelParams, b: f64, profile: Profile, config: &Config) -> Eos {
    let mut engine = HadronsMatter::new(model, b)
        .with_hyperons(config.hyperons)
        .with_limits(0.02, config.mu_max)
        .with_points(config.points);
    if profile == Profile::Bdd && b > 0.0 {
        engine = engine.with_field_profile(FieldProfile::bdd(b));
    }
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();
    let termination = solver.termination();
    Eos {
        magnetization_b: solver.diagnostics().iter().map(|d| d.magnetization_b).collect(),
        stability_pressure: solver.diagnostics().iter().map(|d| d.stability_pressure).collect(),
        termination: termination.map_or("-".into(), |t| format!("{t:?}").replace(',', ";")),
        anomalous: termination.is_some_and(|t| t.is_anomalous()),
        rows,
    }
}

struct Stars {
    m_max: f64,
    r_max: f64,
    nc: f64,
    r14: f64,
    lambda14: f64,
    true_maximum: bool,
}

fn stars(rows: &[Row], pressure: &[f64]) -> Stars {
    let empty = Stars {
        m_max: f64::NAN,
        r_max: f64::NAN,
        nc: f64::NAN,
        r14: f64::NAN,
        lambda14: f64::NAN,
        true_maximum: false,
    };
    if rows.len() < 5 {
        return empty;
    }
    let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let sequence = generate_star_sequence(&eps, pressure, &n, true);
    let branch = stable_branch(&sequence);
    let Some(max) = branch.last() else {
        return empty;
    };
    // Densidade central: menor n_B do núcleo com P >= P_c.
    let nc = rows
        .iter()
        .zip(pressure)
        .filter(|(r, p)| r[0] > 0.0 && **p >= max.central_pressure)
        .map(|(r, _)| r[0])
        .fold(f64::NAN, f64::min);
    Stars {
        m_max: max.mass,
        r_max: max.radius,
        nc,
        r14: interpolate_at_mass(branch, 1.4, |s| s.radius).unwrap_or(f64::NAN),
        lambda14: interpolate_at_mass(branch, 1.4, |s| s.tidal_deformability).unwrap_or(f64::NAN),
        true_maximum: branch.len() + 10 <= sequence.len(),
    }
}

/// Pressão isotrópica a partir da EoS anisotrópica (Maxwell):
/// P_iso = P_par,m - (2/3) 𝓜B + eps_B/3, com P_par,m = P_estab - eps_B.
fn isotropic_pressure(eos: &Eos) -> Vec<f64> {
    eos.rows
        .iter()
        .zip(&eos.magnetization_b)
        .zip(&eos.stability_pressure)
        .map(|((r, mb), stab)| (stab - r[19]) - 2.0 / 3.0 * mb + r[19] / 3.0)
        .collect()
}

/// Linha do stars.csv para um (modelo, perfil, B).
fn star_lines(model_name: &str, profile: Profile, b: f64, eos: &Eos, hyperons: bool) -> Vec<String> {
    let n_max = eos.rows.iter().map(|r| r[0]).fold(0.0, f64::max);
    let mut core: Vec<(f64, f64)> = eos.rows.iter().filter(|r| r[0] > 0.0).map(|r| (r[0], r[2])).collect();
    core.sort_by(|a, b| a.0.total_cmp(&b.0));
    let nonmonotonic = core.windows(2).filter(|w| w[1].1 <= w[0].1).count() as f64 / core.len().max(1) as f64;
    let negative = core.iter().filter(|p| p.1 < 0.0).count();
    let perpendicular: Vec<f64> = eos.rows.iter().map(|r| r[2]).collect();
    let isotropic = isotropic_pressure(eos);
    // Sem núcleo (EoS abaixo de n0) não há estrela de nêutrons: a TOV
    // produziria objetos de crosta e energia de campo.
    let has_core = n_max >= 1.0;
    [("perp", perpendicular), ("iso", isotropic)]
        .into_iter()
        .map(|(topology, pressure)| {
            let s = if has_core { stars(&eos.rows, &pressure) } else { stars(&[], &pressure) };
            format!(
                "{model_name},{},{b:.4e},{hyperons},{},{n_max:.4},{},{},{topology},{:.5},{:.4},{:.4},{:.4},{:.2},{},{nonmonotonic:.4},{negative}",
                profile.label(),
                eos.rows.len(),
                eos.termination,
                eos.anomalous,
                s.m_max,
                s.r_max,
                s.nc,
                s.r14,
                s.lambda14,
                s.true_maximum,
            )
        })
        .collect()
}

/// Ponto da estabilidade local: n_B/n0, mu_n, s com passo delta e delta/2, e
/// dP_perp/dn_B.
struct StabilityPoint {
    n: f64,
    mu: f64,
    s: f64,
    s_half: f64,
    dp_dn: f64,
}

impl StabilityPoint {
    /// Instabilidade magnética robusta: s < 0 com os dois passos.
    fn magnetically_unstable(&self) -> bool {
        self.s < 0.0 && self.s_half < 0.0
    }
}

fn stability_points(model: ModelParams, b: f64, delta: f64, config: &Config) -> Vec<StabilityPoint> {
    let center = solve(model, b, Profile::Constant, config);
    // Cada ponto é resolvido de novo no mesmo mu_n com B(1 ± delta) e
    // B(1 ± delta/2), a partir da solução central. (A malha do solver é
    // adaptativa: execuções separadas não caem nos mesmos mu_n.)
    let engine = |factor: f64| {
        HadronsMatter::new(model, b * factor)
            .with_hyperons(config.hyperons)
            .with_limits(0.02, config.mu_max)
            .with_points(config.points)
    };
    let factors = [1.0 + delta, 1.0 - delta, 1.0 + 0.5 * delta, 1.0 - 0.5 * delta];
    let rows = &center.rows;
    // Os pontos de densidade são independentes: também em paralelo (o rayon
    // divide o trabalho entre este nível e o de B), o que evita que um B
    // lento segure uma única thread no fim da fase.
    let indices: Vec<usize> = (1..rows.len().saturating_sub(1)).filter(|&i| rows[i][0] >= 0.05).collect();
    let points = indices.par_iter().with_max_len(1).filter_map(|&i| {
        let row = &rows[i];
        let guess = [row[18], row[13] / M_NUCLEON, row[14] / M_NUCLEON, row[15] / M_NUCLEON, 0.0];
        // dP_m/dB = 𝓜B/B (erg/cm^3/G) em cada campo.
        let slopes: Vec<f64> = factors
            .iter()
            .map(|&f| {
                let mut e = engine(f);
                e.solve_point(row[17], &guess)?;
                Some(e.magnetization_b * ERG_PER_MEV_FM3 / (b * f))
            })
            .collect::<Option<_>>()?;
        // s = 1 - 4 pi d(dP_m/dB)/dB.
        let s_of = |plus: f64, minus: f64, step: f64| 1.0 - 4.0 * std::f64::consts::PI * (plus - minus) / (2.0 * step * b);
        let dn = rows[i + 1][0] - rows[i - 1][0];
        Some(StabilityPoint {
            n: row[0],
            mu: row[17],
            s: s_of(slopes[0], slopes[1], delta),
            s_half: s_of(slopes[2], slopes[3], 0.5 * delta),
            dp_dn: if dn != 0.0 { (rows[i + 1][2] - rows[i - 1][2]) / dn } else { f64::NAN },
        })
    });
    points.collect()
}

/// Barra de progresso por tarefa concluída.
fn progress(total: usize) -> indicatif::ProgressBar {
    let bar = indicatif::ProgressBar::new(total as u64);
    if let Ok(style) = indicatif::ProgressStyle::with_template("  [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} (resta ~{eta})") {
        bar.set_style(style);
    }
    bar
}

/// `study b [--models GM1] [--bmin 1e14] [--bmax 1e20] [--per-decade 8]
///  [--profiles constante,bdd] [--points 1500] [--mu-max 3.0] [--no-hyperons]
///  [--no-stability] [--delta 1e-4] [--out results/study_b] [--threads N]`
pub fn run(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["no-hyperons", "no-stability"])?;
    let (b_min, b_max) = (args.f64_or("bmin", 1e14)?, args.f64_or("bmax", 1e20)?);
    let per_decade = args.usize_or("per-decade", 8)?;
    let delta = args.f64_or("delta", 1e-4)?;
    if b_min <= 0.0 || b_max <= b_min || per_decade == 0 || !(1e-6..=1e-2).contains(&delta) {
        return Err("use 0 < --bmin < --bmax, --per-decade >= 1 e 1e-6 <= --delta <= 1e-2".into());
    }
    let config = Config {
        points: args.usize_or("points", 1500)?,
        mu_max: args.f64_or("mu-max", 3.0)?,
        hyperons: !args.switch("no-hyperons"),
    };
    let profiles: Vec<Profile> = args
        .value("profiles")
        .unwrap_or("constante,bdd")
        .split(',')
        .map(|p| match p.trim() {
            "constante" | "constant" => Ok(Profile::Constant),
            "bdd" => Ok(Profile::Bdd),
            other => Err(format!("--profiles: perfil desconhecido '{other}' (use constante, bdd)")),
        })
        .collect::<Result<_, _>>()?;
    let out = args.value("out").unwrap_or("results/study_b").to_string();
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads()?).build_global();
    create_dir(&out)?;

    let steps = ((b_max / b_min).log10() * per_decade as f64).round() as usize;
    let mut fields = vec![0.0];
    fields.extend((0..=steps).map(|k| b_min * 10f64.powf(k as f64 / per_decade as f64)));
    let models = args.models(&["GM1"])?;

    println!(
        "study b: {} modelo(s), {} valores de B (0 e {b_min:.1e}..{b_max:.1e} G, {per_decade}/década), perfis: {}, hyperons = {}",
        models.len(),
        fields.len(),
        profiles.iter().map(|p| p.label()).collect::<Vec<_>>().join(","),
        config.hyperons
    );

    // 1. Estrelas e cobertura.
    let stars_path = format!("{out}/stars.csv");
    let mut stars_csv = fs::File::create(&stars_path).map_err(|e| format!("{stars_path}: {e}"))?;
    writeln!(
        stars_csv,
        "model,profile,B_G,hyperons,rows,n_max_over_n0,termination,anomalous,topology,m_max_Msun,r_max_km,nc_over_n0,r14_km,lambda14,true_maximum,pperp_nonmonotonic_fraction,pperp_negative_rows"
    )
    .map_err(|e| e.to_string())?;
    let mut summary = Vec::new();
    for (model_name, model) in &models {
        let jobs: Vec<(Profile, f64)> =
            profiles.iter().flat_map(|&p| fields.iter().map(move |&b| (p, b))).collect();
        println!("\n[{model_name}] estrelas: {} EoS...", jobs.len());
        let bar = progress(jobs.len());
        // with_max_len(1): cada (perfil, B) é uma tarefa; uma thread livre pega
        // a próxima assim que termina (sem blocos contíguos de B).
        let results: Vec<(Profile, f64, Vec<String>, bool, f64)> = jobs
            .par_iter()
            .with_max_len(1)
            .map(|&(profile, b)| {
                let eos = solve(*model, b, profile, &config);
                let lines = star_lines(model_name, profile, b, &eos, config.hyperons);
                let n_max = eos.rows.iter().map(|r| r[0]).fold(0.0, f64::max);
                bar.inc(1);
                (profile, b, lines, eos.anomalous, n_max)
            })
            .collect();
        bar.finish_and_clear();
        for (_, _, lines, _, _) in &results {
            for line in lines {
                writeln!(stars_csv, "{line}").map_err(|e| e.to_string())?;
            }
        }
        // Onde quebra: primeiro B com término anômalo e primeiro B sem estrela.
        for &profile in &profiles {
            let of_profile: Vec<&(Profile, f64, Vec<String>, bool, f64)> =
                results.iter().filter(|r| r.0 == profile).collect();
            let first_anomalous = of_profile.iter().find(|r| r.3).map(|r| (r.1, r.4));
            let first_no_star = of_profile
                .iter()
                .find(|r| r.2.iter().all(|l| l.split(',').nth(9).is_some_and(|m| m == "NaN")))
                .map(|r| r.1);
            summary.push(format!(
                "  {model_name} | {:<9} | primeiro término anômalo: {} | primeiro B sem estrela: {}",
                profile.label(),
                first_anomalous.map_or("nenhum".into(), |(b, n)| format!("B = {b:.2e} G (EoS até {n:.2} n0)")),
                first_no_star.map_or("nenhum".into(), |b| format!("B = {b:.2e} G")),
            ));
        }
    }
    println!("\nResultados estelares em {stars_path}");

    // 2. Estabilidade local (campo constante).
    if !args.switch("no-stability") {
        let path = format!("{out}/stability.csv");
        let mut csv = fs::File::create(&path).map_err(|e| format!("{path}: {e}"))?;
        writeln!(csv, "model,B_G,n_over_n0,mu_n,s_magnetic,s_magnetic_half_step,magnetically_unstable,dpperp_dn")
            .map_err(|e| e.to_string())?;
        for (model_name, model) in &models {
            let positive: Vec<f64> = fields.iter().copied().filter(|&b| b > 0.0).collect();
            println!("\n[{model_name}] estabilidade: {} valores de B, 4 pontos extras por densidade (delta = {delta:e} e delta/2)...", positive.len());
            let bar = progress(positive.len());
            let maps: Vec<(f64, Vec<StabilityPoint>)> = positive
                .par_iter()
                .with_max_len(1)
                .map(|&b| {
                    let points = stability_points(*model, b, delta, &config);
                    bar.inc(1);
                    (b, points)
                })
                .collect();
            bar.finish_and_clear();
            for (b, points) in &maps {
                let unstable = points.iter().filter(|p| p.magnetically_unstable()).count();
                let mechanical = points.iter().filter(|p| p.dp_dn < 0.0).count();
                if unstable + mechanical > 0 {
                    summary.push(format!(
                        "  {model_name} | B = {b:.2e} G | instável: magnético {unstable}/{n}, mecânico (dP_perp/dn < 0) {mechanical}/{n} pontos",
                        n = points.len()
                    ));
                }
                for p in points {
                    writeln!(
                        csv,
                        "{model_name},{b:.4e},{:.5},{:.6},{:.6e},{:.6e},{},{:.4e}",
                        p.n,
                        p.mu,
                        p.s,
                        p.s_half,
                        p.magnetically_unstable(),
                        p.dp_dn
                    )
                    .map_err(|e| e.to_string())?;
                }
            }
        }
        println!("\nEstabilidade local em {path}");
    }

    println!("\nResumo:");
    for line in summary {
        println!("{line}");
    }
    println!("\nFiguras: python plot_scripts/study_b.py {out}");
    Ok(())
}
