// `study b`: varredura completa em B, de B = 0 até o código deixar de produzir
// estrelas, com propriedades estelares, cobertura da EoS e estabilidade local.
//
// Para cada eletrodinâmica (`--nlem`), cada B e cada perfil de campo:
// - estrelas com crosta BPS nas duas pressões da TOV: P_perp (topologia
//   anisotrópica) e a média isotrópica (P_par + 2 P_perp)/3;
// - cobertura da EoS (n_B máximo, motivo do término);
// - estabilidade mecânica: P_perp deve crescer com n_B;
// - estabilidade magnética (só campo constante, física local): convexidade da
//   energia total em B a mu_n fixo,
//       s = 4 pi d^2 eps_campo/dB^2 - 4 pi d^2 P_m/dB^2 = f_vac(B) - 4 pi d(MB/B)/dB
//   (unidades gaussianas; f_vac = 1 em Maxwell, (1 - x)/(1 + x)^2 na
//   eletrodinâmica logarítmica). A NLEM não muda a matéria (acoplamento
//   mínimo), então a curvatura da matéria é calculada uma vez por B e vale
//   para todas as eletrodinâmicas. s < 0: instável à formação de domínios
//   magnéticos. É o complemento de Schur da hessiana de eps(n_B, B), logo o
//   critério a mu fixo já inclui a estabilidade em n_B. A derivada segunda é
//   tomada por diferença central com passo relativo `delta` em B e repetida
//   com delta/2; o ponto só conta como instável se s < 0 nos dois (o passo
//   precisa resolver as oscilações de de Haas-van Alphen e o sinal não pode
//   depender dele).

use std::fs;
use std::io::Write;

use rayon::prelude::*;

use nsrs::constants::{MAX_LANDAU_LIMIT, M_NUCLEON, RESULTS_SIZE};
use nsrs::core::model::ModelParams;
use nsrs::core::observations::{interpolate_at_mass, stable_branch};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::core::magnetic::MagneticStress;
use nsrs::{EngineMode, FieldProfile, HadronsMatter, NlemModel, Solver};

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
    landau_max: usize,
    /// Momentos magnéticos anômalos dos bárions (`--amm`).
    amm: bool,
}

/// Rótulo da eletrodinâmica nos CSVs: maxwell, log:<xi>, modmax:<gamma>.
fn nlem_label(nlem: NlemModel) -> String {
    match nlem {
        NlemModel::Maxwell => "maxwell".into(),
        NlemModel::Log(xi) => format!("log:{xi:.2e}"),
        NlemModel::Modmax(gamma) => format!("modmax:{gamma}"),
    }
}

/// `--nlem maxwell,log:1e17,...` (xi em Gauss).
fn parse_nlem(text: &str) -> Result<Vec<NlemModel>, String> {
    text.split(',')
        .map(|item| {
            let item = item.trim();
            let value = |v: &str| v.parse::<f64>().map_err(|_| format!("--nlem: número inválido em '{item}'"));
            match item.split_once(':') {
                None if item == "maxwell" => Ok(NlemModel::Maxwell),
                Some(("log", v)) => {
                    let xi = value(v)?;
                    if xi > 0.0 { Ok(NlemModel::Log(xi)) } else { Err(format!("--nlem: xi deve ser > 0 em '{item}'")) }
                }
                Some(("modmax", v)) => Ok(NlemModel::Modmax(value(v)?)),
                _ => Err(format!("--nlem: '{item}' (use maxwell, log:<xi em G> ou modmax:<gamma>)")),
            }
        })
        .collect()
}

/// EoS resolvida, com 𝓜B, a pressão sem magnetização e as tensões do campo
/// por linha.
struct Eos {
    rows: Vec<Row>,
    magnetization_b: Vec<f64>,
    stability_pressure: Vec<f64>,
    field_stress: Vec<MagneticStress>,
    termination: String,
    anomalous: bool,
}

fn solve(model: ModelParams, b: f64, profile: Profile, nlem: NlemModel, config: &Config) -> Eos {
    let mut engine = HadronsMatter::new(model, b)
        .with_nlem(nlem)
        .with_hyperons(config.hyperons)
        .with_anomalous_moments(config.amm)
        .with_limits(0.02, config.mu_max)
        .with_points(config.points)
        .with_max_landau_limit(config.landau_max);
    if profile == Profile::Bdd && b > 0.0 {
        engine = engine.with_field_profile(FieldProfile::bdd(b));
    }
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();
    let termination = solver.termination();
    Eos {
        magnetization_b: solver.diagnostics().iter().map(|d| d.magnetization_b).collect(),
        stability_pressure: solver.diagnostics().iter().map(|d| d.stability_pressure).collect(),
        field_stress: solver.diagnostics().iter().map(|d| d.field_stress).collect(),
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

/// Pressão isotrópica a partir da EoS anisotrópica, para qualquer NLEM:
/// P_iso = P_par,m - (2/3) 𝓜B + (P_par,c + 2 P_perp,c)/3, com
/// P_par,m = P_estab - P_perp,c (a pressão de estabilidade é P_par,m + P_perp,c).
/// Em Maxwell, P_par,c = -eps_B e P_perp,c = eps_B: P_iso = P_par,m - (2/3) 𝓜B + eps_B/3.
fn isotropic_pressure(eos: &Eos) -> Vec<f64> {
    eos.magnetization_b
        .iter()
        .zip(&eos.stability_pressure)
        .zip(&eos.field_stress)
        .map(|((mb, stab), f)| {
            (stab - f.p_perpendicular) - 2.0 / 3.0 * mb + (f.p_parallel + 2.0 * f.p_perpendicular) / 3.0
        })
        .collect()
}

/// Linha do stars.csv para um (modelo, NLEM, perfil, B).
fn star_lines(model_name: &str, nlem: &str, profile: Profile, b: f64, eos: &Eos, hyperons: bool) -> Vec<String> {
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
                "{model_name},{nlem},{},{b:.4e},{hyperons},{},{n_max:.4},{},{},{topology},{:.5},{:.4},{:.4},{:.4},{:.2},{},{nonmonotonic:.4},{negative}",
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

/// Ponto da estabilidade local: n_B/n0, mu_n, curvatura da matéria
/// c_m = 4 pi d^2 P_m/dB^2 com passo delta e delta/2, e dP_perp/dn_B para cada
/// eletrodinâmica (mesma ordem de `--nlem`).
struct StabilityPoint {
    n: f64,
    mu: f64,
    c_m: f64,
    c_m_half: f64,
    dp_dn: Vec<f64>,
}

impl StabilityPoint {
    /// s = f_vac - c_m com os dois passos.
    fn s(&self, f_vac: f64) -> (f64, f64) {
        (f_vac - self.c_m, f_vac - self.c_m_half)
    }

    /// Instabilidade magnética robusta: s < 0 com os dois passos.
    fn magnetically_unstable(&self, f_vac: f64) -> bool {
        let (s, s_half) = self.s(f_vac);
        s < 0.0 && s_half < 0.0
    }
}

/// dP_perp/dn_B da EoS no ponto da malha mais próximo de `n` (diferença central).
fn slope_at(rows: &[Row], n: f64) -> f64 {
    if rows.len() < 3 {
        return f64::NAN;
    }
    let i = (1..rows.len() - 1)
        .min_by(|&a, &b| (rows[a][0] - n).abs().total_cmp(&(rows[b][0] - n).abs()))
        .unwrap_or(1);
    let dn = rows[i + 1][0] - rows[i - 1][0];
    if dn != 0.0 { (rows[i + 1][2] - rows[i - 1][2]) / dn } else { f64::NAN }
}

fn stability_points(model: ModelParams, b: f64, delta: f64, nlems: &[NlemModel], config: &Config) -> Vec<StabilityPoint> {
    // A matéria não depende da NLEM: a curvatura c_m vem de uma única EoS
    // (Maxwell). A estabilidade mecânica usa a P_perp de cada eletrodinâmica.
    let center = solve(model, b, Profile::Constant, NlemModel::Maxwell, config);
    let per_nlem: Vec<Option<Eos>> = nlems
        .iter()
        .map(|&nlem| (nlem != NlemModel::Maxwell).then(|| solve(model, b, Profile::Constant, nlem, config)))
        .collect();
    // Cada ponto é resolvido de novo no mesmo mu_n com B(1 ± delta) e
    // B(1 ± delta/2), a partir da solução central. (A malha do solver é
    // adaptativa: execuções separadas não caem nos mesmos mu_n.)
    let engine = |factor: f64| {
        HadronsMatter::new(model, b * factor)
            .with_hyperons(config.hyperons)
            .with_anomalous_moments(config.amm)
            .with_limits(0.02, config.mu_max)
            .with_points(config.points)
            .with_max_landau_limit(config.landau_max)
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
        // c_m = 4 pi d(dP_m/dB)/dB.
        let c_of = |plus: f64, minus: f64, step: f64| 4.0 * std::f64::consts::PI * (plus - minus) / (2.0 * step * b);
        let dn = rows[i + 1][0] - rows[i - 1][0];
        let maxwell_slope = if dn != 0.0 { (rows[i + 1][2] - rows[i - 1][2]) / dn } else { f64::NAN };
        Some(StabilityPoint {
            n: row[0],
            mu: row[17],
            c_m: c_of(slopes[0], slopes[1], delta),
            c_m_half: c_of(slopes[2], slopes[3], 0.5 * delta),
            dp_dn: per_nlem
                .iter()
                .map(|eos| eos.as_ref().map_or(maxwell_slope, |e| slope_at(&e.rows, row[0])))
                .collect(),
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
///  [--no-stability] [--delta 1e-4] [--landau-max 20000] [--amm] [--nlem maxwell,log:1e17]
///  [--out results/study_b] [--threads N]`
pub fn run(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["no-hyperons", "no-stability", "amm"])?;
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
        landau_max: args.usize_or("landau-max", MAX_LANDAU_LIMIT)?,
        amm: args.switch("amm"),
    };
    if config.landau_max == 0 {
        return Err("--landau-max deve ser >= 1".into());
    }
    let nlems = parse_nlem(args.value("nlem").unwrap_or("maxwell"))?;
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
        "study b: {} modelo(s), {} valores de B (0 e {b_min:.1e}..{b_max:.1e} G, {per_decade}/década), perfis: {}, NLEM: {}, hyperons = {}, níveis de Landau <= {}",
        models.len(),
        fields.len(),
        profiles.iter().map(|p| p.label()).collect::<Vec<_>>().join(","),
        nlems.iter().map(|&n| nlem_label(n)).collect::<Vec<_>>().join(","),
        config.hyperons,
        config.landau_max
    );

    // 1. Estrelas e cobertura.
    let stars_path = format!("{out}/stars.csv");
    let mut stars_csv = fs::File::create(&stars_path).map_err(|e| format!("{stars_path}: {e}"))?;
    writeln!(
        stars_csv,
        "model,nlem,profile,B_G,hyperons,rows,n_max_over_n0,termination,anomalous,topology,m_max_Msun,r_max_km,nc_over_n0,r14_km,lambda14,true_maximum,pperp_nonmonotonic_fraction,pperp_negative_rows"
    )
    .map_err(|e| e.to_string())?;
    let mut summary = Vec::new();
    for (model_name, model) in &models {
        let mut jobs: Vec<(NlemModel, Profile, f64)> = Vec::new();
        for &nlem in &nlems {
            for &profile in &profiles {
                jobs.extend(fields.iter().map(|&b| (nlem, profile, b)));
            }
        }
        println!("\n[{model_name}] estrelas: {} EoS...", jobs.len());
        let bar = progress(jobs.len());
        // with_max_len(1): cada (NLEM, perfil, B) é uma tarefa; uma thread livre
        // pega a próxima assim que termina (sem blocos contíguos de B).
        let results: Vec<(NlemModel, Profile, f64, Vec<String>, bool, f64)> = jobs
            .par_iter()
            .with_max_len(1)
            .map(|&(nlem, profile, b)| {
                let eos = solve(*model, b, profile, nlem, &config);
                let lines = star_lines(model_name, &nlem_label(nlem), profile, b, &eos, config.hyperons);
                let n_max = eos.rows.iter().map(|r| r[0]).fold(0.0, f64::max);
                bar.inc(1);
                (nlem, profile, b, lines, eos.anomalous, n_max)
            })
            .collect();
        bar.finish_and_clear();
        for (_, _, _, lines, _, _) in &results {
            for line in lines {
                writeln!(stars_csv, "{line}").map_err(|e| e.to_string())?;
            }
        }
        // Onde quebra: primeiro B com término anômalo e primeiro B sem estrela.
        for (&nlem, &profile) in nlems.iter().flat_map(|n| profiles.iter().map(move |p| (n, p))) {
            let of_profile: Vec<&(NlemModel, Profile, f64, Vec<String>, bool, f64)> =
                results.iter().filter(|r| r.0 == nlem && r.1 == profile).collect();
            let first_anomalous = of_profile.iter().find(|r| r.4).map(|r| (r.2, r.5));
            let first_no_star = of_profile
                .iter()
                .find(|r| r.3.iter().all(|l| l.split(',').nth(10).is_some_and(|m| m == "NaN")))
                .map(|r| r.2);
            summary.push(format!(
                "  {model_name} | {:<12} | {:<9} | primeiro término anômalo: {} | primeiro B sem estrela: {}",
                nlem_label(nlem),
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
        writeln!(
            csv,
            "model,nlem,B_G,n_over_n0,mu_n,f_vac,c_matter,s_magnetic,s_magnetic_half_step,magnetically_unstable,dpperp_dn"
        )
            .map_err(|e| e.to_string())?;
        for (model_name, model) in &models {
            let positive: Vec<f64> = fields.iter().copied().filter(|&b| b > 0.0).collect();
            println!("\n[{model_name}] estabilidade: {} valores de B, 4 pontos extras por densidade (delta = {delta:e} e delta/2)...", positive.len());
            let bar = progress(positive.len());
            let maps: Vec<(f64, Vec<StabilityPoint>)> = positive
                .par_iter()
                .with_max_len(1)
                .map(|&b| {
                    let points = stability_points(*model, b, delta, &nlems, &config);
                    bar.inc(1);
                    (b, points)
                })
                .collect();
            bar.finish_and_clear();
            for (k, &nlem) in nlems.iter().enumerate() {
                let label = nlem_label(nlem);
                for (b, points) in &maps {
                    let f_vac = nlem.vacuum_curvature(*b);
                    let unstable = points.iter().filter(|p| p.magnetically_unstable(f_vac)).count();
                    let mechanical = points.iter().filter(|p| p.dp_dn[k] < 0.0).count();
                    if unstable + mechanical > 0 {
                        summary.push(format!(
                            "  {model_name} | {label:<12} | B = {b:.2e} G | instável: magnético {unstable}/{n}, mecânico (dP_perp/dn < 0) {mechanical}/{n} pontos",
                            n = points.len()
                        ));
                    }
                    for p in points {
                        let (s, s_half) = p.s(f_vac);
                        writeln!(
                            csv,
                            "{model_name},{label},{b:.4e},{:.5},{:.6},{f_vac:.6e},{:.6e},{s:.6e},{s_half:.6e},{},{:.4e}",
                            p.n,
                            p.mu,
                            p.c_m,
                            p.magnetically_unstable(f_vac),
                            p.dp_dn[k]
                        )
                        .map_err(|e| e.to_string())?;
                    }
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
