// Estudos de impacto na estrutura estelar.
//
// `study log`: efeito da eletrodinâmica logarítmica, Log(xi), sobre a estrela.
// A NLEM entra só pela energia e pelas tensões do campo (os níveis de Landau
// usam B), com o campo local do perfil `Constant`, isto é, BDD com
// B_surf = 1e15 G e B0 = bg. Para cada (modelo, topologia, B0) a varredura em
// xi é comparada com os dois limites analíticos: Maxwell (xi -> infinito) e
// sem tensão do campo (xi -> 0, pois eps_B = xi^2 ln(1 + x) -> 0).

use std::fs;
use std::io::Write;

use rayon::prelude::*;

use nsrs::constants::{BDD_ALPHAA, BDD_BETAA};
use nsrs::core::magnetic::{B_SURFACE_G, bdd_field_g, magnetic_stress};
use nsrs::core::model::ModelParams;
use nsrs::core::observations::{interpolate_at_mass, stable_branch};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, HadronsMatter, MagneticTopology, NlemModel, Solver};

use crate::cli::{Args, create_dir, f64_list, format_sci};

#[derive(Clone, Copy, PartialEq)]
enum Case {
    NoStress,
    Maxwell,
    Log(f64),
}

impl Case {
    fn label(self) -> &'static str {
        match self {
            Case::NoStress => "sem_tensao",
            Case::Maxwell => "maxwell",
            Case::Log(_) => "log",
        }
    }

    fn xi(self) -> Option<f64> {
        match self {
            Case::Log(xi) => Some(xi),
            _ => None,
        }
    }

    /// Modelo de tensão do campo; `None` quando a tensão não entra na EoS.
    fn nlem(self) -> Option<NlemModel> {
        match self {
            Case::NoStress => None,
            Case::Maxwell => Some(NlemModel::Maxwell),
            Case::Log(xi) => Some(NlemModel::Log(xi)),
        }
    }
}

struct Job {
    model_name: String,
    model: ModelParams,
    topology: MagneticTopology,
    b0: f64,
    case: Case,
}

/// Resultado de uma EoS e da sequência estelar correspondente.
struct Outcome {
    status: String,
    rows: usize,
    /// Linhas do núcleo (n_B > 0) descartadas por P não crescer com eps.
    core_nonmonotonic: usize,
    m_max: Option<f64>,
    r_at_m_max: Option<f64>,
    nc_over_n0: Option<f64>,
    b_center: Option<f64>,
    p_field_over_p_center: Option<f64>,
    r14: Option<f64>,
    lambda14: Option<f64>,
}

fn topology_label(topology: MagneticTopology) -> &'static str {
    match topology {
        MagneticTopology::Anisotropic => "anisotropica",
        MagneticTopology::Isotropic => "isotropica",
    }
}

/// Campo local do perfil `Constant` (BDD com B_surf = 1e15 G), como em
/// `HadronsMatter::assemble_point`.
fn local_field(b0: f64, nb_over_n0: f64) -> f64 {
    if b0 == 0.0 { 0.0 } else { bdd_field_g(B_SURFACE_G, b0, BDD_BETAA, BDD_ALPHAA, nb_over_n0) }
}

/// Pressão do campo usada na TOV (MeV/fm^3) para o caso e a topologia.
fn field_pressure(case: Case, topology: MagneticTopology, b_gauss: f64) -> f64 {
    case.nlem().map_or(0.0, |nlem| magnetic_stress(nlem, b_gauss).pressure(topology))
}

/// x* = B^2/(2 xi^2) acima do qual a pressão do campo na TOV fica negativa
/// no modelo Log (anisotrópica: P_perp; isotrópica: (P_par + 2 P_perp)/3).
/// Obtido por bisseção sobre `magnetic_stress`, para refletir o código.
pub(crate) fn negative_pressure_threshold(topology: MagneticTopology) -> f64 {
    let b = 1e18;
    let sign = |x: f64| field_pressure(Case::Log(b / (2.0 * x).sqrt()), topology, b);
    let (mut lo, mut hi) = (1e-3_f64, 1e3_f64);
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        if sign(mid) > 0.0 { lo = mid } else { hi = mid }
    }
    (lo * hi).sqrt()
}

/// Menor n_B/n0 em que a pressão do campo fica negativa (Log). `Some(0)`: já
/// negativa na superfície; `None`: nunca, pois B < B_surf + B0.
fn negative_pressure_onset(xi: f64, b0: f64, x_star: f64) -> Option<f64> {
    let b_star = xi * (2.0 * x_star).sqrt();
    if b_star <= B_SURFACE_G {
        return Some(0.0);
    }
    let fraction = (b_star - B_SURFACE_G) / b0;
    if b0 == 0.0 || fraction >= 1.0 {
        return None;
    }
    Some((-(1.0 - fraction).ln() / BDD_BETAA).powf(1.0 / BDD_ALPHAA))
}

fn run(job: &Job, points: usize, hyperons: bool, amm: bool, save_dir: Option<&str>) -> Outcome {
    let mut engine = HadronsMatter::new(job.model, job.b0)
        .with_topology(job.topology)
        .with_hyperons(hyperons)
        .with_anomalous_moments(amm)
        .with_limits(0.02, 3.0)
        .with_points(points);
    engine = match job.case.nlem() {
        Some(nlem) => engine.with_nlem(nlem),
        None => engine.with_field_stress(false),
    };
    if let Some(dir) = save_dir {
        let name = match job.case.xi() {
            Some(xi) => format!("xi_{}", format_sci(xi)),
            None => job.case.label().to_string(),
        };
        engine = engine.with_eos_output(format!("{dir}/{name}.dat"));
    }
    let mut solver = Solver::new(EngineMode::Hadrons(engine));
    let rows = solver.solve();

    let mut outcome = Outcome {
        status: String::new(),
        rows: rows.len(),
        core_nonmonotonic: 0,
        m_max: None,
        r_at_m_max: None,
        nc_over_n0: None,
        b_center: None,
        p_field_over_p_center: None,
        r14: None,
        lambda14: None,
    };
    let mut flags = Vec::new();
    if let Some(t) = solver.termination().filter(|t| t.is_anomalous()) {
        flags.push(format!("truncada({t:?})").replace(',', ";"));
    }

    // Núcleo ordenado por eps, só com P crescente (o que a TOV aproveita).
    let mut core: Vec<(f64, f64, f64)> =
        rows.iter().filter(|r| r[0] > 0.0).map(|r| (r[1], r[2], r[0])).collect();
    core.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut kept: Vec<(f64, f64, f64)> = Vec::with_capacity(core.len());
    for point in core {
        if kept.last().is_none_or(|last| point.1 > last.1) {
            kept.push(point);
        } else {
            outcome.core_nonmonotonic += 1;
        }
    }

    let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
    let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let stars = if rows.len() < 5 { Vec::new() } else { generate_star_sequence(&eps, &p, &n, true) };
    let branch = stable_branch(&stars);
    let Some(max) = branch.last() else {
        flags.insert(0, "vazia".into());
        outcome.status = flags.join("+");
        return outcome;
    };
    if branch.len() + 10 > stars.len() {
        flags.push("max_no_fim".into());
    }
    outcome.m_max = Some(max.mass);
    outcome.r_at_m_max = Some(max.radius);
    outcome.r14 = interpolate_at_mass(branch, 1.4, |s| s.radius);
    outcome.lambda14 = interpolate_at_mass(branch, 1.4, |s| s.tidal_deformability);

    // Densidade central da estrela de massa máxima: n_B em P = P_c.
    let pc = max.central_pressure;
    if let Some(k) = (1..kept.len()).find(|&k| kept[k].1 >= pc) {
        let (a, b) = (kept[k - 1], kept[k]);
        let nc = a.2 + (pc - a.1) / (b.1 - a.1) * (b.2 - a.2);
        let bc = local_field(job.b0, nc);
        outcome.nc_over_n0 = Some(nc);
        outcome.b_center = Some(bc);
        outcome.p_field_over_p_center = Some(field_pressure(job.case, job.topology, bc) / pc);
    }

    outcome.status = if flags.is_empty() { "ok".into() } else { flags.join("+") };
    outcome
}

fn opt(value: Option<f64>, digits: usize) -> String {
    value.map_or(String::new(), |v| format!("{v:.digits$}"))
}

fn opt_sci(value: Option<f64>) -> String {
    value.map_or(String::new(), |v| format!("{v:.4e}"))
}

/// `study log <exp_min> <exp_max> <por_década> <B0...>`
/// Saída: results/study_log/summary.csv (ou --out).
pub fn log(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["no-hyperons", "save-eos", "amm"])?;
    let p = &args.positional;
    if p.len() < 4 {
        return Err("uso: study log <exp_min> <exp_max> <por_década> <B0_1> [B0_2 ...]".into());
    }
    let exp_min: f64 = p[0].parse().map_err(|_| "exp_min inválido")?;
    let exp_max: f64 = p[1].parse().map_err(|_| "exp_max inválido")?;
    let per_decade: usize = p[2].parse().map_err(|_| "por_década inválido")?;
    if per_decade == 0 || exp_max < exp_min {
        return Err("use por_década >= 1 e exp_max >= exp_min".into());
    }
    let fields = f64_list(&p[3..], "B0")?;
    if fields.iter().any(|&b| b <= 0.0) {
        return Err("B0 deve ser positivo".into());
    }
    let steps = ((exp_max - exp_min) * per_decade as f64).round() as usize;
    let xis: Vec<f64> =
        (0..=steps).map(|k| 10f64.powf(exp_min + k as f64 / per_decade as f64)).collect();
    let topologies = match args.value("topology").unwrap_or("ambas") {
        "aniso" | "anisotropica" => vec![MagneticTopology::Anisotropic],
        "iso" | "isotropica" => vec![MagneticTopology::Isotropic],
        "ambas" => vec![MagneticTopology::Anisotropic, MagneticTopology::Isotropic],
        other => return Err(format!("--topology '{other}': use aniso, iso ou ambas")),
    };
    let points = args.usize_or("points", 1500)?;
    let hyperons = !args.switch("no-hyperons");
    let amm = args.switch("amm");
    let out = args.value("out").unwrap_or("results/study_log/summary.csv").to_string();

    let mut jobs = Vec::new();
    for (model_name, model) in args.models(&["GM1"])? {
        for &topology in &topologies {
            for &b0 in &fields {
                let cases = [Case::NoStress, Case::Maxwell]
                    .into_iter()
                    .chain(xis.iter().map(|&xi| Case::Log(xi)));
                for case in cases {
                    jobs.push(Job { model_name: model_name.clone(), model, topology, b0, case });
                }
            }
        }
    }
    let save_dir = |job: &Job| {
        format!(
            "output/study_log{}/{}/{}/B0_{}",
            if amm { "_amm" } else { "" },
            job.model_name,
            topology_label(job.topology),
            format_sci(job.b0)
        )
    };
    if args.switch("save-eos") {
        for job in &jobs {
            create_dir(&save_dir(job))?;
        }
    }

    println!(
        "study log: {} EoS ({} valores de xi, {} B0, {} topologia(s)), hyperons = {hyperons}, amm = {amm}",
        jobs.len(),
        xis.len(),
        fields.len(),
        topologies.len()
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads()?)
        .build()
        .map_err(|e| e.to_string())?;
    let outcomes: Vec<Outcome> = pool.install(|| {
        jobs.par_iter()
            .with_max_len(1)
            .map(|job| {
                let dir = args.switch("save-eos").then(|| save_dir(job));
                run(job, points, hyperons, amm, dir.as_deref())
            })
            .collect()
    });

    if let Some(parent) = std::path::Path::new(&out).parent() {
        create_dir(&parent.to_string_lossy())?;
    }
    let mut csv = fs::File::create(&out).map_err(|e| format!("{out}: {e}"))?;
    let mut line = |text: String| writeln!(csv, "{text}").map_err(|e| e.to_string());
    line(
        [
            "model,topology,hyperons,amm,b0_G,case,xi_G,status,rows,core_nonmonotonic_rows",
            "m_max_Msun,r_at_m_max_km,nc_over_n0,b_center_G,x_center,x_threshold",
            "p_field_over_p_center,n_negative_field_pressure_over_n0",
            "r14_km,lambda14,dm_max_vs_maxwell,dm_max_vs_no_stress,dr14_vs_maxwell_km",
        ]
        .join(","),
    )?;

    let thresholds: Vec<(MagneticTopology, f64)> =
        topologies.iter().map(|&t| (t, negative_pressure_threshold(t))).collect();
    let x_star_of = |t: MagneticTopology| thresholds.iter().find(|(tt, _)| *tt == t).unwrap().1;
    let same_group = |a: &Job, b: &Job| {
        a.model_name == b.model_name && a.topology == b.topology && a.b0 == b.b0
    };

    let mut current_group: Option<usize> = None;
    for (i, (job, outcome)) in jobs.iter().zip(&outcomes).enumerate() {
        let reference = |case: Case| {
            jobs.iter()
                .zip(&outcomes)
                .find(|(j, _)| same_group(j, job) && j.case == case)
                .map(|(_, o)| o)
        };
        let (maxwell, no_stress) = (reference(Case::Maxwell), reference(Case::NoStress));
        let delta = |a: Option<f64>, b: Option<f64>| a.zip(b).map(|(a, b)| a - b);
        let x_star = x_star_of(job.topology);
        let x_center = job.case.xi().zip(outcome.b_center).map(|(xi, b)| b * b / (2.0 * xi * xi));
        let onset = job.case.xi().and_then(|xi| negative_pressure_onset(xi, job.b0, x_star));
        let dm_maxwell = delta(outcome.m_max, maxwell.and_then(|o| o.m_max));
        let dm_no_stress = delta(outcome.m_max, no_stress.and_then(|o| o.m_max));
        let dr14 = delta(outcome.r14, maxwell.and_then(|o| o.r14));

        line(format!(
            "{},{},{},{},{:.4e},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            job.model_name,
            topology_label(job.topology),
            hyperons,
            amm,
            job.b0,
            job.case.label(),
            opt_sci(job.case.xi()),
            outcome.status,
            outcome.rows,
            outcome.core_nonmonotonic,
            opt(outcome.m_max, 5),
            opt(outcome.r_at_m_max, 4),
            opt(outcome.nc_over_n0, 4),
            opt_sci(outcome.b_center),
            opt_sci(x_center),
            if job.case.xi().is_some() { format!("{x_star:.6}") } else { String::new() },
            opt_sci(outcome.p_field_over_p_center),
            opt(onset, 4),
            opt(outcome.r14, 4),
            opt(outcome.lambda14, 2),
            opt(dm_maxwell, 5),
            opt(dm_no_stress, 5),
            opt(dr14, 4),
        ))?;

        let group = (0..=i).find(|&k| same_group(&jobs[k], job)).unwrap();
        if current_group != Some(group) {
            current_group = Some(group);
            println!(
                "\n== {} | {} | B0 = {:.2e} G | x* = {x_star:.3} (pressão do campo < 0 para x > x*)",
                job.model_name,
                topology_label(job.topology),
                job.b0
            );
            println!(
                "  {:<22} {:>8} {:>7} {:>7} {:>10} {:>10} {:>9} {:>8} {:>9} {:>9}  status",
                "caso", "M_max", "R_max", "n_c/n0", "B_c [G]", "x_c", "Pcampo/P", "R_1.4", "dM(Maxw)", "dM(s/t)"
            );
        }
        let name = match job.case.xi() {
            Some(xi) => format!("log xi={xi:.2e}"),
            None => job.case.label().to_string(),
        };
        println!(
            "  {:<22} {:>8} {:>7} {:>7} {:>10} {:>10} {:>9} {:>8} {:>9} {:>9}  {}",
            name,
            opt(outcome.m_max, 4),
            opt(outcome.r_at_m_max, 2),
            opt(outcome.nc_over_n0, 2),
            outcome.b_center.map_or(String::new(), |b| format!("{b:.2e}")),
            x_center.map_or(String::new(), |x| format!("{x:.2e}")),
            outcome.p_field_over_p_center.map_or(String::new(), |r| format!("{r:+.1e}")),
            opt(outcome.r14, 3),
            opt(dm_maxwell, 4),
            opt(dm_no_stress, 4),
            outcome.status
        );
    }
    println!("\nResumo gravado em {out}");
    Ok(())
}
