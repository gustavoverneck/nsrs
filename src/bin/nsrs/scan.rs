// Varreduras de parâmetros (substituem os binários b, nlem_log, nlem_modmax e
// magtop). Caminhos de saída e malhas em mu_n iguais aos originais, para que
// os scripts de plot_scripts/ continuem funcionando.

use std::fs;
use std::io::{BufRead, BufReader, Write};

use nsrs::constants::BCE_G;
use nsrs::core::plotting::Artist;
use nsrs::core::tov_solver::generate_mr_curve;
use nsrs::{EngineMode, HadronsMatter, MagneticTopology, NlemModel, Solver};

use crate::cli::{Args, create_dir, f64_list, format_sci, model};

const ALL_MODELS: [&str; 3] = ["GM1", "GM3", "FSU2"];

/// `scan b`: B = 0 e `points - 1` valores em escala log de bmin (padrão:
/// campo crítico do elétron) até bmax.
/// Saída: output/b/<modelo>/B_<b>/default/eos.dat
pub fn b(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let points = args.usize_or("points", 100)?;
    let (b_min, b_max) = (args.f64_or("bmin", BCE_G)?, args.f64_or("bmax", 3.0e18)?);
    if points < 2 || b_min <= 0.0 || b_max <= b_min {
        return Err("use --points >= 2 e 0 < --bmin < --bmax".into());
    }
    let mut fields = vec![0.0];
    let n = points - 1;
    for i in 0..n {
        let t = if n == 1 { 0.0 } else { i as f64 / (n - 1) as f64 };
        fields.push(10f64.powf(b_min.log10() + t * (b_max.log10() - b_min.log10())));
    }

    for (name, params) in args.models(&ALL_MODELS)? {
        println!("\nModelo={name} | Varrendo {points} valores de B (0, {b_min:.2e}..{b_max:.2e} G)...");
        let mut engines = Vec::new();
        for &b_field in &fields {
            let dir = format!("output/b/{name}/B_{:.2e}/default", b_field);
            create_dir(&dir)?;
            let engine = HadronsMatter::new(params, b_field)
                .with_limits(0.01, 2.0)
                .with_points(1201)
                .with_eos_output(format!("{dir}/eos.dat"));
            engines.push(EngineMode::Hadrons(engine));
        }
        Solver::solve_parallel(engines, args.threads()?);
    }
    println!("\nConcluído. Dados em output/b/");
    Ok(())
}

/// `scan log <exp_min> <exp_max> <pontos_por_década> <B...>`: NLEM
/// logarítmica com xi (Gauss) em escala log.
/// Saída: output/nlem_log/<modelo>/B_<b>/default/csi_<xi>/eos.dat
pub fn log(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let p = &args.positional;
    if p.len() < 4 {
        return Err("uso: scan log <exp_min> <exp_max> <pontos_por_década> <B1> [B2 ...]".into());
    }
    let exp_min: i32 = p[0].parse().map_err(|_| "exp_min inválido")?;
    let exp_max: i32 = p[1].parse().map_err(|_| "exp_max inválido")?;
    let per_decade: usize = p[2].parse().map_err(|_| "pontos_por_década inválido")?;
    if per_decade == 0 || exp_max < exp_min {
        return Err("use pontos_por_década >= 1 e exp_max >= exp_min".into());
    }
    let fields = f64_list(&p[3..], "B")?;
    let mut xis: Vec<f64> = (exp_min..exp_max)
        .flat_map(|e| (0..per_decade).map(move |k| 10f64.powf(e as f64 + k as f64 / per_decade as f64)))
        .collect();
    xis.push(10f64.powi(exp_max));

    for (name, params) in args.models(&ALL_MODELS)? {
        for &b_field in &fields {
            println!("\nModelo={name} | B = {:.2e} G | Varrendo {} valores de ξ...", b_field, xis.len());
            let base = format!("output/nlem_log/{name}/B_{:.2e}/default", b_field);
            let mut engines = Vec::new();
            for &xi in &xis {
                let dir = format!("{base}/csi_{:.2e}", xi);
                create_dir(&dir)?;
                let engine = HadronsMatter::new(params, b_field)
                    .with_nlem(NlemModel::Log(xi))
                    .with_limits(0.01, 2.0)
                    .with_points(2000)
                    .with_eos_output(format!("{dir}/eos.dat"));
                engines.push(EngineMode::Hadrons(engine));
            }
            Solver::solve_parallel(engines, args.threads()?);
        }
    }
    println!("\nConcluído. Dados em output/nlem_log/");
    Ok(())
}

/// `scan modmax <B...>`: ModMax com gamma = {1..9} x 10^{-10..-1} e uma EoS
/// de referência sem NLEM.
/// Saída: output/modmax/<modelo>/B_<b>/{summary.csv, eos_baseline.dat, eos_csi_<gamma>.dat}
pub fn modmax(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    if args.positional.is_empty() {
        return Err("uso: scan modmax <B1> [B2 ...]".into());
    }
    let fields = f64_list(&args.positional, "B")?;
    let gammas: Vec<f64> =
        (-10..=-1).flat_map(|e| (1..=9).map(move |i| i as f64 * 10f64.powi(e))).collect();

    for (name, params) in args.models(&ALL_MODELS)? {
        for &b_field in &fields {
            let base = format!("output/modmax/{name}/B_{}", format_sci(b_field));
            create_dir(&base)?;
            let mut summary = fs::File::create(format!("{base}/summary.csv")).map_err(|e| e.to_string())?;
            writeln!(summary, "label,csi,b_field,eos_file").map_err(|e| e.to_string())?;

            let mut engines = vec![EngineMode::Hadrons(
                HadronsMatter::new(params, b_field)
                    .with_limits(0.02, 2.0)
                    .with_points(1201)
                    .with_eos_output(format!("{base}/eos_baseline.dat")),
            )];
            writeln!(summary, "hadrons_baseline,,,eos_baseline.dat").map_err(|e| e.to_string())?;
            for &gamma in &gammas {
                let file = format!("eos_csi_{}.dat", format_sci(gamma));
                engines.push(EngineMode::Hadrons(
                    HadronsMatter::new(params, b_field)
                        .with_nlem(NlemModel::Modmax(gamma))
                        .with_limits(0.02, 2.0)
                        .with_points(1201)
                        .with_eos_output(format!("{base}/{file}")),
                ));
                writeln!(summary, "modmax,{gamma:.6e},{b_field:.6e},{file}").map_err(|e| e.to_string())?;
            }
            println!("\nModelo={name} | B = {} G | Varrendo {} valores de γ...", format_sci(b_field), gammas.len());
            Solver::solve_parallel(engines, args.threads()?);
        }
    }
    println!("\nConcluído. Dados em output/modmax/");
    Ok(())
}

/// Densidades por espécie (colunas 3..12) de uma linha da EoS.
const SPECIES: [(&str, usize); 10] = [
    ("e-", 3), ("mu-", 4), ("n", 5), ("p", 6), ("L0", 7),
    ("S-", 8), ("S0", 9), ("S+", 10), ("X-", 11), ("X0", 12),
];

/// `scan topology <modelo> <B...> [--prefix TAG] [--plot-only]`: compara as
/// topologias isotrópica e anisotrópica; gera EoS, M-R e populações.
/// Saída: output/magtop/<modelo>/<tag_>B_<b>/{isotropic,anisotropic}/eos.dat,
/// results/magtop/<modelo>/*.svg e results/magtop/<tag_>summary_topology_<modelo>.csv
pub fn topology(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["plot-only"])?;
    if args.positional.len() < 2 {
        return Err("uso: scan topology <GM1|GM3|FSU2> <B1> [B2 ...] [--prefix TAG] [--plot-only]".into());
    }
    let model_name = args.positional[0].clone();
    let params = model(&model_name)?;
    let fields = f64_list(&args.positional[1..], "B")?;
    let prefix = args.value("prefix").map(|p| format!("{p}_")).unwrap_or_default();
    let plots_dir = format!("results/magtop/{model_name}");
    create_dir(&plots_dir)?;
    let mut summary = vec!["b_field,topology,max_mass,radius_at_max".to_string()];

    for &b_field in &fields {
        let b_string = format!("{:.2e}", b_field);
        let base_dir = format!("output/magtop/{model_name}/{prefix}B_{b_string}");
        let tag = b_string.replace('+', "p").replace('-', "m");
        let svg = |kind: &str| format!("{plots_dir}/{prefix}{kind}_{model_name}_B_{tag}.svg");

        if !args.switch("plot-only") {
            println!("Calculando EoS para B = {b_string} G...");
            let engine = |topology| {
                EngineMode::Hadrons(
                    HadronsMatter::new(params, b_field)
                        .with_topology(topology)
                        .with_limits(0.01, 2.2)
                        .with_points(1500),
                )
            };
            let results = Solver::solve_parallel(
                vec![engine(MagneticTopology::Isotropic), engine(MagneticTopology::Anisotropic)],
                args.threads()?,
            );
            for (rows, folder) in results.iter().zip(["isotropic", "anisotropic"]) {
                let dir = format!("{base_dir}/{folder}");
                create_dir(&dir)?;
                Solver::write_eos(rows, &format!("{dir}/eos.dat")).map_err(|e| e.to_string())?;
            }
        }

        let mut eos_plot = Artist::new(&svg("eos_topology"), &format!("Equation of State (Topology) - {model_name} | B = {b_string} G"))
            .with_x_label("Energy Density \u{03B5} [MeV/fm\u{00B3}]")
            .with_y_label("Pressure P [MeV/fm\u{00B3}]")
            .autoscale()
            .with_log_scale();
        let mut mr_plot = Artist::new(&svg("mr_topology"), &format!("Mass-Radius (Topology) - {model_name} | B = {b_string} G"))
            .with_x_label("Radius [km]")
            .with_y_label("Mass [M\u{2299}]")
            .with_x_range(8.0, 16.0);
        let mut any_curve = false;

        for (folder, label, pop_kind) in [("isotropic", "Iso", "pop_iso"), ("anisotropic", "Aniso", "pop_aniso")] {
            let Ok(file) = fs::File::open(format!("{base_dir}/{folder}/eos.dat")) else { continue };
            let rows: Vec<Vec<f64>> = BufReader::new(file)
                .lines()
                .map_while(Result::ok)
                .map(|l| l.split_whitespace().filter_map(|s| s.parse().ok()).collect::<Vec<f64>>())
                .filter(|v| v.len() > 12)
                .collect();

            let physical: Vec<&Vec<f64>> = rows.iter().filter(|r| r[1] > 0.0 && r[2] > 0.0).collect();
            if physical.len() > 10 {
                let eps: Vec<f64> = physical.iter().map(|r| r[1]).collect();
                let p: Vec<f64> = physical.iter().map(|r| r[2]).collect();
                let rho: Vec<f64> = physical.iter().map(|r| r[0]).collect();
                let (masses, radii, _, _) = generate_mr_curve(&eps, &p, &rho, false);
                let curve_label = format!("{b_string} ({label})");
                eos_plot = eos_plot.add_curve(&eps, &p, &curve_label);
                mr_plot = mr_plot.add_curve(&radii, &masses, &curve_label);
                any_curve = true;
                let (m_max, r_max) = masses
                    .iter()
                    .zip(&radii)
                    .filter(|(_, r)| **r < 16.0)
                    .fold((0.0, 0.0), |acc, (&m, &r)| if m > acc.0 { (m, r) } else { acc });
                summary.push(format!("{b_string},{label},{m_max:.4},{r_max:.4}"));
            }

            if rows.len() > 10 {
                let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
                let mut pop = Artist::new(
                    &svg(pop_kind),
                    &format!("Particle Population ({label}) - {model_name} | B = {b_string} G"),
                )
                .with_x_label("n / n0")
                .with_y_label("n_i [fm^-3]")
                .autoscale();
                for (species, col) in SPECIES {
                    let y: Vec<f64> = rows.iter().map(|r| r[col]).collect();
                    pop = pop.add_curve(&n, &y, species);
                }
                pop.plot().ok();
            }
        }

        if any_curve {
            eos_plot.plot().ok();
            mr_plot.plot().ok();
        } else {
            eprintln!("Aviso: nenhum dado válido encontrado para B = {b_string} G");
        }
    }

    let summary_path = format!("results/magtop/{prefix}summary_topology_{model_name}.csv");
    fs::write(&summary_path, summary.join("\n") + "\n").map_err(|e| e.to_string())?;
    println!("\nConcluído. Gráficos em results/magtop/");
    Ok(())
}
