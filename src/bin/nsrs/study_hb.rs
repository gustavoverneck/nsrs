// `study hb`: curva H(B) da matéria a mu_n fixo, para ver onde dH/dB < 0
// (domínios de Condon) com dados do modelo.
//
// Para cada modelo e densidade de referência, escolhe o mu_n da linha da EoS
// (Maxwell, campo constante B_ref) mais próxima de n_B/n0 e varre B numa malha
// logarítmica, resolvendo a matéria a mu_n fixo (continuação do chute de um B
// para o seguinte). Grava a pressão da matéria P_m(mu_n, B), 𝓜B = B dP_m/dB e
// 4 pi 𝓜 (G). O campo total é H = H_vac(B) - 4 pi 𝓜, com H_vac = B (Maxwell)
// ou B/(1 + B^2/2 xi^2) (Log): a parte do vácuo é analítica e fica para o
// script de figura (plot_scripts/hb_domains.py).
//
// Saída: results/study_hb/hb.csv com
// model,n_ref_over_n0,mu_n,B_G,n_over_n0,P_m_MeV_fm3,MB_MeV_fm3,four_pi_M_G

use std::fs;
use std::io::Write;

use rayon::prelude::*;

use nsrs::HadronsMatter;
use nsrs::NlemModel;
use nsrs::constants::M_NUCLEON;

use crate::cli::{Args, create_dir, f64_list, progress, tick};
use crate::study_b::{Config, ERG_PER_MEV_FM3, Profile, solve};

/// `study hb [--models GM1,GM3,FSU2] [--densities 2,4] [--bmin 1e15] [--bmax 1e19]
///  [--per-decade 1000] [--b-ref 1e17] [--no-hyperons] [--amm] [--out results/study_hb]`
pub fn run(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["no-hyperons", "amm"])?;
    let config = Config {
        points: args.usize_or("points", 1500)?,
        mu_max: args.f64_or("mu-max", 3.0)?,
        hyperons: !args.switch("no-hyperons"),
        landau_max: args.usize_or("landau-max", nsrs::constants::MAX_LANDAU_LIMIT)?,
        amm: args.switch("amm"),
        step_tol: crate::study_b::step_tolerance(&args)?,
    };
    let densities = f64_list(
        &args
            .value("densities")
            .unwrap_or("2,4")
            .split(',')
            .map(|s| s.trim().to_string())
            .collect::<Vec<_>>(),
        "--densities",
    )?;
    let (bmin, bmax) = (args.f64_or("bmin", 1e15)?, args.f64_or("bmax", 1e19)?);
    let per_decade = args.usize_or("per-decade", 1000)?;
    let b_ref = args.f64_or("b-ref", 1e17)?;
    if !(bmin > 0.0 && bmax > bmin) {
        return Err("use 0 < --bmin < --bmax".into());
    }
    let models = args.models(&["GM1", "GM3", "FSU2"])?;
    let out = args.value("out").unwrap_or("results/study_hb").to_string();
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads()?)
        .build_global();
    create_dir(&out)?;

    let decades = (bmax / bmin).log10();
    let count = (decades * per_decade as f64).round() as usize + 1;
    let fields: Vec<f64> = (0..count)
        .map(|i| bmin * 10f64.powf(decades * i as f64 / (count - 1) as f64))
        .collect();

    // mu_n e chute inicial de cada densidade, da EoS de Maxwell com B_ref constante.
    let mut jobs = Vec::new();
    for (name, model) in &models {
        let eos = solve(
            *model,
            b_ref,
            Profile::Constant,
            NlemModel::Maxwell,
            &config,
        );
        for &target in &densities {
            let row = eos
                .rows
                .iter()
                .filter(|r| r[0] > 0.0)
                .min_by(|a, b| (a[0] - target).abs().total_cmp(&(b[0] - target).abs()))
                .ok_or(format!("{name}: EoS vazia em B_ref"))?;
            let guess = [
                row[18],
                row[13] / M_NUCLEON,
                row[14] / M_NUCLEON,
                row[15] / M_NUCLEON,
                0.0,
            ];
            jobs.push((name.clone(), *model, target, row[17], guess));
        }
    }
    println!(
        "study hb: {} varredura(s) de {} campos ({bmin:e}..{bmax:e} G), hyperons = {}, amm = {}",
        jobs.len(),
        fields.len(),
        config.hyperons,
        config.amm
    );

    let bar = progress(jobs.len() * fields.len(), "H(B)");
    let sweeps: Vec<Vec<String>> = jobs
        .par_iter()
        .with_max_len(1)
        .map(|(name, model, target, mu, guess)| {
            // Continuação a partir de B_ref nos dois sentidos: o chute vem da linha em B_ref.
            let start = fields.partition_point(|&b| b < b_ref).min(fields.len() - 1);
            let order: Vec<usize> = (start..fields.len()).chain((0..start).rev()).collect();
            let mut lines = vec![String::new(); fields.len()];
            let mut x = *guess;
            for (k, &i) in order.iter().enumerate() {
                if k > 0 && start > 0 && i == start - 1 {
                    x = *guess;
                }
                let b = fields[i];
                let mut engine = HadronsMatter::new(*model, b)
                    .with_hyperons(config.hyperons)
                    .with_anomalous_moments(config.amm)
                    .with_field_stress(false)
                    .with_limits(0.02, config.mu_max)
                    .with_max_landau_limit(config.landau_max);
                if let Some((x_new, result)) = engine.solve_point(*mu, &x) {
                    x = x_new;
                    let p_m = engine.stability_pressure;
                    let mb = engine.magnetization_b;
                    let four_pi_m = 4.0 * std::f64::consts::PI * mb * ERG_PER_MEV_FM3 / b;
                    lines[i] = format!(
                        "{name},{target},{mu:.8},{b:.8e},{:.8},{p_m:.16e},{mb:.10e},{four_pi_m:.10e}",
                        result[0]
                    );
                }
                tick(&bar);
            }
            lines.into_iter().filter(|l| !l.is_empty()).collect()
        })
        .collect();
    bar.finish_and_clear();

    let path = format!("{out}/hb.csv");
    let mut csv = fs::File::create(&path).map_err(|e| format!("{path}: {e}"))?;
    writeln!(
        csv,
        "model,n_ref_over_n0,mu_n,B_G,n_over_n0,P_m_MeV_fm3,MB_MeV_fm3,four_pi_M_G"
    )
    .map_err(|e| e.to_string())?;
    for line in sweeps.iter().flatten() {
        writeln!(csv, "{line}").map_err(|e| e.to_string())?;
    }
    println!("{path}");
    Ok(())
}
