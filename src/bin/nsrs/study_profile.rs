// `study profile`: grandezas magnéticas ponto a ponto da EoS e ao longo do
// raio de estrelas escolhidas, para comparar os perfis de campo constante e
// BDD sob cada eletrodinâmica.
//
// Para cada (modelo, B0, perfil, NLEM) grava:
// - eos.csv: por linha da EoS (topologia anisotrópica, P = P_perp), o campo
//   dos níveis de Landau e o da energia do campo, energia e pressão totais e
//   do campo, 𝓜B, 4 pi 𝓜/B, a permeabilidade B/H, f_vac, a curvatura da
//   matéria c_m e s = dH/dB = f_vac - c_m (critério de `study b`, agora com o
//   campo local de cada ponto: B0 no perfil constante, B(n_B) no BDD);
// - radial.csv: r, m, P e n_B/n0 da estrela de massa máxima e da de 1.4 M_sun
//   (crosta BPS incluída; n_B = NaN na crosta).
//
// H = B (H/B)_vac - 4 pi 𝓜, logo B/H = 1 / ((H/B)_vac - 4 pi 𝓜/B), com
// 4 pi 𝓜/B = 𝓜B / (2 eps_Maxwell) e eps_Maxwell = B^2/8pi.
//
// No perfil constante a energia e as tensões do campo seguem o perfil BDD
// com B0 = B (comportamento padrão do código, `HadronsMatter::assemble_point`);
// só os níveis de Landau usam o B constante. A estabilidade e a
// permeabilidade usam o campo que a matéria sente.

use std::collections::HashMap;
use std::fs;
use std::io::Write;

use rayon::prelude::*;

use nsrs::constants::{BDD_ALPHAA, BDD_BETAA, M_NUCLEON};
use nsrs::core::magnetic::{B_SURFACE_G, bdd_field_g};
use nsrs::core::observations::{interpolate_at_mass, stable_branch};
use nsrs::core::tov_solver::{generate_star_sequence, radial_profile};
use nsrs::HadronsMatter;

use crate::cli::{Args, create_dir, f64_list};
use crate::study_b::{Config, ERG_PER_MEV_FM3, Eos, Profile, Row, nlem_label, parse_nlem, progress, solve};

/// MeV/fm^3 de B^2/8pi para B = 1 G.
const MAXWELL_MEV_FM3_PER_G2: f64 = 1.0 / (8.0 * std::f64::consts::PI * ERG_PER_MEV_FM3);

/// Campo dos níveis de Landau e campo da energia magnética numa linha.
fn fields(profile: Profile, b0: f64, n: f64) -> (f64, f64) {
    let bdd = bdd_field_g(B_SURFACE_G, b0, BDD_BETAA, BDD_ALPHAA, n);
    match profile {
        Profile::Constant => (b0, bdd),
        Profile::Bdd => (bdd, bdd),
    }
}

/// Curvatura da matéria c_m = 4 pi d^2 P_m/dB^2 (gaussiano) no mu_n da linha
/// e no campo `b`, com passos relativos `delta` e `delta/2`.
fn matter_curvature(model: nsrs::core::model::ModelParams, row: &Row, b: f64, delta: f64, config: &Config) -> Option<(f64, f64)> {
    let guess = [row[18], row[13] / M_NUCLEON, row[14] / M_NUCLEON, row[15] / M_NUCLEON, 0.0];
    let slope = |factor: f64| -> Option<f64> {
        let mut e = HadronsMatter::new(model, b * factor)
            .with_hyperons(config.hyperons)
            .with_anomalous_moments(config.amm)
            .with_limits(0.02, config.mu_max)
            .with_points(config.points)
            .with_max_landau_limit(config.landau_max);
        e.solve_point(row[17], &guess)?;
        Some(e.magnetization_b * ERG_PER_MEV_FM3 / (b * factor))
    };
    let c = |step: f64| -> Option<f64> {
        let (plus, minus) = (slope(1.0 + step)?, slope(1.0 - step)?);
        Some(4.0 * std::f64::consts::PI * (plus - minus) / (2.0 * step * b))
    };
    Some((c(delta)?, c(0.5 * delta)?))
}

/// n_B/n0 em função de P no núcleo (ramo em que P cresce com n_B).
fn density_of_pressure(rows: &[Row]) -> impl Fn(f64) -> f64 {
    let mut core: Vec<(f64, f64)> = rows.iter().filter(|r| r[0] > 0.0).map(|r| (r[2], r[0])).collect();
    core.sort_by(|a, b| a.1.total_cmp(&b.1));
    let mut table: Vec<(f64, f64)> = Vec::with_capacity(core.len());
    for (p, n) in core {
        if table.last().is_none_or(|&(lp, _)| p > lp) {
            table.push((p, n));
        }
    }
    move |p: f64| {
        let i = table.partition_point(|&(tp, _)| tp < p);
        if i == 0 || i >= table.len() {
            return f64::NAN;
        }
        let ((p0, n0), (p1, n1)) = (table[i - 1], table[i]);
        n0 + (n1 - n0) * (p - p0) / (p1 - p0)
    }
}

/// Pressões centrais da estrela de massa máxima e da de 1.4 M_sun.
fn chosen_stars(rows: &[Row]) -> Vec<(&'static str, f64, f64)> {
    let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
    let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    if rows.iter().map(|r| r[0]).fold(0.0, f64::max) < 1.0 {
        return Vec::new();
    }
    let sequence = generate_star_sequence(&eps, &p, &n, true);
    let branch = stable_branch(&sequence);
    let mut out = Vec::new();
    if let Some(max) = branch.last() {
        out.push(("max", max.central_pressure, max.mass));
    }
    if let Some(pc) = interpolate_at_mass(branch, 1.4, |s| s.central_pressure) {
        out.push(("1.4", pc, 1.4));
    }
    out
}

/// `study profile [--models GM1,GM3,FSU2] [--b0 1e17,1e18]
///  [--nlem maxwell,log:1e16,log:1e17,log:1e18] [--profiles constante,bdd]
///  [--points 1500] [--delta 1e-4] [--no-hyperons] [--amm] [--out results/study_profile]`
pub fn run(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &["no-hyperons", "amm"])?;
    let delta = args.f64_or("delta", 1e-4)?;
    if !(1e-6..=1e-2).contains(&delta) {
        return Err("use 1e-6 <= --delta <= 1e-2".into());
    }
    let config = Config {
        points: args.usize_or("points", 1500)?,
        mu_max: args.f64_or("mu-max", 3.0)?,
        hyperons: !args.switch("no-hyperons"),
        landau_max: args.usize_or("landau-max", nsrs::constants::MAX_LANDAU_LIMIT)?,
        amm: args.switch("amm"),
    };
    let b0s = f64_list(
        &args.value("b0").unwrap_or("1e17,1e18").split(',').map(|s| s.trim().to_string()).collect::<Vec<_>>(),
        "--b0",
    )?;
    if b0s.iter().any(|&b| b <= 0.0) {
        return Err("--b0: os campos devem ser > 0".into());
    }
    let nlems = parse_nlem(args.value("nlem").unwrap_or("maxwell,log:1e16,log:1e17,log:1e18"))?;
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
    let models = args.models(&["GM1", "GM3", "FSU2"])?;
    let out = args.value("out").unwrap_or("results/study_profile").to_string();
    let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads()?).build_global();
    create_dir(&out)?;

    let mut jobs = Vec::new();
    for (name, model) in &models {
        for &b0 in &b0s {
            for &profile in &profiles {
                for &nlem in &nlems {
                    jobs.push((name.clone(), *model, b0, profile, nlem));
                }
            }
        }
    }
    println!(
        "study profile: {} EoS ({} modelo(s), B0 = {:?} G, perfis {}, NLEM {}), hyperons = {}, amm = {}",
        jobs.len(),
        models.len(),
        b0s,
        profiles.iter().map(|p| p.label()).collect::<Vec<_>>().join(","),
        nlems.iter().map(|&n| nlem_label(n)).collect::<Vec<_>>().join(","),
        config.hyperons,
        config.amm
    );

    // 1. EoS de cada caso.
    println!("\n[1/2] EoS...");
    let bar = progress(jobs.len());
    let eos: Vec<Eos> = jobs
        .par_iter()
        .with_max_len(1)
        .map(|(_, model, b0, profile, nlem)| {
            let e = solve(*model, *b0, *profile, *nlem, &config);
            bar.inc(1);
            e
        })
        .collect();
    bar.finish_and_clear();

    // 2. Curvatura da matéria: não depende da NLEM (acoplamento mínimo), uma
    // vez por (modelo, B0, perfil), nas linhas da EoS de Maxwell (ou da
    // primeira NLEM). As demais NLEMs usam a linha de mesmo mu_n mais próxima.
    let reference: Vec<usize> = (0..jobs.len()).step_by(nlems.len()).collect();
    let eos_ref = &eos;
    let tasks: Vec<(usize, usize)> = reference
        .iter()
        .flat_map(|&j| (1..eos_ref[j].rows.len().saturating_sub(1)).filter(move |&i| eos_ref[j].rows[i][0] >= 0.05).map(move |i| (j, i)))
        .collect();
    println!("[2/2] curvatura da matéria: {} pontos, 4 soluções cada (delta = {delta:e} e delta/2)...", tasks.len());
    let bar = progress(tasks.len());
    let curvatures: HashMap<(usize, usize), (f64, f64)> = tasks
        .par_iter()
        .with_max_len(1)
        .map(|&(j, i)| {
            let (_, model, b0, profile, _) = &jobs[j];
            let row = &eos[j].rows[i];
            let (b_landau, _) = fields(*profile, *b0, row[0]);
            let c = matter_curvature(*model, row, b_landau, delta, &config);
            bar.inc(1);
            c.map(|c| ((j, i), c))
        })
        .flatten()
        .collect();
    bar.finish_and_clear();

    let eos_path = format!("{out}/eos.csv");
    let mut eos_csv = fs::File::create(&eos_path).map_err(|e| format!("{eos_path}: {e}"))?;
    writeln!(
        eos_csv,
        "model,nlem,profile,B0_G,n_over_n0,mu_n,B_landau_G,B_field_G,eps_MeV_fm3,p_perp_MeV_fm3,eps_field_MeV_fm3,p_par_field_MeV_fm3,p_perp_field_MeV_fm3,MB_MeV_fm3,four_pi_M_over_B,h_over_b_vac,permeability,f_vac,c_matter,c_matter_half_step,s,s_half_step,magnetically_unstable"
    )
    .map_err(|e| e.to_string())?;
    let radial_path = format!("{out}/radial.csv");
    let mut radial_csv = fs::File::create(&radial_path).map_err(|e| format!("{radial_path}: {e}"))?;
    writeln!(radial_csv, "model,nlem,profile,B0_G,star,M_Msun,R_km,r_km,m_Msun,P_MeV_fm3,n_over_n0").map_err(|e| e.to_string())?;

    for (j, (name, _, b0, profile, nlem)) in jobs.iter().enumerate() {
        let label = nlem_label(*nlem);
        let reference_job = j - j % nlems.len();
        let reference_rows = &eos[reference_job].rows;
        let curvature_at = |mu: f64| -> Option<(f64, f64)> {
            let i = (0..reference_rows.len()).min_by(|&a, &b| (reference_rows[a][17] - mu).abs().total_cmp(&(reference_rows[b][17] - mu).abs()))?;
            if (reference_rows[i][17] - mu).abs() > 1e-9 * mu.abs().max(1.0) {
                return None;
            }
            curvatures.get(&(reference_job, i)).copied()
        };
        let e = &eos[j];
        for (i, row) in e.rows.iter().enumerate() {
            if row[0] < 0.05 {
                continue;
            }
            let (b_landau, b_field) = fields(*profile, *b0, row[0]);
            let mb = e.magnetization_b[i];
            let stress = e.field_stress[i];
            let four_pi_m_over_b = mb / (2.0 * MAXWELL_MEV_FM3_PER_G2 * b_landau * b_landau);
            let h_over_b = nlem.h_over_b(b_landau);
            let permeability = 1.0 / (h_over_b - four_pi_m_over_b);
            let f_vac = nlem.vacuum_curvature(b_landau);
            let (c, c_half) = curvature_at(row[17]).unwrap_or((f64::NAN, f64::NAN));
            let (s, s_half) = (f_vac - c, f_vac - c_half);
            writeln!(
                eos_csv,
                "{name},{label},{},{b0:.4e},{:.6},{:.8},{b_landau:.6e},{b_field:.6e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{mb:.8e},{four_pi_m_over_b:.8e},{h_over_b:.8e},{permeability:.8e},{f_vac:.8e},{c:.8e},{c_half:.8e},{s:.8e},{s_half:.8e},{}",
                profile.label(),
                row[0],
                row[17],
                row[1],
                row[2],
                stress.energy,
                stress.p_parallel,
                stress.p_perpendicular,
                s < 0.0 && s_half < 0.0
            )
            .map_err(|e| e.to_string())?;
        }

        let eps: Vec<f64> = e.rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = e.rows.iter().map(|r| r[2]).collect();
        let n: Vec<f64> = e.rows.iter().map(|r| r[0]).collect();
        let density = density_of_pressure(&e.rows);
        for (star, pc, mass) in chosen_stars(&e.rows) {
            let points = radial_profile(&eps, &p, &n, pc);
            let Some(surface) = points.last() else { continue };
            for point in &points {
                writeln!(
                    radial_csv,
                    "{name},{label},{},{b0:.4e},{star},{mass:.6},{:.6},{:.6},{:.6e},{:.8e},{:.6}",
                    profile.label(),
                    surface[0],
                    point[0],
                    point[2],
                    point[1],
                    density(point[1])
                )
                .map_err(|e| e.to_string())?;
            }
        }
    }
    println!("\nResultados em {eos_path} e {radial_path}");
    println!("Figuras: python plot_scripts/magnetic_profiles.py {out}");
    Ok(())
}
