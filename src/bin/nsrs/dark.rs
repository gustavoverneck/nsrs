// Setor escuro (substitui os binários darkphotons, darkphotons_base e
// single_darkphotons). Caminhos de saída iguais aos originais.

use std::fs;
use std::io::{self, Write};

use nsrs::constants::{DATA_SIZE, M_NUCLEON, RESULTS_SIZE};
use nsrs::{DarkPhotonsMatter, EngineMode, GM1, HadronsMatter, Solver};

use crate::cli::{Args, create_dir, format_sci};

/// `dark scan`: grade de 10 x 10 x 10 x 10 valores de (epsilon, m_X/M_N,
/// g_D, Y_chi) com m_chi = M_N, mais uma EoS hadrônica de referência.
/// Saída: output/darkphotons_scan/<modelo>/{summary.csv, eos_*.dat}
pub fn scan(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let b_field = args.f64_or("b", 1e17)?;
    let m_chi_mev = M_NUCLEON;
    let eps_values = linspace(1e-6, 1e-3, 10);
    let m_x_values = linspace(0.01, 0.11, 10);
    let g_d_values = linspace(0.35, 1.12, 10);
    let y_chi_values = linspace(0.0, 0.1, 10);

    for (model_name, model) in args.models(&["GM1", "GM3"])? {
        let base_dir = format!("output/darkphotons_scan/{model_name}");
        create_dir(&base_dir)?;
        let mut summary =
            fs::File::create(format!("{base_dir}/summary.csv")).map_err(|e| e.to_string())?;
        let mut line = |text: String| writeln!(summary, "{text}").map_err(|e| e.to_string());
        line("label,epsilon,m_x_over_mN,g_d,y_chi,m_chi_MeV,b_field_G,eos_file".into())?;

        let mut engines = vec![EngineMode::Hadrons(
            HadronsMatter::new(model, b_field)
                .with_limits(0.01, 2.0)
                .with_points(1200)
                .with_eos_output(format!("{base_dir}/eos_hadrons.dat")),
        )];
        line(format!("hadrons,,,,,,{b_field:.6e},eos_hadrons.dat"))?;

        for &epsilon in &eps_values {
            for &m_x in &m_x_values {
                for &g_d in &g_d_values {
                    for &y_chi in &y_chi_values {
                        let file = format!(
                            "eos_eps_{}_mx_{}_gd_{:.3}_ychi_{}.dat",
                            format_sci(epsilon),
                            format_sci(m_x),
                            g_d,
                            format_sci(y_chi)
                        );
                        engines.push(EngineMode::DarkPhotons(
                            DarkPhotonsMatter::new(model, b_field)
                                .with_limits(0.01, 2.0)
                                .with_points(1200)
                                .with_epsilon(epsilon)
                                .with_m_x(m_x)
                                .with_m_chi_mev(m_chi_mev)
                                .with_g_d(g_d)
                                .with_y_chi(y_chi)
                                .with_eos_output(format!("{base_dir}/{file}")),
                        ));
                        line(format!(
                            "darkphotons,{epsilon:.6e},{m_x:.6e},{g_d:.6e},{y_chi:.6e},{m_chi_mev:.6e},{b_field:.6e},{file}"
                        ))?;
                    }
                }
            }
        }
        println!("\nModelo={model_name} | {} EoS...", engines.len());
        Solver::solve_parallel(engines, args.threads()?);
    }
    println!("\nConcluído. Dados em output/darkphotons_scan/");
    Ok(())
}

/// `dark single`: GM1 com B = 1e17 G, uma EoS hadrônica e uma com o setor
/// escuro (epsilon = 1e-4, m_X = 0.1065 M_N, m_chi = M_N, g_D = 0.45,
/// Y_chi = 0.01).
/// Saída: output/darkphotons/GM1/{eos_hadrons.dat, darkphoton.dat}
pub fn single(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let b_field = 1e17;
    let base_dir = "output/darkphotons/GM1";
    create_dir(base_dir)?;
    let engines = vec![
        EngineMode::Hadrons(
            HadronsMatter::new(GM1, b_field)
                .with_limits(0.01, 2.0)
                .with_points(2000)
                .with_eos_output(format!("{base_dir}/eos_hadrons.dat")),
        ),
        EngineMode::DarkPhotons(
            DarkPhotonsMatter::new(GM1, b_field)
                .with_points(2000)
                .with_limits(0.01, 2.0)
                .with_epsilon(1e-4)
                .with_m_x(0.1065)
                .with_m_chi_mev(M_NUCLEON)
                .with_g_d(0.45)
                .with_y_chi(0.01)
                .with_eos_output(format!("{base_dir}/darkphoton.dat")),
        ),
    ];
    Solver::solve_parallel(engines, args.threads()?);
    println!("\nConcluído. Dados em {base_dir}/");
    Ok(())
}

// ============================================================================
// Benchmarks (antigo darkphotons_base)
//
// Four literature-motivated benchmark scenarios.
//
// H0 : pure hadronic matter
// S1 : Kumar et al. Set 1 inspired
// S2 : Kumar et al. Set 2 inspired
// S3 : Kumar et al. Set 3 inspired
//
// - GM1 and GM3 (default) are evaluated for the SAME four scenarios.
// - Microscopic B = 0.
// - epsilon is NOT taken from Kumar et al.; their model is a direct Z' portal.
//   Here epsilon = 1e-4 is a fixed kinetic-mixing benchmark.
// - Y_chi is chosen so that kF_chi = 20 MeV at n0 = 0.153 fm^-3.
// ============================================================================


const B_FIELD_G: f64 = 0.0;

const MU_N_MIN: f64 = 1.00;
const MU_N_MAX: f64 = 2.00;
const EOS_POINTS: usize = 2001;

// Reference saturation density used only to map the literature kF_chi
// benchmark into the Y_chi prescription of the present model.
const N0_REF_FM3: f64 = 0.153;

// Literature benchmark.
const KF_CHI_REF_MEV: f64 = 20.0;

// hbar*c in MeV fm.
const HBARC_MEV_FM: f64 = 197.326_980_4;

// Kinetic-mixing benchmark.
// This is NOT the g_q coupling of Kumar et al.
const EPSILON_BENCHMARK: f64 = 1.0e-4;

// ============================================================================
// Scenario definition
// ============================================================================

#[derive(Clone, Copy)]
struct DarkScenario {
    label: &'static str,
    description: &'static str,

    m_chi_mev: f64,
    m_x_mev: f64,
    g_d: f64,
    epsilon: f64,
    y_chi: f64,
}

/// `dark benchmarks`: cenários H0, S1, S2, S3.
/// Saída: output/darkphotons_benchmarks/<modelo>/{summary.csv, eos_*.dat}
pub fn benchmarks(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let mut all_outputs_complete = true;
    // ------------------------------------------------------------------------
    // Dark fraction corresponding to kF_chi = 20 MeV at n0 = 0.153 fm^-3.
    // ------------------------------------------------------------------------

    let y_chi_benchmark = y_chi_from_reference_kf(KF_CHI_REF_MEV, N0_REF_FM3);

    println!("Reference dark fraction: Y_chi = {:.8e}", y_chi_benchmark);

    println!(
        "By construction: kF_chi(n0={:.3} fm^-3) = {:.1} MeV",
        N0_REF_FM3, KF_CHI_REF_MEV
    );

    // ------------------------------------------------------------------------
    // Literature-inspired dark scenarios.
    //
    // Kumar et al.:
    //
    // Set 1: M_chi = 200 GeV,  m_Z' = 1800 GeV, g_chi = 0.45
    // Set 2: M_chi = 1800 GeV, m_Z' =  900 GeV, g_chi = 0.25
    // Set 3: M_chi = 200 GeV,  m_Z' =  100 MeV, g_chi = 0.45
    //
    // We identify their g_chi only with our DARK coupling g_D.
    //
    // Their visible coupling g_q is NOT identified with epsilon.
    // ------------------------------------------------------------------------

    let dark_scenarios = [
        DarkScenario {
            label: "S1",
            description: "heavy_mediator_200GeV_DM",
            m_chi_mev: 200_000.0,
            m_x_mev: 1_800_000.0,
            g_d: 0.45,
            epsilon: EPSILON_BENCHMARK,
            y_chi: y_chi_benchmark,
        },
        DarkScenario {
            label: "S2",
            description: "heavy_DM_heavy_mediator",
            m_chi_mev: 1_800_000.0,
            m_x_mev: 900_000.0,
            g_d: 0.25,
            epsilon: EPSILON_BENCHMARK,
            y_chi: y_chi_benchmark,
        },
        DarkScenario {
            label: "S3",
            description: "light_mediator_200GeV_DM",
            m_chi_mev: 200_000.0,
            m_x_mev: 100.0,
            g_d: 0.45,
            epsilon: EPSILON_BENCHMARK,
            y_chi: y_chi_benchmark,
        },
    ];

    // ------------------------------------------------------------------------
    // GM1 and GM3 are not additional physical scenarios.
    //
    // They represent two hadronic descriptions under which the same four
    // scenarios H0/S1/S2/S3 are evaluated.
    // ------------------------------------------------------------------------

    let models = args.models(&["GM1", "GM3"])?;

    for (model_name, model) in models {
        let model_name = model_name.as_str();
        let base_dir = format!("output/darkphotons_benchmarks/{}", model_name);

        fs::create_dir_all(&base_dir).map_err(|e| e.to_string())?;

        let summary_tmp = format!("{}/summary.csv.tmp", base_dir);

        let summary_final = format!("{}/summary.csv", base_dir);

        let mut summary = fs::File::create(&summary_tmp).map_err(|e| e.to_string())?;

        writeln!(
            summary,
            concat!(
                "scenario,description,model,",
                "epsilon,m_chi_MeV,m_x_MeV,g_d,y_chi,",
                "kf_chi_ref_MeV,nB_ref_fm-3,",
                "b_field_G,eos_file,status"
            )
        ).map_err(|e| e.to_string())?;

        // ====================================================================
        // H0 -- purely hadronic baseline
        // ====================================================================

        let h0_filename = "eos_H0_hadrons.dat";
        let h0_path = format!("{}/{}", base_dir, h0_filename);

        println!("\n[{}] H0: pure hadronic matter", model_name);

        let hadrons_motor = HadronsMatter::new(model, B_FIELD_G)
            .with_limits(MU_N_MIN, MU_N_MAX)
            .with_points(EOS_POINTS)
            .with_eos_output(&h0_path);

        let mut hadrons_solver = Solver::new(EngineMode::Hadrons(hadrons_motor));

        let _ = hadrons_solver.solve();

        let h0_status = eos_file_status(&h0_path);
        all_outputs_complete &= h0_status == "EOS_WRITTEN";

        writeln!(
            summary,
            concat!(
                "H0,pure_hadronic,{model},",
                "0,0,0,0,0,",
                "0,{n0:.6e},",
                "{b:.6e},{file},{status}"
            ),
            model = model_name,
            n0 = N0_REF_FM3,
            b = B_FIELD_G,
            file = h0_filename,
            status = h0_status,
        ).map_err(|e| e.to_string())?;

        // ====================================================================
        // S1, S2, S3
        // ====================================================================

        for scenario in dark_scenarios {
            let eos_filename = format!("eos_{}_{}.dat", scenario.label, scenario.description);

            let eos_path = format!("{}/{}", base_dir, eos_filename);

            println!(
                "[{}] {}: {}",
                model_name, scenario.label, scenario.description
            );

            println!("    m_chi = {:.6e} MeV", scenario.m_chi_mev);

            println!("    m_X   = {:.6e} MeV", scenario.m_x_mev);

            println!("    g_D   = {:.6e}", scenario.g_d);

            println!("    eps   = {:.6e}", scenario.epsilon);

            println!("    Y_chi = {:.6e}", scenario.y_chi);

            let dark_motor = DarkPhotonsMatter::new(model, B_FIELD_G)
                .with_limits(MU_N_MIN, MU_N_MAX)
                .with_points(EOS_POINTS)
                .with_epsilon(scenario.epsilon)
                .with_m_x_mev(scenario.m_x_mev)
                .with_m_chi_mev(scenario.m_chi_mev)
                .with_g_d(scenario.g_d)
                .with_y_chi(scenario.y_chi)
                .with_eos_output(&eos_path);

            let mut solver = Solver::new(EngineMode::DarkPhotons(dark_motor));

            let _ = solver.solve();

            // The row is written only after the solver has returned
            // and the output file has been checked.
            let status = eos_file_status(&eos_path);
            all_outputs_complete &= status == "EOS_WRITTEN";

            writeln!(
                summary,
                concat!(
                    "{label},{description},{model},",
                    "{epsilon:.8e},",
                    "{mchi:.8e},",
                    "{mx:.8e},",
                    "{gd:.8e},",
                    "{ychi:.8e},",
                    "{kf:.8e},",
                    "{n0:.8e},",
                    "{b:.8e},",
                    "{file},{status}"
                ),
                label = scenario.label,
                description = scenario.description,
                model = model_name,
                epsilon = scenario.epsilon,
                mchi = scenario.m_chi_mev,
                mx = scenario.m_x_mev,
                gd = scenario.g_d,
                ychi = scenario.y_chi,
                kf = KF_CHI_REF_MEV,
                n0 = N0_REF_FM3,
                b = B_FIELD_G,
                file = eos_filename,
                status = status,
            ).map_err(|e| e.to_string())?;
        }

        summary.flush().map_err(|e| e.to_string())?;
        drop(summary);

        // Transactional finalization:
        //
        // summary.csv is only produced once all four scenarios
        // for the current hadronic model have returned.
        fs::rename(&summary_tmp, &summary_final).map_err(|e| e.to_string())?;

        println!("[{}] summary committed: {}", model_name, summary_final);
    }

    if !all_outputs_complete {
        return Err("one or more benchmark EOS files are incomplete; inspect summary.csv".into());
    }

    println!("\nConcluido. Resultados em output/darkphotons_benchmarks/");

    Ok(())
}

// ============================================================================
// Physical conversion
// ============================================================================

/// Returns Y_chi such that the chosen Fermi momentum is reproduced at
/// the reference baryon density:
///
///     n_chi = kF_chi^3 / (3 pi^2)
///     Y_chi = n_chi / n_B
///
/// kF is supplied in MeV, n_B in fm^-3.
fn y_chi_from_reference_kf(kf_mev: f64, n_b_ref_fm3: f64) -> f64 {
    let kf_fm_inv = kf_mev / HBARC_MEV_FM;

    let n_chi_fm3 = kf_fm_inv.powi(3) / (3.0 * std::f64::consts::PI.powi(2));

    n_chi_fm3 / n_b_ref_fm3
}

// ============================================================================
// Minimal file validation
// ============================================================================

fn eos_file_status(path: &str) -> &'static str {
    let content = match fs::read_to_string(path) {
        Ok(content) => content,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return "FAILED_MISSING",
        Err(_) => return "FAILED_READ",
    };

    let mut eos_rows = 0usize;
    let mut valid_mr_rows = 0usize;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let values: Vec<f64> = match line
            .split_whitespace()
            .map(str::parse::<f64>)
            .collect::<Result<_, _>>()
        {
            Ok(values) => values,
            Err(_) => return "FAILED_FORMAT",
        };

        if values.len() != DATA_SIZE
            || values[..RESULTS_SIZE]
                .iter()
                .any(|value| !value.is_finite())
        {
            return "FAILED_FORMAT";
        }

        eos_rows += 1;
        if values[RESULTS_SIZE..]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        {
            valid_mr_rows += 1;
        }
    }

    if eos_rows == 0 {
        "FAILED_NO_EOS_ROWS"
    } else if valid_mr_rows < 3 {
        "FAILED_NO_MR_CURVE"
    } else {
        "EOS_WRITTEN"
    }
}

fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    match n {
        0 => Vec::new(),
        1 => vec![start],
        _ => (0..n).map(|i| start + (end - start) * i as f64 / (n - 1) as f64).collect(),
    }
}
