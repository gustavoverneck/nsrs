// src/bin/darkphotons.rs
//
// Four literature-motivated benchmark scenarios.
//
// H0 : pure hadronic matter
// S1 : Kumar et al. Set 1 inspired
// S2 : Kumar et al. Set 2 inspired
// S3 : Kumar et al. Set 3 inspired
//
// IMPORTANT:
// - GM1 and GM3 are both evaluated for the SAME four physical scenarios.
// - FSU2 remains excluded.
// - Microscopic B = 0.
// - epsilon is NOT taken from Kumar et al.; their model is a direct Z' portal.
//   Here epsilon = 1e-4 is a fixed kinetic-mixing benchmark.
// - Y_chi is chosen so that kF_chi = 20 MeV at n0 = 0.153 fm^-3.

use nsrs::{DarkPhotonsMatter, EngineMode, GM1, GM3, HadronsMatter, Solver};

use std::fs;
use std::io::{self, Write};
use std::path::Path;

// ============================================================================
// Global numerical setup
// ============================================================================

const B_FIELD_G: f64 = 0.0;

const MU_N_MIN: f64 = 0.01;
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

fn main() -> io::Result<()> {
    // ------------------------------------------------------------------------
    // Dark fraction corresponding to kF_chi = 20 MeV at n0 = 0.153 fm^-3.
    // ------------------------------------------------------------------------

    let y_chi_benchmark =
        y_chi_from_reference_kf(KF_CHI_REF_MEV, N0_REF_FM3);

    println!(
        "Reference dark fraction: Y_chi = {:.8e}",
        y_chi_benchmark
    );

    println!(
        "By construction: kF_chi(n0={:.3} fm^-3) = {:.1} MeV",
        N0_REF_FM3,
        KF_CHI_REF_MEV
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

    let models = [
        ("GM1", GM1),
        ("GM3", GM3),
    ];

    for (model_name, model) in models {
        let base_dir =
            format!("output/darkphotons_benchmarks/{}", model_name);

        fs::create_dir_all(&base_dir)?;

        let summary_tmp =
            format!("{}/summary.csv.tmp", base_dir);

        let summary_final =
            format!("{}/summary.csv", base_dir);

        let mut summary =
            fs::File::create(&summary_tmp)?;

        writeln!(
            summary,
            concat!(
                "scenario,description,model,",
                "epsilon,m_chi_MeV,m_x_MeV,g_d,y_chi,",
                "kf_chi_ref_MeV,nB_ref_fm-3,",
                "b_field_G,eos_file,status"
            )
        )?;

        // ====================================================================
        // H0 -- purely hadronic baseline
        // ====================================================================

        let h0_filename = "eos_H0_hadrons.dat";
        let h0_path =
            format!("{}/{}", base_dir, h0_filename);

        println!(
            "\n[{}] H0: pure hadronic matter",
            model_name
        );

        let hadrons_motor =
            HadronsMatter::new(model, B_FIELD_G)
                .with_limits(MU_N_MIN, MU_N_MAX)
                .with_points(EOS_POINTS)
                .with_eos_output(&h0_path);

        let mut hadrons_solver =
            Solver::new(EngineMode::Hadrons(hadrons_motor));

        let _ = hadrons_solver.solve();

        let h0_status =
            eos_file_status(&h0_path);

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
        )?;

        // ====================================================================
        // S1, S2, S3
        // ====================================================================

        for scenario in dark_scenarios {
            let eos_filename =
                format!("eos_{}_{}.dat",
                    scenario.label,
                    scenario.description
                );

            let eos_path =
                format!("{}/{}", base_dir, eos_filename);

            println!(
                "[{}] {}: {}",
                model_name,
                scenario.label,
                scenario.description
            );

            println!(
                "    m_chi = {:.6e} MeV",
                scenario.m_chi_mev
            );

            println!(
                "    m_X   = {:.6e} MeV",
                scenario.m_x_mev
            );

            println!(
                "    g_D   = {:.6e}",
                scenario.g_d
            );

            println!(
                "    eps   = {:.6e}",
                scenario.epsilon
            );

            println!(
                "    Y_chi = {:.6e}",
                scenario.y_chi
            );

            let dark_motor =
                DarkPhotonsMatter::new(model, B_FIELD_G)
                    .with_limits(MU_N_MIN, MU_N_MAX)
                    .with_points(EOS_POINTS)
                    .with_epsilon(scenario.epsilon)
                    .with_m_x_mev(scenario.m_x_mev)
                    .with_m_chi_mev(scenario.m_chi_mev)
                    .with_g_d(scenario.g_d)
                    .with_y_chi(scenario.y_chi)
                    .with_eos_output(&eos_path);

            let mut solver =
                Solver::new(
                    EngineMode::DarkPhotons(dark_motor)
                );

            let _ = solver.solve();

            // The row is written only after the solver has returned
            // and the output file has been checked.
            let status =
                eos_file_status(&eos_path);

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
            )?;
        }

        summary.flush()?;
        drop(summary);

        // Transactional finalization:
        //
        // summary.csv is only produced once all four scenarios
        // for the current hadronic model have returned.
        fs::rename(
            &summary_tmp,
            &summary_final,
        )?;

        println!(
            "[{}] summary committed: {}",
            model_name,
            summary_final
        );
    }

    println!(
        "\nConcluido. Resultados em \
         output/darkphotons_benchmarks/"
    );

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
fn y_chi_from_reference_kf(
    kf_mev: f64,
    n_b_ref_fm3: f64,
) -> f64 {
    let kf_fm_inv =
        kf_mev / HBARC_MEV_FM;

    let n_chi_fm3 =
        kf_fm_inv.powi(3)
            / (3.0 * std::f64::consts::PI.powi(2));

    n_chi_fm3 / n_b_ref_fm3
}

// ============================================================================
// Minimal file validation
// ============================================================================

fn eos_file_status(path: &str) -> &'static str {
    let p = Path::new(path);

    if !p.exists() {
        return "FAILED_MISSING";
    }

    match fs::metadata(p) {
        Ok(metadata) if metadata.len() > 0 => "EOS_WRITTEN",
        Ok(_) => "FAILED_EMPTY",
        Err(_) => "FAILED_METADATA",
    }
}