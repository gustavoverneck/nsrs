// src/bin/darkphotons.rs

use nsrs::constants::M_NUCLEON;
use nsrs::{DarkPhotonsMatter, EngineMode, GM1, GM3, HadronsMatter, Solver};
use std::fs;
use std::io::Write;

fn main() {
    // FSU2 is intentionally excluded until the legacy nonlinear omega/rho
    // normalization documented in docs/PHYSICS.md is corrected globally.
    let models = [("GM1", GM1), ("GM3", GM3)];
    let b_field = 1e17;
    let m_chi_mev = M_NUCLEON;

    let eps_values = linspace(1e-6, 1e-3, 10);
    let m_x_values = linspace(0.01, 0.11, 10);
    let g_d_values = linspace(0.35, 1.12, 10);
    let y_chi_values = linspace(0.0, 0.1, 10);

    for (model_name, model) in models {
        let base_dir = format!("output/darkphotons_scan/{}", model_name);
        if fs::create_dir_all(&base_dir).is_err() {
            continue;
        }

        let summary_path = format!("{}/summary.csv", base_dir);
        let mut summary = match fs::File::create(&summary_path) {
            Ok(file) => file,
            Err(_) => continue,
        };
        let _ = writeln!(
            summary,
            "label,epsilon,m_x_over_mN,g_d,y_chi,m_chi_MeV,b_field_G,eos_file"
        );

        let hadrons_filename = "eos_hadrons.dat";
        let hadrons_path = format!("{}/{}", base_dir, hadrons_filename);
        let hadrons_motor = HadronsMatter::new(model, b_field)
            .with_limits(0.01, 2.0)
            .with_points(1200)
            .with_eos_output(&hadrons_path);
        let mut hadrons_solver = Solver::new(EngineMode::Hadrons(hadrons_motor));
        let _ = hadrons_solver.solve();
        let _ = writeln!(summary, "hadrons,,,,,,{:.6e},{}", b_field, hadrons_filename);

        let mut engines = Vec::new();
        for epsilon in &eps_values {
            for m_x in &m_x_values {
                for g_d in &g_d_values {
                    for y_chi in &y_chi_values {
                        let eos_filename = format!(
                            "eos_eps_{}_mx_{}_gd_{}_ychi_{}.dat",
                            format_sci(*epsilon),
                            format_sci(*m_x),
                            format_fixed(*g_d, 3),
                            format_sci(*y_chi)
                        );
                        let eos_path = format!("{}/{}", base_dir, eos_filename);
                        let dark_motor = DarkPhotonsMatter::new(model, b_field)
                            .with_limits(0.01, 2.0)
                            .with_points(1200)
                            .with_epsilon(*epsilon)
                            .with_m_x(*m_x)
                            .with_m_chi_mev(m_chi_mev)
                            .with_g_d(*g_d)
                            .with_y_chi(*y_chi)
                            .with_eos_output(&eos_path);
                        engines.push(EngineMode::DarkPhotons(dark_motor));
                        let _ = writeln!(
                            summary,
                            "darkphotons,{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{}",
                            epsilon, m_x, g_d, y_chi, m_chi_mev, b_field, eos_filename
                        );
                    }
                }
            }
        }

        let _ = Solver::solve_parallel(engines, 16);
    }

    println!("\nProcesso concluido! Dados salvos em output/darkphotons_scan/");
}

fn format_sci(value: f64) -> String {
    format!("{:.2e}", value).replace('+', "")
}

fn format_fixed(value: f64, decimals: usize) -> String {
    format!("{:.*}", decimals, value)
}

fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![start];
    }
    let step = (end - start) / (n as f64 - 1.0);
    (0..n).map(|i| start + step * i as f64).collect()
}
