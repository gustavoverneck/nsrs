// src/bin/darkphotons.rs

use nsrs::constants::M_NUCLEON;
use nsrs::{DarkPhotonsMatter, EngineMode, GM1, HadronsMatter, Solver};

fn main() {
    let b_field = 1e17;
    let base_dir = "output/darkphotons/GM1";
    let eos_path = format!("{}/darkphoton.dat", base_dir);
    let mut engines = Vec::new();

    let hadrons_filename = "eos_hadrons.dat";
    let hadrons_path = format!("{}/{}", base_dir, hadrons_filename);
    let hadrons_motor = HadronsMatter::new(GM1, b_field)
        .with_limits(0.01, 2.0)
        .with_points(2000)
        .with_eos_output(&hadrons_path);
    // let mut hadrons_solver = Solver::new(EngineMode::Hadrons(hadrons_motor));

    engines.push(EngineMode::Hadrons(hadrons_motor));

    let dark_motor = DarkPhotonsMatter::new(GM1, b_field)
        .with_points(2000)
        .with_limits(0.01, 2.0)
        .with_epsilon(1e-4)
        .with_m_x(0.1065)
        .with_m_chi_mev(M_NUCLEON)
        .with_g_d(0.45)
        .with_y_chi(0.01)
        .with_eos_output(&eos_path);

    engines.push(EngineMode::DarkPhotons(dark_motor));

    Solver::solve_parallel(engines, 2);
    println!("\nProcesso concluido!");
}
