// src/bin/b.rs

use nsrs::{
    generate_mr_curve, write_eos_with_mr, Artist, DarkPhotonsMatter, EngineMode, GM1, GM3,
    Solver,
};
use std::fs;

fn main() {
	let models = [("GM1", GM1), ("GM3", GM3)];

    let mut engines: Vec<EngineMode> = Vec::with_capacity(models.len());
    let mut model_names: Vec<&str> = Vec::with_capacity(models.len());

    for (model_name, model) in models {
        let motor = DarkPhotonsMatter::new(model)
            .with_limits(0.01, 2.0)
            .with_points(1200);
        engines.push(EngineMode::DarkPhotons(motor));
        model_names.push(model_name);
    }

    let all_results = Solver::solve_parallel(engines, 2);

    for (idx, results) in all_results.iter().enumerate() {
        let model_name = model_names[idx];
        let dir_path = format!("output/darkphotons/{}", model_name);
        if fs::create_dir_all(&dir_path).is_err() {
            continue;
        }

        let eps: Vec<f64> = results.iter().map(|r| r[1]).collect();
        let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();
        let rho_arr: Vec<f64> = results.iter().map(|r| r[0]).collect();
        let (masses, radii, b_masses, central_p_list) = generate_mr_curve(&eps, &p_arr, &rho_arr, false);

        let eos_filename = format!("{}/eos.dat", dir_path);
        if write_eos_with_mr(
            results,
            &masses,
            &radii,
            &b_masses,
            &central_p_list,
            &eos_filename,
        )
        .is_err()
        {
            continue;
        }

        let eos_plot = format!("{}/eos.svg", dir_path);
        let mr_plot = format!("{}/mr.svg", dir_path);
        let be_vs_m_plot = format!("{}/be_vs_m.svg", dir_path);
        let m_vs_pc_plot = format!("{}/m_vs_pc.svg", dir_path);

        let _ = Artist::new(&eos_plot, "Equation of State")
            .with_x_label("Energy Density ε [MeV/fm³]")
            .with_y_label("Pressure P [MeV/fm³]")
            .add_curve(&eps, &p_arr, model_name)
            .plot();

        let _ = Artist::new(&mr_plot, "Mass-Radius Relation")
            .with_x_label("Radius [km]")
            .with_y_label("Mass [M⊙]")
            .add_curve(&radii, &masses, model_name)
            .autoscale()
            .plot();

        let n_ratio = masses.len().min(b_masses.len());
        let m_for_be: Vec<f64> = masses.iter().copied().take(n_ratio).collect();
        let binding_energy: Vec<f64> = b_masses
            .iter()
            .zip(masses.iter())
            .take(n_ratio)
            .map(|(&mb, &m)| mb - m)
            .collect();

        let _ = Artist::new(&be_vs_m_plot, "Gravitational Binding Energy")
            .with_x_label("Gravitational Mass M [M⊙]")
            .with_y_label("Binding Energy (Mb - M) [M⊙]")
            .add_curve(&m_for_be, &binding_energy, model_name)
            .plot();

        let n_pc = masses.len().min(central_p_list.len());
        let pc_filtered: Vec<f64> = central_p_list.iter().copied().take(n_pc).collect();
        let m_for_pc: Vec<f64> = masses.iter().copied().take(n_pc).collect();

        let _ = Artist::new(&m_vs_pc_plot, "Stability Curve")
            .with_x_label("Central Pressure Pc [MeV/fm³]")
            .with_y_label("Gravitational Mass [M⊙]")
            .add_curve(&pc_filtered, &m_for_pc, model_name)
            .plot();
    }

    println!("\nProcesso concluído! Dados salvos em output/darkphotons/");
}