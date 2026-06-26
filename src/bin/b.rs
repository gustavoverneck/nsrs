// src/bin/b.rs

use nsrs::{
    EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver, constants::BCE_G, generate_mr_curve,
    write_eos_with_mr,
};
use std::fs;

fn build_b_fields(total_points: usize, b_min: f64, b_max: f64) -> Vec<f64> {
    assert!(total_points >= 2, "total_points deve ser >= 2");
    assert!(b_min > 0.0, "b_min deve ser > 0");
    assert!(b_max > b_min, "b_max deve ser > b_min");

    // Reservamos 1 ponto para B = 0 e distribuímos os demais em escala log
    // de b_min (campo crítico) até b_max.
    let non_zero_points = total_points - 1;
    let log_min = b_min.log10();
    let log_max = b_max.log10();

    let mut b_vals = Vec::with_capacity(total_points);
    b_vals.push(0.0);

    if non_zero_points == 1 {
        b_vals.push(b_min);
        return b_vals;
    }

    for i in 0..non_zero_points {
        let t = i as f64 / (non_zero_points - 1) as f64;
        let log_b = log_min + t * (log_max - log_min);
        b_vals.push(10_f64.powf(log_b));
    }

    b_vals
}

fn main() {
    let total_points = 100;
    let b_min = BCE_G; // campo crítico em Gauss
    let b_max = 3.0e18;

    let b_fields = build_b_fields(total_points, b_min, b_max);
    let models = [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)];

    for (model_name, model_params) in models {
        let engines: Vec<EngineMode> = b_fields
            .iter()
            .map(|&b_field| {
                let motor = HadronsMatter::new(model_params, b_field)
                    .with_limits(0.01, 2.0)
                    .with_points(1201);
                EngineMode::Hadrons(motor)
            })
            .collect();

        println!(
            "\nModelo={} | Topologia=padrão | Varrendo {} valores de B (0, Bc..3e18 G)...",
            model_name, total_points
        );
        let all_results = Solver::solve_parallel(engines, 16);

        for (i, results) in all_results.iter().enumerate() {
            let b_field = b_fields[i];
            let b_string = format!("{:.2e}", b_field);

            let dir_path = format!("output/b/{}/B_{}/default", model_name, b_string);
            if fs::create_dir_all(&dir_path).is_err() {
                continue;
            }

            let eps: Vec<f64> = results.iter().map(|r| r[1]).collect();
            let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();
            let rho_arr: Vec<f64> = results.iter().map(|r| r[0]).collect();
            let (masses, radii, b_masses, central_p_list) =
                generate_mr_curve(&eps, &p_arr, &rho_arr, false);

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
        }
    }

    println!("\nProcesso concluído! Dados salvos em output/b/");
}
