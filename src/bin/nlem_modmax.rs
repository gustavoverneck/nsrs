// src/bin/nlem_modmax.rs

use nsrs::{
    EngineMode, GM1, GM3, HadronsMatter, MagneticTopology, NlemModel, Solver, generate_mr_curve,
    write_eos_with_mr,
};
use std::env;
use std::fs;

fn main() {
    let args: Vec<String> = env::args().collect();

    let csi_vals: Vec<f64> = (0..=10).map(|i| i as f64 * 0.05).collect();
    
    let b_fields: Vec<f64> = if args.len() > 1 {
        args[1..]
            .iter()
            .map(|s| s.parse().expect("B deve ser um número válido"))
            .collect()
    } else {
        eprintln!("Uso: {} <B1> <B2> ...", args[0]);
        eprintln!("Exemplo: {} 1e15 5e16", args[0]);
        return;
    };

    let num_points = csi_vals.len();
    let models = [("GM1", GM1), ("GM3", GM3)];
    let topologies = [
        // ("isotropic", MagneticTopology::Isotropic),
        ("anisotropic", MagneticTopology::Anisotropic),
    ];

    // 2. Loop principal: modelos x campos x topologias
    for (model_name, model_params) in models {
        for &b_field in &b_fields {
            let b_string = format!("{:.2e}", b_field);

            for (topology_name, topology_mode) in topologies {
                let engines: Vec<EngineMode> = csi_vals
                    .iter()
                    .map(|&csi| {
                        let motor = HadronsMatter::new(model_params, b_field)
                            .with_topology(topology_mode)
                            .with_nlem(NlemModel::Modmax(csi))
                            .with_limits(0.02, 2.0)
                            .with_points(1000);

                        EngineMode::Hadrons(motor)
                    })
                    .collect();

                println!(
                    "\nModelo={} | Topologia={} | Varrendo {} valores de \u{03BE} para B = {} G...",
                    model_name, topology_name, num_points, b_string
                );
                let all_results = Solver::solve_parallel(engines, 16);

                let base_dir = format!("output/modmax/{}/B_{}/{}", model_name, b_string, topology_name);

                for (i, results) in all_results.iter().enumerate() {
                    let csi = csi_vals[i];

                    let dir_path = format!("{}/csi_{:.2e}", base_dir, csi);
                    if let Err(_) = fs::create_dir_all(&dir_path) {
                        continue;
                    }

                    let eps: Vec<f64> = results.iter().map(|r| r[1]).collect();
                    let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();

                    let (masses, radii, central_p_list) = generate_mr_curve(&eps, &p_arr, false);

                    let eos_filename = format!("{}/eos.dat", dir_path);
                    if let Err(_) =
                        write_eos_with_mr(results, &masses, &radii, &central_p_list, &eos_filename)
                    {
                        continue;
                    }
                }
            }
        }
    }

    println!("\nProcesso concluído!");
}