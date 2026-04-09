// src/bin/nlem_limits.rs

use nsrs::{
    EngineMode, GM1, GM3, HadronsMatter, MagneticTopology, NlemModel, Solver, generate_mr_curve,
    write_eos_with_mr,
};
use std::env;
use std::fs;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Uso: {} <exp_min> <exp_max> <interval> <B1> <B2> ...", args[0]);
        eprintln!("Exemplo: {} 8 28 2 1e15 5e16", args[0]);
        return;
    }

    let exp_min: i32 = args[1].parse().expect("exp_min inválido (use inteiros como 8)");
    let exp_max: i32 = args[2].parse().expect("exp_max inválido (use inteiros como 28)");
    let interval: usize = args[3]
        .parse()
        .expect("interval inválido (use inteiros como 1, 2, 5 ou 10)");

    if interval == 0 {
        panic!("interval deve ser >= 1");
    }
    
    let b_fields: Vec<f64> = args[4..]
        .iter()
        .map(|s| s.parse().expect("B deve ser um número válido"))
        .collect();

    // 1. Gera os valores de csi com "interval" pontos por década.
    // Ex.: interval=1 -> [10^e]
    //      interval=2 -> [10^e, 10^(e+0.5)]
    //      interval=5 -> [10^(e+k/5), k=0..4]
    //      interval=10 -> [10^(e+k/10), k=0..9]
    let mut csi_vals = Vec::new();

    for exp in exp_min..exp_max {
        for k in 0..interval {
            let log_csi = exp as f64 + (k as f64 / interval as f64);
            csi_vals.push(10_f64.powf(log_csi));
        }
    }

    // adiciona o ponto final 10^exp_max
    let log_csi_last = exp_max as f64;
    csi_vals.push(10_f64.powf(log_csi_last));

    let num_points = csi_vals.len();
    let models = [("GM1", GM1), ("GM3", GM3)];
    let topologies = [
        ("isotropic", MagneticTopology::Isotropic),
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
                            .with_nlem(NlemModel::Log(csi))
                            .with_limits(0.01, 2.0)
                            .with_points(2000);

                        EngineMode::Hadrons(motor)
                    })
                    .collect();

                println!(
                    "\nModelo={} | Topologia={} | Varrendo {} valores de \u{03BE} para B = {} G...",
                    model_name, topology_name, num_points, b_string
                );
                let all_results = Solver::solve_batch(engines);

                let base_dir = format!("output/limits/{}/B_{}/{}", model_name, b_string, topology_name);

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