// src/bin/nlem_modmax.rs

use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, NlemModel, Solver};
use std::env;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = env::args().collect();

    // Gera valores de 1e-10 .. 9e-1 (mantissas 1..9 por década)
    let csi_vals: Vec<f64> = (-10..=-1)
        .flat_map(|exp| {
            let base = 10_f64.powi(exp);
            (1..=9).map(move |i| i as f64 * base)
        })
        .collect();

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

    let models = [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)];

    // Loop principal: modelos x campos
    for (model_name, model_params) in models {
        for &b_field in &b_fields {
            let b_string = format_sci(b_field);

            // Estrutura de diretórios simplificada (sem topologia)
            let base_dir = format!("output/modmax/{}/B_{}", model_name, b_string);
            if fs::create_dir_all(&base_dir).is_err() {
                continue;
            }

            // Criação do arquivo de sumário
            let summary_path = format!("{}/summary.csv", base_dir);
            let mut summary = match fs::File::create(&summary_path) {
                Ok(file) => file,
                Err(_) => continue,
            };
            let _ = writeln!(summary, "label,csi,b_field,eos_file");

            // Baseline sem NLEM Modmax para comparação
            let baseline_filename = "eos_baseline.dat";
            let baseline_path = format!("{}/{}", base_dir, baseline_filename);
            let baseline_motor = HadronsMatter::new(model_params, b_field)
                .with_limits(0.02, 2.0)
                .with_points(1201)
                .with_eos_output(&baseline_path);
            let mut baseline_solver = Solver::new(EngineMode::Hadrons(baseline_motor));
            let _ = baseline_solver.solve();
            let _ = writeln!(summary, "hadrons_baseline,,,{}", baseline_filename);

            // Preparação das engines para varredura do Modmax
            let mut engines = Vec::new();
            for csi in &csi_vals {
                let eos_filename = format!("eos_csi_{}.dat", format_sci(*csi));
                let eos_path = format!("{}/{}", base_dir, eos_filename);

                let motor = HadronsMatter::new(model_params, b_field)
                    .with_nlem(NlemModel::Modmax(*csi))
                    .with_limits(0.02, 2.0)
                    .with_points(1201)
                    .with_eos_output(&eos_path);

                engines.push(EngineMode::Hadrons(motor));

                // Escreve os metadados no csv
                let _ = writeln!(
                    summary,
                    "modmax,{:.6e},{:.6e},{}",
                    csi, b_field, eos_filename
                );
            }

            println!(
                "\nModelo={} | B={} G | Varrendo {} valores de \u{03BE}...",
                model_name,
                b_string,
                engines.len()
            );

            // Roda todas as configurações em paralelo
            let _ = Solver::solve_parallel(engines, 16);
        }
    }

    println!("\nProcesso concluído! Dados salvos em output/modmax/");
}

// Funções utilitárias de formatação
fn format_sci(value: f64) -> String {
    format!("{:.2e}", value).replace('+', "")
}
