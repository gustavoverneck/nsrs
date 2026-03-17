use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::HadronsMatter;
use nsrs::core::solver::{Solver, EngineMode};
use std::env;
use std::fs;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Uso: {} <GM1|GM3> <B1> <B2> ...", args[0]);
        eprintln!("Exemplo: {} GM1 1e16 1e17 1e18", args[0]);
        return;
    }
    
    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1, 
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let b_fields: Vec<f64> = args[2..].iter().map(|s| s.parse().expect("Campo B inválido")).collect();
    
    // Lista de valores de csi para a varredura (0.0 é o Padrão)
    let csi_vals = vec![0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e1, 1e2];

    for &b_field in &b_fields {
        let b_string = format!("{:.2e}", b_field);
        println!("Calculando EoS para B = {} G...", b_string);

        let mut engines = Vec::new();

        // Cria os motores para cada valor de csi
        for &csi in &csi_vals {
            let motor = HadronsMatter::new(model_params, b_field)
                .with_csi(csi)
                .with_limits(0.01, 2.5)
                .with_points(2000);
            
            engines.push(EngineMode::Hadrons(motor));
        }

        // Resolve todos os csi em paralelo para este campo B
        let results = Solver::solve_batch(engines);

        // Exporta os resultados
        for (i, &csi) in csi_vals.iter().enumerate() {
            if results[i].is_empty() {
                println!("  [!] Aviso: A EoS falhou a convergência para B = {} e csi = {}.", b_string, csi);
            }

            // Formata a pasta: 0.0 ou notação científica (ex: 1.0e-5)
            let csi_str = if csi == 0.0 {
                "0.0".to_string()
            } else {
                format!("{:.1e}", csi)
            };

            let dir_path = format!("output/landau/{}/B_{}/csi_{}", model_name, b_string, csi_str);
            fs::create_dir_all(&dir_path).unwrap_or_default();
            
            let eos_path = format!("{}/eos.dat", dir_path);
            Solver::write_eos(&results[i], &eos_path).unwrap();
        }

        println!("  ✔ Dados salvos para B = {} G em output/landau/{}/", b_string, model_name);
    }
}