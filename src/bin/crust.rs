// src/bin/crust.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::PhysicsEngine;
use nsrs::core::solver::Solver;
use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::Artist; 
use nsrs::core::tov_solver::{generate_mr_curve, unify_with_crust}; // <-- Importamos unify_with_crust
use std::env;
use std::fs;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Uso: {} <GM1|GM3> <B1> <B2> ...", args[0]);
        return;
    }
    
    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let b_fields: Vec<f64> = args[2..]
        .iter()
        .map(|s| s.parse().expect("B deve ser um número válido"))
        .collect();

    // 1. Prepara os Engines
    let engines: Vec<PhysicsEngine> = b_fields.iter().map(|&bg| {
        PhysicsEngine::new(model_params, bg)
            .with_limits(0.02, 2.5)
            .with_points(2000)
    }).collect();

    // 2. Executa o cálculo paralelo
    println!("Calculando {} Equações de Estado em paralelo...", b_fields.len());
    let all_results = Solver::solve_batch(engines);

    // Inicializa Artists para EoS e MR
    let mut eos_artist = Artist::new(
        &format!("results/crust_effect_eos_{}.svg", model_name),
        &format!("Crust Effect on EoS - {}", model_name),
    )
    .with_x_label("Energy Density log \u{3B5} [MeV/fm\u{00B3}]")
    .with_y_label("Pressure log P [MeV/fm\u{00B3}]")
    .with_log_scale() 
    .autoscale(); 

    let mut mr_artist = Artist::new(
        &format!("results/crust_effect_mr_{}.svg", model_name),
        &format!("Crust Effect on Mass-Radius - {}", model_name),
    )
    .with_x_label("Radius [km]")
    .with_y_label("Mass [M\u{2299}]") 
    .with_x_range(10.0, 15.0); 

    // 3. Processamento e Exportação Organizada
    for (i, results) in all_results.iter().enumerate() {
        let bg = b_fields[i];
        let b_string = format!("{:.2e}", bg);
        let dir_path = format!("output/{}/{}/", model_name, b_string);
        
        if let Err(e) = fs::create_dir_all(&dir_path) {
            eprintln!("Erro ao criar diretório {}: {}", dir_path, e);
            continue;
        }

        let eos_filename = format!("{}eos.dat", dir_path);

        if let Err(e) = Solver::write_eos(results, &eos_filename) {
            eprintln!("Erro ao salvar EOS em {}: {}", eos_filename, e);
            continue;
        }

        match read_eos_file(&eos_filename) {
            Ok((eps, p)) => {
                let label_nc = format!("B = {} G (Core)", b_string);
                let label_wc = format!("B = {} G (Crust)", b_string);

                // --- PLOT DA EOS ---
                // Adiciona EoS original (Só Núcleo)
                eos_artist = eos_artist.add_curve(&eps, &p, &label_nc);
                
                // Unifica a EoS com a Crosta e adiciona ao plot
                let (unified_eps, unified_p) = unify_with_crust(&eps, &p);
                eos_artist = eos_artist.add_curve(&unified_eps, &unified_p, &label_wc);

                // --- PLOT E CÁLCULO MR ---
                // Integração SEM Crosta
                let (masses_nc, radii_nc) = generate_mr_curve(&eps, &p, false);
                if !masses_nc.is_empty() {
                    let mr_nc_filename = format!("{}mr_no_crust.dat", dir_path);
                    save_mr_data(&radii_nc, &masses_nc, &mr_nc_filename).ok();
                    mr_artist = mr_artist.add_curve(&radii_nc, &masses_nc, &label_nc);
                }

                // Integração COM Crosta
                let (masses_wc, radii_wc) = generate_mr_curve(&eps, &p, true);
                if !masses_wc.is_empty() {
                    let mr_wc_filename = format!("{}mr_with_crust.dat", dir_path);
                    save_mr_data(&radii_wc, &masses_wc, &mr_wc_filename).ok();
                    mr_artist = mr_artist.add_curve(&radii_wc, &masses_wc, &label_wc);
                }
            }
            Err(e) => eprintln!("Erro ao processar dados para B={}: {}", b_string, e),
        }
    }

    // 4. Gera Gráficos
    eos_artist.plot().ok();
    mr_artist.plot().ok();
    
    println!("Exportação concluída!");
    println!("Gráficos salvos em:");
    println!("  - /results/crust_effect_eos_{}.svg", model_name);
    println!("  - /results/crust_effect_mr_{}.svg", model_name);
}

fn save_mr_data(radii: &[f64], masses: &[f64], filename: &str) -> std::io::Result<()> {
    let mut file = fs::File::create(filename)?;
    use std::io::Write;
    for (r, m) in radii.iter().zip(masses.iter()) {
        writeln!(file, "{:12.5e} {:12.5e}", r, m)?;
    }
    Ok(())
}