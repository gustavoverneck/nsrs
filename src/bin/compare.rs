use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::PhysicsEngine;
use nsrs::core::solver::Solver;
use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::Artist; 
use nsrs::core::tov_solver::generate_mr_curve;
use std::env;
use std::fs;
use std::path::Path;

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

    // 1. Prepara os Engines para processamento paralelo
    let engines: Vec<PhysicsEngine> = b_fields.iter().map(|&bg| {
        PhysicsEngine::new(model_params, bg)
            .with_limits(0.01, 2.5)
            .with_points(2200)
    }).collect();

    // 2. Executa o cálculo paralelo
    println!("Calculando {} Equações de Estado em paralelo...", b_fields.len());
    let all_results = Solver::solve_batch(engines);

    // Inicializa Artists para os gráficos comparativos gerais
    let mut eos_artist = Artist::new(
        &format!("results/compare_eos_{}.svg", model_name),
        &format!("Equation of State - {}", model_name),
    )
    .with_x_label("Energy Density \u{3B5} [MeV/fm\u{00B3}]"); // ε [MeV/fm³]


    let mut mr_artist = Artist::new(
        &format!("results/compare_mr_{}.svg", model_name),
        &format!("Mass-Radius Relation - {}", model_name),
    )
    .with_x_label("Radius [km]")
    .with_y_label("Mass [M\u{2299}]") // Símbolo solar ⊙
    .with_x_range(8.0, 14.0); // Fixa o range para estrelas de nêutrons

    // 3. Processamento e Exportação Organizada
    for (i, results) in all_results.iter().enumerate() {
        let bg = b_fields[i];
        let b_string = format!("{:.2e}", bg);
        let label = format!("B = {} G", b_string);

        // Define o caminho da pasta: output/modelo/B/
        let dir_path = format!("output/{}/{}/", model_name, b_string);
        if let Err(e) = fs::create_dir_all(&dir_path) {
            eprintln!("Erro ao criar diretório {}: {}", dir_path, e);
            continue;
        }

        let eos_filename = format!("{}eos.dat", dir_path);
        let mr_filename = format!("{}mr.dat", dir_path);

        // Salva o arquivo EOS
        if let Err(e) = Solver::write_eos(results, &eos_filename) {
            eprintln!("Erro ao salvar EOS em {}: {}", eos_filename, e);
            continue;
        }

        // Processa TOV e salva MR
        match read_eos_file(&eos_filename) {
            Ok((eps, p)) => {
                eos_artist = eos_artist.add_curve(&eps, &p, &label);

                let (masses, radii) = generate_mr_curve(&eps, &p);
                if !masses.is_empty() {
                    // Exporta mr.dat bruto
                    save_mr_data(&radii, &masses, &mr_filename).expect("Falha ao salvar mr.dat");
                    mr_artist = mr_artist.add_curve(&radii, &masses, &label);
                }
            }
            Err(e) => eprintln!("Erro ao processar dados para {}: {}", label, e),
        }
    }

    // 4. Gera Gráficos Comparativos
    eos_artist.plot().ok();
    mr_artist.plot().ok();
    println!("Exportação concluída em /output e /results.");
}

/// Função auxiliar para salvar o arquivo mr.dat formatado
fn save_mr_data(radii: &[f64], masses: &[f64], filename: &str) -> std::io::Result<()> {
    let mut file = fs::File::create(filename)?;
    use std::io::Write;
    for (r, m) in radii.iter().zip(masses.iter()) {
        writeln!(file, "{:12.5e} {:12.5e}", r, m)?;
    }
    Ok(())
}