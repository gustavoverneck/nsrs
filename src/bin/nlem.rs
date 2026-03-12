// src/bin/nlem.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::{PhysicsEngine, NlemModel}; // Importamos NlemModel
use nsrs::core::solver::Solver;
use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::Artist;
use nsrs::core::tov_solver::generate_mr_curve; 
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

    // Parâmetros arbitrários para os modelos NLEM
    // ModMax: csi é adimensional. Log: csi tem unidade de campo (Gauss)
    let csi_modmax = 0.5; 
    let csi_log = 5e16; 

    // 1. Prepara os Engines: 3 modelos para cada campo magnético!
    let mut engines = Vec::new();
    let mut configs = Vec::new();

    for &bg in &b_fields {
        let models_to_test = vec![
            (NlemModel::Maxwell, "Maxwell", format!("Maxwell (B = {:.1e} G)", bg)),
            (NlemModel::Modmax(csi_modmax), "ModMax", format!("ModMax csi={} (B = {:.1e} G)", csi_modmax, bg)),
            (NlemModel::Log(csi_log), "Log", format!("Log csi={:.1e} (B = {:.1e} G)", csi_log, bg)),
        ];

        for (nlem, dir_name, label) in models_to_test {
            engines.push(
                PhysicsEngine::new(model_params, bg)
                    .with_nlem(nlem)
                    .with_limits(0.0, 2.5)
                    .with_points(1201) // Pode aumentar para 2200 para mais precisão
            );
            configs.push((bg, dir_name, label));
        }
    }

    // 2. Executa o cálculo paralelo
    println!("Calculando {} Equações de Estado em paralelo (3 modelos por B)...", engines.len());
    let all_results = Solver::solve_batch(engines);

    // Inicializa Artists para os gráficos comparativos gerais
    let mut eos_artist = Artist::new(
        &format!("results/compare_nlem_eos_{}.svg", model_name),
        &format!("NLEM Effect on EoS - {}", model_name),
    )
    .with_x_label("Energy Density \u{3B5} [MeV/fm\u{00B3}]")
    .with_y_label("Pressure P [MeV/fm\u{00B3}]")
    .with_log_scale()
    .autoscale();

    let mut mr_artist = Artist::new(
        &format!("results/compare_nlem_mr_{}.svg", model_name),
        &format!("NLEM Effect on Mass-Radius - {}", model_name),
    )
    .with_x_label("Radius [km]")
    .with_y_label("Mass [M\u{2299}]") // Símbolo solar ⊙
    .with_x_range(10.0, 16.0);
    // .autoscale();

    // 3. Processamento e Exportação Organizada
    for (i, results) in all_results.iter().enumerate() {
        let (bg, dir_name, label) = &configs[i];
        let b_string = format!("{:.2e}", bg);

        // Define o caminho da pasta: output/modelo/B/NLEM_MODEL/
        let dir_path = format!("output/{}/{}/{}/", model_name, b_string, dir_name);
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
                // Adiciona a curva no EoS
                eos_artist = eos_artist.add_curve(&eps, &p, label);

                // Integra a TOV com a crosta unificada (true)
                let (masses, radii, _) = generate_mr_curve(&eps, &p, true);
                if !masses.is_empty() {
                    // Exporta mr.dat bruto
                    save_mr_data(&radii, &masses, &mr_filename).ok();
                    // Adiciona a curva no gráfico M-R
                    mr_artist = mr_artist.add_curve(&radii, &masses, label);
                }
            }
            Err(e) => eprintln!("Erro ao processar dados para {}: {}", label, e),
        }
    }

    // 4. Gera Gráficos Comparativos
    eos_artist.plot().ok();
    mr_artist.plot().ok();
    println!("Exportação concluída! Verifique as pastas /output e /results.");
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