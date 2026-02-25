use std::error::Error;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use nsrs::solver::{io_utils, plotting, tov_solver};

// Estrutura para armazenar os dados de cada configuração
struct StarSummary {
    model: String,
    b: String,
    csi: String,
    max_mass: f64,
    maxm_radius: f64,
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("--- Iniciando Processamento em Lote TOV ---");
    let total_start_time = Instant::now();

    let base_input_dir = Path::new("output"); 
    let base_output_dir = Path::new("results"); 
    
    // Garante que a pasta results base exista para o summary.csv
    fs::create_dir_all(&base_output_dir)?;

    let mut success_count = 0;
    let mut fail_count = 0;
    
    // Onde guardaremos a nossa lista para o Scatter Plot
    let mut summary_data: Vec<StarSummary> = Vec::new();

    for model_entry in fs::read_dir(base_input_dir)? {
        let model_entry = model_entry?;
        if !model_entry.file_type()?.is_dir() { continue; }
        let model_name = model_entry.file_name().to_string_lossy().to_string();

        for b_entry in fs::read_dir(model_entry.path())? {
            let b_entry = b_entry?;
            if !b_entry.file_type()?.is_dir() { continue; }
            let b_name = b_entry.file_name().to_string_lossy().to_string();

            for csi_entry in fs::read_dir(b_entry.path())? {
                let csi_entry = csi_entry?;
                if !csi_entry.file_type()?.is_dir() { continue; }
                let csi_name = csi_entry.file_name().to_string_lossy().to_string();

                let eos_path = csi_entry.path().join("eos.dat");

                if eos_path.exists() {
                    let target_dir = base_output_dir
                        .join(&model_name)
                        .join(&b_name)
                        .join(&csi_name);
                    
                    fs::create_dir_all(&target_dir)?;
                    let plot_path = target_dir.join("mr_curve.svg"); 

                    print!("Processando: {}/{}... ", &model_name, &b_name);

                    match process_single_eos(&eos_path, &plot_path) {
                        Ok((max_mass, maxm_radius)) => {
                            println!("OK! (M_max: {:.2} M☉, R: {:.2} km)", max_mass, maxm_radius);
                            
                            // Salvando na nossa lista
                            summary_data.push(StarSummary {
                                model: model_name.clone(),
                                b: b_name.clone(),
                                csi: csi_name.clone(),
                                max_mass,
                                maxm_radius,
                            });
                            
                            success_count += 1;
                        },
                        Err(e) => {
                            println!("ERRO: {}", e);
                            fail_count += 1;
                        }
                    }
                }
            }
        }
    }

    // --- Exportando os dados salvos para CSV ---
    let csv_path = base_output_dir.join("summary.csv");
    let mut file = File::create(&csv_path)?;
    
    // Escreve o cabeçalho
    writeln!(file, "model,b,csi,max_mass,max_radius")?;
    
    // Escreve os dados
    for entry in &summary_data {
        // Converte a notação Fortran (d) para a notação Científica Padrão (e)
        let clean_b = entry.b.replace("d", "e").replace("D", "e");
        let clean_csi = entry.csi.replace("d", "e").replace("D", "e");

        writeln!(
            file, 
            "{},{},{},{:.6},{:.6}", 
            entry.model, clean_b, clean_csi, entry.max_mass, entry.maxm_radius
        )?;
    }

    let duration = total_start_time.elapsed();
    println!("\n--- Processamento Finalizado ---");
    println!("Tempo Total: {:.2?}", duration);
    println!("Sucessos: {}", success_count);
    println!("Falhas: {}", fail_count);
    println!("Tabela salva em: {}", csv_path.display());

    Ok(())
}

fn process_single_eos(input_eos_path: &Path, output_plot_path: &Path) -> Result<(f64, f64), Box<dyn Error>> {
    let (eps, p) = io_utils::read_eos_file(input_eos_path)
        .map_err(|e| format!("Falha ao ler {}: {}", input_eos_path.display(), e))?;

    if eps.is_empty() {
        return Err("Arquivo EoS vazio.".into());
    }

    let (masses, radii) = tov_solver::generate_mr_curve(&eps, &p);

    if masses.is_empty() {
        return Err("Nenhuma configuração estelar válida foi gerada.".into());
    }

    // Encontrando a Massa Máxima e o Raio correspondente
    let mut max_mass = 0.0;
    let mut maxm_radius = 0.0;

    for (i, &m) in masses.iter().enumerate() {
        if m > max_mass {
            max_mass = m;
            maxm_radius = radii[i];
        }
    }

    plotting::plot_mr_curve(
        &radii, 
        &masses, 
        output_plot_path.to_str().unwrap()
    ).map_err(|e| format!("Falha ao plotar gráfico: {}", e))?;

    // Retorna a tupla com os dois valores
    Ok((max_mass, maxm_radius))
}