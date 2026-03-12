// src/bin/nlem_limits.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::{PhysicsEngine, NlemModel};
use nsrs::core::solver::Solver;
use nsrs::core::plotting::Artist;
use nsrs::core::tov_solver::generate_mr_curve; 
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};

fn main() {
    let mut args: Vec<String> = env::args().collect();
    
    // Captura a tag --plot-only e a remove dos argumentos
    let plot_only = if let Some(idx) = args.iter().position(|x| x == "--plot-only") {
        args.remove(idx);
        true
    } else {
        false
    };

    if args.len() < 5 {
        eprintln!("Uso: {} [--plot-only] <GM1|GM3> <exp_min> <exp_max> <B1> <B2> ...", args[0]);
        eprintln!("Exemplo: {} --plot-only GM1 8 28 1e15 5e16", args[0]);
        return;
    }
    
    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let exp_min: i32 = args[2].parse().expect("exp_min inválido (use inteiros como 8)");
    let exp_max: i32 = args[3].parse().expect("exp_max inválido (use inteiros como 28)");
    
    let b_fields: Vec<f64> = args[4..]
        .iter()
        .map(|s| s.parse().expect("B deve ser um número válido"))
        .collect();

    // 1. Gera os valores de csi intercalando 1.0eX e 5.0eX
    let mut log_csi_vals = Vec::new();
    let mut csi_vals = Vec::new();

    for exp in exp_min..=exp_max {
        let val1 = 10_f64.powi(exp);
        csi_vals.push(val1);
        log_csi_vals.push(val1.log10());

        if exp < exp_max {
            let val5 = 5.0 * 10_f64.powi(exp);
            csi_vals.push(val5);
            log_csi_vals.push(val5.log10()); 
        }
    }

    let num_points = csi_vals.len();
    let x_ticks_count = (exp_max - exp_min) as usize;

    // 2. Prepara os Artists (Agora com 5 gráficos!)
    let mut artist_m = Artist::new(
        &format!("results/limit_mass_{}_combined.svg", model_name),
        &format!("Maximum Mass Limit - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("Maximum Mass [M\u{2299}]")
    .autoscale().with_x_range(exp_min as f64, exp_max as f64).with_x_labels(x_ticks_count);

    let mut artist_r = Artist::new(
        &format!("results/limit_radius_{}_combined.svg", model_name),
        &format!("Radius at Maximum Mass - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("Radius [km]")
    .autoscale().with_x_range(exp_min as f64, exp_max as f64).with_x_labels(x_ticks_count);

    let mut artist_nc = Artist::new(
        &format!("results/limit_nc_{}_combined.svg", model_name),
        &format!("Central Baryon Density - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("Central Density [n_c / n\u{2080}]")
    .autoscale().with_x_range(exp_min as f64, exp_max as f64).with_x_labels(x_ticks_count);

    let mut artist_meff = Artist::new(
        &format!("results/limit_meff_{}_combined.svg", model_name),
        &format!("Central Effective Mass - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("m* / m_N")
    .autoscale().with_x_range(exp_min as f64, exp_max as f64).with_x_labels(x_ticks_count);

    let mut artist_emag = Artist::new(
        &format!("results/limit_emag_{}_combined.svg", model_name),
        &format!("Central Magnetic Energy - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("\u{03B5}_mag [MeV/fm\u{00B3}]")
    .autoscale().with_x_range(exp_min as f64, exp_max as f64).with_x_labels(x_ticks_count);

    // 3. Loop principal sobre os campos magnéticos
    for &b_field in &b_fields {
        let b_string = format!("{:.2e}", b_field);
        
        let mut max_masses = Vec::new();
        let mut radii_at_max = Vec::new();
        let mut valid_log_csi = Vec::new();
        
        // Novos vetores para a microfísica central
        let mut central_nc_vals = Vec::new();
        let mut central_meff_vals = Vec::new();
        let mut central_emag_vals = Vec::new();

        let csv_path = format!("results/limits_data_{}_B_{}.csv", model_name, b_string);

        if plot_only {
            println!("Lendo dados salvos de {}...", csv_path);
            if let Ok(file) = fs::File::open(&csv_path) {
                let reader = BufReader::new(file);
                for line in reader.lines().skip(1) {
                    if let Ok(content) = line {
                        let parts: Vec<&str> = content.split(',').collect();
                        if parts.len() >= 7 { // Agora lemos 7 colunas!
                            if let (Ok(log_csi), Ok(m_max), Ok(r_max), Ok(nc), Ok(meff), Ok(emag)) = (
                                parts[1].parse::<f64>(),
                                parts[2].parse::<f64>(),
                                parts[3].parse::<f64>(),
                                parts[4].parse::<f64>(),
                                parts[5].parse::<f64>(),
                                parts[6].parse::<f64>()
                            ) {
                                valid_log_csi.push(log_csi);
                                max_masses.push(m_max);
                                radii_at_max.push(r_max);
                                central_nc_vals.push(nc);
                                central_meff_vals.push(meff);
                                central_emag_vals.push(emag);
                            }
                        }
                    }
                }
            } else {
                eprintln!("  -> ERRO: Arquivo {} não encontrado.", csv_path);
                continue;
            }
        } else {
            let engines: Vec<PhysicsEngine> = csi_vals.iter().map(|&csi| {
                PhysicsEngine::new(model_params, b_field)
                    .with_nlem(NlemModel::Log(csi))
                    .with_limits(0.01, 2.5)
                    .with_points(2000)
            }).collect();

            println!("\nVarrendo {} valores de \u{03BE} para B = {} G...", num_points, b_string);
            let all_results = Solver::solve_batch(engines);

            let base_dir = format!("output/limits/{}/B_{}", model_name, b_string);

            for (i, results) in all_results.iter().enumerate() {
                let csi = csi_vals[i];
                let log_csi = log_csi_vals[i];
                
                let dir_path = format!("{}/csi_{:.2e}", base_dir, csi);
                if let Err(_) = fs::create_dir_all(&dir_path) { continue; }
                
                let eos_filename = format!("{}/eos.dat", dir_path);
                if let Err(_) = Solver::write_eos(results, &eos_filename) { continue; }

                // Otimização: Extrai eps e p direto da memória em vez de ler o arquivo!
                let eps: Vec<f64> = results.iter().map(|r| r[1]).collect();
                let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();
                
                let (masses, radii, central_p_list) = generate_mr_curve(&eps, &p_arr, true);

                if !masses.is_empty() {
                    let mut m_max = 0.0;
                    let mut r_max_m = 0.0;
                    let mut best_pc = 0.0; // Armazena a pressão central da estrela de massa máxima
                    
                    for i_star in 0..masses.len() {
                        let r = radii[i_star];
                        let m = masses[i_star];
                        let pc = central_p_list[i_star];

                        if r < 15.0 && m > m_max {
                            m_max = m;
                            r_max_m = r;
                            best_pc = pc;
                        }
                    }
                    
                    if m_max > 0.0 {
                        // 2. Busca o índice exato no array 'results' que corresponde a esta pressão central
                        // Isso garante que você pegue a microfísica do NÚCLEO da estrela correta.
                        let core_idx = results.iter()
                            .position(|r| (r[2] - best_pc).abs() < 1e-8)
                            .unwrap_or(results.len() - 1);

                        max_masses.push(m_max);
                        radii_at_max.push(r_max_m);
                        valid_log_csi.push(log_csi);
                        
                        // 3. Extrai a microfísica usando o índice vinculado à Pc
                        central_nc_vals.push(results[core_idx][0]);  // n/n0
                        central_meff_vals.push(results[core_idx][16]); // m*/mN
                        central_emag_vals.push(results[core_idx][19]); // Energia Magnética
                    }
                }
            }

            // Exporta os Dados expandidos para CSV
            if let Ok(mut file) = fs::File::create(&csv_path) {
                use std::io::Write;
                writeln!(file, "csi,log10_csi,max_mass_msun,radius_at_max_km,central_nc,central_meff,central_emag").ok();
                for i in 0..valid_log_csi.len() {
                    writeln!(file, "{:e},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4e}", 
                        10_f64.powf(valid_log_csi[i]), 
                        valid_log_csi[i], 
                        max_masses[i], 
                        radii_at_max[i],
                        central_nc_vals[i],
                        central_meff_vals[i],
                        central_emag_vals[i]
                    ).ok();
                }
            }
            println!("  -> OK! Dados CSV salvos em: {}", csv_path);
        }

        let label = format!("B = {} G", b_string);
        artist_m = artist_m.add_curve(&valid_log_csi, &max_masses, &label);
        artist_r = artist_r.add_curve(&valid_log_csi, &radii_at_max, &label);
        artist_nc = artist_nc.add_curve(&valid_log_csi, &central_nc_vals, &label);
        artist_meff = artist_meff.add_curve(&valid_log_csi, &central_meff_vals, &label);
        artist_emag = artist_emag.add_curve(&valid_log_csi, &central_emag_vals, &label);
    }

    // 4. Gera os 5 Gráficos Finais Combinados
    artist_m.plot().ok();
    artist_r.plot().ok();
    artist_nc.plot().ok();
    artist_meff.plot().ok();
    artist_emag.plot().ok();

    println!("\nProcesso concluído!");
    println!("Gráficos gerados na pasta results/ limit_mass, limit_radius, limit_nc, limit_meff e limit_emag!");
}