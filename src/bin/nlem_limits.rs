// src/bin/nlem_limits.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::{PhysicsEngine, NlemModel};
use nsrs::core::solver::Solver;
use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::Artist;
use nsrs::core::tov_solver::generate_mr_curve; 
use std::env;
use std::fs;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!("Uso: {} <GM1|GM3> <Campo B> <exp_min> <exp_max>", args[0]);
        eprintln!("Exemplo (10^8 a 10^28): {} GM1 5e16 8 28", args[0]);
        return;
    }
    
    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let b_field: f64 = args[2].parse().expect("Campo B inválido");
    let exp_min: i32 = args[3].parse().expect("exp_min inválido (use inteiros como 8)");
    let exp_max: i32 = args[4].parse().expect("exp_max inválido (use inteiros como 28)");

    // 1. Gera os valores de csi intercalando 1.0eX e 5.0eX
    let mut log_csi_vals = Vec::new();
    let mut csi_vals = Vec::new();

    for exp in exp_min..=exp_max {
        // Valor 1.0 * 10^exp
        let val1 = 10_f64.powi(exp);
        csi_vals.push(val1);
        log_csi_vals.push(val1.log10());

        // Valor 5.0 * 10^exp (adiciona exceto se já estourou o limite máximo desejado)
        if exp < exp_max {
            let val5 = 5.0 * 10_f64.powi(exp);
            csi_vals.push(val5);
            log_csi_vals.push(val5.log10()); // log10(5eX) ≈ X + 0.69897
        }
    }

    let num_points = csi_vals.len();

    // 2. Prepara os Engines para processamento paralelo
    let engines: Vec<PhysicsEngine> = csi_vals.iter().map(|&csi| {
        PhysicsEngine::new(model_params, b_field)
            .with_nlem(NlemModel::Log(csi)) // Fixado no modelo Log
            .with_limits(0.01, 2.5)
            .with_points(1200)
    }).collect();

    // 3. Executa o cálculo paralelo
    println!("Varrendo {} valores de \u{03BE} para B = {:.1e} G...", num_points, b_field);
    let all_results = Solver::solve_batch(engines);

    let mut max_masses = Vec::new();
    let mut radii_at_max = Vec::new();
    let mut valid_log_csi = Vec::new();

    let b_string = format!("{:.2e}", b_field);
    let base_dir = format!("output/limits/{}/B_{}", model_name, b_string);

    // 4. Processamento, Integração TOV e Extração de Limites
    for (i, results) in all_results.iter().enumerate() {
        let csi = csi_vals[i];
        let log_csi = log_csi_vals[i];
        
        let dir_path = format!("{}/csi_{:.2e}", base_dir, csi);
        if let Err(_) = fs::create_dir_all(&dir_path) { continue; }
        
        let eos_filename = format!("{}/eos.dat", dir_path);
        if let Err(_) = Solver::write_eos(results, &eos_filename) { continue; }

        if let Ok((eps, p)) = read_eos_file(&eos_filename) {
            // Gera a curva MR (com crosta ativada para maior precisão de raio)
            let (masses, radii) = generate_mr_curve(&eps, &p, true);
            
            if !masses.is_empty() {
                let mut m_max = 0.0;
                let mut r_max_m = 0.0;
                
                // Encontra a Massa Máxima e o respectivo Raio
                for (r, m) in radii.iter().zip(masses.iter()) {
                    if *m > m_max {
                        m_max = *m;
                        r_max_m = *r;
                    }
                }
                
                max_masses.push(m_max);
                radii_at_max.push(r_max_m);
                valid_log_csi.push(log_csi);
            }
        }
    }

    // A quantidade de intervalos é exatamente a diferença entre o exp máximo e mínimo
    let x_ticks_count = (exp_max - exp_min) as usize;

    let artist_m = Artist::new(
        &format!("results/limit_mass_{}_B_{}.svg", model_name, b_string),
        &format!("Maximum Mass Limit (Log Model) - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("Maximum Mass [M\u{2299}]")
    .autoscale() // Deixa o eixo Y se ajustar sozinho...
    .with_x_range(exp_min as f64, exp_max as f64) // ...mas trava o eixo X aos inputs!
    .with_x_labels(x_ticks_count) // Pede 1 marcação por unidade de expoente
    .add_curve(&valid_log_csi, &max_masses, &format!("B = {} G", b_string));
    
    artist_m.plot().ok();

    let artist_r = Artist::new(
        &format!("results/limit_radius_{}_B_{}.svg", model_name, b_string),
        &format!("Radius at Maximum Mass (Log Model) - {}", model_name),
    )
    .with_x_label("log\u{2081}\u{2080}(\u{03BE}) [G]")
    .with_y_label("Radius [km]")
    .autoscale()
    .with_x_range(exp_min as f64, exp_max as f64) // Trava o eixo X aos inputs
    .with_x_labels(x_ticks_count) // Pede 1 marcação por unidade de expoente
    .add_curve(&valid_log_csi, &radii_at_max, &format!("B = {} G", b_string));
    
    artist_r.plot().ok();

    // 6. Exporta os Dados para CSV
    let csv_path = format!("results/limits_data_{}_B_{}.csv", model_name, b_string);
    if let Ok(mut file) = fs::File::create(&csv_path) {
        use std::io::Write;
        writeln!(file, "csi,log10_csi,max_mass_msun,radius_at_max_km").ok();
        for i in 0..valid_log_csi.len() {
            writeln!(file, "{:e},{:.4},{:.4},{:.4}", 
                10_f64.powf(valid_log_csi[i]), 
                valid_log_csi[i], 
                max_masses[i], 
                radii_at_max[i]
            ).ok();
        }
    }

    println!("Análise de limites concluída!");
    println!("Gráficos gerados:");
    println!("  - /results/limit_mass_{}_B_{}.svg", model_name, b_string);
    println!("  - /results/limit_radius_{}_B_{}.svg", model_name, b_string);
    println!("Dados brutos salvos em: {}", csv_path);
}