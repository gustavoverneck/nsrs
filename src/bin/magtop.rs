// src/bin/magtop.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::{HadronsMatter, NlemModel, MagneticTopology};
use nsrs::core::solver::{Solver, EngineMode};
use nsrs::core::plotting::Artist;
use nsrs::core::tov_solver::generate_mr_curve;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};

fn main() {
    let mut args: Vec<String> = env::args().collect();
    
    let plot_only = if let Some(idx) = args.iter().position(|x| x == "--plot-only") {
        args.remove(idx);
        true
    } else {
        false
    };

    if args.len() < 3 {
        eprintln!("Uso: {} [--plot-only] <GM1|GM3> <B1> <B2> ...", args[0]);
        eprintln!("Exemplo: {} GM1 1e16 5e17 1e18", args[0]);
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

    fs::create_dir_all("results/magtop").unwrap_or_default();

    // Artists Combinados
    let mut artist_eos = Artist::new(
        &format!("results/magtop/eos_topology_{}.svg", model_name),
        &format!("Equation of State (Topology) - {}", model_name),
    )
    .with_x_label("Energy Density \u{03B5} [MeV/fm\u{00B3}]")
    .with_y_label("Pressure P [MeV/fm\u{00B3}]")
    .autoscale()
    .with_log_scale();

    let mut artist_mr = Artist::new(
        &format!("results/magtop/mr_topology_{}.svg", model_name),
        &format!("Mass-Radius (Topology) - {}", model_name),
    )
    .with_x_label("Radius [km]")
    .with_y_label("Mass [M\u{2299}]")
    .autoscale();

    let csv_summary = format!("results/magtop/summary_topology_{}.csv", model_name);
    let mut summary_data = Vec::new();

    for &b_field in &b_fields {
        let b_string = format!("{:.2e}", b_field);
        let base_dir = format!("output/magtop/{}/B_{}", model_name, b_string);
        let dir_iso = format!("{}/isotropic", base_dir);
        let dir_aniso = format!("{}/anisotropic", base_dir);

        if !plot_only {
            println!("Calculando EoS para B = {} G...", b_string);
            
            // Instancia os dois motores com as topologias diferentes
            // Nota: Fixamos NlemModel::Maxwell para isolar APENAS o efeito da topologia
            let engine_iso = EngineMode::Hadrons(
                HadronsMatter::new(model_params, b_field)
                    .with_topology(MagneticTopology::Isotropic)
                    .with_limits(0.01, 2.2)
                    .with_points(1500)
            );

            let engine_aniso = EngineMode::Hadrons(
                HadronsMatter::new(model_params, b_field)
                    .with_topology(MagneticTopology::Anisotropic)
                    .with_limits(0.01, 2.2)
                    .with_points(1500)
            );

            // Resolve as duas em paralelo
            let results = Solver::solve_batch(vec![engine_iso, engine_aniso]);

            if results.len() == 2 {
                fs::create_dir_all(&dir_iso).ok();
                Solver::write_eos(&results[0], &format!("{}/eos.dat", dir_iso)).ok();

                fs::create_dir_all(&dir_aniso).ok();
                Solver::write_eos(&results[1], &format!("{}/eos.dat", dir_aniso)).ok();
            } else {
                eprintln!("Aviso: resultados incompletos para B = {}", b_string);
            }
        }

        // Leitura e Plotagem
        let topologies = [("isotropic", "Iso"), ("anisotropic", "Aniso")];

        for (folder, label_suffix) in topologies.iter() {
            let eos_path = format!("{}/{}/eos.dat", base_dir, folder);
            
            if let Ok(file) = fs::File::open(&eos_path) {
                let reader = BufReader::new(file);
                let mut eps_vec = Vec::new();
                let mut p_vec = Vec::new();

                for line in reader.lines().flatten() {
                    let parts: Vec<f64> = line
                        .split_whitespace()
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect();
                    if parts.len() > 2 {
                        let eps = parts[1];
                        let p = parts[2];
                        if eps > 0.0 && p > 0.0 {
                            eps_vec.push(eps);
                            p_vec.push(p);
                        }
                    }
                }

                if eps_vec.len() > 10 {
                    // Passamos 'false' para não forçar a suavização pesada se não for necessário
                    let (masses, radii, _central_p_list) = generate_mr_curve(&eps_vec, &p_vec, false);
                    let label = format!("{} ({})", b_string, label_suffix);

                    artist_eos = artist_eos.add_curve(&eps_vec, &p_vec, &label);
                    artist_mr = artist_mr.add_curve(&radii, &masses, &label);

                    // Extrair Massa Máxima para o CSV
                    let mut m_max = 0.0;
                    let mut r_at_m_max = 0.0;
                    for i in 0..masses.len() {
                        if masses[i] > m_max && radii[i] < 16.0 {
                            m_max = masses[i];
                            r_at_m_max = radii[i];
                        }
                    }
                    summary_data.push(format!("{},{},{:.4},{:.4}", b_string, label_suffix, m_max, r_at_m_max));
                }
            }
        }
    }

    // Salvar gráficos
    artist_eos.plot().ok();
    artist_mr.plot().ok();

    // Salvar CSV de Resumo
    if let Ok(mut file) = fs::File::create(&csv_summary) {
        use std::io::Write;
        writeln!(file, "b_field,topology,max_mass,radius_at_max").ok();
        for line in summary_data {
            writeln!(file, "{}", line).ok();
        }
    }

    println!("\nAnálise concluída! Gráficos salvos em: results/magtop/");
}