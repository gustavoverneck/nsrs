// src/bin/magtop.rs

use nsrs::core::model::{FSU2, GM1, GM3};
use nsrs::core::physics::{HadronsMatter, MagneticTopology};
use nsrs::core::plotting::Artist;
use nsrs::core::solver::{EngineMode, Solver};
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

    // NOVO: prefixo configurável
    let prefix = if let Some(idx) = args.iter().position(|x| x == "--prefix") {
        if idx + 1 >= args.len() {
            eprintln!("Erro: use --prefix <nome>");
            return;
        }
        let p = args[idx + 1].clone();
        args.drain(idx..=idx + 1);
        p
    } else {
        String::new()
    };

    let file_prefix = if prefix.is_empty() {
        String::new()
    } else {
        format!("{}_", prefix)
    };

    if args.len() < 3 {
        eprintln!(
            "Uso: {} [--plot-only] [--prefix TAG] <GM1|GM3|FSU2> <B1> <B2> ...",
            args[0]
        );
        eprintln!("Exemplo: {} --prefix novo FSU2 1e16 5e17 1e18", args[0]);
        return;
    }

    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        "FSU2" => FSU2,
        _ => panic!("Modelo não reconhecido. Use GM1, GM3 ou FSU2."),
    };

    let b_fields: Vec<f64> = args[2..]
        .iter()
        .map(|s| s.parse().expect("B deve ser um número válido"))
        .collect();

    fs::create_dir_all("results/magtop").unwrap_or_default();
    let plots_dir = format!("results/magtop/{}", model_name);
    fs::create_dir_all(&plots_dir).unwrap_or_default();

    let csv_summary = format!(
        "results/magtop/{}summary_topology_{}.csv",
        file_prefix, model_name
    );
    let mut summary_data = Vec::new();

    for &b_field in &b_fields {
        let b_string = format!("{:.2e}", b_field);

        // NOVO: prefixo também na pasta de saída
        let base_dir = format!("output/magtop/{}/{}B_{}", model_name, file_prefix, b_string);
        let dir_iso = format!("{}/isotropic", base_dir);
        let dir_aniso = format!("{}/anisotropic", base_dir);
        let plot_b_tag = b_string.replace('+', "p").replace('-', "m");

        let mut artist_eos = Artist::new(
            &format!(
                "{}/{}eos_topology_{}_B_{}.svg",
                plots_dir, file_prefix, model_name, plot_b_tag
            ),
            &format!(
                "Equation of State (Topology) - {} | B = {} G",
                model_name, b_string
            ),
        )
        .with_x_label("Energy Density \u{03B5} [MeV/fm\u{00B3}]")
        .with_y_label("Pressure P [MeV/fm\u{00B3}]")
        .autoscale()
        .with_log_scale();

        let mut artist_mr = Artist::new(
            &format!(
                "{}/{}mr_topology_{}_B_{}.svg",
                plots_dir, file_prefix, model_name, plot_b_tag
            ),
            &format!(
                "Mass-Radius (Topology) - {} | B = {} G",
                model_name, b_string
            ),
        )
        .with_x_label("Radius [km]")
        .with_y_label("Mass [M\u{2299}]")
        .with_x_range(8.0, 16.0);

        let mut any_curve_plotted = false;

        // NOVO: artistas de população por topologia
        let mut artist_pop_iso = Artist::new(
            &format!(
                "{}/{}pop_iso_{}_B_{}.svg",
                plots_dir, file_prefix, model_name, plot_b_tag
            ),
            &format!(
                "Particle Population (Iso) - {} | B = {} G",
                model_name, b_string
            ),
        )
        .with_x_label("n / n0")
        .with_y_label("n_i [fm^-3]")
        .autoscale();

        let mut artist_pop_aniso = Artist::new(
            &format!(
                "{}/{}pop_aniso_{}_B_{}.svg",
                plots_dir, file_prefix, model_name, plot_b_tag
            ),
            &format!(
                "Particle Population (Aniso) - {} | B = {} G",
                model_name, b_string
            ),
        )
        .with_x_label("n / n0")
        .with_y_label("n_i [fm^-3]")
        .autoscale();

        let mut any_pop_iso = false;
        let mut any_pop_aniso = false;

        if !plot_only {
            println!("Calculando EoS para B = {} G...", b_string);

            // Instancia os dois motores com as topologias diferentes.
            let engine_iso = EngineMode::Hadrons(
                HadronsMatter::new(model_params, b_field)
                    .with_topology(MagneticTopology::Isotropic)
                    .with_limits(0.01, 2.2)
                    .with_points(1500),
            );

            let engine_aniso = EngineMode::Hadrons(
                HadronsMatter::new(model_params, b_field)
                    .with_topology(MagneticTopology::Anisotropic)
                    .with_limits(0.01, 2.2)
                    .with_points(1500),
            );

            // Resolve as duas em paralelo
            let results = Solver::solve_parallel(vec![engine_iso, engine_aniso], 8);

            if results.len() == 2 {
                fs::create_dir_all(&dir_iso).ok();
                Solver::write_eos(&results[0], &format!("{}/eos.dat", dir_iso)).ok();

                fs::create_dir_all(&dir_aniso).ok();
                Solver::write_eos(&results[1], &format!("{}/eos.dat", dir_aniso)).ok();
            } else {
                eprintln!("Aviso: resultados incompletos para B = {}", b_string);
            }
        }

        let topologies = [("isotropic", "Iso"), ("anisotropic", "Aniso")];

        for (folder, label_suffix) in topologies.iter() {
            let eos_path = format!("{}/{}/eos.dat", base_dir, folder);

            if let Ok(file) = fs::File::open(&eos_path) {
                let reader = BufReader::new(file);
                let mut eps_vec = Vec::new();
                let mut p_vec = Vec::new();
                let mut rho_vec = Vec::new();

                // NOVO: vetores de população
                let mut n_over_n0 = Vec::new();
                let mut n_e = Vec::new();
                let mut n_mu = Vec::new();
                let mut n_n = Vec::new();
                let mut n_p = Vec::new();
                let mut n_l0 = Vec::new();
                let mut n_sm = Vec::new();
                let mut n_s0 = Vec::new();
                let mut n_sp = Vec::new();
                let mut n_xm = Vec::new();
                let mut n_x0 = Vec::new();

                for line in reader.lines().flatten() {
                    let parts: Vec<f64> = line
                        .split_whitespace()
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect();

                    // EOS
                    if parts.len() > 2 {
                        let eps = parts[1];
                        let p = parts[2];
                        let rho = parts.get(0).copied().unwrap_or(0.0);
                        if eps > 0.0 && p > 0.0 {
                            eps_vec.push(eps);
                            p_vec.push(p);
                            rho_vec.push(rho);
                        }
                    }

                    // POPULAÇÃO: colunas 0 e 3..12
                    if parts.len() > 12 {
                        n_over_n0.push(parts[0]);
                        n_e.push(parts[3]);
                        n_mu.push(parts[4]);
                        n_n.push(parts[5]);
                        n_p.push(parts[6]);
                        n_l0.push(parts[7]);
                        n_sm.push(parts[8]);
                        n_s0.push(parts[9]);
                        n_sp.push(parts[10]);
                        n_xm.push(parts[11]);
                        n_x0.push(parts[12]);
                    }
                }

                if eps_vec.len() > 10 {
                    let (masses, radii, _b_masses, _central_p_list) =
                        generate_mr_curve(&eps_vec, &p_vec, &rho_vec, false);
                    let label = format!("{} ({})", b_string, label_suffix);

                    artist_eos = artist_eos.add_curve(&eps_vec, &p_vec, &label);
                    artist_mr = artist_mr.add_curve(&radii, &masses, &label);
                    any_curve_plotted = true;

                    let mut m_max = 0.0;
                    let mut r_at_m_max = 0.0;
                    for i in 0..masses.len() {
                        if masses[i] > m_max && radii[i] < 16.0 {
                            m_max = masses[i];
                            r_at_m_max = radii[i];
                        }
                    }
                    summary_data.push(format!(
                        "{},{},{:.4},{:.4}",
                        b_string, label_suffix, m_max, r_at_m_max
                    ));
                }

                // NOVO: adiciona curvas de população
                if n_over_n0.len() > 10 {
                    if *folder == "isotropic" {
                        artist_pop_iso = artist_pop_iso
                            .add_curve(&n_over_n0, &n_e, "e-")
                            .add_curve(&n_over_n0, &n_mu, "mu-")
                            .add_curve(&n_over_n0, &n_n, "n")
                            .add_curve(&n_over_n0, &n_p, "p")
                            .add_curve(&n_over_n0, &n_l0, "L0")
                            .add_curve(&n_over_n0, &n_sm, "S-")
                            .add_curve(&n_over_n0, &n_s0, "S0")
                            .add_curve(&n_over_n0, &n_sp, "S+")
                            .add_curve(&n_over_n0, &n_xm, "X-")
                            .add_curve(&n_over_n0, &n_x0, "X0");
                        any_pop_iso = true;
                    } else {
                        artist_pop_aniso = artist_pop_aniso
                            .add_curve(&n_over_n0, &n_e, "e-")
                            .add_curve(&n_over_n0, &n_mu, "mu-")
                            .add_curve(&n_over_n0, &n_n, "n")
                            .add_curve(&n_over_n0, &n_p, "p")
                            .add_curve(&n_over_n0, &n_l0, "L0")
                            .add_curve(&n_over_n0, &n_sm, "S-")
                            .add_curve(&n_over_n0, &n_s0, "S0")
                            .add_curve(&n_over_n0, &n_sp, "S+")
                            .add_curve(&n_over_n0, &n_xm, "X-")
                            .add_curve(&n_over_n0, &n_x0, "X0");
                        any_pop_aniso = true;
                    }
                }
            }
        }

        if any_curve_plotted {
            artist_eos.plot().ok();
            artist_mr.plot().ok();
        } else {
            eprintln!(
                "Aviso: nenhum dado válido encontrado para B = {} G",
                b_string
            );
        }

        // NOVO: salva plots de população
        if any_pop_iso {
            artist_pop_iso.plot().ok();
        }
        if any_pop_aniso {
            artist_pop_aniso.plot().ok();
        }
    }

    if let Ok(mut file) = fs::File::create(&csv_summary) {
        use std::io::Write;
        writeln!(file, "b_field,topology,max_mass,radius_at_max").ok();
        for line in summary_data {
            writeln!(file, "{}", line).ok();
        }
    }

    println!("\nAnálise concluída! Gráficos salvos em: results/magtop/");
}
