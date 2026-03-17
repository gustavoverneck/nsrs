// src/bin/lsv.rs

use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::HadronsMatter;
use nsrs::core::solver::{Solver, EngineMode}; 
use nsrs::core::plotting::{Artist, ColorScale, Palette}; 
use nsrs::core::tov_solver::{generate_mr_curve, smooth_array}; 
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

    if args.len() < 6 {
        eprintln!("Uso: {} [--plot-only] <GM1|GM3> <exp_min> <exp_max> <num_pontos> <B1> <B2> ...", args[0]);
        eprintln!("Exemplo (\u{03BE} em MeV\u{207B}\u{00B9}): {} --plot-only GM1 -5 2 800 1e16 1e18", args[0]);
        return;
    }
    
    let model_name = &args[1];
    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let exp_min: f64 = args[2].parse().expect("exp_min inválido");
    let exp_max: f64 = args[3].parse().expect("exp_max inválido");
    let num_points: usize = args[4].parse().expect("num_pontos inválido");
    
    let b_fields: Vec<f64> = args[5..]
        .iter()
        .map(|s| s.parse().expect("B deve ser um número válido"))
        .collect();

    fs::create_dir_all("results/lsv").unwrap_or_default();

    let mut log_csi_vals = Vec::with_capacity(num_points);
    let mut csi_vals = Vec::with_capacity(num_points);

    if num_points == 1 {
        csi_vals.push(10_f64.powf(exp_min));
        log_csi_vals.push(exp_min);
    } else {
        let step = (exp_max - exp_min) / ((num_points - 1) as f64);
        for i in 0..num_points {
            let log_val = exp_min + (i as f64) * step;
            log_csi_vals.push(log_val);
            csi_vals.push(10_f64.powf(log_val));
        }
    }

    let x_ticks_count = (exp_max - exp_min).abs().round() as usize;

    // Artists para LIMITES GLOBAIS
    let mut artist_m = Artist::new(
        &format!("results/lsv/mass_max_{}_combined.svg", model_name),
        &format!("Maximum Mass Limit (LSV) - {}", model_name),
    ).with_x_label("log\u{2081}\u{2080}(\u{03BE}) [MeV\u{207B}\u{00B9}]").with_y_label("Maximum Mass [M\u{2299}]").autoscale().with_x_range(exp_min, exp_max).with_x_labels(x_ticks_count).with_y_range(1.7, 2.1);

    let mut artist_r = Artist::new(
        &format!("results/lsv/radius_max_{}_combined.svg", model_name),
        &format!("Radius at Maximum Mass (LSV) - {}", model_name),
    ).with_x_label("log\u{2081}\u{2080}(\u{03BE}) [MeV\u{207B}\u{00B9}]").with_y_label("Radius [km]").autoscale().with_x_range(exp_min, exp_max).with_x_labels(x_ticks_count).with_y_range(11.5, 12.0);

    let mut artist_nc = Artist::new(
        &format!("results/lsv/nc_central_{}_combined.svg", model_name),
        &format!("Central Baryon Density (LSV) - {}", model_name),
    ).with_x_label("log\u{2081}\u{2080}(\u{03BE}) [MeV\u{207B}\u{00B9}]").with_y_label("Central Density [n_c / n\u{2080}]").autoscale().with_x_range(exp_min, exp_max).with_x_labels(x_ticks_count);

    let mut artist_meff_limit = Artist::new(
        &format!("results/lsv/meff_central_{}_combined.svg", model_name),
        &format!("Central Effective Mass (LSV) - {}", model_name),
    ).with_x_label("log\u{2081}\u{2080}(\u{03BE}) [MeV\u{207B}\u{00B9}]").with_y_label("m* / m_N").autoscale().with_x_range(exp_min, exp_max).with_x_labels(x_ticks_count);

    let color_scale = ColorScale::new(exp_min, exp_max, Palette::Plasma);

    // Closure de Segurança para Frações (Evita o Log(0))
    let clamp = |val: f64| -> f64 { val.max(1e-5) };

    for &b_field in &b_fields {
        let b_string = format!("{:.2e}", b_field);
        
        let mut max_masses = Vec::new();
        let mut radii_at_max = Vec::new();
        let mut valid_log_csi = Vec::new();
        let mut central_nc_vals = Vec::new();
        let mut central_meff_vals = Vec::new();
        
        let csv_path = format!("results/lsv/data_{}_B_{}.csv", model_name, b_string);
        let base_dir = format!("output/lsv/{}/B_{}", model_name, b_string);

        // Prepara os Artists de CORES
        let mut artist_mr_color = Artist::new(&format!("results/lsv/color_mr_{}_B_{}.svg", model_name, b_string), &format!("MR - Varying LSV (\u{03BE}) - B={} G", b_string)).with_x_label("Radius [km]").with_y_label("Mass [M\u{2299}]").with_x_range(10.0, 15.0).with_y_range(0.0, 2.2);
        let mut artist_eos_color = Artist::new(&format!("results/lsv/color_eos_{}_B_{}.svg", model_name, b_string), &format!("EoS - Varying LSV (\u{03BE}) - B={} G", b_string)).with_x_label("\u{03B5} [MeV/fm\u{00B3}]").with_y_label("P [MeV/fm\u{00B3}]").autoscale().with_log_scale();
        let mut artist_meff_color = Artist::new(&format!("results/lsv/color_meff_{}_B_{}.svg", model_name, b_string), &format!("Effective Mass vs Density (\u{03BE}) - B={} G", b_string)).with_x_label("Density [n / n\u{2080}]").with_y_label("m* / m_N").autoscale();
        
        // CORREÇÃO: Escala Logarítmica Adicionada nos Gráficos de Cor Populacionais
        let mut artist_pop_p_color = Artist::new(&format!("results/lsv/color_pop_p_{}_B_{}.svg", model_name, b_string), &format!("Proton Fraction Y_p (\u{03BE}) - B={} G", b_string)).with_x_label("Density [n / n\u{2080}]").with_y_label("Y_p").autoscale().with_log_scale().with_y_range(1e-5, 1.2);
        let mut artist_pop_sm_color = Artist::new(&format!("results/lsv/color_pop_sm_{}_B_{}.svg", model_name, b_string), &format!("\u{03A3}\u{207B} Fraction Y_\u{03A3} (\u{03BE}) - B={} G", b_string)).with_x_label("Density [n / n\u{2080}]").with_y_label("Y_\u{03A3}").autoscale().with_log_scale().with_y_range(1e-5, 1.2);

        if plot_only {
            println!("Lendo dados salvos de {}...", csv_path);
            if let Ok(file) = fs::File::open(&csv_path) {
                let reader = BufReader::new(file);
                for line in reader.lines().skip(1) {
                    if let Ok(content) = line {
                        let parts: Vec<&str> = content.split(',').collect();
                        if parts.len() >= 6 {
                            if let (Ok(log_csi), Ok(m_max), Ok(r_max), Ok(nc), Ok(meff)) = (
                                parts[1].parse::<f64>(), parts[2].parse::<f64>(), parts[3].parse::<f64>(),
                                parts[4].parse::<f64>(), parts[5].parse::<f64>()
                            ) {
                                valid_log_csi.push(log_csi);
                                max_masses.push(m_max);
                                radii_at_max.push(r_max);
                                central_nc_vals.push(nc);
                                central_meff_vals.push(meff);
                            }
                        }
                    }
                }
            } else {
                eprintln!("  -> ERRO: Arquivo {} não encontrado.", csv_path);
                continue;
            }

            for i in 0..csi_vals.len() {
                let csi = csi_vals[i];
                let log_csi = log_csi_vals[i];
                let eos_filename = format!("{}/csi_{:.2e}/eos.dat", base_dir, csi);
                
                let mut n_n0_f = Vec::new(); let mut eps_f = Vec::new(); let mut p_f = Vec::new(); let mut meff_f = Vec::new();
                let mut yp_f = Vec::new(); let mut ysm_f = Vec::new(); let mut yn_f = Vec::new();
                let mut ye_f = Vec::new(); let mut ymu_f = Vec::new(); let mut yl0_f = Vec::new(); let mut yxm_f = Vec::new();
                let mut ys0_f = Vec::new(); let mut ysp_f = Vec::new(); let mut yx0_f = Vec::new();

                if let Ok(file) = fs::File::open(&eos_filename) {
                    let reader = BufReader::new(file);
                    for line in reader.lines().skip(1) {
                        if let Ok(content) = line {
                            let parts: Vec<&str> = content.split_whitespace().collect();
                            if parts.len() >= 17 {
                                let parse = |idx: usize| parts[idx].parse::<f64>().unwrap_or(0.0);
                                let n = parse(0); let e = parse(1); let p = parse(2);
                                let ne = parse(3); let nmu = parse(4); let nn = parse(5); let np = parse(6);
                                let nl0 = parse(7); let nsm = parse(8); let ns0 = parse(9); let nsp = parse(10);
                                let nxm = parse(11); let nx0 = parse(12); let m = parse(16);

                                if e > 0.0 && p > 0.0 {
                                    n_n0_f.push(n); eps_f.push(e); p_f.push(p); meff_f.push(m);
                                    
                                    let nb_tot = nn + np + nl0 + nsm + ns0 + nsp + nxm + nx0;
                                    if nb_tot > 0.0 {
                                        yn_f.push(clamp(nn / nb_tot));
                                        yp_f.push(clamp(np / nb_tot));
                                        ye_f.push(clamp(ne / nb_tot));
                                        ymu_f.push(clamp(nmu / nb_tot));
                                        yl0_f.push(clamp(nl0 / nb_tot));
                                        ysm_f.push(clamp(nsm / nb_tot));
                                        ys0_f.push(clamp(ns0 / nb_tot));
                                        ysp_f.push(clamp(nsp / nb_tot));
                                        yxm_f.push(clamp(nxm / nb_tot));
                                        yx0_f.push(clamp(nx0 / nb_tot));
                                    }
                                }
                            }
                        }
                    }
                    if eps_f.len() > 10 {
                        let window_size = 30; 
                        // let eps_smooth = smooth_array(&eps_f, window_size);
                        // let p_smooth = smooth_array(&p_f, window_size);

                        // Passa os vetores suavizados para o solver TOV
                        let (masses, radii, central_p_list) = generate_mr_curve(&eps_f, &p_f, false);
                        let label = format!("log(\u{03BE}) = {:.2}", log_csi);
                        let (r, g, b) = color_scale.get_color(log_csi);
                        
                        // Altere para plotar todas as curvas, mas apenas algumas na legenda
                        if i % 40 == 0 || i == csi_vals.len() - 1 {
                            artist_mr_color = artist_mr_color.add_curve_color(&radii, &masses, &label, r, g, b);
                        } else {
                            // Se o seu Artist suportar, adicione a curva sem label para não poluir a legenda
                            artist_mr_color = artist_mr_color.add_curve_color(&radii, &masses, "", r, g, b);
                        }
                        
                        artist_mr_color = artist_mr_color.add_curve_color(&radii, &masses, &label, r, g, b);
                        artist_eos_color = artist_eos_color.add_curve_color(&eps_f, &p_f, &label, r, g, b);
                        artist_meff_color = artist_meff_color.add_curve_color(&n_n0_f, &meff_f, &label, r, g, b);
                        artist_pop_p_color = artist_pop_p_color.add_curve_color(&n_n0_f, &yp_f, &label, r, g, b);
                        artist_pop_sm_color = artist_pop_sm_color.add_curve_color(&n_n0_f, &ysm_f, &label, r, g, b);

                        // Gráfico Clássico de Populações
                        if i == 0 || i == csi_vals.len() - 1 {
                            let tag = if i == 0 { "min" } else { "max" };
                            let mut artist_pop_all = Artist::new(
                                &format!("results/lsv/pop_all_{}_B_{}_{}_csi.svg", model_name, b_string, tag),
                                &format!("Particle Fractions (\u{03BE} = {:.1e}) - B={} G", csi, b_string),
                            ).with_x_label("Density [n / n\u{2080}]").with_y_label("Fraction Y_i (Log)").autoscale().with_log_scale().with_y_range(1e-5, 1.5);

                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yn_f, "n");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yp_f, "p");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ye_f, "e\u{207B}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ymu_f, "\u{03BC}\u{207B}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yl0_f, "\u{039B}\u{2070}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ysm_f, "\u{03A3}\u{207B}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ys0_f, "\u{03A3}\u{2070}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ysp_f, "\u{03A3}\u{207A}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yxm_f, "\u{039E}\u{207B}");
                            artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yx0_f, "\u{039E}\u{2070}");
                            artist_pop_all.plot().ok();
                        }
                    }
                }
            }

        } else {
            let engines: Vec<EngineMode> = csi_vals.iter().map(|&csi| {
                let motor = HadronsMatter::new(model_params, b_field).with_csi(csi).with_limits(0.01, 2.5).with_points(2000);
                EngineMode::Hadrons(motor)
            }).collect();

            println!("\nVarrendo {} valores de \u{03BE} para B = {} G...", num_points, b_string);
            let all_results = Solver::solve_batch(engines);

            for (i, results) in all_results.iter().enumerate() {
                let csi = csi_vals[i];
                let log_csi = log_csi_vals[i];
                
                let dir_path = format!("{}/csi_{:.2e}", base_dir, csi);
                if let Err(_) = fs::create_dir_all(&dir_path) { continue; }
                
                let eos_filename = format!("{}/eos.dat", dir_path);
                if let Err(_) = Solver::write_eos(results, &eos_filename) { continue; }

                let mut n_n0_f = Vec::new(); let mut eps_f = Vec::new(); let mut p_f = Vec::new(); let mut meff_f = Vec::new();
                let mut yp_f = Vec::new(); let mut ysm_f = Vec::new(); let mut yn_f = Vec::new();
                let mut ye_f = Vec::new(); let mut ymu_f = Vec::new(); let mut yl0_f = Vec::new(); let mut yxm_f = Vec::new();
                let mut ys0_f = Vec::new(); let mut ysp_f = Vec::new(); let mut yx0_f = Vec::new();

                for res in results.iter() {
                    if res[1] > 0.0 && res[2] > 0.0 {
                        n_n0_f.push(res[0]); eps_f.push(res[1]); p_f.push(res[2]); meff_f.push(res[16]);
                        
                        let nb_tot = res[5] + res[6] + res[7] + res[8] + res[9] + res[10] + res[11] + res[12];
                        if nb_tot > 0.0 {
                            yn_f.push(clamp(res[5] / nb_tot));
                            yp_f.push(clamp(res[6] / nb_tot));
                            ye_f.push(clamp(res[3] / nb_tot));
                            ymu_f.push(clamp(res[4] / nb_tot));
                            yl0_f.push(clamp(res[7] / nb_tot));
                            ysm_f.push(clamp(res[8] / nb_tot));
                            ys0_f.push(clamp(res[9] / nb_tot));
                            ysp_f.push(clamp(res[10] / nb_tot));
                            yxm_f.push(clamp(res[11] / nb_tot));
                            yx0_f.push(clamp(res[12] / nb_tot));
                        }
                    }
                }
                
                if eps_f.len() > 10 {
                    let (masses, radii, central_p_list) = generate_mr_curve(&eps_f, &p_f, false);
                    let label = format!("log(\u{03BE}) = {:.2}", log_csi);
                    let (cr, cg, cb) = color_scale.get_color(log_csi);
                    
                    artist_mr_color = artist_mr_color.add_curve_color(&radii, &masses, &label, cr, cg, cb);
                    artist_eos_color = artist_eos_color.add_curve_color(&eps_f, &p_f, &label, cr, cg, cb);
                    artist_meff_color = artist_meff_color.add_curve_color(&n_n0_f, &meff_f, &label, cr, cg, cb);
                    artist_pop_p_color = artist_pop_p_color.add_curve_color(&n_n0_f, &yp_f, &label, cr, cg, cb);
                    artist_pop_sm_color = artist_pop_sm_color.add_curve_color(&n_n0_f, &ysm_f, &label, cr, cg, cb);

                    if i == 0 || i == csi_vals.len() - 1 {
                        let tag = if i == 0 { "min" } else { "max" };
                        let mut artist_pop_all = Artist::new(
                            &format!("results/lsv/pop_all_{}_B_{}_{}_csi.svg", model_name, b_string, tag),
                            &format!("Particle Fractions (\u{03BE} = {:.1e}) - B={} G", csi, b_string),
                        ).with_x_label("Density [n / n\u{2080}]").with_y_label("Fraction Y_i (Log)").autoscale().with_log_scale().with_y_range(1e-5, 1.5);

                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yn_f, "n");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yp_f, "p");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ye_f, "e\u{207B}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ymu_f, "\u{03BC}\u{207B}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yl0_f, "\u{039B}\u{2070}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ysm_f, "\u{03A3}\u{207B}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ys0_f, "\u{03A3}\u{2070}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &ysp_f, "\u{03A3}\u{207A}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yxm_f, "\u{039E}\u{207B}");
                        artist_pop_all = artist_pop_all.add_curve(&n_n0_f, &yx0_f, "\u{039E}\u{2070}");
                        artist_pop_all.plot().ok();
                    }

                    let mut m_max = 0.0; let mut r_max_m = 0.0; let mut best_pc = 0.0; 
                    
                    for idx_star in 0..masses.len() {
                        let r = radii[idx_star]; let m = masses[idx_star]; let pc = central_p_list[idx_star];
                        if r < 15.0 && m > m_max { m_max = m; r_max_m = r; best_pc = pc; }
                    }
                    
                    if m_max > 0.0 {
                        let core_idx = results.iter().position(|r| (r[2] - best_pc).abs() < 1e-8).unwrap_or(results.len() - 1);
                        max_masses.push(m_max); radii_at_max.push(r_max_m); valid_log_csi.push(log_csi);
                        central_nc_vals.push(results[core_idx][0]); central_meff_vals.push(results[core_idx][16]); 
                    }
                }
            }

            if let Ok(mut file) = fs::File::create(&csv_path) {
                use std::io::Write;
                writeln!(file, "csi,log10_csi,max_mass_msun,radius_at_max_km,central_nc,central_meff").ok();
                for i in 0..valid_log_csi.len() {
                    writeln!(file, "{:e},{:.4},{:.4},{:.4},{:.4},{:.4}", 
                        10_f64.powf(valid_log_csi[i]), valid_log_csi[i], max_masses[i], radii_at_max[i], 
                        central_nc_vals[i], central_meff_vals[i]
                    ).ok();
                }
            }
        }

        artist_mr_color.plot().ok();
        artist_eos_color.plot().ok();
        artist_meff_color.plot().ok();
        artist_pop_p_color.plot().ok();
        artist_pop_sm_color.plot().ok();
        println!("  -> Gráficos de degradê e populações gerados em 'results/lsv/' para B = {} G.", b_string);

        let label = format!("B = {} G", b_string);
        artist_m = artist_m.add_curve(&valid_log_csi, &max_masses, &label);
        artist_r = artist_r.add_curve(&valid_log_csi, &radii_at_max, &label);
        artist_nc = artist_nc.add_curve(&valid_log_csi, &central_nc_vals, &label);
        artist_meff_limit = artist_meff_limit.add_curve(&valid_log_csi, &central_meff_vals, &label);
    }

    artist_m.plot().ok();
    artist_r.plot().ok();
    artist_nc.plot().ok();
    artist_meff_limit.plot().ok();

    println!("\nProcesso concluído com sucesso!");
}