// Ferramentas sobre arquivos de EoS (substitui os binários tov e bdd).

use std::fs;

use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::plot_mr_curve;
use nsrs::core::tov_solver::generate_mr_curve;
use plotters::prelude::*;

use crate::cli::create_dir;

/// `tov <eos.dat>`: curva M-R de uma EoS em arquivo -> results/mr_<nome>.svg
pub fn tov(args: &[String]) -> Result<(), String> {
    let [eos_path] = args else {
        return Err("uso: nsrs tov <eos.dat>".into());
    };
    println!("Lendo arquivo EOS: {eos_path}");
    let (eps, p, rho) =
        read_eos_file(eos_path).map_err(|e| format!("falha ao ler '{eos_path}': {e}"))?;
    println!("EoS lida com sucesso. Pontos carregados: {}", eps.len());
    println!("Integrando equações TOV...");
    let (masses, radii, _b_masses, _) = generate_mr_curve(&eps, &p, &rho, false);
    if masses.is_empty() || radii.is_empty() {
        return Err(
            "a curva M-R retornou vazia; verifique se as pressões da EoS suportam uma estrela".into(),
        );
    }
    let output_name = std::path::Path::new(eos_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("eos");
    create_dir("results")?;
    let plot_path = format!("results/mr_{output_name}.svg");
    println!("Gerando gráfico em: {plot_path}");
    plot_mr_curve(&radii, &masses, &plot_path)
        .map_err(|e| format!("erro ao salvar o gráfico SVG: {e}"))?;
    let max_mass = masses.iter().copied().fold(f64::NAN, f64::max);
    println!("Massa máxima da estrela: {max_mass:.2} M_sol");
    Ok(())
}

/// `plot <eos.dat> <saida.png> <col_x> <col_y>`: uma coluna contra outra
/// (índices a partir de 0).
pub fn plot(args: &[String]) -> Result<(), String> {
    let [input, output, x_col, y_col] = args else {
        return Err("uso: nsrs plot <eos.dat> <saida.png> <col_x> <col_y>".into());
    };
    let column = |s: &str| s.parse::<usize>().map_err(|_| format!("coluna inválida: '{s}'"));
    let (x_col, y_col) = (column(x_col)?, column(y_col)?);
    let output = output.as_str();
    let content = fs::read_to_string(input).map_err(|e| format!("falha ao ler '{input}': {e}"))?;

    let data = parse_data(&content, x_col, y_col);
    if data.is_empty() {
        return Err("nenhum dado válido; confira os índices das colunas e o formato".into());
    }

    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for (x, y) in &data {
        x_min = x_min.min(*x);
        x_max = x_max.max(*x);
        y_min = y_min.min(*y);
        y_max = y_max.max(*y);
    }

    if (x_max - x_min).abs() < 1e-14 {
        x_min -= 1.0;
        x_max += 1.0;
    }
    if (y_max - y_min).abs() < 1e-14 {
        y_min -= 1.0;
        y_max += 1.0;
    }

    let root = BitMapBackend::new(output, (1024, 768)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| e.to_string())?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("EOS column {y_col} vs column {x_col}"),
            ("sans-serif", 32),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(|e| e.to_string())?;

    chart
        .configure_mesh()
        .x_desc(format!("EOS column {x_col}"))
        .y_desc(format!("EOS column {y_col}"))
        .y_label_formatter(&|v| format!("{:.3e}", v))
        .draw()
        .map_err(|e| e.to_string())?;

    chart
        .draw_series(LineSeries::new(data.clone(), &BLUE))
        .map_err(|e| e.to_string())?
        .label(format!("column {y_col}"))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 25, y)], BLUE));

    chart
        .draw_series(data.into_iter().map(|p| Circle::new(p, 2, BLUE.filled())))
        .map_err(|e| e.to_string())?;

    chart
        .configure_series_labels()
        .border_style(BLACK)
        .draw()
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())?;
    println!("Plot written to {}", output);
    Ok(())
}

fn parse_data(content: &str, x_col: usize, y_col: usize) -> Vec<(f64, f64)> {
    let mut data = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Accept whitespace-separated and tolerate commas
        let clean = line.replace(',', " ");
        let cols: Vec<&str> = clean.split_whitespace().collect();

        if cols.len() <= x_col || cols.len() <= y_col {
            continue;
        }

        let x = cols[x_col].parse::<f64>();
        let y = cols[y_col].parse::<f64>();

        if let (Ok(xv), Ok(yv)) = (x, y) {
            if xv.is_finite() && yv.is_finite() {
                data.push((xv, yv));
            }
        }
    }

    data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    data
}

