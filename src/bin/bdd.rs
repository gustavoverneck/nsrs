use plotters::prelude::*;
use std::{env, fs, process};

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

fn main() {
    // Usage:
    // cargo run --bin bdd -- <input.dat> [output.png] [x_col] [y_col]
    // Defaults: output=bdd_vs_density.png, x_col=0, y_col=20
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Usage: {} <input.dat> [output.png] [x_col] [y_col]",
            args[0]
        );
        process::exit(1);
    }

    let input = &args[1];
    let output = args
        .get(2)
        .map(String::as_str)
        .unwrap_or("bdd_vs_density.png");
    let x_col = args
        .get(3)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let y_col = args
        .get(4)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20);

    let content = fs::read_to_string(input).unwrap_or_else(|e| {
        eprintln!("Failed to read '{}': {}", input, e);
        process::exit(1);
    });

    let data = parse_data(&content, x_col, y_col);
    if data.is_empty() {
        eprintln!("No valid data found. Check column indexes and file format.");
        process::exit(1);
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
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("Bdd vs Density", ("sans-serif", 32))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("Density (n/n0)")
        .y_desc("Bdd (Tesla)")
        .y_label_formatter(&|v| format!("{:.3e}", v))
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(data.clone(), &BLUE))
        .unwrap()
        .label("Bdd")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 25, y)], BLUE));

    chart
        .draw_series(data.into_iter().map(|p| Circle::new(p, 2, BLUE.filled())))
        .unwrap();

    chart
        .configure_series_labels()
        .border_style(BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();
    println!("Plot written to {}", output);
}
