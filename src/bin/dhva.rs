// src/bin/dhva.rs

use nsrs::core::constants::{ME, M_NUCLEON, QE, BCE};
use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::HadronsMatter;
use nsrs::core::solver::{Solver, EngineMode};
use nsrs::core::plotting::{Artist, ColorScale, Palette};

use std::env;

fn calculate_nu_max(mu_e_norm: f64, b_gauss: f64, xi: f64) -> f64 {
    let m = ME / M_NUCLEON;
    let b0 = b_gauss / 4.41e13;
    let b = b0 * BCE;

    let lsv_shift = 2.0 * xi * m * b + (xi * b).powi(2);

    if mu_e_norm.powi(2) > (m.powi(2) + lsv_shift) && b > 0.0 {
        ((mu_e_norm.powi(2) - m.powi(2) - lsv_shift) / (2.0 * QE * b)).floor()
    } else {
        0.0
    }
}

fn calculate_frequency(mu_e_norm: f64) -> f64 {
    let mu = mu_e_norm * M_NUCLEON;
    if mu > ME {
        mu.powi(2) - ME.powi(2)
    } else {
        0.0
    }
}

fn derivative(x: &Vec<f64>, y: &Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    if x.len() < 2 {
        return (Vec::new(), Vec::new());
    }

    let mut dx = Vec::new();
    let mut dy = Vec::new();

    for i in 1..x.len() {
        let ddx = x[i] - x[i - 1];
        if ddx != 0.0 {
            dx.push(x[i]);
            dy.push((y[i] - y[i - 1]) / ddx);
        }
    }

    (dx, dy)
}

fn find_closest_density(data: &Vec<[f64; 20]>, target: f64) -> Option<[f64; 20]> {
    data.iter()
        .min_by(|a, b| {
            (a[0] - target)
                .abs()
                .partial_cmp(&(b[0] - target).abs())
                .unwrap()
        })
        .copied()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 9 {
        eprintln!("Uso: {} <GM1|GM3> <n/n0> <log_xi_min> <log_xi_max> <num_xi> <Bmin> <Bmax> <num_B>", args[0]);
        return;
    }

    let model_name = &args[1];
    let target_density: f64 = args[2].parse().unwrap();

    let exp_min: f64 = args[3].parse().unwrap();
    let exp_max: f64 = args[4].parse().unwrap();
    let num_xi: usize = args[5].parse().unwrap();

    let b_min: f64 = args[6].parse().unwrap();
    let b_max: f64 = args[7].parse().unwrap();
    let num_b: usize = args[8].parse().unwrap();

    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo inválido"),
    };

    let xi_vals: Vec<f64> = (0..num_xi)
        .map(|i| {
            let log_val = exp_min
                + (i as f64) * (exp_max - exp_min)
                    / (num_xi.saturating_sub(1) as f64).max(1.0);
            10_f64.powf(log_val)
        })
        .collect();

    let b_vals: Vec<f64> = (0..num_b)
        .map(|i| {
            b_min
                + (i as f64) * (b_max - b_min)
                    / (num_b.saturating_sub(1) as f64).max(1.0)
        })
        .collect();

    let color_scale = ColorScale::new(exp_min, exp_max, Palette::Plasma);

    // ===== ARTISTS =====
    let mut a_nu_invb = Artist::new("results/dhva_scan/nu_vs_invB.svg","ν vs 1/B")
        .with_x_label("1/B").with_y_label("ν_max").autoscale();

    let mut a_freq_invb = Artist::new("results/dhva_scan/freq_vs_invB.svg","F vs 1/B")
        .with_x_label("1/B").with_y_label("F").autoscale();

    let mut a_osc = Artist::new("results/dhva_scan/osc_vs_invB.svg","dF/d(1/B)")
        .with_x_label("1/B").with_y_label("dF/d(1/B)").autoscale();

    let mut a_nu_b = Artist::new("results/dhva_scan/nu_vs_B.svg","ν vs B")
        .with_x_label("B").with_y_label("ν_max").autoscale();

    let mut a_freq_b = Artist::new("results/dhva_scan/freq_vs_B.svg","F vs B")
        .with_x_label("B").with_y_label("F").autoscale();

    let mut a_df_db = Artist::new("results/dhva_scan/df_db.svg","dF/dB")
        .with_x_label("B").with_y_label("dF/dB").autoscale();

    let mut a_dnu = Artist::new("results/dhva_scan/dnu.svg","Δν (transições)")
        .with_x_label("1/B").with_y_label("Δν").autoscale();

    // ===== LOOP =====
    for &xi in &xi_vals {
        let mut x_invb = Vec::new();
        let mut x_b = Vec::new();
        let mut y_nu = Vec::new();
        let mut y_f = Vec::new();

        for &b in &b_vals {
            let engine = EngineMode::Hadrons(
                HadronsMatter::new(model_params, b)
                    .with_csi(xi)
                    .with_limits(0.01, 3.0)
                    .with_points(800),
            );

            let mut solver = Solver::new(engine);
            let result = solver.solve();

            if let Some(row) = find_closest_density(&result, target_density) {
                let mu_e = row[18];

                x_invb.push(1.0 / b);
                x_b.push(b);
                y_nu.push(calculate_nu_max(mu_e, b, xi));
                y_f.push(calculate_frequency(mu_e));
            }
        }

        let (x_dinvb, y_dinvb) = derivative(&x_invb, &y_f);
        let (x_db, y_db) = derivative(&x_b, &y_f);
        let (_, dnu) = derivative(&x_invb, &y_nu);

        let (r,g,b) = color_scale.get_color(xi.log10());
        let label = format!("log ξ={:.1}", xi.log10());

        a_nu_invb = a_nu_invb.add_curve_color(&x_invb,&y_nu,&label,r,g,b);
        a_freq_invb = a_freq_invb.add_curve_color(&x_invb,&y_f,&label,r,g,b);

        if !x_dinvb.is_empty() {
            a_osc = a_osc.add_curve_color(&x_dinvb,&y_dinvb,&label,r,g,b);
            a_dnu = a_dnu.add_curve_color(&x_dinvb,&dnu,&label,r,g,b);
        }

        a_nu_b = a_nu_b.add_curve_color(&x_b,&y_nu,&label,r,g,b);
        a_freq_b = a_freq_b.add_curve_color(&x_b,&y_f,&label,r,g,b);

        if !x_db.is_empty() {
            a_df_db = a_df_db.add_curve_color(&x_db,&y_db,&label,r,g,b);
        }
    }

    a_nu_invb.plot().ok();
    a_freq_invb.plot().ok();
    a_osc.plot().ok();
    a_nu_b.plot().ok();
    a_freq_b.plot().ok();
    a_df_db.plot().ok();
    a_dnu.plot().ok();

    println!("✅ Suite completa dHvA gerada em results/dhva_scan/");
}