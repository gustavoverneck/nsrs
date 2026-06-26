// src/bin/nlem_log_limits.rs
// Find the limits of magnetic fields in neutron stars using nlem_log
// Tests increasingly absurd magnetic fields until physical limits are reached

use nsrs::{
    EngineMode, FSU2, GM1, GM3, HadronsMatter, NlemModel, Solver, core::model::ModelParams,
    generate_mr_curve, write_eos_with_mr,
};
use std::fs;

#[derive(Debug, Clone)]
struct FieldTestResult {
    b_field: f64,
    csi: f64,
    max_mass: Option<f64>,
    max_radius: Option<f64>,
    success: bool,
    error_msg: String,
}

fn compute_csi_for_effective_field(b_bare: f64, b_target: f64) -> Option<f64> {
    // For Log model: B_eff = B / (1 + B^2 / (2*csi^2))
    // Solving for csi: csi = sqrt(B_target * B^2 / (2 * (B - B_target)))

    if b_target >= b_bare || b_target <= 0.0 || b_bare <= 0.0 {
        return None;
    }

    let numerator = b_target * b_bare * b_bare;
    let denominator = 2.0 * (b_bare - b_target);

    if denominator <= 0.0 {
        return None;
    }

    let csi_sq = numerator / denominator;
    if csi_sq <= 0.0 {
        return None;
    }

    Some(csi_sq.sqrt())
}

fn generate_csi_values_for_field(b_bare: f64, num_steps: usize) -> Vec<f64> {
    // Generate csi values that produce effective fields from small to nearly b_bare
    let mut csi_vals = Vec::new();

    // Start from a small fraction of b_bare and go up to nearly b_bare (logarithmic scale)
    let b_min = b_bare * 0.01; // 1% of bare field
    let b_max = b_bare * 0.99; // 99% of bare field

    let log_min = b_min.log10();
    let log_max = b_max.log10();

    for i in 0..=num_steps {
        let t = i as f64 / num_steps as f64;
        let log_b_eff = log_min + t * (log_max - log_min);
        let b_eff = 10_f64.powf(log_b_eff);

        if let Some(csi) = compute_csi_for_effective_field(b_bare, b_eff) {
            csi_vals.push(csi);
        }
    }

    csi_vals
}

fn test_magnetic_field_no_nlem(
    model_params: ModelParams,
    b_field: f64,
    model_name: &str,
) -> FieldTestResult {
    let b_string = format!("{:.2e}", b_field);

    let engines: Vec<EngineMode> = vec![{
        let motor = HadronsMatter::new(model_params, b_field)
            .with_limits(0.01, 2.0)
            .with_points(2000);
        EngineMode::Hadrons(motor)
    }];

    // Safely attempt solver execution (use single thread to avoid GSL thread safety issues)
    let all_results = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        Solver::solve_parallel(engines, 1)
    })) {
        Ok(results) => results,
        Err(e) => {
            let err_msg = if let Some(s) = e.downcast_ref::<String>() {
                format!("Solver error: {}", s)
            } else if let Some(&s) = e.downcast_ref::<&str>() {
                format!("Solver error: {}", s)
            } else {
                "Solver panic or GSL error".to_string()
            };
            eprintln!("    ⚠ {} for B = {} G (no NLEM)", err_msg, b_string);
            return FieldTestResult {
                b_field,
                csi: 0.0,
                max_mass: None,
                max_radius: None,
                success: false,
                error_msg: err_msg,
            };
        }
    };

    let results_eos = &all_results[0];

    match (|| {
        if results_eos.is_empty() {
            return Err("Empty EOS");
        }

        let eps: Vec<f64> = results_eos.iter().map(|r| r[1]).collect();
        let p_arr: Vec<f64> = results_eos.iter().map(|r| r[2]).collect();

        if eps.is_empty() || p_arr.is_empty() {
            return Err("Empty eps or p_arr");
        }

        // Check for minimum data points required by GSL interpolation (akima needs at least 5 points)
        if eps.len() < 5 {
            return Err("Insufficient data points for interpolation");
        }

        let rho_arr: Vec<f64> = results_eos.iter().map(|r| r[0]).collect();
        let (masses, radii, b_masses, central_p_list) =
            generate_mr_curve(&eps, &p_arr, &rho_arr, false);

        if masses.is_empty() || radii.is_empty() {
            return Err("Empty M-R curve");
        }

        let max_mass = masses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !max_mass.is_finite() {
            return Err("Invalid max mass");
        }

        let max_radius = radii
            .iter()
            .zip(masses.iter())
            .filter(|(_, m)| (**m - max_mass).abs() < 0.01)
            .map(|(r, _)| *r)
            .next();

        // Try to save results
        let base_dir = format!(
            "output/nlem_log_limits/{}/B_{}/no_nlem",
            model_name, b_string
        );
        let _ = fs::create_dir_all(&base_dir);
        let eos_filename = format!("{}/eos.dat", base_dir);
        let _ = write_eos_with_mr(
            &results_eos,
            &masses,
            &radii,
            &b_masses,
            &central_p_list,
            &eos_filename,
        );

        Ok((max_mass, max_radius))
    })() {
        Ok((max_mass, max_radius)) => FieldTestResult {
            b_field,
            csi: 0.0,
            max_mass: Some(max_mass),
            max_radius,
            success: true,
            error_msg: String::new(),
        },
        Err(err_msg) => FieldTestResult {
            b_field,
            csi: 0.0,
            max_mass: None,
            max_radius: None,
            success: false,
            error_msg: err_msg.to_string(),
        },
    }
}

fn test_magnetic_field(
    model_params: ModelParams,
    b_field: f64,
    csi_vals: &[f64],
    model_name: &str,
) -> Vec<FieldTestResult> {
    let b_string = format!("{:.2e}", b_field);
    let mut results = Vec::new();

    let engines: Vec<EngineMode> = csi_vals
        .iter()
        .map(|&csi| {
            let motor = HadronsMatter::new(model_params, b_field)
                .with_nlem(NlemModel::Log(csi))
                .with_limits(0.01, 2.0)
                .with_points(2000);

            EngineMode::Hadrons(motor)
        })
        .collect();

    println!(
        "  Testing B = {} G with {} csi values...",
        b_string,
        csi_vals.len()
    );

    // Safely attempt solver execution (use single thread to avoid GSL thread safety issues)
    let all_results = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        Solver::solve_parallel(engines, 1)
    })) {
        Ok(results) => results,
        Err(e) => {
            let err_msg = if let Some(s) = e.downcast_ref::<String>() {
                format!("Solver error: {}", s)
            } else if let Some(&s) = e.downcast_ref::<&str>() {
                format!("Solver error: {}", s)
            } else {
                "Solver panic or GSL error".to_string()
            };
            eprintln!("    ⚠ {} for B = {} G", err_msg, b_string);
            return vec![FieldTestResult {
                b_field,
                csi: csi_vals.first().copied().unwrap_or(1e10),
                max_mass: None,
                max_radius: None,
                success: false,
                error_msg: err_msg,
            }];
        }
    };

    for (i, results_eos) in all_results.iter().enumerate() {
        let csi = csi_vals[i];

        // Try to extract M-R curve with comprehensive error handling
        match (|| {
            if results_eos.is_empty() {
                return Err("Empty EOS");
            }

            let eps: Vec<f64> = results_eos.iter().map(|r| r[1]).collect();
            let p_arr: Vec<f64> = results_eos.iter().map(|r| r[2]).collect();
            let rho_arr: Vec<f64> = results_eos.iter().map(|r| r[0]).collect();

            if eps.is_empty() || p_arr.is_empty() {
                return Err("Empty eps or p_arr");
            }

            // Check for minimum data points required by GSL interpolation (akima needs at least 5 points)
            if eps.len() < 5 {
                return Err("Insufficient data points for interpolation");
            }

            let (masses, radii, b_masses, central_p_list) =
                generate_mr_curve(&eps, &p_arr, &rho_arr, false);

            if masses.is_empty() || radii.is_empty() {
                return Err("Empty M-R curve");
            }

            let max_mass = masses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if !max_mass.is_finite() {
                return Err("Invalid max mass");
            }

            let max_radius = radii
                .iter()
                .zip(masses.iter())
                .filter(|(_, m)| (**m - max_mass).abs() < 0.01)
                .map(|(r, _)| *r)
                .next();

            // Try to save results, but don't fail if it doesn't work
            let base_dir = format!("output/nlem_log_limits/{}/B_{}", model_name, b_string);
            let dir_path = format!("{}/csi_{:.2e}", base_dir, csi);
            let _ = fs::create_dir_all(&dir_path);
            let eos_filename = format!("{}/eos.dat", dir_path);
            let _ = write_eos_with_mr(
                &results_eos,
                &masses,
                &radii,
                &b_masses,
                &central_p_list,
                &eos_filename,
            );

            Ok((max_mass, max_radius))
        })() {
            Ok((max_mass, max_radius)) => {
                results.push(FieldTestResult {
                    b_field,
                    csi,
                    max_mass: Some(max_mass),
                    max_radius,
                    success: true,
                    error_msg: String::new(),
                });
            }
            Err(err_msg) => {
                results.push(FieldTestResult {
                    b_field,
                    csi,
                    max_mass: None,
                    max_radius: None,
                    success: false,
                    error_msg: err_msg.to_string(),
                });
            }
        }
    }

    results
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // Parse arguments
    let (exp_start, exp_end) = if args.len() >= 3 {
        (
            args[1]
                .parse::<i32>()
                .expect("exp_start must be an integer"),
            args[2].parse::<i32>().expect("exp_end must be an integer"),
        )
    } else {
        // Default: test from 1e14 to 1e24
        eprintln!("Usage: {} <exp_start> <exp_end>", args[0]);
        eprintln!("Example: {} 14 24", args[0]);
        eprintln!("\nUsing defaults: B from 1e14 to 1e24");
        (14, 24)
    };

    println!("=== Neutron Star Magnetic Field Limits ===");
    println!("Testing B fields from 10^{} to 10^{} G", exp_start, exp_end);
    println!("Using adaptive csi values based on effective field expressions");

    let models = [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)];
    let mut all_field_results: Vec<(String, f64, Vec<FieldTestResult>)> = Vec::new();
    let mut no_nlem_limits: Vec<(String, f64)> = Vec::new();

    // PHASE 1: Find limits WITHOUT NLEM
    println!("\n╔═══════════════════════════════════════════════════════════╗");
    println!("║ PHASE 1: Finding Limits WITHOUT NLEM (Baseline)          ║");

    println!("╚═══════════════════════════════════════════════════════════╝\n");

    for (model_name, model_params) in &models {
        println!("\n=== Model: {} (No NLEM) ===", model_name);
        let mut last_success_b = 0.0;

        // Generate magnetic field values logarithmically
        for exp in exp_start..=exp_end {
            let b_field = 10_f64.powf(exp as f64);
            let result = test_magnetic_field_no_nlem(*model_params, b_field, model_name);

            let b_str = format!("{:.2e}", b_field);
            if result.success {
                if let Some(mass) = result.max_mass {
                    println!("  B = {} G | SUCCESS | M_max = {:.3} M☉", b_str, mass);
                    last_success_b = b_field;
                }
            } else {
                println!("  B = {} G | FAILURE | {}", b_str, result.error_msg);
                if last_success_b > 0.0 {
                    println!(
                        "  → Limit found: between 10^{:.1} and 10^{:.1} G",
                        last_success_b.log10(),
                        b_field.log10()
                    );
                    no_nlem_limits.push((model_name.to_string(), last_success_b));
                    break;
                }
            }
        }
    }

    // PHASE 2: Test with NLEM from baseline
    println!("\n╔═══════════════════════════════════════════════════════════╗");
    println!("║ PHASE 2: Testing with NLEM (adaptive csi per field)      ║");
    println!("╚═══════════════════════════════════════════════════════════╝\n");

    // Test each model
    for (model_name, model_params) in &models {
        println!("\n=== Model: {} (With NLEM) ===", model_name);

        // Determine starting exponent for NLEM tests
        let base_limit = no_nlem_limits
            .iter()
            .find(|(m, _)| m == model_name)
            .map(|(_, b)| b.log10().ceil() as i32)
            .unwrap_or(exp_start);

        // Generate magnetic field values logarithmically, starting from base limit
        for exp in base_limit..=exp_end {
            let b_field = 10_f64.powf(exp as f64);

            // For this bare B field, generate csi values that produce effective fields
            // from 1% to 99% of the bare field (logarithmic spacing)
            let csi_vals = generate_csi_values_for_field(b_field, 25);

            if csi_vals.is_empty() {
                println!(
                    "  B = {:.2e} G | WARNING: Could not generate valid csi values",
                    b_field
                );
                continue;
            }

            let results = test_magnetic_field(*model_params, b_field, &csi_vals, model_name);
            all_field_results.push((model_name.to_string(), b_field, results));
        }
    }

    // Analyze and report results
    println!("\n=== COMPARATIVE ANALYSIS: NLEM vs NO-NLEM ===\n");

    for (model_name, b_field, results) in &all_field_results {
        let success_count = results.iter().filter(|r| r.success).count();
        let b_str = format!("{:.2e}", b_field);

        if success_count == 0 {
            println!(
                "Model {} | B = {} G | FAILURE: 0/{} csi values",
                model_name,
                b_str,
                results.len()
            );
        } else if success_count == results.len() {
            // Find max mass across successful results
            let max_mass = results
                .iter()
                .filter(|r| r.success)
                .filter_map(|r| r.max_mass)
                .fold(f64::NEG_INFINITY, f64::max);

            println!(
                "Model {} | B = {} G | SUCCESS: {}/{} | Max M = {:.3} M☉",
                model_name,
                b_str,
                success_count,
                results.len(),
                max_mass
            );
        } else {
            println!(
                "Model {} | B = {} G | PARTIAL: {}/{}",
                model_name,
                b_str,
                success_count,
                results.len()
            );
        }
    }

    // Compare limits
    println!("\n=== NO-NLEM BASELINE LIMITS ===\n");
    for (model_name, limit_b) in &no_nlem_limits {
        println!("Model {} | Limit: {:.2e} G", model_name, limit_b);
    }

    println!("\n=== NLEM ENHANCEMENT ===\n");
    for (model_name, _) in &models {
        if let Some((_, baseline_b)) = no_nlem_limits.iter().find(|(m, _)| m == model_name) {
            let nlem_results: Vec<_> = all_field_results
                .iter()
                .filter(|(m, _, _)| m == model_name)
                .collect();

            if let Some((_, highest_b, res)) = nlem_results.iter().max_by(|a, b| {
                let a_success = a.2.iter().filter(|r| r.success).count();
                let b_success = b.2.iter().filter(|r| r.success).count();
                a_success.cmp(&b_success)
            }) {
                let success_count = res.iter().filter(|r| r.success).count();
                if success_count > 0 {
                    let enhancement = (highest_b / baseline_b) - 1.0;
                    println!(
                        "Model {} | Baseline: {:.2e} G → NLEM extends to: {:.2e} G (+{:.1}%)",
                        model_name,
                        baseline_b,
                        highest_b,
                        enhancement * 100.0
                    );
                }
            }
        }
    }

    // Find breaking points for each model
    println!("\n=== BREAKING POINTS ===\n");

    for (model_name, _) in &models {
        let model_results: Vec<_> = all_field_results
            .iter()
            .filter(|(m, _, _)| m == model_name)
            .collect();

        if let Some((_, b_before, _res_before)) = model_results
            .iter()
            .rev()
            .find(|(_, _, r)| r.iter().filter(|x| x.success).count() > r.len() / 2)
        {
            if let Some(idx) = model_results
                .iter()
                .position(|(_, b, _)| (b - b_before).abs() < 1.0 && b > b_before)
            {
                if idx + 1 < model_results.len() {
                    let (_, b_after, res_after) = model_results[idx + 1];
                    if res_after.iter().filter(|x| x.success).count() == 0 {
                        println!(
                            "Model {} breaks between 10^{:.1} G and 10^{:.1} G",
                            model_name,
                            b_before.log10(),
                            b_after.log10()
                        );
                    }
                }
            }
        }
    }

    // Save detailed report
    let report_path = "output/nlem_log_limits/field_limits_report.txt";
    if let Err(e) = fs::create_dir_all("output/nlem_log_limits") {
        eprintln!("Warning: Could not create output directory: {}", e);
    } else {
        let mut report = String::new();
        report.push_str("=== NEUTRON STAR MAGNETIC FIELD LIMITS REPORT ===\n\n");

        for (model_name, b_field, results) in &all_field_results {
            report.push_str(&format!("Model: {} | B = {:.2e} G\n", model_name, b_field));
            let success = results.iter().filter(|r| r.success).count();
            report.push_str(&format!(
                "  Status: {}/{} successful\n",
                success,
                results.len()
            ));

            for result in results.iter().take(5) {
                if let Some(mass) = result.max_mass {
                    report.push_str(&format!(
                        "    csi={:.2e}: M_max={:.3} M☉, R={:.2}km\n",
                        result.csi,
                        mass,
                        result.max_radius.unwrap_or(0.0)
                    ));
                } else {
                    report.push_str(&format!(
                        "    csi={:.2e}: FAILED ({})\n",
                        result.csi, result.error_msg
                    ));
                }
            }
            report.push_str("\n");
        }

        if let Err(e) = fs::write(report_path, report) {
            eprintln!("Warning: Could not write report: {}", e);
        } else {
            println!("\nDetailed report saved to {}", report_path);
        }
    }

    println!("\n✓ Magnetic field limit exploration complete!");
}
