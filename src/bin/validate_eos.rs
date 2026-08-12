// src/bin/validate_eos.rs

use std::env;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use nsrs::core::constants::{HBAR_C, N0, RESULTS_SIZE};

const COL_NB_OVER_N0: usize = 0;
const COL_EPS: usize = 1;
const COL_PRESS: usize = 2;
const COL_NE: usize = 3;
const COL_NMU: usize = 4;
const COL_NP: usize = 6;
const COL_NSM: usize = 8;
const COL_NSP: usize = 10;
const COL_NXM: usize = 11;
const COL_MU_TOTAL: usize = 20;
const COL_N_CHI: usize = 21;
const COL_Y_CHI: usize = 22;
const COL_M_CHI: usize = 23;
const COL_M_X: usize = 24;
const COL_EPSILON: usize = 25;
const COL_G_D: usize = 26;
const COL_X0: usize = 27;
const COL_KF_CHI: usize = 28;
const COL_MU_CHI: usize = 29;
const COL_EPS_CHI_KIN: usize = 30;
const COL_P_CHI_KIN: usize = 31;
const COL_EPS_X: usize = 32;
const COL_P_X: usize = 33;
const COL_MR_MASS: usize = RESULTS_SIZE;
const COL_MR_RADIUS: usize = RESULTS_SIZE + 1;
const COL_MR_BARYONIC: usize = RESULTS_SIZE + 2;

const MIN_BASE_COLUMNS: usize = RESULTS_SIZE;
const MIN_MR_COLUMNS: usize = COL_MR_BARYONIC + 1;
const EPS_DIFF_TOL: f64 = 1.0e-10;
const PRESS_DIFF_TOL: f64 = 1.0e-10;
const CHARGE_ABS_TOL: f64 = 1.0e-5;
const CHARGE_REL_TOL: f64 = 1.0e-3;
const MAX_VALID_MASS_MSUN: f64 = 3.0;
const MAX_VALID_RADIUS_KM: f64 = 20.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Severity {
    Pass,
    Warn,
    Fail,
}

impl Severity {
    fn as_str(self) -> &'static str {
        match self {
            Severity::Pass => "PASS",
            Severity::Warn => "WARN",
            Severity::Fail => "FAIL",
        }
    }

    fn combine(self, other: Severity) -> Severity {
        self.max(other)
    }
}

#[derive(Clone, Debug)]
struct Config {
    inputs: Vec<PathBuf>,
    csv_path: Option<PathBuf>,
    strict: bool,
    require_mr: bool,
    min_max_mass: f64,
    heavy_mass_warning: f64,
    r14_min: f64,
    r14_max: f64,
    r208_min: f64,
    r208_max: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            csv_path: None,
            strict: false,
            require_mr: false,
            min_max_mass: 2.01,
            heavy_mass_warning: 2.35,
            r14_min: 10.5,
            r14_max: 13.8,
            r208_min: 10.0,
            r208_max: 15.0,
        }
    }
}

#[derive(Clone, Debug)]
struct EosData {
    rows: Vec<Vec<f64>>,
    max_cols: usize,
}

#[derive(Clone, Debug)]
struct Check {
    name: &'static str,
    severity: Severity,
    detail: String,
}

#[derive(Clone, Debug)]
struct FileReport {
    path: PathBuf,
    row_count: usize,
    has_mr: bool,
    overall: Severity,
    checks: Vec<Check>,
}

fn main() {
    let config = match parse_args() {
        Ok(config) => config,
        Err(message) => {
            eprintln!("{message}");
            print_usage();
            std::process::exit(2);
        }
    };

    let files = collect_inputs(&config.inputs);
    if files.is_empty() {
        eprintln!("No .dat EoS files found.");
        std::process::exit(2);
    }

    let mut reports = Vec::new();
    for path in files {
        reports.push(validate_file(&path, &config));
    }

    print_reports(&reports);

    if let Some(csv_path) = &config.csv_path {
        if let Err(err) = write_csv(&reports, csv_path) {
            eprintln!("Could not write CSV {}: {err}", csv_path.display());
            std::process::exit(2);
        }
        println!("\nCSV report written to {}", csv_path.display());
    }

    let failed = reports.iter().any(|r| r.overall == Severity::Fail);
    let warned = reports.iter().any(|r| r.overall == Severity::Warn);
    if failed || (config.strict && warned) {
        std::process::exit(1);
    }
}

fn parse_args() -> Result<Config, String> {
    let mut config = Config::default();
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            "--csv" => {
                config.csv_path = Some(PathBuf::from(next_value(&mut args, "--csv")?));
            }
            "--strict" => {
                config.strict = true;
            }
            "--require-mr" => {
                config.require_mr = true;
            }
            "--min-max-mass" => {
                config.min_max_mass =
                    parse_f64(next_value(&mut args, "--min-max-mass")?, "--min-max-mass")?;
            }
            "--heavy-mass-warning" => {
                config.heavy_mass_warning = parse_f64(
                    next_value(&mut args, "--heavy-mass-warning")?,
                    "--heavy-mass-warning",
                )?;
            }
            "--r14-min" => {
                config.r14_min = parse_f64(next_value(&mut args, "--r14-min")?, "--r14-min")?;
            }
            "--r14-max" => {
                config.r14_max = parse_f64(next_value(&mut args, "--r14-max")?, "--r14-max")?;
            }
            "--r208-min" => {
                config.r208_min = parse_f64(next_value(&mut args, "--r208-min")?, "--r208-min")?;
            }
            "--r208-max" => {
                config.r208_max = parse_f64(next_value(&mut args, "--r208-max")?, "--r208-max")?;
            }
            _ if arg.starts_with('-') => return Err(format!("Unknown option: {arg}")),
            _ => config.inputs.push(PathBuf::from(arg)),
        }
    }

    if config.inputs.is_empty() {
        return Err("At least one EoS file or directory is required.".to_string());
    }
    if config.r14_min >= config.r14_max {
        return Err("--r14-min must be smaller than --r14-max".to_string());
    }
    if config.r208_min >= config.r208_max {
        return Err("--r208-min must be smaller than --r208-max".to_string());
    }

    Ok(config)
}

fn next_value(args: &mut impl Iterator<Item = String>, opt: &str) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("Missing value after {opt}"))
}

fn parse_f64(text: String, opt: &str) -> Result<f64, String> {
    let value = text
        .parse::<f64>()
        .map_err(|_| format!("Invalid numeric value for {opt}: {text}"))?;
    if !value.is_finite() {
        return Err(format!("Value for {opt} must be finite"));
    }
    Ok(value)
}

fn print_usage() {
    eprintln!(
        "Usage: cargo run --bin validate_eos -- [OPTIONS] <eos.dat|directory>...\n\
\n\
Options:\n\
  --csv <path>                 Write machine-readable validation report\n\
  --strict                     Exit non-zero on warnings as well as failures\n\
  --require-mr                 Fail files without M-R columns\n\
  --min-max-mass <Msun>        Hard lower maximum-mass threshold [default: 2.01]\n\
  --heavy-mass-warning <Msun>  Soft massive-pulsar benchmark [default: 2.35]\n\
  --r14-min <km>               Soft lower R(1.4 Msun) band [default: 10.5]\n\
  --r14-max <km>               Soft upper R(1.4 Msun) band [default: 13.8]\n\
  --r208-min <km>              Soft lower R(2.08 Msun) band [default: 10.0]\n\
  --r208-max <km>              Soft upper R(2.08 Msun) band [default: 15.0]\n\
\n\
Checks include finite data, positive EoS, monotonic P(epsilon), 0 <= c_s^2 <= 1,\n\
dominant-energy sanity P <= epsilon, charge neutrality, non-negative populations,\n\
M-R maximum mass, stable branch monotonicity, compactness, and radius bands."
    );
}

fn collect_inputs(inputs: &[PathBuf]) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for input in inputs {
        collect_one(input, &mut files);
    }
    files.sort();
    files.dedup();
    files
}

fn collect_one(path: &Path, files: &mut Vec<PathBuf>) {
    if path.is_file() {
        if is_dat_file(path) {
            files.push(path.to_path_buf());
        }
        return;
    }

    if !path.is_dir() {
        return;
    }

    let entries = match fs::read_dir(path) {
        Ok(entries) => entries,
        Err(_) => return,
    };

    for entry in entries.flatten() {
        collect_one(&entry.path(), files);
    }
}

fn is_dat_file(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("dat"))
        .unwrap_or(false)
}

fn validate_file(path: &Path, config: &Config) -> FileReport {
    let mut checks = Vec::new();
    let data = match load_eos(path) {
        Ok(data) => data,
        Err(err) => {
            checks.push(Check {
                name: "load",
                severity: Severity::Fail,
                detail: err,
            });
            return finalize_report(path, 0, false, checks);
        }
    };

    let row_count = data.rows.len();
    let has_mr = has_complete_mr(&data);

    checks.extend(check_shape(&data, config));
    checks.extend(check_thermodynamics(&data));
    checks.extend(check_composition(&data));
    checks.extend(check_dark_diagnostics(&data));
    checks.extend(check_total_chemical_potential(&data));
    checks.extend(check_mass_radius(&data, config));

    finalize_report(path, row_count, has_mr, checks)
}

fn finalize_report(path: &Path, row_count: usize, has_mr: bool, checks: Vec<Check>) -> FileReport {
    let overall = checks
        .iter()
        .fold(Severity::Pass, |acc, check| acc.combine(check.severity));
    FileReport {
        path: path.to_path_buf(),
        row_count,
        has_mr,
        overall,
        checks,
    }
}

fn load_eos(path: &Path) -> Result<EosData, String> {
    let file = fs::File::open(path).map_err(|err| format!("open failed: {err}"))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    let mut max_cols = 0;

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|err| format!("read line {} failed: {err}", line_no + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let mut row = Vec::new();
        for token in trimmed.split_whitespace() {
            let value = token
                .parse::<f64>()
                .map_err(|_| format!("line {} has non-numeric token '{token}'", line_no + 1))?;
            row.push(value);
        }
        max_cols = max_cols.max(row.len());
        rows.push(row);
    }

    if rows.is_empty() {
        return Err("no numeric rows found".to_string());
    }

    Ok(EosData { rows, max_cols })
}

fn check_shape(data: &EosData, config: &Config) -> Vec<Check> {
    let mut checks = Vec::new();
    let malformed_rows = data
        .rows
        .iter()
        .filter(|row| row.len() != MIN_BASE_COLUMNS && row.len() != MIN_MR_COLUMNS)
        .count();
    let mixed_widths = data
        .rows
        .first()
        .map(|first| data.rows.iter().any(|row| row.len() != first.len()))
        .unwrap_or(false);

    checks.push(Check {
        name: "columns",
        severity: if malformed_rows == 0 && !mixed_widths {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "rows={}, max_cols={}, invalid_width_rows={}, mixed_widths={}",
            data.rows.len(),
            data.max_cols,
            malformed_rows,
            mixed_widths
        ),
    });

    let has_mr = has_complete_mr(data);
    checks.push(Check {
        name: "mr_columns",
        severity: if has_mr {
            Severity::Pass
        } else if config.require_mr {
            Severity::Fail
        } else {
            Severity::Warn
        },
        detail: if has_mr {
            "M-R columns present".to_string()
        } else {
            "M-R columns absent; observational M-R checks skipped".to_string()
        },
    });

    checks
}

fn check_thermodynamics(data: &EosData) -> Vec<Check> {
    let mut checks = Vec::new();
    let usable: Vec<&Vec<f64>> = data
        .rows
        .iter()
        .filter(|row| row.len() >= MIN_BASE_COLUMNS)
        .collect();

    let finite_violations = usable
        .iter()
        .flat_map(|row| row.iter().take(MIN_BASE_COLUMNS))
        .filter(|value| !value.is_finite())
        .count();

    checks.push(Check {
        name: "finite",
        severity: if finite_violations == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!("non_finite_values={finite_violations}"),
    });

    let nonpositive_eps = usable.iter().filter(|row| row[COL_EPS] <= 0.0).count();
    let negative_press = usable.iter().filter(|row| row[COL_PRESS] < 0.0).count();
    checks.push(Check {
        name: "positive_eos",
        severity: if nonpositive_eps == 0 && negative_press == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!("eps<=0 rows={nonpositive_eps}, pressure<0 rows={negative_press}"),
    });

    let mut eps_p: Vec<(f64, f64)> = usable
        .iter()
        .filter_map(|row| {
            let eps = row[COL_EPS];
            let press = row[COL_PRESS];
            if eps.is_finite() && press.is_finite() {
                Some((eps, press))
            } else {
                None
            }
        })
        .collect();
    eps_p.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut eps_nonmonotonic = 0;
    let mut press_nonmonotonic = 0;
    let mut cs2_negative = 0;
    let mut cs2_superluminal = 0;
    let mut cs2_values = Vec::new();
    let mut dominant_energy_violations = 0;

    for row in &usable {
        let eps = row[COL_EPS];
        let press = row[COL_PRESS];
        if eps.is_finite() && press.is_finite() && press > eps + PRESS_DIFF_TOL {
            dominant_energy_violations += 1;
        }
    }

    for pair in eps_p.windows(2) {
        let (eps0, p0) = pair[0];
        let (eps1, p1) = pair[1];
        let de = eps1 - eps0;
        let dp = p1 - p0;

        if de <= EPS_DIFF_TOL {
            eps_nonmonotonic += 1;
            continue;
        }
        if dp < -PRESS_DIFF_TOL {
            press_nonmonotonic += 1;
        }

        let cs2 = dp / de;
        if cs2.is_finite() {
            if cs2 < -1.0e-8 {
                cs2_negative += 1;
            }
            if cs2 > 1.0 + 1.0e-8 {
                cs2_superluminal += 1;
            }
            cs2_values.push(cs2);
        }
    }

    let cs2_min = finite_min(&cs2_values);
    let cs2_max = finite_max(&cs2_values);

    checks.push(Check {
        name: "monotonic_eos",
        severity: if eps_nonmonotonic == 0 && press_nonmonotonic == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "eps_nonmonotonic_pairs={eps_nonmonotonic}, p_nonmonotonic_pairs={press_nonmonotonic}"
        ),
    });

    checks.push(Check {
        name: "causality_cs2",
        severity: if cs2_negative == 0 && cs2_superluminal == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "cs2_min={}, cs2_max={}, negative_pairs={}, superluminal_pairs={}",
            fmt_opt(cs2_min),
            fmt_opt(cs2_max),
            cs2_negative,
            cs2_superluminal
        ),
    });

    checks.push(Check {
        name: "dominant_energy",
        severity: if dominant_energy_violations == 0 {
            Severity::Pass
        } else {
            Severity::Warn
        },
        detail: format!("rows_with_pressure_greater_than_energy={dominant_energy_violations}"),
    });

    checks
}

fn check_composition(data: &EosData) -> Vec<Check> {
    let usable: Vec<&Vec<f64>> = data
        .rows
        .iter()
        .filter(|row| row.len() >= MIN_BASE_COLUMNS)
        .collect();

    let mut checks = Vec::new();
    let mut negative_population_rows = 0;
    let mut max_charge_abs = 0.0_f64;
    let mut max_charge_rel = 0.0_f64;

    for row in &usable {
        if (3..=12).any(|col| row[col].is_finite() && row[col] < -1.0e-12) {
            negative_population_rows += 1;
        }

        let charge =
            row[COL_NP] + row[COL_NSP] - row[COL_NSM] - row[COL_NXM] - row[COL_NE] - row[COL_NMU];
        let nb = (row[COL_NB_OVER_N0] * N0).abs().max(1.0e-30);
        max_charge_abs = max_charge_abs.max(charge.abs());
        max_charge_rel = max_charge_rel.max(charge.abs() / nb);
    }

    checks.push(Check {
        name: "nonnegative_populations",
        severity: if negative_population_rows == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!("rows_with_negative_population={negative_population_rows}"),
    });

    checks.push(Check {
        name: "charge_neutrality",
        severity: if max_charge_abs <= CHARGE_ABS_TOL || max_charge_rel <= CHARGE_REL_TOL {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "max_abs_charge={:.6e} fm^-3, max_abs_charge_over_nB={:.6e}",
            max_charge_abs, max_charge_rel
        ),
    });

    checks
}

fn check_dark_diagnostics(data: &EosData) -> Vec<Check> {
    let mut dark_rows = 0;
    let mut invalid_parameters = 0;
    let mut negative_thermodynamics = 0;
    let mut fraction_violations = 0;
    let mut fermi_momentum_violations = 0;
    let mut chemical_potential_violations = 0;
    let mut proca_violations = 0;
    let mut vector_pressure_violations = 0;

    for row in data.rows.iter().filter(|row| row.len() >= MIN_BASE_COLUMNS) {
        let is_dark = row[COL_N_CHI..=COL_P_X]
            .iter()
            .any(|value| !value.is_finite() || *value != 0.0);
        if !is_dark {
            continue;
        }
        dark_rows += 1;

        let n_chi = row[COL_N_CHI];
        let y_chi = row[COL_Y_CHI];
        let m_chi = row[COL_M_CHI];
        let m_x = row[COL_M_X];
        let epsilon = row[COL_EPSILON];
        let g_d = row[COL_G_D];
        let x0 = row[COL_X0];
        let kf_chi = row[COL_KF_CHI];
        let mu_chi = row[COL_MU_CHI];
        let eps_chi = row[COL_EPS_CHI_KIN];
        let p_chi = row[COL_P_CHI_KIN];
        let eps_x = row[COL_EPS_X];
        let p_x = row[COL_P_X];

        let parameters_valid = n_chi.is_finite()
            && n_chi >= 0.0
            && y_chi.is_finite()
            && y_chi >= 0.0
            && m_chi.is_finite()
            && m_chi > 0.0
            && m_x.is_finite()
            && m_x > 0.0
            && epsilon.is_finite()
            && epsilon.abs() < 1.0
            && g_d.is_finite()
            && x0.is_finite()
            && kf_chi.is_finite()
            && kf_chi >= 0.0
            && mu_chi.is_finite();
        if !parameters_valid {
            invalid_parameters += 1;
            continue;
        }

        if !eps_chi.is_finite()
            || !p_chi.is_finite()
            || !eps_x.is_finite()
            || !p_x.is_finite()
            || eps_chi < -1.0e-10
            || p_chi < -1.0e-10
            || eps_x < -1.0e-10
            || p_x < -1.0e-10
        {
            negative_thermodynamics += 1;
        }

        let expected_n_chi = y_chi * row[COL_NB_OVER_N0] * N0;
        if !approximately_equal(n_chi, expected_n_chi, 2.0e-10, 2.0e-5) {
            fraction_violations += 1;
        }

        let expected_n_from_kf =
            kf_chi.powi(3) / (3.0 * std::f64::consts::PI.powi(2) * HBAR_C.powi(3));
        if !approximately_equal(n_chi, expected_n_from_kf, 2.0e-10, 3.0e-5) {
            fermi_momentum_violations += 1;
        }

        if n_chi > 1.0e-14 {
            let norm = (1.0 - epsilon.powi(2)).sqrt();
            let expected_mu = kf_chi.hypot(m_chi) + g_d * x0 / norm;
            if !approximately_equal(mu_chi, expected_mu, 2.0e-4, 3.0e-5) {
                chemical_potential_violations += 1;
            }

            let proca_left = m_x.powi(2) * x0;
            let proca_right = g_d * n_chi * HBAR_C.powi(3) / norm;
            if !approximately_equal(proca_left, proca_right, 2.0e-4, 5.0e-5) {
                proca_violations += 1;
            }
        } else if kf_chi.abs() > 1.0e-10 || mu_chi.abs() > 1.0e-10 {
            chemical_potential_violations += 1;
        }

        let expected_eps_x = 0.5 * m_x.powi(2) * x0.powi(2) / HBAR_C.powi(3);
        if !approximately_equal(eps_x, expected_eps_x, 2.0e-10, 5.0e-5)
            || !approximately_equal(p_x, eps_x, 2.0e-10, 2.0e-5)
        {
            vector_pressure_violations += 1;
        }
    }

    let violations = invalid_parameters
        + negative_thermodynamics
        + fraction_violations
        + fermi_momentum_violations
        + chemical_potential_violations
        + proca_violations
        + vector_pressure_violations;

    vec![Check {
        name: "dark_diagnostics",
        severity: if violations == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "dark_rows={dark_rows}, invalid_parameters={invalid_parameters}, negative_thermodynamics={negative_thermodynamics}, fraction={fraction_violations}, kf={fermi_momentum_violations}, mu_chi={chemical_potential_violations}, proca={proca_violations}, vector_pressure={vector_pressure_violations}"
        ),
    }]
}

fn approximately_equal(left: f64, right: f64, abs_tol: f64, rel_tol: f64) -> bool {
    (left - right).abs() <= abs_tol.max(rel_tol * left.abs().max(right.abs()))
}

fn check_total_chemical_potential(data: &EosData) -> Vec<Check> {
    let usable: Vec<&Vec<f64>> = data
        .rows
        .iter()
        .filter(|row| row.len() >= MIN_BASE_COLUMNS)
        .collect();

    let mut nonfinite = 0;
    let mut nonpositive = 0;
    let mut min_mu = f64::INFINITY;
    let mut max_mu = f64::NEG_INFINITY;

    for row in usable {
        let mu_total = row[COL_MU_TOTAL];
        if !mu_total.is_finite() {
            nonfinite += 1;
            continue;
        }
        if mu_total <= 0.0 {
            nonpositive += 1;
        }
        min_mu = min_mu.min(mu_total);
        max_mu = max_mu.max(mu_total);
    }

    vec![Check {
        name: "total_chemical_potential",
        severity: if nonfinite == 0 && nonpositive == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "mu_total_min={}, mu_total_max={}, nonfinite={}, nonpositive={}",
            fmt_opt(finite_value(min_mu)),
            fmt_opt(finite_value(max_mu)),
            nonfinite,
            nonpositive
        ),
    }]
}

fn check_mass_radius(data: &EosData, config: &Config) -> Vec<Check> {
    let mut checks = Vec::new();
    if !has_complete_mr(data) {
        return checks;
    }

    let mut mr: Vec<(f64, f64, f64)> = data
        .rows
        .iter()
        .filter(|row| row.len() >= MIN_MR_COLUMNS)
        .filter_map(|row| {
            let mass = row[COL_MR_MASS];
            let radius = row[COL_MR_RADIUS];
            let central_p = row[COL_PRESS];
            if is_valid_mr_point(mass, radius) && central_p.is_finite() {
                Some((central_p, mass, radius))
            } else {
                None
            }
        })
        .collect();

    mr.sort_by(|a, b| a.0.total_cmp(&b.0));

    if mr.len() < 3 {
        checks.push(Check {
            name: "mr_valid_points",
            severity: Severity::Fail,
            detail: format!("valid_mr_points={}", mr.len()),
        });
        return checks;
    }

    checks.push(Check {
        name: "mr_valid_points",
        severity: Severity::Pass,
        detail: format!("valid_mr_points={}", mr.len()),
    });

    let (max_idx, max_mass, radius_at_max) = max_mass_point(&mr);
    checks.push(Check {
        name: "max_mass",
        severity: if max_mass >= config.min_max_mass {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!(
            "Mmax={:.6} Msun, R_at_Mmax={:.6} km, required_Mmax>={:.3} Msun",
            max_mass, radius_at_max, config.min_max_mass
        ),
    });

    checks.push(Check {
        name: "heavy_pulsar_benchmark",
        severity: if max_mass >= config.heavy_mass_warning {
            Severity::Pass
        } else {
            Severity::Warn
        },
        detail: format!(
            "Mmax={:.6} Msun compared with high-mass benchmark {:.3} Msun",
            max_mass, config.heavy_mass_warning
        ),
    });

    let stable_violations = stable_branch_violations(&mr, max_idx);
    checks.push(Check {
        name: "stable_branch",
        severity: if stable_violations == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!("pre_Mmax_dM_dPc_negative_pairs={stable_violations}"),
    });

    let compactness_violations = compactness_violations(&mr);
    checks.push(Check {
        name: "compactness",
        severity: if compactness_violations == 0 {
            Severity::Pass
        } else {
            Severity::Fail
        },
        detail: format!("points_inside_schwarzschild_radius={compactness_violations}"),
    });

    checks.push(radius_band_check(
        "radius_1p4",
        &mr,
        1.4,
        config.r14_min,
        config.r14_max,
    ));
    checks.push(radius_band_check(
        "radius_2p08",
        &mr,
        2.08,
        config.r208_min,
        config.r208_max,
    ));

    checks
}

fn has_complete_mr(data: &EosData) -> bool {
    !data.rows.is_empty() && data.rows.iter().all(|row| row.len() == MIN_MR_COLUMNS)
}

fn is_valid_mr_point(mass: f64, radius: f64) -> bool {
    mass.is_finite()
        && radius.is_finite()
        && mass > 0.0
        && radius > 0.0
        && mass <= MAX_VALID_MASS_MSUN
        && radius <= MAX_VALID_RADIUS_KM
}

fn max_mass_point(mr: &[(f64, f64, f64)]) -> (usize, f64, f64) {
    let mut max_idx = 0;
    let mut max_mass = f64::NEG_INFINITY;
    let mut radius_at_max = f64::NAN;
    for (idx, &(_pc, mass, radius)) in mr.iter().enumerate() {
        if mass > max_mass {
            max_mass = mass;
            radius_at_max = radius;
            max_idx = idx;
        }
    }
    (max_idx, max_mass, radius_at_max)
}

fn stable_branch_violations(mr: &[(f64, f64, f64)], max_idx: usize) -> usize {
    mr[..=max_idx]
        .windows(2)
        .filter(|pair| pair[1].1 + 1.0e-8 < pair[0].1)
        .count()
}

fn compactness_violations(mr: &[(f64, f64, f64)]) -> usize {
    const GM_SUN_OVER_C2_KM: f64 = 1.47662504;
    mr.iter()
        .filter(|&&(_pc, mass, radius)| radius <= 2.0 * GM_SUN_OVER_C2_KM * mass)
        .count()
}

fn radius_band_check(
    name: &'static str,
    mr: &[(f64, f64, f64)],
    target_mass: f64,
    min_radius: f64,
    max_radius: f64,
) -> Check {
    match interpolate_radius_at_mass(mr, target_mass) {
        Some(radius) => Check {
            name,
            severity: if radius >= min_radius && radius <= max_radius {
                Severity::Pass
            } else {
                Severity::Warn
            },
            detail: format!(
                "R({:.2} Msun)={:.6} km, soft_band=[{:.3}, {:.3}] km",
                target_mass, radius, min_radius, max_radius
            ),
        },
        None => Check {
            name,
            severity: Severity::Warn,
            detail: format!(
                "target mass {:.2} Msun not bracketed by M-R curve",
                target_mass
            ),
        },
    }
}

fn interpolate_radius_at_mass(mr: &[(f64, f64, f64)], target_mass: f64) -> Option<f64> {
    for pair in mr.windows(2) {
        let (_pc0, m0, r0) = pair[0];
        let (_pc1, m1, r1) = pair[1];
        if (target_mass - m0) * (target_mass - m1) <= 0.0 && (m1 - m0).abs() > 1.0e-12 {
            let t = (target_mass - m0) / (m1 - m0);
            return Some(r0 + t * (r1 - r0));
        }
    }
    None
}

fn finite_min(values: &[f64]) -> Option<f64> {
    values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .min_by(|a, b| a.total_cmp(b))
}

fn finite_max(values: &[f64]) -> Option<f64> {
    values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .max_by(|a, b| a.total_cmp(b))
}

fn finite_value(value: f64) -> Option<f64> {
    if value.is_finite() { Some(value) } else { None }
}

fn fmt_opt(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{v:.6e}"),
        None => "n/a".to_string(),
    }
}

fn print_reports(reports: &[FileReport]) {
    let mut pass = 0;
    let mut warn = 0;
    let mut fail = 0;

    for report in reports {
        match report.overall {
            Severity::Pass => pass += 1,
            Severity::Warn => warn += 1,
            Severity::Fail => fail += 1,
        }

        println!(
            "\n[{}] {} (rows={}, mr={})",
            report.overall.as_str(),
            report.path.display(),
            report.row_count,
            report.has_mr
        );
        for check in &report.checks {
            println!(
                "  {:>4} {:<24} {}",
                check.severity.as_str(),
                check.name,
                check.detail
            );
        }
    }

    println!(
        "\nSummary: files={}, pass={}, warn={}, fail={}",
        reports.len(),
        pass,
        warn,
        fail
    );
}

fn write_csv(reports: &[FileReport], csv_path: &Path) -> std::io::Result<()> {
    if let Some(parent) = csv_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }

    let mut file = fs::File::create(csv_path)?;
    writeln!(file, "path,overall,rows,has_mr,check,severity,detail")?;
    for report in reports {
        for check in &report.checks {
            writeln!(
                file,
                "{},{},{},{},{},{},{}",
                csv_escape(&report.path.display().to_string()),
                report.overall.as_str(),
                report.row_count,
                report.has_mr,
                check.name,
                check.severity.as_str(),
                csv_escape(&check.detail)
            )?;
        }
    }
    Ok(())
}

fn csv_escape(text: &str) -> String {
    if text.contains(',') || text.contains('"') || text.contains('\n') {
        format!("\"{}\"", text.replace('"', "\"\""))
    } else {
        text.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_dark_row() -> Vec<f64> {
        let mut row = vec![0.0; RESULTS_SIZE];
        row[COL_NB_OVER_N0] = 1.0;
        row[COL_Y_CHI] = 0.02;
        row[COL_N_CHI] = row[COL_Y_CHI] * N0;
        row[COL_M_CHI] = 750.0;
        row[COL_M_X] = 100.0;
        row[COL_EPSILON] = 0.08;
        row[COL_G_D] = 0.35;

        let norm = (1.0 - row[COL_EPSILON].powi(2)).sqrt();
        row[COL_KF_CHI] =
            (3.0 * std::f64::consts::PI.powi(2) * row[COL_N_CHI] * HBAR_C.powi(3)).cbrt();
        row[COL_X0] =
            row[COL_G_D] * row[COL_N_CHI] * HBAR_C.powi(3) / (row[COL_M_X].powi(2) * norm);
        row[COL_MU_CHI] = row[COL_KF_CHI].hypot(row[COL_M_CHI]) + row[COL_G_D] * row[COL_X0] / norm;
        row[COL_EPS_CHI_KIN] = 2.5;
        row[COL_P_CHI_KIN] = 0.2;
        row[COL_EPS_X] = 0.5 * row[COL_M_X].powi(2) * row[COL_X0].powi(2) / HBAR_C.powi(3);
        row[COL_P_X] = row[COL_EPS_X];
        row
    }

    #[test]
    fn dark_diagnostic_validator_accepts_consistent_exported_units() {
        let data = EosData {
            rows: vec![valid_dark_row()],
            max_cols: RESULTS_SIZE,
        };
        let check = check_dark_diagnostics(&data).remove(0);
        assert_eq!(check.severity, Severity::Pass, "{}", check.detail);
    }

    #[test]
    fn dark_diagnostic_validator_rejects_broken_fraction_relation() {
        let mut row = valid_dark_row();
        row[COL_N_CHI] *= 2.0;
        let data = EosData {
            rows: vec![row],
            max_cols: RESULTS_SIZE,
        };
        let check = check_dark_diagnostics(&data).remove(0);
        assert_eq!(check.severity, Severity::Fail);
        assert!(check.detail.contains("fraction=1"));
    }
}
