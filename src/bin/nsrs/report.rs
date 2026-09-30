// Relatórios de validação (substitui os binários properties e observations).

use std::io::Write;

use nsrs::core::io_utils::derived_diagnostics;
use nsrs::core::nuclear::saturation_properties;
use nsrs::core::observations::{Status, assess, load_constraints};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver};

use crate::cli::{Args, create_dir};

const HYPERONS: [&str; 6] = ["Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi-", "Xi0"];

/// `report properties`: propriedades nucleares e estelares dos modelos (B = 0):
/// saturação (n0, E/A, K, J, L, M*/M), sequência estelar com crosta
/// (M_max, R_1.4, Lambda_1.4, I_1.4, z_1.4), limiar de URCA direto, início
/// dos hyperons e c_s^2 máximo no ramo estável.
pub fn properties(_raw: &[String]) -> Result<(), String> {
    for (name, model) in [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)] {
        println!("== {name}");
        match saturation_properties(model) {
            Some(s) => println!(
                "  saturação: n0 = {:.4} fm^-3, E/A = {:.2} MeV, K = {:.1} MeV, J = {:.2} MeV, L = {:.1} MeV, M*/M = {:.3}",
                s.n0,
                s.energy_per_nucleon,
                s.incompressibility,
                s.symmetry_energy,
                s.symmetry_slope,
                s.effective_mass
            ),
            None => println!("  saturação: não encontrada"),
        }

        for hyperons in [true, false] {
            println!(
                "  -- {}",
                if hyperons {
                    "com hyperons"
                } else {
                    "só núcleons (npe mu)"
                }
            );
            // Malha em mu até 3 M_N: com o padrão (1.8 M_N) EoS rígidas terminam
            // antes do centro da estrela de massa máxima.
            let engine = HadronsMatter::new(model, 0.0)
                .with_hyperons(hyperons)
                .with_limits(0.02, 3.0)
                .with_points(2000);
            let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
            let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
            let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
            let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
            let stars = generate_star_sequence(&eps, &p, &n, true);
            let Some(i_max) =
                (0..stars.len()).max_by(|&a, &b| stars[a].mass.total_cmp(&stars[b].mass))
            else {
                println!("  sem sequência estelar");
                continue;
            };
            let max = stars[i_max];
            if i_max + 10 >= stars.len() {
                println!(
                    "  AVISO: massa máxima no fim da sequência (EoS curta); M_max é um limite inferior"
                );
            }
            println!(
                "  estrelas (com crosta): M_max = {:.3} Msun, R(M_max) = {:.2} km, P_c = {:.1} MeV/fm^3",
                max.mass, max.radius, max.central_pressure
            );
            if let Some(k) = (1..=i_max).find(|&k| stars[k].mass >= 1.4) {
                let (a, b) = (stars[k - 1], stars[k]);
                let t = (1.4 - a.mass) / (b.mass - a.mass);
                let lerp = |x: f64, y: f64| x + t * (y - x);
                println!(
                    "  1.4 Msun: R = {:.2} km, Lambda = {:.0}, k2 = {:.4}, I = {:.3}e45 g cm^2, z = {:.3}",
                    lerp(a.radius, b.radius),
                    lerp(a.tidal_deformability, b.tidal_deformability),
                    lerp(a.love_k2, b.love_k2),
                    lerp(a.moment_of_inertia, b.moment_of_inertia),
                    lerp(a.redshift, b.redshift)
                );
            }

            let derived = derived_diagnostics(&rows);
            let onset = |pred: &dyn Fn(usize) -> bool| {
                (0..rows.len()).find(|&i| rows[i][0] > 1e-3 && pred(i))
            };
            match onset(&|i| derived[i].direct_urca_electron) {
                Some(i) => println!(
                    "  URCA direto (e): n_B = {:.2} n0, Y_p = {:.3}",
                    rows[i][0], derived[i].proton_fraction
                ),
                None => println!("  URCA direto (e): não ocorre"),
            }
            let onsets: Vec<String> = HYPERONS
                .iter()
                .enumerate()
                .filter_map(|(h, label)| {
                    onset(&|i| rows[i][7 + h] > 1e-6).map(|i| format!("{label} {:.2}", rows[i][0]))
                })
                .collect();
            println!("  início dos hyperons (n_B/n0): {}", onsets.join(", "));
            let p_c_max = max.central_pressure;
            let cs2_max = rows
                .iter()
                .zip(&derived)
                .filter(|(r, _)| r[2] <= p_c_max)
                .map(|(_, d)| d.sound_speed_squared)
                .fold(0.0, f64::max);
            println!("  c_s^2 máximo até o centro de M_max: {cs2_max:.3}");
        }
    }
    Ok(())
}

/// `report observations [constraints.csv]`: Nível 3, confronto de GM1, GM3 e
/// FSU2 (com e sem hyperons, B = 0) com os vínculos observacionais. Imprime a
/// tabela e grava results/observations_report.csv.
pub fn observations(raw: &[String]) -> Result<(), String> {
    let args = Args::parse(raw, &[])?;
    let path = args
        .positional
        .first()
        .map(String::as_str)
        .unwrap_or("input/observations/constraints.csv");
    let constraints = load_constraints(path).map_err(|e| format!("falha ao ler '{path}': {e}"))?;

    create_dir("results")?;
    let mut csv = std::fs::File::create("results/observations_report.csv")
        .map_err(|e| format!("results/observations_report.csv: {e}"))?;
    writeln!(csv, "model,hyperons,constraint,model_value,value,err_minus,err_plus,credibility,distance_sigma,status,reference,doi")
        .unwrap();

    println!("Critério: d <= 1 compatível, 1 < d <= 2 tensão, d > 2 excluído (d em desvios-padrão).");
    for (name, model) in [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)] {
        let saturation = saturation_properties(model);
        for hyperons in [true, false] {
            let engine = HadronsMatter::new(model, 0.0)
                .with_hyperons(hyperons)
                .with_limits(0.02, 3.0)
                .with_points(2000);
            let rows = Solver::new(EngineMode::Hadrons(engine)).solve();
            let e: Vec<f64> = rows.iter().map(|r| r[1]).collect();
            let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
            let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
            let stars = generate_star_sequence(&e, &p, &n, true);

            let tag = if hyperons { "com hyperons" } else { "só núcleons" };
            println!("\n== {name} ({tag})");
            let (mut ok, mut tension, mut excluded) = (0, 0, 0);
            for c in &constraints {
                match assess(c, &stars, saturation.as_ref()) {
                    Some(a) => {
                        match a.status {
                            Status::Compatible => ok += 1,
                            Status::Tension => tension += 1,
                            Status::Excluded => excluded += 1,
                        }
                        println!(
                            "  {:<34} modelo {:>9.4} | obs {:>8.3} (-{}, +{}) {:<10} d = {:>5.2}  {}",
                            a.label, a.model_value, c.value, c.err_minus, c.err_plus,
                            c.credibility.replace(" sigma", "σ"), a.distance, a.status.label()
                        );
                        writeln!(
                            csv,
                            "{name},{hyperons},\"{}\",{:.6},{},{},{},\"{}\",{:.4},{},\"{}\",{}",
                            a.label, a.model_value, c.value, c.err_minus, c.err_plus,
                            c.credibility, a.distance, a.status.label(), c.reference, c.doi
                        )
                        .unwrap();
                    }
                    None => println!("  {:<34} não comparável (massa acima da máxima do modelo)", c.label),
                }
            }
            println!("  resumo: {ok} compatíveis, {tension} em tensão, {excluded} excluídos");
        }
    }
    println!("\nRelatório gravado em results/observations_report.csv");
    Ok(())
}
