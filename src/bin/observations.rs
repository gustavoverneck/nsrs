// src/bin/observations.rs
//
// Nível 3: confronto de GM1, GM3 e FSU2 (com e sem hyperons, B = 0) com os
// vínculos de input/observations/constraints.csv. Imprime a tabela e grava
// results/observations_report.csv.
//
// Uso: cargo run --release --bin observations [constraints.csv]

use std::io::Write;

use nsrs::core::nuclear::saturation_properties;
use nsrs::core::observations::{Status, assess, load_constraints};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "input/observations/constraints.csv".to_string());
    let constraints = load_constraints(&path).unwrap_or_else(|e| {
        eprintln!("falha ao ler '{path}': {e}");
        std::process::exit(1);
    });

    std::fs::create_dir_all("results").expect("cannot create results/");
    let mut csv = std::fs::File::create("results/observations_report.csv").expect("report file");
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
}
