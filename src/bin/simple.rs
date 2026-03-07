use nsrs::core::model::{GM1, GM3};
use nsrs::core::physics::PhysicsEngine;
use nsrs::core::solver::Solver;
use nsrs::core::io_utils::read_eos_file;
use nsrs::core::plotting::plot_mr_curve;
use nsrs::core::tov_solver::generate_mr_curve;
use std::env;

// usage:  cargo run --release --bin simple_ns -- <B> <model>
// example: cargo run --release --bin simple_ns -- 1e15 GM1

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Uso: {} <GM1|GM3> <B>", args[0]);
        return;
    }
    let model_name = &args[1];
    let bg: f64 = args[2].parse().expect("B deve ser número");

    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let engine = PhysicsEngine::new(model_params, bg)
        .with_limits(0.02, 5.0)
        .with_points(2003);

    let mut solver = Solver::new(engine);
    let results = solver.solve();

    let eos_path = "eos.dat";
    if let Err(e) = Solver::write_eos(&results, eos_path) {
        eprintln!("Erro ao escrever arquivo: {}", e);
        return;
    }

    // Gera curva Massa-Raio a partir da EOS recém-calculada
    match read_eos_file(eos_path) {
        Ok((eps, p)) => {
            let (masses, radii) = generate_mr_curve(&eps, &p, false);
            if masses.is_empty() || radii.is_empty() {
                eprintln!("Curva M-R vazia; verifique a EOS gerada.");
            } else {
                let mr_path = format!("results/mr_{}_B{:.3e}.svg", model_name, bg);
                if let Err(e) = plot_mr_curve(&radii, &masses, &mr_path) {
                    eprintln!("Erro ao plotar M-R: {}", e);
                } else {
                    println!("Curva M-R salva em {}", mr_path);
                }
            }
        }
        Err(e) => eprintln!("Erro ao ler EOS para TOV: {}", e),
    }
    
}