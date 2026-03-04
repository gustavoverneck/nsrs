mod core {
    pub mod constants;
    pub mod model;
    pub mod physics;
    pub mod particles;
    pub mod eos;
    pub mod solver;
    pub mod tov_solver;
    pub mod io_utils;
    pub mod plotting;
}

use core::model::{GM1, GM3};
use core::physics::PhysicsEngine;
use core::solver::Solver;
use core::io_utils::read_eos_file;
use core::plotting::plot_mr_curve;
use core::tov_solver::generate_mr_curve;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Uso: {} <B> <GM1|GM3>", args[0]);
        return;
    }
    let bg: f64 = args[1].parse().expect("B deve ser número");
    let model_name = &args[2];

    let model_params = match model_name.as_str() {
        "GM1" => GM1,
        "GM3" => GM3,
        _ => panic!("Modelo não reconhecido. Use GM1 ou GM3."),
    };

    let engine = PhysicsEngine::new(model_params, bg)
        .with_limits(0.01, 2.5)
        .with_points(4000);

    let mut solver = Solver::new(engine);
    let results = solver.solve();

    let eos_path = "output/eos.dat";
    if let Err(e) = Solver::write_eos(&results, eos_path) {
        eprintln!("Erro ao escrever arquivo: {}", e);
        return;
    }

    // Gera curva Massa-Raio a partir da EOS recém-calculada
    match read_eos_file(eos_path) {
        Ok((eps, p)) => {
            let (masses, radii) = generate_mr_curve(&eps, &p);
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