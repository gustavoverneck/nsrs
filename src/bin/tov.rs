use std::env;

use nsrs::core::io_utils::read_eos_file;
use nsrs::core::tov_solver::generate_mr_curve;
use nsrs::core::plotting::plot_mr_curve;

fn main() {
    // 1. Captura os argumentos de linha de comando
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Uso: {} <caminho_para_arquivo_eos>", args[0]);
        eprintln!("Exemplo: cargo run --bin tov eos.dat");
        std::process::exit(1);
    }

    let eos_path = &args[1];
    println!("Lendo arquivo EOS: {}", eos_path);

    // 2. Lê os dados da EoS (Densidade de Energia e Pressão)
    match read_eos_file(eos_path) {
        Ok((eps, p)) => {
            println!("EoS lida com sucesso. Pontos carregados: {}", eps.len());
            println!("Integrando equações TOV...");

            // 3. Gera a curva Massa-Raio
            let (masses, radii) = generate_mr_curve(&eps, &p);

            if masses.is_empty() || radii.is_empty() {
                eprintln!("Erro: A curva M-R retornou vazia. Verifique se as pressões da EoS suportam uma estrela.");
            } else {
                // 4. Salva o gráfico extraindo o nome base do arquivo EoS
                let output_name = eos_path
                    .split('/')
                    .last()
                    .unwrap_or("eos")
                    .replace(".dat", "");
                
                let plot_path = format!("results/mr_{}.svg", output_name);

                println!("Gerando gráfico em: {}", plot_path);
                if let Err(e) = plot_mr_curve(&radii, &masses, &plot_path) {
                    eprintln!("Erro ao salvar o gráfico SVG: {}", e);
                } else {
                    println!("Concluído! Curva M-R gerada com sucesso.");
                    
                    // Exibe a massa máxima atingida no terminal
                    let max_mass = masses.iter().copied().fold(f64::NAN, f64::max);
                    println!("Massa máxima da estrela: {:.2} M_sol", max_mass);
                }
            }
        }
        Err(e) => {
            eprintln!("Falha ao ler o arquivo '{}': {}", eos_path, e);
        }
    }
}