use std::error::Error;
use std::time::Instant;
use nsrs::solver::{io_utils, tov_solver, plotting};

// O 'Result<(), Box<dyn Error>>' permite usar o operador '?' para tratamento de erro fácil no main
fn main() -> Result<(), Box<dyn Error>> {
    println!("--- Iniciando Resolvedor TOV em Rust ---");

    // 1. Caminhos dos arquivos
    let input_eos_path = "input/eos.dat";
    let output_plot_path = "output/mr_curve.svg";

    // 2. Leitura dos Dados
    println!("Lendo EoS de '{}'...", input_eos_path);
    // O operador '?' vai encerrar o programa e mostrar o erro se falhar aqui
    let (eps_data, p_data) = io_utils::read_eos_file(input_eos_path)
        .map_err(|e| format!("Falha ao ler arquivo EoS: {}", e))?;

    // 3. Cálculo da Curva M-R
    println!("Iniciando integração das equações TOV...");
    let start_time = Instant::now();
    
    let (masses, radii) = tov_solver::generate_mr_curve(&eps_data, &p_data);
    
    let duration = start_time.elapsed();
    println!("Cálculo concluído em {:.2?}s.", duration);
    println!("Modelos estelares gerados: {}", masses.len());

    if masses.is_empty() {
        return Err("Nenhuma estrela estável foi gerada. Verifique sua EoS.".into());
    }

    // Encontra a massa máxima aproximada
    let max_mass = masses.iter().fold(0./0., |a: f64, b| a.max(*b));
    println!("Massa Máxima aproximada: {:.2} M☉", max_mass);

    // 4. Geração do Gráfico Científico
    println!("Gerando gráfico em '{}'...", output_plot_path);
    plotting::plot_mr_curve(&radii, &masses, output_plot_path)
        .map_err(|e| format!("Falha ao plotar gráfico: {}", e))?;

    println!("--- Processo Finalizado com Sucesso ---");
    Ok(())
}