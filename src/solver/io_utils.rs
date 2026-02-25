use std::error::Error;
use std::fs;
use std::path::Path;

pub const HBAR_C: f64 = 197.3269804; // MeV * fm

pub fn read_eos_file<P: AsRef<Path>>(path: P) -> Result<(Vec<f64>, Vec<f64>), Box<dyn Error>> {
    let content = fs::read_to_string(path)?;
    
    let mut eps_data = Vec::new();
    let mut p_data = Vec::new();

    for line in content.lines() {
        let trimmed_line = line.trim();
        
        if trimmed_line.is_empty() || trimmed_line.starts_with('#') {
            continue;
        }

        let columns: Vec<&str> = trimmed_line.split_whitespace().collect();
        
        // Agora esperamos pelo menos 3 colunas: n_B (0), eps (1), P (2)
        if columns.len() >= 3 {
            // Lemos direto das colunas 1 e 2. A coluna 0 (n_B) é ignorada.
            let eps_fm4: f64 = columns[1].parse().map_err(|_| format!("Erro ao ler: {}", columns[1]))?;
            let p_fm4: f64 = columns[2].parse().map_err(|_| format!("Erro ao ler: {}", columns[2]))?;
            
            // Conversão de unidades [fm^-4] -> [MeV/fm^3]
            let eps = eps_fm4 * HBAR_C;
            let p = p_fm4 * HBAR_C;

            // Filtro básico de sanidade
            if eps > 0.0 && p >= 0.0 {
                eps_data.push(eps);
                p_data.push(p);
            }
        }
    }

    if eps_data.is_empty() {
        return Err("Arquivo vazio ou formato incorreto (precisa de 3 colunas).".into());
    }
    
    // O seu script Python verifica se precisa inverter a ordem.
    // Nossa função sort_eos_data garante que sempre ficará em ordem crescente de Pressão,
    // o que é a maneira mais segura e robusta de preparar os dados para interpolação.
    sort_eos_data(&mut eps_data, &mut p_data);

    Ok((eps_data, p_data))
}

fn sort_eos_data(eps: &mut Vec<f64>, p: &mut Vec<f64>) {
    let mut combined: Vec<(f64, f64)> = p.iter().cloned().zip(eps.iter().cloned()).collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    *p = combined.iter().map(|x| x.0).collect();
    *eps = combined.iter().map(|x| x.1).collect();
}