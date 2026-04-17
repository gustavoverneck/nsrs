use std::error::Error;
use std::fs;
use std::io::Write;
use std::path::Path;

use crate::core::constants::RESULTS_SIZE;

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
            let eps: f64 = columns[1].parse().map_err(|_| format!("Erro ao ler: {}", columns[1]))?;
            let p: f64 = columns[2].parse().map_err(|_| format!("Erro ao ler: {}", columns[2]))?;
            
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
    
    sort_eos_data(&mut eps_data, &mut p_data);

    Ok((eps_data, p_data))
}

fn sort_eos_data(eps: &mut Vec<f64>, p: &mut Vec<f64>) {
    let mut combined: Vec<(f64, f64)> = p.iter().cloned().zip(eps.iter().cloned()).collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    *p = combined.iter().map(|x| x.0).collect();
    *eps = combined.iter().map(|x| x.1).collect();
}

pub fn write_eos_with_mr<P: AsRef<Path>>(
    results: &[[f64; RESULTS_SIZE]],
    masses: &[f64],
    radii: &[f64],
    central_pressures: &[f64],
    path: P,
) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;

    writeln!(
        file,
        "# nB eps p ... mr_mass_msun mr_radius_km"
    )?;

    let n_mr = masses.len().min(radii.len()).min(central_pressures.len());

    for row in results.iter() {
        let p_row = row[2];

        let mut matched_mass = f64::NAN;
        let mut matched_radius = f64::NAN;
        let mut best_diff = f64::INFINITY;

        for i in 0..n_mr {
            let p_c = central_pressures[i];
            let diff = (p_c - p_row).abs();
            let tol = p_row.abs().max(1.0) * 1e-8;

            if diff <= tol && diff < best_diff {
                best_diff = diff;
                matched_mass = masses[i];
                matched_radius = radii[i];
            }
        }

        // Skip rows with NaN values
        if matched_mass.is_nan() || matched_radius.is_nan() || row.iter().any(|v| v.is_nan()) {
            continue;
        }

        let base_line = row
            .iter()
            .map(|val| format!("{:12.5e}", val))
            .collect::<Vec<String>>()
            .join(" ");

        writeln!(
            file,
            "{} {:12.5e} {:12.5e}",
            base_line,
            matched_mass,
            matched_radius
        )?;
    }

    Ok(())
}