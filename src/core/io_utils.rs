use std::error::Error;
use std::fs;
use std::io::Write;
use std::path::Path;

use crate::core::constants::RESULTS_SIZE;

pub fn read_eos_file<P: AsRef<Path>>(
    path: P,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), Box<dyn Error>> {
    let content = fs::read_to_string(path)?;

    let mut rho_data = Vec::new();
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
            let rho: f64 = columns[0]
                .parse()
                .map_err(|_| format!("Erro ao ler: {}", columns[0]))?;
            let eps: f64 = columns[1]
                .parse()
                .map_err(|_| format!("Erro ao ler: {}", columns[1]))?;
            let p: f64 = columns[2]
                .parse()
                .map_err(|_| format!("Erro ao ler: {}", columns[2]))?;

            // Filtro básico de sanidade
            if eps > 0.0 && p >= 0.0 {
                rho_data.push(rho);
                eps_data.push(eps);
                p_data.push(p);
            }
        }
    }

    if eps_data.is_empty() {
        return Err("Arquivo vazio ou formato incorreto (precisa de 3 colunas).".into());
    }

    sort_eos_data(&mut rho_data, &mut eps_data, &mut p_data);

    Ok((eps_data, p_data, rho_data))
}

fn sort_eos_data(rho: &mut Vec<f64>, eps: &mut Vec<f64>, p: &mut Vec<f64>) {
    let mut combined: Vec<(f64, f64, f64)> = p
        .iter()
        .cloned()
        .zip(eps.iter().cloned())
        .zip(rho.iter().cloned())
        .map(|((p_val, eps_val), rho_val)| (p_val, eps_val, rho_val))
        .collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    *p = combined.iter().map(|x| x.0).collect();
    *eps = combined.iter().map(|x| x.1).collect();
    *rho = combined.iter().map(|x| x.2).collect();
}

pub fn write_eos_with_mr<P: AsRef<Path>>(
    results: &[[f64; RESULTS_SIZE]],
    masses: &[f64],
    radii: &[f64],
    baryonic_masses: &[f64],
    central_pressures: &[f64],
    path: P,
) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;

    writeln!(
        file,
        "# 0:nB_over_n0 1:eps_MeV_fm3 2:p_MeV_fm3 3:ne_fm^-3 4:nmu_fm^-3 \
5:nn_fm^-3 6:np_fm^-3 7:nL0_fm^-3 8:nS-_fm^-3 9:nS0_fm^-3 \
10:nS+_fm^-3 11:nX-_fm^-3 12:nX0_fm^-3 13:g_sigma_sigma_MeV \
14:g_omega_omega0_MeV 15:g_rho_rho03_MeV 16:mstar_over_mN \
17:mu_n_over_mN 18:mu_e_over_mN \
19:eps_mag_MeV_fm3 20:mu_total_over_mN 21:n_chi_fm^-3 22:Y_chi \
23:m_chi_MeV 24:m_x_MeV 25:epsilon 26:g_d 27:X0_MeV 28:kF_chi_MeV \
29:mu_chi_MeV 30:eps_chi_kin_MeV_fm3 31:P_chi_kin_MeV_fm3 \
32:eps_X_MeV_fm3 33:P_X_MeV_fm3 \
34:mr_mass_msun 35:mr_radius_km 36:mr_baryonic_mass_msun"
    )?;

    let n_mr = masses
        .len()
        .min(radii.len())
        .min(baryonic_masses.len())
        .min(central_pressures.len());

    for row in results.iter() {
        let p_row = row[2];

        let mut matched_mass = f64::NAN;
        let mut matched_radius = f64::NAN;
        let mut matched_baryonic_mass = f64::NAN;
        let mut best_diff = f64::INFINITY;

        for i in 0..n_mr {
            let p_c = central_pressures[i];
            let diff = (p_c - p_row).abs();
            let tol = p_row.abs().max(1.0) * 1e-8;

            if diff <= tol && diff < best_diff {
                best_diff = diff;
                matched_mass = masses[i];
                matched_radius = radii[i];
                matched_baryonic_mass = baryonic_masses[i];
            }
        }

        // Skip rows with NaN values
        if matched_mass.is_nan()
            || matched_radius.is_nan()
            || matched_baryonic_mass.is_nan()
            || row.iter().any(|v| v.is_nan())
        {
            continue;
        }

        let base_line = row
            .iter()
            .map(|val| format!("{:12.5e}", val))
            .collect::<Vec<String>>()
            .join(" ");

        writeln!(
            file,
            "{} {:12.5e} {:12.5e} {:12.5e}",
            base_line, matched_mass, matched_radius, matched_baryonic_mass
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::constants::DATA_SIZE;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn eos_writer_preserves_the_34_plus_3_column_contract() {
        let mut row = [0.0; RESULTS_SIZE];
        row[0] = 1.0;
        row[1] = 150.0;
        row[2] = 25.0;
        row[21] = 0.001;
        row[22] = 0.01;
        row[23] = 500.0;
        row[24] = 25.0;

        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock must follow Unix epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "nsrs_eos_contract_{}_{}.dat",
            std::process::id(),
            unique
        ));

        write_eos_with_mr(&[row], &[1.4], &[12.0], &[1.5], &[25.0], &path)
            .expect("contract fixture must be writable");
        let content = fs::read_to_string(&path).expect("contract fixture must be readable");
        let _ = fs::remove_file(&path);

        let mut lines = content.lines();
        let header = lines.next().expect("writer must emit a header");
        assert!(header.contains("33:P_X_MeV_fm3"));
        assert!(header.contains("34:mr_mass_msun"));
        assert!(header.contains("35:mr_radius_km"));
        assert!(header.contains("36:mr_baryonic_mass_msun"));

        let values: Vec<f64> = lines
            .next()
            .expect("writer must emit the matched EOS row")
            .split_whitespace()
            .map(|token| token.parse().expect("every EOS token must be numeric"))
            .collect();
        assert_eq!(values.len(), DATA_SIZE);
        assert_eq!(values[RESULTS_SIZE], 1.4);
        assert_eq!(values[RESULTS_SIZE + 1], 12.0);
        assert_eq!(values[RESULTS_SIZE + 2], 1.5);
    }
}
