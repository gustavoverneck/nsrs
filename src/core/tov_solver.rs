// src/core/tov_solver.rs
use std::f64::consts::PI;
use rgsl::{Interp, InterpAccel, InterpType, interpolation}; // Importado 'interpolation'

// Constantes de Conversão
const MEV_FM3_TO_MSUN_KM3: f64 = 8.9653e-7;
const G_C2: f64 = 1.4766; 

// Equações Diferenciais TOV
fn tov_derivatives(
    r: f64, 
    p: f64, 
    m: f64, 
    p_array: &[f64], 
    eps_array: &[f64],
    spline: &Interp,
    accel: &mut InterpAccel
) -> (f64, f64) {
    
    // Interpolação via GSL usando a função do módulo interpolation
    let eps = if p <= p_array[0] {
        eps_array[0]
    } else if p >= p_array[p_array.len() - 1] {
        eps_array[p_array.len() - 1]
    } else {
        // Na rgsl, usa-se a função interpolation::eval diretamente
        interpolation::eval(spline, p_array, eps_array, p, accel)
    };

    let num = (eps + p) * (m + 4.0 * PI * r.powi(3) * p);
    let den = r * (r - 2.0 * G_C2 * m);
    
    if den <= 0.0 {
        return (f64::NEG_INFINITY, 0.0);
    }
    
    let dp_dr = -G_C2 * num / den; 
    let dm_dr = 4.0 * PI * r.powi(2) * eps; 
    
    (dp_dr, dm_dr)
}

pub fn integrate_star(
    pc_mev: f64, 
    eps_array: &[f64], 
    p_array: &[f64]
) -> Option<(f64, f64)> {
    
    let eps_tov: Vec<f64> = eps_array.iter().map(|&e| e * MEV_FM3_TO_MSUN_KM3).collect();
    let p_tov: Vec<f64> = p_array.iter().map(|&p| p * MEV_FM3_TO_MSUN_KM3).collect();
    let pc_tov = pc_mev * MEV_FM3_TO_MSUN_KM3;
    let p_min = p_tov[0];

    // Inicializa GSL: removido o .ok() pois Interp::new já retorna Option
    let mut accel = InterpAccel::new();
    let mut spline = Interp::new(InterpType::akima(), p_tov.len())?;
    spline.init(&p_tov, &eps_tov);

    let dr = 0.01; 
    let mut r = 1e-5;
    let mut p = pc_tov;
    let mut m = 0.0;
    let r_end = 30000.0;

    while r < r_end {
        if p.is_nan() || p <= p_min || r <= 2.0 * G_C2 * m {
            break;
        }

        let (k1_p, k1_m) = tov_derivatives(r, p, m, &p_tov, &eps_tov, &spline, &mut accel);
        let (k2_p, k2_m) = tov_derivatives(r + dr / 2.0, p + k1_p * dr / 2.0, m + k1_m * dr / 2.0, &p_tov, &eps_tov, &spline, &mut accel);
        let (k3_p, k3_m) = tov_derivatives(r + dr / 2.0, p + k2_p * dr / 2.0, m + k2_m * dr / 2.0, &p_tov, &eps_tov, &spline, &mut accel);
        let (k4_p, k4_m) = tov_derivatives(r + dr, p + k3_p * dr, m + k3_m * dr, &p_tov, &eps_tov, &spline, &mut accel);

        p += (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) * dr / 6.0;
        m += (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m) * dr / 6.0;
        r += dr;
    }

    Some((m, r))
}

pub fn generate_mr_curve(eps_array: &[f64], p_array: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut masses = Vec::new();
    let mut radii = Vec::new();

    // 1. Limpamos e ordenamos a EoS antes de começar qualquer integração
    let (clean_eps, clean_p) = clean_eos(eps_array, p_array);

    if clean_p.len() < 3 {
        return (masses, radii); // Retorna vazio se a EoS for insuficiente
    }

    // 2. Iteramos sobre as pressões centrais da EoS limpa
    for &pc in &clean_p {
        if let Some((m, r)) = integrate_star(pc, &clean_eps, &clean_p) {
            if m > 0.05 && r > 2.0 && m.is_finite() && r.is_finite() {
                masses.push(m);
                radii.push(r);
            }
        }
    }

    (masses, radii)
}

/// Garante que a pressão seja estritamente crescente para a GSL
fn clean_eos(eps: &[f64], p: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut combined: Vec<(f64, f64)> = p.iter().cloned().zip(eps.iter().cloned()).collect();

    // Ordena por pressão crescente
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut safe_p = Vec::with_capacity(combined.len());
    let mut safe_eps = Vec::with_capacity(combined.len());

    let mut last_p = -f64::INFINITY;

    for (pres, energy) in combined {
        // Regra de Ouro da GSL: pres deve ser estritamente maior que o anterior
        // Usamos um epsilon de 1e-18 para evitar ruídos de ponto flutuante
        if pres > last_p + 1e-18 && pres.is_finite() && energy.is_finite() {
            safe_p.push(pres);
            safe_eps.push(energy);
            last_p = pres;
        }
    }

    (safe_eps, safe_p)
}