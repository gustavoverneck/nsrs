// src/tov_solver.rs
use std::f64::consts::PI;

// Constantes de Conversão
const MEV_FM3_TO_MSUN_KM3: f64 = 8.9653e-7;
const G_C2: f64 = 1.4766; 

// A mesma função, pura e rápida (equivalente exato ao scipy interp1d com extrapolate)
fn interpolate_eps(p: f64, p_array: &[f64], eps_array: &[f64]) -> f64 {
    if p <= p_array[0] {
        let slope = (eps_array[1] - eps_array[0]) / (p_array[1] - p_array[0]);
        return eps_array[0] + slope * (p - p_array[0]);
    }
    
    let n = p_array.len();
    if p >= p_array[n - 1] {
        let slope = (eps_array[n - 1] - eps_array[n - 2]) / (p_array[n - 1] - p_array[n - 2]);
        return eps_array[n - 1] + slope * (p - p_array[n - 1]);
    }

    let idx = p_array.binary_search_by(|v| v.partial_cmp(&p).unwrap()).unwrap_or_else(|x| x);
    let i = idx.max(1) - 1;
    
    let t = (p - p_array[i]) / (p_array[i + 1] - p_array[i]);
    eps_array[i] + t * (eps_array[i + 1] - eps_array[i])
}

// Equações Diferenciais TOV
fn tov_derivatives(r: f64, p: f64, m: f64, p_array: &[f64], eps_array: &[f64]) -> (f64, f64) {
    let eps = interpolate_eps(p, p_array, eps_array);
    let num = (eps + p) * (m + 4.0 * PI * r.powi(3) * p);
    let den = r * (r - 2.0 * G_C2 * m);
    
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

    // Passo fixo de 10 metros (0.01 km) para garantir precisão absurda na superfície
    let dr = 0.01; 
    let mut r = 1e-5;
    let mut p = pc_tov;
    let mut m = 0.0;

    let r_end = 30000.0;

    while r < r_end {
        if p <= p_min {
            break;
        }

        // Integrador Runge-Kutta de 4ª Ordem (Manual e Transparente)
        let (k1_p, k1_m) = tov_derivatives(r, p, m, &p_tov, &eps_tov);
        let (k2_p, k2_m) = tov_derivatives(r + dr / 2.0, p + k1_p * dr / 2.0, m + k1_m * dr / 2.0, &p_tov, &eps_tov);
        let (k3_p, k3_m) = tov_derivatives(r + dr / 2.0, p + k2_p * dr / 2.0, m + k2_m * dr / 2.0, &p_tov, &eps_tov);
        let (k4_p, k4_m) = tov_derivatives(r + dr, p + k3_p * dr, m + k3_m * dr, &p_tov, &eps_tov);

        p += (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) * dr / 6.0;
        m += (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m) * dr / 6.0;
        r += dr;
    }

    Some((m, r))
}

pub fn generate_mr_curve(eps_array: &[f64], p_array: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut masses = Vec::new();
    let mut radii = Vec::new();
    
    let central_pressures = if p_array.len() > 5 { &p_array[5..] } else { p_array };

    for &pc in central_pressures {
        if let Some((m, r)) = integrate_star(pc, eps_array, p_array) {
            if m > 0.05 && r > 2.0 {
                masses.push(m);
                radii.push(r);
            }
        }
    }

    (masses, radii)
}