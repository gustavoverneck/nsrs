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
    let _ = spline.init(&p_tov, &eps_tov);

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

/// Unifica a crosta personalizada (1/fm⁴) com a EoS do núcleo, descartando dados inválidos
pub fn unify_with_crust(core_eps: &[f64], core_p: &[f64]) -> (Vec<f64>, Vec<f64>) {
    // Constante de conversão de 1/fm⁴ para MeV/fm³
    const HBARC: f64 = 197.3269804;

    // Dados da crosta em 1/fm⁴ - from https://github.com/mrpelicer/nuclear_physics
    const CRUST_P_FM4: &[f64] = &[
        1.212e-11, 8.236e-11, 2.764e-10, 5.152e-10, 1.593e-09, 4.023e-09, 1.380e-08,
        3.315e-08, 1.077e-07, 2.559e-07, 3.479e-07, 4.729e-07, 6.430e-07, 9.147e-07,
        1.041e-06, 1.840e-06, 2.469e-06, 2.496e-06, 2.642e-06, 2.878e-06, 3.110e-06,
        3.425e-06, 3.852e-06, 4.425e-06, 5.181e-06, 6.168e-06, 8.198e-06, 1.109e-05,
        1.509e-05, 2.050e-05, 2.767e-05, 3.701e-05, 5.361e-05
    ];

    const CRUST_E_FM4: &[f64] = &[
        9.387e-08, 3.738e-07, 9.392e-07, 1.489e-06, 3.741e-06, 7.465e-06, 1.877e-05,
        3.747e-05, 9.418e-05, 1.881e-04, 2.369e-04, 2.982e-04, 3.758e-04, 5.242e-04,
        5.958e-04, 9.452e-04, 1.222e-03, 1.268e-03, 1.486e-03, 1.879e-03, 2.264e-03,
        2.765e-03, 3.400e-03, 4.182e-03, 5.131e-03, 6.260e-03, 8.329e-03, 1.090e-02,
        1.402e-02, 1.776e-02, 2.218e-02, 2.732e-02, 3.542e-02
    ];

    let mut raw_eps = Vec::with_capacity(CRUST_P_FM4.len() + core_p.len());
    let mut raw_p = Vec::with_capacity(CRUST_P_FM4.len() + core_p.len());

    if core_p.is_empty() {
        return (raw_eps, raw_p);
    }

    // 1. Inserir a Crosta (aplicando a conversão para MeV/fm³)
    for i in 0..CRUST_P_FM4.len() {
        raw_p.push(CRUST_P_FM4[i] * HBARC);
        raw_eps.push(CRUST_E_FM4[i] * HBARC);
    }

    // O ponto de transição agora é dinamicamente o último valor da sua crosta
    let p_transition = raw_p.last().copied().unwrap_or(0.0);
    let e_transition = raw_eps.last().copied().unwrap_or(0.0);

    // 2. Inserir o Núcleo (GM1/GM3)
    for i in 0..core_p.len() {
        // A costura só ocorre quando o núcleo supera tanto a pressão quanto a 
        // densidade de energia máximas da crosta. Isso preserva (e_c, p_c) como a fronteira absoluta.
        if core_p[i] > p_transition && core_eps[i] > e_transition {
            raw_p.push(core_p[i]);
            raw_eps.push(core_eps[i]);
        }
    }

    // 3. FILTRO DE MONOTONIA ESTRITA (Garante compatibilidade com a GSL)
    let mut combined: Vec<(f64, f64)> = raw_p.into_iter().zip(raw_eps.into_iter()).collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut final_eps = Vec::with_capacity(combined.len());
    let mut final_p = Vec::with_capacity(combined.len());
    let mut last_p = -1.0;
    let mut last_eps = -1.0;

    for (p, eps) in combined {
        // A física exige que P e EPS cresçam juntos estritamente
        if p > last_p && eps > last_eps {
            final_p.push(p);
            final_eps.push(eps);
            last_p = p;
            last_eps = eps;
        }
    }

    (final_eps, final_p)
}

pub fn generate_mr_curve(eps_array: &[f64], p_array: &[f64], with_crust: bool) -> (Vec<f64>, Vec<f64>) {
    let mut masses = Vec::new();
    let mut radii = Vec::new();

    // 1. Costura a crosta APENAS se a flag for verdadeira
    let (working_eps, working_p) = if with_crust {
        unify_with_crust(eps_array, p_array)
    } else {
        // Se falso, usa apenas a EoS original do núcleo (GM1/GM3)
        (eps_array.to_vec(), p_array.to_vec())
    };

    // 2. Limpa e ordena para satisfazer a GSL
    let (clean_eps, clean_p) = clean_eos(&working_eps, &working_p);

    if clean_p.len() < 3 {
        return (masses, radii); 
    }

    // 3. Define onde começar a iterar as pressões centrais
    // Se tiver crosta, pulamos os pontos de baixa pressão para não criar estrelas "ocas".
    let core_start_idx = if with_crust && !p_array.is_empty() {
        clean_p.iter().position(|&p| p >= p_array[0]).unwrap_or(0)
    } else {
        0
    };

    for &pc in &clean_p[core_start_idx..] {
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