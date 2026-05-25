// src/core/tov_solver.rs
use std::f64::consts::PI;

use crate::constants::{M_NUCLEON, RESULTS_SIZE, MEV_FM3_TO_MSUN_KM3, G_C2};

fn log_linear_interp(x: &[f64], y: &[f64], xval: f64) -> f64 {
    let n = x.len();
    if n < 2 || y.len() != n || !xval.is_finite() {
        return f64::NAN;
    }

    if xval <= x[0] {
        return y[0];
    }
    if xval >= x[n - 1] {
        return y[n - 1];
    }

    let mut klo = 0usize;
    let mut khi = n - 1;
    while khi - klo > 1 {
        let k = (khi + klo) >> 1;
        if x[k] > xval {
            khi = k;
        } else {
            klo = k;
        }
    }

    let x0 = x[klo];
    let x1 = x[khi];
    let h = x1 - x0;

    if h == 0.0 {
        return f64::NAN;
    }

    let y0 = y[klo];
    let y1 = y[khi];

    if !y0.is_finite() || !y1.is_finite() {
        return f64::NAN;
    }

    // Use log-log interpolation only when all values are strictly positive.
    // This avoids NaN from ln(0) or ln(negative) and keeps the GSL input safe.
    if x0 > 0.0 && x1 > 0.0 && xval > 0.0 && y0 > 0.0 && y1 > 0.0 {
        let log_x0 = x0.ln();
        let log_x1 = x1.ln();
        let log_xval = xval.ln();

        let log_y0 = y0.ln();
        let log_y1 = y1.ln();

        let t_log = (log_xval - log_x0) / (log_x1 - log_x0);
        (log_y0 + t_log * (log_y1 - log_y0)).exp()
    } else {
        // Fallback: safe linear interpolation for zero/negative values,
        // such as near the stellar surface.
        let t = (xval - x0) / h;
        y0 + t * (y1 - y0)
    }
}

// Equações Diferenciais TOV (usando interpolação log-linear local)
fn tov_derivatives(
    r: f64,
    p: f64,
    m: f64,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> (f64, f64, f64) {

    let eps = log_linear_interp(p_array, eps_array, p);
    let rho = log_linear_interp(p_array, rho_array, p);

    let num = (eps + p) * (m + 4.0 * PI * r.powi(3) * p);
    let den = r * (r - 2.0 * G_C2 * m);

    let metric_term = 1.0 - 2.0 * G_C2 * m / r;

    if den <= 0.0 || metric_term <= 0.0 {
        return (f64::NEG_INFINITY, 0.0, 0.0);
    }

    let dp_dr = -G_C2 * num / den;
    let dm_dr = 4.0 * PI * r.powi(2) * eps;
    let dmb_dr = 4.0 * PI * r.powi(2) * rho / metric_term.sqrt();

    (dp_dr, dm_dr, dmb_dr)
}

fn rkck_step(
    r: f64,
    y: [f64; 3],
    h: f64,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> ([f64; 3], [f64; 3]) {
    let (k1_p, k1_m, k1_mb) = tov_derivatives(
        r,
        y[0],
        y[1],
        p_array,
        eps_array,
        rho_array,
    );

    let a2 = 0.2;
    let a3 = 0.3;
    let a4 = 0.6;
    let a5 = 1.0;
    let a6 = 0.875;

    let b21 = 0.2;
    let b31 = 3.0 / 40.0;
    let b32 = 9.0 / 40.0;
    let b41 = 0.3;
    let b42 = -0.9;
    let b43 = 1.2;
    let b51 = -11.0 / 54.0;
    let b52 = 2.5;
    let b53 = -70.0 / 27.0;
    let b54 = 35.0 / 27.0;
    let b61 = 1631.0 / 55296.0;
    let b62 = 175.0 / 512.0;
    let b63 = 575.0 / 13824.0;
    let b64 = 44275.0 / 110592.0;
    let b65 = 253.0 / 4096.0;

    let c1 = 37.0 / 378.0;
    let c3 = 250.0 / 621.0;
    let c4 = 125.0 / 594.0;
    let c6 = 512.0 / 1771.0;

    let dc1 = c1 - 2825.0 / 27648.0;
    let dc3 = c3 - 18575.0 / 48384.0;
    let dc4 = c4 - 13525.0 / 55296.0;
    let dc5 = -277.0 / 14336.0;
    let dc6 = c6 - 0.25;

    let y2 = [
        y[0] + h * b21 * k1_p,
        y[1] + h * b21 * k1_m,
        y[2] + h * b21 * k1_mb,
    ];
    let (k2_p, k2_m, k2_mb) = tov_derivatives(
        r + a2 * h,
        y2[0],
        y2[1],
        p_array,
        eps_array,
        rho_array,
    );

    let y3 = [
        y[0] + h * (b31 * k1_p + b32 * k2_p),
        y[1] + h * (b31 * k1_m + b32 * k2_m),
        y[2] + h * (b31 * k1_mb + b32 * k2_mb),
    ];
    let (k3_p, k3_m, k3_mb) = tov_derivatives(
        r + a3 * h,
        y3[0],
        y3[1],
        p_array,
        eps_array,
        rho_array,
    );

    let y4 = [
        y[0] + h * (b41 * k1_p + b42 * k2_p + b43 * k3_p),
        y[1] + h * (b41 * k1_m + b42 * k2_m + b43 * k3_m),
        y[2] + h * (b41 * k1_mb + b42 * k2_mb + b43 * k3_mb),
    ];
    let (k4_p, k4_m, k4_mb) = tov_derivatives(
        r + a4 * h,
        y4[0],
        y4[1],
        p_array,
        eps_array,
        rho_array,
    );

    let y5 = [
        y[0] + h * (b51 * k1_p + b52 * k2_p + b53 * k3_p + b54 * k4_p),
        y[1] + h * (b51 * k1_m + b52 * k2_m + b53 * k3_m + b54 * k4_m),
        y[2] + h * (b51 * k1_mb + b52 * k2_mb + b53 * k3_mb + b54 * k4_mb),
    ];
    let (k5_p, k5_m, k5_mb) = tov_derivatives(
        r + a5 * h,
        y5[0],
        y5[1],
        p_array,
        eps_array,
        rho_array,
    );

    let y6 = [
        y[0]
            + h
                * (b61 * k1_p + b62 * k2_p + b63 * k3_p + b64 * k4_p + b65 * k5_p),
        y[1]
            + h
                * (b61 * k1_m + b62 * k2_m + b63 * k3_m + b64 * k4_m + b65 * k5_m),
        y[2]
            + h
                * (b61 * k1_mb
                    + b62 * k2_mb
                    + b63 * k3_mb
                    + b64 * k4_mb
                    + b65 * k5_mb),
    ];
    let (k6_p, k6_m, k6_mb) = tov_derivatives(
        r + a6 * h,
        y6[0],
        y6[1],
        p_array,
        eps_array,
        rho_array,
    );

    let yout = [
        y[0] + h * (c1 * k1_p + c3 * k3_p + c4 * k4_p + c6 * k6_p),
        y[1] + h * (c1 * k1_m + c3 * k3_m + c4 * k4_m + c6 * k6_m),
        y[2] + h * (c1 * k1_mb + c3 * k3_mb + c4 * k4_mb + c6 * k6_mb),
    ];

    let yerr = [
        h * (dc1 * k1_p + dc3 * k3_p + dc4 * k4_p + dc5 * k5_p + dc6 * k6_p),
        h * (dc1 * k1_m + dc3 * k3_m + dc4 * k4_m + dc5 * k5_m + dc6 * k6_m),
        h * (dc1 * k1_mb + dc3 * k3_mb + dc4 * k4_mb + dc5 * k5_mb + dc6 * k6_mb),
    ];

    (yout, yerr)
}

fn rkqs_step(
    r: f64,
    y: [f64; 3],
    htry: f64,
    eps: f64,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> Option<([f64; 3], f64, f64)> {
    let safety = 0.9;
    let pgrow = -0.2;
    let pshrink = -0.25;
    let errcon = 1.89e-4;
    let tiny = 1.0e-30;

    let mut h = htry;
    loop {
        let (yout, yerr) = rkck_step(
            r,
            y,
            h,
            p_array,
            eps_array,
            rho_array,
        );

        let yscal = [
            y[0].abs() + (h * yerr[0]).abs() + tiny,
            y[1].abs() + (h * yerr[1]).abs() + tiny,
            y[2].abs() + (h * yerr[2]).abs() + tiny,
        ];
        let mut errmax: f64 = 0.0;
        for i in 0..3 {
            errmax = errmax.max((yerr[i] / yscal[i]).abs());
        }
        errmax /= eps;

        if errmax > 1.0 {
            let htemp = safety * h * errmax.powf(pshrink);
            let hnew = h.signum() * htemp.abs().max(0.1 * h.abs());
            if r + hnew == r {
                return None;
            }
            h = hnew;
            continue;
        }

        let hnext = if errmax > errcon {
            safety * h * errmax.powf(pgrow)
        } else {
            5.0 * h
        };

        return Some((yout, h, hnext));
    }
}

pub fn integrate_star(
    pc_tov: f64,
    p_min: f64,
    p_tov: &[f64],
    eps_tov: &[f64],
    rho_tov: &[f64],
) -> Option<(f64, f64, f64, f64)> {
    if p_tov.len() < 5 || eps_tov.len() < 5 || rho_tov.len() < 5 {
        return None;
    }
    if p_tov.len() != eps_tov.len() || p_tov.len() != rho_tov.len() {
        return None;
    }
    let mut r = 1e-5;
    let mut y = [pc_tov, 0.0, 0.0];
    let mut h = 1.0e-2;
    let hmin = 0.0;
    let r_end = 30000.0;
    let eps = 3.0e-16;
    let mut steps = 0u64;
    let max_steps = 90000u64;

    while r < r_end && steps < max_steps {
        if y[0].is_nan() || y[0] <= p_min || r <= 2.0 * G_C2 * y[1] {
            break;
        }

        if r + h > r_end {
            h = r_end - r;
        }

        let (ynew, hdid, hnext) = match rkqs_step(
            r,
            y,
            h,
            eps,
            p_tov,
            &eps_tov,
            &rho_tov,
        ) {
            Some(result) => result,
            None => break,
        };

        if ynew[0] <= p_min {
            y = ynew;
            r += hdid;
            break;
        }

        y = ynew;
        r += hdid;

        if hnext.abs() < hmin {
            break;
        }
        h = hnext;
        steps += 1;
    }

    Some((y[1], r, y[2], pc_tov / MEV_FM3_TO_MSUN_KM3))
}

/// Unifica a crosta personalizada (1/fm⁴) com a EoS do núcleo, descartando dados inválidos
pub fn unify_with_crust(core_eps: &[f64], core_p: &[f64], core_rho: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // Constante de conversão de 1/fm⁴ para MeV/fm³
    const HBARC: f64 = 197.3269804;

    // Dados da crosta em 1/fm⁴- Baym-Pethick-Sutherland - from https://github.com/mrpelicer/nuclear_physics
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

    // Fallback: use crust energy density as a rho proxy until a rho table is available.
    const CRUST_RHO_FM4: &[f64] = CRUST_E_FM4;

    let mut raw_eps = Vec::with_capacity(CRUST_P_FM4.len() + core_p.len());
    let mut raw_p = Vec::with_capacity(CRUST_P_FM4.len() + core_p.len());
    let mut raw_rho = Vec::with_capacity(CRUST_P_FM4.len() + core_p.len());

    if core_p.is_empty() || core_eps.len() != core_p.len() || core_rho.len() != core_p.len() {
        return (raw_eps, raw_p, raw_rho);
    }

    // 1. Inserir a Crosta (aplicando a conversão para MeV/fm³)
    for i in 0..CRUST_P_FM4.len() {
        raw_p.push(CRUST_P_FM4[i] * HBARC);
        raw_eps.push(CRUST_E_FM4[i] * HBARC);
        raw_rho.push(CRUST_RHO_FM4[i] * HBARC);
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
            
            // CONVERSÃO APLICADA AQUI: De 1/fm³ para MeV/fm³
            raw_rho.push(core_rho[i] * M_NUCLEON); 
        }
    }

    // 3. FILTRO DE MONOTONIA ESTRITA (Garante compatibilidade com a GSL)
    let mut combined: Vec<(f64, f64, f64)> = raw_p
        .into_iter()
        .zip(raw_eps.into_iter())
        .zip(raw_rho.into_iter())
        .map(|((p, eps), rho)| (p, eps, rho))
        .collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut final_eps = Vec::with_capacity(combined.len());
    let mut final_p = Vec::with_capacity(combined.len());
    let mut final_rho = Vec::with_capacity(combined.len());
    let mut last_p = -1.0;
    let mut last_eps = -1.0;

    for (p, eps, rho) in combined {
        // A física exige que P e EPS cresçam juntos estritamente
        if p > last_p && eps > last_eps && eps.is_finite() && rho.is_finite() {
            final_p.push(p);
            final_eps.push(eps);
            final_rho.push(rho);
            last_p = p;
            last_eps = eps;
        }
    }

    (final_eps, final_p, final_rho)
}

pub fn generate_mr_curve(eps_array: &[f64], p_array: &[f64], rho_array: &[f64], with_crust: bool) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut masses = Vec::new();
    let mut radii = Vec::new();
    let mut baryonic_masses = Vec::new();
    let mut central_pressures = Vec::new();

    // 1. Costura a crosta APENAS se a flag for verdadeira
    let (clean_eps, clean_p, clean_rho) = if with_crust {
        // A função unify_with_crust já faz a conversão do núcleo internamente agora
        unify_with_crust(eps_array, p_array, rho_array) 
    } else {
        // Se NÃO usar a crosta, precisamos converter a EoS original de 1/fm³ para MeV/fm³
        const NUCLEON_MASS_MEV: f64 = M_NUCLEON;
        let converted_rho: Vec<f64> = rho_array
            .iter()
            .map(|&nb| nb * NUCLEON_MASS_MEV)
            .collect();
        
        (eps_array.to_vec(), p_array.to_vec(), converted_rho)
    };

    // 2. Limpa e ordena para satisfazer a GSL
    let (clean_eps, clean_p, clean_rho) = clean_eos_with_rho(&clean_eps, &clean_p, &clean_rho);

    if clean_p.len() < 5 {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new()); 
    }

    // 3. Converte unidades uma unica vez
    let eps_tov: Vec<f64> = clean_eps.iter().map(|&e| e * MEV_FM3_TO_MSUN_KM3).collect();
    let p_tov: Vec<f64> = clean_p.iter().map(|&p| p * MEV_FM3_TO_MSUN_KM3).collect();
    let rho_tov: Vec<f64> = clean_rho.iter().map(|&rho| rho * MEV_FM3_TO_MSUN_KM3).collect();
    let p_min = p_tov[0];

    // 4. Define onde começar a iterar as pressões centrais
    // Se tiver crosta, pulamos os pontos de baixa pressão para não criar estrelas "ocas".
    let core_start_idx = if with_crust && !p_array.is_empty() {
        clean_p.iter().position(|&p| p >= p_array[0]).unwrap_or(0)
    } else {
        0
    };

    for &pc_mev in &clean_p[core_start_idx..] {
        let pc_tov = pc_mev * MEV_FM3_TO_MSUN_KM3;
        if let Some((m, r, mb, pc_final)) = integrate_star(
            pc_tov,
            p_min,
            &p_tov,
            &eps_tov,
            &rho_tov,
        ) {
            if m > 0.05 && r > 2.0 && m.is_finite() && r.is_finite() {
                masses.push(m);
                radii.push(r);
                baryonic_masses.push(mb);
                central_pressures.push(pc_final);
            }
        }
    }

    (masses, radii, baryonic_masses, central_pressures)
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

/// Procura e interpola os dados de quarks para um mun específico
fn find_quark_point(quark_eos: &[[f64; RESULTS_SIZE]], target_mun: f64) -> Option<[f64; RESULTS_SIZE]> {
    // Busca binária para encontrar a posição do mun (índice 17)
    let pos = quark_eos.binary_search_by(|row| {
        row[17].partial_cmp(&target_mun).unwrap_or(std::cmp::Ordering::Equal)
    });

    match pos {
        // Encontrou o valor exato
        Ok(idx) => Some(quark_eos[idx]),
        
        // Não encontrou valor exato, tenta interpolar entre idx-1 e idx
        Err(idx) => {
            if idx == 0 || idx >= quark_eos.len() {
                return None; // mun fora do range da tabela de quarks
            }
            
            let q_low = &quark_eos[idx - 1];
            let q_high = &quark_eos[idx];
            
            // Fator de interpolação linear
            let factor = (target_mun - q_low[17]) / (q_high[17] - q_low[17]);
            
            let mut interpolated = [0.0; RESULTS_SIZE];
            for i in 0..RESULTS_SIZE {
                interpolated[i] = q_low[i] + factor * (q_high[i] - q_low[i]);
            }
            Some(interpolated)
        }
    }
}

pub fn unify_hybrid_eos(hadron_eos: &[[f64; RESULTS_SIZE]], quark_eos: &[[f64; RESULTS_SIZE]]) -> (Vec<f64>, Vec<f64>) {
    let mut final_eps = Vec::new();
    let mut final_p = Vec::new();

    for h_row in hadron_eos {
        let mun = h_row[17]; 
        let p_had = h_row[2];

        // Agora a função existe e retorna os dados interpolados
        if let Some(q_row) = find_quark_point(quark_eos, mun) {
            let p_qrk = q_row[2];

            if p_qrk > p_had {
                // Construção de Maxwell: Fase de Quarks é mais estável
                final_eps.push(q_row[1]);
                final_p.push(q_row[2]);
            } else {
                // Fase Hadrónica ainda é mais estável
                final_eps.push(h_row[1]);
                final_p.push(h_row[2]);
            }
        } else {
            // Se o mun for muito baixo (ex: na crosta), a fase de quarks nem existe
            final_eps.push(h_row[1]);
            final_p.push(h_row[2]);
        }
    }
    (final_eps, final_p)
}

/// Garante que a pressão seja estritamente crescente para a GSL (com rho)
fn clean_eos_with_rho(eps: &[f64], p: &[f64], rho: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut combined: Vec<(f64, f64, f64)> = p
        .iter()
        .cloned()
        .zip(eps.iter().cloned())
        .zip(rho.iter().cloned())
        .map(|((p, eps), rho)| (p, eps, rho))
        .collect();

    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut safe_p = Vec::with_capacity(combined.len());
    let mut safe_eps = Vec::with_capacity(combined.len());
    let mut safe_rho = Vec::with_capacity(combined.len());

    let mut last_p = -f64::INFINITY;

    for (pres, energy, rho_val) in combined {
        if pres > last_p + 1e-18 && pres.is_finite() && energy.is_finite() && rho_val.is_finite() {
            safe_p.push(pres);
            safe_eps.push(energy);
            safe_rho.push(rho_val);
            last_p = pres;
        }
    }

    (safe_eps, safe_p, safe_rho)
}