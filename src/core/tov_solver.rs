// src/core/tov_solver.rs
use std::f64::consts::PI;

use crate::constants::{G_C2, M_NUCLEON, MEV_FM3_TO_MSUN_KM3, N0, RESULTS_SIZE};

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

/// Interpolação log-log (ou linear, fora do domínio positivo) que devolve
/// também a derivada dy/dx no segmento.
fn log_linear_interp_with_slope(x: &[f64], y: &[f64], xval: f64) -> (f64, f64) {
    let n = x.len();
    if n < 2 || !xval.is_finite() || xval <= x[0] || xval >= x[n - 1] {
        return (log_linear_interp(x, y, xval), 0.0);
    }
    let k = x.partition_point(|&v| v <= xval).clamp(1, n - 1);
    let (x0, x1, y0, y1) = (x[k - 1], x[k], y[k - 1], y[k]);
    let value = log_linear_interp(x, y, xval);
    let slope = if x0 > 0.0 && x1 > 0.0 && y0 > 0.0 && y1 > 0.0 && xval > 0.0 {
        (y1 / y0).ln() / (x1 / x0).ln() * value / xval
    } else {
        (y1 - y0) / (x1 - x0)
    };
    (value, slope)
}

/// Estado integrado: [P, m, m_B, y (maré), omega_bar, d omega_bar/dr].
/// P e densidades em M_sun/km^3, massas em M_sun, r em km.
type TovState = [f64; 6];

/// Equações de estrutura. Além da TOV (P, m, m_B):
/// - maré (Hinderer, ApJ 677, 1216 (2008); Postnikov, Prakash & Lattimer,
///   PRD 82, 024016 (2010)): r dy/dr = -y^2 - y F - r^2 Q;
/// - rotação lenta (Hartle, ApJ 150, 1005 (1967)):
///   (1/r^4) d/dr(r^4 j dw/dr) + (4/r)(dj/dr) w = 0, j = e^{-nu/2} sqrt(1-2m/r).
/// Em unidades geometrizadas, m_g = G m, eps_g = G eps, p_g = G p (km^-2).
fn tov_derivatives(
    r: f64,
    s: &TovState,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> TovState {
    let (p, m) = (s[0], s[1]);
    let (eps, deps_dp) = log_linear_interp_with_slope(p_array, eps_array, p);
    let rho = log_linear_interp(p_array, rho_array, p);

    let num = (eps + p) * (m + 4.0 * PI * r.powi(3) * p);
    let den = r * (r - 2.0 * G_C2 * m);
    let metric_term = 1.0 - 2.0 * G_C2 * m / r;

    if den <= 0.0 || metric_term <= 0.0 {
        return [f64::NEG_INFINITY, 0.0, 0.0, 0.0, 0.0, 0.0];
    }

    let dp_dr = -G_C2 * num / den;
    let dm_dr = 4.0 * PI * r.powi(2) * eps;
    let dmb_dr = 4.0 * PI * r.powi(2) * rho / metric_term.sqrt();

    // Geometrizado.
    let (m_g, eps_g, p_g) = (G_C2 * m, G_C2 * eps, G_C2 * p);
    let f = metric_term;
    let mass_term = (m_g + 4.0 * PI * r.powi(3) * p_g) / (r * r * f);

    // Maré: (eps + p)/c_s^2 = (eps + p) deps/dp.
    let big_f = (1.0 - 4.0 * PI * r * r * (eps_g - p_g)) / f;
    let big_q = 4.0 * PI * (5.0 * eps_g + 9.0 * p_g + (eps_g + p_g) * deps_dp) / f
        - 6.0 / (r * r * f)
        - 4.0 * mass_term * mass_term;
    let y = s[3];
    let dy_dr = -(y * y + y * big_f + r * r * big_q) / r;

    // Rotação lenta: j'/j = -nu'/2 + (1/2) d ln f / dr.
    let dnu_dr = 2.0 * mass_term;
    let dlnf_dr = (-2.0 * 4.0 * PI * r * r * eps_g / r + 2.0 * m_g / (r * r)) / f;
    let dlnj_dr = -0.5 * dnu_dr + 0.5 * dlnf_dr;
    let (w, psi) = (s[4], s[5]);
    let dpsi_dr = -(4.0 / r) * psi - dlnj_dr * (psi + 4.0 * w / r);

    [dp_dr, dm_dr, dmb_dr, dy_dr, psi, dpsi_dr]
}

fn axpy(s: &TovState, h: f64, terms: &[(f64, &TovState)]) -> TovState {
    let mut out = *s;
    for i in 0..6 {
        let mut acc = 0.0;
        for (c, k) in terms {
            acc += c * k[i];
        }
        out[i] += h * acc;
    }
    out
}

/// Passo de Cash-Karp. O erro devolvido considera apenas (P, m, m_B), de
/// modo que os passos (e a curva M-R) são os mesmos da TOV pura.
fn rkck_step(
    r: f64,
    y: &TovState,
    h: f64,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> (TovState, [f64; 3]) {
    let d = |r: f64, s: &TovState| tov_derivatives(r, s, p_array, eps_array, rho_array);

    let (a2, a3, a4, a5, a6) = (0.2, 0.3, 0.6, 1.0, 0.875);
    let b21 = 0.2;
    let (b31, b32) = (3.0 / 40.0, 9.0 / 40.0);
    let (b41, b42, b43) = (0.3, -0.9, 1.2);
    let (b51, b52, b53, b54) = (-11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0);
    let (b61, b62, b63, b64, b65) = (
        1631.0 / 55296.0,
        175.0 / 512.0,
        575.0 / 13824.0,
        44275.0 / 110592.0,
        253.0 / 4096.0,
    );
    let (c1, c3, c4, c6) = (37.0 / 378.0, 250.0 / 621.0, 125.0 / 594.0, 512.0 / 1771.0);
    let dc1 = c1 - 2825.0 / 27648.0;
    let dc3 = c3 - 18575.0 / 48384.0;
    let dc4 = c4 - 13525.0 / 55296.0;
    let dc5 = -277.0 / 14336.0;
    let dc6 = c6 - 0.25;

    let k1 = d(r, y);
    let k2 = d(r + a2 * h, &axpy(y, h, &[(b21, &k1)]));
    let k3 = d(r + a3 * h, &axpy(y, h, &[(b31, &k1), (b32, &k2)]));
    let k4 = d(r + a4 * h, &axpy(y, h, &[(b41, &k1), (b42, &k2), (b43, &k3)]));
    let k5 = d(
        r + a5 * h,
        &axpy(y, h, &[(b51, &k1), (b52, &k2), (b53, &k3), (b54, &k4)]),
    );
    let k6 = d(
        r + a6 * h,
        &axpy(y, h, &[(b61, &k1), (b62, &k2), (b63, &k3), (b64, &k4), (b65, &k5)]),
    );

    let yout = axpy(y, h, &[(c1, &k1), (c3, &k3), (c4, &k4), (c6, &k6)]);
    let mut yerr = [0.0; 3];
    for i in 0..3 {
        yerr[i] = h * (dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + dc6 * k6[i]);
    }
    (yout, yerr)
}

fn rkqs_step(
    r: f64,
    y: &TovState,
    htry: f64,
    eps: f64,
    p_array: &[f64],
    eps_array: &[f64],
    rho_array: &[f64],
) -> Option<(TovState, f64, f64)> {
    let safety = 0.9;
    let pgrow = -0.2;
    let pshrink = -0.25;
    let errcon = 1.89e-4;
    let tiny = 1.0e-30;

    let mut h = htry;
    loop {
        let (yout, yerr) = rkck_step(r, y, h, p_array, eps_array, rho_array);

        let mut errmax: f64 = 0.0;
        for i in 0..3 {
            let yscal = y[i].abs() + (h * yerr[i]).abs() + tiny;
            errmax = errmax.max((yerr[i] / yscal).abs());
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

/// Propriedades de uma estrela da sequência.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarProperties {
    /// Massa gravitacional (M_sun).
    pub mass: f64,
    /// Raio (km).
    pub radius: f64,
    /// Massa bariônica (M_sun).
    pub baryonic_mass: f64,
    /// Pressão central (MeV/fm^3).
    pub central_pressure: f64,
    /// Compacidade G M / (R c^2).
    pub compactness: f64,
    /// Redshift gravitacional de superfície, (1 - 2C)^{-1/2} - 1.
    pub redshift: f64,
    /// Número de Love de maré l = 2.
    pub love_k2: f64,
    /// Deformabilidade de maré adimensional, (2/3) k2 C^{-5}.
    pub tidal_deformability: f64,
    /// Momento de inércia (10^45 g cm^2), rotação lenta.
    pub moment_of_inertia: f64,
    /// Momento de inércia adimensional I c^4 / (G^2 M^3).
    pub moment_of_inertia_bar: f64,
}

/// k2 a partir de C e y_R (já corrigido pela descontinuidade de superfície).
/// Abaixo de C = 5e-3 a expressão relativística perde precisão (numerador e
/// denominador são O(C^5)); usa-se o limite newtoniano (2 - y)/(2(3 + y)),
/// com erro O(C) < 1%.
pub fn love_number_k2(c: f64, y: f64) -> f64 {
    if c < 5e-3 {
        return (2.0 - y) / (2.0 * (3.0 + y));
    }
    let num = 8.0 / 5.0 * c.powi(5) * (1.0 - 2.0 * c).powi(2) * (2.0 + 2.0 * c * (y - 1.0) - y);
    let den = 2.0 * c * (6.0 - 3.0 * y + 3.0 * c * (5.0 * y - 8.0))
        + 4.0 * c.powi(3) * (13.0 - 11.0 * y + c * (3.0 * y - 2.0) + 2.0 * c * c * (1.0 + y))
        + 3.0 * (1.0 - 2.0 * c).powi(2) * (2.0 - y + 2.0 * c * (y - 1.0)) * (1.0 - 2.0 * c).ln();
    num / den
}

/// 1 M_sun km^2 em unidades de 10^45 g cm^2.
const MSUN_KM2_IN_1E45_G_CM2: f64 = 1.98847e33 * 1e10 / 1e45;

pub fn integrate_star_properties(
    pc_tov: f64,
    p_min: f64,
    p_tov: &[f64],
    eps_tov: &[f64],
    rho_tov: &[f64],
) -> Option<StarProperties> {
    if p_tov.len() < 5 || eps_tov.len() < 5 || rho_tov.len() < 5 {
        return None;
    }
    if p_tov.len() != eps_tov.len() || p_tov.len() != rho_tov.len() {
        return None;
    }
    if !pc_tov.is_finite() || !p_min.is_finite() || pc_tov <= p_min {
        return None;
    }
    let mut r = 1e-5;
    // Condições regulares no centro: y = 2, omega_bar = 1, omega_bar' = 0.
    let mut y: TovState = [pc_tov, 0.0, 0.0, 2.0, 1.0, 0.0];
    let mut h = 1.0e-2;
    let r_end = 30000.0;
    // Tolerância relativa do passo adaptativo. Com a superfície localizada
    // por bisseção abaixo, 1e-10 muda M e R de EoS realistas em < 1e-7
    // relativo e é ~8x mais rápido que 3e-16 (ver tests/verification.rs).
    let eps = 1.0e-10;
    let mut steps = 0u64;
    let max_steps = 90000u64;

    while r < r_end && steps < max_steps {
        if !y.iter().all(|value| value.is_finite()) || y[0] <= p_min || r <= 2.0 * G_C2 * y[1] {
            return None;
        }

        if r + h > r_end {
            h = r_end - r;
        }

        let (ynew, hdid, hnext) = rkqs_step(r, &y, h, eps, p_tov, eps_tov, rho_tov)?;

        if !ynew.iter().all(|value| value.is_finite())
            || !hdid.is_finite()
            || hdid <= 0.0
            || !hnext.is_finite()
        {
            return None;
        }

        if ynew[0] <= p_min {
            // Locate the surface inside the accepted step instead of
            // returning the overshot state.  A TOV point is successful only
            // when this pressure event is actually reached.  The step size
            // to the event P = p_min is found by bisection, re-integrating
            // from (r, y) with the same Cash-Karp step, so the surface is as
            // accurate as the integration itself.
            let pressure_drop = y[0] - ynew[0];
            if !pressure_drop.is_finite() || pressure_drop <= 0.0 {
                return None;
            }
            let (mut h_lo, mut h_hi) = (0.0, hdid);
            let mut y_surface = y;
            for _ in 0..60 {
                let h_mid = 0.5 * (h_lo + h_hi);
                let (y_mid, _) = rkck_step(r, &y, h_mid, p_tov, eps_tov, rho_tov);
                if !y_mid.iter().all(|value| value.is_finite()) {
                    return None;
                }
                if y_mid[0] > p_min {
                    h_lo = h_mid;
                    y_surface = y_mid;
                } else {
                    h_hi = h_mid;
                }
            }
            let radius = r + h_lo;
            let (mass, baryonic_mass) = (y_surface[1], y_surface[2]);
            if !(radius.is_finite() && mass.is_finite() && baryonic_mass.is_finite()) {
                return None;
            }

            let compactness = G_C2 * mass / radius;
            // Descontinuidade de densidade na superfície (Damour & Nagar 2009;
            // Postnikov et al. 2010): y_R -> y_R - 3 eps_s / <eps>.
            let eps_surface = log_linear_interp(p_tov, eps_tov, p_min);
            let mean_density = 3.0 * mass / (4.0 * PI * radius.powi(3));
            let y_r = y_surface[3] - 3.0 * eps_surface / mean_density;
            let love_k2 = love_number_k2(compactness, y_r);

            // Exterior: omega_bar = Omega - 2J/r^3, logo J = R^4 omega_bar'/6.
            let j = radius.powi(4) * y_surface[5] / 6.0;
            let omega = y_surface[4] + 2.0 * j / radius.powi(3);
            let inertia_geo = j / omega; // km^3

            return Some(StarProperties {
                mass,
                radius,
                baryonic_mass,
                central_pressure: pc_tov / MEV_FM3_TO_MSUN_KM3,
                compactness,
                redshift: 1.0 / (1.0 - 2.0 * compactness).sqrt() - 1.0,
                love_k2,
                tidal_deformability: 2.0 / 3.0 * love_k2 / compactness.powi(5),
                moment_of_inertia: inertia_geo / G_C2 * MSUN_KM2_IN_1E45_G_CM2,
                moment_of_inertia_bar: inertia_geo / (G_C2 * mass).powi(3),
            });
        }

        y = ynew;
        r += hdid;

        h = hnext;
        steps += 1;
    }

    // Reaching the radial/step budget is not reaching the stellar surface.
    None
}

pub fn integrate_star(
    pc_tov: f64,
    p_min: f64,
    p_tov: &[f64],
    eps_tov: &[f64],
    rho_tov: &[f64],
) -> Option<(f64, f64, f64, f64)> {
    integrate_star_properties(pc_tov, p_min, p_tov, eps_tov, rho_tov)
        .map(|s| (s.mass, s.radius, s.baryonic_mass, s.central_pressure))
}

/// Unifica a crosta personalizada (1/fm⁴) com a EoS do núcleo, descartando dados inválidos
pub fn unify_with_crust(
    core_eps: &[f64],
    core_p: &[f64],
    core_rho: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // Constante de conversão de 1/fm⁴ para MeV/fm³
    const HBARC: f64 = 197.3269804;

    // Dados da crosta em 1/fm⁴- Baym-Pethick-Sutherland - from https://github.com/mrpelicer/nuclear_physics
    const CRUST_P_FM4: &[f64] = &[
        1.212e-11, 8.236e-11, 2.764e-10, 5.152e-10, 1.593e-09, 4.023e-09, 1.380e-08, 3.315e-08,
        1.077e-07, 2.559e-07, 3.479e-07, 4.729e-07, 6.430e-07, 9.147e-07, 1.041e-06, 1.840e-06,
        2.469e-06, 2.496e-06, 2.642e-06, 2.878e-06, 3.110e-06, 3.425e-06, 3.852e-06, 4.425e-06,
        5.181e-06, 6.168e-06, 8.198e-06, 1.109e-05, 1.509e-05, 2.050e-05, 2.767e-05, 3.701e-05,
        5.361e-05,
    ];

    const CRUST_E_FM4: &[f64] = &[
        9.387e-08, 3.738e-07, 9.392e-07, 1.489e-06, 3.741e-06, 7.465e-06, 1.877e-05, 3.747e-05,
        9.418e-05, 1.881e-04, 2.369e-04, 2.982e-04, 3.758e-04, 5.242e-04, 5.958e-04, 9.452e-04,
        1.222e-03, 1.268e-03, 1.486e-03, 1.879e-03, 2.264e-03, 2.765e-03, 3.400e-03, 4.182e-03,
        5.131e-03, 6.260e-03, 8.329e-03, 1.090e-02, 1.402e-02, 1.776e-02, 2.218e-02, 2.732e-02,
        3.542e-02,
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

            // A coluna EOS é n_B/n_0. Converta primeiro para fm^-3 e
            // depois para uma densidade de massa-energia em MeV/fm^3.
            raw_rho.push(core_rho[i] * N0 * M_NUCLEON);
        }
    }

    // 3. FILTRO DE MONOTONIA ESTRITA (Garante compatibilidade com a GSL)
    let mut combined: Vec<(f64, f64, f64)> = raw_p
        .into_iter()
        .zip(raw_eps.into_iter())
        .zip(raw_rho.into_iter())
        .map(|((p, eps), rho)| (p, eps, rho))
        .collect();
    // Ordena pela densidade de energia, que cresce com a densidade em qualquer
    // EoS; pontos em que P não cresce (p.ex. P_perp com magnetização na
    // matéria diluída) são descartados pelo filtro abaixo. Para uma EoS
    // monótona o resultado é idêntico ao da ordenação por P.
    combined.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

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

/// Sequência de estrelas com M, R, M_B, maré, momento de inércia e redshift.
pub fn generate_star_sequence(
    eps_array: &[f64],
    p_array: &[f64],
    rho_array: &[f64],
    with_crust: bool,
) -> Vec<StarProperties> {
    let mut stars = Vec::new();

    // 1. Costura a crosta APENAS se a flag for verdadeira
    let (clean_eps, clean_p, clean_rho) = if with_crust {
        // A função unify_with_crust já faz a conversão do núcleo internamente agora
        unify_with_crust(eps_array, p_array, rho_array)
    } else {
        // ``rho_array`` segue o contrato da coluna 0 da EOS: n_B/n_0.
        let converted_rho: Vec<f64> = rho_array
            .iter()
            .map(|&nb_over_n0| nb_over_n0 * N0 * M_NUCLEON)
            .collect();

        (eps_array.to_vec(), p_array.to_vec(), converted_rho)
    };

    // 2. Limpa e ordena para satisfazer a GSL
    let (clean_eps, clean_p, clean_rho) = clean_eos_with_rho(&clean_eps, &clean_p, &clean_rho);

    if clean_p.len() < 5 {
        return stars;
    }

    // 3. Converte unidades uma unica vez
    let eps_tov: Vec<f64> = clean_eps.iter().map(|&e| e * MEV_FM3_TO_MSUN_KM3).collect();
    let p_tov: Vec<f64> = clean_p.iter().map(|&p| p * MEV_FM3_TO_MSUN_KM3).collect();
    let rho_tov: Vec<f64> = clean_rho
        .iter()
        .map(|&rho| rho * MEV_FM3_TO_MSUN_KM3)
        .collect();
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
        if let Some(star) = integrate_star_properties(pc_tov, p_min, &p_tov, &eps_tov, &rho_tov) {
            if star.mass > 0.05 && star.radius > 2.0 {
                stars.push(star);
            }
        }
    }

    stars
}

/// Perfil radial da estrela de pressão central `pc_mev` (MeV/fm^3), com a
/// mesma EoS (crosta BPS incluída) e o mesmo integrador de
/// `generate_star_sequence`: pontos [r (km), P (MeV/fm^3), m (M_sun)] em cada
/// passo aceito, do centro até a superfície (P = P_min, interpolada no último
/// passo). Vazio se a integração falhar.
pub fn radial_profile(eps_array: &[f64], p_array: &[f64], rho_array: &[f64], pc_mev: f64) -> Vec<[f64; 3]> {
    let (eps, p, rho) = unify_with_crust(eps_array, p_array, rho_array);
    let (eps, p, rho) = clean_eos_with_rho(&eps, &p, &rho);
    if p.len() < 5 || !pc_mev.is_finite() {
        return Vec::new();
    }
    let eps_tov: Vec<f64> = eps.iter().map(|&e| e * MEV_FM3_TO_MSUN_KM3).collect();
    let p_tov: Vec<f64> = p.iter().map(|&v| v * MEV_FM3_TO_MSUN_KM3).collect();
    let rho_tov: Vec<f64> = rho.iter().map(|&v| v * MEV_FM3_TO_MSUN_KM3).collect();
    let p_min = p_tov[0];
    let pc_tov = pc_mev * MEV_FM3_TO_MSUN_KM3;
    if pc_tov <= p_min {
        return Vec::new();
    }

    let (mut r, mut h) = (1e-5, 1.0e-2);
    let mut y: TovState = [pc_tov, 0.0, 0.0, 2.0, 1.0, 0.0];
    let mut points = vec![[0.0, pc_mev, 0.0]];
    for _ in 0..90000 {
        let Some((ynew, hdid, hnext)) = rkqs_step(r, &y, h, 1.0e-10, &p_tov, &eps_tov, &rho_tov) else {
            return Vec::new();
        };
        if !ynew.iter().all(|v| v.is_finite()) || hdid <= 0.0 {
            return Vec::new();
        }
        if ynew[0] <= p_min {
            let t = (y[0] - p_min) / (y[0] - ynew[0]);
            let m = y[1] + t * (ynew[1] - y[1]);
            points.push([r + t * hdid, p_min / MEV_FM3_TO_MSUN_KM3, m]);
            return points;
        }
        y = ynew;
        r += hdid;
        h = hnext;
        points.push([r, y[0] / MEV_FM3_TO_MSUN_KM3, y[1]]);
    }
    Vec::new()
}

pub fn generate_mr_curve(
    eps_array: &[f64],
    p_array: &[f64],
    rho_array: &[f64],
    with_crust: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let stars = generate_star_sequence(eps_array, p_array, rho_array, with_crust);
    (
        stars.iter().map(|s| s.mass).collect(),
        stars.iter().map(|s| s.radius).collect(),
        stars.iter().map(|s| s.baryonic_mass).collect(),
        stars.iter().map(|s| s.central_pressure).collect(),
    )
}

/// Procura e interpola os dados de quarks para um mun específico
fn find_quark_point(
    quark_eos: &[[f64; RESULTS_SIZE]],
    target_mun: f64,
) -> Option<[f64; RESULTS_SIZE]> {
    // Busca binária para encontrar a posição do mun (índice 17)
    let pos = quark_eos.binary_search_by(|row| {
        row[17]
            .partial_cmp(&target_mun)
            .unwrap_or(std::cmp::Ordering::Equal)
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

pub fn unify_hybrid_eos(
    hadron_eos: &[[f64; RESULTS_SIZE]],
    quark_eos: &[[f64; RESULTS_SIZE]],
) -> (Vec<f64>, Vec<f64>) {
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

    // Ordena pela densidade de energia, que cresce com a densidade em qualquer
    // EoS; pontos em que P não cresce (p.ex. P_perp com magnetização na
    // matéria diluída) são descartados pelo filtro abaixo. Para uma EoS
    // monótona o resultado é idêntico ao da ordenação por P.
    combined.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut safe_p = Vec::with_capacity(combined.len());
    let mut safe_eps = Vec::with_capacity(combined.len());
    let mut safe_rho = Vec::with_capacity(combined.len());

    let mut last_p = -f64::INFINITY;

    for (pres, energy, rho_val) in combined {
        // Repeated vacuum points are continuation aids for the microscopic
        // solver, not material through which a star should be integrated.
        if pres >= 0.0
            && energy > 0.0
            && rho_val > 0.0
            && pres > last_p + 1e-18
            && pres.is_finite()
            && energy.is_finite()
            && rho_val.is_finite()
        {
            safe_p.push(pres);
            safe_eps.push(energy);
            safe_rho.push(rho_val);
            last_p = pres;
        }
    }

    (safe_eps, safe_p, safe_rho)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vacuum_prefix_does_not_change_the_mr_curve_or_baryonic_mass_scale() {
        let eps: Vec<f64> = (0..101).map(|i| 10.0 + i as f64 * 10.0).collect();
        let pressure: Vec<f64> = eps.iter().map(|energy| 0.25 * (energy - 10.0)).collect();
        let density: Vec<f64> = eps.iter().map(|energy| energy / (M_NUCLEON * N0)).collect();

        let reference = generate_mr_curve(&eps, &pressure, &density, false);
        assert!(!reference.0.is_empty());

        let mut prefixed_eps = vec![0.0; 20];
        let mut prefixed_pressure = vec![0.0; 20];
        let mut prefixed_density = vec![0.0; 20];
        prefixed_eps.extend_from_slice(&eps);
        prefixed_pressure.extend_from_slice(&pressure);
        prefixed_density.extend_from_slice(&density);
        let prefixed =
            generate_mr_curve(&prefixed_eps, &prefixed_pressure, &prefixed_density, false);

        assert_eq!(prefixed, reference);
        let max_index = reference
            .0
            .iter()
            .enumerate()
            .max_by(|(_, lhs), (_, rhs)| lhs.total_cmp(rhs))
            .map(|(index, _)| index)
            .unwrap();
        let mass = reference.0[max_index];
        let baryonic_mass = reference.2[max_index];
        assert!(baryonic_mass > mass);
        assert!(baryonic_mass < 1.5 * mass);
    }
}
