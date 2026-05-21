// src/core/darkphotons.rs

use crate::constants::PI2;
use crate::core::constants::{
    M_NUCLEON, MB, ML, N0, QE, RESULTS_SIZE
};
use crate::core::model::ModelParams;
use nalgebra::{Matrix4, Vector4};

#[derive(Clone)]
pub struct DarkPhotonsMatter {
    // Parâmetros fixos
    pub model: ModelParams,
    pub m_nuc: f64,
    pub qe: f64,
    pub ml: [f64; 2],
    pub mb: [f64; 8],
    pub m_eff: [f64; 8],
    pub mu_b: [f64; 8],
    pub charges_b: [f64; 8],
    pub xs: f64,

    // Limites do loop (podem ser ajustados)
    pub mun_inf: f64,
    pub mun_sup: f64,
    pub n_points: usize,

    // --- Estado mutável para o ponto atual ---
    // Potenciais químicos
    pub mun: f64,
    pub mue: f64,
    pub mup: f64,

    // Densidades
    pub nb: [f64; 8],
    pub nl: [f64; 2],
    pub nbt: f64,               // densidade bariônica total

    // Densidades escalares
    pub rhosb: f64,
    pub rhos_b: [f64; 8],
    pub rhos_l: [f64; 2],

    // Energias de Fermi e momentos (para EOS)
    pub ef_b: [f64; 8],
    pub ef_l: [f64; 2],      // Energias de Fermi: [0]=e, [1]=mu

    // Acoplamentos (xv_v para omega, xv_r para rho)
    pub xv_v: [f64; 8], // g_wB / g_wN
    pub xv_r: [f64; 8], // g_rB / g_rN

    // Momentos de Fermi por nível de Landau (Vetorizados)
    pub kf_b_up: [Vec<f64>; 8],   // [Barião][Nível nu]
    pub kf_b_down: [Vec<f64>; 8],
    pub f_l: [Vec<f64>; 2],  // Momentos de Fermi: [0]=fe, [1]=fmu

    pub isospin_factor: [f64; 8],
}

impl DarkPhotonsMatter {
    // constante estática para acoplamento sigma
    const X_SIGMA: [f64; 8] = [1.0, 1.0, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];

    pub fn new(model: ModelParams) -> Self {
        let m_nuc = M_NUCLEON;
        let qe = QE;
        let ml = ML;
        let mb = MB;
        let xs = 0.7;

        let m_eff = [0.0; 8];
        let mu_b = [0.0; 8];
        let charges_b = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0];

        let xv_v = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];
        let xv_r = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];

        let kf_b_up = std::array::from_fn(|_| Vec::new());
        let kf_b_down = std::array::from_fn(|_| Vec::new());
        let f_l = std::array::from_fn(|_| Vec::new());

        let ef_l = [0.0; 2];

        let isospin_factor = [-0.5, 0.5, 0.0, -1.0, 0.0, 1.0, -0.5, 0.5];

        DarkPhotonsMatter {
            model,
            m_nuc,
            qe,
            ml,
            mb,
            m_eff,
            mu_b,
            charges_b,
            xs,
            mun_inf: 0.02,
            mun_sup: 1.80,
            n_points: 1201,

            mun: 0.0,
            mue: 0.0,
            mup: 0.0,

            nb: [0.0; 8],
            nl: [0.0; 2],
            nbt: 0.0,
            
            rhosb: 0.0,
            rhos_b: [0.0; 8],
            rhos_l: [0.0; 2],

            ef_b: [0.0; 8],
            ef_l,

            xv_v,
            xv_r,

            kf_b_up,
            kf_b_down,
            f_l,

            isospin_factor: isospin_factor,
        }
    }

     // Métodos builder
    pub fn with_limits(mut self, inf: f64, sup: f64) -> Self {
        self.mun_inf = inf;
        self.mun_sup = sup;
        self
    }

    pub fn with_points(mut self, n: usize) -> Self {
        self.n_points = n;
        self
    }

    // Mapeamento das variáveis (vindo do solver)
    pub fn mapping(&self, x: &[f64]) -> (f64, f64, f64, f64) {
        let mue = x[0];
        let vsigma = x[1]; // Removido o .sin().powi(2) que destruía o Jacobiano
        let vomega = x[2];
        let vrho = x[3];
        (mue, vsigma, vomega, vrho)
    }

    // Função de resíduo (chamada pelo solver numérico)
    pub fn funcv(&mut self, x: &[f64]) -> [f64; 4] {
        let (mue, vsigma, vomega, vrho) = self.mapping(x);

        self.mue = mue;
        self.mup = self.mun - mue;
        self.mu_b[0] = self.mun;

        // massas efetivas
        for i in 0..8 {
            self.m_eff[i] = self.mb[i] - Self::X_SIGMA[i] * vsigma;
        }

        // potenciais químicos de todas as outras partículas
        for i in 1..8 {
            self.mu_b[i] = self.mu_b[0] - self.charges_b[i] * mue;
        }

        // calcular densidades
        calculate_all_densities(self, vomega, vrho);

        let fsigma = self.equation_sigma(vsigma);
        let fomega = self.equation_omega(vomega);
        let frho = self.equation_rho(vrho);
        let charge_neutral = self.charge_neutrality();

        [fsigma, fomega, frho, charge_neutral]
    }

    fn equation_sigma(&self, vsigma: f64) -> f64 {
        let gs2 = self.model.gs.powi(2);
        gs2 * (self.rhosb - self.model.rb * vsigma.powi(2) - self.model.rc * vsigma.powi(3)) - vsigma
    }

    // Equações de Campo Vetorizadas para suportar as partículas com total precisão
    fn equation_omega(&self, vomega: f64) -> f64 {
        let mut sum_baryon = 0.0;
        for i in 0..8 {
            sum_baryon += self.nb[i] * self.xv_v[i];
        }
        self.model.gv.powi(2) * sum_baryon - vomega
    }

    fn equation_rho(&self, vrho: f64) -> f64 {
        let mut sum_source = 0.0;
        for i in 0..8 {
            // A fonte para o rho é baseada no negativo do isospin
            sum_source += self.isospin_factor[i] * self.nb[i] * self.xv_r[i];
        }
        self.model.gr.powi(2) * sum_source - vrho
    }
    
    fn charge_neutrality(&self) -> f64 {
        let charge_baryons: f64 = self.nb
            .iter()
            .zip(self.charges_b.iter())
            .map(|(n, q)| n * q)
            .sum();

        // Leptons: e⁻ and μ⁻ have charge -1
        let charge_leptons: f64 = self.nl.iter().map(|n| -n).sum();

        charge_baryons + charge_leptons
    }

    // Resolve para um dado mun e chute inicial, retorna solução e resultado
    pub fn solve_point(&mut self, mun: f64, initial_x: &[f64]) -> Option<([f64; 4], [f64; RESULTS_SIZE])> {
        self.mun = mun;

        let mut x = Vector4::from_column_slice(initial_x);
        let tolerance = 1e-13;
        let max_iterations = 100;
        let mut converged = false;

        for _ in 0..max_iterations {
            let f_val_arr = self.funcv(x.as_slice());
            let f_val = Vector4::from_column_slice(&f_val_arr);
            let f_norm = f_val.norm();

            if f_norm < tolerance {
                converged = true;
                break;
            }

            let mut j_matrix = Matrix4::zeros();

            for i in 0..4 {
                let h = 1e-7 * (x[i].abs() + 1e-2);
                let mut x_temp = x;
                x_temp[i] += h;

                let f_temp_arr = self.funcv(x_temp.as_slice());
                let f_temp = Vector4::from_column_slice(&f_temp_arr);

                let column_derivative = (f_temp - f_val) / h;
                j_matrix.set_column(i, &column_derivative);
            }

            let delta_x = match j_matrix.lu().solve(&(-f_val)) {
                Some(step) => step,
                None => break,
            };

            let mut alpha = 1.0;
            let mut step_accepted = false;

            for _ in 0..15 {
                let x_try = x + alpha * delta_x;
                let f_new_arr = self.funcv(x_try.as_slice());
                let f_new = Vector4::from_column_slice(&f_new_arr);
                let f_new_norm = f_new.norm();

                if f_new_norm.is_nan() {
                    alpha *= 0.5;
                    continue;
                }

                if f_new_norm < f_norm {
                    x = x_try;
                    step_accepted = true;
                    break;
                }

                alpha *= 0.5;
            }

            if !step_accepted {
                x += 0.001 * delta_x;
                let _ = self.funcv(x.as_slice()); // mantém estado interno consistente com x
            }
        }

        if !converged {
            return None;
        }

        // garante estado físico final consistente
        let _ = self.funcv(x.as_slice());

        let x_final = [x[0], x[1], x[2], x[3]];
        let (mue, vsigma, vomega, vrho) = self.mapping(&x_final);
        let (ener, press) = compute_eos(self, mue, vsigma, vomega, vrho);

        let nb_total = self.nb.iter().sum::<f64>();
        let nbtd = nb_total * (self.m_nuc / 197.32).powi(3);

        let factor_mev_fm3 = self.m_nuc * (self.m_nuc / 197.32).powi(3);
        let ener_conv = ener * factor_mev_fm3;
        let press_conv = press * factor_mev_fm3;

        let ener_final = ener_conv;
        let press_final = press_conv;

        if ener_final >= 0.0 && press_final >= 0.0 {
            let result = [
                nbtd / N0, //  0: n/n0 (adimensional)
                ener_final,   //  1: Energia Total [MeV/fm^3]
                press_final,  //  2: Pressão Total [MeV/fm^3]
                self.nl[0],   //  3: e- [fm^-3]
                self.nl[1],   //  4: mu- [fm^-3]
                self.nb[0],   //  5: n [fm^-3]
                self.nb[1],   //  6: p [fm^-3]
                self.nb[2],   //  7: L0 [fm^-3]
                self.nb[3],   //  8: S- [fm^-3]
                self.nb[4],   //  9: S0 [fm^-3]
                self.nb[5],   // 10: S+ [fm^-3]
                self.nb[6],   // 11: X- [fm^-3]
                self.nb[7],   // 12: X0 [fm^-3]
                vsigma,       // 13: Campo Sigma [MeV]
                vomega,       // 14: Campo Omega [MeV]
                vrho,         // 15: Campo Rho [MeV]
                self.m_eff[0] / self.m_nuc, // 16: m*/mN (adimensional)
                self.mun,     // 17: mu_n [adimensional no código atual]
                mue,          // 18: mu_e [adimensional no código atual]
                0.0,
                0.0,
            ];
            Some((x_final, result))
        } else {
            None
        }
    }

}


pub fn calculate_all_densities(engine: &mut DarkPhotonsMatter, vomega: f64, vrho: f64) {
    for i in 0..8 {
        let (rs, rb) = density_baryon(engine, i, vomega, vrho);
        engine.rhos_b[i] = rs;
        engine.nb[i] = rb;
    }

    for i in 0..2 {
        let (rl, nl) = density_lepton(engine, i);
        engine.rhos_l[i] = rl;
        engine.nl[i] = nl;
    }

    let rhos_n_p = engine.rhos_b[0] + engine.rhos_b[1];
    let rhos_h = engine.rhos_b[2..8].iter().sum::<f64>();
    
    engine.rhosb = rhos_n_p + engine.xs * rhos_h; 
    engine.nbt = engine.nb.iter().sum();
}

fn density_baryon(engine: &mut DarkPhotonsMatter, idx: usize, vomega: f64, vrho: f64) -> (f64, f64) {
    let m = engine.m_eff[idx];
    
    let ef = engine.mu_b[idx] 
             - (engine.xv_v[idx] * vomega) 
             - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx]); 
    
    engine.ef_b[idx] = ef;
    
    if ef <= 0.0 { return (0.0, 0.0); }

    // TRATAMENTO PARA O CASO ISOTRÓPICO (B=0)
    let kf2 = ef.powi(2) - m.powi(2);
    if kf2 > 0.0 {
        let kf = kf2.sqrt();
        let dens = kf.powi(3) / (3.0 * PI2);
        let m_safe = m.abs().max(1e-15);
        let rhos = (m / (2.0 * PI2)) * (ef * kf - m.powi(2) * ((kf + ef) / m_safe).ln());
        
        return (rhos, dens);
    }
    return (0.0, 0.0);

}

pub fn density_lepton(engine: &mut DarkPhotonsMatter, idx: usize) -> (f64, f64) {
    let mue = engine.mue;      
       
    if mue <= 0.0 { return (0.0, 0.0); }

    let m = engine.ml[idx];

    let kf2 = mue.powi(2) - m.powi(2);
    if kf2 > 0.0 {
        let kf = kf2.sqrt();
        let dens_val = kf.powi(3) / (3.0 * PI2);
        let m_safe = m.abs().max(1e-15);
        let rhos_val = (m / (2.0 * PI2)) * (mue * kf - m.powi(2) * ((kf + mue) / m_safe).ln());
        
        return (rhos_val, dens_val);
    }
    return (0.0, 0.0)
}

pub fn compute_eos(engine: &DarkPhotonsMatter, mue: f64, vsigma: f64, vomega: f64, vrho: f64) -> (f64, f64) {
    // 1. Energia dos mésons (Potenciais de campo)
    // Inclui termos de massa e auto-interações (rb, rc para sigma e rxi para omega)
    let enerf = (vsigma / engine.model.gs).powi(2) / 2.0
        + (vomega / engine.model.gv).powi(2) / 2.0
        + (vrho / engine.model.gr).powi(2) / 2.0
        + engine.model.rb * vsigma.powi(3) / 3.0
        + engine.model.rc * vsigma.powi(4) / 4.0
        + engine.model.rxi * vomega.powi(4) / 4.0;

    let mut enerbar = 0.0;

    // --- LOOP DE BARIÕES (0:n, 1:p, 2:L0, 3:S-, 4:S0, 5:S+, 6:X-, 7:X0) ---
    for i in 0..8 {
        let ef = engine.ef_b[i];
        if ef <= 0.0 { continue; }

        // --- Partículas Neutras (n, L0, S0, X0) ---
        // O AMM desdobra a partícula em 2 estados de spin (Up e Down)
        for &kf in [engine.kf_b_up[i].first().copied(), engine.kf_b_down[i].first().copied()]
            .iter()
            .flatten()
        {
            if kf > 0.0 {
                let m_spin = (ef.powi(2) - kf.powi(2)).max(0.0).sqrt();
                let m_safe = m_spin.max(1e-15);
                
                // Fórmula para um único estado de spin (g=1), fator 1/4pi^2
                enerbar += (1.0 / (4.0 * PI2)) * (
                    ef.powi(3) * kf / 2.0
                    - (m_spin / 4.0) * (m_spin * kf * ef + m_spin.powi(3) * ((kf + ef) / m_safe.abs()).ln())
                );
            }
        }
    }

        // --- LOOP DE LÉPTONS (0:e-, 1:mu-) ---
    let mut enerlep = 0.0;
    for i in 0..2 {
        let ef = engine.ef_l[i];
        if ef <= 0.0 { continue; }

        let kf = engine.f_l[i].first().copied().unwrap_or(0.0);
        if kf > 0.0 {
            let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
            enerlep += 2.0 * (1.0 / (4.0 * PI2)) * (
                ef.powi(3) * kf / 2.0
                - (m_spin / 4.0) * (m_spin * kf * ef + m_spin.powi(3) * ((kf + ef) / m_spin.abs()).ln())
            );
        }
    }

    // Energia Total (Mésons + Bariões + Léptons)
    let ener = enerf + enerbar + enerlep;

    // Pressão via relação termodinâmica: P = sum(mu_i * n_i) - epsilon
    let mut press_sum = 0.0;
    for i in 0..8 {
        press_sum += engine.mu_b[i] * engine.nb[i];
    }
    for i in 0..2 {
        press_sum += mue * engine.nl[i]; // mu_e = mu_mu = mue
    }
    
    let press = press_sum - ener;

    (ener, press)
}