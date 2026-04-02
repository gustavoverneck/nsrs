// solver/physics.rs
#![allow(unused)]

use crate::core::constants::{
    AMML0, AMMN, AMMP, AMMS0, AMMSM, AMMSP, AMMX0, AMMXM, BCE, BCE_G, BDD_ALPHAA, BDD_BETAA, M_NUCLEON, MB, ML, N0, QE, RESULTS_SIZE
};
use crate::core::model::ModelParams;
use nalgebra::{Matrix4, Vector4};
use rgsl::exponential::exp;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NlemModel {
    Maxwell,        // Eletromagnetismo Clássico (Linear)
    Modmax(f64),    // Eletrodinâmica ModMax (recebe parâmetro csi)
    Log(f64),       // Eletrodinâmica Logarítmica (recebe parâmetro csi)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MagneticTopology {
    Isotropic,   // Campo caótico (P_mag = 1/3 eps_mag) -> Obrigatório para TOV 1D
    Anisotropic, // Campo alinhado (P_mag = eps_mag) -> Para solvers 2D futuros
}

impl NlemModel {
    /// Recebe o campo magnético original (bg) e retorna o campo magnético 
    /// EFETIVO alterado pelo modelo não-linear.
    pub fn effective_bg(&self, bg: f64) -> f64 {
        match self {
            NlemModel::Maxwell => bg,
            
            NlemModel::Modmax(csi) => {
                // Fórmula: bg * exp(-csi)
                bg * exp(-csi)
            },
            
            NlemModel::Log(csi) => {
                // Modelo de Soleng Padrão: bg / (1 + bg^2 / 2*csi^2)
                // O sinal de '+' garante que nunca haverá divisão por zero!
                let denom = 1.0 + bg.powi(2) / (2.0 * csi.powi(2));
                bg / denom
            }
        }
    }

    /// Calcula a Densidade de Energia Magnética Macroscópica da NLEM
    /// Essa é a energia real que curva o espaço-tempo na equação TOV.
    pub fn magnetic_energy(&self, bg: f64, ebsi_maxwell: f64) -> f64 {
        match self {
            NlemModel::Maxwell => ebsi_maxwell,
            
            NlemModel::Modmax(csi) => {
                // ε_MM = ε_Max * e^(-csi)
                ebsi_maxwell * exp(-csi)
            },
            
            NlemModel::Log(csi) => {
                let bg_sq = bg.powi(2);
                let csi_sq = csi.powi(2);
                
                // x = B^2 / (2 * csi^2)
                let x = bg_sq / (2.0 * csi_sq);
                
                // Limite de segurança matemático para x muito próximo de zero
                if x < 1e-15 {
                    ebsi_maxwell
                } else {
                    // Usamos ln_1p(x) que calcula ln(1 + x)
                    ebsi_maxwell * (x.ln_1p() / x)
                }
            }
        }
    }
}

#[derive(Clone)]
pub struct HadronsMatter {
    // Parâmetros fixos
    pub model: ModelParams,
    pub nlem: NlemModel,
    pub topology: MagneticTopology,
    pub bg: f64,
    pub b: f64,
    pub m_nuc: f64,
    pub qe: f64,
    pub ml: [f64; 2],
    pub mb: [f64; 8],
    pub m_eff: [f64; 8],
    pub mu_b: [f64; 8],
    pub charges_b: [f64; 8],
    pub amm_b: [f64; 8],
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

    // Contadores de níveis por spin
    pub n_b_up: [usize; 8],
    pub n_b_down: [usize; 8],
    pub n_l: [usize; 2],     // Contadores: [0]=ne, [1]=nu

    pub max_landau_limit: usize,

    pub isospin_factor: [f64; 8],
}

impl HadronsMatter {
    pub fn new(model: ModelParams, bg: f64) -> Self {
        let m_nuc = M_NUCLEON;
        let qe = QE;
        let ml = ML;
        let mb = MB;
        let xs = 0.7;
        let b0 = bg / BCE_G;
        let b = b0 * BCE;

        let m_eff = [0.0; 8];
        let mu_b = [0.0; 8];
        let charges_b = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0];

        let amm_b = [AMMN, AMMP, AMML0, AMMSM, AMMS0, AMMSP, AMMXM, AMMX0];

        let xv_v = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];
        let xv_r = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];
        let max_landau_limit = 100_000;

        let kf_b_up = std::array::from_fn(|_| vec![0.0; max_landau_limit]);
        let kf_b_down = std::array::from_fn(|_| vec![0.0; max_landau_limit]);
        let f_l = std::array::from_fn(|_| vec![0.0; max_landau_limit]);

        let ef_l = [0.0; 2];
        let n_l = [0; 2];

        let isospin_factor = [-0.5, 0.5, 0.0, -1.0, 0.0, 1.0, -0.5, 0.5];

        HadronsMatter {
            model,
            nlem: NlemModel::Maxwell,
            topology: MagneticTopology::Isotropic,
            bg,
            b,
            m_nuc,
            qe,
            ml,
            mb,
            m_eff,
            mu_b,
            charges_b,
            amm_b,
            xs,
            mun_inf: 0.92,
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

            n_b_up: [0; 8],
            n_b_down: [0; 8],
            n_l,

            max_landau_limit: max_landau_limit,

            isospin_factor: isospin_factor,
        }
    }
    /// Define a topologia das linhas de campo magnético
    pub fn with_topology(mut self, top: MagneticTopology) -> Self {
        self.topology = top;
        self
    }

    /// Builder para acoplar o Eletromagnetismo Não-Linear
    pub fn with_nlem(mut self, nlem: NlemModel) -> Self {
        self.nlem = nlem;
        
        // 1. Calcula o campo macroscópico efetivo usando o Enum
        let bg_effective = self.nlem.effective_bg(self.bg);
        
        // 2. Recalcula o 'b' que vai para os Níveis de Landau usando o novo bg
        let b0 = bg_effective / BCE_G;
        self.b = b0 * BCE;
        
        self
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
    pub fn funcv(&mut self, x: &[f64]) -> Vec<f64> {
        let (mue, vsigma, vomega, vrho) = self.mapping(x);

        let x_sigma = [1.0, 1.0, self.xs, self.xs, self.xs, self.xs, self.xs, self.xs];

        self.mue = mue;
        self.mup = self.mun - mue;

        self.mu_b[0] = self.mun;

        // 1) Massas efetivas
        for i in 0..8 {
            self.m_eff[i] = self.mb[i] - x_sigma[i] * vsigma;
        }

        // 1b) Potenciais químicos dos outros bárions
        for i in 1..8 {
            self.mu_b[i] = self.mu_b[0] - self.charges_b[i] * mue;
        }

        // 2) Atualização dinâmica do campo magnético local
        let nb_total: f64 = self.nb.iter().sum();
        let nbtd = nb_total * (self.m_nuc / 197.32).powi(3); // fm^-3

        let bsurf = 1e11; // Tesla
        let btsl = self.bg * 1e-4; // Gauss -> Tesla
        let bdd = bsurf + btsl * (1.0 - (-BDD_BETAA * (nbtd / N0).powf(BDD_ALPHAA)).exp());

        let bdd_nlem = self.nlem.effective_bg(bdd);
        
        let bdd_nlem_gauss = bdd_nlem * 1e4;
        let b0 = bdd_nlem_gauss / BCE_G; 
        
        // Atualiza o engine com as unidades naturais (MeV^2 ou fm^-2) do projeto
        self.b = b0 * BCE;

        // 3) Densidades com B local atualizado
        crate::core::particles::calculate_all_densities(self, vomega, vrho);

        // 4) Equações de campo
        let fsigma = self.equation_sigma(vsigma);
        let fomega = self.equation_omega(vomega);
        let frho = self.equation_rho(vrho);
        let charge_neutral = self.charge_neutrality();

        vec![fsigma, fomega, frho, charge_neutral]
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
        let tolerance = 1e-10;
        let max_iterations = 100;
        let mut converged = false;

        for _ in 0..max_iterations {
            // Estado congelado para Jacobiano consistente
            let base_state = self.clone();

            let mut eval_state = base_state.clone();
            let f_val_vec = eval_state.funcv(x.as_slice());
            let f_val = Vector4::from_column_slice(&f_val_vec);
            let f_norm = f_val.norm();

            if f_norm < tolerance {
                *self = eval_state;
                converged = true;
                break;
            }

            // Jacobiana por diferenças finitas
            let mut j_matrix = Matrix4::zeros();
            for i in 0..4 {
                let h = 1e-8 * (x[i].abs() + 1e-2);
                let mut x_temp = x;
                x_temp[i] += h;

                let mut temp_state = base_state.clone();
                let f_temp_vec = temp_state.funcv(x_temp.as_slice());
                let f_temp = Vector4::from_column_slice(&f_temp_vec);

                let col = (f_temp - f_val) / h;
                j_matrix.set_column(i, &col);
            }

            let delta_x = match j_matrix.lu().solve(&(-f_val)) {
                Some(step) => step,
                None => break,
            };

            // Backtracking line-search com estado congelado
            let mut alpha = 1.0;
            let mut step_accepted = false;

            for _ in 0..15 {
                let x_try = x + alpha * delta_x;

                let mut trial_state = base_state.clone();
                let f_new_vec = trial_state.funcv(x_try.as_slice());
                let f_new = Vector4::from_column_slice(&f_new_vec);
                let f_new_norm = f_new.norm();

                if f_new_norm.is_nan() {
                    alpha *= 0.5;
                    continue;
                }

                if f_new_norm < f_norm {
                    x = x_try;
                    *self = trial_state;
                    step_accepted = true;
                    break;
                }

                alpha *= 0.5;
            }

            if !step_accepted {
                x += 0.001 * delta_x;
                let mut forced_state = base_state.clone();
                let _ = forced_state.funcv(x.as_slice());
                *self = forced_state;
            }
        }

        if !converged {
            return None;
        }

        let x_final = [x[0], x[1], x[2], x[3]];

        let (mue, vsigma, vomega, vrho) = self.mapping(&x_final);
        let (ener, press) = crate::core::eos::compute(self, mue, vsigma, vomega, vrho);

        let nb_total = self.nb.iter().sum::<f64>();
        let nbtd = nb_total * (self.m_nuc / 197.32).powi(3);

        let factor_mev_fm3 = self.m_nuc * (self.m_nuc / 197.32).powi(3);
        let ener_conv = ener * factor_mev_fm3;
        let press_conv = press * factor_mev_fm3;

        let bsurf = 1e11;
        let btsl = self.bg * 1e-4;
        let bdd = bsurf + btsl * (1.0 - (-BDD_BETAA * (nbtd / N0).powf(BDD_ALPHAA)).exp());

        let ebsi_maxwell = bdd.powi(2) / (8.0 * std::f64::consts::PI * 1e-7);
        let ebsi_nlem = self.nlem.magnetic_energy(bdd, ebsi_maxwell);
        let ebsd = ebsi_nlem / 1.602176634e32;

        let pmag_effective = match self.topology {
            MagneticTopology::Isotropic => ebsd / 3.0,
            MagneticTopology::Anisotropic => ebsd,
        };

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + pmag_effective;

        if ener_final >= 0.0 && press_final >= 0.0 {
            let result = [
                nbtd / 0.153, //  0: n/n0 (adimensional)
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
                ebsd,         // 19: Densidade de energia magnética [MeV/fm^3]
                bdd,          // 20: Campo magnético local B(n) [T]
            ];
            Some((x_final, result))
        } else {
            None
        }
    }

}