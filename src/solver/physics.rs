// solver/physics.rs

use crate::solver::constants::{
    HBAR_C, M_NUCLEON, QE, ML, MB, AMMP, AMMN, BCE,
};
use crate::solver::model::ModelParams;
use nalgebra::{Matrix4, Vector4};

pub struct PhysicsEngine {
    // Parâmetros fixos
    pub model: ModelParams,
    pub bg: f64,
    pub b: f64,
    pub m_nuc: f64,
    pub hc: f64,
    pub qe: f64,
    pub ml: [f64; 2],
    pub mb: [f64; 8],
    pub xs: f64,
    pub xv: f64,
    pub ammp: f64,
    pub ammn: f64,
    pub bce: f64,

    // Limites do loop (podem ser ajustados)
    pub mun_inf: f64,
    pub mun_sup: f64,
    pub n_points: usize,

    // --- Estado mutável para o ponto atual ---
    // Potenciais químicos
    pub mun: f64,
    pub mue: f64,
    pub mup: f64,

    // Massas efetivas
    pub mns: f64,
    pub ml0s: f64,
    pub msms: f64,
    pub ms0s: f64,
    pub msps: f64,
    pub mxms: f64,
    pub mx0s: f64,

    // Potenciais químicos dos híperons
    pub mul0: f64,
    pub musm: f64,
    pub mus0: f64,
    pub musp: f64,
    pub muxm: f64,
    pub mux0: f64,

    // Densidades
    pub nb: [f64; 8],
    pub nl: [f64; 2],
    pub nbt: f64,               // densidade bariônica total

    // Densidade escalar total (para equação sigma)
    pub rhosb: f64,

    // Energias de Fermi e momentos (para EOS)
    pub efn: f64,
    pub efp: f64,
    pub efe: f64,
    pub efmu: f64,
    pub efsm: f64,
    pub efsp: f64,
    pub efxm: f64,
    pub rkfl0: f64,
    pub efl0: f64,
    pub rkfs0: f64,
    pub efs0: f64,
    pub rkfx0: f64,
    pub efx0: f64,

    // Arrays para níveis de Landau
    pub fpu: Vec<f64>,
    pub fpd: Vec<f64>,
    pub fe: Vec<f64>,
    pub fmu: Vec<f64>,
    pub fsm: Vec<f64>,
    pub fsp: Vec<f64>,
    pub fxm: Vec<f64>,

    // Contadores de níveis
    pub npu: usize,
    pub npd: usize,
    pub ne: usize,
    pub nu: usize,
    pub nsm: usize,
    pub nsp: usize,
    pub nxm: usize,
}

impl PhysicsEngine {
    pub fn new(model: ModelParams, bg: f64) -> Self {
        let m_nuc = M_NUCLEON;
        let hc = HBAR_C;
        let qe = QE;
        let ml = ML;
        let mb = MB;
        let xs = 0.7;
        let xv = 0.783;
        let ammp = AMMP;
        let ammn = AMMN;
        let bce = BCE;

        let b0 = bg / 4.41e13;
        let b = b0 * bce;

        let max_landau = 19999;

        PhysicsEngine {
            model,
            bg,
            b,
            m_nuc,
            hc,
            qe,
            ml,
            mb,
            xs,
            xv,
            ammp,
            ammn,
            bce,
            mun_inf: 0.92,
            mun_sup: 1.80,
            n_points: 1001,

            mun: 0.0,
            mue: 0.0,
            mup: 0.0,

            mns: 0.0,
            ml0s: 0.0,
            msms: 0.0,
            ms0s: 0.0,
            msps: 0.0,
            mxms: 0.0,
            mx0s: 0.0,

            mul0: 0.0,
            musm: 0.0,
            mus0: 0.0,
            musp: 0.0,
            muxm: 0.0,
            mux0: 0.0,

            nb: [0.0; 8],
            nl: [0.0; 2],
            nbt: 0.0,

            rhosb: 0.0,

            efn: 0.0,
            efp: 0.0,
            efe: 0.0,
            efmu: 0.0,
            efsm: 0.0,
            efsp: 0.0,
            efxm: 0.0,
            rkfl0: 0.0,
            efl0: 0.0,
            rkfs0: 0.0,
            efs0: 0.0,
            rkfx0: 0.0,
            efx0: 0.0,

            fpu: vec![0.0; max_landau],
            fpd: vec![0.0; max_landau],
            fe: vec![0.0; max_landau],
            fmu: vec![0.0; max_landau],
            fsm: vec![0.0; max_landau],
            fsp: vec![0.0; max_landau],
            fxm: vec![0.0; max_landau],

            npu: 0,
            npd: 0,
            ne: 0,
            nu: 0,
            nsm: 0,
            nsp: 0,
            nxm: 0,
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
        let vsigma = x[1].sin().powi(2);
        let vomega = x[2];
        let vrho = x[3];
        (mue, vsigma, vomega, vrho)
    }

    // Função de resíduo (chamada pelo solver numérico)
    pub fn funcv(&mut self, x: &[f64]) -> Vec<f64> {
        let (mue, vsigma, vomega, vrho) = self.mapping(x);
        self.mue = mue;
        self.mup = self.mun - mue;

        // massas efetivas
        self.mns = 1.0 - vsigma;
        self.ml0s = (1116.0 / M_NUCLEON) - self.xs * vsigma;
        self.msms = (1193.0 / M_NUCLEON) - self.xs * vsigma;
        self.ms0s = (1193.0 / M_NUCLEON) - self.xs * vsigma;
        self.msps = (1193.0 / M_NUCLEON) - self.xs * vsigma;
        self.mxms = (1318.0 / M_NUCLEON) - self.xs * vsigma;
        self.mx0s = (1318.0 / M_NUCLEON) - self.xs * vsigma;

        // potenciais químicos dos híperons
        self.mul0 = self.mun;
        self.musm = self.mun + mue;
        self.mus0 = self.mun;
        self.musp = self.mun - mue;
        self.muxm = self.mun + mue;
        self.mux0 = self.mun;

        // calcular densidades (chama funções em particles.rs)
        crate::solver::particles::calculate_all_densities(self, vomega, vrho);

        // equações de movimento
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

    fn equation_omega(&self, vomega: f64) -> f64 {
        let ddd = self.nb[0] + self.nb[1]
            + self.xv * (self.nb[2] + self.nb[3] + self.nb[4] + self.nb[5] + self.nb[6] + self.nb[7]);
        self.model.gv.powi(2) * ddd - vomega
    }

    fn equation_rho(&self, vrho: f64) -> f64 {
        let frhnuc = (self.nb[1] - self.nb[0]) / 2.0;
        let fhh = -self.nb[3] + self.nb[5] - self.nb[6] / 2.0 + self.nb[7] / 2.0;
        self.model.gr.powi(2) * (frhnuc + self.xv * fhh) - vrho
    }

    fn charge_neutrality(&self) -> f64 {
        let charge_had = self.nb[1] - self.nb[3] + self.nb[5] - self.nb[6];
        charge_had - self.nl[0] - self.nl[1]
    }

    // Resolve para um dado mun e chute inicial, retorna solução e resultado
    pub fn solve_point(&mut self, mun: f64, initial_x: &[f64]) -> Option<([f64; 4], [f64; 13])> {
        self.mun = mun;

        let mut x = Vector4::from_column_slice(initial_x);
        let tolerance = 1e-10;
        let max_iterations = 1000; // Newton puro converge muito rápido (geralmente < 10 passos)
        let mut converged = false;

        for _ in 0..max_iterations {
            let f_val_vec = self.funcv(x.as_slice());
            let f_val = Vector4::from_column_slice(&f_val_vec);

            if f_val.norm() < tolerance {
                converged = true;
                break;
            }

            // 1. Calcula a Matriz Jacobiana EXATA por diferenças finitas em CADA iteração
            let mut j_matrix = Matrix4::zeros();
            
            for i in 0..4 {
                // Passo dinâmico para evitar erros de ponto flutuante em variáveis pequenas
                let h = 1e-8 * (x[i].abs() + 1e-2); 
                let mut x_temp = x;
                x_temp[i] += h;
                
                let f_temp_vec = self.funcv(x_temp.as_slice());
                let f_temp = Vector4::from_column_slice(&f_temp_vec);
                
                let column_derivative = (f_temp - f_val) / h;
                j_matrix.set_column(i, &column_derivative);
            }

            // 2. Resolve o sistema linear J * Δx = -F usando Decomposição LU
            let delta_x = match j_matrix.lu().solve(&(-f_val)) {
                Some(step) => step,
                None => {
                    // Se cair aqui, a derivada é zero (matriz singular). 
                    // O Newton não consegue prosseguir.
                    break;
                }
            };

            // 3. Backtracking Line Search (Mecanismo de Segurança)
            let mut alpha = 1.0;
            let mut step_accepted = false;

            for _ in 0..15 {
                let x_try = x + alpha * delta_x;
                let f_new_vec = self.funcv(x_try.as_slice());
                let f_new = Vector4::from_column_slice(&f_new_vec);

                // Rejeita o passo se ele atirou o solver para uma área não física (gerando NaN)
                if f_new.norm().is_nan() {
                    alpha *= 0.5;
                    continue;
                }

                // Condição de Armijo simplificada (se o erro diminuiu, aceitamos o passo)
                if f_new.norm() < f_val.norm() {
                    x = x_try;
                    step_accepted = true;
                    break;
                }
                
                alpha *= 0.5;
            }

            // Se as 15 tentativas de corte de passo falharam, damos um micro-passo forçado
            if !step_accepted {
                x += 0.001 * delta_x;
            }
        }

        if !converged {
            return None;
        }

        let x_final = [x[0], x[1], x[2], x[3]];

        // Mapeamento e computação física usando o resultado convergido
        let (mue, vsigma, vomega, vrho) = self.mapping(&x_final);
        let (ener, press) = crate::solver::eos::compute(self, mue, vsigma, vomega, vrho);

        let nb_total = self.nb.iter().sum::<f64>();
        let nbtd = nb_total * (self.m_nuc / 197.32).powi(3);
        
        let ener_conv = ener * self.m_nuc * (self.m_nuc / 197.32).powi(3) / self.hc;
        let press_conv = press * self.m_nuc * (self.m_nuc / 197.32).powi(3) / self.hc;

        // Adição da Pressão/Energia do Campo Magnético
        let bsurf = 1e11; 
        let btsl = self.bg * 1e-4; 
        let betaa = 1e-2;
        let alphaa = 3.0;

        let bdd = bsurf + btsl * (1.0 - (-betaa * (nbtd / 0.153).powf(alphaa)).exp());
        let ebsi = bdd.powi(2) / (8.0 * std::f64::consts::PI * 1e-7); 
        let ebsd = ebsi / 3.161e34; 

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + ebsd;

        // Limiar da crosta: retorna None se a pressão ficar negativa
        if ener_final >= 0.0 && press_final >= 0.0 {
            let result = [
                nbtd / 0.153,
                ener_final,
                press_final,
                self.nl[0],
                self.nl[1],
                self.nb[0],
                self.nb[1],
                self.nb[2],
                self.nb[3],
                self.nb[4],
                self.nb[5],
                self.nb[6],
                self.nb[7],
            ];
            Some((x_final, result))
        } else {
            None
        }
    }
}

impl Clone for PhysicsEngine {
    fn clone(&self) -> Self {
        // Clona os arrays
        PhysicsEngine {
            fpu: self.fpu.clone(),
            fpd: self.fpd.clone(),
            fe: self.fe.clone(),
            fmu: self.fmu.clone(),
            fsm: self.fsm.clone(),
            fsp: self.fsp.clone(),
            fxm: self.fxm.clone(),
            ..*self
        }
    }
}