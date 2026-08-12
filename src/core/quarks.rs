// src/core/quarks.rs

use crate::core::constants::{BDD_ALPHAA, BDD_BETAA, HBAR_C, M_NUCLEON, MV, N0, RESULTS_SIZE};
use crate::core::physics::NlemModel;
use nalgebra::{Matrix2, Vector2};
use std::f64::consts::PI;

/// Calcula a termodinâmica de um gás de Fermi ideal livre
/// Retorna: (Densidade [fm^-3], Energia [MeV/fm^3], Pressão [MeV/fm^3])
fn fermi_gas_thermo(kf_mev: f64, m_mev: f64, g: f64) -> (f64, f64, f64) {
    let hc3 = HBAR_C.powi(3); // Fator de conversão de MeV^3 para fm^-3

    if kf_mev <= 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let n_fm3 = g * kf_mev.powi(3) / (6.0 * PI.powi(2) * hc3);

    if m_mev == 0.0 {
        let ener = g * kf_mev.powi(4) / (8.0 * PI.powi(2) * hc3);
        return (n_fm3, ener, ener / 3.0);
    }

    let x = kf_mev / m_mev;
    let x2 = x * x;
    let x_sqrt = (x2 + 1.0).sqrt();

    let coef = g * m_mev.powi(4) / (16.0 * PI.powi(2) * hc3);
    let ener = coef * (x * x_sqrt * (2.0 * x2 + 1.0) - (x + x_sqrt).ln());
    let press = (coef / 3.0) * (x * x_sqrt * (2.0 * x2 - 3.0) + 3.0 * (x + x_sqrt).ln());

    (n_fm3, ener, press)
}

#[derive(Clone)]
pub struct QuarksMatter {
    pub bag_constant: f64, // MeV/fm^3
    pub gv: f64,           // Acoplamento Vetorial
    pub mv: f64,           // Massa do Méson Vetorial (MeV)
    pub bg: f64,
    pub nlem: NlemModel,
    pub mun_inf: f64, // Agora em unidades normalizadas (ex: 1.0)
    pub mun_sup: f64, // Agora em unidades normalizadas (ex: 2.5)
    pub n_points: usize,

    pub eos_output: Option<String>,

    // Variáveis para o Newton-Raphson interno
    pub last_mue: f64,
    pub last_v0: f64,
}

impl QuarksMatter {
    pub fn new(bag_constant: f64, bg: f64) -> Self {
        let nlem = NlemModel::Maxwell;

        Self {
            bag_constant,
            gv: 0.0, // (2.67357438647 / HBAR_C) * M_NUCLEON
            mv: MV * M_NUCLEON,
            bg,
            nlem,
            // Valores normalizados (mu_n / M_N)
            mun_inf: 1.0,
            mun_sup: 2.5,
            n_points: 1201,
            eos_output: None,
            last_mue: 10.0,
            last_v0: 0.0,
        }
    }

    pub fn with_nlem(mut self, nlem: NlemModel) -> Self {
        self.nlem = nlem;
        self
    }

    pub fn with_limits(mut self, inf: f64, sup: f64) -> Self {
        self.mun_inf = inf;
        self.mun_sup = sup;
        self
    }

    pub fn with_points(mut self, n: usize) -> Self {
        self.n_points = n;
        self
    }

    pub fn with_eos_output<P: Into<String>>(mut self, path: P) -> Self {
        self.eos_output = Some(path.into());
        self
    }

    /// Resolve o ponto para um mun NORMALIZADO
    pub fn solve_point(&mut self, mun_norm: f64) -> Option<[f64; RESULTS_SIZE]> {
        let mun_mev = mun_norm * M_NUCLEON;

        // Massas baseadas no program7.f do Fortran
        let m_u = 5.0;
        let m_d = 7.0;
        let m_s = 150.0;
        let m_e = 0.51099907;
        let m_mu = 105.66;
        let hc3 = HBAR_C.powi(3);

        // 1. Newton-Raphson Robusto com Clamping para g_v alto
        let guesses = [self.last_mue, 10.0, 50.0, 100.0];
        let mut final_x = Vector2::new(0.0, 0.0);
        let mut converged = false;

        for &guess_mue in guesses.iter() {
            // Chute inicial conservador para V0
            let guess_v0 = if self.gv > 0.0 { mun_mev * 0.05 } else { 0.0 };
            let mut x = Vector2::new(guess_mue, guess_v0);

            for _ in 0..200 {
                // Aumentamos as iterações devido ao amortecimento
                let f_val = self.residuals(mun_mev, x[0], x[1], m_u, m_d, m_s, m_e, m_mu);
                if f_val.norm() < 1e-10 {
                    final_x = x;
                    converged = true;
                    break;
                }

                let mut j_matrix = Matrix2::zeros();
                for i in 0..2 {
                    let h = 1e-7 * (x[i].abs() + 1e-3);
                    let mut x_temp = x;
                    x_temp[i] += h;
                    let f_temp =
                        self.residuals(mun_mev, x_temp[0], x_temp[1], m_u, m_d, m_s, m_e, m_mu);
                    j_matrix.set_column(i, &((f_temp - f_val) / h));
                }

                if let Some(mut delta) = j_matrix.lu().solve(&(-f_val)) {
                    // PROTEÇÃO 1: Evita saltos absurdos da derivada (Clipping)
                    if delta[0].abs() > 20.0 {
                        delta[0] = 20.0 * delta[0].signum();
                    }
                    if delta[1].abs() > 20.0 {
                        delta[1] = 20.0 * delta[1].signum();
                    }

                    // PROTEÇÃO 2: Amortecimento agressivo para acoplamentos fortes
                    let damping = if self.gv > 2.0 {
                        0.05
                    } else if self.gv > 0.0 {
                        0.2
                    } else {
                        1.0
                    };
                    x += delta * damping;

                    // PROTEÇÃO 3: V0 nunca pode ser negativo (a repulsão é positiva)
                    if x[1] < 0.0 {
                        x[1] = 0.0;
                    }

                    // PROTEÇÃO 4 (A MAIS IMPORTANTE): V0 não pode aniquilar o plasma
                    // Garante que (mu - gv*V0) > massa do quark, senão a densidade zera e o solver quebra
                    if self.gv > 0.0 {
                        let mu_u_base = mun_mev / 3.0 - (2.0 / 3.0) * x[0];
                        if mu_u_base > m_u {
                            let max_v0 = (mu_u_base - m_u) / self.gv;
                            if x[1] >= max_v0 {
                                x[1] = max_v0 * 0.98; // Trava o V0 logo antes de matar a física
                            }
                        } else {
                            x[1] = 0.0;
                        }
                    }
                } else {
                    break; // Matriz singular, aborta e tenta o próximo chute
                }
            }
            if converged {
                break;
            }
        }

        if !converged {
            return None;
        }

        self.last_mue = final_x[0];
        self.last_v0 = final_x[1];

        // --- CÁLCULO FÍSICO ---
        let get_kf = |mu_total: f64, mass: f64, gv_v0: f64| {
            let mu_eff = mu_total - gv_v0;
            if mu_eff > mass {
                (mu_eff.powi(2) - mass.powi(2)).sqrt()
            } else {
                0.0
            }
        };

        let kf_u = get_kf(
            mun_mev / 3.0 - (2.0 / 3.0) * final_x[0],
            m_u,
            self.gv * final_x[1],
        );
        let kf_d = get_kf(
            mun_mev / 3.0 + (1.0 / 3.0) * final_x[0],
            m_d,
            self.gv * final_x[1],
        );
        let kf_s = get_kf(
            mun_mev / 3.0 + (1.0 / 3.0) * final_x[0],
            m_s,
            self.gv * final_x[1],
        );
        let kf_e = if final_x[0] > m_e {
            (final_x[0].powi(2) - m_e.powi(2)).sqrt()
        } else {
            0.0
        };
        let kf_mu = if final_x[0] > m_mu {
            (final_x[0].powi(2) - m_mu.powi(2)).sqrt()
        } else {
            0.0
        };

        let mu_u = mun_mev / 3.0 - (2.0 / 3.0) * final_x[0];
        let mu_d = mun_mev / 3.0 + (1.0 / 3.0) * final_x[0];
        let mu_s = mu_d;
        let mu_l = final_x[0];

        let (n_u, e_u, p_u) = fermi_gas_thermo(kf_u, m_u, 6.0);
        let (n_d, e_d, p_d) = fermi_gas_thermo(kf_d, m_d, 6.0);
        let (n_s, e_s, p_s) = fermi_gas_thermo(kf_s, m_s, 6.0);
        let (n_e, e_e, p_e) = fermi_gas_thermo(kf_e, m_e, 2.0);
        let (n_mu, e_mu, p_mu) = fermi_gas_thermo(kf_mu, m_mu, 2.0);

        let nb_total = (n_u + n_d + n_s) / 3.0;
        let meson_energy = 0.5 * self.mv.powi(2) * final_x[1].powi(2) / hc3;

        let ener_conv = e_u + e_d + e_s + e_e + e_mu + self.bag_constant + meson_energy;
        let press_conv = p_u + p_d + p_s + p_e + p_mu - self.bag_constant + meson_energy;

        // --- CONTRIBUIÇÃO MAGNÉTICA (NLEM) ---
        let bsurf = 1e11;
        let btsl = self.bg * 1e-4;

        let bdd = if self.bg == 0.0 {
            0.0
        } else {
            bsurf + btsl * (1.0 - (-BDD_BETAA * (nb_total / N0).powf(BDD_ALPHAA)).exp())
        };
        let ebsi_maxwell = bdd.powi(2) / (8.0 * PI * 1e-7);
        let ebsi_nlem = self.nlem.magnetic_energy(bdd, ebsi_maxwell);
        let ebsd = ebsi_nlem / 1.602176634e32;

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + ebsd;

        // Tratamento para evitar pressões negativas
        if press_final <= 0.0 {
            // Enquanto a pressão do gás não vencer a Constante de Sacola, o estado é de vácuo.
            let mut vac = [0.0; RESULTS_SIZE];
            vac[17] = mun_norm;
            return Some(vac);
        }

        if ener_final >= 0.0 && press_final >= 0.0 {
            let mut result = [0.0; RESULTS_SIZE];
            let fermion_mu_density = mu_u * n_u + mu_d * n_d + mu_s * n_s + mu_l * (n_e + n_mu);
            let mu_total_per_baryon = if nb_total > 0.0 {
                fermion_mu_density / nb_total / M_NUCLEON
            } else {
                0.0
            };

            result[0] = nb_total / N0;
            result[1] = ener_final;
            result[2] = press_final;
            result[3] = n_e;
            result[4] = n_mu;
            result[5] = n_u;
            result[6] = n_d;
            result[7] = n_s;
            result[14] = final_x[1];
            result[17] = mun_norm;
            result[18] = final_x[0];
            result[19] = ebsd;
            result[20] = mu_total_per_baryon;
            Some(result)
        } else {
            None
        }
    }

    fn residuals(
        &self,
        mun_mev: f64,
        mue: f64,
        v0: f64,
        mu: f64,
        md: f64,
        ms: f64,
        me: f64,
        mmu: f64,
    ) -> Vector2<f64> {
        let hc3 = HBAR_C.powi(3); // Atualizado para usar a constante global

        // Função auxiliar para densidade com proteção física mu_eff > mass
        let get_n = |mu_total: f64, mass: f64, gv_v0: f64, g: f64| {
            let mu_eff = mu_total - gv_v0;
            // IMPORTANTE: mu_eff deve ser positivo e maior que a massa
            if mu_eff > mass {
                let kf = (mu_eff.powi(2) - mass.powi(2)).sqrt();
                g * kf.powi(3) / (6.0 * PI.powi(2) * hc3)
            } else {
                0.0
            }
        };

        let n_u = get_n(mun_mev / 3.0 - 2.0 / 3.0 * mue, mu, self.gv * v0, 6.0);
        let n_d = get_n(mun_mev / 3.0 + 1.0 / 3.0 * mue, md, self.gv * v0, 6.0);
        let n_s = get_n(mun_mev / 3.0 + 1.0 / 3.0 * mue, ms, self.gv * v0, 6.0);

        let n_e = if mue > me {
            let kf = (mue.powi(2) - me.powi(2)).sqrt();
            2.0 * kf.powi(3) / (6.0 * PI.powi(2) * hc3)
        } else {
            0.0
        };

        let n_mu = if mue > mmu {
            let kf = (mue.powi(2) - mmu.powi(2)).sqrt();
            2.0 * kf.powi(3) / (6.0 * PI.powi(2) * hc3)
        } else {
            0.0
        };

        // 1. Carga do Plasma = 0
        let charge = (2.0 / 3.0) * n_u - (1.0 / 3.0) * n_d - (1.0 / 3.0) * n_s - n_e - n_mu;

        // 2. Equação do Campo Vetorial NORMALIZADA
        // Dividimos por mv^2 para que o resíduo seja o próprio erro em MeV.
        let eom_v = (self.gv * (n_u + n_d + n_s) * hc3) / self.mv.powi(2) - v0;

        Vector2::new(charge, eom_v)
    }
}
