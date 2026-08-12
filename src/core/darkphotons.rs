// src/core/darkphotons.rs

use crate::core::constants::{
    AMML0, AMMN, AMMP, AMMS0, AMMSM, AMMSP, AMMX0, AMMXM, BCE, BCE_G, BDD_ALPHAA, BDD_BETAA,
    HBAR_C, M_NUCLEON, MAX_LANDAU_LIMIT, MB, ML, N0, PI2, QE, RESULTS_SIZE,
};
use crate::core::model::ModelParams;
use crate::core::physics::MagneticTopology;
use nalgebra::{SMatrix, SVector};

type DarkStateVector = SVector<f64, 5>;
type DarkJacobian = SMatrix<f64, 5, 5>;

#[derive(Clone)]
pub struct DarkPhotonsMatter {
    // Parâmetros fixos
    pub model: ModelParams,
    pub bg: f64,
    pub b: f64,
    /// Kinetic-mixing parameter (dimensionless, |epsilon| < 1).
    pub epsilon: f64,
    /// Physical dark-photon mass divided by M_NUCLEON.
    pub m_x: f64,
    /// Dark U(1) coupling (dimensionless).
    pub g_d: f64,
    /// Physical Dirac-fermion mass divided by M_NUCLEON.
    pub m_chi: f64,
    /// Imposed number fraction n_chi / n_B (dimensionless).
    pub y_chi: f64,

    // Estado do gás de Fermi escuro (unidades naturais normalizadas por M_NUCLEON)
    /// Dark number density in internal natural units (M_NUCLEON^3).
    pub n_chi: f64,
    /// Dark Fermi momentum divided by M_NUCLEON.
    pub kf_chi: f64,
    /// Dark effective Fermi energy divided by M_NUCLEON.
    pub ef_chi: f64,
    /// Full dark chemical potential divided by M_NUCLEON.
    pub mu_chi: f64,
    /// Dark kinetic energy density in internal units (M_NUCLEON^4).
    pub ener_chi_kin: f64,
    /// Dark kinetic pressure in internal units (M_NUCLEON^4).
    pub press_chi_kin: f64,
    /// Mean dark-photon potential divided by M_NUCLEON.
    pub v_x0: f64,
    pub topology: MagneticTopology,
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
    pub nbt: f64, // densidade bariônica total

    // Densidades escalares
    pub rhosb: f64,
    pub rhos_b: [f64; 8],
    pub rhos_l: [f64; 2],

    // Energias de Fermi e momentos (para EOS)
    pub ef_b: [f64; 8],
    pub ef_l: [f64; 2], // Energias de Fermi: [0]=e, [1]=mu

    // Acoplamentos (xv_v para omega, xv_r para rho)
    pub xv_v: [f64; 8], // g_wB / g_wN
    pub xv_r: [f64; 8], // g_rB / g_rN

    // Momentos de Fermi por nível de Landau (Vetorizados)
    pub kf_b_up: [Vec<f64>; 8], // [Barião][Nível nu]
    pub kf_b_down: [Vec<f64>; 8],
    pub f_l: [Vec<f64>; 2], // Momentos de Fermi: [0]=fe, [1]=fmu

    // Contadores de níveis por spin
    pub n_b_up: [usize; 8],
    pub n_b_down: [usize; 8],
    pub n_l: [usize; 2], // Contadores: [0]=ne, [1]=nu

    pub max_landau_limit: usize,

    pub isospin_factor: [f64; 8],

    pub eos_output: Option<String>,
}

impl DarkPhotonsMatter {
    // constante estática para acoplamento sigma
    const X_SIGMA: [f64; 8] = [1.0, 1.0, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];

    fn density_natural_to_fm3(density: f64) -> f64 {
        density * (M_NUCLEON / HBAR_C).powi(3)
    }

    fn kinetic_mixing_norm(epsilon: f64) -> f64 {
        assert!(
            epsilon.is_finite() && epsilon.abs() < 1.0,
            "dark photon kinetic mixing epsilon must be finite and satisfy |epsilon| < 1"
        );
        (1.0 - epsilon.powi(2)).sqrt()
    }

    fn dark_shift_for_charge(&self, charge_units: f64) -> f64 {
        self.epsilon * charge_units * self.qe * self.v_x0 / Self::kinetic_mixing_norm(self.epsilon)
    }

    pub fn new(model: ModelParams, bg: f64) -> Self {
        let m_nuc = M_NUCLEON;
        let qe = QE;
        let ml = ML;
        let mb = MB;
        let xs = 0.7;
        let b0 = bg / BCE_G;
        let b = b0 * BCE;
        let epsilon = 0.0;
        let m_x = 1.0;
        let g_d = 0.0;
        let m_chi = 1.0;
        let y_chi = 0.0;
        let n_chi = 0.0;
        let kf_chi = 0.0;
        let ef_chi = 0.0;
        let mu_chi = 0.0;
        let ener_chi_kin = 0.0;
        let press_chi_kin = 0.0;
        let v_x0 = 0.0;
        let topology = MagneticTopology::Anisotropic;

        let m_eff = [0.0; 8];
        let mu_b = [0.0; 8];
        let charges_b = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0];

        let amm_b = [AMMN, AMMP, AMML0, AMMSM, AMMS0, AMMSP, AMMXM, AMMX0];

        let xv_v = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];
        let xv_r = [1.0, 1.0, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783];

        let kf_b_up = std::array::from_fn(|_| Vec::new());
        let kf_b_down = std::array::from_fn(|_| Vec::new());
        let f_l = std::array::from_fn(|_| Vec::new());

        let ef_l = [0.0; 2];
        let n_l = [0; 2];

        let isospin_factor = [-0.5, 0.5, 0.0, -1.0, 0.0, 1.0, -0.5, 0.5];

        DarkPhotonsMatter {
            model,
            bg,
            b,
            epsilon,
            m_x,
            g_d,
            m_chi,
            y_chi,
            n_chi,
            kf_chi,
            ef_chi,
            mu_chi,
            ener_chi_kin,
            press_chi_kin,
            v_x0,
            topology,
            m_nuc,
            qe,
            ml,
            mb,
            m_eff,
            mu_b,
            charges_b,
            amm_b,
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

            n_b_up: [0; 8],
            n_b_down: [0; 8],
            n_l,

            max_landau_limit: MAX_LANDAU_LIMIT,

            isospin_factor,
            eos_output: None,
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

    pub fn with_eos_output<P: Into<String>>(mut self, path: P) -> Self {
        self.eos_output = Some(path.into());
        self
    }

    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        // dimensionless kinetic mixing
        Self::kinetic_mixing_norm(epsilon);
        self.epsilon = epsilon;
        self
    }

    pub fn with_m_x(mut self, m_x: f64) -> Self {
        // dark photon mass in units of M_NUCLEON (same convention as mb/ml)
        assert!(
            m_x.is_finite() && m_x > 0.0,
            "dark photon mass m_x must be finite and positive"
        );
        self.m_x = m_x;
        self
    }

    pub fn with_m_x_mev(self, m_x_mev: f64) -> Self {
        assert!(
            m_x_mev.is_finite() && m_x_mev > 0.0,
            "dark photon mass m_x must be finite and positive"
        );
        self.with_m_x(m_x_mev / M_NUCLEON)
    }

    pub fn with_g_d(mut self, g_d: f64) -> Self {
        // dimensionless dark coupling
        assert!(g_d.is_finite(), "dark coupling g_d must be finite");
        self.g_d = g_d;
        self
    }

    pub fn with_m_chi(mut self, m_chi: f64) -> Self {
        // Dirac dark-fermion mass in units of M_NUCLEON.
        assert!(
            m_chi.is_finite() && m_chi > 0.0,
            "dark fermion mass m_chi must be finite and positive"
        );
        self.m_chi = m_chi;
        self
    }

    pub fn with_m_chi_mev(self, m_chi_mev: f64) -> Self {
        assert!(
            m_chi_mev.is_finite() && m_chi_mev > 0.0,
            "dark fermion mass m_chi must be finite and positive"
        );
        self.with_m_chi(m_chi_mev / M_NUCLEON)
    }

    pub fn with_y_chi(mut self, y_chi: f64) -> Self {
        assert!(
            y_chi.is_finite() && y_chi >= 0.0,
            "dark number fraction y_chi must be finite and non-negative"
        );
        self.y_chi = y_chi;
        self
    }

    pub fn with_topology(mut self, top: MagneticTopology) -> Self {
        self.topology = top;
        self
    }

    // Mapeamento das variáveis (vindo do solver)
    fn mapping(&self, x: &[f64]) -> (f64, f64, f64, f64, f64) {
        assert!(
            x.len() >= 5,
            "dark-sector solver state must have five variables"
        );
        let mue = x[0];
        let vsigma = x[1]; // Removido o .sin().powi(2) que destruía o Jacobiano
        let vomega = x[2];
        let vrho = x[3];
        let v_x0 = x[4];
        (mue, vsigma, vomega, vrho, v_x0)
    }

    // Função de resíduo (chamada pelo solver numérico)
    fn funcv(&mut self, x: &[f64]) -> [f64; 5] {
        let (mue, vsigma, vomega, vrho, v_x0) = self.mapping(x);

        self.mue = mue;
        self.v_x0 = v_x0;
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

        self.n_chi = self.y_chi * self.nbt;
        self.update_dark_fermion_state();

        let fsigma = self.equation_sigma(vsigma);
        let fomega = self.equation_omega(vomega, vrho);
        let frho = self.equation_rho(vrho, vomega);
        let charge_neutral = self.charge_neutrality();
        // Scale the Proca row by m_X^2 so the global Newton tolerance controls
        // the absolute error in X_0 even for light mediators.
        let fdark = self.dark_photon_residual(charge_neutral) / self.m_x.powi(2);

        [fsigma, fomega, frho, charge_neutral, fdark]
    }

    fn update_dark_fermion_state(&mut self) {
        if self.n_chi <= 0.0 {
            self.n_chi = 0.0;
            self.kf_chi = 0.0;
            self.ef_chi = 0.0;
            self.mu_chi = 0.0;
            self.ener_chi_kin = 0.0;
            self.press_chi_kin = 0.0;
            return;
        }

        self.kf_chi = dark_fermion_kf_from_density(self.n_chi);
        self.ef_chi = (self.kf_chi.powi(2) + self.m_chi.powi(2)).sqrt();
        self.mu_chi = self.ef_chi + self.g_d * self.v_x0 / Self::kinetic_mixing_norm(self.epsilon);
        self.ener_chi_kin = dark_fermion_energy_density(self.kf_chi, self.m_chi);
        self.press_chi_kin = dark_fermion_pressure(self.kf_chi, self.m_chi);
    }

    fn dark_photon_residual(&self, charge_density: f64) -> f64 {
        self.m_x.powi(2) * self.v_x0
            - (self.g_d * self.n_chi + self.epsilon * self.qe * charge_density)
                / Self::kinetic_mixing_norm(self.epsilon)
    }

    fn dark_vector_energy_density(&self) -> f64 {
        0.5 * self.m_x.powi(2) * self.v_x0.powi(2)
    }

    fn equation_sigma(&self, vsigma: f64) -> f64 {
        let gs2 = self.model.gs.powi(2);
        gs2 * (self.rhosb - self.model.rb * vsigma.powi(2) - self.model.rc * vsigma.powi(3))
            - vsigma
    }

    // Equações de Campo Vetorizadas para suportar as partículas com total precisão
    fn equation_omega(&self, vomega: f64, vrho: f64) -> f64 {
        // Legacy RMF normalization retained for parity with HadronsMatter.
        // For nonzero rxi/lambda_v (FSU2), these nonlinear terms are not the
        // exact derivative of the energy functional used by compute(); see
        // docs/PHYSICS.md before changing both visible-matter paths together.
        let mut sum_baryon = 0.0;
        for i in 0..8 {
            sum_baryon += self.nb[i] * self.xv_v[i];
        }
        self.model.gv.powi(2) * sum_baryon
            - vomega
            - self.model.rxi * vomega.powi(3)
            - 2.0 * self.model.lambda_v * vomega * vrho.powi(2)
    }

    fn equation_rho(&self, vrho: f64, vomega: f64) -> f64 {
        let mut sum_source = 0.0;
        for i in 0..8 {
            // A fonte para o rho é baseada no negativo do isospin
            sum_source += self.isospin_factor[i] * self.nb[i] * self.xv_r[i];
        }
        self.model.gr.powi(2) * sum_source
            - vrho
            - 2.0 * self.model.lambda_v * vrho * vomega.powi(2)
    }

    fn charge_neutrality(&self) -> f64 {
        let charge_baryons: f64 = self
            .nb
            .iter()
            .zip(self.charges_b.iter())
            .map(|(n, q)| n * q)
            .sum();

        // Leptons: e⁻ and μ⁻ have charge -1
        let charge_leptons: f64 = self.nl.iter().map(|n| -n).sum();

        charge_baryons + charge_leptons
    }

    // Resolve para um dado mun e chute inicial, retorna solução e resultado
    pub fn solve_point(
        &mut self,
        mun: f64,
        initial_x: &[f64; 5],
    ) -> Option<([f64; 5], [f64; RESULTS_SIZE])> {
        self.mun = mun;

        let mut x = DarkStateVector::from_column_slice(initial_x);
        let tolerance = 1e-10;
        let max_iterations = 100;
        let mut converged = false;

        for _ in 0..max_iterations {
            let f_val_arr = self.funcv(x.as_slice());
            let f_val = DarkStateVector::from_column_slice(&f_val_arr);
            let f_norm = f_val.norm();

            if f_norm.is_finite() && f_norm < tolerance {
                converged = true;
                break;
            }

            if !f_norm.is_finite() {
                break;
            }

            let mut j_matrix = DarkJacobian::zeros();

            for i in 0..5 {
                let h = f64::EPSILON.sqrt() * (1.0 + x[i].abs());
                let mut x_temp = x;
                x_temp[i] += h;

                let f_temp_arr = self.funcv(x_temp.as_slice());
                let f_temp = DarkStateVector::from_column_slice(&f_temp_arr);

                let column_derivative = (f_temp - f_val) / h;
                j_matrix.set_column(i, &column_derivative);
            }

            let delta_x = match j_matrix.lu().solve(&(-f_val)) {
                Some(step) if step.iter().all(|value| value.is_finite()) => step,
                None => break,
                Some(_) => break,
            };

            let mut alpha = 1.0;
            let mut step_accepted = false;

            for _ in 0..15 {
                let x_try = x + alpha * delta_x;
                let f_new_arr = self.funcv(x_try.as_slice());
                let f_new = DarkStateVector::from_column_slice(&f_new_arr);
                let f_new_norm = f_new.norm();

                if !f_new_norm.is_finite() {
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
                break;
            }
        }

        // A root reached by the final permitted Newton step must not be
        // discarded merely because the next loop-head check cannot run.
        if !converged {
            let final_residual = DarkStateVector::from_column_slice(&self.funcv(x.as_slice()));
            let final_norm = final_residual.norm();
            converged = final_norm.is_finite() && final_norm < tolerance;
        }

        if !converged {
            return None;
        }

        // garante estado físico final consistente
        let _ = self.funcv(x.as_slice());

        let x_final = [x[0], x[1], x[2], x[3], x[4]];
        let (mue, vsigma, vomega, vrho, _) = self.mapping(&x_final);
        let (ener, press) = compute(self, mue, vsigma, vomega, vrho);

        let nb_total = self.nb.iter().sum::<f64>();
        let nbtd = nb_total * (self.m_nuc / HBAR_C).powi(3);

        let factor_mev_fm3 = self.m_nuc * (self.m_nuc / HBAR_C).powi(3);
        let ener_conv = ener * factor_mev_fm3;
        let press_conv = press * factor_mev_fm3;

        let bsurf = 1e11;
        let btsl = self.bg * 1e-4;
        let bdd = bsurf + btsl * (1.0 - (-BDD_BETAA * (nbtd / N0).powf(BDD_ALPHAA)).exp());

        let ebsi_maxwell = bdd.powi(2) / (8.0 * std::f64::consts::PI * 1e-7);
        let ebsd = ebsi_maxwell / 1.602176634e32;

        let pmag_effective = match self.topology {
            MagneticTopology::Isotropic => ebsd / 3.0,
            MagneticTopology::Anisotropic => ebsd,
        };

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + pmag_effective;

        if ener_final >= 0.0 && press_final >= 0.0 {
            let fermion_mu_density = self
                .mu_b
                .iter()
                .zip(self.nb.iter())
                .map(|(mu, n)| mu * n)
                .sum::<f64>()
                + mue * self.nl.iter().sum::<f64>()
                + self.mu_chi * self.n_chi;
            let mu_total_per_baryon = if self.nbt > 0.0 {
                fermion_mu_density / self.nbt
            } else {
                0.0
            };

            let density_factor = (self.m_nuc / HBAR_C).powi(3);
            let dark_vector_energy = self.dark_vector_energy_density();
            let mut result = [0.0; RESULTS_SIZE];
            result[0] = nbtd / N0;
            result[1] = ener_final;
            result[2] = press_final;
            result[3] = self.nl[0] * density_factor;
            result[4] = self.nl[1] * density_factor;
            for i in 0..8 {
                result[5 + i] = self.nb[i] * density_factor;
            }
            result[13] = vsigma * self.m_nuc;
            result[14] = vomega * self.m_nuc;
            result[15] = vrho * self.m_nuc;
            result[16] = self.m_eff[0];
            result[17] = self.mun;
            result[18] = mue;
            result[19] = ebsd;
            result[20] = mu_total_per_baryon;
            result[21] = Self::density_natural_to_fm3(self.n_chi);
            result[22] = self.y_chi;
            result[23] = self.m_chi * self.m_nuc;
            result[24] = self.m_x * self.m_nuc;
            result[25] = self.epsilon;
            result[26] = self.g_d;
            result[27] = self.v_x0 * self.m_nuc;
            result[28] = self.kf_chi * self.m_nuc;
            result[29] = self.mu_chi * self.m_nuc;
            result[30] = self.ener_chi_kin * factor_mev_fm3;
            result[31] = self.press_chi_kin * factor_mev_fm3;
            result[32] = dark_vector_energy * factor_mev_fm3;
            result[33] = dark_vector_energy * factor_mev_fm3;
            Some((x_final, result))
        } else {
            None
        }
    }
}

/// Number density of a zero-temperature spin-1/2 Dirac gas in internal natural units.
pub fn dark_fermion_number_density_from_kf(kf: f64) -> f64 {
    if kf <= 0.0 {
        0.0
    } else {
        kf.powi(3) / (3.0 * PI2)
    }
}

/// Fermi momentum of a zero-temperature spin-1/2 Dirac gas in internal natural units.
pub fn dark_fermion_kf_from_density(density: f64) -> f64 {
    if density <= 0.0 {
        0.0
    } else {
        (3.0 * PI2 * density).cbrt()
    }
}

/// Kinetic (including rest mass) energy density of a free Dirac gas.
pub fn dark_fermion_energy_density(kf: f64, mass: f64) -> f64 {
    assert!(
        mass.is_finite() && mass > 0.0,
        "dark fermion mass must be finite and positive"
    );
    if kf <= 0.0 {
        return 0.0;
    }

    let x = kf / mass;
    let dimensionless = if x.abs() < 1e-2 {
        // Series avoids cancellation between the algebraic and asinh terms.
        (8.0 / 3.0) * x.powi(3) + (4.0 / 5.0) * x.powi(5) - x.powi(7) / 7.0 + x.powi(9) / 18.0
    } else {
        let ef_over_m = x.hypot(1.0);
        x * ef_over_m * (2.0 * x.powi(2) + 1.0) - x.asinh()
    };

    mass.powi(4) * dimensionless / (8.0 * PI2)
}

/// Kinetic pressure of a zero-temperature spin-1/2 Dirac gas.
pub fn dark_fermion_pressure(kf: f64, mass: f64) -> f64 {
    assert!(
        mass.is_finite() && mass > 0.0,
        "dark fermion mass must be finite and positive"
    );
    if kf <= 0.0 {
        return 0.0;
    }

    let x = kf / mass;
    let dimensionless = if x.abs() < 1e-2 {
        (8.0 / 5.0) * x.powi(5) - (4.0 / 7.0) * x.powi(7) + x.powi(9) / 3.0
    } else {
        let ef_over_m = x.hypot(1.0);
        x * ef_over_m * (2.0 * x.powi(2) - 3.0) + 3.0 * x.asinh()
    };

    mass.powi(4) * dimensionless / (24.0 * PI2)
}

pub fn calculate_all_densities(engine: &mut DarkPhotonsMatter, vomega: f64, vrho: f64) {
    for i in 0..8 {
        let (rs, rb) = if engine.charges_b[i] == 0.0 {
            density_baryon_neutral(engine, i, vomega, vrho)
        } else {
            density_baryon_charged(engine, i, vomega, vrho)
        };
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

pub fn density_baryon_neutral(
    engine: &mut DarkPhotonsMatter,
    idx: usize,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    let ef = engine.mu_b[idx]
        - (engine.xv_v[idx] * vomega)
        - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx]);

    engine.ef_b[idx] = ef;

    // ZERA OS MOMENTOS PARA EVITAR "FANTASMAS" DO NEWTON-RAPHSON
    engine.kf_b_up[idx].clear();
    engine.kf_b_down[idx].clear();

    if ef <= 0.0 {
        return (0.0, 0.0);
    }

    let m_star = engine.m_eff[idx];
    let amm = engine.amm_b[idx];
    let b = engine.b;

    let m_up = m_star - amm * b;
    let m_down = m_star + amm * b;

    let mut rhos_total = 0.0;
    let mut dens_total = 0.0;

    // Quando B=0, m_up == m_down
    let spins = [m_up, m_down];
    for (spin_idx, &m_spin) in spins.iter().enumerate() {
        let kf2 = ef.powi(2) - m_spin.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let m_safe = m_spin.abs().max(1e-15);

            rhos_total +=
                (m_spin / (4.0 * PI2)) * (ef * kf - m_spin.powi(2) * ((kf + ef) / m_safe).ln());

            dens_total += kf.powi(3) / (6.0 * PI2);

            // Grava o kf no buffer correto usando push()
            if spin_idx == 0 {
                engine.kf_b_up[idx].push(kf);
            } else {
                engine.kf_b_down[idx].push(kf);
            }
        }
    }
    (rhos_total, dens_total)
}

fn density_baryon_charged(
    engine: &mut DarkPhotonsMatter,
    idx: usize,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    let q = engine.charges_b[idx].abs() * engine.qe;
    let b = engine.b;
    let m = engine.m_eff[idx];
    let amm = engine.amm_b[idx];
    let dark_shift = engine.dark_shift_for_charge(engine.charges_b[idx]);

    let ef = engine.mu_b[idx]
        - (engine.xv_v[idx] * vomega)
        - (engine.xv_r[idx] * vrho * engine.isospin_factor[idx])
        - dark_shift;

    engine.ef_b[idx] = ef;

    // ZERA OS MOMENTOS PARA EVITAR FANTASMAS
    engine.n_b_up[idx] = 0;
    engine.n_b_down[idx] = 0;
    engine.kf_b_up[idx].clear();
    engine.kf_b_down[idx].clear();

    if ef <= 0.0 {
        return (0.0, 0.0);
    }

    // TRATAMENTO PARA O CASO ISOTRÓPICO (B=0)
    if b == 0.0 {
        let kf2 = ef.powi(2) - m.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let dens = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos = (m / (2.0 * PI2)) * (ef * kf - m.powi(2) * ((kf + ef) / m_safe).ln());

            engine.kf_b_up[idx].push(kf);
            engine.kf_b_down[idx].push(kf);
            engine.n_b_up[idx] = 1;
            engine.n_b_down[idx] = 1;

            return (rhos, dens);
        }
        return (0.0, 0.0);
    }

    let nu_max_approx_up = ((ef + amm * b).powi(2) - m.powi(2)) / (2.0 * q * b);
    let nu_max_approx_down = ((ef - amm * b).powi(2) - m.powi(2)) / (2.0 * q * b);

    let nu_max = if nu_max_approx_up > 0.0 || nu_max_approx_down > 0.0 {
        let max_nu = nu_max_approx_up.max(nu_max_approx_down);
        (max_nu.floor() as usize + 1).min(engine.max_landau_limit)
    } else {
        0
    };

    let q_sign = engine.charges_b[idx].signum();
    let (nu_start_up, nu_start_down) = if q_sign > 0.0 {
        (0, 1) // Cargas positivas: Spin UP tem nu=0
    } else {
        (1, 0) // Cargas negativas: Spin DOWN tem nu=0
    };

    let (mut rhos, mut dens) = (0.0, 0.0);
    let mut n_up = 0;
    let mut n_down = 0;

    engine.kf_b_up[idx].reserve(nu_max);
    engine.kf_b_down[idx].reserve(nu_max);

    // --- Spin UP ---
    for nu in nu_start_up..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau - amm * b;

        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        engine.kf_b_up[idx].push(kf);
        n_up += 1;

        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15);

        rhos +=
            (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    // --- Spin DOWN ---
    for nu in nu_start_down..nu_max {
        let m_landau = (m.powi(2) + 2.0 * q * b * nu as f64).sqrt();
        let m_eff_spin = m_landau + amm * b;

        let kf2 = ef.powi(2) - m_eff_spin.powi(2);
        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        engine.kf_b_down[idx].push(kf);
        n_down += 1;

        let m_safe = m_eff_spin.abs().max(1e-15);
        let m_landau_safe = m_landau.max(1e-15);

        rhos +=
            (q * b / (2.0 * PI2)) * m * (m_eff_spin / m_landau_safe) * ((kf + ef) / m_safe).ln();
        dens += (q * b / (2.0 * PI2)) * kf;
    }

    engine.n_b_up[idx] = n_up;
    engine.n_b_down[idx] = n_down;

    (rhos, dens)
}

pub fn density_lepton(engine: &mut DarkPhotonsMatter, idx: usize) -> (f64, f64) {
    let mue = engine.mue;

    // Carga do lépton (elétron/múon) é -1 em unidades de carga fundamental e
    let q_lepton = -1.0;

    // Cálculo do acoplamento do lépton com o fóton escuro devido à mistura cinética
    let dark_shift = engine.dark_shift_for_charge(q_lepton);

    // O potencial químico efetivo (energia de Fermi) que rege o preenchimento dos estados
    let ef_lepton = mue - dark_shift;

    // ZERA OS ESTADOS PARA EVITAR FANTASMAS
    engine.n_l[idx] = 0;
    engine.f_l[idx].clear();
    engine.ef_l[idx] = ef_lepton;

    if ef_lepton <= 0.0 {
        return (0.0, 0.0);
    }

    let b = engine.b;
    let q = engine.qe;
    let m = engine.ml[idx];

    let mut rhos = 0.0;
    let mut dens = 0.0;
    let mut n_occupied = 0;

    if b == 0.0 {
        // Usa ef_lepton para o cálculo do momento de Fermi padrão
        let kf2 = ef_lepton.powi(2) - m.powi(2);
        if kf2 > 0.0 {
            let kf = kf2.sqrt();
            let dens_val = kf.powi(3) / (3.0 * PI2);
            let m_safe = m.abs().max(1e-15);
            let rhos_val =
                (m / (2.0 * PI2)) * (ef_lepton * kf - m.powi(2) * ((kf + ef_lepton) / m_safe).ln());

            engine.f_l[idx].push(kf);
            engine.n_l[idx] = 1;

            return (rhos_val, dens_val);
        }
        return (0.0, 0.0);
    }

    // Usa ef_lepton para estimar o nível máximo de Landau ocupável
    let nu_max_approx = (ef_lepton.powi(2) - m.powi(2)) / (2.0 * q * b);
    let nu_max = if nu_max_approx > 0.0 {
        (nu_max_approx.floor() as usize + 1).min(engine.max_landau_limit)
    } else {
        0
    };

    engine.f_l[idx].reserve(nu_max);

    for nu in 0..nu_max {
        let m_landau_2 = m.powi(2) + 2.0 * q * b * nu as f64;
        // O corte do preenchimento de cada nível depende do potencial efetivo ef_lepton
        let kf2 = ef_lepton.powi(2) - m_landau_2;

        if kf2 <= 0.0 {
            break;
        }

        let kf = kf2.sqrt();
        let g = if nu == 0 { 1.0 } else { 2.0 };

        engine.f_l[idx].push(kf);

        let factor = (g * q * b) / (2.0 * PI2);
        let m_safe = m_landau_2.sqrt().max(1e-15);

        // Integração termodinâmica atualizada com ef_lepton
        rhos += factor * m * ((kf + ef_lepton) / m_safe).ln();
        dens += factor * kf;

        n_occupied += 1;
    }

    engine.n_l[idx] = n_occupied;

    (rhos, dens)
}

fn compute(
    engine: &DarkPhotonsMatter,
    mue: f64,
    vsigma: f64,
    vomega: f64,
    vrho: f64,
) -> (f64, f64) {
    // 1. Energia dos mésons (Potenciais de campo)
    // Inclui termos de massa e auto-interações (rb, rc para sigma e rxi para omega)
    let ener_mesons = (vsigma / engine.model.gs).powi(2) / 2.0
        + (vomega / engine.model.gv).powi(2) / 2.0
        + (vrho / engine.model.gr).powi(2) / 2.0
        + engine.model.rb * vsigma.powi(3) / 3.0
        + engine.model.rc * vsigma.powi(4) / 4.0
        + engine.model.rxi * vomega.powi(4) / 4.0
        + engine.model.lambda_v * vomega.powi(2) * vrho.powi(2);

    let mut enerbar = 0.0;

    // --- LOOP DE BARIÕES (0:n, 1:p, 2:L0, 3:S-, 4:S0, 5:S+, 6:X-, 7:X0) ---
    for i in 0..8 {
        let ef = engine.ef_b[i];
        if ef <= 0.0 {
            continue;
        }

        // Se a partícula for neutra OU B=0, usa a fórmula contínua!
        if engine.charges_b[i] == 0.0 || engine.b == 0.0 {
            // --- Partículas Neutras (n, L0, S0, X0) ---
            // O AMM desdobra a partícula em 2 estados de spin (Up e Down)
            for &kf in [
                engine.kf_b_up[i].first().copied(),
                engine.kf_b_down[i].first().copied(),
            ]
            .iter()
            .flatten()
            {
                if kf > 0.0 {
                    let m_spin = (ef.powi(2) - kf.powi(2)).max(0.0).sqrt();
                    let m_safe = m_spin.max(1e-15);

                    // Fórmula para um único estado de spin (g=1), fator 1/4pi^2
                    enerbar += (1.0 / (4.0 * PI2))
                        * (ef.powi(3) * kf / 2.0
                            - (m_spin / 4.0)
                                * (m_spin * kf * ef
                                    + m_spin.powi(3) * ((kf + ef) / m_safe.abs()).ln()));
                }
            }
        } else {
            // --- Partículas Carregadas (p, S-, S+, X-) ---
            // Soma sobre os níveis de Landau (nu) para ambos os spins
            let qb = engine.charges_b[i].abs() * engine.qe * engine.b;
            let factor = qb / (4.0 * PI2);

            // Contribuição Spin Up
            for nu in 0..engine.n_b_up[i] {
                let kf = engine.kf_b_up[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerbar += factor * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }

            // Contribuição Spin Down
            for nu in 0..engine.n_b_down[i] {
                let kf = engine.kf_b_down[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerbar += factor * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }
        }
    }

    // --- LOOP DE LÉPTONS (0:e-, 1:mu-) ---
    let mut enerlep = 0.0;
    for i in 0..2 {
        let ef = engine.ef_l[i];
        if ef <= 0.0 {
            continue;
        }

        if engine.b == 0.0 {
            // Fórmula isotrópica para léptons se B=0 (com fator 2 para spin-up e spin-down)
            let kf = engine.f_l[i].first().copied().unwrap_or(0.0);
            if kf > 0.0 {
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                enerlep += 2.0
                    * (1.0 / (4.0 * PI2))
                    * (ef.powi(3) * kf / 2.0
                        - (m_spin / 4.0)
                            * (m_spin * kf * ef
                                + m_spin.powi(3) * ((kf + ef) / m_spin.abs()).ln()));
            }
        } else {
            // Léptons sob Efeito de Landau (B > 0)
            let qb = engine.qe * engine.b;

            for nu in 0..engine.n_l[i] {
                let kf = engine.f_l[i][nu];
                let m_spin = (ef.powi(2) - kf.powi(2)).sqrt();
                let g = if nu == 0 { 1.0 } else { 2.0 }; // Degenerescência de Landau para Dirac

                enerlep += (g * qb / (4.0 * PI2))
                    * (ef * kf + m_spin.powi(2) * ((kf + ef) / m_spin.abs()).ln());
            }
        }
    }

    let ener_dark_vector = engine.dark_vector_energy_density();
    let ener = ener_mesons + enerbar + enerlep + engine.ener_chi_kin + ener_dark_vector;

    // Pressão via relação termodinâmica: P = sum(mu_i * n_i) - epsilon
    let mut press_sum = 0.0;
    for i in 0..8 {
        press_sum += engine.mu_b[i] * engine.nb[i];
    }
    for i in 0..2 {
        press_sum += mue * engine.nl[i]; // mu_e = mu_mu = mue
    }
    press_sum += engine.mu_chi * engine.n_chi;

    let press = press_sum - ener;

    (ener, press)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::model::{FSU2, GM1, GM3};
    use crate::core::physics::HadronsMatter;
    use crate::core::solver::{EngineMode, Solver};

    const TEST_MUN: f64 = 1.15;
    const TEST_Y_CHI: f64 = 0.02;
    const TEST_M_CHI: f64 = 0.8;
    const TEST_M_X: f64 = 0.4;

    fn assert_close(left: f64, right: f64, rel: f64, abs: f64) {
        let error = (left - right).abs();
        let scale = left.abs().max(right.abs());
        assert!(
            error <= abs.max(rel * scale),
            "left={left:.16e}, right={right:.16e}, error={error:.3e}"
        );
    }

    fn solve(engine: &mut DarkPhotonsMatter, mun: f64) -> ([f64; 5], [f64; RESULTS_SIZE]) {
        engine
            .solve_point(mun, &[0.0; 5])
            .unwrap_or_else(|| panic!("dark-sector point failed to converge at mun={mun}"))
    }

    fn visible_state(engine: &DarkPhotonsMatter) -> Vec<f64> {
        engine
            .nb
            .iter()
            .chain(engine.nl.iter())
            .copied()
            .chain([engine.mue])
            .collect()
    }

    fn assert_same_visible_state(left: &DarkPhotonsMatter, right: &DarkPhotonsMatter) {
        for (a, b) in visible_state(left).iter().zip(visible_state(right).iter()) {
            assert_close(*a, *b, 2e-8, 2e-11);
        }
    }

    #[test]
    fn dirac_gas_density_and_thermodynamic_identity() {
        for ratio in [0.0, 1e-5, 9.9e-3, 1.01e-2, 0.1, 1.0, 10.0, 100.0] {
            let mass = 0.73;
            let kf = ratio * mass;
            let density = dark_fermion_number_density_from_kf(kf);
            let energy = dark_fermion_energy_density(kf, mass);
            let pressure = dark_fermion_pressure(kf, mass);
            let ef = kf.hypot(mass);

            assert_close(dark_fermion_kf_from_density(density), kf, 2e-13, 1e-18);
            assert_close(energy + pressure, density * ef, 5e-10, 1e-27);
            assert!(energy >= 0.0 && pressure >= 0.0);
        }
    }

    #[test]
    fn builders_validate_and_convert_dark_parameters() {
        let engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi_mev(500.0)
            .with_m_x_mev(25.0)
            .with_g_d(0.4)
            .with_epsilon(1e-3)
            .with_y_chi(0.03);

        assert_close(engine.m_chi, 500.0 / M_NUCLEON, 0.0, f64::EPSILON);
        assert_close(engine.m_x, 25.0 / M_NUCLEON, 0.0, f64::EPSILON);
        assert_eq!(engine.g_d, 0.4);
        assert_eq!(engine.epsilon, 1e-3);
        assert_eq!(engine.y_chi, 0.03);

        assert!(
            std::panic::catch_unwind(|| { DarkPhotonsMatter::new(GM1, 0.0).with_m_x(f64::NAN) })
                .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                DarkPhotonsMatter::new(GM1, 0.0).with_g_d(f64::INFINITY)
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| { DarkPhotonsMatter::new(GM1, 0.0).with_epsilon(1.0) })
                .is_err()
        );
    }

    #[test]
    #[should_panic(expected = "dark fermion mass")]
    fn builder_rejects_nonpositive_dark_fermion_mass() {
        let _ = DarkPhotonsMatter::new(GM1, 0.0).with_m_chi(0.0);
    }

    #[test]
    #[should_panic(expected = "dark number fraction")]
    fn builder_rejects_negative_dark_fraction() {
        let _ = DarkPhotonsMatter::new(GM1, 0.0).with_y_chi(-1e-3);
    }

    #[test]
    fn zero_fraction_recovers_visible_dark_engine() {
        for model in [GM1, GM3, FSU2] {
            let mut baseline = HadronsMatter::new(model, 0.0);
            let mut zero_fraction = DarkPhotonsMatter::new(model, 0.0)
                .with_m_chi(TEST_M_CHI)
                .with_m_x(TEST_M_X)
                .with_g_d(0.8)
                .with_epsilon(0.2)
                .with_y_chi(0.0);

            let mut visible_guess = [0.0; 4];
            let mut dark_guess = [0.0; 5];
            for mun in [1.12, 1.15, 1.18] {
                let (next_visible, baseline_result) = baseline
                    .solve_point(mun, &visible_guess)
                    .expect("hadronic reference point must converge");
                let (next_dark, result) = zero_fraction
                    .solve_point(mun, &dark_guess)
                    .expect("zero-fraction dark point must converge");
                visible_guess = next_visible;
                dark_guess = next_dark;

                for i in 0..=20 {
                    assert_close(result[i], baseline_result[i], 2e-8, 2e-9);
                }
                assert_eq!(zero_fraction.n_chi, 0.0);
                assert_close(zero_fraction.v_x0, 0.0, 0.0, 1e-15);
                assert_eq!(zero_fraction.ener_chi_kin, 0.0);
                assert_eq!(zero_fraction.press_chi_kin, 0.0);
            }
        }
    }

    #[test]
    fn free_dark_gas_changes_only_kinetic_thermodynamics() {
        let mut baseline = DarkPhotonsMatter::new(GM1, 0.0);
        let (_, baseline_result) = solve(&mut baseline, TEST_MUN);

        let mut dark = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_y_chi(TEST_Y_CHI);
        let (_, result) = solve(&mut dark, TEST_MUN);

        assert_same_visible_state(&baseline, &dark);
        assert_close(dark.v_x0, 0.0, 0.0, 1e-14);
        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        assert_close(
            result[1] - baseline_result[1],
            dark.ener_chi_kin * conversion,
            2e-8,
            2e-8,
        );
        assert_close(
            result[2] - baseline_result[2],
            dark.press_chi_kin * conversion,
            2e-8,
            2e-8,
        );
    }

    #[test]
    fn dark_self_interaction_matches_analytic_proca_solution() {
        let mut free = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_y_chi(TEST_Y_CHI);
        let (_, free_result) = solve(&mut free, TEST_MUN);

        let mut interacting = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let (_, result) = solve(&mut interacting, TEST_MUN);

        assert_same_visible_state(&free, &interacting);
        let expected_x0 = interacting.g_d * interacting.n_chi / interacting.m_x.powi(2);
        assert_close(interacting.v_x0, expected_x0, 2e-9, 2e-11);

        let vector_energy = interacting.dark_vector_energy_density();
        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        assert_close(
            result[1] - free_result[1],
            vector_energy * conversion,
            2e-8,
            2e-8,
        );
        assert_close(
            result[2] - free_result[2],
            vector_energy * conversion,
            2e-8,
            2e-8,
        );

        let vector_gibbs = interacting.g_d * interacting.n_chi * interacting.v_x0;
        assert_close(vector_gibbs, 2.0 * vector_energy, 2e-9, 2e-12);
        assert_close(
            interacting.mu_chi * interacting.n_chi - interacting.ener_chi_kin - vector_energy,
            interacting.press_chi_kin + vector_energy,
            2e-9,
            2e-12,
        );
    }

    #[test]
    fn kinetic_mixing_satisfies_neutrality_proca_and_sign_convention() {
        let epsilon = 0.08;
        let mut no_portal = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let _ = solve(&mut no_portal, TEST_MUN);

        let mut portal = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(epsilon)
            .with_y_chi(TEST_Y_CHI);
        let (solution, _) = solve(&mut portal, TEST_MUN);

        let charge = portal.charge_neutrality();
        assert!(charge.abs() < 1e-10, "charge residual={charge:.3e}");
        assert!(
            portal.dark_photon_residual(charge).abs() < 1e-10,
            "Proca residual={:.3e}",
            portal.dark_photon_residual(charge)
        );
        let norm = DarkPhotonsMatter::kinetic_mixing_norm(epsilon);
        let expected_x0 = portal.g_d * portal.n_chi / (portal.m_x.powi(2) * norm);
        assert_close(portal.v_x0, expected_x0, 2e-8, 2e-10);
        assert_close(solution[4], portal.v_x0, 0.0, 1e-15);
        assert_close(portal.n_chi, portal.y_chi * portal.nbt, 1e-14, 1e-15);
        assert_close(
            portal.mu_chi,
            portal.ef_chi + portal.g_d * portal.v_x0 / norm,
            1e-14,
            1e-14,
        );
        assert_close(
            portal.g_d * portal.n_chi * portal.v_x0 / norm,
            2.0 * portal.dark_vector_energy_density(),
            2e-8,
            2e-11,
        );

        let delta_x = epsilon * portal.qe * portal.v_x0 / norm;
        assert_close(portal.ef_l[0], portal.mue + delta_x, 2e-12, 2e-13);
        assert_close(portal.mue + delta_x, no_portal.mue, 3e-5, 5e-8);
        for (left, right) in portal
            .nb
            .iter()
            .chain(portal.nl.iter())
            .zip(no_portal.nb.iter().chain(no_portal.nl.iter()))
        {
            assert_close(*left, *right, 3e-5, 5e-9);
        }
    }

    #[test]
    fn total_pressure_obeys_gibbs_relation() {
        let mut baseline = HadronsMatter::new(GM1, 0.0);
        let (_, baseline_result) = baseline
            .solve_point(TEST_MUN, &[0.0; 4])
            .expect("hadronic reference point must converge");

        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.04)
            .with_y_chi(TEST_Y_CHI);
        let (solution, result) = solve(&mut engine, TEST_MUN);
        let (energy, pressure) =
            compute(&engine, solution[0], solution[1], solution[2], solution[3]);
        let gibbs = engine
            .mu_b
            .iter()
            .zip(engine.nb.iter())
            .map(|(mu, n)| mu * n)
            .sum::<f64>()
            + engine.mue * engine.nl.iter().sum::<f64>()
            + engine.mu_chi * engine.n_chi
            - energy;
        assert_close(pressure, gibbs, 1e-13, 1e-14);

        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        let expected_dark_pressure =
            (engine.press_chi_kin + engine.dark_vector_energy_density()) * conversion;
        assert_close(
            result[2] - baseline_result[2],
            expected_dark_pressure,
            5e-6,
            5e-7,
        );
    }

    #[test]
    fn proca_residual_includes_off_shell_visible_charge_source() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.2)
            .with_y_chi(TEST_Y_CHI);
        engine.n_chi = 3e-4;
        engine.v_x0 = 0.07;
        let n_q = -2.5e-3;
        let norm = DarkPhotonsMatter::kinetic_mixing_norm(engine.epsilon);
        let expected = engine.m_x.powi(2) * engine.v_x0
            - (engine.g_d * engine.n_chi + engine.epsilon * engine.qe * n_q) / norm;

        assert_close(engine.dark_photon_residual(n_q), expected, 1e-14, 1e-14);
        assert_ne!(
            engine.dark_photon_residual(n_q),
            engine.dark_photon_residual(0.0)
        );
    }

    #[test]
    fn dark_sector_has_vacuum_limit() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.05)
            .with_y_chi(TEST_Y_CHI);
        engine.mun = 0.1;
        let residual = engine.funcv(&[0.0; 5]);

        assert!(residual.iter().all(|value| value.abs() < 1e-15));
        assert_eq!(engine.nbt, 0.0);
        assert_eq!(engine.n_chi, 0.0);
        assert_eq!(engine.v_x0, 0.0);
        assert_eq!(engine.ener_chi_kin, 0.0);
        assert_eq!(engine.press_chi_kin, 0.0);

        let mut previous_energy = f64::INFINITY;
        let mut previous_pressure = f64::INFINITY;
        for n_b in [1e-4, 1e-6, 1e-8, 1e-10] {
            engine.n_chi = engine.y_chi * n_b;
            engine.v_x0 = engine.g_d * engine.n_chi
                / (engine.m_x.powi(2) * DarkPhotonsMatter::kinetic_mixing_norm(engine.epsilon));
            engine.update_dark_fermion_state();
            let energy = engine.ener_chi_kin + engine.dark_vector_energy_density();
            let pressure = engine.press_chi_kin + engine.dark_vector_energy_density();
            assert_close(engine.n_chi / n_b, engine.y_chi, 2e-14, 1e-15);
            assert!(energy < previous_energy);
            assert!(pressure < previous_pressure);
            previous_energy = energy;
            previous_pressure = pressure;
        }
        assert!(previous_energy < 1e-10);
        assert!(previous_pressure < 1e-15);
    }

    #[test]
    fn linear_gm1_vector_residual_matches_documented_normalization() {
        // This locks the linear GM1 normalization only. FSU2 is deliberately
        // documented as a known legacy inconsistency until both visible-matter
        // paths are revised together.
        let mut gm1 = DarkPhotonsMatter::new(GM1, 0.0);
        gm1.nb[0] = 0.01;
        let vomega = 0.12;
        let derivative = vomega / GM1.gv.powi(2);
        assert_close(
            gm1.equation_omega(vomega, 0.0),
            GM1.gv.powi(2) * (gm1.nb[0] - derivative),
            1e-14,
            1e-14,
        );
        assert_ne!(FSU2.rxi, 0.0);
        assert_ne!(FSU2.lambda_v, 0.0);
    }

    #[test]
    fn eos_columns_use_documented_physical_units() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let (solution, result) = solve(&mut engine, TEST_MUN);
        let density_factor = (M_NUCLEON / HBAR_C).powi(3);
        let energy_factor = M_NUCLEON * density_factor;

        assert_eq!(result.len(), 34);
        assert_close(result[3], engine.nl[0] * density_factor, 1e-14, 1e-14);
        assert_close(result[5], engine.nb[0] * density_factor, 1e-14, 1e-14);
        assert_close(result[13], solution[1] * M_NUCLEON, 1e-14, 1e-14);
        assert_close(result[16], engine.m_eff[0], 1e-14, 1e-14);
        assert_close(result[21], engine.n_chi * density_factor, 1e-14, 1e-14);
        assert_close(result[23], engine.m_chi * M_NUCLEON, 1e-14, 1e-14);
        assert_close(result[27], engine.v_x0 * M_NUCLEON, 1e-14, 1e-14);
        assert_close(
            result[30],
            engine.ener_chi_kin * energy_factor,
            1e-14,
            1e-14,
        );
        assert_close(result[32], result[33], 0.0, 0.0);
    }

    #[test]
    fn previous_dark_solution_is_a_valid_five_dimensional_guess() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(0.01)
            .with_g_d(0.8)
            .with_epsilon(1e-3)
            .with_y_chi(0.05);
        let (first, _) = solve(&mut engine, 1.12);
        assert_ne!(first[4], 0.0);

        let (second, _) = engine
            .solve_point(1.14, &first)
            .expect("next EOS point must accept the previous five-dimensional solution");
        let residual = engine.funcv(&second);
        assert!(
            DarkStateVector::from_column_slice(&residual).norm() < 1e-10,
            "residual at continued point is {:?}",
            residual
        );
    }

    #[test]
    fn generic_solver_carries_the_dark_field_between_eos_points() {
        let engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_limits(1.12, 1.14)
            .with_points(3)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let mut solver = Solver::new(EngineMode::DarkPhotons(engine));
        let rows = solver.solve();

        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|row| row[21] > 0.0));
        assert!(rows.iter().all(|row| row[27] > 0.0));
    }
}
