// solver/physics.rs
#![allow(unused)]

use crate::core::constants::{
    AMML0, AMMN, AMMP, AMMS0, AMMSM, AMMSP, AMMX0, AMMXM, BCE, BCE_G, BDD_ALPHAA, BDD_BETAA,
    HBAR_C, M_NUCLEON, MAX_LANDAU_LIMIT, MB, ML, N0, QE, RESULTS_SIZE,
};
use crate::core::magnetic::{FieldProfile, magnetic_stress};
use crate::core::model::ModelParams;
use nalgebra::{Matrix4, Vector4};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NlemModel {
    Maxwell,     // Eletromagnetismo Clássico (Linear)
    Modmax(f64), // Eletrodinâmica ModMax (recebe parâmetro csi)
    Log(f64),    // Eletrodinâmica Logarítmica (recebe parâmetro csi)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MagneticTopology {
    Isotropic,   // campo emaranhado: P_mag = (P_par + 2 P_perp)/3 (Maxwell: eps/3)
    Anisotropic, // P_mag = P_perp (Maxwell: eps)
}

impl NlemModel {
    /// Recebe o campo magnético original (bg, em Gauss) e retorna o campo
    /// EFETIVO dos níveis de Landau. Para `Log`, csi (xi) também é em Gauss.
    pub fn effective_bg(&self, bg: f64) -> f64 {
        match self {
            NlemModel::Maxwell => bg,

            NlemModel::Modmax(csi) => {
                // Fórmula: bg * exp(-csi)
                bg * (-csi).exp()
            }

            NlemModel::Log(csi) => bg * (1.0 + bg.powi(2) / (2.0 * csi.powi(2))),
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

    /// Perfil do campo local (ver `core::magnetic`). Padrão: `Constant`.
    pub field_profile: FieldProfile,
    /// Campo local (Gauss) usado no último ponto resolvido.
    pub local_field_g: f64,
    /// n_B/n0 do último ponto aceito: chute da iteração do perfil BDD.
    last_nb_over_n0: f64,
}

impl HadronsMatter {
    // constante estática para acoplamento sigma
    const X_SIGMA: [f64; 8] = [1.0, 1.0, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];

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

        let kf_b_up = std::array::from_fn(|_| Vec::new());
        let kf_b_down = std::array::from_fn(|_| Vec::new());
        let f_l = std::array::from_fn(|_| Vec::new());

        let ef_l = [0.0; 2];
        let n_l = [0; 2];

        let isospin_factor = [-0.5, 0.5, 0.0, -1.0, 0.0, 1.0, -0.5, 0.5];

        HadronsMatter {
            model,
            nlem: NlemModel::Maxwell,
            topology: MagneticTopology::Anisotropic,
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

            isospin_factor: isospin_factor,
            eos_output: None,

            field_profile: FieldProfile::Constant,
            local_field_g: bg,
            last_nb_over_n0: 0.0,
        }
    }

    /// Define o perfil do campo magnético local. Com `Bdd` e `Dexheimer2017`
    /// o campo local entra nos níveis de Landau e na energia magnética; o
    /// argumento `bg` de `new()` deixa de ser usado.
    pub fn with_field_profile(mut self, profile: FieldProfile) -> Self {
        self.field_profile = profile;
        self
    }

    /// Campo dos níveis de Landau a partir do campo local em Gauss.
    fn set_landau_field(&mut self, b_gauss: f64) {
        self.local_field_g = b_gauss;
        self.b = self.nlem.effective_bg(b_gauss) / BCE_G * BCE;
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

    pub fn with_eos_output<P: Into<String>>(mut self, path: P) -> Self {
        self.eos_output = Some(path.into());
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
        crate::core::particles::calculate_all_densities(self, vomega, vrho);

        let fsigma = self.equation_sigma(vsigma);
        let fomega = self.equation_omega(vomega, vrho);
        let frho = self.equation_rho(vrho, vomega);
        let charge_neutral = self.charge_neutrality();

        [fsigma, fomega, frho, charge_neutral]
    }

    fn equation_sigma(&self, vsigma: f64) -> f64 {
        let gs2 = self.model.gs.powi(2);
        gs2 * (self.rhosb - self.model.rb * vsigma.powi(2) - self.model.rc * vsigma.powi(3))
            - vsigma
    }

    // Equações de Campo Vetorizadas para suportar as partículas com total precisão
    fn equation_omega(&self, vomega: f64, vrho: f64) -> f64 {
        let mut sum_baryon = 0.0;
        for i in 0..8 {
            sum_baryon += self.nb[i] * self.xv_v[i];
        }
        self.model.gv.powi(2)
            * (sum_baryon
                - self.model.rxi * vomega.powi(3)
                - 2.0 * self.model.lambda_v * vomega * vrho.powi(2))
            - vomega
    }

    fn equation_rho(&self, vrho: f64, vomega: f64) -> f64 {
        let mut sum_source = 0.0;
        for i in 0..8 {
            // A fonte para o rho é baseada no negativo do isospin
            sum_source += self.isospin_factor[i] * self.nb[i] * self.xv_r[i];
        }
        self.model.gr.powi(2)
            * (sum_source - 2.0 * self.model.lambda_v * vrho * vomega.powi(2))
            - vrho
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

    /// Resolve para um dado mun e chute inicial, retorna solução e resultado.
    pub fn solve_point(
        &mut self,
        mun: f64,
        initial_x: &[f64],
    ) -> Option<([f64; 4], [f64; RESULTS_SIZE])> {
        let profile = self.field_profile;
        let mu_b_mev = mun * self.m_nuc;
        let solution = if !profile.depends_on_density() {
            let field = profile.local_field_g(0.0, mu_b_mev);
            if let Some(b_gauss) = field {
                self.set_landau_field(b_gauss);
            }
            self.solve_point_with_field(mun, initial_x, field)
        } else {
            // B depende de n_B, que depende de B: iteração de ponto fixo,
            // partindo da densidade do último ponto aceito.
            let mut nb_over_n0 = self.last_nb_over_n0;
            let mut x = [initial_x[0], initial_x[1], initial_x[2], initial_x[3]];
            let mut converged = None;
            for _ in 0..200 {
                let b_gauss = profile.local_field_g(nb_over_n0, mu_b_mev)?;
                self.set_landau_field(b_gauss);
                let (x_new, result) = self.solve_point_with_field(mun, &x, Some(b_gauss))?;
                let change = (result[0] - nb_over_n0).abs();
                x = x_new;
                nb_over_n0 = result[0];
                if change <= 1e-10 * nb_over_n0.max(1e-6) {
                    converged = Some((x_new, result));
                    break;
                }
            }
            converged
        };
        if let Some((_, result)) = &solution {
            self.last_nb_over_n0 = result[0];
        }
        solution
    }

    /// Resolve um ponto com o campo dos níveis de Landau já definido em
    /// `self.b`. `field_g` é o campo local usado na energia magnética; `None`
    /// usa o perfil legado do modo `Constant`.
    pub fn solve_point_with_field(
        &mut self,
        mun: f64,
        initial_x: &[f64],
        field_g: Option<f64>,
    ) -> Option<([f64; 4], [f64; RESULTS_SIZE])> {
        self.mun = mun;

        let mut x = Vector4::from_column_slice(initial_x);
        let tolerance = 1e-10;
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
                let h = 1e-8 * (x[i].abs() + 1e-2);
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
        let (ener, press) = crate::core::eos::compute(self, mue, vsigma, vomega, vrho);

        let nb_total = self.nb.iter().sum::<f64>();
        let nbtd = nb_total * (self.m_nuc / HBAR_C).powi(3);

        let factor_mev_fm3 = self.m_nuc * (self.m_nuc / HBAR_C).powi(3);
        let ener_conv = ener * factor_mev_fm3;
        let press_conv = press * factor_mev_fm3;

        // Campo local da energia magnética. No modo `Constant` (legado) é o
        // perfil BDD com B_surf = 1e15 G e B0 = bg, nulo se bg = 0.
        let b_local_g = field_g.unwrap_or_else(|| {
            if self.bg == 0.0 {
                0.0
            } else {
                crate::core::magnetic::bdd_field_g(
                    crate::core::magnetic::B_SURFACE_G,
                    self.bg,
                    BDD_BETAA,
                    BDD_ALPHAA,
                    nbtd / N0,
                )
            }
        });
        let stress = magnetic_stress(self.nlem, b_local_g);
        let ebsd = stress.energy;
        let pmag_effective = stress.pressure(self.topology);

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + pmag_effective;

        if ener_final >= 0.0 && press_final >= 0.0 {
            let fermion_mu_density = self
                .mu_b
                .iter()
                .zip(self.nb.iter())
                .map(|(mu, n)| mu * n)
                .sum::<f64>()
                + mue * self.nl.iter().sum::<f64>();
            let mu_total_per_baryon = if self.nbt > 0.0 {
                fermion_mu_density / self.nbt
            } else {
                0.0
            };

            let density_factor = (self.m_nuc / HBAR_C).powi(3);
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
            Some((x_final, result))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::model::{FSU2, GM1};
    use nalgebra::{Matrix3, Vector3};

    /// Solves the three meson equations at fixed (mu_n, mu_e) with B = 0,
    /// i.e. without imposing charge neutrality. Returns the fields and
    /// (n_B, eps, P) in units of M_N^3 and M_N^4.
    fn solve_fields(
        engine: &mut HadronsMatter,
        mun: f64,
        mue: f64,
        guess: [f64; 3],
    ) -> Option<([f64; 3], f64, f64, f64)> {
        engine.mun = mun;
        let residual = |engine: &mut HadronsMatter, y: &Vector3<f64>| {
            let f = engine.funcv(&[mue, y[0], y[1], y[2]]);
            Vector3::new(f[0], f[1], f[2])
        };
        let mut y = Vector3::from(guess);
        for _ in 0..200 {
            let f = residual(engine, &y);
            if f.norm() < 1e-13 {
                let _ = engine.funcv(&[mue, y[0], y[1], y[2]]);
                let (ener, press) = crate::core::eos::compute(engine, mue, y[0], y[1], y[2]);
                return Some(([y[0], y[1], y[2]], engine.nbt, ener, press));
            }
            let mut jac = Matrix3::zeros();
            for i in 0..3 {
                let h = 1e-7 * (y[i].abs() + 1e-3);
                let mut yh = y;
                yh[i] += h;
                jac.set_column(i, &((residual(engine, &yh) - f) / h));
            }
            let step = jac.lu().solve(&(-f))?;
            let mut alpha = 1.0;
            while alpha > 1e-6 && residual(engine, &(y + alpha * step)).norm() >= f.norm() {
                alpha *= 0.5;
            }
            y += alpha * step;
        }
        None
    }

    /// Locates the symmetric-matter saturation point as the zero of P(mu_n)
    /// on the dense branch. Returns (n0 [fm^-3], E/A [MeV], M*/M).
    fn saturation_point(model: ModelParams) -> (f64, f64, f64) {
        let mut engine = HadronsMatter::new(model, 0.0);
        // Walk down the dense branch until the pressure changes sign.
        let (mut mu_hi, mut g_hi) = (1.2, [0.0; 3]);
        let mut mu = 1.2;
        let mut guess = [0.4, 0.4, 0.0];
        let mut mu_lo = loop {
            let (fields, nb, _, press) =
                solve_fields(&mut engine, mu, 0.0, guess).expect("dense branch must converge");
            assert!(nb > 1e-5, "continuation fell onto the vacuum branch");
            guess = fields;
            if press < 0.0 {
                break mu;
            }
            (mu_hi, g_hi) = (mu, fields);
            mu -= if mu > 1.0 { 5e-3 } else { 5e-4 };
        };
        for _ in 0..50 {
            let mid = 0.5 * (mu_hi + mu_lo);
            let (fields, _, _, press) =
                solve_fields(&mut engine, mid, 0.0, g_hi).expect("bisection point must converge");
            if press > 0.0 {
                (mu_hi, g_hi) = (mid, fields);
            } else {
                mu_lo = mid;
            }
        }
        let (fields, nb, _, _) = solve_fields(&mut engine, mu_hi, 0.0, g_hi).unwrap();
        let n0 = nb * (M_NUCLEON / HBAR_C).powi(3);
        // At P = 0 the Gibbs relation gives eps/n_B = mu_n; the binding energy
        // is measured from the mean bare nucleon mass (m_n != m_p here).
        let m_avg = 0.5 * (MB[0] + MB[1]);
        (n0, (mu_hi - m_avg) * M_NUCLEON, 1.0 - fields[0])
    }

    #[test]
    fn gm1_saturation_matches_glendenning_moszkowski() {
        let (n0, ea, mstar) = saturation_point(GM1);
        assert!((n0 - 0.153).abs() < 0.002, "n0 = {n0}");
        assert!((ea + 16.3).abs() < 0.1, "E/A = {ea}");
        assert!((mstar - 0.70).abs() < 0.005, "M*/M = {mstar}");
    }

    #[test]
    fn fsu2_saturation_matches_chen_piekarewicz() {
        // Phys. Rev. C 90, 044305 (2014): n0 = 0.1505 fm^-3,
        // E/A = -16.28 MeV, M*/M = 0.593.
        let (n0, ea, mstar) = saturation_point(FSU2);
        assert!((n0 - 0.1505).abs() < 0.002, "n0 = {n0}");
        assert!((ea + 16.28).abs() < 0.1, "E/A = {ea}");
        assert!((mstar - 0.593).abs() < 0.005, "M*/M = {mstar}");
    }

    #[test]
    fn fsu2_pressure_is_consistent_with_field_equations() {
        // At fixed mu_e, dP/dmu_n = n_B only if the meson field equations are
        // the stationarity conditions of the energy functional in compute().
        // A nonzero mu_e makes the matter asymmetric, so rho and the
        // omega-rho coupling Lambda_v are exercised.
        let mut engine = HadronsMatter::new(FSU2, 0.0);
        let (mun, mue, dmu) = (1.08, 0.12, 1e-5);
        let mut guess = [0.4, 0.4, -0.05];
        for mu in [1.20, 1.15, 1.10, mun] {
            guess = solve_fields(&mut engine, mu, mue, guess).expect("continuation").0;
        }
        let (fields, nb, _, _) = solve_fields(&mut engine, mun, mue, guess).unwrap();
        assert!(fields[2].abs() > 1e-3, "rho field must be active");
        let (_, _, _, p_plus) = solve_fields(&mut engine, mun + dmu, mue, fields).unwrap();
        let (_, _, _, p_minus) = solve_fields(&mut engine, mun - dmu, mue, fields).unwrap();
        let dp_dmu = (p_plus - p_minus) / (2.0 * dmu);
        assert!(
            ((dp_dmu - nb) / nb).abs() < 1e-5,
            "dP/dmu_n = {dp_dmu}, n_B = {nb}"
        );
    }
}
