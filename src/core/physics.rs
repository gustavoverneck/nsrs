// solver/physics.rs
#![allow(unused)]

use crate::core::constants::{
    BCE, BCE_G, BDD_ALPHAA, BDD_BETAA, KAPPA_B, RNCM,
    HBAR_C, M_NUCLEON, MAX_LANDAU_LIMIT, MB, ML, N0, QE, RESULTS_SIZE,
};
use crate::core::magnetic::{FieldProfile, magnetic_stress};
use crate::core::model::ModelParams;
use nalgebra::{SMatrix, SVector};

/// Estado do Newton: (mu_e, g_s sigma, g_v omega, g_rho rho, X_0), todos
/// divididos por M_N. Sem setor escuro, X_0 = 0.
type StateVector = SVector<f64, 5>;
type Jacobian = SMatrix<f64, 5, 5>;

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
    /// Razão H/B do campo auxiliar, H = d eps_B / dB, para um campo
    /// magnético estático puro B (Gauss). Maxwell: 1; ModMax: e^{-gamma};
    /// Log: 1/(1 + B^2/(2 xi^2)), com xi em Gauss (Soleng 1995).
    ///
    /// As partículas acoplam ao potencial vetor, portanto os níveis de Landau
    /// usam B; a eletrodinâmica não linear entra apenas na energia e nas
    /// tensões do campo (ver `core::magnetic::magnetic_stress`).
    pub fn h_over_b(&self, b_gauss: f64) -> f64 {
        match *self {
            NlemModel::Maxwell => 1.0,
            NlemModel::Modmax(gamma) => (-gamma).exp(),
            NlemModel::Log(xi) => 1.0 / (1.0 + b_gauss * b_gauss / (2.0 * xi * xi)),
        }
    }

    /// Curvatura do vácuo f_vac = dH/dB = 4 pi d^2 eps_B/dB^2 (gaussiano), o
    /// termo do campo no critério de estabilidade s = f_vac - 4 pi d^2 P_m/dB^2.
    /// Maxwell: 1; ModMax: e^{-gamma}; Log: (1 - x)/(1 + x)^2, x = B^2/(2 xi^2),
    /// negativa para B > sqrt(2) xi.
    pub fn vacuum_curvature(&self, b_gauss: f64) -> f64 {
        match *self {
            NlemModel::Maxwell => 1.0,
            NlemModel::Modmax(gamma) => (-gamma).exp(),
            NlemModel::Log(xi) => {
                let x = b_gauss * b_gauss / (2.0 * xi * xi);
                (1.0 - x) / ((1.0 + x) * (1.0 + x))
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

    // --- Setor escuro opcional (férmion de Dirac chi + fóton escuro X) ---
    // Ver docs/PHYSICS.md, "Setor escuro fermiônico". Todos os parâmetros
    // são nulos por padrão; nesse caso X_0 = 0 e a matéria é a hadrônica.
    /// Mistura cinética (adimensional, |epsilon| < 1).
    pub epsilon: f64,
    /// Massa do fóton escuro / M_N.
    pub m_x: f64,
    /// Acoplamento U(1) escuro (adimensional).
    pub g_d: f64,
    /// Massa do férmion escuro / M_N.
    pub m_chi: f64,
    /// Fração imposta n_chi / n_B.
    pub y_chi: f64,
    /// Densidade escura (unidades de M_N^3).
    pub n_chi: f64,
    /// Momento de Fermi escuro / M_N.
    pub kf_chi: f64,
    /// Energia de Fermi escura / M_N.
    pub ef_chi: f64,
    /// Potencial químico escuro completo / M_N.
    pub mu_chi: f64,
    /// Energia cinética escura (M_N^4).
    pub ener_chi_kin: f64,
    /// Pressão cinética escura (M_N^4).
    pub press_chi_kin: f64,
    /// Potencial médio do fóton escuro / M_N.
    pub v_x0: f64,
    /// Se algum parâmetro escuro foi definido (controla as colunas 21-33).
    pub dark_enabled: bool,

    /// Inclui o octeto de hyperons (padrão). `false`: matéria npe(mu).
    pub include_hyperons: bool,

    /// Perfil do campo local (ver `core::magnetic`). Padrão: `Constant`.
    pub field_profile: FieldProfile,
    /// Campo local (Gauss) usado no último ponto resolvido.
    pub local_field_g: f64,
    /// n_B/n0 do último ponto aceito: chute da iteração do perfil BDD.
    last_nb_over_n0: f64,
    /// M B = B dP/dB|_mu (MeV/fm^3) do último ponto resolvido.
    pub magnetization_b: f64,
    /// Pressão total sem o termo de magnetização (MeV/fm^3) do último ponto:
    /// a pressão termodinâmica, monótona em mu, usada nos critérios de
    /// validade da varredura.
    pub stability_pressure: f64,
    /// Energia e tensões do campo (MeV/fm^3) somadas à EoS no último ponto;
    /// zero quando o campo não entra na EoS.
    pub field_stress: crate::core::magnetic::MagneticStress,
    /// Se a energia e as tensões do próprio campo entram na EoS da TOV.
    /// `None`: padrão do perfil (`FieldProfile::field_stress_in_eos`).
    field_stress_override: Option<bool>,
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

        // Sem momentos anômalos por padrão; ver `with_anomalous_moments`.
        let amm_b = [0.0; 8];

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

            epsilon: 0.0,
            m_x: 1.0,
            g_d: 0.0,
            m_chi: 1.0,
            y_chi: 0.0,
            n_chi: 0.0,
            kf_chi: 0.0,
            ef_chi: 0.0,
            mu_chi: 0.0,
            ener_chi_kin: 0.0,
            press_chi_kin: 0.0,
            v_x0: 0.0,
            dark_enabled: false,
            include_hyperons: true,

            field_profile: FieldProfile::Constant,
            local_field_g: bg,
            last_nb_over_n0: 0.0,
            magnetization_b: 0.0,
            stability_pressure: 0.0,
            field_stress: crate::core::magnetic::MagneticStress::default(),
            field_stress_override: None,
        }
    }

    /// Liga (ou desliga, o padrão) os momentos magnéticos anômalos dos bárions
    /// (`KAPPA_B`): acoplamento de Pauli a = s kappa mu_N B nos níveis de Landau dos
    /// carregados e no espectro anisotrópico dos neutros.
    pub fn with_anomalous_moments(mut self, include: bool) -> Self {
        self.amm_b = if include { KAPPA_B.map(|k| k * RNCM) } else { [0.0; 8] };
        self
    }

    /// Inclui (padrão) ou exclui os hyperons; sem eles a matéria é npe(mu),
    /// como nas parametrizações originais calibradas só com núcleons.
    pub fn with_hyperons(mut self, include: bool) -> Self {
        self.include_hyperons = include;
        self
    }

    // --- Builders do setor escuro ---

    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        Self::kinetic_mixing_norm(epsilon);
        self.epsilon = epsilon;
        self.dark_enabled = true;
        self
    }

    /// Massa do fóton escuro em unidades de M_N.
    pub fn with_m_x(mut self, m_x: f64) -> Self {
        assert!(
            m_x.is_finite() && m_x > 0.0,
            "dark photon mass m_x must be finite and positive"
        );
        self.m_x = m_x;
        self.dark_enabled = true;
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
        assert!(g_d.is_finite(), "dark coupling g_d must be finite");
        self.g_d = g_d;
        self.dark_enabled = true;
        self
    }

    /// Massa do férmion escuro em unidades de M_N.
    pub fn with_m_chi(mut self, m_chi: f64) -> Self {
        assert!(
            m_chi.is_finite() && m_chi > 0.0,
            "dark fermion mass m_chi must be finite and positive"
        );
        self.m_chi = m_chi;
        self.dark_enabled = true;
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
        self.dark_enabled = true;
        self
    }

    pub(crate) fn kinetic_mixing_norm(epsilon: f64) -> f64 {
        assert!(
            epsilon.is_finite() && epsilon.abs() < 1.0,
            "dark photon kinetic mixing epsilon must be finite and satisfy |epsilon| < 1"
        );
        (1.0 - epsilon.powi(2)).sqrt()
    }

    /// Deslocamento da energia de Fermi de um férmion visível de carga
    /// `charge_units` (em unidades de e) pelo fóton escuro. Nulo sem mistura.
    pub(crate) fn dark_shift_for_charge(&self, charge_units: f64) -> f64 {
        self.epsilon * charge_units * self.qe * self.v_x0 / Self::kinetic_mixing_norm(self.epsilon)
    }

    pub(crate) fn update_dark_fermion_state(&mut self) {
        use crate::core::darkphotons::{
            dark_fermion_energy_density, dark_fermion_kf_from_density, dark_fermion_pressure,
        };
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

    /// Equação de Proca para X_0 com fontes escura e visível (off-shell).
    pub(crate) fn dark_photon_residual(&self, charge_density: f64) -> f64 {
        self.m_x.powi(2) * self.v_x0
            - (self.g_d * self.n_chi + self.epsilon * self.qe * charge_density)
                / Self::kinetic_mixing_norm(self.epsilon)
    }

    pub(crate) fn dark_vector_energy_density(&self) -> f64 {
        0.5 * self.m_x.powi(2) * self.v_x0.powi(2)
    }

    /// Define o perfil do campo magnético local. Com `Bdd` e `Dexheimer2017`
    /// o campo local entra nos níveis de Landau e na energia magnética; o
    /// argumento `bg` de `new()` deixa de ser usado.
    pub fn with_field_profile(mut self, profile: FieldProfile) -> Self {
        self.field_profile = profile;
        self
    }

    /// Inclui (ou não) a energia e as tensões do próprio campo, B^2/8pi e
    /// suas generalizações NLEM, na EoS usada pela TOV. Sem esta chamada vale
    /// o padrão do perfil. A matéria (Landau e magnetização) usa o campo em
    /// qualquer caso.
    pub fn with_field_stress(mut self, include: bool) -> Self {
        self.field_stress_override = Some(include);
        self
    }

    /// Se a energia/tensão do campo entra na EoS deste motor.
    pub fn field_stress_in_eos(&self) -> bool {
        self.field_stress_override
            .unwrap_or_else(|| self.field_profile.field_stress_in_eos())
    }

    /// Campo dos níveis de Landau a partir do campo local em Gauss.
    fn set_landau_field(&mut self, b_gauss: f64) {
        self.local_field_g = b_gauss;
        self.b = b_gauss / BCE_G * BCE;
    }
    /// Define a topologia das linhas de campo magnético
    pub fn with_topology(mut self, top: MagneticTopology) -> Self {
        self.topology = top;
        self
    }

    /// Builder para acoplar o Eletromagnetismo Não-Linear. O campo dos
    /// níveis de Landau continua sendo B (acoplamento mínimo); a NLEM altera a
    /// energia e as tensões do campo.
    pub fn with_nlem(mut self, nlem: NlemModel) -> Self {
        self.nlem = nlem;
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

    /// Número máximo de níveis de Landau somados por espécie (padrão
    /// `MAX_LANDAU_LIMIT`). Em B baixo e densidade alta a soma precisa de mais
    /// níveis; o corte trunca a soma e distorce a magnetização.
    pub fn with_max_landau_limit(mut self, n: usize) -> Self {
        assert!(n > 0, "max_landau_limit must be positive");
        self.max_landau_limit = n;
        self
    }

    pub fn with_eos_output<P: Into<String>>(mut self, path: P) -> Self {
        self.eos_output = Some(path.into());
        self
    }

    // Mapeamento das variáveis (vindo do solver)
    /// (mu_e, vsigma, vomega, vrho, X_0). Aceita estados de 4 entradas
    /// (sem setor escuro), com X_0 = 0.
    pub fn mapping(&self, x: &[f64]) -> (f64, f64, f64, f64, f64) {
        let mue = x[0];
        let vsigma = x[1];
        let vomega = x[2];
        let vrho = x[3];
        let v_x0 = x.get(4).copied().unwrap_or(0.0);
        (mue, vsigma, vomega, vrho, v_x0)
    }

    // Função de resíduo (chamada pelo solver numérico)
    pub fn funcv(&mut self, x: &[f64]) -> [f64; 5] {
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
        crate::core::particles::calculate_all_densities(self, vomega, vrho);

        self.n_chi = self.y_chi * self.nbt;
        self.update_dark_fermion_state();

        let fsigma = self.equation_sigma(vsigma);
        let fomega = self.equation_omega(vomega, vrho);
        let frho = self.equation_rho(vrho, vomega);
        let charge_neutral = self.charge_neutrality();
        // Linha de Proca escalada por m_X^2, para que a tolerância global
        // controle o erro absoluto em X_0 mesmo para mediadores leves. Sem
        // setor escuro (m_X = 1, g_d = epsilon = 0) ela se reduz a X_0 = 0.
        let fdark = self.dark_photon_residual(charge_neutral) / self.m_x.powi(2);

        [fsigma, fomega, frho, charge_neutral, fdark]
    }

    pub(crate) fn equation_sigma(&self, vsigma: f64) -> f64 {
        let gs2 = self.model.gs.powi(2);
        gs2 * (self.rhosb - self.model.rb * vsigma.powi(2) - self.model.rc * vsigma.powi(3))
            - vsigma
    }

    // Equações de Campo Vetorizadas para suportar as partículas com total precisão
    pub(crate) fn equation_omega(&self, vomega: f64, vrho: f64) -> f64 {
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

    pub(crate) fn equation_rho(&self, vrho: f64, vomega: f64) -> f64 {
        let mut sum_source = 0.0;
        for i in 0..8 {
            // A fonte para o rho é baseada no negativo do isospin
            sum_source += self.isospin_factor[i] * self.nb[i] * self.xv_r[i];
        }
        self.model.gr.powi(2)
            * (sum_source - 2.0 * self.model.lambda_v * vrho * vomega.powi(2))
            - vrho
    }

    pub(crate) fn charge_neutrality(&self) -> f64 {
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
    ) -> Option<([f64; 5], [f64; RESULTS_SIZE])> {
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
            let (mue, vs, vw, vr, x0) = self.mapping(initial_x);
            let mut x = [mue, vs, vw, vr, x0];
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
    ) -> Option<([f64; 5], [f64; RESULTS_SIZE])> {
        let x = self.newton(mun, initial_x)?;
        let x_final = [x[0], x[1], x[2], x[3], x[4]];

        // Magnetização: M B = B dP/dB a mu fixo (Landau), calculada com duas
        // soluções em B (1 +- delta) partindo da solução convergida. Entra na
        // pressão perpendicular da matéria, P_perp = P_par - M B (Ferrer et
        // al. 2010; Strickland, Dexheimer & Menezes 2012).
        self.magnetization_b = if self.b > 0.0 {
            self.magnetization_times_field(mun, &x_final)?
        } else {
            0.0
        };

        // garante estado físico final consistente
        let _ = self.funcv(&x_final);
        self.assemble_point(x_final, field_g)
    }

    /// Pressão termodinâmica da matéria P = -Omega (MeV/fm^3) no estado atual.
    fn matter_pressure_mev_fm3(&self, x: &[f64; 5]) -> f64 {
        let (mue, vsigma, vomega, vrho, _) = self.mapping(x);
        let (_, press) = crate::core::eos::compute(self, mue, vsigma, vomega, vrho);
        press * self.m_nuc * (self.m_nuc / HBAR_C).powi(3)
    }

    fn magnetization_times_field(&mut self, mun: f64, x: &[f64; 5]) -> Option<f64> {
        let (b, delta) = (self.b, 1e-5);
        let mut pressure_at = |engine: &mut Self, scale: f64| -> Option<f64> {
            engine.b = b * scale;
            let xs = engine.newton(mun, x)?;
            let xs = [xs[0], xs[1], xs[2], xs[3], xs[4]];
            let _ = engine.funcv(&xs);
            Some(engine.matter_pressure_mev_fm3(&xs))
        };
        let p_plus = pressure_at(self, 1.0 + delta);
        let p_minus = pressure_at(self, 1.0 - delta);
        self.b = b;
        Some((p_plus? - p_minus?) / (2.0 * delta))
    }

    /// Newton amortecido para (mu_e, sigma, omega, rho, X_0) a mu_n fixo.
    fn newton(&mut self, mun: f64, initial_x: &[f64]) -> Option<StateVector> {
        self.mun = mun;

        let (mue0, vs0, vw0, vr0, x00) = self.mapping(initial_x);
        let mut x = StateVector::from_column_slice(&[mue0, vs0, vw0, vr0, x00]);
        let tolerance = 1e-10;
        let max_iterations = 100;
        let mut converged = false;

        for _ in 0..max_iterations {
            let f_val = StateVector::from_column_slice(&self.funcv(x.as_slice()));
            let f_norm = f_val.norm();

            if f_norm.is_finite() && f_norm < tolerance {
                converged = true;
                break;
            }
            if !f_norm.is_finite() {
                break;
            }

            let mut j_matrix = Jacobian::zeros();
            for i in 0..5 {
                // Passo pequeno: perto de limiares de partículas, uma sonda
                // maior pode atravessar dois ramos de densidade.
                let h = 1e-8 * (x[i].abs() + 1e-2);
                let mut x_temp = x;
                x_temp[i] += h;
                let f_temp = StateVector::from_column_slice(&self.funcv(x_temp.as_slice()));
                j_matrix.set_column(i, &((f_temp - f_val) / h));
            }

            let delta_x = match j_matrix.lu().solve(&(-f_val)) {
                Some(step) if step.iter().all(|value| value.is_finite()) => step,
                _ => break,
            };

            let mut alpha = 1.0;
            let mut step_accepted = false;
            for _ in 0..15 {
                let x_try = x + alpha * delta_x;
                let f_new_norm = StateVector::from_column_slice(&self.funcv(x_try.as_slice())).norm();
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
                // Perto do limiar vácuo-matéria as densidades são apenas
                // diferenciáveis por partes; um passo amortecido evita que a
                // continuação fique presa na raiz de vácuo.
                let x_fallback = x + 0.001 * delta_x;
                if !x_fallback.iter().all(|value| value.is_finite()) {
                    break;
                }
                x = x_fallback;
                let _ = self.funcv(x.as_slice());
            }
        }

        // Uma raiz alcançada pelo último passo permitido não é descartada.
        if !converged {
            let final_norm = StateVector::from_column_slice(&self.funcv(x.as_slice())).norm();
            converged = final_norm.is_finite() && final_norm < tolerance;
        }

        if converged { Some(x) } else { None }
    }

    /// Monta a linha de saída a partir do estado convergido (já aplicado).
    fn assemble_point(
        &mut self,
        x_final: [f64; 5],
        field_g: Option<f64>,
    ) -> Option<([f64; 5], [f64; RESULTS_SIZE])> {
        let (mue, vsigma, vomega, vrho, _) = self.mapping(&x_final);
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
        let stress = if self.field_stress_in_eos() {
            magnetic_stress(self.nlem, b_local_g)
        } else {
            crate::core::magnetic::MagneticStress {
                energy: 0.0,
                p_parallel: 0.0,
                p_perpendicular: 0.0,
            }
        };
        self.field_stress = stress;
        let ebsd = stress.energy;
        // Pressão do campo + termo de magnetização da matéria: anisotrópica,
        // P_perp = P - M B; isotrópica (campo emaranhado), P - (2/3) M B.
        let magnetization_weight = match self.topology {
            MagneticTopology::Anisotropic => 1.0,
            MagneticTopology::Isotropic => 2.0 / 3.0,
        };
        let pmag_effective =
            stress.pressure(self.topology) - magnetization_weight * self.magnetization_b;

        let ener_final = ener_conv + ebsd;
        let press_final = press_conv + pmag_effective;

        // P_perp pode ficar negativa na matéria diluída em campo forte
        // (M B > P); a validade do ponto é julgada pela pressão sem o termo de
        // magnetização. A TOV descarta os trechos em que P_perp não cresce.
        self.stability_pressure = press_final + magnetization_weight * self.magnetization_b;
        // Tolerância absoluta de 1e-12 MeV/fm^3: no limiar vácuo-matéria o
        // cancelamento em P = sum(mu n) - eps pode dar ~-1e-22.
        if ener_final >= -1e-12 && self.stability_pressure >= -1e-12 {
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
            if self.dark_enabled {
                let dark_vector_energy = self.dark_vector_energy_density();
                result[21] = self.n_chi * density_factor;
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
            }
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
    fn anomalous_moments_keep_pressure_consistent() {
        // Com AMM e B = 1e18 G (Landau + Pauli nos prótons, espectro anisotrópico
        // nos nêutrons), dP/dmu_n = n_B a mu_e fixo continua valendo; e o spin com
        // momento ao longo de B (kappa_n < 0: spin "down") é mais populado.
        let mut engine = HadronsMatter::new(GM1, 1e18).with_anomalous_moments(true);
        let (mun, mue, dmu) = (1.08, 0.12, 1e-5);
        let mut guess = [0.4, 0.4, -0.05];
        for mu in [1.20, 1.15, 1.10, mun] {
            guess = solve_fields(&mut engine, mu, mue, guess).expect("continuation").0;
        }
        let (fields, nb, _, _) = solve_fields(&mut engine, mun, mue, guess).unwrap();
        let (a, m_star, ef) = (engine.amm_b[0] * engine.b, engine.m_eff[0], engine.ef_b[0]);
        let up = crate::core::particles::neutral_amm_spin(m_star, ef, a).unwrap().density;
        let down = crate::core::particles::neutral_amm_spin(m_star, ef, -a).unwrap().density;
        assert!(down > up, "n_down = {down}, n_up = {up}");
        let (_, _, _, p_plus) = solve_fields(&mut engine, mun + dmu, mue, fields).unwrap();
        let (_, _, _, p_minus) = solve_fields(&mut engine, mun - dmu, mue, fields).unwrap();
        let dp_dmu = (p_plus - p_minus) / (2.0 * dmu);
        assert!(((dp_dmu - nb) / nb).abs() < 1e-5, "dP/dmu_n = {dp_dmu}, n_B = {nb}");
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
