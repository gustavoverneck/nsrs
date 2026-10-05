// src/core/magnetic.rs
//
// Perfis do campo magnético local e tensor de tensões do campo para
// eletrodinâmica linear (Maxwell) e não linear (ModMax, logarítmica).
// Todas as intensidades de campo são em Gauss.

use crate::core::constants::{BDD_ALPHAA, BDD_BETAA, M_NUCLEON};
use crate::core::physics::{MagneticTopology, NlemModel};

/// Campo crítico do elétron, m_e^2 c^3 / (e hbar), em Gauss.
pub const B_CRIT_ELECTRON_G: f64 = 4.414e13;

/// Conversão de erg/cm^3 (unidades gaussianas) para MeV/fm^3.
const ERG_CM3_PER_MEV_FM3: f64 = 1.602176634e33;

/// Campo de superfície usado nos perfis dependentes de densidade: 1e15 G,
/// a maior intensidade de superfície observada em magnetares, adotada nos
/// trabalhos recentes (ver Dexheimer et al., Astron. Nachr. 338, 1052 (2017)).
pub const B_SURFACE_G: f64 = 1e15;

/// Coeficientes (a, b, c) do ajuste quadrático de Dexheimer et al.,
/// Phys. Lett. B 773, 487 (2017), Tabela 2, para B*(mu_B) = (a + b mu_B +
/// c mu_B^2) mu / B_c^2, com mu_B em MeV e o dipolo mu em A m^2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DexheimerFit {
    /// Estrela de massa bariônica 2.2 M_sun (~2 M_sun gravitacional).
    BaryonMass22,
    /// Estrela de massa bariônica 1.6 M_sun (~1.4 M_sun gravitacional).
    BaryonMass16,
}

impl DexheimerFit {
    fn coefficients(self) -> (f64, f64, f64) {
        match self {
            DexheimerFit::BaryonMass22 => (-7.69e-1, 1.20e-3, -3.46e-7),
            DexheimerFit::BaryonMass16 => (-1.02, 1.58e-3, -4.85e-7),
        }
    }

    /// Vértice do polinômio (c < 0): acima dele o ajuste decresceria, o que
    /// é um artefato da extrapolação; o campo é mantido constante ali.
    pub fn vertex_mev(self) -> f64 {
        let (_, b, c) = self.coefficients();
        -b / (2.0 * c)
    }
}

/// Como a intensidade local do campo varia dentro da estrela.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FieldProfile {
    /// Comportamento legado: os níveis de Landau usam o campo central `bg`
    /// em todas as densidades; a energia magnética usa o perfil BDD com
    /// B_surf = 1e15 G.
    Constant,
    /// Perfil dependente de densidade de Bandyopadhyay, Chakrabarty & Pal,
    /// PRL 79, 2176 (1997): B = B_surf + B0 [1 - exp(-beta (n_B/n0)^gamma)].
    /// O mesmo campo local entra nos níveis de Landau e na energia.
    Bdd {
        b_surf_g: f64,
        b0_g: f64,
        beta: f64,
        gamma: f64,
    },
    /// Perfil polar ajustado a soluções de Einstein-Maxwell por Dexheimer
    /// et al., Phys. Lett. B 773, 487 (2017), Eq. (1):
    /// B = (a + b mu_B + c mu_B^2) mu / B_c.
    Dexheimer2017 { dipole_am2: f64, fit: DexheimerFit },
}

impl FieldProfile {
    /// Perfil BDD com os parâmetros usuais (beta = 0.01, gamma = 3) e
    /// B_surf = 1e15 G.
    pub fn bdd(b0_g: f64) -> Self {
        FieldProfile::Bdd {
            b_surf_g: B_SURFACE_G,
            b0_g,
            beta: BDD_BETAA,
            gamma: BDD_ALPHAA,
        }
    }

    /// Campo local em Gauss, dado n_B/n0 e mu_B (MeV). `None` para o perfil
    /// `Constant`, cujo tratamento é específico de cada engine.
    pub fn local_field_g(&self, nb_over_n0: f64, mu_b_mev: f64) -> Option<f64> {
        match *self {
            FieldProfile::Constant => None,
            FieldProfile::Bdd {
                b_surf_g,
                b0_g,
                beta,
                gamma,
            } => Some(bdd_field_g(b_surf_g, b0_g, beta, gamma, nb_over_n0)),
            FieldProfile::Dexheimer2017 { dipole_am2, fit } => {
                // Abaixo da superfície (vácuo) usa-se mu_B = m_N; acima do
                // vértice o campo é congelado (ver `vertex_mev`).
                let mu = mu_b_mev.clamp(M_NUCLEON, fit.vertex_mev());
                let (a, b, c) = fit.coefficients();
                Some(((a + b * mu + c * mu * mu) * dipole_am2 / B_CRIT_ELECTRON_G).max(0.0))
            }
        }
    }

    /// Se a energia e as tensões do próprio campo entram, por padrão, na EoS
    /// usada pela TOV. `Constant` e `Bdd`: sim (prática da literatura com o
    /// perfil BDD). `Dexheimer2017`: não; o ajuste vem de soluções de
    /// Einstein-Maxwell, nas quais o campo é tratado na estrutura, e é
    /// destinado apenas à EoS microscópica ("to be used as input in
    /// microscopic calculations"). Com B(m_N) ~ 4e17 G na superfície, somar
    /// B^2/8pi ~ 3 MeV/fm^3 à EoS criaria um envelope sem matéria.
    pub fn field_stress_in_eos(&self) -> bool {
        !matches!(self, FieldProfile::Dexheimer2017 { .. })
    }

    /// Se o campo depende de n_B e precisa de iteração em cada ponto.
    pub fn depends_on_density(&self) -> bool {
        matches!(self, FieldProfile::Bdd { .. })
    }
}

pub fn bdd_field_g(b_surf_g: f64, b0_g: f64, beta: f64, gamma: f64, nb_over_n0: f64) -> f64 {
    b_surf_g + b0_g * (1.0 - (-beta * nb_over_n0.max(0.0).powf(gamma)).exp())
}

/// Densidade de energia e pressões do campo magnético puro, em MeV/fm^3.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MagneticStress {
    pub energy: f64,
    /// Pressão ao longo das linhas de campo.
    pub p_parallel: f64,
    /// Pressão perpendicular às linhas de campo.
    pub p_perpendicular: f64,
}

impl MagneticStress {
    /// Pressão efetiva usada na TOV. Anisotrópica: P_perp (para Maxwell,
    /// P_perp = eps). Isotrópica (campo emaranhado): média sobre direções,
    /// (P_par + 2 P_perp)/3 (para Maxwell, eps/3).
    pub fn pressure(&self, topology: MagneticTopology) -> f64 {
        match topology {
            MagneticTopology::Anisotropic => self.p_perpendicular,
            MagneticTopology::Isotropic => (self.p_parallel + 2.0 * self.p_perpendicular) / 3.0,
        }
    }
}

/// Tensor de tensões de um campo magnético estático puro B (Gauss).
///
/// Para uma Lagrangiana L(B), a densidade de energia é eps = -L e
/// H = d eps / dB. As tensões são sigma_ij = H_i B_j - delta_ij (H B - eps),
/// logo P_par = -eps e P_perp = H B - eps (Soleng, PRD 52, 6178 (1995),
/// Eq. 3, para o caso logarítmico). Para Maxwell e ModMax, eps é quadrático
/// em B, H B = 2 eps e P_perp = eps.
pub fn magnetic_stress(nlem: NlemModel, b_gauss: f64) -> MagneticStress {
    let eps_maxwell = b_gauss * b_gauss / (8.0 * std::f64::consts::PI) / ERG_CM3_PER_MEV_FM3;
    // eps/eps_Maxwell; H B/eps_Maxwell = 2 H/B.
    let energy_ratio = match nlem {
        NlemModel::Maxwell => 1.0,
        NlemModel::Modmax(gamma) => (-gamma).exp(),
        NlemModel::Log(xi_gauss) => {
            // eps = xi^2 ln(1 + x), com x = B^2 / (2 xi^2).
            let x = b_gauss * b_gauss / (2.0 * xi_gauss * xi_gauss);
            if x < 1e-12 { 1.0 - 0.5 * x } else { x.ln_1p() / x }
        }
    };
    let hb_ratio = 2.0 * nlem.h_over_b(b_gauss);
    let energy = eps_maxwell * energy_ratio;
    MagneticStress {
        energy,
        p_parallel: -energy,
        p_perpendicular: eps_maxwell * hb_ratio - energy,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * a.abs().max(b.abs())
    }

    #[test]
    fn maxwell_stress_reproduces_legacy_topologies() {
        let s = magnetic_stress(NlemModel::Maxwell, 1e18);
        // B^2/(8 pi) for 1e18 G in MeV/fm^3.
        assert!(close(s.energy, 1e36 / (8.0 * std::f64::consts::PI) / 1.602176634e33, 1e-14));
        assert!(close(s.pressure(MagneticTopology::Anisotropic), s.energy, 1e-14));
        assert!(close(s.pressure(MagneticTopology::Isotropic), s.energy / 3.0, 1e-14));
        assert!(close(s.p_parallel, -s.energy, 1e-14));
    }

    #[test]
    fn modmax_stress_is_rescaled_maxwell() {
        let m = magnetic_stress(NlemModel::Maxwell, 3e17);
        let s = magnetic_stress(NlemModel::Modmax(0.7), 3e17);
        let f = (-0.7f64).exp();
        assert!(close(s.energy, f * m.energy, 1e-14));
        assert!(close(s.p_perpendicular, f * m.p_perpendicular, 1e-14));
    }

    #[test]
    fn log_stress_follows_from_the_energy_density() {
        let xi = 2e17;
        let eps = |b: f64| magnetic_stress(NlemModel::Log(xi), b).energy;
        for b in [1e16, 1e17, 2e17, 5e17, 3e18] {
            let s = magnetic_stress(NlemModel::Log(xi), b);
            // H B = B d eps / dB (derivada numérica).
            let h = 1e-6 * b;
            let hb = b * (eps(b + h) - eps(b - h)) / (2.0 * h);
            assert!(close(s.p_perpendicular, hb - s.energy, 1e-7), "B = {b:e}");
            assert!(close(s.p_parallel, -s.energy, 1e-14));
        }
        // xi = B (mesmas unidades, Gauss): x = 1/2 e eps = 2 ln(3/2) eps_Maxwell.
        let m = magnetic_stress(NlemModel::Maxwell, xi);
        let s = magnetic_stress(NlemModel::Log(xi), xi);
        assert!(close(s.energy, 2.0 * 1.5f64.ln() * m.energy, 1e-14));
    }

    #[test]
    fn log_stress_recovers_maxwell_for_weak_fields() {
        let m = magnetic_stress(NlemModel::Maxwell, 1e15);
        let s = magnetic_stress(NlemModel::Log(1e25), 1e15);
        assert!(close(s.energy, m.energy, 1e-15));
        assert!(close(s.p_perpendicular, m.p_perpendicular, 1e-15));
    }

    #[test]
    fn dexheimer_profile_matches_eq1_units() {
        // Tabela 2 (M_B = 2.2): mu_B = 1000 MeV, mu = 3e32 A m^2.
        let p = FieldProfile::Dexheimer2017 {
            dipole_am2: 3e32,
            fit: DexheimerFit::BaryonMass22,
        };
        let expected = (-7.69e-1 + 1.20e-3 * 1000.0 - 3.46e-7 * 1e6) * 3e32 / 4.414e13;
        assert!(close(p.local_field_g(0.0, 1000.0).unwrap(), expected, 1e-14));
        assert!(close(expected, 5.777e17, 1e-3));
        // Congelado no vértice e na superfície.
        let v = DexheimerFit::BaryonMass22.vertex_mev();
        assert_eq!(p.local_field_g(0.0, v + 100.0), p.local_field_g(0.0, v));
        assert_eq!(p.local_field_g(0.0, 20.0), p.local_field_g(0.0, M_NUCLEON));
    }

    #[test]
    fn bdd_profile_limits() {
        let p = FieldProfile::bdd(1e18);
        assert_eq!(p.local_field_g(0.0, 0.0), Some(B_SURFACE_G));
        assert!(close(p.local_field_g(100.0, 0.0).unwrap(), B_SURFACE_G + 1e18, 1e-12));
    }
}
