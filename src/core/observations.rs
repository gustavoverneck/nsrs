// src/core/observations.rs
//
// Nível 3: confronto de modelos com vínculos observacionais e empíricos
// (input/observations/constraints.csv). Cada vínculo é convertido numa
// distância d em desvios-padrão (barras assimétricas; intervalos de 90% são
// convertidos para 1 sigma dividindo por 1.645). Classificação:
// d <= 1 compatível, 1 < d <= 2 tensão, d > 2 excluído.

use crate::core::nuclear::SaturationProperties;
use crate::core::tov_solver::StarProperties;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstraintKind {
    /// Massa de pulsar: a massa máxima do modelo deve alcançá-la.
    MaxMass,
    /// Ponto (M, R) de NICER.
    MassRadius,
    /// Raio a uma massa fixa.
    RadiusAtMass,
    /// Deformabilidade de maré a uma massa fixa.
    TidalAtMass,
    NuclearN0,
    NuclearEnergyPerNucleon,
    NuclearK,
    NuclearJ,
    NuclearL,
}

impl ConstraintKind {
    fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "max_mass" => Self::MaxMass,
            "mass_radius" => Self::MassRadius,
            "radius_at_mass" => Self::RadiusAtMass,
            "tidal_at_mass" => Self::TidalAtMass,
            "nuclear_n0" => Self::NuclearN0,
            "nuclear_ea" => Self::NuclearEnergyPerNucleon,
            "nuclear_k" => Self::NuclearK,
            "nuclear_j" => Self::NuclearJ,
            "nuclear_l" => Self::NuclearL,
            _ => return None,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Constraint {
    pub kind: ConstraintKind,
    pub label: String,
    pub value: f64,
    pub err_minus: f64,
    pub err_plus: f64,
    /// Massa em que o vínculo vale (R ou Lambda a massa fixa).
    pub at_mass: Option<f64>,
    /// Massa observada com erros (pontos M-R): (M, -, +).
    pub mass: Option<(f64, f64, f64)>,
    pub credibility: String,
    pub reference: String,
    pub doi: String,
    pub arxiv: String,
}

impl Constraint {
    /// Fator que converte as barras publicadas para 1 sigma.
    fn sigma_scale(&self) -> f64 {
        if self.credibility.starts_with("90%") { 1.0 / 1.645 } else { 1.0 }
    }

    /// Desvio assimétrico (x - value)/sigma do lado correspondente.
    fn z(&self, x: f64) -> f64 {
        let err = if x >= self.value { self.err_plus } else { self.err_minus };
        (x - self.value) / (err * self.sigma_scale())
    }
}

pub fn load_constraints<P: AsRef<Path>>(path: P) -> Result<Vec<Constraint>, String> {
    let mut reader = csv::Reader::from_path(path.as_ref()).map_err(|e| e.to_string())?;
    let headers = reader.headers().map_err(|e| e.to_string())?.clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| format!("missing column '{name}'"))
    };
    let (c_kind, c_label, c_value, c_em, c_ep) =
        (col("kind")?, col("label")?, col("value")?, col("err_minus")?, col("err_plus")?);
    let (c_at, c_m, c_mm, c_mp) =
        (col("at_mass_msun")?, col("mass_msun")?, col("mass_err_minus")?, col("mass_err_plus")?);
    let (c_cred, c_ref, c_doi, c_arxiv) =
        (col("credibility")?, col("reference")?, col("doi")?, col("arxiv")?);

    let mut constraints = Vec::new();
    for (line, record) in reader.records().enumerate() {
        let record = record.map_err(|e| e.to_string())?;
        let field = |i: usize| record.get(i).unwrap_or("").trim();
        let number = |i: usize| -> Result<f64, String> {
            field(i).parse::<f64>().map_err(|_| format!("row {}: bad number '{}'", line + 2, field(i)))
        };
        let optional = |i: usize| -> Result<Option<f64>, String> {
            if field(i).is_empty() { Ok(None) } else { number(i).map(Some) }
        };
        let kind = ConstraintKind::parse(field(c_kind))
            .ok_or_else(|| format!("row {}: unknown kind '{}'", line + 2, field(c_kind)))?;
        let mass = match optional(c_m)? {
            Some(m) => Some((m, number(c_mm)?, number(c_mp)?)),
            None => None,
        };
        constraints.push(Constraint {
            kind,
            label: field(c_label).to_string(),
            value: number(c_value)?,
            err_minus: number(c_em)?,
            err_plus: number(c_ep)?,
            at_mass: optional(c_at)?,
            mass,
            credibility: field(c_cred).to_string(),
            reference: field(c_ref).to_string(),
            doi: field(c_doi).to_string(),
            arxiv: field(c_arxiv).to_string(),
        });
    }
    Ok(constraints)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Compatible,
    Tension,
    Excluded,
}

impl Status {
    pub fn from_distance(d: f64) -> Self {
        if d <= 1.0 {
            Status::Compatible
        } else if d <= 2.0 {
            Status::Tension
        } else {
            Status::Excluded
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Status::Compatible => "compatível",
            Status::Tension => "tensão",
            Status::Excluded => "excluído",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Assessment {
    pub label: String,
    /// Valor do modelo comparado com o vínculo.
    pub model_value: f64,
    /// Distância em desvios-padrão (>= 0).
    pub distance: f64,
    pub status: Status,
}

/// Ramo estável: estrelas até a massa máxima.
pub fn stable_branch(stars: &[StarProperties]) -> &[StarProperties] {
    match (0..stars.len()).max_by(|&a, &b| stars[a].mass.total_cmp(&stars[b].mass)) {
        Some(i) => &stars[..=i],
        None => &[],
    }
}

/// Interpola uma grandeza do ramo estável na massa `mass`.
pub fn interpolate_at_mass(
    branch: &[StarProperties],
    mass: f64,
    value: impl Fn(&StarProperties) -> f64,
) -> Option<f64> {
    let k = (1..branch.len()).find(|&k| branch[k - 1].mass <= mass && branch[k].mass >= mass)?;
    let (a, b) = (&branch[k - 1], &branch[k]);
    let t = (mass - a.mass) / (b.mass - a.mass);
    Some(value(a) + t * (value(b) - value(a)))
}

/// Avalia um vínculo. `None` se o modelo não permite a comparação (p.ex.
/// massa pedida acima da massa máxima: nesse caso o vínculo de massa
/// máxima correspondente já registra a exclusão).
pub fn assess(
    constraint: &Constraint,
    stars: &[StarProperties],
    saturation: Option<&SaturationProperties>,
) -> Option<Assessment> {
    let branch = stable_branch(stars);
    let (model_value, distance) = match constraint.kind {
        ConstraintKind::MaxMass => {
            let m_max = branch.last()?.mass;
            (m_max, (-constraint.z(m_max)).max(0.0))
        }
        ConstraintKind::MassRadius => {
            let (m_obs, m_minus, m_plus) = constraint.mass?;
            let best = branch
                .iter()
                .map(|s| {
                    let zm = (s.mass - m_obs) / if s.mass >= m_obs { m_plus } else { m_minus };
                    let zr = constraint.z(s.radius);
                    ((zm * zm + zr * zr).sqrt(), s.radius)
                })
                .min_by(|a, b| a.0.total_cmp(&b.0))?;
            (best.1, best.0)
        }
        ConstraintKind::RadiusAtMass => {
            let r = interpolate_at_mass(branch, constraint.at_mass?, |s| s.radius)?;
            (r, constraint.z(r).abs())
        }
        ConstraintKind::TidalAtMass => {
            let l = interpolate_at_mass(branch, constraint.at_mass?, |s| s.tidal_deformability)?;
            (l, constraint.z(l).abs())
        }
        ConstraintKind::NuclearN0
        | ConstraintKind::NuclearEnergyPerNucleon
        | ConstraintKind::NuclearK
        | ConstraintKind::NuclearJ
        | ConstraintKind::NuclearL => {
            let s = saturation?;
            let x = match constraint.kind {
                ConstraintKind::NuclearN0 => s.n0,
                ConstraintKind::NuclearEnergyPerNucleon => s.energy_per_nucleon,
                ConstraintKind::NuclearK => s.incompressibility,
                ConstraintKind::NuclearJ => s.symmetry_energy,
                _ => s.symmetry_slope,
            };
            (x, constraint.z(x).abs())
        }
    };
    Some(Assessment {
        label: constraint.label.clone(),
        model_value,
        distance,
        status: Status::from_distance(distance),
    })
}
