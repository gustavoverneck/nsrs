use crate::core::constants::{M_NUCLEON, N0, RESULTS_SIZE};
// src/core/solver.rs
use crate::core::darkphotons::DarkPhotonsMatter;
use crate::core::hybrid::HybridMatter;
use crate::core::io_utils::write_eos_with_mr;
use crate::core::magnetic::MagneticStress;
use crate::core::physics::HadronsMatter;
use crate::core::quarks::QuarksMatter;
use crate::core::tov_solver::generate_mr_curve;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

// 1. Enum que define o modelo físico a ser resolvido
pub enum EngineMode {
    Hadrons(HadronsMatter),
    Quarks(QuarksMatter),
    Hybrid(HybridMatter),
    DarkPhotons(DarkPhotonsMatter),
}

/// Motivo pelo qual a varredura em mu_n terminou. `mu_n_mev` é o potencial
/// químico do último ponto aceito e `nb_over_n0` a densidade correspondente.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EosTermination {
    /// A malha chegou a `mun_sup`.
    ReachedUpperLimit,
    /// Corte intencional em densidade (n_B > 11 n0).
    DensityCap,
    /// Corte intencional por causalidade (c_s^2 = dP/deps > 1).
    Acausal { mu_n_mev: f64, nb_over_n0: f64 },
    /// Corte intencional: a massa efetiva do nêutron (coluna 16) chegou a
    /// M*/M <= 0, fora da validade física do modelo de campo médio.
    NonPositiveEffectiveMass { mu_n_mev: f64, nb_over_n0: f64 },
    /// Queda resolvida de eps ou P entre pontos consecutivos.
    NonMonotonic { mu_n_mev: f64, nb_over_n0: f64 },
    /// O Newton não convergiu nem com o passo mínimo; a EoS está truncada.
    ConvergenceFailure { mu_n_mev: f64, nb_over_n0: f64 },
}

impl EosTermination {
    /// Terminações que truncam a EoS sem que isso tenha sido pedido.
    pub fn is_anomalous(&self) -> bool {
        matches!(
            self,
            EosTermination::NonMonotonic { .. } | EosTermination::ConvergenceFailure { .. }
        )
    }
}

/// Grandezas por ponto da EoS que não cabem no formato de 34 colunas.
/// Exportadas em `<saída>_diag.txt` (ver `io_utils::write_diagnostics`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointDiagnostics {
    /// M B = B dP/dB a mu fixo (MeV/fm^3); já descontado da pressão
    /// perpendicular exportada na coluna 2.
    pub magnetization_b: f64,
    /// Pressão sem o termo de magnetização (MeV/fm^3); igual à coluna 2 para
    /// motores sem campo na matéria. Usada nos critérios de validade.
    pub stability_pressure: f64,
    /// Energia e tensões do campo somadas à EoS (MeV/fm^3); zero para motores
    /// sem campo. Permitem montar a pressão isotrópica para qualquer NLEM.
    pub field_stress: MagneticStress,
}

/// Passos do modo adaptativo, em mu_n/M_N: mínimo (~0.09 MeV), máximo na
/// matéria (~9.4 MeV) e no vácuo (~19 MeV).
const ADAPTIVE_H_MIN: f64 = 1e-4;
const ADAPTIVE_H_MAX: f64 = 1e-2;
const ADAPTIVE_H_VACUUM: f64 = 2e-2;
/// Sondagens adiante após falha do Newton no passo mínimo: até 2^7 H_MIN.
const ADAPTIVE_PROBES: usize = 7;

/// Controle do passo adaptativo em mu_n. `a` é o último ponto aceito, `older`
/// o penúltimo com o passo entre os dois, `h` o passo até `new`. Devolve se
/// `new` é aceito e o próximo passo.
///
/// No vácuo o passo é `ADAPTIVE_H_VACUUM`. O primeiro ponto com matéria só é
/// aceito com passo mínimo: o limiar é cercado pelo lado do vácuo. Na matéria
/// valem, para n_B e eps, a variação relativa por passo (até 20 tol) e o
/// desvio da extrapolação linear dos dois pontos anteriores (até tol), que
/// mede o erro da interpolação linear usada pela TOV e afina o passo nas
/// quinas (híperons, níveis de Landau).
fn adaptive_step(
    tol: f64,
    h: f64,
    a: &[f64; RESULTS_SIZE],
    older: Option<(&[f64; RESULTS_SIZE], f64)>,
    new: &[f64; RESULTS_SIZE],
) -> (bool, f64) {
    let shrinkable = h > ADAPTIVE_H_MIN * (1.0 + 1e-9);
    if a[0] <= 0.0 && new[0] <= 0.0 {
        return (true, ADAPTIVE_H_VACUUM);
    }
    if a[0] <= 0.0 {
        return if shrinkable { (false, (0.25 * h).max(ADAPTIVE_H_MIN)) } else { (true, ADAPTIVE_H_MIN) };
    }
    // Escalas mínimas: n_B/n0 e eps (MeV/fm^3) perto do limiar.
    let mut worst: f64 = 0.0;
    let mut factor: f64 = 2.0;
    for (k, floor) in [(0, 1e-2), (1, 1.0)] {
        let scale = new[k].abs().max(floor);
        let change = (new[k] - a[k]).abs() / scale / (20.0 * tol);
        worst = worst.max(change);
        factor = factor.min(0.9 / change.max(1e-12));
        if let Some((b, h_ab)) = older.filter(|(b, h_ab)| b[0] > 0.0 && *h_ab > 0.0) {
            let predicted = a[k] + (a[k] - b[k]) * h / h_ab;
            let curvature = (new[k] - predicted).abs() / scale / tol;
            worst = worst.max(curvature);
            factor = factor.min(0.9 / curvature.max(1e-12).sqrt());
        }
    }
    let h_next = (h * factor.clamp(0.2, 2.0)).clamp(ADAPTIVE_H_MIN, ADAPTIVE_H_MAX);
    if worst > 1.0 && shrinkable {
        (false, h_next.min(0.9 * h))
    } else {
        (true, h_next)
    }
}

pub struct Solver {
    engine: EngineMode,
    termination: Option<EosTermination>,
    diagnostics: Vec<PointDiagnostics>,
}

impl EngineMode {
    /// Arquivo em que a EoS é gravada, se houver.
    pub fn eos_output(&self) -> Option<&str> {
        match self {
            EngineMode::Hadrons(h) | EngineMode::DarkPhotons(h) => h.eos_output.as_deref(),
            EngineMode::Quarks(q) => q.eos_output.as_deref(),
            EngineMode::Hybrid(h) => h.eos_output.as_deref(),
        }
    }
}

impl Solver {
    pub fn new(engine: EngineMode) -> Self {
        Solver {
            engine,
            termination: None,
            diagnostics: Vec::new(),
        }
    }

    /// Motivo do término da última chamada a `solve()`.
    pub fn termination(&self) -> Option<EosTermination> {
        self.termination
    }

    /// Diagnósticos por linha da última chamada a `solve()` (mesma ordem das
    /// linhas da EoS).
    pub fn diagnostics(&self) -> &[PointDiagnostics] {
        &self.diagnostics
    }

    pub fn solve(&mut self) -> Vec<[f64; RESULTS_SIZE]> {
        // Obtém limites e metadados conforme o modo da engine
        let (mun_inf, mun_sup, n, bg_val) = match &self.engine {
            EngineMode::Hadrons(h) => (h.mun_inf, h.mun_sup, h.n_points, h.bg),
            EngineMode::Quarks(q) => (q.mun_inf, q.mun_sup, q.n_points, q.bg),
            EngineMode::Hybrid(hyb) => (
                hyb.hadrons.mun_inf,
                hyb.hadrons.mun_sup,
                hyb.hadrons.n_points,
                hyb.hadrons.bg,
            ),
            EngineMode::DarkPhotons(d) => (d.mun_inf, d.mun_sup, d.n_points, 0.0),
        };

        let initial_dmub = (mun_sup - mun_inf) / (n - 1) as f64;
        let adaptive = match &self.engine {
            EngineMode::Hadrons(h) => h.step_tolerance,
            _ => None,
        };
        let mut dmub = initial_dmub;
        // mu_n do penúltimo ponto aceito (passo adaptativo).
        let mut prev_mun = f64::NAN;
        // Sondagens adiante desde a última falha (passo adaptativo).
        let mut probe = 0usize;
        // O próximo ponto que convergir é aceito sem o controle de passo (vem
        // de uma sondagem ou de um salto de ramo).
        let mut force_accept = false;
        // Maior passo no vácuo depois de um ponto com matéria rejeitado.
        let mut onset_cap = f64::INFINITY;
        let min_dmub = 1e-6; // passo mínimo aceitável

        let mut results: Vec<[f64; RESULTS_SIZE]> = Vec::with_capacity(n);
        let mut diagnostics: Vec<PointDiagnostics> = Vec::with_capacity(n);
        let mut last_visible_x = [0.0; 5];
        let mut last_dark_x = [0.0; 5];
        let mut last_mun = mun_inf; // último mun que convergiu
        // mu_n do último salto de ramo (no máximo um a partir de cada ponto).
        let mut jumped_from: Option<f64> = None;
        let mut mun = mun_inf;
        let mut termination = EosTermination::ReachedUpperLimit;
        // Apenas estas engines exportam M*/M do nêutron na coluna 16.
        let exports_effective_mass = matches!(
            self.engine,
            EngineMode::Hadrons(_) | EngineMode::DarkPhotons(_)
        );
        // Último ponto aceito, em (MeV, n_B/n0), para os diagnósticos.
        let last_point = |results: &[[f64; RESULTS_SIZE]]| {
            results
                .last()
                .map_or((mun_inf * M_NUCLEON, 0.0), |r| (r[17] * M_NUCLEON, r[0]))
        };

        while mun <= mun_sup + 1e-9 {
            let x_before = last_visible_x;
            // Tenta resolver o ponto com o mun atual
            let point_data =
                match &mut self.engine {
                    EngineMode::Hadrons(h_engine) => h_engine
                        .solve_point(mun, &last_visible_x)
                        .map(|(x, result)| {
                            last_visible_x = x;
                            result
                        }),
                    EngineMode::Quarks(q_engine) => q_engine.solve_point(mun),
                    EngineMode::Hybrid(hyb_engine) => hyb_engine
                        .solve_point(mun, &last_visible_x)
                        .map(|(x, result)| {
                            last_visible_x = x;
                            result
                        }),
                    EngineMode::DarkPhotons(d_engine) => {
                        d_engine.solve_point(mun, &last_dark_x).map(|(x, result)| {
                            last_dark_x = x;
                            result
                        })
                    }
                };

            if let Some(point_result) = point_data {
                // Passo adaptativo: um ponto que muda a EoS demais é refeito
                // com passo menor, a partir do último ponto aceito.
                let mut next_step = None;
                if let Some(tol) = adaptive {
                    let previous = results.last().map(|a| {
                        let older = (results.len() >= 2).then(|| (&results[results.len() - 2], last_mun - prev_mun));
                        (a, older)
                    });
                    let (accept, h_next) = match previous {
                        _ if force_accept => (true, ADAPTIVE_H_MIN),
                        Some((a, older)) => adaptive_step(tol, mun - last_mun, a, older, &point_result),
                        None => (true, if point_result[0] > 0.0 { ADAPTIVE_H_MIN } else { ADAPTIVE_H_VACUUM }),
                    };
                    // No vácuo, depois de cercar o limiar, o passo não volta a crescer.
                    let h_next = if results.last().is_some_and(|a| a[0] <= 0.0) {
                        if !accept {
                            onset_cap = h_next;
                        }
                        h_next.min(onset_cap)
                    } else {
                        h_next
                    };
                    if !accept {
                        last_visible_x = x_before;
                        if let EngineMode::Hadrons(h) = &mut self.engine {
                            h.last_nb_over_n0 = results.last().map_or(0.0, |r| r[0]);
                        }
                        dmub = h_next;
                        mun = last_mun + dmub;
                        continue;
                    }
                    next_step = Some(h_next);
                    probe = 0;
                    force_accept = false;
                }
                prev_mun = last_mun;
                last_mun = mun;
                let point_diagnostics = match &self.engine {
                    EngineMode::Hadrons(h) | EngineMode::DarkPhotons(h) => PointDiagnostics {
                        magnetization_b: h.magnetization_b,
                        stability_pressure: h.stability_pressure,
                        field_stress: h.field_stress,
                    },
                    // Na fase hadrônica da híbrida, P_perp contém -M B; o teste
                    // de monotonia usa a pressão sem magnetização, como acima.
                    EngineMode::Hybrid(hyb) if hyb.hadronic_phase => PointDiagnostics {
                        magnetization_b: hyb.hadrons.magnetization_b,
                        stability_pressure: hyb.hadrons.stability_pressure,
                        field_stress: hyb.hadrons.field_stress,
                    },
                    _ => PointDiagnostics {
                        magnetization_b: 0.0,
                        stability_pressure: point_result[2],
                        field_stress: MagneticStress::default(),
                    },
                };

                // Para a integração se a densidade bariônica ultrapassar 11 N0.
                if point_result[0] * N0 > 11.0 * N0 {
                    termination = EosTermination::DensityCap;
                    break;
                }

                if exports_effective_mass && point_result[16] <= 0.0 {
                    let (mu_n_mev, nb_over_n0) = last_point(&results);
                    termination = EosTermination::NonPositiveEffectiveMass {
                        mu_n_mev,
                        nb_over_n0,
                    };
                    break;
                }

                // Verificação de estabilidade (dP/dE). Perto do limiar
                // vácuo-matéria, diferenças no nível do arredondamento não
                // representam uma instabilidade física e não devem encerrar a
                // continuação da EOS.
                if !results.is_empty() {
                    // Pressão sem magnetização: monótona em mu para matéria
                    // estável (dP/dmu = n_B), ao contrário de P_perp.
                    let prev = results.last().unwrap();
                    let prev_pressure = diagnostics.last().unwrap().stability_pressure;
                    let point_pressure = point_diagnostics.stability_pressure;
                    let de = point_result[1] - prev[1];
                    let dp = point_pressure - prev_pressure;
                    let resolved_matter = prev[0] > 1e-6 && point_result[0] > 1e-6;

                    let de_tol = 1e-10 * point_result[1].abs().max(prev[1].abs()).max(1.0);
                    let dp_tol = 1e-10 * point_pressure.abs().max(prev_pressure.abs()).max(1.0);

                    if resolved_matter && de > de_tol && dp > dp_tol {
                        let cs2 = dp / de;
                        if cs2 > 1.0 {
                            let (mu_n_mev, nb_over_n0) = last_point(&results);
                            termination = EosTermination::Acausal {
                                mu_n_mev,
                                nb_over_n0,
                            };
                            break;
                        }
                    } else if resolved_matter && (de < -de_tol || dp < -dp_tol) {
                        // Queda resolvida de energia ou pressão: esta sim é
                        // incompatível com a ramificação monótona usada pela TOV.
                        let (mu_n_mev, nb_over_n0) = last_point(&results);
                        termination = EosTermination::NonMonotonic {
                            mu_n_mev,
                            nb_over_n0,
                        };
                        break;
                    }
                }

                results.push(point_result);
                diagnostics.push(point_diagnostics);

                if let Some(h_next) = next_step {
                    dmub = h_next;
                    mun += dmub;
                } else {
                    // Avança para o próximo mun
                    mun += dmub;

                    // Se o passo estava reduzido, tenta aumentá-lo gradualmente
                    if dmub < initial_dmub {
                        dmub = (dmub * 1.5).min(initial_dmub);
                    }
                }
            } else {
                // Falha na convergência: reduz o passo e tenta novamente a partir do último sucesso
                dmub *= 0.5;
                // Modo adaptativo: abaixo do passo mínimo, sonda adiante com
                // passos dobrados a partir do último ponto aceito (como a
                // malha uniforme, que salta a faixa em que o Newton falha).
                if adaptive.is_some() && (dmub < ADAPTIVE_H_MIN || probe > 0) && !results.is_empty() {
                    probe += 1;
                    if probe <= ADAPTIVE_PROBES {
                        dmub = ADAPTIVE_H_MIN * 2f64.powi(probe as i32);
                        mun = last_mun + dmub;
                        force_accept = true;
                        continue;
                    }
                    probe = 0;
                    force_accept = false;
                    dmub = 0.0;
                }
                if dmub < min_dmub {
                    // Antes de desistir, tenta saltar para o ramo denso: sob
                    // campo forte o ramo de baixa densidade pode dobrar logo
                    // depois do limiar (transição de primeira ordem em mu).
                    let last_nb = results.last().map_or(0.0, |r| r[0]);
                    let jump = (jumped_from != Some(last_mun))
                        .then(|| {
                            self.jump_to_dense_branch(last_mun, initial_dmub, &last_visible_x, last_nb)
                        })
                        .flatten();
                    if let Some((mu_jump, x)) = jump {
                        jumped_from = Some(last_mun);
                        last_visible_x = x;
                        mun = mu_jump;
                        dmub = if adaptive.is_some() { ADAPTIVE_H_MIN } else { initial_dmub };
                        force_accept = adaptive.is_some();
                        continue;
                    }
                    let (mu_n_mev, nb_over_n0) = last_point(&results);
                    termination = EosTermination::ConvergenceFailure {
                        mu_n_mev,
                        nb_over_n0,
                    };
                    break;
                }
                mun = if results.is_empty() {
                    mun_inf
                } else {
                    last_mun + dmub
                };
            }
        }

        let output_path = self.engine.eos_output().map(str::to_string);

        self.termination = Some(termination);
        if termination.is_anomalous() {
            eprintln!(
                "Aviso: EoS truncada antes de mun_sup = {:.2} MeV (B = {:.3e} G, saída = {}): {:?}",
                mun_sup * M_NUCLEON,
                bg_val,
                output_path.as_deref().unwrap_or("-"),
                termination
            );
        }

        if let Some(path) = &output_path {
            if let Err(error) =
                crate::core::io_utils::write_diagnostics(&results, &diagnostics, path)
            {
                eprintln!("failed to write diagnostics for '{}': {error}", path);
            }
        }
        self.diagnostics = diagnostics;

        if let Some(path) = output_path {
            let eps_arr: Vec<f64> = results.iter().map(|r| r[1]).collect();
            let p_arr: Vec<f64> = results.iter().map(|r| r[2]).collect();
            let rho_arr: Vec<f64> = results.iter().map(|r| r[0]).collect();
            // Sequência com crosta, maré, momento de inércia e redshift.
            let stars = crate::core::tov_solver::generate_star_sequence(
                &eps_arr, &p_arr, &rho_arr, true,
            );
            if let Err(error) = crate::core::io_utils::write_stars(&stars, &path) {
                eprintln!("failed to write star sequence for '{}': {error}", path);
            }
            let (masses, radii, b_masses, pc_list) =
                generate_mr_curve(&eps_arr, &p_arr, &rho_arr, false);
            if let Err(error) =
                write_eos_with_mr(&results, &masses, &radii, &b_masses, &pc_list, &path)
            {
                eprintln!("failed to write EOS output '{}': {error}", path);
            }
        }

        results
    }

    /// Procura, a partir de `last_mun`, uma raiz no ramo denso quando a
    /// continuação trava: avança mu_n em passos de `step` e parte de chutes
    /// com os campos mesônicos ampliados. Só para hádrons e só enquanto a
    /// matéria ainda é diluída (n_B < n0); devolve (mu_n, estado) da primeira
    /// raiz com n_B pelo menos 1,5 vez maior que o último ponto aceito.
    /// Ainda no vácuo, aceita a primeira raiz com matéria: com AMM e campo
    /// forte há uma faixa estreita de mu_n logo acima do limiar (~1 MeV com
    /// GM1 a 1e18 G) em que o Newton não converge, e a malha que cai nela
    /// encerrava a EoS sem matéria.
    fn jump_to_dense_branch(
        &mut self,
        last_mun: f64,
        step: f64,
        last_x: &[f64; 5],
        last_nb: f64,
    ) -> Option<(f64, [f64; 5])> {
        let EngineMode::Hadrons(engine) = &mut self.engine else { return None };
        if last_nb >= 1.0 {
            return None;
        }
        let saved_nb = engine.last_nb_over_n0;
        if last_nb <= 0.0 {
            for j in 1..=20 {
                let mu_try = last_mun + j as f64 * step;
                engine.last_nb_over_n0 = saved_nb;
                if let Some((x, row)) = engine.solve_point(mu_try, last_x) {
                    if row[0] > 0.0 {
                        engine.last_nb_over_n0 = saved_nb;
                        return Some((mu_try, x));
                    }
                }
            }
            engine.last_nb_over_n0 = saved_nb;
            return None;
        }
        if last_nb <= 1e-6 {
            return None;
        }
        for j in 1..=20 {
            let mu_try = last_mun + j as f64 * step;
            for scale in [2.0, 4.0, 8.0, 16.0, 32.0] {
                let mut guess = *last_x;
                for value in &mut guess[1..4] {
                    *value *= scale;
                }
                engine.last_nb_over_n0 = saved_nb;
                if let Some((x, row)) = engine.solve_point(mu_try, &guess) {
                    if row[0] >= 1.5 * last_nb {
                        engine.last_nb_over_n0 = saved_nb;
                        return Some((mu_try, x));
                    }
                }
            }
        }
        engine.last_nb_over_n0 = saved_nb;
        None
    }

    /// Resolve múltiplas EoS de forma paralela usando Rayon.
    pub fn solve_parallel(
        engines: Vec<EngineMode>,
        num_threads: usize,
    ) -> Vec<Vec<[f64; RESULTS_SIZE]>> {
        let pb = ProgressBar::new(engines.len() as u64);
        let style = ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar());
        pb.set_style(style);
        pb.set_message("NSRS");

        // 1. Criamos um construtor de pool de threads personalizado
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("Falha ao criar o ThreadPool do Rayon");

        // 2. Executamos o processamento paralelo dentro deste pool específico
        let results = pool.install(|| {
            // with_max_len(1): cada engine é uma tarefa; uma thread livre pega a
            // próxima assim que termina (sem blocos contíguos de custo desigual).
            engines
                .into_par_iter()
                .with_max_len(1)
                .map(|engine| {
                    let written = engine.eos_output().is_some();
                    let mut solver = Solver::new(engine);
                    let result = solver.solve();
                    pb.inc(1);
                    // Varreduras de produção gravam cada EoS em disco e não
                    // mantêm todas em memória.
                    if written { Vec::new() } else { result }
                })
                .collect()
        });

        pb.finish_and_clear();
        results
    }

    /// Exporta os resultados da EoS para um arquivo formatado.
    pub fn write_eos(results: &[[f64; RESULTS_SIZE]], filename: &str) -> std::io::Result<()> {
        use std::io::Write; // Garante que o trait Write está no escopo
        let mut file = std::fs::File::create(filename)?;

        for data in results.iter() {
            // Transforma todas as colunas EOS em uma String separada por espaços.
            let line = data
                .iter()
                .map(|val| format!("{:12.5e}", val))
                .collect::<Vec<String>>()
                .join(" ");

            // 2. Escreve a linha inteira no arquivo
            writeln!(file, "{}", line)?;
        }
        Ok(())
    }
}
