// src/bin/hybrid.rs

use nsrs::{
    Artist, EngineMode, GM1, HadronsMatter, HybridMatter, QuarksMatter, Solver,
    constants::RESULTS_SIZE, generate_mr_curve,
};
fn main() {
    let bg = 1e15;
    let bag_constant = 85.38; // MeV/fm³

    // 1. Instancia os componentes base
    let hadron_engine = HadronsMatter::new(GM1, bg)
        .with_limits(0.0, 2.5)
        .with_points(1500);

    let mut quark_engine = QuarksMatter::new(bag_constant, 0.0)
        .with_limits(1.5, 5.0)
        .with_points(3000);

    quark_engine.gv = 0.0;

    // 2. Resolve a Estrela de Hádrons Pura
    println!("Resolvendo EoS Hadrônica (GM1)...");
    let mut solver_h = Solver::new(EngineMode::Hadrons(hadron_engine.clone()));
    let res_h = solver_h.solve();
    let (eps_h, p_h, rho_h) = extract_eos(&res_h);
    let (m_h, r_h, _b_h, _) = generate_mr_curve(&eps_h, &p_h, &rho_h, true);

    // 3. Resolve a Estrela de Quarks Pura (MIT Bag)
    println!("Resolvendo EoS de Quarks (MIT Bag)...");
    let mut solver_q = Solver::new(EngineMode::Quarks(quark_engine.clone()));
    let res_q = solver_q.solve();
    let (eps_q, p_q, rho_q) = extract_eos(&res_q);
    let (m_q, r_q, _b_q, _) = generate_mr_curve(&eps_q, &p_q, &rho_q, true);

    // 4. Resolve a Estrela Híbrida (Maxwell Construction)
    println!("Resolvendo EoS Híbrida (Maxwell)...");
    let hybrid_engine = HybridMatter::new(hadron_engine, quark_engine);
    let mut solver_hyb = Solver::new(EngineMode::Hybrid(hybrid_engine));
    let res_hyb = solver_hyb.solve();
    let (eps_hyb, p_hyb, rho_hyb) = extract_eos(&res_hyb);
    let (m_hyb, r_hyb, _b_hyb, _) = generate_mr_curve(&eps_hyb, &p_hyb, &rho_hyb, true);

    // 5. Geração do Gráfico de Comparação M-R
    println!("Gerando gráficos de comparação...");
    let mut artist_mr = Artist::new("results/comparison_mr.svg", "Mass-Radius Comparison")
        .with_x_label("Radius [km]")
        .with_y_label("Mass [M\u{2299}]")
        .autoscale();

    artist_mr = artist_mr.add_curve(&r_h, &m_h, "Pure Hadron (GM1)");
    artist_mr = artist_mr.add_curve(&r_q, &m_q, "Pure Quark (MIT Bag)");
    artist_mr = artist_mr.add_curve(&r_hyb, &m_hyb, "Hybrid (Maxwell)");
    artist_mr.plot().ok();

    // 6. Geração do Gráfico de Comparação da EoS (P vs epsilon)
    let mut artist_eos = Artist::new("results/comparison_eos.svg", "Equation of State Comparison")
        .with_x_label("\u{03B5} [MeV/fm\u{00B3}]")
        .with_y_label("P [MeV/fm\u{00B3}]")
        .autoscale()
        .with_log_scale();

    artist_eos = artist_eos.add_curve(&eps_h, &p_h, "Hadron");
    artist_eos = artist_eos.add_curve(&eps_q, &p_q, "Quark");
    artist_eos = artist_eos.add_curve(&eps_hyb, &p_hyb, "Hybrid");
    artist_eos.plot().ok();

    println!("Sucesso! Gráficos gerados em results/comparison_mr.svg e results/comparison_eos.svg");
}

/// Função auxiliar para extrair Eps e P dos resultados do solver
fn extract_eos(results: &[[f64; RESULTS_SIZE]]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let eps = results.iter().map(|r| r[1]).collect();
    let p = results.iter().map(|r| r[2]).collect();
    let rho = results.iter().map(|r| r[0]).collect();
    (eps, p, rho)
}
