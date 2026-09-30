// Matéria de quarks (substitui os binários bag_model e hybrid).
// Em desenvolvimento; os gráficos vão para results/.

use nsrs::core::constants::RESULTS_SIZE;
use nsrs::core::plotting::{Artist, ColorScale, Palette};
use nsrs::core::quarks::QuarksMatter;
use nsrs::core::solver::{EngineMode, Solver};
use nsrs::core::tov_solver::generate_mr_curve;
use nsrs::{GM1, HadronsMatter, HybridMatter};

use crate::cli::create_dir;

/// `quarks bag`: estrelas de quarks (MIT bag), variando a constante de
/// sacola (gv = 0) e o acoplamento vetorial gv (B = 85.38 MeV/fm^3).
pub fn bag(_raw: &[String]) -> Result<(), String> {
    create_dir("results")?;
    let bg = 0.0; // Sem campo magnético para o MIT Bag puro

    // ====================================================================
    // CENÁRIO 1: Variando bag_constant (gv fixo em 0.0)
    // ====================================================================
    println!("======================================================");
    println!("[1/2] Cenário 1: Variação da Constante de Sacola (B)");
    println!("======================================================");

    // Lista de constantes de sacola para testar (em MeV/fm³)
    let bag_values = [57.5, 65.0, 75.0, 85.38, 100.0];

    let mut artist_mr_bag = Artist::new("results/bag_var_mr.svg", "MR - Varying Bag Constant")
        .with_x_label("Radius [km]")
        .with_y_label("Mass [M\u{2299}]")
        .with_x_range(6.0, 13.0);

    let mut artist_eos_bag = Artist::new("results/bag_var_eos.svg", "EoS - Varying Bag Constant")
        .with_x_label("\u{03B5} [MeV/fm\u{00B3}]")
        .with_y_label("P [MeV/fm\u{00B3}]")
        .autoscale()
        .with_log_scale();

    for &bag in &bag_values {
        println!("  -> Calculando bag_constant = {:.2}...", bag);

        // Cria a matéria de quarks e fixa gv = 0.0
        let mut quark_engine = QuarksMatter::new(bag, bg)
            .with_limits(1.0, 3.5)
            .with_points(2000);
        quark_engine.gv = 0.0;

        let mut solver_q = Solver::new(EngineMode::Quarks(quark_engine));
        let res_q = solver_q.solve();

        let (eps_q, p_q, rho_q) = extract_eos_filtered(&res_q);

        // Garante que o TOV só rode se houver matéria física
        if eps_q.len() > 10 {
            let (m_q, r_q, _b_q, _) = generate_mr_curve(&eps_q, &p_q, &rho_q, false);
            let label = format!("B = {:.2}", bag);

            artist_mr_bag = artist_mr_bag.add_curve(&r_q, &m_q, &label);
            artist_eos_bag = artist_eos_bag.add_curve(&eps_q, &p_q, &label);
        }
    }

    artist_mr_bag.plot().ok();
    artist_eos_bag.plot().ok();
    println!("Gráficos do Cenário 1 gerados: bag_var_mr.svg e bag_var_eos.svg");

    // ====================================================================
    // CENÁRIO 2: Variando gv de 0.0 a 2.0 (bag_constant fixa em 85.38)
    // ====================================================================
    println!("\n======================================================");
    println!("[2/2] Cenário 2: Variação da Repulsão Vetorial (gv)");
    println!("======================================================");

    let bag_fixa = 85.38;

    // DEFINIÇÃO DA ESCALA (Vai de gv=0.0 até gv=2.0)
    let gv_scale = ColorScale::new(0.0, 2.0, Palette::Plasma);

    let mut artist_mr_gv =
        Artist::new("results/gv_var_mr.svg", "MR - Varying Vector Coupling (gv)")
            .with_x_label("Radius [km]")
            .with_y_label("Mass [M\u{2299}]")
            .with_x_range(6.0, 13.0);

    let mut artist_eos_gv = Artist::new(
        "results/gv_var_eos.svg",
        "EoS - Varying Vector Coupling (gv)",
    )
    .with_x_label("\u{03B5} [MeV/fm\u{00B3}]")
    .with_y_label("P [MeV/fm\u{00B3}]")
    .autoscale()
    .with_log_scale();

    // Loop de 0 a 20 gerando (0.0, 0.1, 0.2 ... 2.0)
    for i in 0..=20 {
        let gv = (i as f64) * 0.1;
        println!("  -> Calculando gv = {:.1}...", gv);

        let mut quark_engine = QuarksMatter::new(bag_fixa, bg)
            .with_limits(1.0, 3.5)
            .with_points(1500);
        quark_engine.gv = gv;

        let mut solver_q = Solver::new(EngineMode::Quarks(quark_engine));
        let res_q = solver_q.solve();
        let (eps_q, p_q, rho_q) = extract_eos_filtered(&res_q);

        if eps_q.len() > 10 {
            let (m_q, r_q, _b_q, _) = generate_mr_curve(&eps_q, &p_q, &rho_q, true);
            let label = format!("gv = {:.1}", gv);

            // A mágica acontece aqui: pegamos a cor direto do ColorScale
            let (r, g, b) = gv_scale.get_color(gv);

            artist_mr_gv = artist_mr_gv.add_curve_color(&r_q, &m_q, &label, r, g, b);
            artist_eos_gv = artist_eos_gv.add_curve_color(&eps_q, &p_q, &label, r, g, b);
        }
    }

    artist_mr_gv.plot().ok();
    artist_eos_gv.plot().ok();
    println!("Gráficos do Cenário 2 gerados: gv_var_mr.svg e gv_var_eos.svg");
    println!("\nProcesso finalizado com sucesso!");
    Ok(())
}

/// `quarks hybrid`: GM1 (B = 1e15 G) x MIT bag x híbrida (Maxwell).
pub fn hybrid(_raw: &[String]) -> Result<(), String> {
    create_dir("results")?;
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
    Ok(())
}

/// Função auxiliar para extrair Eps e P dos resultados do solver.
/// Diferente do hybrid.rs, filtramos as energias e pressões > 0
/// para garantir que o TOV não falhe ao lidar com o vácuo autoligado das Strange Stars.
fn extract_eos_filtered(results: &[[f64; RESULTS_SIZE]]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut eps = Vec::new();
    let mut p = Vec::new();
    let mut rho = Vec::new();
    for r in results {
        // Apenas pressões e densidades estritamente positivas entram no TOV
        if r[1] > 0.0 && r[2] > 0.0 {
            eps.push(r[1]);
            p.push(r[2]);
            rho.push(r[0]);
        }
    }
    (eps, p, rho)
}

/// Função auxiliar para extrair Eps e P dos resultados do solver
fn extract_eos(results: &[[f64; RESULTS_SIZE]]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let eps = results.iter().map(|r| r[1]).collect();
    let p = results.iter().map(|r| r[2]).collect();
    let rho = results.iter().map(|r| r[0]).collect();
    (eps, p, rho)
}
