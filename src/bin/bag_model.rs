// src/bin/bag_model.rs

use nsrs::core::constants::RESULTS_SIZE;
use nsrs::core::quarks::QuarksMatter;
use nsrs::core::solver::{Solver, EngineMode};
use nsrs::core::tov_solver::generate_mr_curve;
use nsrs::core::plotting::{Artist, ColorScale, Palette};

fn main() {
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
    
    let mut artist_mr_gv = Artist::new("results/gv_var_mr.svg", "MR - Varying Vector Coupling (gv)")
        .with_x_label("Radius [km]")
        .with_y_label("Mass [M\u{2299}]")
        .with_x_range(6.0, 13.0);
        
    let mut artist_eos_gv = Artist::new("results/gv_var_eos.svg", "EoS - Varying Vector Coupling (gv)")
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