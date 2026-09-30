// src/bin/properties.rs
//
// Relatório de propriedades nucleares e estelares dos modelos (B = 0):
// saturação (n0, E/A, K, J, L, M*/M), sequência estelar com crosta
// (M_max, R_1.4, Lambda_1.4, I_1.4, z_1.4), limiar de URCA direto, início
// dos hyperons e c_s^2 máximo no ramo estável.
//
// Uso: cargo run --release --bin properties

use nsrs::core::io_utils::derived_diagnostics;
use nsrs::core::nuclear::saturation_properties;
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, FSU2, GM1, GM3, HadronsMatter, Solver};

const HYPERONS: [&str; 6] = ["Lambda", "Sigma-", "Sigma0", "Sigma+", "Xi-", "Xi0"];

fn main() {
    for (name, model) in [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)] {
        println!("== {name}");
        match saturation_properties(model) {
            Some(s) => println!(
                "  saturação: n0 = {:.4} fm^-3, E/A = {:.2} MeV, K = {:.1} MeV, J = {:.2} MeV, L = {:.1} MeV, M*/M = {:.3}",
                s.n0, s.energy_per_nucleon, s.incompressibility, s.symmetry_energy, s.symmetry_slope, s.effective_mass
            ),
            None => println!("  saturação: não encontrada"),
        }

        let rows = Solver::new(EngineMode::Hadrons(HadronsMatter::new(model, 0.0))).solve();
        let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
        let n: Vec<f64> = rows.iter().map(|r| r[0]).collect();
        let stars = generate_star_sequence(&eps, &p, &n, true);
        let Some(i_max) = (0..stars.len()).max_by(|&a, &b| stars[a].mass.total_cmp(&stars[b].mass))
        else {
            println!("  sem sequência estelar");
            continue;
        };
        let max = stars[i_max];
        println!(
            "  estrelas (com crosta): M_max = {:.3} Msun, R(M_max) = {:.2} km, P_c = {:.1} MeV/fm^3",
            max.mass, max.radius, max.central_pressure
        );
        if let Some(k) = (1..=i_max).find(|&k| stars[k].mass >= 1.4) {
            let (a, b) = (stars[k - 1], stars[k]);
            let t = (1.4 - a.mass) / (b.mass - a.mass);
            let lerp = |x: f64, y: f64| x + t * (y - x);
            println!(
                "  1.4 Msun: R = {:.2} km, Lambda = {:.0}, k2 = {:.4}, I = {:.3}e45 g cm^2, z = {:.3}",
                lerp(a.radius, b.radius),
                lerp(a.tidal_deformability, b.tidal_deformability),
                lerp(a.love_k2, b.love_k2),
                lerp(a.moment_of_inertia, b.moment_of_inertia),
                lerp(a.redshift, b.redshift)
            );
        }

        let derived = derived_diagnostics(&rows);
        let onset = |pred: &dyn Fn(usize) -> bool| (0..rows.len()).find(|&i| rows[i][0] > 1e-3 && pred(i));
        match onset(&|i| derived[i].direct_urca_electron) {
            Some(i) => println!(
                "  URCA direto (e): n_B = {:.2} n0, Y_p = {:.3}",
                rows[i][0], derived[i].proton_fraction
            ),
            None => println!("  URCA direto (e): não ocorre"),
        }
        let onsets: Vec<String> = HYPERONS
            .iter()
            .enumerate()
            .filter_map(|(h, label)| {
                onset(&|i| rows[i][7 + h] > 1e-6).map(|i| format!("{label} {:.2}", rows[i][0]))
            })
            .collect();
        println!("  início dos hyperons (n_B/n0): {}", onsets.join(", "));
        let p_c_max = max.central_pressure;
        let cs2_max = rows
            .iter()
            .zip(&derived)
            .filter(|(r, _)| r[2] <= p_c_max)
            .map(|(_, d)| d.sound_speed_squared)
            .fold(0.0, f64::max);
        println!("  c_s^2 máximo até o centro de M_max: {cs2_max:.3}");
    }
}
