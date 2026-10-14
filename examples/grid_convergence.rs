// Teste de convergência da malha em mu_n: uniforme x adaptativa.
use nsrs::core::observations::{interpolate_at_mass, stable_branch};
use nsrs::core::tov_solver::generate_star_sequence;
use nsrs::{EngineMode, FSU2, FieldProfile, GM1, GM3, HadronsMatter, Solver};
use rayon::prelude::*;

fn main() {
    let models = [("GM1", GM1), ("GM3", GM3), ("FSU2", FSU2)];
    let fields = [(0.0, false), (1e17, false), (1e18, false), (1e18, true)];
    let grids: Vec<(String, usize, Option<f64>)> = vec![
        ("u6000".into(), 6000, None),
        ("u1500".into(), 1500, None),
        ("a3e-3".into(), 1500, Some(3e-3)),
        ("a1e-3".into(), 1500, Some(1e-3)),
        ("a3e-4".into(), 1500, Some(3e-4)),
    ];
    let mut jobs = Vec::new();
    for (mname, m) in models { for (b, bdd) in fields { for amm in [false, true] { for g in &grids {
        jobs.push((mname, m, b, bdd, amm, g.clone()));
    }}}}
    let out: Vec<String> = jobs.par_iter().with_max_len(1).map(|(mname, m, b, bdd, amm, (gname, n, tol))| {
        let mut e = HadronsMatter::new(*m, *b).with_limits(0.02, 3.0).with_points(*n).with_anomalous_moments(*amm);
        if *bdd { e = e.with_field_profile(FieldProfile::bdd(*b)); }
        if let Some(t) = tol { e = e.with_step_tolerance(*t); }
        let t0 = std::time::Instant::now();
        let mut s = Solver::new(EngineMode::Hadrons(e));
        let rows = s.solve();
        let t_eos = t0.elapsed().as_secs_f64();
        let eps: Vec<f64> = rows.iter().map(|r| r[1]).collect();
        let p: Vec<f64> = rows.iter().map(|r| r[2]).collect();
        let nb: Vec<f64> = rows.iter().map(|r| r[0]).collect();
        let t1 = std::time::Instant::now();
        let stars = generate_star_sequence(&eps, &p, &nb, true);
        let t_tov = t1.elapsed().as_secs_f64();
        let br = stable_branch(&stars);
        let mmax = br.last().map_or(f64::NAN, |s| s.mass);
        let r14 = interpolate_at_mass(br, 1.4, |s| s.radius).unwrap_or(f64::NAN);
        let l14 = interpolate_at_mass(br, 1.4, |s| s.tidal_deformability).unwrap_or(f64::NAN);
        let matter = rows.iter().filter(|r| r[0] > 0.0).count();
        format!("{mname},{b:e},{},{amm},{gname},{},{matter},{t_eos:.3},{t_tov:.3},{mmax:.6},{r14:.5},{l14:.3},{:?}",
            if *bdd {"bdd"} else {"const"}, rows.len(), s.termination().map(|t| t.is_anomalous()))
    }).collect();
    println!("model,B,profile,amm,grid,rows,matter_rows,t_eos_s,t_tov_s,m_max,r14,lambda14,anomalous");
    for l in out { println!("{l}"); }
}
