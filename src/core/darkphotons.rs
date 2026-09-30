// src/core/darkphotons.rs
//
// Setor escuro fermiônico: férmion de Dirac chi e fóton escuro X com mistura
// cinética. A física do setor está em `HadronsMatter` (builders `with_y_chi`,
// `with_m_chi`, `with_m_x`, `with_g_d`, `with_epsilon`); este módulo mantém
// a termodinâmica do gás de Dirac escuro e os testes do setor.

use crate::core::constants::PI2;

/// Engine com setor escuro: é a mesma `HadronsMatter`, com os parâmetros
/// escuros definidos pelos builders. Mantido para compatibilidade.
pub type DarkPhotonsMatter = crate::core::physics::HadronsMatter;

/// Number density of a zero-temperature spin-1/2 Dirac gas in internal natural units.
pub fn dark_fermion_number_density_from_kf(kf: f64) -> f64 {
    if kf <= 0.0 {
        0.0
    } else {
        kf.powi(3) / (3.0 * PI2)
    }
}

/// Fermi momentum of a zero-temperature spin-1/2 Dirac gas in internal natural units.
pub fn dark_fermion_kf_from_density(density: f64) -> f64 {
    if density <= 0.0 {
        0.0
    } else {
        (3.0 * PI2 * density).cbrt()
    }
}

/// Kinetic (including rest mass) energy density of a free Dirac gas.
pub fn dark_fermion_energy_density(kf: f64, mass: f64) -> f64 {
    assert!(
        mass.is_finite() && mass > 0.0,
        "dark fermion mass must be finite and positive"
    );
    if kf <= 0.0 {
        return 0.0;
    }

    let x = kf / mass;
    let dimensionless = if x.abs() < 1e-2 {
        // Series avoids cancellation between the algebraic and asinh terms.
        (8.0 / 3.0) * x.powi(3) + (4.0 / 5.0) * x.powi(5) - x.powi(7) / 7.0 + x.powi(9) / 18.0
    } else {
        let ef_over_m = x.hypot(1.0);
        x * ef_over_m * (2.0 * x.powi(2) + 1.0) - x.asinh()
    };

    mass.powi(4) * dimensionless / (8.0 * PI2)
}

/// Kinetic pressure of a zero-temperature spin-1/2 Dirac gas.
pub fn dark_fermion_pressure(kf: f64, mass: f64) -> f64 {
    assert!(
        mass.is_finite() && mass > 0.0,
        "dark fermion mass must be finite and positive"
    );
    if kf <= 0.0 {
        return 0.0;
    }

    let x = kf / mass;
    let dimensionless = if x.abs() < 1e-2 {
        (8.0 / 5.0) * x.powi(5) - (4.0 / 7.0) * x.powi(7) + x.powi(9) / 3.0
    } else {
        let ef_over_m = x.hypot(1.0);
        x * ef_over_m * (2.0 * x.powi(2) - 3.0) + 3.0 * x.asinh()
    };

    mass.powi(4) * dimensionless / (24.0 * PI2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::constants::{HBAR_C, M_NUCLEON, RESULTS_SIZE};
    use crate::core::eos::compute;
    use crate::core::model::{FSU2, GM1, GM3};
    use crate::core::physics::HadronsMatter;
    use crate::core::solver::{EngineMode, Solver};

    type DarkStateVector = nalgebra::SVector<f64, 5>;

    const TEST_MUN: f64 = 1.15;
    const TEST_Y_CHI: f64 = 0.02;
    const TEST_M_CHI: f64 = 0.8;
    const TEST_M_X: f64 = 0.4;

    fn assert_close(left: f64, right: f64, rel: f64, abs: f64) {
        let error = (left - right).abs();
        let scale = left.abs().max(right.abs());
        assert!(
            error <= abs.max(rel * scale),
            "left={left:.16e}, right={right:.16e}, error={error:.3e}"
        );
    }

    fn solve(engine: &mut DarkPhotonsMatter, mun: f64) -> ([f64; 5], [f64; RESULTS_SIZE]) {
        engine
            .solve_point(mun, &[0.0; 5])
            .unwrap_or_else(|| panic!("dark-sector point failed to converge at mun={mun}"))
    }

    fn visible_state(engine: &DarkPhotonsMatter) -> Vec<f64> {
        engine
            .nb
            .iter()
            .chain(engine.nl.iter())
            .copied()
            .chain([engine.mue])
            .collect()
    }

    fn assert_same_visible_state(left: &DarkPhotonsMatter, right: &DarkPhotonsMatter) {
        for (a, b) in visible_state(left).iter().zip(visible_state(right).iter()) {
            assert_close(*a, *b, 2e-8, 2e-11);
        }
    }

    #[test]
    fn dirac_gas_density_and_thermodynamic_identity() {
        for ratio in [0.0, 1e-5, 9.9e-3, 1.01e-2, 0.1, 1.0, 10.0, 100.0] {
            let mass = 0.73;
            let kf = ratio * mass;
            let density = dark_fermion_number_density_from_kf(kf);
            let energy = dark_fermion_energy_density(kf, mass);
            let pressure = dark_fermion_pressure(kf, mass);
            let ef = kf.hypot(mass);

            assert_close(dark_fermion_kf_from_density(density), kf, 2e-13, 1e-18);
            assert_close(energy + pressure, density * ef, 5e-10, 1e-27);
            assert!(energy >= 0.0 && pressure >= 0.0);
        }
    }

    #[test]
    fn builders_validate_and_convert_dark_parameters() {
        let engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi_mev(500.0)
            .with_m_x_mev(25.0)
            .with_g_d(0.4)
            .with_epsilon(1e-3)
            .with_y_chi(0.03);

        assert_close(engine.m_chi, 500.0 / M_NUCLEON, 0.0, f64::EPSILON);
        assert_close(engine.m_x, 25.0 / M_NUCLEON, 0.0, f64::EPSILON);
        assert_eq!(engine.g_d, 0.4);
        assert_eq!(engine.epsilon, 1e-3);
        assert_eq!(engine.y_chi, 0.03);

        assert!(
            std::panic::catch_unwind(|| { DarkPhotonsMatter::new(GM1, 0.0).with_m_x(f64::NAN) })
                .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                DarkPhotonsMatter::new(GM1, 0.0).with_g_d(f64::INFINITY)
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| { DarkPhotonsMatter::new(GM1, 0.0).with_epsilon(1.0) })
                .is_err()
        );
    }

    #[test]
    #[should_panic(expected = "dark fermion mass")]
    fn builder_rejects_nonpositive_dark_fermion_mass() {
        let _ = DarkPhotonsMatter::new(GM1, 0.0).with_m_chi(0.0);
    }

    #[test]
    #[should_panic(expected = "dark number fraction")]
    fn builder_rejects_negative_dark_fraction() {
        let _ = DarkPhotonsMatter::new(GM1, 0.0).with_y_chi(-1e-3);
    }

    #[test]
    fn zero_fraction_recovers_visible_dark_engine() {
        for model in [GM1, GM3, FSU2] {
            let mut baseline = HadronsMatter::new(model, 0.0);
            let mut zero_fraction = DarkPhotonsMatter::new(model, 0.0)
                .with_m_chi(TEST_M_CHI)
                .with_m_x(TEST_M_X)
                .with_g_d(0.8)
                .with_epsilon(0.2)
                .with_y_chi(0.0);

            let mut visible_guess = [0.0; 5];
            let mut dark_guess = [0.0; 5];
            for mun in [1.12, 1.15, 1.18] {
                let (next_visible, baseline_result) = baseline
                    .solve_point(mun, &visible_guess)
                    .expect("hadronic reference point must converge");
                let (next_dark, result) = zero_fraction
                    .solve_point(mun, &dark_guess)
                    .expect("zero-fraction dark point must converge");
                visible_guess = next_visible;
                dark_guess = next_dark;

                for i in 0..=20 {
                    assert_close(result[i], baseline_result[i], 2e-8, 2e-9);
                }
                assert_eq!(zero_fraction.n_chi, 0.0);
                assert_close(zero_fraction.v_x0, 0.0, 0.0, 1e-15);
                assert_eq!(zero_fraction.ener_chi_kin, 0.0);
                assert_eq!(zero_fraction.press_chi_kin, 0.0);
            }
        }
    }

    #[test]
    fn free_dark_gas_changes_only_kinetic_thermodynamics() {
        let mut baseline = DarkPhotonsMatter::new(GM1, 0.0);
        let (_, baseline_result) = solve(&mut baseline, TEST_MUN);

        let mut dark = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_y_chi(TEST_Y_CHI);
        let (_, result) = solve(&mut dark, TEST_MUN);

        assert_same_visible_state(&baseline, &dark);
        assert_close(dark.v_x0, 0.0, 0.0, 1e-14);
        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        assert_close(
            result[1] - baseline_result[1],
            dark.ener_chi_kin * conversion,
            2e-8,
            2e-8,
        );
        assert_close(
            result[2] - baseline_result[2],
            dark.press_chi_kin * conversion,
            2e-8,
            2e-8,
        );
    }

    #[test]
    fn dark_self_interaction_matches_analytic_proca_solution() {
        let mut free = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_y_chi(TEST_Y_CHI);
        let (_, free_result) = solve(&mut free, TEST_MUN);

        let mut interacting = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let (_, result) = solve(&mut interacting, TEST_MUN);

        assert_same_visible_state(&free, &interacting);
        let expected_x0 = interacting.g_d * interacting.n_chi / interacting.m_x.powi(2);
        assert_close(interacting.v_x0, expected_x0, 2e-9, 2e-11);

        let vector_energy = interacting.dark_vector_energy_density();
        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        assert_close(
            result[1] - free_result[1],
            vector_energy * conversion,
            2e-8,
            2e-8,
        );
        assert_close(
            result[2] - free_result[2],
            vector_energy * conversion,
            2e-8,
            2e-8,
        );

        let vector_gibbs = interacting.g_d * interacting.n_chi * interacting.v_x0;
        assert_close(vector_gibbs, 2.0 * vector_energy, 2e-9, 2e-12);
        assert_close(
            interacting.mu_chi * interacting.n_chi - interacting.ener_chi_kin - vector_energy,
            interacting.press_chi_kin + vector_energy,
            2e-9,
            2e-12,
        );
    }

    #[test]
    fn kinetic_mixing_satisfies_neutrality_proca_and_sign_convention() {
        let epsilon = 0.08;
        let mut no_portal = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let _ = solve(&mut no_portal, TEST_MUN);

        let mut portal = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(epsilon)
            .with_y_chi(TEST_Y_CHI);
        let (solution, _) = solve(&mut portal, TEST_MUN);

        let charge = portal.charge_neutrality();
        assert!(charge.abs() < 1e-10, "charge residual={charge:.3e}");
        assert!(
            portal.dark_photon_residual(charge).abs() < 1e-10,
            "Proca residual={:.3e}",
            portal.dark_photon_residual(charge)
        );
        let norm = DarkPhotonsMatter::kinetic_mixing_norm(epsilon);
        let expected_x0 = portal.g_d * portal.n_chi / (portal.m_x.powi(2) * norm);
        assert_close(portal.v_x0, expected_x0, 2e-8, 2e-10);
        assert_close(solution[4], portal.v_x0, 0.0, 1e-15);
        assert_close(portal.n_chi, portal.y_chi * portal.nbt, 1e-14, 1e-15);
        assert_close(
            portal.mu_chi,
            portal.ef_chi + portal.g_d * portal.v_x0 / norm,
            1e-14,
            1e-14,
        );
        assert_close(
            portal.g_d * portal.n_chi * portal.v_x0 / norm,
            2.0 * portal.dark_vector_energy_density(),
            2e-8,
            2e-11,
        );

        let delta_x = epsilon * portal.qe * portal.v_x0 / norm;
        assert_close(portal.ef_l[0], portal.mue + delta_x, 2e-12, 2e-13);
        assert_close(portal.mue + delta_x, no_portal.mue, 3e-5, 5e-8);
        for (left, right) in portal
            .nb
            .iter()
            .chain(portal.nl.iter())
            .zip(no_portal.nb.iter().chain(no_portal.nl.iter()))
        {
            assert_close(*left, *right, 3e-5, 5e-9);
        }
    }

    #[test]
    fn total_pressure_obeys_gibbs_relation() {
        let mut baseline = HadronsMatter::new(GM1, 0.0);
        let (_, baseline_result) = baseline
            .solve_point(TEST_MUN, &[0.0; 4])
            .expect("hadronic reference point must converge");

        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.04)
            .with_y_chi(TEST_Y_CHI);
        let (solution, result) = solve(&mut engine, TEST_MUN);
        let (energy, pressure) =
            compute(&engine, solution[0], solution[1], solution[2], solution[3]);
        let gibbs = engine
            .mu_b
            .iter()
            .zip(engine.nb.iter())
            .map(|(mu, n)| mu * n)
            .sum::<f64>()
            + engine.mue * engine.nl.iter().sum::<f64>()
            + engine.mu_chi * engine.n_chi
            - energy;
        assert_close(pressure, gibbs, 1e-13, 1e-14);

        let conversion = M_NUCLEON * (M_NUCLEON / HBAR_C).powi(3);
        let expected_dark_pressure =
            (engine.press_chi_kin + engine.dark_vector_energy_density()) * conversion;
        assert_close(
            result[2] - baseline_result[2],
            expected_dark_pressure,
            5e-6,
            5e-7,
        );
    }

    #[test]
    fn proca_residual_includes_off_shell_visible_charge_source() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.2)
            .with_y_chi(TEST_Y_CHI);
        engine.n_chi = 3e-4;
        engine.v_x0 = 0.07;
        let n_q = -2.5e-3;
        let norm = DarkPhotonsMatter::kinetic_mixing_norm(engine.epsilon);
        let expected = engine.m_x.powi(2) * engine.v_x0
            - (engine.g_d * engine.n_chi + engine.epsilon * engine.qe * n_q) / norm;

        assert_close(engine.dark_photon_residual(n_q), expected, 1e-14, 1e-14);
        assert_ne!(
            engine.dark_photon_residual(n_q),
            engine.dark_photon_residual(0.0)
        );
    }

    #[test]
    fn dark_sector_has_vacuum_limit() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_epsilon(0.05)
            .with_y_chi(TEST_Y_CHI);
        engine.mun = 0.1;
        let residual = engine.funcv(&[0.0; 5]);

        assert!(residual.iter().all(|value| value.abs() < 1e-15));
        assert_eq!(engine.nbt, 0.0);
        assert_eq!(engine.n_chi, 0.0);
        assert_eq!(engine.v_x0, 0.0);
        assert_eq!(engine.ener_chi_kin, 0.0);
        assert_eq!(engine.press_chi_kin, 0.0);

        let mut previous_energy = f64::INFINITY;
        let mut previous_pressure = f64::INFINITY;
        for n_b in [1e-4, 1e-6, 1e-8, 1e-10] {
            engine.n_chi = engine.y_chi * n_b;
            engine.v_x0 = engine.g_d * engine.n_chi
                / (engine.m_x.powi(2) * DarkPhotonsMatter::kinetic_mixing_norm(engine.epsilon));
            engine.update_dark_fermion_state();
            let energy = engine.ener_chi_kin + engine.dark_vector_energy_density();
            let pressure = engine.press_chi_kin + engine.dark_vector_energy_density();
            assert_close(engine.n_chi / n_b, engine.y_chi, 2e-14, 1e-15);
            assert!(energy < previous_energy);
            assert!(pressure < previous_pressure);
            previous_energy = energy;
            previous_pressure = pressure;
        }
        assert!(previous_energy < 1e-10);
        assert!(previous_pressure < 1e-15);
    }

    #[test]
    fn linear_gm1_vector_residual_matches_documented_normalization() {
        // Locks the linear GM1 normalization of the omega equation.
        let mut gm1 = DarkPhotonsMatter::new(GM1, 0.0);
        gm1.nb[0] = 0.01;
        let vomega = 0.12;
        let derivative = vomega / GM1.gv.powi(2);
        assert_close(
            gm1.equation_omega(vomega, 0.0),
            GM1.gv.powi(2) * (gm1.nb[0] - derivative),
            1e-14,
            1e-14,
        );
        assert_ne!(FSU2.rxi, 0.0);
        assert_ne!(FSU2.lambda_v, 0.0);
    }

    #[test]
    fn eos_columns_use_documented_physical_units() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let (solution, result) = solve(&mut engine, TEST_MUN);
        let density_factor = (M_NUCLEON / HBAR_C).powi(3);
        let energy_factor = M_NUCLEON * density_factor;

        assert_eq!(result.len(), 34);
        assert_close(result[3], engine.nl[0] * density_factor, 1e-14, 1e-14);
        assert_close(result[5], engine.nb[0] * density_factor, 1e-14, 1e-14);
        assert_close(result[13], solution[1] * M_NUCLEON, 1e-14, 1e-14);
        assert_close(result[16], engine.m_eff[0], 1e-14, 1e-14);
        assert_close(result[21], engine.n_chi * density_factor, 1e-14, 1e-14);
        assert_close(result[23], engine.m_chi * M_NUCLEON, 1e-14, 1e-14);
        assert_close(result[27], engine.v_x0 * M_NUCLEON, 1e-14, 1e-14);
        assert_close(
            result[30],
            engine.ener_chi_kin * energy_factor,
            1e-14,
            1e-14,
        );
        assert_close(result[32], result[33], 0.0, 0.0);
    }

    #[test]
    fn previous_dark_solution_is_a_valid_five_dimensional_guess() {
        let mut engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(0.01)
            .with_g_d(0.8)
            .with_epsilon(1e-3)
            .with_y_chi(0.05);
        let (first, _) = solve(&mut engine, 1.12);
        assert_ne!(first[4], 0.0);

        let (second, _) = engine
            .solve_point(1.14, &first)
            .expect("next EOS point must accept the previous five-dimensional solution");
        let residual = engine.funcv(&second);
        assert!(
            DarkStateVector::from_column_slice(&residual).norm() < 1e-10,
            "residual at continued point is {:?}",
            residual
        );
    }

    #[test]
    fn generic_solver_carries_the_dark_field_between_eos_points() {
        let engine = DarkPhotonsMatter::new(GM1, 0.0)
            .with_limits(1.12, 1.14)
            .with_points(3)
            .with_m_chi(TEST_M_CHI)
            .with_m_x(TEST_M_X)
            .with_g_d(0.35)
            .with_y_chi(TEST_Y_CHI);
        let mut solver = Solver::new(EngineMode::DarkPhotons(engine));
        let rows = solver.solve();

        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|row| row[21] > 0.0));
        assert!(rows.iter().all(|row| row[27] > 0.0));
    }

    #[test]
    fn zero_background_has_no_magnetic_vacuum_floor() {
        for model in [GM1, GM3] {
            let mut visible = HadronsMatter::new(model, 0.0);
            let (_, visible_row) = visible
                .solve_point(0.99, &[0.0; 4])
                .expect("visible vacuum point must converge");

            let mut dark = DarkPhotonsMatter::new(model, 0.0)
                .with_m_chi(TEST_M_CHI)
                .with_m_x(TEST_M_X)
                .with_g_d(0.35)
                .with_epsilon(1e-4)
                .with_y_chi(TEST_Y_CHI);
            let (_, dark_row) = dark
                .solve_point(0.99, &[0.0; 5])
                .expect("dark vacuum point must converge");

            for row in [visible_row, dark_row] {
                assert_eq!(row[0], 0.0);
                assert_eq!(row[1], 0.0);
                assert_eq!(row[2], 0.0);
                assert_eq!(row[19], 0.0);
            }
        }
    }

    #[test]
    fn dark_continuation_crosses_the_vacuum_matter_onset() {
        for model in [GM1, GM3] {
            let engine = DarkPhotonsMatter::new(model, 0.0)
                .with_limits(1.0, 1.02)
                .with_points(41)
                .with_m_chi_mev(200_000.0)
                .with_m_x_mev(100.0)
                .with_g_d(0.45)
                .with_epsilon(1e-4)
                .with_y_chi(2.298_360_48e-4);
            let rows = Solver::new(EngineMode::DarkPhotons(engine)).solve();

            assert!(
                rows.iter().any(|row| row[0] > 1e-4),
                "dark continuation remained trapped in the {model:?} vacuum branch"
            );
        }
    }
}
