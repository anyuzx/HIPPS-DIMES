"""Tests for the globalized noise-aware covariance-cone solver."""

import numpy as np
import pytest

import hipps_dimes.covariance_geometry as geometry
import hipps_dimes.covariance_solver as solver
import hipps_dimes.numerics as numerics


def _rouse_squared_distance_target(n=7, spring_constant=1.7):
    connectivity = numerics.construct_connectivity_matrix_rouse(
        n, spring_constant
    )
    basis = numerics._centered_orthonormal_basis(n)
    reduced_gram, _ = numerics._reduced_gram_from_connectivity(
        connectivity, basis
    )
    full_gram = basis @ reduced_gram @ basis.T
    return numerics._squared_distances_from_gram(full_gram)


def _symmetric_noise(n, scale, seed=20260826):
    rng = np.random.default_rng(seed)
    noise = rng.normal(0.0, scale, size=(n, n))
    noise = 0.5 * (noise + noise.T)
    np.fill_diagonal(noise, 0.0)
    return noise


def test_globalized_covariance_solver_is_installed():
    solver.install_covariance_solver()
    assert (
        numerics.fit_gaussian_noise_covariance
        is solver.fit_gaussian_noise_covariance
    )
    assert hasattr(numerics, "_legacy_fit_gaussian_noise_covariance")


def test_anchored_data_gradient_and_hessian_match_finite_differences():
    rng = np.random.default_rng(7)
    n_modes = 4
    n = n_modes + 1
    gram = np.diag(np.linspace(0.8, 1.4, n_modes))
    gram += 0.03 * np.ones_like(gram)
    target = np.zeros((n, n))
    weights = np.zeros((n, n))
    pair_i, pair_j = np.triu_indices(n, k=1)
    target_values = np.linspace(0.7, 2.0, pair_i.size)
    pair_weights = np.linspace(1.0, 3.0, pair_i.size)
    target[pair_i, pair_j] = target_values
    target[pair_j, pair_i] = target_values
    weights[pair_i, pair_j] = pair_weights
    weights[pair_j, pair_i] = pair_weights

    objective, gradient, _ = geometry.anchored_data_objective_gradient(
        gram,
        target,
        weights,
        pair_i,
        pair_j,
        array_module=np,
    )
    assert np.isfinite(objective)
    direction = rng.normal(size=(n_modes, n_modes))
    direction = 0.5 * (direction + direction.T)
    epsilon = 1e-6
    forward = geometry.anchored_data_objective_gradient(
        gram + epsilon * direction,
        target,
        weights,
        pair_i,
        pair_j,
        array_module=np,
    )[0]
    backward = geometry.anchored_data_objective_gradient(
        gram - epsilon * direction,
        target,
        weights,
        pair_i,
        pair_j,
        array_module=np,
    )[0]
    finite_difference = (forward - backward) / (2.0 * epsilon)
    assert finite_difference == pytest.approx(
        np.sum(gradient * direction), rel=1e-8, abs=1e-8
    )

    analytic_action = geometry.anchored_data_hessian_action(
        direction, weights, array_module=np
    )
    forward_gradient = geometry.anchored_data_objective_gradient(
        gram + epsilon * direction,
        target,
        weights,
        pair_i,
        pair_j,
        array_module=np,
    )[1]
    backward_gradient = geometry.anchored_data_objective_gradient(
        gram - epsilon * direction,
        target,
        weights,
        pair_i,
        pair_j,
        array_module=np,
    )[1]
    finite_difference_action = (
        forward_gradient - backward_gradient
    ) / (2.0 * epsilon)
    assert np.allclose(
        analytic_action,
        finite_difference_action,
        rtol=1e-8,
        atol=1e-8,
    )


def test_svec_preconditioner_diagonal_matches_exact_coordinate_curvatures():
    n_modes = 4
    n = n_modes + 1
    gram = np.diag(np.linspace(0.7, 1.5, n_modes))
    gram += 0.04 * np.ones_like(gram)
    inverse_gram = np.linalg.inv(gram)
    weights = np.zeros((n, n))
    pair_i, pair_j = np.triu_indices(n, k=1)
    pair_weights = np.linspace(0.5, 2.5, pair_i.size)
    weights[pair_i, pair_j] = pair_weights
    weights[pair_j, pair_i] = pair_weights

    preconditioner = geometry.anchored_data_svec_diagonal(
        weights, np
    ) + geometry.entropy_svec_diagonal(inverse_gram, np)

    for i in range(n_modes):
        for j in range(i, n_modes):
            direction = np.zeros((n_modes, n_modes))
            if i == j:
                direction[i, i] = 1.0
            else:
                direction[i, j] = 1.0 / np.sqrt(2.0)
                direction[j, i] = 1.0 / np.sqrt(2.0)
            action = geometry.anchored_data_hessian_action(
                direction, weights, array_module=np
            )
            action += 1.5 * inverse_gram @ direction @ inverse_gram
            curvature = np.sum(direction * action)
            assert preconditioner[i, j] == pytest.approx(
                curvature, rel=1e-12, abs=1e-12
            )


def test_anchored_and_centered_parameterizations_reach_same_optimum():
    target = _rouse_squared_distance_target(n=7, spring_constant=1.8)
    target = target + _symmetric_noise(len(target), 0.015)
    variance = 0.0004
    options = {
        "max_iterations": 80,
        "relative_tolerance": 1e-8,
        "proximal_warmup_iterations": 5,
        "continuation_factors": (0.1, 1.0),
    }
    anchored = solver.fit_gaussian_noise_covariance(
        target,
        variance,
        coordinate_parameterization="anchored",
        **options,
    )
    centered = solver.fit_gaussian_noise_covariance(
        target,
        variance,
        coordinate_parameterization="centered",
        **options,
    )

    assert anchored[3]["converged"]
    assert centered[3]["converged"]
    assert anchored[3]["coordinate_parameterization"] == "anchored"
    assert centered[3]["coordinate_parameterization"] == "centered"
    assert np.allclose(anchored[0], centered[0], rtol=1e-8, atol=1e-9)
    assert np.allclose(anchored[2], centered[2], rtol=1e-8, atol=1e-9)
    assert anchored[3]["objective"] == pytest.approx(
        centered[3]["objective"], rel=1e-10, abs=1e-10
    )
    assert np.max(np.abs(np.sum(anchored[1], axis=1))) <= 1e-11
    assert np.max(np.abs(np.sum(anchored[2], axis=1))) <= 1e-11


def test_globalization_records_continuation_cg_and_line_search_diagnostics():
    target = _rouse_squared_distance_target(n=7, spring_constant=1.5)
    target = target + _symmetric_noise(len(target), 0.025, seed=11)
    _, gram, connectivity, info = solver.fit_gaussian_noise_covariance(
        target,
        0.0009,
        max_iterations=100,
        relative_tolerance=1e-7,
        continuation_factors=(0.1, 1.0),
        continuation_activation_relative_gradient=1e-12,
        proximal_warmup_iterations=5,
    )

    assert info["converged"]
    assert info["continuation"]["schedule"] == (0.1, 1.0)
    assert info["continuation"]["final_stage_reached"]
    history = info["history"]
    assert set(history["phase"]) <= {"proximal_warmup", "newton"}
    for values in history.values():
        assert len(values) == info["iterations"]
    assert np.all(history["line_search_backtracks"] >= 0)
    assert np.any(history["phase"] == "newton")
    newton = history["phase"] == "newton"
    assert np.all(np.isfinite(history["cg_forcing_tolerance"][newton]))
    assert np.all(
        history["cg_forcing_tolerance"][newton]
        >= info["globalization"]["cg_minimum_tolerance"]
    )
    assert np.all(
        history["cg_forcing_tolerance"][newton]
        <= info["globalization"]["cg_maximum_forcing_tolerance"]
    )
    assert info["preconditioner_kind"] == "svec_diagonal"
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-11
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-11


def test_two_locus_solution_matches_analytic_covariance_optimum():
    target_distance = 2.0
    variance = 0.3
    target = np.array(
        [[0.0, target_distance], [target_distance, 0.0]]
    )
    fitted, gram, connectivity, info = solver.fit_gaussian_noise_covariance(
        target,
        variance,
        max_iterations=20,
        relative_tolerance=1e-11,
    )
    expected_distance = 0.5 * (
        target_distance
        + np.sqrt(target_distance**2 + 6.0 * variance)
    )
    assert info["converged"]
    assert fitted[0, 1] == pytest.approx(expected_distance, rel=1e-12)
    assert info["relative_stationarity_residual"] <= 1e-12
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-12
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
