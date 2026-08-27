"""Tests for the exact noise-aware dual proximal-gradient solver."""

import numpy as np
import pytest

import HippsDimes
from hipps_dimes import proximal_fista as pf


def _rouse_squared_distances(n, spring_constant=1.0):
    indices = np.arange(n)
    return 3.0 * np.abs(indices[:, None] - indices[None, :]) / spring_constant


def test_exact_smooth_gradient_matches_finite_difference():
    n = 4
    target = _rouse_squared_distances(n, spring_constant=1.3)
    variance = 0.2
    (
        _,
        pair_i,
        pair_j,
        target_pairs,
        pair_variance,
        _,
        _,
    ) = pf._validate_inputs(target, variance, None)
    initial, _ = pf._prepare_initial_connectivity(
        n, pair_i, pair_j, target_pairs, None
    )
    fixed = initial.copy()
    np.fill_diagonal(fixed, 0.0)
    theta = -0.5 * initial[pair_i, pair_j]
    basis = pf._centered_orthonormal_basis(n)

    def evaluate(value):
        return pf._evaluate_state(
            value,
            fixed,
            basis,
            pair_i,
            pair_j,
            target_pairs,
            pair_variance,
            np,
        )

    state = evaluate(theta)
    assert state is not None
    direction = np.zeros_like(theta)
    direction[2] = 1.0
    epsilon = 1e-6
    plus = evaluate(theta + epsilon * direction)
    minus = evaluate(theta - epsilon * direction)
    assert plus is not None and minus is not None
    finite_difference = (
        plus["smooth_objective"] - minus["smooth_objective"]
    ) / (2.0 * epsilon)
    analytic = float(np.dot(state["smooth_gradient"], direction))
    assert finite_difference == pytest.approx(analytic, rel=2e-6, abs=2e-8)


def test_two_locus_solution_matches_analytic_stationary_point():
    observed_distance = 3.4
    variance = 0.6
    target = np.array([[0.0, observed_distance], [observed_distance, 0.0]])

    fitted, gram, connectivity, info = (
        HippsDimes.fit_gaussian_noise_dual_fista(
            target,
            noise_variance=variance,
            max_iterations=3000,
            relative_tolerance=1e-9,
            absolute_tolerance=1e-12,
        )
    )

    expected_spring = (
        -observed_distance
        + np.sqrt(observed_distance**2 + 6.0 * variance)
    ) / variance
    assert info["converged"]
    assert connectivity[0, 1] == pytest.approx(
        expected_spring, rel=2e-7, abs=2e-9
    )
    assert fitted[0, 1] - observed_distance == pytest.approx(
        0.5 * variance * connectivity[0, 1],
        rel=2e-7,
        abs=2e-9,
    )
    assert np.linalg.matrix_rank(gram, tol=1e-10) == 1
    assert info["relative_kkt_residual"] <= 1e-9


def test_monotone_fista_preserves_validity_and_decreases_objective():
    n = 6
    target = _rouse_squared_distances(n)
    target = target.copy()
    target[0, 5] = target[5, 0] = 0.35 * target[0, 5]
    target[1, 4] = target[4, 1] = 1.65 * target[1, 4]

    _, _, connectivity, info = HippsDimes.fit_gaussian_noise_dual_fista(
        target,
        noise_variance=0.4,
        max_iterations=2500,
        relative_tolerance=2e-6,
        absolute_tolerance=1e-10,
        accelerated=True,
        monotone=True,
    )

    objective = info["history"]["objective"]
    assert objective.size > 1
    assert np.all(
        np.diff(objective)
        <= 2e-10 * np.maximum(1.0, np.abs(objective[:-1]))
    )
    assert np.all(
        info["history"]["minimum_internal_precision_eigenvalue"] > 0.0
    )
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 2e-10
    assert info["relative_kkt_residual"] <= 2e-6


def test_heteroskedastic_variance_obeys_pair_stationarity():
    n = 5
    target = _rouse_squared_distances(n, spring_constant=1.4)
    pair_i, pair_j = np.triu_indices(n, k=1)
    variance = np.zeros((n, n), dtype=float)
    pair_values = 0.03 + 0.01 * (pair_j - pair_i)
    variance[pair_i, pair_j] = pair_values
    variance[pair_j, pair_i] = pair_values

    fitted, _, connectivity, info = (
        HippsDimes.fit_gaussian_noise_connectivity_fista(
            target,
            noise_variance=variance,
            max_iterations=3000,
            relative_tolerance=2e-7,
            absolute_tolerance=1e-11,
        )
    )

    assert info["converged"]
    assert info["noise_model"] == "heteroskedastic_variance_matrix"
    stationarity = (
        (fitted - target)[pair_i, pair_j]
        - 0.5
        * variance[pair_i, pair_j]
        * connectivity[pair_i, pair_j]
    )
    scale = max(
        1.0,
        np.linalg.norm((fitted - target)[pair_i, pair_j]),
        np.linalg.norm(
            0.5
            * variance[pair_i, pair_j]
            * connectivity[pair_i, pair_j]
        ),
    )
    assert np.linalg.norm(stationarity) / scale <= 2e-7
    assert np.max(np.abs(stationarity)) <= 5e-8


def test_accelerated_and_plain_proximal_gradient_reach_same_solution():
    n = 4
    target = _rouse_squared_distances(n, spring_constant=0.9)
    target[0, 3] *= 0.8
    target[3, 0] = target[0, 3]

    accelerated = HippsDimes.fit_gaussian_noise_dual_fista(
        target,
        noise_variance=0.2,
        max_iterations=2500,
        relative_tolerance=1e-6,
        accelerated=True,
    )
    plain = HippsDimes.fit_gaussian_noise_dual_fista(
        target,
        noise_variance=0.2,
        max_iterations=5000,
        relative_tolerance=1e-6,
        accelerated=False,
    )

    assert accelerated[3]["relative_kkt_residual"] <= 1e-6
    assert plain[3]["relative_kkt_residual"] <= 1e-6
    assert accelerated[3]["objective"] == pytest.approx(
        plain[3]["objective"], rel=2e-7, abs=2e-8
    )
    assert np.allclose(accelerated[2], plain[2], rtol=2e-5, atol=2e-7)
