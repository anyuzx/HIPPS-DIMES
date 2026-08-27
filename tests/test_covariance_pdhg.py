"""Tests for the primal-dual covariance-cone optimizer."""

import numpy as np
import pytest

import HippsDimes
from hipps_dimes import covariance_pdhg as pdhg


def _rouse_squared_distance_target(n=6, spring_constant=1.0):
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(
        n, spring_constant
    )
    eigenvalues, eigenvectors = np.linalg.eigh(-connectivity)
    inverse = np.zeros_like(eigenvalues)
    inverse[eigenvalues > 1e-12] = 1.0 / eigenvalues[eigenvalues > 1e-12]
    gram = 3.0 * (eigenvectors * inverse) @ eigenvectors.T
    diagonal = np.diag(gram)
    target = diagonal[:, None] + diagonal - 2.0 * gram
    target = 0.5 * (target + target.T)
    np.fill_diagonal(target, 0.0)
    return target


def test_distance_operator_and_adjoint_are_exact_duals():
    rng = np.random.default_rng(20260827)
    n = 7
    pair_i, pair_j = np.triu_indices(n, k=1)
    selector = np.arange(len(pair_i)) % 3 != 0
    pair_i = pair_i[selector]
    pair_j = pair_j[selector]

    matrix = rng.normal(size=(n, n))
    matrix = pdhg._center(0.5 * (matrix + matrix.T), np)
    dual = rng.normal(size=len(pair_i))

    left = np.dot(
        pdhg._distance_operator(matrix, pair_i, pair_j, np), dual
    )
    right = np.sum(
        matrix
        * pdhg._distance_adjoint(dual, n, pair_i, pair_j, np)
    )
    assert left == pytest.approx(right, rel=1e-13, abs=1e-13)


def test_logdet_proximal_map_is_centered_spd_and_stationary():
    rng = np.random.default_rng(20260828)
    n = 8
    matrix = rng.normal(size=(n, n))
    matrix = pdhg._center(0.5 * (matrix + matrix.T), np)
    householder = pdhg._householder_vector(n, np)
    step = 0.17

    fitted, inverse, _, eigenvalues = (
        pdhg._prox_centered_negative_logdet(
            matrix, step, householder, np
        )
    )

    assert np.max(np.abs(np.sum(fitted, axis=1))) <= 1e-12
    assert np.min(eigenvalues) > 0.0
    residual = (
        (fitted - matrix) / step
        - pdhg._ENTROPY_COEFFICIENT * inverse
    )
    assert np.linalg.norm(residual, ord="fro") <= 1e-11


def test_pdhg_matches_analytic_two_locus_solution():
    target_value = 2.7
    variance = 0.14
    target = np.array(
        [[0.0, target_value], [target_value, 0.0]], dtype=float
    )

    fitted, gram, connectivity, info = (
        HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            variance,
            max_iterations=500,
            relative_tolerance=1e-10,
            absolute_tolerance=1e-12,
        )
    )

    expected_distance = 0.5 * (
        target_value + np.sqrt(target_value**2 + 6.0 * variance)
    )
    assert info["converged"]
    assert fitted[0, 1] == pytest.approx(
        expected_distance, rel=1e-9, abs=1e-10
    )
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-12
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    assert fitted[0, 1] - target_value == pytest.approx(
        0.5 * variance * connectivity[0, 1],
        rel=1e-8,
        abs=1e-10,
    )
    assert info["relative_stationarity_residual"] <= 1e-8


def test_pdhg_reaches_same_small_system_optimum_as_newton():
    target = _rouse_squared_distance_target(6)
    relative_std = 0.05

    reference = HippsDimes.fit_gaussian_noise_covariance(
        target,
        relative_noise_std=relative_std,
        max_iterations=100,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
        initialization="nearest_edm",
    )
    fitted, gram, connectivity, info = (
        HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            relative_noise_std=relative_std,
            max_iterations=2000,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
        )
    )

    assert reference[3]["converged"]
    assert info["converged"]
    assert info["algorithm"] == "pdhg"
    assert np.allclose(fitted, reference[0], rtol=2e-7, atol=2e-8)
    assert np.allclose(gram, reference[1], rtol=2e-7, atol=2e-8)
    assert np.allclose(
        connectivity, reference[2], rtol=2e-6, atol=2e-7
    )
    assert info["relative_stationarity_residual"] <= 1e-7
    assert np.all(np.isfinite(info["history"]["objective"]))
    assert len(info["history"]["objective"]) == info["iterations"]
