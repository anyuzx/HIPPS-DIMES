"""Tests for the primal-dual covariance-cone optimizer."""

import inspect

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


def _independent_relative_kkt(target, gram, variance):
    """Recompute the eliminated KKT residual without PDHG internals."""
    target = np.asarray(target, dtype=np.float64)
    gram = 0.5 * (
        np.asarray(gram, dtype=np.float64)
        + np.asarray(gram, dtype=np.float64).T
    )
    n = target.shape[0]
    observed = np.isfinite(target) & ~np.eye(n, dtype=bool)
    pair_i, pair_j = np.where(np.triu(observed, k=1))
    target_pairs = target[pair_i, pair_j]
    if np.isscalar(variance):
        pair_variance = np.full(len(pair_i), variance, dtype=np.float64)
    else:
        pair_variance = np.asarray(variance, dtype=np.float64)

    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    zero_index = int(np.argmin(np.abs(eigenvalues)))
    internal_mask = np.ones(n, dtype=bool)
    internal_mask[zero_index] = False
    internal_eigenvalues = eigenvalues[internal_mask]
    internal_vectors = eigenvectors[:, internal_mask]
    gram_pseudoinverse = (
        internal_vectors * (1.0 / internal_eigenvalues)
    ) @ internal_vectors.T

    diagonal = np.diag(gram)
    fitted_pairs = (
        diagonal[pair_i]
        + diagonal[pair_j]
        - 2.0 * gram[pair_i, pair_j]
    )
    eliminated_dual = (fitted_pairs - target_pairs) / pair_variance
    pair_matrix = np.zeros((n, n), dtype=np.float64)
    pair_matrix[pair_i, pair_j] = eliminated_dual
    pair_matrix[pair_j, pair_i] = eliminated_dual
    data_gradient = -pair_matrix
    np.fill_diagonal(data_gradient, np.sum(pair_matrix, axis=1))
    entropy_gradient = 1.5 * gram_pseudoinverse
    residual = data_gradient - entropy_gradient
    scale = max(
        np.linalg.norm(data_gradient),
        np.linalg.norm(entropy_gradient),
        np.finfo(np.float64).tiny,
    )
    return float(np.linalg.norm(residual) / scale)


def test_pdhg_default_relative_tolerance_is_production_value():
    parameter = inspect.signature(
        HippsDimes.fit_gaussian_noise_covariance_pdhg
    ).parameters["relative_tolerance"]

    assert parameter.default == 1e-5


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


def test_pdhg_reports_fresh_eliminated_kkt_certificate():
    target = _rouse_squared_distance_target(6)
    variance = 0.1

    _, gram, _, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target,
        variance,
        max_iterations=500,
    )
    independent = _independent_relative_kkt(target, gram, variance)

    assert info["converged"]
    assert info["independent_kkt_recomputed_from_returned_gram"]
    assert info["independent_kkt_converged"]
    assert independent == pytest.approx(
        info["relative_eliminated_kkt_residual"], rel=1e-12, abs=1e-14
    )
    assert info["relative_gradient_norm"] == pytest.approx(independent)
    assert info["relative_stationarity_residual"] == pytest.approx(independent)
    assert independent <= 1e-5
    assert np.allclose(
        info["history"]["relative_gradient_norm"],
        info["history"]["relative_eliminated_kkt_residual"],
    )


def test_pdhg_independent_certificate_controls_converged(monkeypatch):
    target = _rouse_squared_distance_target(4)
    original = pdhg._independent_eliminated_kkt_residuals

    def force_certificate_failure(*args, **kwargs):
        norm, scale, _, fitted, residual = original(*args, **kwargs)
        return max(norm, scale), scale, 1.0, fitted, residual

    monkeypatch.setattr(
        pdhg,
        "_independent_eliminated_kkt_residuals",
        force_certificate_failure,
    )

    with pytest.warns(RuntimeWarning, match="independent_kkt_failed"):
        _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            0.1,
            max_iterations=500,
        )

    assert info["termination_internal_kkt_converged"]
    assert not info["independent_kkt_converged"]
    assert not info["converged"]
    assert info["status"] == "independent_kkt_failed"


def test_pdhg_supports_heteroskedastic_noise_with_missing_pairs():
    target = _rouse_squared_distance_target(6)
    target[0, 4] = np.nan
    target[4, 0] = np.nan
    relative_std = 0.08

    _, gram, _, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target,
        relative_noise_std=relative_std,
        max_iterations=2000,
    )
    observed_pairs = target[np.triu(np.isfinite(target), k=1)]
    variance = np.square(relative_std * observed_pairs)
    independent = _independent_relative_kkt(target, gram, variance)

    assert info["converged"]
    assert info["noise_model"] == "heteroskedastic_relative_std"
    assert info["observed_pair_count"] == 14
    assert independent <= 1e-5


@pytest.mark.skipif(
    not HippsDimes.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
def test_pdhg_cpu_gpu_results_agree():
    target = _rouse_squared_distance_target(5)
    keyword_arguments = {
        "noise_variance": 0.1,
        "max_iterations": 1000,
        "relative_tolerance": 1e-7,
        "absolute_tolerance": 1e-10,
    }

    cpu = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target, use_gpu=False, **keyword_arguments
    )
    gpu = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target, use_gpu=True, **keyword_arguments
    )

    assert cpu[3]["converged"]
    assert gpu[3]["converged"]
    assert np.allclose(gpu[0], cpu[0], rtol=2e-6, atol=2e-8)
    assert np.allclose(gpu[1], cpu[1], rtol=2e-6, atol=2e-8)
    assert np.allclose(gpu[2], cpu[2], rtol=2e-5, atol=2e-7)
