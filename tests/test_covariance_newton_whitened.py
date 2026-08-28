"""Tests for entropy-whitened Newton-PCG and its PDHG handoff."""

import numpy as np
import pytest

import HippsDimes
from hipps_dimes import covariance_newton_whitened as newton
from hipps_dimes import numerics


def _rouse_target(n=6, spring_constant=1.0):
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(
        n, spring_constant
    )
    eigenvalues, eigenvectors = np.linalg.eigh(-connectivity)
    inverse = np.zeros_like(eigenvalues)
    inverse[eigenvalues > 1e-12] = 1.0 / eigenvalues[eigenvalues > 1e-12]
    gram = 3.0 * (eigenvectors * inverse) @ eigenvectors.T
    return numerics._squared_distances_from_gram(gram), gram


def test_svec_roundtrip_and_frobenius_isometry():
    rng = np.random.default_rng(20260829)
    size = 7
    matrix = rng.normal(size=(size, size))
    matrix = 0.5 * (matrix + matrix.T)
    diagonal, off_i, off_j = newton._svec_index_arrays(size, np)
    vector = newton._svec(matrix, diagonal, off_i, off_j, np)
    recovered = newton._smat(
        vector, size, diagonal, off_i, off_j, np
    )

    assert np.allclose(recovered, matrix, rtol=0.0, atol=0.0)
    assert np.dot(vector, vector) == pytest.approx(
        np.sum(matrix * matrix), rel=1e-14
    )


def test_inexact_newton_forcing_is_not_capped_by_final_tolerance():
    assert newton._forcing_tolerance(0.177, 1e-6, 0.1, 0.1) == pytest.approx(
        0.1 * np.sqrt(0.177)
    )
    assert newton._forcing_tolerance(1e-3, 1e-6, 0.1, 0.1) == pytest.approx(
        0.1 * np.sqrt(1e-3)
    )
    assert newton._forcing_tolerance(1e-12, 1e-6, 0.1, 0.1) == 1e-6


def test_entropy_whitened_hessian_matches_congruence_transform():
    rng = np.random.default_rng(20260830)
    n = 6
    target, _ = _rouse_target(n, 1.2)
    pair_i, pair_j = np.triu_indices(n, k=1)
    inverse_variance = 1.0 / (
        0.05 + 0.02 * (pair_j - pair_i)
    )
    _, inverse_variance_matrix = numerics._gaussian_covariance_pair_matrices(
        n,
        pair_i,
        pair_j,
        target[pair_i, pair_j],
        inverse_variance,
    )
    basis = numerics._centered_orthonormal_basis(n)
    raw = rng.normal(size=(n - 1, n - 1))
    reduced_gram = raw @ raw.T + np.eye(n - 1)
    inverse_gram = np.linalg.inv(reduced_gram)
    cholesky = np.linalg.cholesky(reduced_gram)
    whitened_basis = basis @ cholesky
    direction = rng.normal(size=(n - 1, n - 1))
    direction = 0.5 * (direction + direction.T)

    whitened_action = (
        numerics._gaussian_covariance_data_hessian_action(
            direction,
            whitened_basis,
            inverse_variance_matrix,
        )
        + 1.5 * direction
    )
    physical_direction = cholesky @ direction @ cholesky.T
    direct_action = numerics._gaussian_covariance_data_hessian_action(
        physical_direction,
        basis,
        inverse_variance_matrix,
    ) + 1.5 * inverse_gram @ physical_direction @ inverse_gram
    transformed_direct = cholesky.T @ direct_action @ cholesky

    assert np.allclose(
        whitened_action,
        transformed_direct,
        rtol=2e-12,
        atol=2e-12,
    )


def test_exact_whitened_svec_preconditioner_matches_operator_diagonal():
    rng = np.random.default_rng(20260831)
    n = 5
    modes = n - 1
    pair_i, pair_j = np.triu_indices(n, k=1)
    weights = 0.5 + rng.random(len(pair_i))
    inverse_variance_matrix = np.zeros((n, n))
    inverse_variance_matrix[pair_i, pair_j] = weights
    inverse_variance_matrix[pair_j, pair_i] = weights
    whitened_basis = rng.normal(size=(n, modes))
    data_diagonal, _ = numerics._gaussian_covariance_data_preconditioner_diagonal(
        whitened_basis,
        pair_i,
        pair_j,
        weights,
    )
    diagonal, off_i, off_j = newton._svec_index_arrays(modes, np)
    preconditioner = newton._whitened_svec_preconditioner(
        data_diagonal,
        diagonal,
        off_i,
        off_j,
        np,
    )

    explicit = []
    for index in range(len(preconditioner)):
        vector = np.zeros(len(preconditioner))
        vector[index] = 1.0
        matrix = newton._smat(
            vector, modes, diagonal, off_i, off_j, np
        )
        action = numerics._gaussian_covariance_data_hessian_action(
            matrix,
            whitened_basis,
            inverse_variance_matrix,
        ) + 1.5 * matrix
        explicit.append(
            np.dot(
                vector,
                newton._svec(action, diagonal, off_i, off_j, np),
            )
        )

    assert np.allclose(preconditioner, explicit, rtol=2e-13, atol=2e-13)


def test_pcg_reports_progress_and_obeys_hard_cap():
    matrix = np.diag([1.0, 100.0, 10000.0])
    events = []

    def operator(vector):
        return matrix @ vector

    _, info = newton._pcg_svec(
        operator,
        np.ones(3),
        np.ones(3),
        1e-14,
        1,
        xp=np,
        progress_callback=events.append,
        progress_interval=1,
        true_residual_interval=0,
    )

    assert info["iterations"] == 1
    assert not info["converged"]
    assert info["termination_reason"] == "maximum_iterations"
    assert events[0]["iteration"] == 0
    assert events[-1]["iteration"] == 1
    assert "cg_relative_residual" in events[-1]


def test_direct_gram_handoff_skips_scalar_calibration():
    target, gram = _rouse_target(5, 1.1)
    with pytest.warns(RuntimeWarning, match="without satisfying"):
        _, returned_gram, _, info = (
            HippsDimes.fit_gaussian_noise_covariance_newton_whitened(
                target,
                noise_variance=0.2,
                initial_gram=gram,
                max_iterations=1,
                readiness_probe=False,
                cg_max_iterations=20,
                relative_tolerance=1e-14,
            )
        )

    initialization = info["initialization"]
    assert initialization["kind"] == "provided_gram"
    assert not initialization["scalar_calibration_applied"]
    assert initialization["scalar_calibration"]["scale_factor"] == 1.0
    assert np.max(np.abs(np.sum(returned_gram, axis=1))) <= 1e-10


def test_failed_readiness_direction_returns_untouched_handoff(monkeypatch):
    target, gram = _rouse_target(5, 1.0)

    def zero_pcg(
        operator,
        right_hand_side,
        preconditioner,
        tolerance,
        maximum,
        **kwargs,
    ):
        del operator, preconditioner, tolerance, maximum, kwargs
        return np.zeros_like(right_hand_side), {
            "iterations": 1,
            "relative_residual": 1.0,
            "best_relative_residual": 1.0,
            "converged": False,
            "termination_reason": "maximum_iterations",
            "hessian_actions": 1,
            "wall_seconds": 0.0,
        }

    monkeypatch.setattr(newton, "_pcg_svec", zero_pcg)
    with pytest.warns(RuntimeWarning, match="newton_readiness_failed"):
        _, returned_gram, _, info = (
            HippsDimes.fit_gaussian_noise_covariance_newton_whitened(
                target,
                noise_variance=0.2,
                initial_gram=gram,
                max_iterations=3,
                readiness_probe=True,
                readiness_cg_max_iterations=1,
            )
        )

    assert info["status"] == "newton_readiness_failed"
    assert info["iterations"] == 0
    assert np.allclose(returned_gram, gram, rtol=2e-12, atol=2e-12)


def test_whitened_newton_never_accepts_a_kkt_increase():
    target, gram = _rouse_target(6, 1.2)
    target = target.copy()
    target[0, 5] *= 0.85
    target[5, 0] = target[0, 5]

    with pytest.warns(RuntimeWarning, match="without satisfying"):
        _, _, _, info = (
            HippsDimes.fit_gaussian_noise_covariance_newton_whitened(
                target,
                relative_noise_std=0.1,
                initial_gram=gram,
                max_iterations=4,
                relative_tolerance=1e-12,
                readiness_probe=True,
                readiness_cg_max_iterations=50,
                cg_max_iterations=100,
            )
        )

    history = info["history"]
    if len(history["iteration"]):
        assert np.all(history["kkt_after"] <= history["kkt_before"])
        assert info["relative_eliminated_kkt_residual"] <= (
            info["initialization"]["relative_kkt"]
        )
    else:
        assert info["status"] == "newton_readiness_failed"
