"""Tests for variance-whitened, inverse-free PDHG."""

import numpy as np
import pytest

import HippsDimes
from hipps_dimes import covariance_pdhg as base
from hipps_dimes import covariance_pdhg_whitened as white
from hipps_dimes import numerics


def _rouse_target(n=6, spring_constant=1.0):
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(
        n, spring_constant
    )
    eigenvalues, eigenvectors = np.linalg.eigh(-connectivity)
    inverse = np.zeros_like(eigenvalues)
    inverse[eigenvalues > 1e-12] = 1.0 / eigenvalues[eigenvalues > 1e-12]
    gram = 3.0 * (eigenvectors * inverse) @ eigenvectors.T
    return numerics._squared_distances_from_gram(gram)


def test_weighted_distance_operator_and_adjoint_are_exact_duals():
    rng = np.random.default_rng(20260828)
    n = 7
    pair_i, pair_j = np.triu_indices(n, k=1)
    keep = np.arange(len(pair_i)) % 4 != 0
    pair_i = pair_i[keep]
    pair_j = pair_j[keep]
    inverse_std = np.exp(rng.normal(scale=0.7, size=len(pair_i)))
    matrix = rng.normal(size=(n, n))
    matrix = base._center(0.5 * (matrix + matrix.T), np)
    dual = rng.normal(size=len(pair_i))

    left = np.dot(
        white._weighted_distance_operator(
            matrix, inverse_std, pair_i, pair_j, np
        ),
        dual,
    )
    right = np.sum(
        matrix
        * white._weighted_distance_adjoint(
            dual, inverse_std, n, pair_i, pair_j, np
        )
    )
    assert left == pytest.approx(right, rel=1e-13, abs=1e-13)


def test_weighted_operator_power_iteration_matches_explicit_small_system():
    rng = np.random.default_rng(20260829)
    n = 5
    pair_i, pair_j = np.triu_indices(n, k=1)
    inverse_std = np.exp(rng.normal(scale=0.8, size=len(pair_i)))
    q = numerics._centered_orthonormal_basis(n)
    m = n - 1
    columns = []
    for i in range(m):
        for j in range(i, m):
            internal = np.zeros((m, m), dtype=float)
            if i == j:
                internal[i, j] = 1.0
            else:
                internal[i, j] = 1.0 / np.sqrt(2.0)
                internal[j, i] = 1.0 / np.sqrt(2.0)
            full = q @ internal @ q.T
            columns.append(
                white._weighted_distance_operator(
                    full, inverse_std, pair_i, pair_j, np
                )
            )
    design = np.column_stack(columns)
    exact = float(np.linalg.eigvalsh(design.T @ design)[-1])
    initial = q @ np.diag(np.linspace(1.0, 2.0, m)) @ q.T
    estimated, info = white._estimate_weighted_operator_norm_squared(
        initial,
        inverse_std,
        n,
        pair_i,
        pair_j,
        np,
        max_iterations=80,
        relative_tolerance=1e-10,
        safety_factor=1.05,
    )

    assert info["rayleigh_quotient"] == pytest.approx(
        exact, rel=2e-7
    )
    assert estimated > exact
    assert estimated < 1.08 * exact


def test_inverse_free_eliminated_kkt_matches_explicit_inverse():
    rng = np.random.default_rng(20260830)
    n = 6
    pair_i, pair_j = np.triu_indices(n, k=1)
    variance = 0.05 + rng.random(len(pair_i))
    inverse_std = 1.0 / np.sqrt(variance)
    target = 1.0 + rng.random(len(pair_i))
    whitened_target = target * inverse_std
    householder = base._householder_vector(n, np)
    q = numerics._centered_orthonormal_basis(n)
    internal = np.diag(np.linspace(0.8, 1.7, n - 1))
    previous = q @ internal @ q.T
    dual = rng.normal(scale=0.1, size=len(pair_i))
    tau = 0.03
    dual_adjoint = white._weighted_distance_adjoint(
        dual, inverse_std, n, pair_i, pair_j, np
    )
    gram, _, eigenvalues, eigenvectors = (
        white._prox_centered_negative_logdet_inverse_free(
            previous - tau * dual_adjoint,
            tau,
            householder,
            np,
        )
    )
    fitted = white._weighted_distance_operator(
        gram, inverse_std, pair_i, pair_j, np
    )
    residuals = white._inverse_free_residuals(
        previous,
        gram,
        dual_adjoint,
        dual,
        fitted,
        whitened_target,
        inverse_std,
        tau,
        n,
        pair_i,
        pair_j,
        np,
    )
    inverse = white._inverse_from_internal_spectrum(
        eigenvalues, eigenvectors, householder, np
    )
    explicit = white._weighted_distance_adjoint(
        fitted - whitened_target,
        inverse_std,
        n,
        pair_i,
        pair_j,
        np,
    ) - 1.5 * inverse

    assert np.allclose(
        residuals["eliminated_residual"],
        explicit,
        rtol=2e-11,
        atol=2e-11,
    )


def test_whitened_pdhg_matches_original_objective_on_heteroskedastic_case():
    target = _rouse_target(6, spring_constant=1.3)
    target[0, 5] *= 0.82
    target[5, 0] = target[0, 5]

    with pytest.warns(RuntimeWarning, match="without satisfying"):
        old = HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            relative_noise_std=0.12,
            max_iterations=120,
            relative_tolerance=1e-10,
            adaptive_steps=True,
        )
    new = HippsDimes.fit_gaussian_noise_covariance_pdhg_whitened(
        target,
        relative_noise_std=0.12,
        max_iterations=2000,
        relative_tolerance=2e-5,
        operator_norm_power_iterations=60,
        operator_norm_tolerance=1e-8,
    )

    assert new[3]["converged"]
    assert new[3]["relative_eliminated_kkt_residual"] <= 2e-5
    assert new[3]["objective"] <= old[3]["objective"] + 1e-5
    assert new[3]["pdhg"]["variance_whitened"]
    assert new[3]["pdhg"]["inverse_free_runtime_kkt"]
    assert new[3]["pdhg"]["weighted_residual_balancing"]


def test_inverse_is_not_reconstructed_at_every_iteration():
    target = _rouse_target(5, spring_constant=1.1)
    with pytest.warns(RuntimeWarning, match="without satisfying"):
        _, _, _, info = (
            HippsDimes.fit_gaussian_noise_covariance_pdhg_whitened(
                target,
                relative_noise_std=0.1,
                max_iterations=8,
                relative_tolerance=1e-14,
                save_steps=[4],
                adaptive_steps=False,
                operator_norm_power_iterations=12,
            )
        )

    assert info["inverse_reconstruction_count"] == 4
    assert not info["inverse_reconstructed_each_iteration"]
    assert len(info["history"]["iteration"]) == 8
    assert np.isfinite(
        info["history"]["relative_eliminated_kkt_residual"]
    ).all()
