"""Tests for the primal-dual covariance-cone optimizer."""

import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pytest

import HippsDimes
from hipps_dimes import covariance_pdhg as pdhg
from hipps_dimes import numerics

_REAL_DATA_DIRECTORY = (
    Path(__file__).parent / "data" / "gm12878_hic_chr1_31mb_41mb_25kb"
)


def _load_real_chromosome_case():
    metadata_path = _REAL_DATA_DIRECTORY / "metadata.json"
    contact_path = _REAL_DATA_DIRECTORY / "contact_hipps_ready.npy"
    reference_path = _REAL_DATA_DIRECTORY / "cov_pdhg_reference.npz"
    metadata = json.loads(metadata_path.read_text())
    contact = np.load(contact_path)
    with np.load(reference_path) as reference:
        gram = reference["gram_matrix"].copy()

    with np.errstate(divide="ignore"):
        distance = numerics.cmap2dmap_missing_data(
            contact,
            metadata["contact_to_squared_distance"]["alpha"],
            metadata["contact_to_squared_distance"]["not_normalize"],
        )
    target = (3.0 * np.pi / 8.0) * np.square(distance)
    target[np.isinf(target)] = np.nan
    return metadata, contact, target, gram


def _rouse_squared_distance_target(n=6, spring_constant=1.0):
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(n, spring_constant)
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
        np.asarray(gram, dtype=np.float64) + np.asarray(gram, dtype=np.float64).T
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
    fitted_pairs = diagonal[pair_i] + diagonal[pair_j] - 2.0 * gram[pair_i, pair_j]
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
    parameters = inspect.signature(
        HippsDimes.fit_gaussian_noise_covariance_pdhg
    ).parameters

    assert parameters["relative_tolerance"].default == 1e-5
    assert parameters["dual_initialization"].default == "residual"


@pytest.mark.parametrize(
    "name",
    [
        "fit_gaussian_noise_covariance_pdhg_whitened",
        "fit_gaussian_noise_covariance_preconditioned_pdhg",
        "fit_gaussian_noise_covariance_hybrid_whitened",
    ],
)
def test_pdhg_has_one_canonical_public_surface(name):
    assert not hasattr(HippsDimes, name)


def test_direct_b_fista_is_module_only():
    assert hasattr(pdhg, "fit_gaussian_noise_covariance_fista")
    assert not hasattr(HippsDimes, "fit_gaussian_noise_covariance_fista")
    assert not hasattr(HippsDimes, "fit_gaussian_noise_covariance")
    assert not hasattr(numerics, "fit_gaussian_noise_covariance")


def test_hybrid_defaults_to_fixed_production_handoff():
    signature = inspect.signature(HippsDimes.fit_gaussian_noise_covariance_hybrid)

    assert signature.parameters["relative_tolerance"].default == 1e-5
    assert signature.parameters["handoff_relative_tolerance"].default == 1e-2


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

    left = np.dot(pdhg._distance_operator(matrix, pair_i, pair_j, np), dual)
    right = np.sum(matrix * pdhg._distance_adjoint(dual, n, pair_i, pair_j, np))
    assert left == pytest.approx(right, rel=1e-13, abs=1e-13)


def test_weighted_distance_operator_and_adjoint_are_exact_duals():
    rng = np.random.default_rng(20260828)
    n = 7
    pair_i, pair_j = np.triu_indices(n, k=1)
    keep = np.arange(len(pair_i)) % 4 != 0
    pair_i = pair_i[keep]
    pair_j = pair_j[keep]
    inverse_std = np.exp(rng.normal(scale=0.7, size=len(pair_i)))
    matrix = rng.normal(size=(n, n))
    matrix = pdhg._center(0.5 * (matrix + matrix.T), np)
    dual = rng.normal(size=len(pair_i))

    left = np.dot(
        pdhg._weighted_distance_operator(matrix, inverse_std, pair_i, pair_j, np),
        dual,
    )
    right = np.sum(
        matrix
        * pdhg._weighted_distance_adjoint(dual, inverse_std, n, pair_i, pair_j, np)
    )
    assert left == pytest.approx(right, rel=1e-13, abs=1e-13)


def test_weighted_operator_norm_certificate_is_safe_on_small_system():
    rng = np.random.default_rng(20260829)
    n = 5
    pair_i, pair_j = np.triu_indices(n, k=1)
    inverse_std = np.exp(rng.normal(scale=0.8, size=len(pair_i)))
    basis = numerics._centered_orthonormal_basis(n)
    n_modes = n - 1
    columns = []
    for i in range(n_modes):
        for j in range(i, n_modes):
            internal = np.zeros((n_modes, n_modes), dtype=float)
            if i == j:
                internal[i, j] = 1.0
            else:
                internal[i, j] = 1.0 / np.sqrt(2.0)
                internal[j, i] = 1.0 / np.sqrt(2.0)
            columns.append(
                pdhg._weighted_distance_operator(
                    basis @ internal @ basis.T,
                    inverse_std,
                    pair_i,
                    pair_j,
                    np,
                )
            )
    design = np.column_stack(columns)
    exact = float(np.linalg.eigvalsh(design.T @ design)[-1])
    safe_bound, info = pdhg._certify_weighted_operator_norm_squared(
        inverse_std,
        n,
        pair_i,
        pair_j,
        np,
        max_iterations=80,
        relative_tolerance=1e-10,
        safety_factor=1.05,
        phase="pdhg",
    )

    assert info["method"] == "edge_space_collatz"
    assert info["bound_tolerance_reached"]
    assert info["lower_bound"] <= exact <= info["certified_upper_bound"]
    assert info["certified_upper_bound"] == pytest.approx(exact, rel=2e-7)
    assert safe_bound == pytest.approx(1.05 * info["certified_upper_bound"])


def test_weighted_operator_norm_certificate_handles_adversarial_weights():
    """The certified bound must protect PDHG on a connected weighted graph."""
    n = 6
    pair_i = np.array([0, 1, 2, 3, 4], dtype=np.int64)
    pair_j = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    inverse_std = np.sqrt(np.array([1.3, 1e-8, 1.0, 1e-8, 1e-8]))
    basis = numerics._centered_orthonormal_basis(n)
    columns = []
    for i in range(n - 1):
        for j in range(i, n - 1):
            internal = np.zeros((n - 1, n - 1), dtype=float)
            if i == j:
                internal[i, j] = 1.0
            else:
                internal[i, j] = 1.0 / np.sqrt(2.0)
                internal[j, i] = 1.0 / np.sqrt(2.0)
            columns.append(
                pdhg._weighted_distance_operator(
                    basis @ internal @ basis.T,
                    inverse_std,
                    pair_i,
                    pair_j,
                    np,
                )
            )
    design = np.column_stack(columns)
    exact = float(np.linalg.eigvalsh(design.T @ design)[-1])

    safe_bound, info = pdhg._certify_weighted_operator_norm_squared(
        inverse_std,
        n,
        pair_i,
        pair_j,
        np,
        max_iterations=25,
        relative_tolerance=1e-4,
        safety_factor=1.10,
        phase="pdhg",
    )

    assert exact <= info["certified_upper_bound"]
    assert 0.95**2 * exact / safe_bound < 1.0


def test_covariance_rejects_disconnected_observation_graph():
    target = np.full((4, 4), np.nan)
    np.fill_diagonal(target, 0.0)
    target[0, 1] = target[1, 0] = 1.0
    target[2, 3] = target[3, 2] = 1.0

    with pytest.raises(
        ValueError,
        match=r"must be connected; found 2 components with sizes \[2, 2\]",
    ):
        HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            noise_variance=0.1,
        )


def test_covariance_accepts_sparse_connected_observation_graph():
    target = np.full((5, 5), np.nan)
    np.fill_diagonal(target, 0.0)
    indices = np.arange(4)
    target[indices, indices + 1] = 1.0
    target[indices + 1, indices] = 1.0

    validated = numerics._validate_gaussian_covariance_inputs(
        target,
        noise_variance=0.1,
        relative_noise_std=None,
    )

    assert len(validated[1]) == 4


def test_inverse_free_eliminated_kkt_matches_explicit_inverse():
    rng = np.random.default_rng(20260830)
    n = 6
    pair_i, pair_j = np.triu_indices(n, k=1)
    variance = 0.05 + rng.random(len(pair_i))
    inverse_std = 1.0 / np.sqrt(variance)
    target = 1.0 + rng.random(len(pair_i))
    whitened_target = target * inverse_std
    householder = pdhg._householder_vector(n, np)
    basis = numerics._centered_orthonormal_basis(n)
    internal = np.diag(np.linspace(0.8, 1.7, n - 1))
    previous = basis @ internal @ basis.T
    dual = rng.normal(scale=0.1, size=len(pair_i))
    step = 0.03
    dual_adjoint = pdhg._weighted_distance_adjoint(
        dual, inverse_std, n, pair_i, pair_j, np
    )
    gram, _, eigenvalues, eigenvectors = (
        pdhg._prox_centered_negative_logdet_inverse_free(
            previous - step * dual_adjoint,
            step,
            householder,
            np,
        )
    )
    fitted = pdhg._weighted_distance_operator(gram, inverse_std, pair_i, pair_j, np)
    residuals = pdhg._inverse_free_residuals(
        previous,
        gram,
        dual_adjoint,
        dual,
        fitted,
        whitened_target,
        inverse_std,
        step,
        n,
        pair_i,
        pair_j,
        np,
    )
    inverse = pdhg._inverse_from_internal_spectrum(
        eigenvalues, eigenvectors, householder, np
    )
    explicit = (
        pdhg._weighted_distance_adjoint(
            fitted - whitened_target,
            inverse_std,
            n,
            pair_i,
            pair_j,
            np,
        )
        - 1.5 * inverse
    )

    assert np.allclose(
        residuals["eliminated_residual"],
        explicit,
        rtol=2e-11,
        atol=2e-11,
    )


def test_logdet_proximal_map_is_centered_spd_and_stationary():
    rng = np.random.default_rng(20260828)
    n = 8
    matrix = rng.normal(size=(n, n))
    matrix = pdhg._center(0.5 * (matrix + matrix.T), np)
    householder = pdhg._householder_vector(n, np)
    step = 0.17

    fitted, _, eigenvalues, eigenvectors = (
        pdhg._prox_centered_negative_logdet_inverse_free(matrix, step, householder, np)
    )
    inverse = pdhg._inverse_from_internal_spectrum(
        eigenvalues, eigenvectors, householder, np
    )

    assert np.max(np.abs(np.sum(fitted, axis=1))) <= 1e-12
    assert np.min(eigenvalues) > 0.0
    residual = (fitted - matrix) / step - pdhg._ENTROPY_COEFFICIENT * inverse
    assert np.linalg.norm(residual, ord="fro") <= 1e-11


def test_logdet_proximal_map_stably_handles_large_negative_eigenvalue():
    n = 3
    eigenvalues = np.array([-1e9, 2.0])
    internal = np.diag(eigenvalues)
    householder = pdhg._householder_vector(n, np)
    matrix = pdhg._full_centered_from_internal(internal, householder, np)
    step = 1.0

    fitted, logdet, updated, eigenvectors = (
        pdhg._prox_centered_negative_logdet_inverse_free(matrix, step, householder, np)
    )
    inverse = pdhg._inverse_from_internal_spectrum(
        updated, eigenvectors, householder, np
    )

    assert updated[0] == pytest.approx(1.5e-9, rel=1e-14)
    assert np.all(np.isfinite(updated))
    assert np.all(updated > 0.0)
    assert np.all(np.isfinite(fitted))
    assert np.all(np.isfinite(inverse))
    assert np.isfinite(logdet)
    assert np.allclose(
        updated * (updated - eigenvalues),
        pdhg._ENTROPY_COEFFICIENT * step,
        rtol=1e-14,
        atol=1e-14,
    )


@pytest.mark.parametrize(
    "use_gpu",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.skipif(
                not HippsDimes.is_gpu_available(),
                reason="requires CuPy and an accessible CUDA GPU",
            ),
        ),
    ],
)
def test_direct_b_fista_matches_analytic_two_locus_solution(use_gpu):
    target_value = 2.7
    variance = 0.14
    target = np.array([[0.0, target_value], [target_value, 0.0]], dtype=float)
    centering = np.eye(2) - np.ones((2, 2)) / 2.0
    initial_gram = -0.5 * centering @ target @ centering

    fitted, gram, connectivity, info = pdhg.fit_gaussian_noise_covariance_fista(
        target,
        variance,
        initial_gram=initial_gram,
        use_gpu=use_gpu,
        max_iterations=100,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )

    expected_distance = 0.5 * (target_value + np.sqrt(target_value**2 + 6.0 * variance))
    assert info["converged"]
    assert info["algorithm"] == "direct_b_monotone_fista"
    assert info["backend"] == ("gpu" if use_gpu else "cpu")
    assert fitted[0, 1] == pytest.approx(expected_distance, rel=1e-9, abs=1e-10)
    assert info["relative_eliminated_kkt_residual"] <= 1e-9
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-12
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12


def test_direct_b_fista_preserves_converged_input_gram_scale():
    target = _rouse_squared_distance_target(5)
    centering = np.eye(5) - np.ones((5, 5)) / 5.0
    initial_gram = -0.5 * centering @ target @ centering

    _, returned_gram, _, info = pdhg.fit_gaussian_noise_covariance_fista(
        target,
        noise_variance=0.1,
        initial_gram=initial_gram,
        max_iterations=10,
        relative_tolerance=2.0,
    )

    assert info["iterations"] == 0
    assert info["initialization"]["kind"] == "provided_gram"
    assert info["initialization"]["scalar_calibration"] is None
    assert info["initialization"]["physical_gram_used_directly"]
    assert np.allclose(returned_gram, initial_gram, rtol=1e-14, atol=1e-14)


def test_direct_b_fista_refines_missing_pair_handoff():
    target = _rouse_squared_distance_target(6)
    target[1, 4] = np.nan
    target[4, 1] = np.nan

    _, handoff_gram, _, handoff_info = pdhg.fit_gaussian_noise_covariance_pdhg(
        target,
        relative_noise_std=0.1,
        max_iterations=1000,
        relative_tolerance=1e-2,
        return_best=False,
    )
    fitted, gram, connectivity, info = pdhg.fit_gaussian_noise_covariance_fista(
        target,
        relative_noise_std=0.1,
        initial_gram=handoff_gram,
        max_iterations=1000,
        relative_tolerance=1e-7,
    )

    objective = info["history"]["objective"]
    objective_scale = max(1.0, float(np.max(np.abs(objective))))
    assert handoff_info["relative_eliminated_kkt_residual"] <= 1e-2
    assert info["converged"]
    assert info["relative_eliminated_kkt_residual"] <= 1e-7
    assert np.all(np.diff(objective) <= 1e-12 * objective_scale)
    assert np.all(info["history"]["minimum_internal_gram_eigenvalue"] > 0.0)
    assert np.all(np.isfinite(fitted))
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-12
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12


def test_direct_b_fista_rejects_uncentered_initial_gram():
    target = _rouse_squared_distance_target(4)
    initial_gram = np.eye(4)

    with pytest.raises(ValueError, match="must be centered"):
        pdhg.fit_gaussian_noise_covariance_fista(
            target,
            noise_variance=0.1,
            initial_gram=initial_gram,
        )


@pytest.mark.parametrize(
    "use_gpu",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.skipif(
                not HippsDimes.is_gpu_available(),
                reason="requires CuPy and an accessible CUDA GPU",
            ),
        ),
    ],
)
def test_pdhg_small_variance_inconsistent_target_stays_finite(use_gpu):
    target = np.array(
        [
            [0.0, 0.19, 54.47],
            [0.19, 0.0, 5.0],
            [54.47, 5.0, 0.0],
        ]
    )

    with pytest.warns(RuntimeWarning, match="without satisfying"):
        fitted, gram, connectivity, info = (
            HippsDimes.fit_gaussian_noise_covariance_pdhg(
                target,
                noise_variance=1e-8,
                max_iterations=20,
                use_gpu=use_gpu,
            )
        )

    assert not info["converged"]
    assert info["status"] == "independent_kkt_failed"
    assert info["termination_internal_kkt_converged"]
    assert not info["independent_kkt_converged"]
    assert np.all(np.isfinite(fitted))
    assert np.all(np.isfinite(gram))
    assert np.all(np.isfinite(connectivity))
    assert np.all(np.isfinite(info["history"]["objective"]))
    assert np.all(np.isfinite(info["history"]["gram_condition_number"]))
    assert np.all(info["history"]["minimum_internal_gram_eigenvalue"] > 0.0)


def test_pdhg_matches_analytic_two_locus_solution():
    target_value = 2.7
    variance = 0.14
    target = np.array([[0.0, target_value], [target_value, 0.0]], dtype=float)

    fitted, gram, connectivity, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target,
        variance,
        max_iterations=500,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )

    expected_distance = 0.5 * (target_value + np.sqrt(target_value**2 + 6.0 * variance))
    assert info["converged"]
    assert fitted[0, 1] == pytest.approx(expected_distance, rel=1e-9, abs=1e-10)
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-12
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    assert fitted[0, 1] - target_value == pytest.approx(
        0.5 * variance * connectivity[0, 1],
        rel=1e-8,
        abs=1e-10,
    )
    assert info["relative_stationarity_residual"] <= 1e-8


def test_pdhg_and_hybrid_reach_same_small_system_optimum():
    target = _rouse_squared_distance_target(6)
    relative_std = 0.05

    reference = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        relative_noise_std=relative_std,
        max_iterations=1000,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    fitted, gram, connectivity, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target,
        relative_noise_std=relative_std,
        max_iterations=2000,
        relative_tolerance=1e-8,
        absolute_tolerance=1e-10,
    )

    assert reference[3]["converged"]
    assert reference[3]["phase_iterations"]["fista"] > 0
    assert info["converged"]
    assert info["algorithm"] == "pdhg"
    assert info["pdhg"]["variance_whitened"]
    assert info["pdhg"]["inverse_free_runtime_kkt"]
    assert info["pdhg"]["weighted_residual_balancing"]
    assert np.allclose(fitted, reference[0], rtol=2e-7, atol=2e-8)
    assert np.allclose(gram, reference[1], rtol=3e-6, atol=2e-7)
    assert np.allclose(connectivity, reference[2], rtol=2e-6, atol=2e-7)
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


def test_pdhg_does_not_reconstruct_inverse_at_every_iteration():
    target = _rouse_squared_distance_target(5, spring_constant=1.1)
    with pytest.warns(RuntimeWarning, match="without satisfying"):
        _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            relative_noise_std=0.1,
            max_iterations=8,
            relative_tolerance=1e-14,
            save_steps=[4],
            adaptive_steps=False,
            operator_norm_power_iterations=12,
        )

    assert info["inverse_reconstruction_count"] == 4
    assert not info["inverse_reconstructed_each_iteration"]
    assert len(info["history"]["iteration"]) == 8
    assert np.isfinite(info["history"]["relative_eliminated_kkt_residual"]).all()


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


def test_hybrid_hands_off_to_fista_and_matches_cov_optimum():
    target = _rouse_squared_distance_target(6)
    variance = 0.1

    reference = HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target,
        variance,
        max_iterations=3000,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    fitted, gram, connectivity, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        variance,
        max_iterations=500,
    )

    assert reference[3]["converged"]
    assert info["converged"]
    assert info["algorithm"] == "hybrid"
    assert info["pdhg"]["variance_whitened"]
    assert info["pdhg"]["inverse_free_runtime_kkt"]
    assert info["handoff"]["reached"]
    assert info["handoff"]["relative_eliminated_kkt_residual"] <= 1e-2
    assert info["handoff"]["physical_gram_used_directly"]
    assert info["handoff"]["scalar_recalibration"] is None
    assert info["phase_iterations"]["pdhg"] > 0
    assert info["phase_iterations"]["fista"] > 0
    assert sum(info["phase_iterations"].values()) == info["iterations"]
    assert "fista_weighted_operator_norm" not in info
    assert info["relative_eliminated_kkt_residual"] <= 1e-5
    assert info["independent_kkt_recomputed_from_returned_gram"]
    assert np.allclose(fitted, reference[0], rtol=2e-7, atol=2e-8)
    assert np.allclose(gram, reference[1], rtol=3e-6, atol=2e-7)
    assert np.allclose(connectivity, reference[2], rtol=2e-6, atol=2e-7)

    history = info["history"]
    assert len(history["iteration"]) == info["iterations"]
    assert np.array_equal(history["iteration"], np.arange(1, info["iterations"] + 1))
    assert set(history["phase"]) == {"pdhg", "fista"}
    assert history["phase"][-1] == "fista"
    assert np.all(np.isfinite(history["objective"]))
    assert np.isfinite(history["relative_eliminated_kkt_residual"][-1])


def test_hybrid_handoff_does_not_wait_for_split_residuals(monkeypatch):
    target = _rouse_squared_distance_target(6)
    original = pdhg._inverse_free_residuals

    def force_split_residuals_to_lag(*args, **kwargs):
        residuals = dict(original(*args, **kwargs))
        residuals["primal_relative"] = 1.0
        residuals["dual_relative"] = 1.0
        return residuals

    monkeypatch.setattr(
        pdhg,
        "_inverse_free_residuals",
        force_split_residuals_to_lag,
    )

    _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        0.1,
        max_iterations=500,
    )

    assert info["converged"]
    assert info["handoff"]["reached"]
    assert info["handoff"]["relative_eliminated_kkt_residual"] <= 1e-2
    assert info["phase_iterations"]["pdhg"] < 500
    assert info["phase_iterations"]["fista"] > 0
    pdhg_end = info["phase_iterations"]["pdhg"] - 1
    assert info["history"]["relative_primal_kkt_residual"][pdhg_end] == 1.0
    assert info["history"]["relative_dual_kkt_residual"][pdhg_end] == 1.0


def test_hybrid_reuses_phase_certificates(monkeypatch):
    target = _rouse_squared_distance_target(6)
    original = pdhg._independent_eliminated_kkt_residuals
    certificate_calls = 0

    def count_certificates(*args, **kwargs):
        nonlocal certificate_calls
        certificate_calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        pdhg,
        "_independent_eliminated_kkt_residuals",
        count_certificates,
    )

    _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        0.1,
        max_iterations=500,
    )

    assert info["converged"]
    assert info["phase_iterations"]["fista"] > 0
    assert certificate_calls == 2


def test_hybrid_rejects_handoff_stricter_than_final_tolerance():
    target = _rouse_squared_distance_target(4)

    with pytest.raises(ValueError, match="must not be smaller"):
        HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            0.1,
            relative_tolerance=1e-5,
            handoff_relative_tolerance=1e-6,
        )


def test_hybrid_reports_when_total_budget_ends_before_handoff():
    target = _rouse_squared_distance_target(6)

    with pytest.warns(RuntimeWarning):
        _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            0.1,
            max_iterations=1,
        )

    assert not info["converged"]
    assert info["status"] == "pdhg_handoff_not_reached"
    assert info["phase_iterations"] == {"pdhg": 1, "fista": 0}
    assert info["iterations"] == 1


def test_hybrid_accepts_final_certificate_at_iteration_budget():
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        0.1,
        max_iterations=1,
        relative_tolerance=1e-12,
        handoff_relative_tolerance=1e-2,
    )

    assert info["converged"]
    assert info["independent_kkt_converged"]
    assert info["status"] == "optimality_tolerance"
    assert info["phase_iterations"] == {"pdhg": 1, "fista": 0}
    assert info["relative_eliminated_kkt_residual"] <= 1e-12


def test_hybrid_distinguishes_handoff_from_final_tolerance_at_budget():
    target = _rouse_squared_distance_target(6)

    with pytest.warns(RuntimeWarning, match="iteration_budget_exhausted_at_handoff"):
        _, _, _, info = HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            0.1,
            max_iterations=1,
            relative_tolerance=1e-12,
            absolute_tolerance=1e-15,
            handoff_relative_tolerance=10.0,
        )

    assert info["handoff"]["reached"]
    assert not info["converged"]
    assert not info["independent_kkt_converged"]
    assert info["status"] == "iteration_budget_exhausted_at_handoff"


def test_real_n400_chromosome_reference_satisfies_cov_objective_and_kkt():
    """Protect the COV objective with an experimental chromosome fixture."""
    metadata, contact, target, gram = _load_real_chromosome_case()
    contact_path = _REAL_DATA_DIRECTORY / "contact_hipps_ready.npy"
    reference_path = _REAL_DATA_DIRECTORY / "cov_pdhg_reference.npz"
    contact_sha256 = hashlib.sha256(contact_path.read_bytes()).hexdigest()
    reference_sha256 = hashlib.sha256(reference_path.read_bytes()).hexdigest()

    n = metadata["region"]["n_bins"]
    pair_mask = np.triu(np.isfinite(target), k=1)
    pair_i, pair_j = np.where(pair_mask)
    target_pairs = target[pair_i, pair_j]
    relative_std = metadata["covariance_model"]["relative_noise_std"]
    pair_variance = np.square(relative_std * target_pairs)

    eigenvalues = np.linalg.eigvalsh(gram)
    zero_index = int(np.argmin(np.abs(eigenvalues)))
    internal_eigenvalues = np.delete(eigenvalues, zero_index)
    diagonal = np.diag(gram)
    fitted_pairs = diagonal[pair_i] + diagonal[pair_j] - 2.0 * gram[pair_i, pair_j]
    objective = -1.5 * np.sum(np.log(internal_eigenvalues))
    objective += 0.5 * np.sum(np.square(fitted_pairs - target_pairs) / pair_variance)
    relative_kkt = _independent_relative_kkt(target, gram, pair_variance)
    expected = metadata["reference_solution"]

    assert contact_sha256 == metadata["files"]["contact_map_sha256"]
    assert reference_sha256 == metadata["files"]["reference_sha256"]
    assert contact.shape == (n, n) == gram.shape
    assert np.array_equal(contact, contact.T)
    assert len(target_pairs) == expected["observed_pair_count"]
    assert n * (n - 1) // 2 - len(target_pairs) == 1660
    assert np.max(np.abs(np.sum(gram, axis=1))) <= 1e-10
    assert np.min(internal_eigenvalues) > 0.0
    assert objective == pytest.approx(expected["objective"], rel=1e-11, abs=1e-6)
    assert objective == pytest.approx(
        expected["trusted_reference_objective"], rel=1e-11, abs=1e-6
    )
    assert relative_kkt == pytest.approx(
        expected["independent_relative_eliminated_kkt_residual"],
        rel=1e-7,
        abs=1e-10,
    )
    assert relative_kkt <= metadata["solver"]["relative_tolerance"]


@pytest.mark.real_data
@pytest.mark.filterwarnings("ignore:divide by zero encountered:RuntimeWarning")
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
@pytest.mark.skipif(
    not HippsDimes.is_gpu_available(),
    reason="full N=400 convergence regression requires a CUDA GPU",
)
def test_default_hybrid_converges_from_rouse_on_real_n400_contact_map():
    """Run and certify the default public COV workflow on the real fixture."""
    metadata, contact, _, reference_gram = _load_real_chromosome_case()
    solver = metadata["solver"]
    expected = metadata["reference_solution"]

    results = HippsDimes.run_optimization(
        input_matrix=contact,
        input_type="cmap",
        input_format="npy",
        method="COV",
        gaussian_noise_relative_std=(
            metadata["covariance_model"]["relative_noise_std"]
        ),
        covariance_relative_tolerance=solver["relative_tolerance"],
        covariance_absolute_tolerance=solver["absolute_tolerance"],
        iteration=solver["maximum_iterations"],
        use_gpu=True,
        ignore_missing_data=True,
        repair_fully_missing_loci=True,
        not_normalize=True,
        no_log=True,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )
    info = results["covariance_optimization"]

    assert info["converged"]
    assert info["status"] == "optimality_tolerance"
    assert info["algorithm"] == "hybrid"
    assert info["pdhg"]["variance_whitened"]
    assert info["pdhg"]["inverse_free_runtime_kkt"]
    assert info["handoff_relative_tolerance"] == pytest.approx(
        solver["handoff_relative_tolerance"]
    )
    assert info["handoff"]["reached"]
    assert info["handoff"]["relative_eliminated_kkt_residual"] <= (
        solver["handoff_relative_tolerance"]
    )
    assert info["phase_iterations"]["pdhg"] > 0
    assert info["phase_iterations"]["pdhg"] <= 1500
    assert info["phase_iterations"]["fista"] > 0
    assert info["independent_kkt_converged"]
    assert info["iterations"] <= solver["maximum_iterations"]
    assert info["relative_eliminated_kkt_residual"] <= (solver["relative_tolerance"])
    assert info["objective"] == pytest.approx(
        expected["trusted_reference_objective"], rel=1e-8, abs=1e-6
    )

    fitted_gram = results["gram_matrix"]
    fitted_connectivity = results["connectivity_matrix"]
    fitted_distance = results["dmap_final"]
    fitted_squared_distance = (3.0 * np.pi / 8.0) * np.square(fitted_distance)
    reference_diagonal = np.diag(reference_gram)
    reference_squared_distance = (
        reference_diagonal[:, None] + reference_diagonal - 2.0 * reference_gram
    )
    reference_distance = np.sqrt(
        np.maximum(0.0, 8.0 * reference_squared_distance / (3.0 * np.pi))
    )
    eigenvalues, eigenvectors = np.linalg.eigh(reference_gram)
    internal = eigenvalues > 1e-12
    reference_inverse = (
        eigenvectors[:, internal] * (1.0 / eigenvalues[internal])
    ) @ eigenvectors[:, internal].T
    reference_connectivity = HippsDimes.a2a(-3.0 * reference_inverse)

    def relative_l2(observed, reference):
        return np.linalg.norm(observed - reference) / np.linalg.norm(reference)

    assert relative_l2(fitted_squared_distance, reference_squared_distance) <= 2e-5
    assert relative_l2(fitted_distance, reference_distance) <= 1e-5
    assert relative_l2(fitted_gram, reference_gram) <= 5e-5
    assert relative_l2(fitted_connectivity, reference_connectivity) <= 3e-5

    fitted_eigenvalues = np.linalg.eigvalsh(fitted_gram)
    zero_index = int(np.argmin(np.abs(fitted_eigenvalues)))
    assert np.min(np.delete(fitted_eigenvalues, zero_index)) > 0.0
    assert np.max(np.abs(fitted_gram - fitted_gram.T)) <= 1e-10
    assert np.max(np.abs(np.sum(fitted_gram, axis=1))) <= 1e-10
    assert np.max(np.abs(fitted_connectivity - fitted_connectivity.T)) <= 1e-10
    assert np.max(np.abs(np.sum(fitted_connectivity, axis=1))) <= 1e-10
    assert info["wall_seconds"] > 0.0
    assert info["phase_wall_seconds"]["pdhg"] > 0.0
    assert info["phase_wall_seconds"]["fista"] > 0.0


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
