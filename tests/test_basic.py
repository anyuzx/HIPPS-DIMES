"""Basic tests for HippsDimes package."""

import numpy as np
import pytest

# Import the main module
import HippsDimes
import hipps_dimes
import hipps_dimes.numerics as numerics


@pytest.fixture
def working_cupy():
    """Require a CuPy runtime that can execute kernels and dense eigensolvers."""
    if not numerics.is_gpu_available():
        pytest.skip("requires CuPy and an accessible CUDA GPU")
    cp = numerics.cp
    try:
        values = cp.arange(4, dtype=cp.float64)
        assert float(cp.sum(values * values).item()) == pytest.approx(14.0)
        with numerics.cupyx.errstate(linalg="raise"):
            cp.linalg.eigh(cp.eye(2, dtype=cp.float64))
        cp.cuda.get_current_stream().synchronize()
    except Exception as error:
        pytest.skip(f"CuPy kernel or eigensolver execution is unavailable: {error}")
    return cp


def test_import():
    """Test that the package can be imported."""
    assert HippsDimes is not None


def test_construct_connectivity_matrix_rouse():
    """Test the Rouse chain connectivity matrix construction."""
    n = 5
    k = 1.0
    A = HippsDimes.construct_connectivity_matrix_rouse(n, k)
    
    # Check shape
    assert A.shape == (n, n)
    
    # Check symmetry
    assert np.allclose(A, A.T)
    
    # Check row sums (should be zero for Laplacian)
    assert np.allclose(np.sum(A, axis=1), 0.0)


def test_a2a_enforces_nonnegative_spring_constants():
    """Negative off-diagonals are removed while positive springs are retained."""
    source = np.array(
        [
            [-9.0, 2.0, -3.0],
            [2.0, -9.0, 4.0],
            [-3.0, 4.0, -9.0],
        ]
    )
    source_before = source.copy()

    result = HippsDimes.a2a(source, fill_negative=True)

    expected = np.array(
        [
            [-2.0, 2.0, 0.0],
            [2.0, -6.0, 4.0],
            [0.0, 4.0, -4.0],
        ]
    )
    assert np.array_equal(source, source_before)
    assert np.array_equal(result, expected)
    assert np.allclose(np.sum(result, axis=1), 0.0)
    assert np.all(np.linalg.eigvalsh(result) <= 1e-12)


@pytest.mark.skipif(
    not HippsDimes.is_gpu_available(),
    reason="CuPy GPU is not available",
)
def test_a2a_enforces_nonnegative_spring_constants_on_gpu(working_cupy):
    """GPU normalization should preserve the backend and match the CPU result."""
    source_cpu = np.array(
        [
            [-9.0, 2.0, -3.0],
            [2.0, -9.0, 4.0],
            [-3.0, 4.0, -9.0],
        ]
    )
    source_gpu = working_cupy.asarray(source_cpu)

    result_gpu = HippsDimes.a2a(source_gpu, fill_negative=True)

    expected = HippsDimes.a2a(source_cpu, fill_negative=True)
    assert isinstance(result_gpu, working_cupy.ndarray)
    assert np.array_equal(working_cupy.asnumpy(source_gpu), source_cpu)
    assert np.array_equal(working_cupy.asnumpy(result_gpu), expected)


def test_a2dmap_theory():
    """Test the theoretical distance map calculation."""
    n = 4
    k = 1.0
    A = HippsDimes.construct_connectivity_matrix_rouse(n, k)
    
    dmap = HippsDimes.a2dmap_theory(A)
    
    # Check shape
    assert dmap.shape == (n, n)
    
    # Check symmetry
    assert np.allclose(dmap, dmap.T)
    
    # Check diagonal (should be zero)
    assert np.allclose(np.diag(dmap), 0.0)


def test_centered_psd_projection_constraints_and_idempotence():
    """The EDM Gram projection should be centered, PSD, and idempotent."""
    source = np.array(
        [
            [2.0, -3.0, 1.0],
            [-3.0, -1.0, 4.0],
            [1.0, 4.0, 0.5],
        ]
    )

    projected = numerics._project_centered_psd(source)
    projected_twice = numerics._project_centered_psd(projected)

    assert np.allclose(projected, projected.T, atol=1e-12)
    assert np.allclose(np.sum(projected, axis=1), 0.0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(projected)) >= -1e-12
    assert np.allclose(projected_twice, projected, rtol=1e-11, atol=1e-12)


def test_centered_psd_projection_enforces_internal_eigenvalue_floor():
    """A positive floor should preserve COM while bounding all internal modes."""
    source = np.array(
        [
            [2.0, -3.0, 1.0, 0.5],
            [-3.0, -1.0, 4.0, -2.0],
            [1.0, 4.0, 0.5, 1.5],
            [0.5, -2.0, 1.5, -0.5],
        ]
    )
    floor = 0.25
    n = len(source)
    centering = np.eye(n) - np.ones((n, n)) / n

    projected = numerics._project_centered_psd(source, floor)
    projected_twice = numerics._project_centered_psd(projected, floor)
    eigenvalues = np.linalg.eigvalsh(projected)

    assert np.allclose(projected, projected.T, atol=1e-12)
    assert np.allclose(np.sum(projected, axis=1), 0.0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(projected - floor * centering)) >= -1e-12
    assert np.count_nonzero(eigenvalues > 0.5 * floor) == n - 1
    assert np.allclose(projected_twice, projected, rtol=1e-11, atol=1e-12)


def test_nearest_edm_gradient_matches_finite_difference():
    """The weighted unique-pair objective should have the implemented gradient."""
    target = np.array(
        [
            [0.0, 1.0, 4.0],
            [1.0, 0.0, 2.0],
            [4.0, 2.0, 0.0],
        ]
    )
    weights = np.array(
        [
            [0.0, 1.0, 0.5],
            [1.0, 0.0, 2.0],
            [0.5, 2.0, 0.0],
        ]
    )
    source = np.array(
        [
            [1.0, -0.2, 0.3],
            [-0.2, 0.5, -0.1],
            [0.3, -0.1, 0.8],
        ]
    )
    gram = numerics._project_centered_psd(source)
    direction = np.array(
        [
            [0.2, -0.5, 0.7],
            [-0.5, 0.4, 0.1],
            [0.7, 0.1, -0.3],
        ]
    )

    _, gradient, _ = numerics._nearest_edm_objective_gradient(
        gram, target, weights
    )
    epsilon = 1e-6
    forward, _, _ = numerics._nearest_edm_objective_gradient(
        gram + epsilon * direction, target, weights
    )
    backward, _, _ = numerics._nearest_edm_objective_gradient(
        gram - epsilon * direction, target, weights
    )

    finite_difference = (forward - backward) / (2.0 * epsilon)
    assert finite_difference == pytest.approx(
        np.sum(gradient * direction), rel=1e-8, abs=1e-9
    )


def test_nearest_edm_matches_analytic_three_point_solution_and_exports():
    """The closest invalid three-point EDM has a closed-form boundary solution."""
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )
    expected = np.array(
        [
            [0.0, 19.0 / 9.0, 19.0 / 9.0],
            [19.0 / 9.0, 0.0, 76.0 / 9.0],
            [19.0 / 9.0, 76.0 / 9.0, 0.0],
        ]
    )

    progress_events = []
    fitted, gram, info = hipps_dimes.nearest_edm(
        target, progress_callback=progress_events.append
    )

    assert info["converged"]
    assert info["status"] == "optimality_tolerance"
    assert np.allclose(fitted, expected, rtol=1e-8, atol=1e-9)
    assert np.allclose(np.sum(gram, axis=1), 0.0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(gram)) >= -1e-12
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-12)
    assert info["backend"] == "cpu"
    assert info["dtype"] == "float64"
    assert info["gpu_device"] is None
    assert info["cupy_version"] is None
    assert info["wall_seconds"] >= 0.0
    assert info["projection_count"] == (
        1
        + info["line_search_projection_count"]
        + info["certificate_projection_count"]
    )
    assert info["certificate_projection_count"] == info["iterations"]
    assert info["gpu_memory_pool_baseline_used_bytes"] is None
    assert info["gpu_memory_pool_maximum_used_bytes"] is None
    assert info["gpu_memory_pool_baseline_total_bytes"] is None
    assert info["gpu_memory_pool_maximum_total_bytes"] is None
    assert len(progress_events) == info["iterations"]
    assert all(
        event["stage"] == "nearest_edm_initialization"
        and event["use_gpu"] is False
        for event in progress_events
    )
    assert HippsDimes.nearest_edm is hipps_dimes.nearest_edm


def test_nearest_edm_leaves_valid_squared_distances_unchanged():
    coordinates = np.array(
        [
            [-1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, -1.0],
        ]
    )
    differences = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    target = np.sum(differences * differences, axis=-1)

    fitted, gram, info = hipps_dimes.nearest_edm(target)

    assert info["converged"]
    assert info["objective"] == pytest.approx(0.0, abs=1e-20)
    assert np.allclose(fitted, target, rtol=1e-10, atol=1e-11)
    assert np.allclose(np.sum(gram, axis=1), 0.0, atol=1e-12)


def test_nearest_edm_floor_has_analytic_zero_target_solution():
    """The smallest feasible Gram matrix should fit an all-zero target best."""
    n = 4
    floor = 0.2
    target = np.zeros((n, n))
    centering = np.eye(n) - np.ones((n, n)) / n
    expected_gram = floor * centering
    expected_fitted = np.full((n, n), 2.0 * floor)
    np.fill_diagonal(expected_fitted, 0.0)

    fitted, gram, info = hipps_dimes.nearest_edm(
        target, gram_eigenvalue_floor=floor
    )

    assert info["converged"]
    assert info["gram_eigenvalue_floor"] == pytest.approx(floor)
    assert np.allclose(gram, expected_gram, rtol=1e-11, atol=1e-12)
    assert np.allclose(fitted, expected_fitted, rtol=1e-11, atol=1e-12)
    assert np.count_nonzero(np.linalg.eigvalsh(gram) > 0.5 * floor) == n - 1


def test_nearest_edm_honors_weights_and_is_invariant_to_weight_scale():
    """Pair weights should move the optimum without depending on global scale."""
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )
    leg_weight = 2.0
    base_weight = 0.25
    weights = np.array(
        [
            [0.0, leg_weight, leg_weight],
            [leg_weight, 0.0, base_weight],
            [leg_weight, base_weight, 0.0],
        ]
    )
    expected_leg = (leg_weight + 18.0 * base_weight) / (
        leg_weight + 8.0 * base_weight
    )
    expected = np.array(
        [
            [0.0, expected_leg, expected_leg],
            [expected_leg, 0.0, 4.0 * expected_leg],
            [expected_leg, 4.0 * expected_leg, 0.0],
        ]
    )

    fitted, _, info = hipps_dimes.nearest_edm(target, weights)
    scaled_fitted, _, scaled_info = hipps_dimes.nearest_edm(
        target, 7.0 * weights
    )

    assert info["converged"]
    assert scaled_info["converged"]
    assert np.allclose(fitted, expected, rtol=1e-8, atol=1e-9)
    assert np.allclose(scaled_fitted, fitted, rtol=1e-10, atol=1e-11)


def test_nearest_edm_completes_missing_pairs_without_imputation():
    """NaNs should remove pairs from the objective while retaining a valid EDM."""
    coordinates = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ]
    )
    differences = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    target = np.sum(differences * differences, axis=-1)
    target[0, 3] = np.nan
    target[3, 0] = np.nan

    fitted, gram, info = hipps_dimes.nearest_edm(target)
    observed = np.isfinite(target) & ~np.eye(len(target), dtype=bool)

    assert info["converged"]
    assert np.allclose(fitted[observed], target[observed], atol=1e-7)
    assert np.allclose(np.diag(fitted), 0.0)
    assert np.min(fitted) >= -1e-12
    assert np.min(np.linalg.eigvalsh(gram)) >= -1e-12


@pytest.mark.parametrize(
    ("target", "weights", "message"),
    [
        (np.ones((2, 3)), None, "square"),
        (np.array([[0.0, np.inf], [np.inf, 0.0]]), None, "infinite"),
        (np.array([[0.0, 1.0], [2.0, 0.0]]), None, "symmetric"),
        (
            np.array([[0.0, 1.0], [1.0, 0.0]]),
            np.array([[0.0, -1.0], [-1.0, 0.0]]),
            "nonnegative",
        ),
        (
            np.array([[0.0, np.nan], [np.nan, 0.0]]),
            np.array([[0.0, 1.0], [1.0, 0.0]]),
            "finite distance",
        ),
    ],
)
def test_nearest_edm_rejects_invalid_inputs(target, weights, message):
    with pytest.raises(ValueError, match=message):
        hipps_dimes.nearest_edm(target, weights)


def test_nearest_edm_reports_iteration_limit():
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )

    with pytest.warns(RuntimeWarning, match="max_iterations"):
        _, _, info = hipps_dimes.nearest_edm(
            target,
            max_iterations=1,
            relative_tolerance=1e-15,
            absolute_tolerance=0.0,
        )

    assert not info["converged"]
    assert info["status"] == "max_iterations"
    assert info["iterations"] == 1


@pytest.mark.parametrize("floor", [-0.1, np.nan, np.inf, True, [0.1]])
def test_nearest_edm_rejects_invalid_gram_eigenvalue_floor(floor):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="finite nonnegative scalar"):
        hipps_dimes.nearest_edm(target, gram_eigenvalue_floor=floor)


def test_nearest_edm_gpu_request_fails_without_available_gpu(monkeypatch):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])
    monkeypatch.setattr(numerics, "_CUPY_AVAILABLE", False)

    with pytest.raises(RuntimeError, match="CuPy.*accessible CUDA GPU"):
        hipps_dimes.nearest_edm(target, use_gpu=True)


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
def test_nearest_edm_gpu_helpers_match_cpu(working_cupy):
    """The shared float64 projection and objective helpers should match."""
    source = np.array(
        [
            [1.0, -0.4, 0.3, 0.2],
            [-0.4, 0.8, -0.1, 0.4],
            [0.3, -0.1, 0.5, -0.2],
            [0.2, 0.4, -0.2, 1.2],
        ]
    )
    target = np.array(
        [
            [0.0, 1.0, 4.0, 2.0],
            [1.0, 0.0, 2.5, 3.0],
            [4.0, 2.5, 0.0, 1.5],
            [2.0, 3.0, 1.5, 0.0],
        ]
    )
    weights = np.array(
        [
            [0.0, 1.0, 0.5, 0.0],
            [1.0, 0.0, 2.0, 0.75],
            [0.5, 2.0, 0.0, 1.25],
            [0.0, 0.75, 1.25, 0.0],
        ]
    )
    floor = 0.05

    cpu_projection = numerics._project_centered_psd(source, floor)
    gpu_projection = numerics._project_centered_psd(
        working_cupy.asarray(source),
        floor,
        array_module=working_cupy,
    )
    cpu_objective, cpu_gradient, cpu_fitted = (
        numerics._nearest_edm_objective_gradient(
            cpu_projection, target, weights
        )
    )
    gpu_objective, gpu_gradient, gpu_fitted = (
        numerics._nearest_edm_objective_gradient(
            gpu_projection,
            working_cupy.asarray(target),
            working_cupy.asarray(weights),
            array_module=working_cupy,
        )
    )

    assert np.allclose(
        working_cupy.asnumpy(gpu_projection),
        cpu_projection,
        rtol=1e-12,
        atol=1e-12,
    )
    assert gpu_objective == pytest.approx(cpu_objective, rel=1e-12, abs=1e-12)
    assert np.allclose(
        working_cupy.asnumpy(gpu_gradient), cpu_gradient, rtol=1e-12, atol=1e-12
    )
    assert np.allclose(
        working_cupy.asnumpy(gpu_fitted), cpu_fitted, rtol=1e-12, atol=1e-12
    )


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
@pytest.mark.parametrize("case", ["weighted_floor", "missing_pair"])
def test_nearest_edm_gpu_matches_cpu_end_to_end(case, working_cupy):
    """GPU nearest EDM should preserve weights, floors, and missing-pair logic."""
    if case == "weighted_floor":
        target = np.array(
            [
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 9.0],
                [1.0, 9.0, 0.0],
            ]
        )
        weights = np.array(
            [
                [0.0, 2.0, 2.0],
                [2.0, 0.0, 0.25],
                [2.0, 0.25, 0.0],
            ]
        )
        floor = 1e-3
    else:
        coordinates = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
            ]
        )
        differences = (
            coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        )
        target = np.sum(differences * differences, axis=-1)
        target[0, 3] = np.nan
        target[3, 0] = np.nan
        weights = None
        floor = 0.0

    options = {
        "gram_eigenvalue_floor": floor,
        "max_iterations": 1000,
        "relative_tolerance": 1e-8,
        "absolute_tolerance": 1e-10,
    }
    cpu = hipps_dimes.nearest_edm(target, weights, **options)
    progress_events = []
    gpu = hipps_dimes.nearest_edm(
        target,
        weights,
        use_gpu=True,
        progress_callback=progress_events.append,
        **options,
    )

    assert cpu[2]["converged"]
    assert gpu[2]["converged"]
    assert isinstance(gpu[0], np.ndarray)
    assert isinstance(gpu[1], np.ndarray)
    assert gpu[2]["backend"] == "gpu"
    assert gpu[2]["dtype"] == "float64"
    assert gpu[2]["gpu_device"] == numerics.get_gpu_name()
    assert gpu[2]["cupy_version"] == working_cupy.__version__
    assert gpu[2]["wall_seconds"] >= 0.0
    assert gpu[2]["gpu_memory_pool_maximum_used_bytes"] >= gpu[2][
        "gpu_memory_pool_baseline_used_bytes"
    ]
    assert gpu[2]["gpu_memory_pool_maximum_total_bytes"] >= gpu[2][
        "gpu_memory_pool_baseline_total_bytes"
    ]
    assert gpu[2]["certificate_projection_count"] == gpu[2]["iterations"]
    assert len(progress_events) == gpu[2]["iterations"]
    assert all(event["use_gpu"] for event in progress_events)
    assert np.allclose(gpu[0], cpu[0], rtol=1e-7, atol=1e-8)
    assert np.allclose(gpu[1], cpu[1], rtol=1e-7, atol=1e-8)


def test_gaussian_covariance_objective_gradient_matches_finite_difference():
    """The calibrated COV objective should expose its true matrix gradient."""
    n = 5
    basis = numerics._centered_orthonormal_basis(n)
    pair_i, pair_j = np.triu_indices(n, k=1)
    pair_vectors = basis[pair_i] - basis[pair_j]
    reduced_gram = np.diag(np.linspace(0.8, 2.0, n - 1))
    target = np.einsum(
        "ij,ij->i", pair_vectors @ reduced_gram, pair_vectors, optimize=True
    )
    target *= np.linspace(0.95, 1.05, len(target))
    inverse_variance = np.linspace(2.0, 5.0, len(target))
    target_matrix, inverse_variance_matrix = (
        numerics._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target, inverse_variance
        )
    )
    direction = np.arange((n - 1) ** 2, dtype=float).reshape(n - 1, n - 1)
    direction = 0.5 * (direction + direction.T)
    direction /= np.linalg.norm(direction)

    objective, gradient, *_ = numerics._gaussian_covariance_objective_gradient(
        reduced_gram,
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
    )
    step = 1e-6
    forward = numerics._gaussian_covariance_objective_gradient(
        reduced_gram + step * direction,
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
    )[0]
    backward = numerics._gaussian_covariance_objective_gradient(
        reduced_gram - step * direction,
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
    )[0]

    assert np.isfinite(objective)
    assert (forward - backward) / (2.0 * step) == pytest.approx(
        np.sum(gradient * direction), rel=1e-7, abs=1e-8
    )


def test_covariance_structured_hessian_and_block_preconditioner_match_pairs():
    """Structured and blockwise COV kernels should equal explicit pair formulas."""
    rng = np.random.default_rng(20260826)
    n = 7
    basis = numerics._centered_orthonormal_basis(n)
    all_i, all_j = np.triu_indices(n, k=1)
    selector = np.arange(len(all_i)) % 4 != 0
    pair_i = all_i[selector]
    pair_j = all_j[selector]
    pair_vectors = basis[pair_i] - basis[pair_j]
    inverse_variance = np.linspace(1.5, 4.5, len(pair_i))
    target_pairs = np.linspace(0.9, 1.3, len(pair_i))
    _, inverse_variance_matrix = numerics._gaussian_covariance_pair_matrices(
        n, pair_i, pair_j, target_pairs, inverse_variance
    )
    direction = rng.normal(size=(n - 1, n - 1))
    direction = 0.5 * (direction + direction.T)

    pair_quadratic = np.einsum(
        "ij,ij->i", pair_vectors @ direction, pair_vectors, optimize=True
    )
    expected_hessian = pair_vectors.T @ (
        (inverse_variance * pair_quadratic)[:, np.newaxis] * pair_vectors
    )
    expected_hessian = 0.5 * (expected_hessian + expected_hessian.T)
    observed_hessian = numerics._gaussian_covariance_data_hessian_action(
        direction, basis, inverse_variance_matrix
    )

    pair_vectors_squared = np.square(pair_vectors)
    expected_diagonal = pair_vectors_squared.T @ (
        inverse_variance[:, np.newaxis] * pair_vectors_squared
    )
    observed_diagonal, _ = (
        numerics._gaussian_covariance_data_preconditioner_diagonal(
            basis,
            pair_i,
            pair_j,
            inverse_variance,
            block_size=3,
        )
    )

    assert np.allclose(
        observed_hessian, expected_hessian, rtol=1e-12, atol=1e-12
    )
    assert np.allclose(
        observed_diagonal, expected_diagonal, rtol=1e-12, atol=1e-12
    )


def _rouse_squared_distance_target(n=6, spring_constant=1.0):
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, spring_constant)
    mean_distance = HippsDimes.a2dmap_theory(
        truth, force_positive_definite=True
    )
    return (3.0 * np.pi / 8.0) * np.square(mean_distance)


@pytest.mark.parametrize(
    "missing_pairs",
    [(), ((0, 6), (1, 3), (2, 5))],
)
def test_rouse_initialization_matches_observed_pair_mean(missing_pairs):
    """Rouse k0 should use exactly the observed COV pair set."""
    n = 7
    true_spring_constant = 2.5
    target = _rouse_squared_distance_target(n, true_spring_constant)
    for i, j in missing_pairs:
        target[i, j] = np.nan
        target[j, i] = np.nan
    pair_i, pair_j = np.where(np.triu(np.isfinite(target), k=1))

    (
        connectivity,
        observed_pair_mean,
        unit_spring_rouse_pair_mean,
        spring_constant,
    ) = numerics._rouse_initial_connectivity(target, pair_i, pair_j)
    basis = numerics._centered_orthonormal_basis(n)
    reduced_gram, _ = numerics._reduced_gram_from_connectivity(
        connectivity, basis
    )
    gram = basis @ reduced_gram @ basis.T
    initialized = numerics._squared_distances_from_gram(gram)

    assert observed_pair_mean == pytest.approx(
        np.mean(target[pair_i, pair_j])
    )
    assert unit_spring_rouse_pair_mean == pytest.approx(
        3.0 * np.mean(pair_j - pair_i)
    )
    assert spring_constant == pytest.approx(true_spring_constant, rel=1e-12)
    assert np.mean(initialized[pair_i, pair_j]) == pytest.approx(
        observed_pair_mean, rel=1e-12
    )


@pytest.mark.parametrize("noise_model", ["homoskedastic", "heteroskedastic"])
def test_gaussian_covariance_initial_scalar_calibration_is_exact(noise_model):
    """Every initial shape should be optimally scaled along its Gram ray."""
    n = 5
    basis = numerics._centered_orthonormal_basis(n)
    reduced_gram = np.diag(np.linspace(0.7, 1.6, n - 1))
    reduced_gram += 0.05 * np.ones((n - 1, n - 1))
    gram = basis @ reduced_gram @ basis.T
    initial_distances = numerics._squared_distances_from_gram(gram)
    all_i, all_j = np.triu_indices(n, k=1)
    observed = np.arange(len(all_i)) % 3 != 0
    pair_i = all_i[observed]
    pair_j = all_j[observed]
    initial_pairs = initial_distances[pair_i, pair_j]
    target_pairs = initial_pairs * np.linspace(0.75, 1.25, len(pair_i))
    if noise_model == "homoskedastic":
        inverse_variance = np.full(len(pair_i), 2.5)
    else:
        inverse_variance = 1.0 / np.square(0.2 * target_pairs)

    scaled, calibration = (
        numerics._calibrate_gaussian_covariance_initial_scale(
            reduced_gram,
            basis,
            pair_i,
            pair_j,
            target_pairs,
            inverse_variance,
        )
    )

    a = np.sum(inverse_variance * np.square(initial_pairs))
    c = np.sum(inverse_variance * initial_pairs * target_pairs)
    expected_scale = (
        c + np.sqrt(c * c + 6.0 * (n - 1) * a)
    ) / (2.0 * a)
    target_matrix, inverse_variance_matrix = (
        numerics._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )

    def objective(candidate):
        return numerics._gaussian_covariance_objective_gradient(
            candidate,
            basis,
            target_matrix,
            inverse_variance_matrix,
            pair_i,
            pair_j,
        )[0]

    assert calibration["scale_factor"] == pytest.approx(
        expected_scale, rel=1e-14
    )
    assert np.allclose(scaled, expected_scale * reduced_gram)
    assert calibration["objective_before"] == pytest.approx(
        objective(reduced_gram), rel=1e-13, abs=1e-13
    )
    assert calibration["objective_after"] == pytest.approx(
        objective(scaled), rel=1e-13, abs=1e-13
    )
    assert calibration["objective_after"] < objective(0.9 * scaled)
    assert calibration["objective_after"] < objective(1.1 * scaled)
    assert calibration["objective_reduction"] > 0.0
    assert calibration["relative_derivative_residual"] <= 1e-14
    assert calibration["backend"] == "cpu"
    assert calibration["wall_seconds"] >= 0.0


def test_gaussian_covariance_relative_noise_uses_rouse_default_and_stationarity():
    """Relative-noise COV should converge from the calibrated Rouse start."""
    n = 6
    target = _rouse_squared_distance_target(n)
    relative_std = 0.05

    fitted, gram, connectivity, info = (
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            relative_noise_std=relative_std,
            max_iterations=30,
            relative_tolerance=1e-8,
            save_steps=[1],
        )
    )

    variance = np.square(relative_std * target)
    upper = np.triu_indices(n, k=1)
    assert info["converged"]
    assert info["noise_model"] == "heteroskedastic_relative_std"
    assert info["initialization"]["kind"] == "rouse"
    assert info["initialization"]["backend"] == "cpu"
    assert info["initialization"]["spring_constant"] == pytest.approx(1.0)
    scalar_calibration = info["initialization"]["scalar_calibration"]
    assert scalar_calibration["scale_factor"] > 0.0
    assert scalar_calibration["objective_after"] < scalar_calibration[
        "objective_before"
    ]
    assert scalar_calibration["relative_derivative_residual"] <= 1e-14
    assert info["initialization"]["effective_spring_constant"] == pytest.approx(
        scalar_calibration["connectivity_scale_factor"]
    )
    assert info["initialization"][
        "observed_pair_mean_squared_distance"
    ] == pytest.approx(np.mean(target[upper]))
    assert info["backend"] == "cpu"
    assert info["gpu_device"] is None
    assert info["preconditioner_data_setup_count"] == 1
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-10)
    assert np.count_nonzero(np.linalg.eigvalsh(gram) > 1e-12) == n - 1
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    assert np.allclose(
        (fitted - target)[upper],
        0.5 * variance[upper] * connectivity[upper],
        rtol=1e-7,
        atol=1e-9,
    )
    assert info["relative_stationarity_residual"] <= 1e-8
    assert sorted(info["connectivity_matrix_at_steps"]) == [1]


def test_gaussian_covariance_absolute_noise_initializations_match():
    """All scalar-calibrated starts should reach the homoskedastic optimum."""
    target = _rouse_squared_distance_target()
    variance = 0.02
    solver_options = {
        "max_iterations": 40,
        "relative_tolerance": 2e-8,
    }

    rouse = HippsDimes.fit_gaussian_noise_covariance(
        target,
        variance,
        **solver_options,
    )
    nearest = HippsDimes.fit_gaussian_noise_covariance(
        target,
        variance,
        initialization="nearest_edm",
        **solver_options,
    )
    provided = HippsDimes.fit_gaussian_noise_covariance(
        target,
        variance,
        initial_connectivity=HippsDimes.construct_connectivity_matrix_rouse(
            len(target), 0.4
        ),
        **solver_options,
    )

    for result in (rouse, nearest, provided):
        assert result[3]["converged"]
        scalar_calibration = result[3]["initialization"]["scalar_calibration"]
        assert scalar_calibration["scale_factor"] > 0.0
        assert scalar_calibration["objective_after"] <= scalar_calibration[
            "objective_before"
        ]
    assert rouse[3]["noise_model"] == "homoskedastic_absolute_variance"
    assert nearest[3]["initialization"]["kind"] == "weighted_nearest_edm"
    assert provided[3]["initialization"]["kind"] == "provided_connectivity"
    assert nearest[3]["initialization"]["backend"] == "cpu"
    assert nearest[3]["initialization"]["nearest_edm_wall_seconds"] >= 0.0
    assert nearest[3]["initialization"]["nearest_edm_projection_count"] > 0
    assert nearest[3]["initialization"]["wall_seconds"] >= nearest[3][
        "initialization"
    ]["nearest_edm_wall_seconds"]
    assert np.allclose(rouse[0], nearest[0], rtol=1e-8, atol=1e-9)
    assert np.allclose(rouse[2], nearest[2], rtol=1e-8, atol=1e-9)
    assert np.allclose(rouse[0], provided[0], rtol=1e-8, atol=1e-9)
    assert np.allclose(rouse[2], provided[2], rtol=1e-8, atol=1e-9)


def test_gaussian_covariance_relative_noise_supports_missing_pairs():
    target = _rouse_squared_distance_target()
    missing_pairs = ((0, 3), (1, 4))
    for i, j in missing_pairs:
        target[i, j] = np.nan
        target[j, i] = np.nan

    fitted, _, connectivity, info = HippsDimes.fit_gaussian_noise_covariance(
        target,
        relative_noise_std=0.05,
        max_iterations=40,
        relative_tolerance=1e-8,
    )

    assert info["converged"]
    assert info["observed_pair_count"] == 13
    assert info["initialization"]["scalar_calibration"]["scale_factor"] > 0.0
    assert np.all(np.isfinite(fitted))
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12


@pytest.mark.parametrize(
    ("noise_variance", "relative_noise_std", "message"),
    [
        (None, None, "exactly one"),
        (0.1, 0.2, "exactly one"),
        (0.0, None, "positive finite scalar"),
        (None, -0.1, "positive finite scalar"),
        (np.ones((2, 2)), None, "positive finite scalar"),
        (None, 1e-300, "produce positive finite pair variances"),
    ],
)
def test_gaussian_covariance_rejects_invalid_noise_contract(
    noise_variance, relative_noise_std, message
):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match=message):
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            noise_variance,
            relative_noise_std=relative_noise_std,
        )


def test_gaussian_covariance_rejects_nearest_edm_with_connectivity():
    target = _rouse_squared_distance_target(4)
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(4, 1.0)

    with pytest.raises(ValueError, match="cannot be combined"):
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            0.1,
            initialization="nearest_edm",
            initial_connectivity=connectivity,
        )


def test_gaussian_covariance_gpu_request_fails_without_available_gpu(monkeypatch):
    target = _rouse_squared_distance_target(4)
    monkeypatch.setattr(numerics, "_CUPY_AVAILABLE", False)

    with pytest.raises(RuntimeError, match="CuPy.*accessible CUDA GPU"):
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            0.1,
            use_gpu=True,
        )


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
def test_gaussian_covariance_gpu_kernels_match_cpu(working_cupy):
    """Float64 GPU COV kernels and block preconditioner should match CPU."""
    cp = working_cupy
    rng = np.random.default_rng(20260827)
    n = 7
    basis = numerics._centered_orthonormal_basis(n)
    all_i, all_j = np.triu_indices(n, k=1)
    selector = np.arange(len(all_i)) % 4 != 0
    pair_i = all_i[selector]
    pair_j = all_j[selector]
    reduced_gram = np.diag(np.linspace(0.8, 2.0, n - 1))
    pair_vectors = basis[pair_i] - basis[pair_j]
    target_pairs = np.einsum(
        "ij,ij->i", pair_vectors @ reduced_gram, pair_vectors, optimize=True
    )
    target_pairs *= np.linspace(0.95, 1.05, len(target_pairs))
    inverse_variance = np.linspace(1.5, 4.5, len(pair_i))
    target_matrix, inverse_variance_matrix = (
        numerics._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    direction = rng.normal(size=(n - 1, n - 1))
    direction = 0.5 * (direction + direction.T)

    cpu_state = numerics._gaussian_covariance_objective_gradient(
        reduced_gram,
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
    )
    gpu_state = numerics._gaussian_covariance_objective_gradient(
        cp.asarray(reduced_gram),
        cp.asarray(basis),
        cp.asarray(target_matrix),
        cp.asarray(inverse_variance_matrix),
        cp.asarray(pair_i, dtype=cp.int32),
        cp.asarray(pair_j, dtype=cp.int32),
        array_module=cp,
    )
    cpu_hessian = numerics._gaussian_covariance_data_hessian_action(
        direction, basis, inverse_variance_matrix
    )
    gpu_hessian = numerics._gaussian_covariance_data_hessian_action(
        cp.asarray(direction),
        cp.asarray(basis),
        cp.asarray(inverse_variance_matrix),
        array_module=cp,
    )
    cpu_diagonal, _ = (
        numerics._gaussian_covariance_data_preconditioner_diagonal(
            basis,
            pair_i,
            pair_j,
            inverse_variance,
            block_size=3,
        )
    )
    gpu_diagonal, _ = (
        numerics._gaussian_covariance_data_preconditioner_diagonal(
            cp.asarray(basis),
            cp.asarray(pair_i, dtype=cp.int32),
            cp.asarray(pair_j, dtype=cp.int32),
            cp.asarray(inverse_variance),
            block_size=3,
            array_module=cp,
            use_gpu=True,
        )
    )

    assert gpu_state[0] == pytest.approx(cpu_state[0], rel=1e-13, abs=1e-13)
    assert np.allclose(cp.asnumpy(gpu_state[1]), cpu_state[1], rtol=1e-12, atol=1e-12)
    assert np.allclose(cp.asnumpy(gpu_state[2]), cpu_state[2], rtol=1e-12, atol=1e-12)
    assert np.allclose(cp.asnumpy(gpu_state[3]), cpu_state[3], rtol=1e-12, atol=1e-12)
    assert np.allclose(cp.asnumpy(gpu_hessian), cpu_hessian, rtol=1e-12, atol=1e-12)
    assert np.allclose(cp.asnumpy(gpu_diagonal), cpu_diagonal, rtol=1e-12, atol=1e-12)


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
@pytest.mark.parametrize(
    "noise_options",
    [
        {"relative_noise_std": 0.05},
        {"noise_variance": 0.02, "initialization": "nearest_edm"},
    ],
)
def test_gaussian_covariance_gpu_matches_cpu_end_to_end(
    noise_options, working_cupy
):
    """GPU COV should match CPU for both noise models and both initializers."""
    target = _rouse_squared_distance_target(6)
    solver_options = {
        "max_iterations": 30,
        "relative_tolerance": 2e-8,
        "save_steps": [1],
        **noise_options,
    }

    cpu = HippsDimes.fit_gaussian_noise_covariance(target, **solver_options)
    progress_events = []
    gpu = HippsDimes.fit_gaussian_noise_covariance(
        target,
        use_gpu=True,
        progress_callback=progress_events.append,
        **solver_options,
    )

    assert cpu[3]["converged"]
    assert gpu[3]["converged"]
    assert gpu[3]["backend"] == "gpu"
    assert gpu[3]["dtype"] == "float64"
    assert gpu[3]["gpu_device"] == numerics.get_gpu_name()
    assert gpu[3]["cupy_version"] == numerics.cp.__version__
    assert gpu[3]["preconditioner_setup_seconds"] >= 0.0
    assert sorted(gpu[3]["connectivity_matrix_at_steps"]) == [1]
    assert np.allclose(gpu[0], cpu[0], rtol=1e-8, atol=1e-9)
    assert np.allclose(gpu[1], cpu[1], rtol=1e-8, atol=1e-9)
    assert np.allclose(gpu[2], cpu[2], rtol=1e-8, atol=1e-9)
    assert gpu[3]["relative_stationarity_residual"] <= 1e-8
    cpu_calibration = cpu[3]["initialization"]["scalar_calibration"]
    gpu_calibration = gpu[3]["initialization"]["scalar_calibration"]
    assert cpu_calibration["backend"] == "cpu"
    assert gpu_calibration["backend"] == "gpu"
    assert gpu_calibration["scale_factor"] == pytest.approx(
        cpu_calibration["scale_factor"], rel=1e-12
    )
    assert gpu_calibration["objective_after"] == pytest.approx(
        cpu_calibration["objective_after"], rel=1e-12, abs=1e-12
    )
    expected_stages = {
        "covariance_preconditioner",
        "covariance_optimization",
    }
    if noise_options.get("initialization") == "nearest_edm":
        expected_stages.add("nearest_edm_initialization")
        assert gpu[3]["initialization"]["backend"] == "gpu"
        assert gpu[3]["initialization"]["gpu_device"] == numerics.get_gpu_name()
        assert gpu[3]["initialization"]["nearest_edm_wall_seconds"] >= 0.0
        assert gpu[3]["initialization"]["nearest_edm_projection_count"] > 0
    else:
        assert gpu[3]["initialization"]["backend"] == "cpu"
    assert {event["stage"] for event in progress_events} == expected_stages
    assert all(event["use_gpu"] for event in progress_events)


def test_compute_modulus():
    """Moduli should use tau_p / 2 for stress-mode relaxation."""
    A = HippsDimes.construct_connectivity_matrix_rouse(2, 1.0)
    zeta = 1.0

    internal_eigenvalue = np.linalg.eigvalsh(A)[0]
    tau_p = -zeta / internal_eigenvalue
    freq = np.array([2.0 / tau_p])

    storage_mod, loss_mod = HippsDimes.compute_modulus(A, freq, zeta)

    # Check shapes
    assert storage_mod.shape == (len(freq), 2)
    assert loss_mod.shape == (len(freq), 2)

    # Check that frequencies match
    assert np.allclose(storage_mod[:, 0], freq)
    assert np.allclose(loss_mod[:, 0], freq)

    # At omega * tau_p / 2 = 1, one unnormalized Maxwell mode has G' = G'' = 1/2.
    assert storage_mod[0, 1] == pytest.approx(0.5)
    assert loss_mod[0, 1] == pytest.approx(0.5)


def test_compute_monomer_mechanical_susceptibility():
    """Monomer susceptibility should use tau_p directly with a 1/zeta prefactor."""
    A = HippsDimes.construct_connectivity_matrix_rouse(2, 1.0)
    zeta = 2.0
    internal_eigenvalue = np.linalg.eigvalsh(A)[0]
    tau_p = -zeta / internal_eigenvalue
    freq = np.array([0.0, 1.0 / tau_p, 2.0 / tau_p])

    freq_out, chi_prime, chi_double_prime = (
        HippsDimes.compute_monomer_mechanical_susceptibility(A, freq, zeta)
    )

    mode_weight = 0.5
    denominator = 1.0 + (freq * tau_p) ** 2
    expected_prime = mode_weight * tau_p / (zeta * denominator)
    expected_double_prime = (
        mode_weight * freq * tau_p**2 / (zeta * denominator)
    )
    expected_prime = np.column_stack((expected_prime, expected_prime))
    expected_double_prime = np.column_stack(
        (expected_double_prime, expected_double_prime)
    )

    assert np.allclose(freq_out, freq)
    assert np.allclose(chi_prime, expected_prime)
    assert np.allclose(chi_double_prime, expected_double_prime)
    assert chi_prime[1, 0] == pytest.approx(chi_double_prime[1, 0])
    assert (
        hipps_dimes.compute_monomer_mechanical_susceptibility
        is HippsDimes.compute_monomer_mechanical_susceptibility
    )
    assert not hasattr(HippsDimes, "compute_monomer_modulus")
    assert not hasattr(hipps_dimes, "compute_monomer_modulus")


def test_compute_stress_relaxation():
    """Stress relaxation should normalize the non-zero mode sum by chain length."""
    A = HippsDimes.construct_connectivity_matrix_rouse(4, 1.0)
    t = [0.0, 0.1, 1.0, 10.0]

    relaxation = HippsDimes.compute_stress_relaxation(A, t, zeta=1.0)

    eigvals = np.linalg.eigvalsh(A)
    lam_nz = eigvals[np.abs(eigvals) > 1e-12]
    expected = np.sum(
        np.exp(-2.0 * np.asarray(t)[:, None] / (-1.0 / lam_nz)[None, :]),
        axis=1,
    ) / len(A)

    assert relaxation.shape == (len(t), 2)
    assert np.allclose(relaxation[:, 0], t)
    assert relaxation[0, 1] == pytest.approx(len(lam_nz) / len(A))
    assert np.all(np.diff(relaxation[:, 1]) <= 1e-12)
    assert np.allclose(relaxation[:, 1], expected)


def test_compute_stress_relaxation_stretched_exponent():
    """Stress relaxation should support stretched-exponential mode decay."""
    A = HippsDimes.construct_connectivity_matrix_rouse(4, 1.0)
    t = np.array([0.0, 0.1, 1.0, 10.0])
    stretched_exponent = 0.5

    relaxation = HippsDimes.compute_stress_relaxation(
        A,
        t,
        zeta=1.0,
        stretched_exponent=stretched_exponent,
    )

    eigvals = np.linalg.eigvalsh(A)
    lam_nz = eigvals[np.abs(eigvals) > 1e-12]
    tau_p = -1.0 / lam_nz
    expected = np.sum(
        np.exp(-2.0 * np.power(t[:, None] / tau_p[None, :], stretched_exponent)),
        axis=1,
    ) / len(A)

    assert relaxation.shape == (len(t), 2)
    assert np.allclose(relaxation[:, 0], t)
    assert relaxation[0, 1] == pytest.approx(len(lam_nz) / len(A))
    assert np.all(np.diff(relaxation[:, 1]) <= 1e-12)
    assert np.allclose(relaxation[:, 1], expected)


def test_dynamics_pairwise_distances():
    """Dynamics should compute pairwise Euclidean distances correctly."""
    dyn = HippsDimes.Dynamics(3, k=1.0, model="rouse")
    dyn.xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
        ]
    )
    distances = dyn._calculate_pairwise_distances()

    expected = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, np.sqrt(5.0)],
            [2.0, np.sqrt(5.0), 0.0],
        ]
    )
    assert np.allclose(distances, expected)


def test_a2dmap_theory_with_force_applied_smoke():
    """Force-applied dmap computation should return a valid symmetric map."""
    A = HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    force = {"loci": [0, 4], "amplitude": 1.0, "direction": [1.0, 0.0, 0.0]}
    dmap = HippsDimes.a2dmap_theory_with_force_applied(A, force)

    assert dmap.shape == (5, 5)
    assert np.allclose(dmap, dmap.T)
    assert np.allclose(np.diag(dmap), 0.0)


def test_restore_matrix_with_nans_reinserts_removed_loci():
    """Reduced matrices should be restored with NaN-filled rows/cols at removed loci."""
    small = np.array(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 3.0],
            [2.0, 3.0, 0.0],
        ]
    )

    restored = HippsDimes.restore_matrix_with_nans(small, removed_idx=[1, 3], original_size=5)

    expected = np.array(
        [
            [0.0, np.nan, 1.0, np.nan, 2.0],
            [np.nan, np.nan, np.nan, np.nan, np.nan],
            [1.0, np.nan, 0.0, np.nan, 3.0],
            [np.nan, np.nan, np.nan, np.nan, np.nan],
            [2.0, np.nan, 3.0, np.nan, 0.0],
        ]
    )

    assert restored.shape == (5, 5)
    assert np.array_equal(np.isnan(restored), np.isnan(expected))
    assert np.allclose(np.nan_to_num(restored), np.nan_to_num(expected))


def test_restore_matrix_with_nans_validates_shape():
    """Restoration should reject reduced matrices whose shape is inconsistent with removed_idx."""
    small = np.eye(3)

    with pytest.raises(ValueError, match="incompatible shape"):
        HippsDimes.restore_matrix_with_nans(small, removed_idx=[1], original_size=5)


def test_a2xyz_sample_conditioned_pair_distance_enforces_pair_norm():
    """Conditioned ensemble sampling should enforce the requested pair distance for each sample."""
    np.random.seed(12345)
    A = HippsDimes.construct_connectivity_matrix_rouse(6, 1.0)
    pair = (1, 4)
    target_distance = 2.5

    xyzs = HippsDimes.a2xyz_sample_conditioned_pair_distance(
        A,
        pair=pair,
        b_scalar=target_distance,
        ensemble=24,
    )

    pair_vectors = xyzs[:, pair[0], :] - xyzs[:, pair[1], :]
    pair_distances = np.linalg.norm(pair_vectors, axis=1)

    assert xyzs.shape == (24, 6, 3)
    assert np.allclose(pair_distances, target_distance, rtol=1e-10, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__]) 
