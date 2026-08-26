"""Basic tests for HippsDimes package."""

import numpy as np
import pytest

# Import the main module
import HippsDimes
import hipps_dimes
import hipps_dimes.numerics as numerics


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
def test_a2a_enforces_nonnegative_spring_constants_on_gpu():
    """GPU normalization should preserve the backend and match the CPU result."""
    source_cpu = np.array(
        [
            [-9.0, 2.0, -3.0],
            [2.0, -9.0, 4.0],
            [-3.0, 4.0, -9.0],
        ]
    )
    source_gpu = numerics.cp.asarray(source_cpu)
    try:
        source_gpu.copy()
    except RuntimeError as error:
        pytest.skip(f"CuPy kernel execution is unavailable: {error}")

    result_gpu = HippsDimes.a2a(source_gpu, fill_negative=True)

    expected = HippsDimes.a2a(source_cpu, fill_negative=True)
    assert isinstance(result_gpu, numerics.cp.ndarray)
    assert np.array_equal(numerics.cp.asnumpy(source_gpu), source_cpu)
    assert np.array_equal(numerics.cp.asnumpy(result_gpu), expected)


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
    J = np.eye(n) - np.ones((n, n)) / n

    projected = numerics._project_centered_psd(source, floor)
    projected_twice = numerics._project_centered_psd(projected, floor)
    eigenvalues = np.linalg.eigvalsh(projected)

    assert np.allclose(projected, projected.T, atol=1e-12)
    assert np.allclose(np.sum(projected, axis=1), 0.0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(projected - floor * J)) >= -1e-12
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
    B = numerics._project_centered_psd(source)
    direction = np.array(
        [
            [0.2, -0.5, 0.7],
            [-0.5, 0.4, 0.1],
            [0.7, 0.1, -0.3],
        ]
    )

    _, gradient, _ = numerics._nearest_edm_objective_gradient(B, target, weights)
    epsilon = 1e-6
    forward, _, _ = numerics._nearest_edm_objective_gradient(
        B + epsilon * direction, target, weights
    )
    backward, _, _ = numerics._nearest_edm_objective_gradient(
        B - epsilon * direction, target, weights
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

    fitted, gram, info = hipps_dimes.nearest_edm(target)

    assert info["converged"]
    assert info["status"] == "optimality_tolerance"
    assert np.allclose(fitted, expected, rtol=1e-8, atol=1e-9)
    assert np.allclose(np.sum(gram, axis=1), 0.0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(gram)) >= -1e-12
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-12)
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
    J = np.eye(n) - np.ones((n, n)) / n
    expected_gram = floor * J
    expected_fitted = np.full((n, n), 2.0 * floor)
    np.fill_diagonal(expected_fitted, 0.0)

    fitted, gram, info = hipps_dimes.nearest_edm(target, gram_eigenvalue_floor=floor)

    assert info["converged"]
    assert info["gram_eigenvalue_floor"] == pytest.approx(floor)
    assert np.allclose(gram, expected_gram, rtol=1e-11, atol=1e-12)
    assert np.allclose(fitted, expected_fitted, rtol=1e-11, atol=1e-12)
    assert np.count_nonzero(np.linalg.eigvalsh(gram) > 0.5 * floor) == n - 1


def test_nearest_edm_floor_leaves_strictly_feasible_target_unchanged():
    """A target whose Gram spectrum clears the floor should retain zero loss."""
    n = 4
    floor = 0.1
    J = np.eye(n) - np.ones((n, n)) / n
    source_gram = 2.0 * floor * J
    target = numerics._squared_distances_from_gram(source_gram)

    fitted, gram, info = hipps_dimes.nearest_edm(target, gram_eigenvalue_floor=floor)

    assert info["converged"]
    assert info["objective"] == pytest.approx(0.0, abs=1e-20)
    assert np.allclose(fitted, target, rtol=1e-10, atol=1e-11)
    assert np.allclose(gram, source_gram, rtol=1e-10, atol=1e-11)


def test_nearest_edm_zero_floor_is_backward_compatible():
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )

    fitted_default, gram_default, info_default = hipps_dimes.nearest_edm(target)
    fitted_zero, gram_zero, info_zero = hipps_dimes.nearest_edm(
        target, gram_eigenvalue_floor=0.0
    )

    assert np.array_equal(fitted_zero, fitted_default)
    assert np.array_equal(gram_zero, gram_default)
    assert info_default["gram_eigenvalue_floor"] == 0.0
    assert info_zero["gram_eigenvalue_floor"] == 0.0


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
    expected_leg = (leg_weight + 18.0 * base_weight) / (leg_weight + 8.0 * base_weight)
    expected = np.array(
        [
            [0.0, expected_leg, expected_leg],
            [expected_leg, 0.0, 4.0 * expected_leg],
            [expected_leg, 4.0 * expected_leg, 0.0],
        ]
    )

    fitted, _, info = hipps_dimes.nearest_edm(target, weights)
    scaled_fitted, _, scaled_info = hipps_dimes.nearest_edm(target, 7.0 * weights)

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


def test_nearest_edm_floor_supports_weights_and_missing_pairs():
    """The floored feasible set should retain weighted/missing-data semantics."""
    n = 4
    floor = 0.1
    J = np.eye(n) - np.ones((n, n)) / n
    target = numerics._squared_distances_from_gram(0.4 * J)
    target[0, 3] = np.nan
    target[3, 0] = np.nan
    weights = np.ones((n, n)) - np.eye(n)
    weights[0, 3] = 0.0
    weights[3, 0] = 0.0
    weights[1, 2] = 3.0
    weights[2, 1] = 3.0
    observed = weights > 0.0

    fitted, gram, info = hipps_dimes.nearest_edm(
        target,
        weights,
        gram_eigenvalue_floor=floor,
        absolute_tolerance=1e-8,
    )

    assert info["converged"]
    assert np.allclose(fitted[observed], target[observed], atol=1e-7)
    assert np.min(np.linalg.eigvalsh(gram - floor * J)) >= -1e-12
    assert np.allclose(np.sum(gram, axis=1), 0.0, atol=1e-12)


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


@pytest.mark.parametrize(
    "floor",
    [-0.1, np.nan, np.inf, True, [0.1]],
)
def test_nearest_edm_rejects_invalid_gram_eigenvalue_floor(floor):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="finite nonnegative scalar"):
        hipps_dimes.nearest_edm(target, gram_eigenvalue_floor=floor)


def test_centered_orthonormal_basis_has_expected_geometry():
    basis = numerics._centered_orthonormal_basis(7)

    assert basis.shape == (7, 6)
    assert np.allclose(basis.T @ basis, np.eye(6), atol=1e-14)
    assert np.allclose(np.sum(basis, axis=0), 0.0, atol=1e-14)


def test_rouse_kl_prox_is_positive_and_satisfies_scalar_optimality():
    trial = np.diag([-2.0, 0.5, 3.0])
    reference_inverse = np.diag([0.25, 2.0, 4.0])
    step_size = 0.2
    coefficient = 0.3

    proximal = numerics._rouse_kl_prox(trial, step_size, coefficient, reference_inverse)
    eigenvalues = np.diag(proximal)
    residual = (
        eigenvalues
        - np.diag(trial)
        + step_size * coefficient * (np.diag(reference_inverse) - 1.0 / eigenvalues)
    )

    assert np.all(eigenvalues > 0.0)
    assert np.allclose(residual, 0.0, atol=1e-13)


def test_nearest_edm_rouse_prior_recovers_exact_rouse_reference():
    n = 6
    spring_constant = 2.0
    reference, _, _ = numerics._rouse_reference_gram(n, spring_constant)
    target = numerics._squared_distances_from_gram(reference)

    fitted, gram, info = hipps_dimes.nearest_edm(
        target,
        rouse_prior_weight=0.1,
        rouse_spring_constant=spring_constant,
    )

    assert info["converged"]
    assert info["mean_rouse_kl"] == pytest.approx(0.0, abs=1e-13)
    assert info["normalized_data_objective"] == pytest.approx(0.0, abs=1e-26)
    assert np.allclose(fitted, target, atol=1e-13)
    assert np.allclose(gram, reference, atol=1e-13)
    assert np.all(np.linalg.eigvalsh(gram)[1:] > 0.0)


def test_nearest_edm_stronger_rouse_prior_moves_fit_toward_reference():
    target = np.array(
        [
            [0.0, 1.0, 1.0, 2.0],
            [1.0, 0.0, 9.0, 3.0],
            [1.0, 9.0, 0.0, 1.0],
            [2.0, 3.0, 1.0, 0.0],
        ]
    )
    weights = np.zeros_like(target)
    off_diagonal = ~np.eye(len(target), dtype=bool)
    weights[off_diagonal] = 1.0 / np.square(target[off_diagonal])

    _, weak_gram, weak = hipps_dimes.nearest_edm(
        target,
        weights,
        rouse_prior_weight=1e-4,
        rouse_spring_constant=3.0,
        max_iterations=1000,
        relative_tolerance=1e-7,
    )
    _, strong_gram, strong = hipps_dimes.nearest_edm(
        target,
        weights,
        rouse_prior_weight=1.0,
        rouse_spring_constant=3.0,
        max_iterations=1000,
        relative_tolerance=1e-7,
    )

    assert weak["converged"] and strong["converged"]
    assert strong["mean_rouse_kl"] < weak["mean_rouse_kl"]
    assert strong["normalized_data_objective"] > weak["normalized_data_objective"]
    assert np.all(np.linalg.eigvalsh(weak_gram)[1:] > 0.0)
    assert np.all(np.linalg.eigvalsh(strong_gram)[1:] > 0.0)
    assert np.all(np.diff(weak["history"]["total_objective"]) <= 1e-12)
    assert np.all(np.diff(strong["history"]["total_objective"]) <= 1e-12)


def test_nearest_edm_rouse_prior_supports_missing_pairs_and_infers_scale():
    target = np.array(
        [
            [0.0, 1.0, 4.0, np.nan],
            [1.0, 0.0, 2.0, 5.0],
            [4.0, 2.0, 0.0, 3.0],
            [np.nan, 5.0, 3.0, 0.0],
        ]
    )

    fitted, gram, info = hipps_dimes.nearest_edm(
        target,
        rouse_prior_weight=0.01,
        max_iterations=1000,
        relative_tolerance=1e-6,
    )

    assert info["converged"]
    assert info["observed_pair_count"] == 5
    assert info["rouse_spring_constant"] == pytest.approx(3.0 / 2.0)
    assert np.isfinite(fitted).all()
    assert np.all(np.linalg.eigvalsh(gram)[1:] > 0.0)
    assert np.allclose(np.sum(gram, axis=1), 0.0, atol=1e-12)


def test_nearest_edm_rouse_prior_relative_objective_is_scale_equivariant():
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )
    off_diagonal = ~np.eye(len(target), dtype=bool)
    weights = np.zeros_like(target)
    weights[off_diagonal] = 1.0 / np.square(target[off_diagonal])
    options = {
        "rouse_prior_weight": 0.1,
        "rouse_spring_constant": 3.0,
        "max_iterations": 1000,
        "relative_tolerance": 1e-7,
    }

    fitted, gram, info = hipps_dimes.nearest_edm(target, weights, **options)
    scale = 7.0
    scaled_fitted, scaled_gram, scaled_info = hipps_dimes.nearest_edm(
        scale * target,
        weights / scale**2,
        **{
            **options,
            "rouse_spring_constant": options["rouse_spring_constant"] / scale,
        },
    )

    assert info["converged"] and scaled_info["converged"]
    assert np.allclose(scaled_fitted / scale, fitted, rtol=1e-8, atol=1e-9)
    assert np.allclose(scaled_gram / scale, gram, rtol=1e-8, atol=1e-9)
    assert scaled_info["total_objective"] == pytest.approx(
        info["total_objective"], rel=1e-10
    )


def test_nearest_edm_zero_rouse_prior_is_bitwise_backward_compatible():
    target = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 9.0],
            [1.0, 9.0, 0.0],
        ]
    )

    default = hipps_dimes.nearest_edm(target)
    explicit = hipps_dimes.nearest_edm(
        target, rouse_prior_weight=0.0, rouse_spring_constant=2.0
    )

    assert np.array_equal(default[0], explicit[0])
    assert np.array_equal(default[1], explicit[1])
    assert default[2].keys() == explicit[2].keys()
    for key in default[2]["history"]:
        assert np.array_equal(default[2]["history"][key], explicit[2]["history"][key])


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"rouse_prior_weight": -0.1}, "finite nonnegative"),
        ({"rouse_prior_weight": np.nan}, "finite nonnegative"),
        (
            {"rouse_prior_weight": 0.1, "rouse_spring_constant": 0.0},
            "finite positive",
        ),
        (
            {"rouse_prior_weight": 0.1, "gram_eigenvalue_floor": 0.1},
            "cannot both be positive",
        ),
    ],
)
def test_nearest_edm_rejects_invalid_rouse_prior_options(options, message):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match=message):
        hipps_dimes.nearest_edm(target, **options)


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


def test_gaussian_covariance_objective_gradient_matches_finite_difference():
    """The exact soft-constraint objective should expose its true matrix gradient."""
    n = 5
    basis = numerics._centered_orthonormal_basis(n)
    upper = np.triu_indices(n, k=1)
    pair_vectors = basis[upper[0]] - basis[upper[1]]
    reduced_gram = np.diag(np.linspace(0.8, 2.0, n - 1))
    target = np.einsum(
        "ij,ij->i", pair_vectors @ reduced_gram, pair_vectors, optimize=True
    )
    target *= np.linspace(0.95, 1.05, len(target))
    inverse_variance = np.linspace(2.0, 5.0, len(target))
    target_matrix, inverse_variance_matrix = (
        numerics._gaussian_covariance_pair_matrices(
            n, upper[0], upper[1], target, inverse_variance
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
        upper[0],
        upper[1],
    )
    step = 1e-6
    forward = numerics._gaussian_covariance_objective_gradient(
        reduced_gram + step * direction,
        basis,
        target_matrix,
        inverse_variance_matrix,
        upper[0],
        upper[1],
    )[0]
    backward = numerics._gaussian_covariance_objective_gradient(
        reduced_gram - step * direction,
        basis,
        target_matrix,
        inverse_variance_matrix,
        upper[0],
        upper[1],
    )[0]

    assert np.isfinite(objective)
    assert (forward - backward) / (2.0 * step) == pytest.approx(
        np.sum(gradient * direction), rel=1e-7, abs=1e-8
    )


def test_structured_gaussian_covariance_operators_match_pair_vectors():
    """Structured EDM kernels should equal the explicit pair-vector formulas."""
    rng = np.random.default_rng(20260826)
    n = 7
    basis = numerics._centered_orthonormal_basis(n)
    upper = np.triu_indices(n, k=1)
    reduced_gram = np.diag(np.linspace(0.7, 1.9, n - 1))
    reduced_gram += 0.03 * np.ones_like(reduced_gram)
    direction = rng.normal(size=(n - 1, n - 1))
    direction = 0.5 * (direction + direction.T)

    cases = (
        (np.ones(len(upper[0]), dtype=bool), 0.0, False),
        (np.arange(len(upper[0])) % 4 != 0, 0.0, True),
        (np.ones(len(upper[0]), dtype=bool), 0.04, True),
    )
    for selector, connectivity_l1, heteroscedastic in cases:
        pair_i = upper[0][selector]
        pair_j = upper[1][selector]
        pair_vectors = basis[pair_i] - basis[pair_j]
        fitted_reference = np.einsum(
            "ij,ij->i",
            pair_vectors @ reduced_gram,
            pair_vectors,
            optimize=True,
        )
        target_pairs = fitted_reference * np.linspace(
            0.94, 1.06, len(fitted_reference)
        )
        inverse_variance = (
            np.linspace(1.5, 4.5, len(target_pairs))
            if heteroscedastic
            else np.full(len(target_pairs), 2.5)
        )
        target_matrix, inverse_variance_matrix = (
            numerics._gaussian_covariance_pair_matrices(
                n,
                pair_i,
                pair_j,
                target_pairs,
                inverse_variance,
            )
        )

        residual_reference = fitted_reference - target_pairs
        effective_reference = (
            residual_reference
            if connectivity_l1 == 0.0
            else numerics._soft_threshold(
                residual_reference, 2.0 * connectivity_l1
            )
        )
        weighted_reference = inverse_variance * effective_reference
        objective_reference = 0.5 * np.dot(
            weighted_reference, effective_reference
        )
        gradient_reference = pair_vectors.T @ (
            weighted_reference[:, np.newaxis] * pair_vectors
        )
        gradient_reference = 0.5 * (
            gradient_reference + gradient_reference.T
        )

        objective, gradient, fitted, residual = (
            numerics._gaussian_covariance_data_objective_gradient(
                reduced_gram,
                basis,
                target_matrix,
                inverse_variance_matrix,
                pair_i,
                pair_j,
                connectivity_l1,
            )
        )
        assert objective == pytest.approx(
            objective_reference, rel=1e-12, abs=1e-13
        )
        assert np.allclose(gradient, gradient_reference, rtol=1e-12, atol=1e-12)
        assert np.allclose(fitted, fitted_reference, rtol=1e-12, atol=1e-12)
        assert np.allclose(residual, residual_reference, rtol=1e-12, atol=1e-12)

        pair_quadratic = np.einsum(
            "ij,ij->i", pair_vectors @ direction, pair_vectors, optimize=True
        )
        hessian_reference = pair_vectors.T @ (
            (inverse_variance * pair_quadratic)[:, np.newaxis] * pair_vectors
        )
        hessian_reference = 0.5 * (
            hessian_reference + hessian_reference.T
        )
        hessian_action = numerics._gaussian_covariance_data_hessian_action(
            direction, basis, inverse_variance_matrix
        )
        assert np.allclose(
            hessian_action, hessian_reference, rtol=1e-12, atol=1e-12
        )


def test_gaussian_negative_entropy_proximal_map_is_positive_and_exact():
    """The log-determinant proximal eigenvalues should solve their quadratics."""
    trial_eigenvalues = np.array([-1e12, -2.0, 0.0, 3.0])
    trial = np.diag(trial_eigenvalues)
    step_size = 0.2

    proximal, observed_eigenvalues, updated_eigenvalues = (
        numerics._proximal_gaussian_negative_entropy(trial, step_size)
    )

    assert np.array_equal(observed_eigenvalues, trial_eigenvalues)
    assert np.all(updated_eigenvalues > 0.0)
    assert np.allclose(np.linalg.eigvalsh(proximal), updated_eigenvalues)
    assert np.allclose(
        np.square(updated_eigenvalues)
        - observed_eigenvalues * updated_eigenvalues,
        1.5 * step_size,
        rtol=1e-13,
        atol=1e-14,
    )


def test_gaussian_covariance_fit_stays_physical_and_satisfies_stationarity():
    """The COV solver should converge inside the cone for heteroscedastic noise."""
    n = 6
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    mean_distance = HippsDimes.a2dmap_theory(
        truth, force_positive_definite=True
    )
    target = (3.0 * np.pi / 8.0) * np.square(mean_distance)
    variance = np.square(0.05 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)

    fitted, gram, connectivity, info = HippsDimes.fit_gaussian_noise_covariance(
        target,
        variance,
        max_iterations=30,
        relative_tolerance=1e-9,
        save_steps=[1],
    )

    assert info["converged"]
    assert info["status"] == "optimality_tolerance"
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-10)
    internal_gram = np.linalg.eigvalsh(gram)[1:]
    assert np.all(internal_gram > 0.0)
    connectivity_eigenvalues = np.linalg.eigvalsh(connectivity)
    assert connectivity_eigenvalues[-1] <= 1e-10
    assert np.count_nonzero(connectivity_eigenvalues < -1e-10) == n - 1
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    upper = np.triu_indices(n, k=1)
    # With Hamiltonian 1/2 sum k_ij r_ij^2, the exact first-order condition is
    # D_fit-D_obs = variance_ij * k_ij / 2.
    assert np.allclose(
        (fitted - target)[upper],
        0.5 * variance[upper] * connectivity[upper],
        rtol=1e-7,
        atol=1e-10,
    )
    assert sorted(info["connectivity_matrix_at_steps"]) == [1]


def test_gaussian_covariance_fit_supports_missing_heteroscedastic_pairs():
    """The structured COV operator should ignore unobserved pair entries."""
    n = 6
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    mean_distance = HippsDimes.a2dmap_theory(
        truth, force_positive_definite=True
    )
    target = (3.0 * np.pi / 8.0) * np.square(mean_distance)
    variance = np.square(0.05 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)
    missing_pairs = ((0, 3), (1, 4))
    for i, j in missing_pairs:
        target[i, j] = np.nan
        target[j, i] = np.nan

    fitted, _, connectivity, info = HippsDimes.fit_gaussian_noise_covariance(
        target,
        variance,
        max_iterations=40,
        relative_tolerance=1e-9,
    )

    observed = np.isfinite(target) & ~np.eye(n, dtype=bool)
    assert info["converged"]
    assert info["observed_pair_count"] == n * (n - 1) // 2 - len(missing_pairs)
    assert np.all(np.isfinite(fitted))
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    assert np.allclose(
        (fitted - target)[observed],
        (0.5 * variance * connectivity)[observed],
        rtol=1e-7,
        atol=1e-10,
    )


def test_covariance_fista_matches_covariance_solver_with_signed_couplings():
    """FISTA should monotonically converge to the calibrated COV solution."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    truth = numerics._connectivity_from_grounded_precision(grounded_precision)
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    variance = np.square(0.02 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)

    fitted_cov, gram_cov, connectivity_cov, info_cov = (
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            variance,
            max_iterations=40,
            relative_tolerance=1e-10,
        )
    )
    fitted_fista, gram_fista, connectivity_fista, info_fista = (
        HippsDimes.fit_gaussian_noise_covariance_fista(
            target,
            variance,
            max_iterations=300,
            relative_tolerance=1e-8,
            save_steps=[1],
        )
    )

    assert info_cov["converged"]
    assert info_fista["converged"]
    assert info_fista["status"] == "optimality_tolerance"
    assert info_fista["accelerated"]
    assert info_fista["monotone_restart"]
    assert info_fista["allows_signed_offdiagonal_connectivity"]
    assert info_fista["backtracking_reductions"] > 0
    assert sorted(info_fista["connectivity_matrix_at_steps"]) == [1]
    objective = info_fista["history"]["objective"]
    assert np.all(np.diff(objective) <= 1e-10)
    assert np.all(
        info_fista["history"]["minimum_internal_gram_eigenvalue"] > 0.0
    )
    assert np.any(truth[np.triu_indices(n, k=1)] < 0.0)
    assert np.any(connectivity_fista[np.triu_indices(n, k=1)] < 0.0)
    assert np.allclose(fitted_fista, fitted_cov, rtol=2e-8, atol=1e-9)
    assert np.allclose(gram_fista, gram_cov, rtol=2e-8, atol=1e-9)
    assert np.allclose(
        connectivity_fista, connectivity_cov, rtol=2e-8, atol=1e-9
    )
    assert np.max(np.abs(np.sum(connectivity_fista, axis=1))) <= 1e-12


def test_covariance_fista_handles_invalid_target_edm_without_projection():
    """The analytic proximal map should keep every iterate SPD for an invalid EDM."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    perturbed_value = 0.2 * target[0, n - 1]
    target[0, n - 1] = perturbed_value
    target[n - 1, 0] = perturbed_value
    variance = np.square(0.03 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)
    assert np.linalg.eigvalsh(
        numerics._centered_gram_from_squared_distances(target)
    )[0] < -1e-3

    _, gram, connectivity, info = (
        HippsDimes.fit_gaussian_noise_covariance_fista(
            target,
            variance,
            max_iterations=500,
            relative_tolerance=1e-7,
        )
    )

    assert info["converged"]
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-10)
    assert np.all(info["history"]["minimum_internal_gram_eigenvalue"] > 0.0)
    assert np.linalg.eigvalsh(gram)[1] > 0.0
    assert np.linalg.eigvalsh(-connectivity)[1] > 0.0


def test_gaussian_connectivity_cholesky_gradient_matches_finite_difference():
    """The signed connectivity dual should expose its exact packed gradient."""
    dimension = 4
    lower = np.tril_indices(dimension)
    diagonal_mask = lower[0] == lower[1]
    factor = np.eye(dimension)
    factor[lower] = np.linspace(-0.15, 0.20, len(lower[0]))
    factor[np.diag_indices(dimension)] = [0.8, 1.1, 0.9, 1.2]
    parameters = numerics._pack_log_cholesky(factor, lower, diagonal_mask)
    inverse_sqrt_reference = np.diag([0.7, 1.2, 0.9, 1.1])
    whitened_observed_covariance = np.array(
        [
            [1.0, 0.10, -0.03, 0.02],
            [0.10, 0.9, 0.04, 0.0],
            [-0.03, 0.04, 1.2, 0.08],
            [0.02, 0.0, 0.08, 0.8],
        ]
    )
    internal_variance = np.ones((dimension, dimension))
    np.fill_diagonal(internal_variance, 0.0)
    pinned_variance = np.linspace(0.5, 1.1, dimension)
    direction = np.linspace(-0.3, 0.4, len(parameters))
    direction /= np.linalg.norm(direction)

    objective, gradient, _ = (
        numerics._gaussian_connectivity_cholesky_objective_gradient(
            parameters,
            lower,
            diagonal_mask,
            inverse_sqrt_reference,
            whitened_observed_covariance,
            internal_variance,
            pinned_variance,
        )
    )
    step = 1e-6
    forward = numerics._gaussian_connectivity_cholesky_objective_gradient(
        parameters + step * direction,
        lower,
        diagonal_mask,
        inverse_sqrt_reference,
        whitened_observed_covariance,
        internal_variance,
        pinned_variance,
    )[0]
    backward = numerics._gaussian_connectivity_cholesky_objective_gradient(
        parameters - step * direction,
        lower,
        diagonal_mask,
        inverse_sqrt_reference,
        whitened_observed_covariance,
        internal_variance,
        pinned_variance,
    )[0]

    assert np.isfinite(objective)
    assert (forward - backward) / (2.0 * step) == pytest.approx(
        np.dot(gradient, direction), rel=1e-7, abs=1e-8
    )


def test_connectivity_cholesky_matches_covariance_solver_with_signed_couplings():
    """CHOL and COV should solve the same calibrated Gaussian model."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    truth = numerics._connectivity_from_grounded_precision(grounded_precision)
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    variance = np.square(0.02 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)

    fitted_cov, gram_cov, connectivity_cov, info_cov = (
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            variance,
            max_iterations=40,
            relative_tolerance=1e-10,
        )
    )
    fitted_chol, gram_chol, connectivity_chol, info_chol = (
        HippsDimes.fit_gaussian_noise_connectivity_cholesky(
            target,
            variance,
            max_iterations=300,
            gradient_tolerance=1e-10,
            function_tolerance=1e-15,
            stationarity_tolerance=1e-8,
            save_steps=[1],
        )
    )

    assert info_cov["converged"]
    assert info_chol["converged"]
    assert info_chol["allows_signed_offdiagonal_connectivity"]
    assert info_chol["maximum_function_evaluations"] == 300 * 51 + 1
    assert info_chol["optimizer_function_evaluations"] <= info_chol[
        "maximum_function_evaluations"
    ]
    assert "0.125" in info_chol["objective_definition"]
    assert sorted(info_chol["connectivity_matrix_at_steps"]) == [1]
    assert np.any(truth[np.triu_indices(n, k=1)] < 0.0)
    assert np.any(connectivity_chol[np.triu_indices(n, k=1)] < 0.0)
    assert np.allclose(fitted_chol, fitted_cov, rtol=2e-7, atol=2e-8)
    assert np.allclose(gram_chol, gram_cov, rtol=2e-7, atol=2e-8)
    assert np.allclose(connectivity_chol, connectivity_cov, rtol=2e-7, atol=2e-8)
    assert np.max(np.abs(np.sum(connectivity_chol, axis=1))) <= 1e-12
    eigenvalues = np.linalg.eigvalsh(connectivity_chol)
    assert eigenvalues[-1] <= 1e-10
    assert np.count_nonzero(eigenvalues < -1e-10) == n - 1
    upper = np.triu_indices(n, k=1)
    assert np.allclose(
        (fitted_chol - target)[upper],
        0.5 * variance[upper] * connectivity_chol[upper],
        rtol=2e-5,
        atol=2e-8,
    )


def test_gaussian_connectivity_coordinate_minimum_is_exact_and_stable():
    """A CIS edge step should attain its line optimum inside the SPD domain."""
    current_connectivity = -0.2
    current_squared_distance = 2.4
    observed_squared_distance = 2.0
    variance = 0.3

    delta, determinant_ratio, objective_change = (
        numerics._gaussian_connectivity_coordinate_minimum(
            current_connectivity,
            current_squared_distance,
            observed_squared_distance,
            variance,
        )
    )
    updated_connectivity = current_connectivity + delta
    updated_squared_distance = current_squared_distance / determinant_ratio

    assert determinant_ratio > 0.0
    assert determinant_ratio == pytest.approx(
        1.0 + delta * current_squared_distance / 3.0
    )
    assert objective_change < 0.0
    assert updated_squared_distance - observed_squared_distance == pytest.approx(
        0.5 * variance * updated_connectivity,
        rel=1e-13,
        abs=1e-13,
    )


def test_gaussian_l1_coordinate_minimum_selects_exact_zero_kink():
    """The exact L1 edge update should select zero when its KKT interval holds."""
    current_connectivity = 0.2
    current_squared_distance = 2.4
    rho = current_squared_distance / 3.0
    zero_determinant_ratio = 1.0 - rho * current_connectivity
    observed_squared_distance = (
        current_squared_distance / zero_determinant_ratio
    )

    delta, determinant_ratio, objective_change = (
        numerics._gaussian_connectivity_coordinate_minimum(
            current_connectivity,
            current_squared_distance,
            observed_squared_distance,
            0.3,
            0.1,
        )
    )

    assert current_connectivity + delta == 0.0
    assert determinant_ratio == pytest.approx(zero_determinant_ratio)
    assert objective_change < 0.0
    assert (
        current_squared_distance / determinant_ratio
        - observed_squared_distance
    ) == pytest.approx(0.0, abs=1e-14)


def test_l1_fista_and_coordinate_descent_match_with_exact_sparsity():
    """L1-FISTA and L1-CIS should attain the same unique signed sparse fit."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    variance = np.square(0.05 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)
    connectivity_l1 = 0.1

    fitted_fista, gram_fista, connectivity_fista, info_fista = (
        HippsDimes.fit_gaussian_noise_covariance_fista(
            target,
            variance,
            connectivity_l1=connectivity_l1,
            max_iterations=500,
            relative_tolerance=1e-9,
        )
    )
    fitted_cis, gram_cis, connectivity_cis, info_cis = (
        HippsDimes.fit_gaussian_noise_connectivity_coordinate_descent(
            target,
            variance,
            connectivity_l1=connectivity_l1,
            max_iterations=200,
            stationarity_tolerance=1e-9,
        )
    )

    assert info_fista["converged"]
    assert info_cis["converged"]
    assert info_fista["status"] == "stationarity_tolerance"
    assert info_fista["convergence_definition"] == (
        "physical L1 KKT stationarity residual norm"
    )
    assert info_fista["relative_gradient_norm"] > 1e-9
    assert info_fista["connectivity_l1"] == connectivity_l1
    assert info_cis["connectivity_l1"] == connectivity_l1
    assert info_fista["relative_stationarity_residual"] <= 1e-9
    assert info_cis["relative_stationarity_residual"] <= 1e-9
    assert np.count_nonzero(
        connectivity_cis[np.triu_indices(n, k=1)] == 0.0
    ) > 0
    assert np.all(np.diff(info_fista["history"]["objective"]) <= 1e-10)
    assert np.all(np.diff(info_cis["history"]["objective"]) <= 1e-10)
    assert np.linalg.norm(fitted_fista - fitted_cis) / np.linalg.norm(
        fitted_cis
    ) <= 1e-7
    assert np.linalg.norm(gram_fista - gram_cis) / np.linalg.norm(
        gram_cis
    ) <= 1e-7
    assert np.linalg.norm(connectivity_fista - connectivity_cis) / np.linalg.norm(
        connectivity_cis
    ) <= 2e-7
    assert info_fista["connectivity_dual_objective"] == pytest.approx(
        info_cis["connectivity_dual_objective"], rel=2e-9, abs=2e-8
    )


@pytest.mark.parametrize(
    "solver",
    [
        HippsDimes.fit_gaussian_noise_covariance_fista,
        HippsDimes.fit_gaussian_noise_connectivity_coordinate_descent,
    ],
)
@pytest.mark.parametrize("invalid_l1", [-1.0, np.inf, np.nan, True])
def test_exact_l1_solvers_reject_invalid_coefficient(solver, invalid_l1):
    """Both exact L1 solvers should reject invalid sparsity coefficients."""
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="connectivity_l1"):
        solver(target, 0.1, connectivity_l1=invalid_l1, max_iterations=1)


@pytest.mark.parametrize(
    "solver, solver_options",
    [
        (
            HippsDimes.fit_gaussian_noise_covariance_fista,
            {"max_iterations": 200, "relative_tolerance": 1e-8},
        ),
        (
            HippsDimes.fit_gaussian_noise_connectivity_coordinate_descent,
            {"max_iterations": 100, "stationarity_tolerance": 1e-8},
        ),
    ],
)
def test_exact_l1_solvers_preserve_zero_coefficient_path(solver, solver_options):
    """An explicit zero L1 coefficient should preserve the original result."""
    n = 4
    grounded_precision = np.eye(n - 1) + 0.1 * np.ones((n - 1, n - 1))
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )

    baseline = solver(target, 0.01, **solver_options)
    explicit_zero = solver(
        target, 0.01, connectivity_l1=0.0, **solver_options
    )

    for baseline_matrix, zero_matrix in zip(baseline[:3], explicit_zero[:3]):
        assert np.array_equal(baseline_matrix, zero_matrix)
    assert baseline[3]["iterations"] == explicit_zero[3]["iterations"]
    assert baseline[3]["objective"] == explicit_zero[3]["objective"]


def test_connectivity_coordinate_descent_matches_covariance_solver_signed_case():
    """CIS should converge monotonically to the calibrated COV solution."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    truth = numerics._connectivity_from_grounded_precision(grounded_precision)
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    variance = np.square(0.02 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)

    fitted_cov, gram_cov, connectivity_cov, info_cov = (
        HippsDimes.fit_gaussian_noise_covariance(
            target,
            variance,
            max_iterations=40,
            relative_tolerance=1e-10,
        )
    )
    fitted_cis, gram_cis, connectivity_cis, info_cis = (
        HippsDimes.fit_gaussian_noise_connectivity_coordinate_descent(
            target,
            variance,
            max_iterations=120,
            stationarity_tolerance=1e-10,
            save_steps=[1],
        )
    )

    assert info_cov["converged"]
    assert info_cis["converged"]
    assert info_cis["status"] == "stationarity_tolerance"
    assert info_cis["relative_stationarity_residual"] <= 1e-10
    assert info_cis["allows_signed_offdiagonal_connectivity"]
    assert info_cis["coordinate_updates"] == (
        info_cis["iterations"] * info_cis["pairs_per_sweep"]
    )
    assert sorted(info_cis["connectivity_matrix_at_steps"]) == [1]
    assert np.all(np.diff(info_cis["history"]["objective"]) <= 1e-11)
    assert np.all(info_cis["history"]["minimum_coordinate_determinant_ratio"] > 0.0)
    assert np.any(truth[np.triu_indices(n, k=1)] < 0.0)
    assert np.any(connectivity_cis[np.triu_indices(n, k=1)] < 0.0)
    assert np.allclose(fitted_cis, fitted_cov, rtol=3e-6, atol=1e-7)
    assert np.allclose(gram_cis, gram_cov, rtol=3e-6, atol=1e-7)
    assert np.allclose(connectivity_cis, connectivity_cov, rtol=3e-6, atol=3e-7)
    assert np.max(np.abs(np.sum(connectivity_cis, axis=1))) <= 1e-12
    eigenvalues = np.linalg.eigvalsh(connectivity_cis)
    assert eigenvalues[-1] <= 1e-10
    assert np.count_nonzero(eigenvalues < -1e-10) == n - 1
    upper = np.triu_indices(n, k=1)
    assert np.allclose(
        (fitted_cis - target)[upper],
        0.5 * variance[upper] * connectivity_cis[upper],
        rtol=1e-6,
        atol=2e-9,
    )


def test_coordinate_descent_safeguards_invalid_edm_boundary_start():
    """CIS should damp and retry instead of losing SPD to boundary roundoff."""
    n = 5
    grounded_precision = np.eye(n - 1) + 0.15 * np.ones((n - 1, n - 1))
    target = numerics._squared_distances_from_grounded_precision(
        grounded_precision
    )
    perturbed_value = 0.2 * target[0, n - 1]
    target[0, n - 1] = perturbed_value
    target[n - 1, 0] = perturbed_value
    variance = np.square(0.03 * np.maximum(target, 1e-8))
    np.fill_diagonal(variance, 0.0)
    assert np.linalg.eigvalsh(
        numerics._centered_gram_from_squared_distances(target)
    )[0] < -1e-3

    with pytest.warns(RuntimeWarning, match="max_iterations"):
        _, _, connectivity, info = (
            HippsDimes.fit_gaussian_noise_connectivity_coordinate_descent(
                target,
                variance,
                max_iterations=5,
                initial_gram_floor_relative=1e-8,
            )
        )

    assert info["sweep_retries"] > 0
    assert info["damped_sweeps"] > 0
    assert info["minimum_accepted_sweep_relaxation"] < 1.0
    assert np.all(info["history"]["sweep_relaxation"] > 0.0)
    assert np.all(info["history"]["sweep_relaxation"] <= 1.0)
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-8)
    assert np.linalg.eigvalsh(-connectivity[:-1, :-1])[0] > 0.0


def test_heteroscedastic_iterative_scaling_uses_elementwise_proximal_shrink():
    """A variance matrix should generalize the historical scalar shrink pairwise."""
    n = 5
    initial = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    mean_distance = HippsDimes.a2dmap_theory(
        initial, force_positive_definite=True
    )
    target = (3.0 * np.pi / 8.0) * np.square(mean_distance)
    variance = np.zeros((n, n), dtype=float)
    upper = np.triu_indices(n, k=1)
    variance[upper] = np.linspace(0.01, 0.10, len(upper[0]))
    variance[(upper[1], upper[0])] = variance[upper]
    learning_rate = 2.0
    optimizer = HippsDimes.Optimize(target, connectivity_matrix=initial.copy())

    _, _, _, observed, _ = optimizer.run_noisy(
        1,
        variance,
        learning_rate=learning_rate,
        method="IS",
        show_progress=False,
    )

    expected = initial / (1.0 + learning_rate * variance)
    expected = HippsDimes.a2a(expected)
    assert np.allclose(observed, expected, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    pytest.main([__file__]) 
