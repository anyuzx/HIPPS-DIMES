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

    def objective(candidate):
        _, log_determinant = np.linalg.slogdet(candidate)
        candidate_gram = basis @ candidate @ basis.T
        candidate_distances = numerics._squared_distances_from_gram(
            candidate_gram
        )
        residual = candidate_distances[pair_i, pair_j] - target_pairs
        return (
            -1.5 * log_determinant
            + 0.5 * np.sum(inverse_variance * np.square(residual))
        )

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
        HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            relative_noise_std=relative_std,
            max_iterations=500,
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
    assert info["phase_iterations"]["fista"] > 0
    assert np.all(np.diff(info["history"]["objective"]) <= 1e-10)
    assert np.count_nonzero(np.linalg.eigvalsh(gram) > 1e-12) == n - 1
    assert np.max(np.abs(np.sum(connectivity, axis=1))) <= 1e-12
    assert np.allclose(
        (fitted - target)[upper],
        0.5 * variance[upper] * connectivity[upper],
        rtol=3e-6,
        atol=2e-8,
    )
    assert info["relative_stationarity_residual"] <= 2e-8
    assert sorted(info["connectivity_matrix_at_steps"]) == [1]


def test_gaussian_covariance_absolute_noise_initializations_match():
    """Rouse and restart starts should reach the homoskedastic optimum."""
    target = _rouse_squared_distance_target()
    variance = 0.02
    solver_options = {
        "max_iterations": 1000,
        "relative_tolerance": 2e-8,
    }

    rouse = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        variance,
        **solver_options,
    )
    provided = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target,
        variance,
        initial_connectivity=HippsDimes.construct_connectivity_matrix_rouse(
            len(target), 0.4
        ),
        **solver_options,
    )

    for result in (rouse, provided):
        assert result[3]["converged"]
        scalar_calibration = result[3]["initialization"]["scalar_calibration"]
        assert scalar_calibration["scale_factor"] > 0.0
        assert scalar_calibration["objective_after"] <= scalar_calibration[
            "objective_before"
        ]
    assert rouse[3]["noise_model"] == "homoskedastic_absolute_variance"
    assert provided[3]["initialization"]["kind"] == "provided_connectivity"
    assert np.allclose(rouse[0], provided[0], rtol=1e-8, atol=1e-9)
    assert np.allclose(rouse[2], provided[2], rtol=1e-8, atol=1e-9)


def test_gaussian_covariance_relative_noise_supports_missing_pairs():
    target = _rouse_squared_distance_target()
    missing_pairs = ((0, 3), (1, 4))
    for i, j in missing_pairs:
        target[i, j] = np.nan
        target[j, i] = np.nan

    fitted, _, connectivity, info = (
        HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            relative_noise_std=0.05,
            max_iterations=1000,
            relative_tolerance=1e-8,
        )
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
        HippsDimes.fit_gaussian_noise_covariance_pdhg(
            target,
            noise_variance,
            relative_noise_std=relative_noise_std,
        )


def test_gaussian_covariance_gpu_request_fails_without_available_gpu(monkeypatch):
    target = _rouse_squared_distance_target(4)
    monkeypatch.setattr(numerics, "_CUPY_AVAILABLE", False)

    with pytest.raises(RuntimeError, match="CuPy.*accessible CUDA GPU"):
        HippsDimes.fit_gaussian_noise_covariance_hybrid(
            target,
            0.1,
            use_gpu=True,
        )


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
@pytest.mark.parametrize(
    "noise_options",
    [
        {"relative_noise_std": 0.05},
        {"noise_variance": 0.02},
    ],
)
def test_gaussian_covariance_gpu_matches_cpu_end_to_end(
    noise_options, working_cupy
):
    """GPU COV should match CPU for both noise models."""
    target = _rouse_squared_distance_target(6)
    solver_options = {
        "max_iterations": 1000,
        "relative_tolerance": 2e-8,
        "save_steps": [1],
        **noise_options,
    }

    cpu = HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target, **solver_options
    )
    progress_events = []
    gpu = HippsDimes.fit_gaussian_noise_covariance_hybrid(
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
    assert gpu[3]["phase_iterations"]["fista"] > 0
    assert sorted(gpu[3]["connectivity_matrix_at_steps"]) == [1]
    assert np.allclose(gpu[0], cpu[0], rtol=2e-6, atol=2e-8)
    assert np.allclose(gpu[1], cpu[1], rtol=2e-6, atol=2e-8)
    assert np.allclose(gpu[2], cpu[2], rtol=2e-5, atol=2e-7)
    assert gpu[3]["relative_stationarity_residual"] <= 2e-8
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
        "covariance_operator_norm",
        "covariance_optimization",
    }
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
