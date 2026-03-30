"""Basic tests for HippsDimes package."""

import numpy as np
import pytest

# Import the main module
import HippsDimes


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


def test_compute_modulus():
    """Test the modulus computation."""
    n = 4
    k = 1.0
    A = HippsDimes.construct_connectivity_matrix_rouse(n, k)
    freq = np.logspace(-2, 2, 10)
    zeta = 1.0
    
    storage_mod, loss_mod = HippsDimes.compute_modulus(A, freq, zeta)
    
    # Check shapes
    assert storage_mod.shape == (len(freq), 2)
    assert loss_mod.shape == (len(freq), 2)
    
    # Check that frequencies match
    assert np.allclose(storage_mod[:, 0], freq)
    assert np.allclose(loss_mod[:, 0], freq)


def test_compute_stress_relaxation():
    """Stress relaxation should normalize the non-zero mode sum by chain length."""
    A = HippsDimes.construct_connectivity_matrix_rouse(4, 1.0)
    t = [0.0, 0.1, 1.0, 10.0]

    relaxation = HippsDimes.compute_stress_relaxation(A, t, zeta=1.0)

    eigvals = np.linalg.eigvalsh(A)
    lam_nz = eigvals[np.abs(eigvals) > 1e-12]
    expected = np.sum(
        np.exp(-np.asarray(t)[:, None] / (-1.0 / lam_nz)[None, :]),
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
