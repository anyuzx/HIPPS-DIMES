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


if __name__ == "__main__":
    pytest.main([__file__]) 