"""Additional tests for entropy and library API behavior."""

import numpy as np

import HippsDimes


def test_compute_entropy_from_random_connectivity_matrix():
    """Entropy should match the direct eigenvalue-based calculation."""
    rng = np.random.default_rng(12345)
    n = 8

    # Build a connected Laplacian-like connectivity matrix:
    # off-diagonal entries are nonnegative couplings; diagonal enforces row-sum zero.
    offdiag = rng.uniform(0.05, 1.0, size=(n, n))
    offdiag = 0.5 * (offdiag + offdiag.T)
    np.fill_diagonal(offdiag, 0.0)

    A = offdiag.copy()
    np.fill_diagonal(A, -np.sum(offdiag, axis=1))

    entropy = HippsDimes.compute_entropy_from_A(A)
    eigvals = np.linalg.eigvalsh(-A)
    expected = np.sum(-np.log(eigvals[eigvals > 1e-12]))

    assert np.isfinite(entropy)
    assert np.isclose(entropy, expected, rtol=1e-10, atol=1e-10)


def test_run_optimization_smoke_dmap():
    """Basic run_optimization call should execute and return expected outputs."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
    )

    assert "log" in results
    assert "dmap_final" in results
    assert "connectivity_matrix" in results
    assert "xyzs" not in results

    log_df = results["log"]
    assert list(log_df.columns) == ["iteration", "loss", "entropy"]
    assert len(log_df) == 5

    A_est = results["connectivity_matrix"]
    assert A_est.shape == (n, n)
    assert np.allclose(A_est, A_est.T)
    assert np.allclose(np.sum(A_est, axis=1), 0.0, atol=1e-8)
