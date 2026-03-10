"""Additional tests for entropy and library API behavior."""

import numpy as np
import pandas as pd

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


def test_run_optimization_smoke_dmap_npy_input(tmp_path):
    """run_optimization should accept dmap input from a .npy file."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    dmap_path = tmp_path / "dmap.npy"
    np.save(dmap_path, dmap_target)

    results = HippsDimes.run_optimization(
        input_path=str(dmap_path),
        input_type="dmap",
        input_format="npy",
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
    )

    assert "dmap_final" in results
    assert "connectivity_matrix" in results
    assert results["connectivity_matrix"].shape == (n, n)


def test_run_optimization_smoke_cmap_npy_input(tmp_path):
    """run_optimization should accept cmap input from a .npy file."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    cmap_target = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    cmap_path = tmp_path / "cmap.npy"
    np.save(cmap_path, cmap_target)

    results = HippsDimes.run_optimization(
        input_path=str(cmap_path),
        input_type="cmap",
        input_format="npy",
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
    )

    assert "dmap_final" in results
    assert "connectivity_matrix" in results
    assert "cmap_final" in results
    assert "rc_optimal" in results


def test_run_optimization_writes_default_log_files(tmp_path):
    """When output_prefix is provided, both log files should be written by default."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    output_prefix = tmp_path / "run"

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        output_prefix=str(output_prefix),
        input_type="dmap",
        method="IS",
        iteration=3,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
    )

    iteration_series_path = tmp_path / "run_iteration_series.csv"
    run_parameters_path = tmp_path / "run_run_parameters.csv"

    assert "iteration_series" in results
    assert "run_parameters" in results
    assert iteration_series_path.exists()
    assert run_parameters_path.exists()

    iteration_series_df = pd.read_csv(iteration_series_path)
    assert list(iteration_series_df.columns) == ["iteration", "loss", "entropy"]
    assert len(iteration_series_df) == 3

    run_parameters_df = pd.read_csv(run_parameters_path)
    assert list(run_parameters_df.columns) == ["parameter", "value"]
    assert "method" in set(run_parameters_df["parameter"])


def test_run_optimization_no_log_disables_both_log_files(tmp_path):
    """no_log should suppress both log-file writes while leaving main outputs intact."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    output_prefix = tmp_path / "run"

    HippsDimes.run_optimization(
        input_matrix=dmap_target,
        output_prefix=str(output_prefix),
        input_type="dmap",
        method="IS",
        iteration=3,
        learning_rate=5.0,
        no_log=True,
        no_xyzs=True,
        verbose=False,
    )

    assert not (tmp_path / "run_iteration_series.csv").exists()
    assert not (tmp_path / "run_run_parameters.csv").exists()
    assert (tmp_path / "run_dmap_final.txt").exists()
    assert (tmp_path / "run_connectivity_matrix.txt").exists()


def test_run_optimization_saves_target_cmap_for_cooler_or_hic_formats(tmp_path):
    """The internal target cmap should be saved for cooler/hic cmap inputs."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    cmap_target = HippsDimes.a2cmap_theory(A_true, rc=1.0)

    for input_format in ("cooler", "hic"):
        output_prefix = tmp_path / f"run_{input_format}"
        HippsDimes.run_optimization(
            input_matrix=cmap_target,
            output_prefix=str(output_prefix),
            input_type="cmap",
            input_format=input_format,
            method="IS",
            iteration=3,
            learning_rate=5.0,
            no_xyzs=True,
            verbose=False,
        )

        saved_target_cmap = np.loadtxt(tmp_path / f"run_{input_format}_cmap_target.txt")
        assert np.allclose(saved_target_cmap, cmap_target)


def test_run_optimization_does_not_save_target_cmap_for_text_input(tmp_path):
    """Text cmap inputs should not write the internal target cmap output file."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    cmap_target = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    output_prefix = tmp_path / "run_text"

    HippsDimes.run_optimization(
        input_matrix=cmap_target,
        output_prefix=str(output_prefix),
        input_type="cmap",
        input_format="text",
        method="IS",
        iteration=3,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
    )

    assert not (tmp_path / "run_text_cmap_target.txt").exists()
