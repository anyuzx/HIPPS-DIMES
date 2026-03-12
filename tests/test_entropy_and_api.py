"""Additional tests for entropy and library API behavior."""

import os
import pickle

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import HippsDimes
from hipps_dimes.cli import main as cli_main


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


@pytest.mark.parametrize(
    ("input_type", "input_format", "expected_message"),
    [
        ("dmap", "text", "Distance map must be a square 2D array"),
        ("ddmap", "npy", "Squared distance map must be a square 2D array"),
        ("cmap", "text", "Contact map must be a square 2D array"),
    ],
)
def test_run_optimization_rejects_nonsquare_file_inputs(
    tmp_path, input_type, input_format, expected_message
):
    """File-backed inputs should fail fast with a clear square-matrix error."""
    matrix = np.arange(6.0).reshape(2, 3)
    path = tmp_path / ("matrix.npy" if input_format == "npy" else "matrix.txt")

    if input_format == "npy":
        np.save(path, matrix)
    else:
        np.savetxt(path, matrix)

    with pytest.raises(ValueError, match=expected_message):
        HippsDimes.run_optimization(
            input_path=str(path),
            input_type=input_type,
            input_format=input_format,
            method="IS",
            iteration=1,
            learning_rate=5.0,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


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


def test_run_optimization_save_pickle_writes_only_pickle_file(tmp_path):
    """save_pickle should suppress default artifact files and write only the pickle."""
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
        ensemble=2,
        save_steps=[2],
        save_pickle=True,
        verbose=False,
        show_progress=False,
    )

    pickle_path = tmp_path / "run_HIPPS_DIMES_results.pkl"
    assert pickle_path.exists()
    assert not (tmp_path / "run_iteration_series.csv").exists()
    assert not (tmp_path / "run_run_parameters.csv").exists()
    assert not (tmp_path / "run_dmap_final.txt").exists()
    assert not (tmp_path / "run_connectivity_matrix.txt").exists()
    assert not (tmp_path / "run_connectivity_matrix_iter2.txt").exists()
    assert not (tmp_path / "run.xyz").exists()

    with pickle_path.open("rb") as fin:
        saved_results = pickle.load(fin)

    assert saved_results.keys() == results.keys()
    assert "xyzs" in saved_results
    assert saved_results["xyzs"].shape[0] == 2
    assert sorted(saved_results["connectivity_matrix_at_steps"]) == [2]
    pd.testing.assert_frame_equal(saved_results["iteration_series"], results["iteration_series"])
    pd.testing.assert_frame_equal(saved_results["run_parameters"], results["run_parameters"])
    assert np.allclose(saved_results["dmap_final"], results["dmap_final"])
    assert np.allclose(saved_results["connectivity_matrix"], results["connectivity_matrix"])


def test_run_optimization_save_pickle_requires_output_prefix():
    """save_pickle needs an output prefix to derive the pickle filename."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    with pytest.raises(ValueError, match="output_prefix must be provided when save_pickle=True"):
        HippsDimes.run_optimization(
            input_matrix=dmap_target,
            input_type="dmap",
            method="IS",
            iteration=1,
            learning_rate=5.0,
            save_pickle=True,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


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


def test_run_optimization_progress_callback_receives_structured_updates():
    """Progress callbacks should receive monotonic structured updates."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    updates = []

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
        progress_callback=updates.append,
    )

    assert [update["iteration"] for update in updates] == [1, 2, 3, 4, 5]
    assert all(update["total"] == 5 for update in updates)
    assert all(update["stage"] == "optimization" for update in updates)
    assert all(update["method"] == "IS" for update in updates)
    assert all(update["general_method"] == "optimization" for update in updates)
    assert np.isclose(updates[-1]["entropy"], results["iteration_series"]["entropy"].iloc[-1])


def test_run_optimization_noisy_save_steps_return_snapshots_without_output_prefix():
    """Noisy optimization should return save-step snapshots even in library mode."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="IS",
        iteration=3,
        learning_rate=5.0,
        gaussian_noise_variance=0.1,
        save_steps=[1, 3],
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    assert "connectivity_matrix_at_steps" in results
    assert sorted(results["connectivity_matrix_at_steps"]) == [1, 3]
    assert results["connectivity_matrix_at_steps"][1].shape == (n, n)
    assert results["connectivity_matrix_at_steps"][3].shape == (n, n)


def test_run_optimization_applies_eigh_threads_in_library(monkeypatch):
    """Library callers should get the same thread-setting behavior as the CLI."""
    n = 4
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    thread_env_vars = [
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ]

    for name in thread_env_vars:
        monkeypatch.delenv(name, raising=False)

    HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="IS",
        iteration=1,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
        eigh_threads=1,
    )

    assert {name: os.environ.get(name) for name in thread_env_vars} == {
        name: "1" for name in thread_env_vars
    }


def test_cli_save_pickle_writes_pickle_only(tmp_path):
    """CLI --save-pickle should route through to pickle-only output behavior."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)
    dmap_path = tmp_path / "dmap.npy"
    output_prefix = tmp_path / "cli_run"
    np.save(dmap_path, dmap_target)

    runner = CliRunner()
    result = runner.invoke(
        cli_main,
        [
            str(dmap_path),
            str(output_prefix),
            "--input-type",
            "dmap",
            "--input-format",
            "npy",
            "--iteration",
            "2",
            "--learning-rate",
            "5.0",
            "--no-xyzs",
            "--save-pickle",
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "cli_run_HIPPS_DIMES_results.pkl").exists()
    assert not (tmp_path / "cli_run_dmap_final.txt").exists()
    assert not (tmp_path / "cli_run_connectivity_matrix.txt").exists()
    assert not (tmp_path / "cli_run_iteration_series.csv").exists()



def test_run_optimization_show_progress_false_suppresses_progress_bar_output(capsys):
    """show_progress=False should prevent tqdm progress output."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="IS",
        iteration=3,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    captured = capsys.readouterr()
    combined_output = captured.out + captured.err
    assert "Performing optimization" not in combined_output
    assert "iteration/s" not in combined_output
