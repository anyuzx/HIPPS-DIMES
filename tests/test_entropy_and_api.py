"""Additional tests for entropy and library API behavior."""

import json
import os
import pickle

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import HippsDimes
import hipps_dimes.api as api
import hipps_dimes.covariance_pdhg as covariance_pdhg
import hipps_dimes.numerics as numerics
from hipps_dimes.api import _repair_fully_missing_loci_nearest_neighbors, _summarize_missing_data
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


def test_contact_threshold_fit_uses_normalized_observed_pairs_only():
    """Threshold fitting should ignore zero/NaN pairs and recover known rc."""
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(8, 1.0)
    model_dmap = HippsDimes.a2dmap_theory(connectivity)
    expected_rc = 0.7
    experimental_cmap = 4.0 * HippsDimes.dmap2cmap(
        model_dmap, expected_rc
    )
    experimental_cmap[1, 5] = experimental_cmap[5, 1] = 0.0
    experimental_cmap[2, 6] = experimental_cmap[6, 2] = np.nan

    fitted_rc = numerics._optimize_contact_threshold(
        model_dmap, experimental_cmap
    )

    assert fitted_rc == pytest.approx(expected_rc, rel=1e-6)


@pytest.mark.parametrize("method", ["IS", "GD", "DI", "COV"])
@pytest.mark.filterwarnings(
    "ignore:fit_gaussian_noise_covariance.*stopped:RuntimeWarning"
)
def test_contact_threshold_postprocessing_reuses_any_solver_dmap(
    monkeypatch, method
):
    """All contact-map solvers should use their returned distance map for rc."""
    connectivity = HippsDimes.construct_connectivity_matrix_rouse(6, 1.0)
    cmap_target = HippsDimes.a2cmap_theory(connectivity, rc=1.0)
    captured = {}

    def capture_threshold_fit(model_dmap, experimental_cmap):
        captured["model_dmap"] = np.asarray(model_dmap).copy()
        captured["experimental_cmap"] = np.asarray(experimental_cmap).copy()
        return 0.8

    monkeypatch.setattr(
        api, "_optimize_contact_threshold", capture_threshold_fit
    )
    results = HippsDimes.run_optimization(
        input_matrix=cmap_target,
        input_type="cmap",
        input_format="text",
        method=method,
        gaussian_noise_relative_std=0.1 if method == "COV" else None,
        iteration=2,
        learning_rate=1e-4 if method == "GD" else 5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    assert np.allclose(captured["model_dmap"], results["dmap_final"])
    assert np.allclose(captured["experimental_cmap"], cmap_target)
    assert np.allclose(
        results["cmap_final"],
        HippsDimes.dmap2cmap(results["dmap_final"], 0.8),
    )


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


@pytest.mark.parametrize(
    "optimizer, expected_descriptions, expected_stages",
    [
        (
            "hybrid",
            [
                "COV operator norm",
                "COV PDHG optimization",
                "COV FISTA optimization",
            ],
            {
                "covariance_operator_norm",
                "covariance_optimization",
            },
        ),
        (
            "pdhg",
            ["COV operator norm", "COV PDHG optimization"],
            {"covariance_operator_norm", "covariance_optimization"},
        ),
    ],
)
def test_cov_progress_distinguishes_pdhg_and_fista(
    monkeypatch,
    optimizer,
    expected_descriptions,
    expected_stages,
):
    """COV should render solver-specific stages and preserve callbacks."""

    class FakeProgressBar:
        instances = []

        def __init__(self, *, total, desc, unit):
            self.total = total
            self.desc = desc
            self.unit = unit
            self.n = 0
            self.postfix = {}
            self.closed = False
            self.instances.append(self)

        def set_postfix(self, **values):
            values.pop("refresh", None)
            self.postfix = values

        def update(self, amount):
            self.n += amount

        def close(self):
            self.closed = True

    monkeypatch.setattr(api, "tqdm", FakeProgressBar)

    n = 6
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(truth)
    updates = []

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="COV",
        covariance_optimizer=optimizer,
        gaussian_noise_variance=0.1,
        iteration=200,
        no_xyzs=True,
        verbose=False,
        show_progress=True,
        progress_callback=updates.append,
    )

    bars = FakeProgressBar.instances
    assert [bar.desc for bar in bars] == expected_descriptions
    assert all(bar.closed and bar.n > 0 for bar in bars)
    assert {update["stage"] for update in updates} == expected_stages
    optimization_bar = bars[-1]
    assert "cg_iterations" not in optimization_bar.postfix
    assert "relative_kkt" in optimization_bar.postfix
    assert results["covariance_optimization"]["converged"]


def test_run_optimization_cov_save_steps_return_snapshots_without_output_prefix():
    """COV should return save-step snapshots and calibrated diagnostics."""
    n = 6
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="COV",
        iteration=200,
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
    assert results["covariance_optimization"]["converged"]
    assert results["covariance_optimization"]["algorithm"] == "hybrid"
    assert results["covariance_optimization"]["pdhg"]["variance_whitened"]
    assert results["covariance_optimization"]["pdhg"][
        "inverse_free_runtime_kkt"
    ]
    assert results["covariance_optimization"]["handoff"]["reached"]
    assert results["covariance_optimization"]["phase_iterations"]["pdhg"] > 0
    assert results["covariance_optimization"]["phase_iterations"]["fista"] > 0
    assert (
        results["covariance_optimization"][
            "relative_eliminated_kkt_residual"
        ]
        <= 1e-5
    )
    assert results["covariance_optimization"]["initialization"]["kind"] == "rouse"
    assert results["gram_matrix"].shape == (n, n)
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )
    assert parameters["covariance_optimizer_requested"] == "hybrid"
    assert parameters["covariance_optimizer_resolved"] == "hybrid"
    assert parameters["covariance_relative_tolerance"] == pytest.approx(1e-5)
    assert parameters["covariance_absolute_tolerance"] == pytest.approx(1e-10)
    assert parameters[
        "covariance_handoff_relative_tolerance"
    ] == pytest.approx(1e-2)
    assert json.loads(parameters["covariance_phase_iterations"])["fista"] > 0
    assert bool(parameters["covariance_independent_kkt_converged"])
    assert parameters["covariance_fista_initial_step_size"] > 0.0
    assert parameters["covariance_fista_final_step_size"] > 0.0


def test_run_optimization_cov_relative_noise_is_applied_after_dmap_conversion():
    n = 5
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(truth)
    relative_std = 0.08

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="COV",
        gaussian_noise_relative_std=relative_std,
        iteration=200,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    ddmap_target = (3.0 * np.pi / 8.0) * np.square(dmap_target)
    expected_variance = np.square(
        relative_std * ddmap_target[np.triu_indices(n, k=1)]
    )
    info = results["covariance_optimization"]
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert info["converged"]
    assert info["noise_model"] == "heteroskedastic_relative_std"
    assert info["noise_variance_minimum"] == pytest.approx(
        np.min(expected_variance)
    )
    assert info["noise_variance_median"] == pytest.approx(
        np.median(expected_variance)
    )
    assert info["noise_variance_maximum"] == pytest.approx(
        np.max(expected_variance)
    )
    assert parameters["gaussian_noise_relative_std"] == relative_std
    assert parameters["covariance_initialization_resolved"] == "rouse"


def test_run_optimization_cov_cmap_rouse_ignores_missing_pairs(tmp_path):
    """Direct contact input should preserve missing pairs through COV setup."""
    n = 6
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    cmap_target = HippsDimes.a2cmap_theory(truth, rc=1.0)
    missing_pairs = ((0, 4), (1, 5))
    for i, j in missing_pairs:
        cmap_target[i, j] = 0.0
        cmap_target[j, i] = 0.0
    cmap_path = tmp_path / "cmap_missing.npy"
    np.save(cmap_path, cmap_target)

    results = HippsDimes.run_optimization(
        input_path=str(cmap_path),
        input_type="cmap",
        input_format="npy",
        method="COV",
        gaussian_noise_relative_std=0.1,
        iteration=300,
        ignore_missing_data=True,
        not_normalize=True,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    info = results["covariance_optimization"]
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert info["converged"]
    assert info["observed_pair_count"] == n * (n - 1) // 2 - len(missing_pairs)
    assert info["initialization"]["kind"] == "rouse"
    assert parameters["missing_pairs"] == len(missing_pairs)
    assert bool(parameters["ignore_missing_data"]) is True


@pytest.mark.skipif(
    not numerics.is_gpu_available(),
    reason="requires CuPy and an accessible CUDA GPU",
)
def test_run_optimization_cov_uses_gpu_backend():
    n = 5
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(truth)

    results = HippsDimes.run_optimization(
        input_matrix=dmap_target,
        input_type="dmap",
        input_format="text",
        method="COV",
        covariance_optimizer="hybrid",
        gaussian_noise_variance=0.05,
        iteration=500,
        use_gpu=True,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert results["covariance_optimization"]["converged"]
    assert results["covariance_optimization"]["backend"] == "gpu"
    assert parameters["use_gpu_requested"] is True
    assert parameters["use_gpu_enabled"] is True
    assert parameters["covariance_initialization_backend"] == "cpu"
    assert parameters["covariance_initialization_wall_seconds"] >= 0.0
    assert parameters["covariance_backend"] == "gpu"
    assert parameters["covariance_dtype"] == "float64"
    assert parameters["covariance_gpu_device"] == numerics.get_gpu_name()
    assert parameters["covariance_cupy_version"] == numerics.cp.__version__
    assert parameters["covariance_optimizer_resolved"] == "hybrid"
    assert parameters["covariance_fista_initial_step_size"] > 0.0
    assert parameters["covariance_fista_final_step_size"] > 0.0


@pytest.mark.parametrize("method", ["IS", "GD", "DI"])
def test_run_optimization_rejects_gaussian_noise_outside_cov(method):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="only with method='COV'"):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="ddmap",
            method=method,
            gaussian_noise_variance=0.1,
            no_xyzs=True,
            verbose=False,
        )


@pytest.mark.parametrize(
    "options, message",
    [
        ({}, "requires exactly one"),
        (
            {
                "gaussian_noise_variance": 0.1,
                "gaussian_noise_relative_std": 0.1,
            },
            "requires exactly one",
        ),
        ({"gaussian_noise_variance": 0.1, "lamd": 0.2}, "lamd"),
        (
            {
                "gaussian_noise_variance": 0.1,
                "enforce_nonnegative_connectivity_matrix": True,
            },
            "nonnegative-connectivity",
        ),
        (
            {"gaussian_noise_variance": 0.1, "gpu_float32": True},
            "float64 only",
        ),
        (
            {"gaussian_noise_variance": 0.1, "covariance_optimizer": "bad"},
            "covariance_optimizer",
        ),
        (
            {
                "gaussian_noise_variance": 0.1,
                "covariance_optimizer": "newton",
            },
            "covariance_optimizer",
        ),
        (
            {
                "gaussian_noise_variance": 0.1,
                "covariance_optimizer": "fista",
            },
            "covariance_optimizer",
        ),
        (
            {
                "gaussian_noise_variance": 0.1,
                "covariance_relative_tolerance": 0.0,
                "covariance_absolute_tolerance": 0.0,
            },
            "at least one covariance convergence tolerance",
        ),
        (
            {
                "gaussian_noise_variance": 0.1,
                "covariance_handoff_relative_tolerance": 0.0,
            },
            "covariance_handoff_relative_tolerance",
        ),
        (
            {
                "gaussian_noise_variance": 0.1,
                "covariance_handoff_relative_tolerance": 1e-6,
            },
            "must not be smaller",
        ),
    ],
)
def test_run_optimization_rejects_invalid_cov_combinations(options, message):
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match=message):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="ddmap",
            method="COV",
            no_xyzs=True,
            verbose=False,
            **options,
        )


def test_run_optimization_rejects_covariance_solver_options_outside_cov():
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(ValueError, match="only with method='COV'"):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="ddmap",
            method="IS",
            covariance_optimizer="pdhg",
            no_xyzs=True,
            verbose=False,
        )


def test_legacy_optimize_no_longer_exposes_noisy_solver():
    target = np.array([[0.0, 1.0], [1.0, 0.0]])

    model = HippsDimes.Optimize(target)

    assert not hasattr(model, "run_noisy")


def test_cli_help_exposes_cov_solver_and_noise_contract():
    result = CliRunner().invoke(cli_main, ["--help"])
    parameter_names = {parameter.name for parameter in cli_main.params}

    assert result.exit_code == 0
    assert not hasattr(HippsDimes, "nearest_edm")
    assert not hasattr(numerics, "nearest_edm")
    assert "COV" in result.output
    assert "gaussian_noise_relative_std" in parameter_names
    assert "covariance_initialization" not in parameter_names
    assert "covariance_optimizer" in parameter_names
    assert "covariance_relative_tolerance" in parameter_names
    assert "covariance_absolute_tolerance" in parameter_names
    assert "covariance_handoff_relative_tolerance" in parameter_names
    assert "repair_fully_missing_loci" in parameter_names
    relative_noise_parameter = next(
        parameter
        for parameter in cli_main.params
        if parameter.name == "gaussian_noise_relative_std"
    )
    assert relative_noise_parameter.type.min_open
    optimizer_parameter = next(
        parameter
        for parameter in cli_main.params
        if parameter.name == "covariance_optimizer"
    )
    assert optimizer_parameter.default == "hybrid"
    assert set(optimizer_parameter.type.choices) == {"hybrid", "pdhg"}
    handoff_parameter = next(
        parameter
        for parameter in cli_main.params
        if parameter.name == "covariance_handoff_relative_tolerance"
    )
    assert handoff_parameter.default == pytest.approx(1e-2)


def test_cli_routes_cov_to_pdhg_with_configurable_tolerance(tmp_path):
    n = 5
    truth = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(truth)
    input_path = tmp_path / "target.npy"
    output_prefix = tmp_path / "cov_pdhg"
    np.save(input_path, dmap_target)

    result = CliRunner().invoke(
        cli_main,
        [
            str(input_path),
            str(output_prefix),
            "--input-type",
            "dmap",
            "--input-format",
            "npy",
            "--method",
            "COV",
            "--gaussian-noise-variance",
            "0.1",
            "--covariance-optimizer",
            "pdhg",
            "--covariance-relative-tolerance",
            "1e-5",
            "--iteration",
            "200",
            "--no-xyzs",
            "--save-pickle",
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with open(
        tmp_path / "cov_pdhg_HIPPS_DIMES_results.pkl", "rb"
    ) as fin:
        payload = pickle.load(fin)
    info = payload["covariance_optimization"]
    assert info["algorithm"] == "pdhg"
    assert info["pdhg"]["variance_whitened"]
    assert info["pdhg"]["inverse_free_runtime_kkt"]
    assert info["converged"]
    assert info["relative_eliminated_kkt_residual"] <= 1e-5


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


def test_ignore_missing_data_l2_stays_bounded_over_long_iterations():
    """L2 with missing-data constraints should not drift into late-iteration blow-up."""
    A_true = HippsDimes.construct_connectivity_matrix_rouse(10, 1.0)
    cmap = HippsDimes.a2cmap_theory(A_true, rc=1.0)

    missing_mask = np.ones_like(cmap, dtype=bool)
    for i in range(cmap.shape[0]):
        for j in range(cmap.shape[1]):
            if abs(i - j) >= 2 and (i + 2 * j) % 3 != 0:
                missing_mask[i, j] = False
    missing_mask = np.logical_or(missing_mask, missing_mask.T)
    np.fill_diagonal(missing_mask, True)

    cmap_missing = cmap.copy()
    cmap_missing[~missing_mask] = 0.0
    cmap_missing = 0.5 * (cmap_missing + cmap_missing.T)
    np.fill_diagonal(cmap_missing, 1.0)

    results = HippsDimes.run_optimization(
        input_matrix=cmap_missing,
        input_type="cmap",
        input_format="text",
        ignore_missing_data=True,
        method="IS",
        iteration=1000,
        learning_rate=10.0,
        lamd=0.01,
        reg="L2",
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    loss = results["iteration_series"]["loss"].to_numpy()
    assert np.isfinite(loss).all()
    assert np.max(loss) < 0.7


def test_repair_fully_missing_contact_locus_repairs_only_nearest_neighbors():
    """Nearest-neighbor repair should only fill i-1/i+1 pairs for fully missing loci."""
    A_true = HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    cmap = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    cmap_missing = cmap.copy()
    cmap_missing[2, :] = 0.0
    cmap_missing[:, 2] = 0.0
    np.fill_diagonal(cmap_missing, 1.0)

    summary = _summarize_missing_data(cmap_missing, "cmap")
    assert summary["fully_missing_loci"] == [2]

    repaired, repair_info = _repair_fully_missing_loci_nearest_neighbors(
        cmap_missing,
        "cmap",
        summary["missing_mask"],
        summary["fully_missing_loci"],
    )

    repaired_summary = _summarize_missing_data(repaired, "cmap")
    assert repair_info["nearest_neighbor_repaired_pairs"] == [(1, 2), (2, 3)]
    assert repaired[2, 1] > 0.0
    assert repaired[2, 3] > 0.0
    assert repaired[2, 0] == 0.0
    assert repaired[2, 4] == 0.0
    assert repaired_summary["fully_missing_loci"] == []


@pytest.mark.parametrize("matrix_kind", ["dmap", "ddmap"])
def test_repair_fully_missing_distance_locus_excludes_diagonal_source(
    matrix_kind,
):
    """Distance repair must not copy the zero diagonal off diagonal."""
    dmap = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    )
    target = dmap if matrix_kind == "dmap" else np.square(dmap)
    target[2, :] = np.nan
    target[:, 2] = np.nan
    target[2, 2] = 0.0
    summary = _summarize_missing_data(target, matrix_kind)

    repaired, repair_info = _repair_fully_missing_loci_nearest_neighbors(
        target,
        matrix_kind,
        summary["missing_mask"],
        summary["fully_missing_loci"],
    )

    assert repair_info["nearest_neighbor_repaired_pairs"] == [(1, 2), (2, 3)]
    assert repaired[1, 2] > 0.0
    assert repaired[2, 3] > 0.0
    assert np.isnan(repaired[0, 2])
    assert np.isnan(repaired[2, 4])


def test_run_optimization_records_missing_data_analysis_and_repair():
    """run_optimization should record missing-data analysis and nearest-neighbor repairs."""
    A_true = HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    cmap = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    cmap_missing = cmap.copy()
    cmap_missing[2, :] = 0.0
    cmap_missing[:, 2] = 0.0
    np.fill_diagonal(cmap_missing, 1.0)

    results = HippsDimes.run_optimization(
        input_matrix=cmap_missing,
        input_type="cmap",
        input_format="text",
        ignore_missing_data=True,
        repair_fully_missing_loci=True,
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    run_parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert run_parameters["missing_pairs"] == 4
    assert np.isclose(float(run_parameters["missing_pair_fraction"]), 0.4)
    assert json.loads(run_parameters["fully_missing_loci"]) == [2]
    assert bool(run_parameters["repair_fully_missing_loci"]) is True
    assert run_parameters["nearest_neighbor_repair_count"] == 2
    assert json.loads(run_parameters["nearest_neighbor_repaired_pairs"]) == [[1, 2], [2, 3]]
    assert json.loads(run_parameters["remaining_fully_missing_loci"]) == []


def test_cov_repair_is_explicit_and_enters_final_noise_model():
    target = HippsDimes.a2cmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(5, 1.0),
        rc=1.0,
    )
    target[2, :] = 0.0
    target[:, 2] = 0.0
    target[2, 2] = 1.0

    results = HippsDimes.run_optimization(
        input_matrix=target,
        input_type="cmap",
        input_format="text",
        method="COV",
        gaussian_noise_relative_std=0.1,
        ignore_missing_data=True,
        repair_fully_missing_loci=True,
        iteration=500,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )
    info = results["covariance_optimization"]
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert info["converged"]
    assert info["observed_pair_count"] == 8
    assert parameters["nearest_neighbor_repair_count"] == 2
    assert bool(parameters["repair_fully_missing_loci"]) is True


def test_cov_returned_iterate_is_consistent_across_reports(
    monkeypatch,
    capsys,
):
    target = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(4, 1.0)
    )
    original = covariance_pdhg._inverse_free_residuals
    call_count = 0

    def prefer_first_runtime_iterate(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        residuals = dict(original(*args, **kwargs))
        residuals["eliminated_relative"] = float(call_count - 1)
        return residuals

    monkeypatch.setattr(
        covariance_pdhg,
        "_inverse_free_residuals",
        prefer_first_runtime_iterate,
    )
    with pytest.warns(RuntimeWarning, match="without satisfying"):
        results = HippsDimes.run_optimization(
            input_matrix=target,
            input_type="dmap",
            input_format="text",
            method="COV",
            covariance_optimizer="pdhg",
            gaussian_noise_variance=0.1,
            covariance_relative_tolerance=0.0,
            covariance_absolute_tolerance=1e-30,
            iteration=3,
            no_xyzs=True,
            verbose=True,
            show_progress=False,
        )

    info = results["covariance_optimization"]
    series = results["iteration_series"]
    parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )
    output = capsys.readouterr().out

    assert info["iterations"] == 3
    assert info["returned_iteration"] == 1
    assert series["is_returned_iterate"].sum() == 1
    assert bool(series.loc[0, "is_returned_iterate"])
    assert info["loss"] == pytest.approx(series.loc[0, "loss"])
    assert info["entropy"] == pytest.approx(series.loc[0, "entropy"])
    assert info["objective"] == pytest.approx(series.loc[0, "objective"])
    assert parameters["covariance_iterations_executed"] == 3
    assert parameters["covariance_returned_iteration"] == 1
    assert parameters["covariance_returned_loss"] == pytest.approx(info["loss"])
    assert "Returned loss at iteration 1" in output


def test_cli_cov_nonconvergence_exits_nonzero_after_writing_outputs(tmp_path):
    target = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(6, 1.0)
    )
    input_path = tmp_path / "target.npy"
    output_prefix = tmp_path / "nonconverged"
    np.save(input_path, target)

    result = CliRunner().invoke(
        cli_main,
        [
            str(input_path),
            str(output_prefix),
            "--input-type",
            "dmap",
            "--input-format",
            "npy",
            "--method",
            "COV",
            "--gaussian-noise-variance",
            "0.1",
            "--iteration",
            "1",
            "--no-xyzs",
            "--quiet",
        ],
    )

    assert result.exit_code == 1
    assert "COV did not converge" in result.output
    assert "nonconverged partial results" in result.output
    assert (tmp_path / "nonconverged_run_parameters.csv").exists()
    assert (tmp_path / "nonconverged_iteration_series.csv").exists()
    assert (tmp_path / "nonconverged_dmap_final.txt").exists()
    assert (tmp_path / "nonconverged_connectivity_matrix.txt").exists()
    parameters = pd.read_csv(
        tmp_path / "nonconverged_run_parameters.csv"
    ).set_index("parameter")["value"]
    assert parameters["covariance_converged"] == "False"


def test_remove_fully_missing_loci_requires_ignore_missing_data():
    """Removing fully missing loci is only valid when ignore_missing_data is enabled."""
    n = 5
    A_true = HippsDimes.construct_connectivity_matrix_rouse(n, 1.0)
    dmap_target = HippsDimes.a2dmap_theory(A_true)

    with pytest.raises(ValueError, match="remove_fully_missing_loci=True requires ignore_missing_data=True"):
        HippsDimes.run_optimization(
            input_matrix=dmap_target,
            input_type="dmap",
            input_format="text",
            remove_fully_missing_loci=True,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


def test_repair_fully_missing_loci_requires_ignore_missing_data():
    target = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    )

    with pytest.raises(
        ValueError,
        match="repair_fully_missing_loci=True requires ignore_missing_data=True",
    ):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="dmap",
            input_format="text",
            repair_fully_missing_loci=True,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


def test_fully_missing_locus_policy_must_be_explicit():
    target = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    )
    target[2, :] = np.nan
    target[:, 2] = np.nan
    target[2, 2] = 0.0

    with pytest.raises(ValueError, match="explicitly set"):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="dmap",
            input_format="text",
            ignore_missing_data=True,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


def test_fully_missing_locus_repair_and_removal_are_mutually_exclusive():
    target = HippsDimes.a2dmap_theory(
        HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    )

    with pytest.raises(ValueError, match="mutually exclusive"):
        HippsDimes.run_optimization(
            input_matrix=target,
            input_type="dmap",
            input_format="text",
            ignore_missing_data=True,
            repair_fully_missing_loci=True,
            remove_fully_missing_loci=True,
            no_xyzs=True,
            verbose=False,
            show_progress=False,
        )


def test_run_optimization_remove_fully_missing_loci_reduces_problem_and_records_mapping():
    """Removing fully missing loci should reduce the optimization problem and keep index mapping."""
    A_true = HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    cmap = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    cmap_missing = cmap.copy()
    cmap_missing[2, :] = 0.0
    cmap_missing[:, 2] = 0.0
    np.fill_diagonal(cmap_missing, 1.0)

    results = HippsDimes.run_optimization(
        input_matrix=cmap_missing,
        connectivity_matrix=A_true,
        input_type="cmap",
        input_format="text",
        ignore_missing_data=True,
        remove_fully_missing_loci=True,
        method="IS",
        iteration=5,
        learning_rate=5.0,
        no_xyzs=True,
        verbose=False,
        show_progress=False,
    )

    run_parameters = dict(
        zip(
            results["run_parameters"]["parameter"],
            results["run_parameters"]["value"],
        )
    )

    assert results["connectivity_matrix"].shape == (4, 4)
    assert results["dmap_final"].shape == (4, 4)
    assert results["cmap_final"].shape == (4, 4)
    assert results["kept_loci"] == [0, 1, 3, 4]
    assert results["removed_fully_missing_loci"] == [2]
    assert run_parameters["matrix_rows_original"] == 5
    assert run_parameters["matrix_rows"] == 4
    assert bool(run_parameters["remove_fully_missing_loci"]) is True
    assert json.loads(run_parameters["removed_fully_missing_loci"]) == [2]
    assert run_parameters["removed_fully_missing_loci_count"] == 1
    assert run_parameters["nearest_neighbor_repair_count"] == 0
    assert json.loads(run_parameters["remaining_fully_missing_loci"]) == []


def test_cli_remove_fully_missing_loci_writes_reduced_results_pickle(tmp_path):
    """CLI should expose fully missing locus removal through the same API path."""
    A_true = HippsDimes.construct_connectivity_matrix_rouse(5, 1.0)
    cmap = HippsDimes.a2cmap_theory(A_true, rc=1.0)
    cmap_missing = cmap.copy()
    cmap_missing[2, :] = 0.0
    cmap_missing[:, 2] = 0.0
    np.fill_diagonal(cmap_missing, 1.0)

    cmap_path = tmp_path / "cmap_missing.npy"
    output_prefix = tmp_path / "remove_missing_cli"
    np.save(cmap_path, cmap_missing)

    runner = CliRunner()
    result = runner.invoke(
        cli_main,
        [
            str(cmap_path),
            str(output_prefix),
            "--input-type",
            "cmap",
            "--input-format",
            "npy",
            "--iteration",
            "5",
            "--learning-rate",
            "5.0",
            "--no-xyzs",
            "--ignore-missing-data",
            "--remove-fully-missing-loci",
            "--save-pickle",
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output

    with open(tmp_path / "remove_missing_cli_HIPPS_DIMES_results.pkl", "rb") as fin:
        payload = pickle.load(fin)

    run_parameters = dict(
        zip(
            payload["run_parameters"]["parameter"],
            payload["run_parameters"]["value"],
        )
    )

    assert payload["connectivity_matrix"].shape == (4, 4)
    assert payload["kept_loci"] == [0, 1, 3, 4]
    assert payload["removed_fully_missing_loci"] == [2]
    assert json.loads(run_parameters["removed_fully_missing_loci"]) == [2]
    assert run_parameters["removed_fully_missing_loci_count"] == 1


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
