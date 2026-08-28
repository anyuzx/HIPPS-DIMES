"""Scalable PDHG-to-entropy-whitened-Newton hybrid COV solver."""

from __future__ import annotations

import time
import warnings

import numpy as np

from .covariance_newton_whitened import (
    fit_gaussian_noise_covariance_newton_whitened,
)
from .covariance_pdhg import (
    _combine_hybrid_histories,
    fit_gaussian_noise_covariance_pdhg,
)


def fit_gaussian_noise_covariance_hybrid_whitened_newton(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=20_000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    handoff_relative_tolerance=1e-3,
    save_steps=None,
    progress_callback=None,
    pdhg_options=None,
    newton_options=None,
):
    """Run variance-whitened PDHG and polish with whitened Newton-PCG.

    The PDHG physical Gram matrix is passed directly to Newton. The handoff
    skips scalar calibration and verifies that the Newton objective reproduces
    the PDHG objective before any Newton update. A capped Newton-readiness
    probe returns the untouched PDHG solution if the linear solve or KKT-safe
    trial step is not yet useful.
    """
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if (
        not np.isfinite(handoff_relative_tolerance)
        or handoff_relative_tolerance <= 0.0
    ):
        raise ValueError(
            "handoff_relative_tolerance must be positive and finite"
        )
    if (
        relative_tolerance > 0.0
        and handoff_relative_tolerance < relative_tolerance
    ):
        raise ValueError(
            "handoff_relative_tolerance must not be below relative_tolerance"
        )

    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer))
        or step < 1
        or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    pdhg_options = dict(pdhg_options or {})
    newton_options = dict(newton_options or {})
    for reserved in (
        "max_iterations",
        "relative_tolerance",
        "absolute_tolerance",
        "save_steps",
        "progress_callback",
        "initial_connectivity",
        "use_gpu",
    ):
        if reserved in pdhg_options:
            raise ValueError(
                f"pdhg_options must not override reserved key {reserved!r}"
            )
    for reserved in (
        "max_iterations",
        "relative_tolerance",
        "absolute_tolerance",
        "save_steps",
        "progress_callback",
        "initial_gram",
        "initial_connectivity",
        "calibrate_initial_scale",
        "expected_initial_objective",
        "expected_initial_relative_kkt",
        "use_gpu",
    ):
        if reserved in newton_options:
            raise ValueError(
                f"newton_options must not override reserved key {reserved!r}"
            )

    def phase_callback(phase, offset):
        if progress_callback is None:
            return None

        def report(event):
            event = dict(event)
            event["optimizer"] = "hybrid_whitened_newton"
            event["hybrid_phase"] = phase
            if event.get("stage") == "covariance_optimization":
                event["phase_iteration"] = int(event["iteration"])
                event["phase_total"] = int(event["total"])
                event["iteration"] = offset + int(event["iteration"])
                event["total"] = int(max_iterations)
                event["method"] = "COV-hybrid-whitened-newton"
            elif event.get("stage") == "covariance_newton_cg":
                event["global_outer_iteration"] = (
                    offset + int(event.get("outer_iteration", 0))
                )
            progress_callback(event)

        return report

    start_time = time.perf_counter()
    pdhg_start = time.perf_counter()
    pdhg_result = fit_gaussian_noise_covariance_pdhg(
        squared_distances,
        noise_variance,
        relative_noise_std=relative_noise_std,
        initial_connectivity=initial_connectivity,
        use_gpu=use_gpu,
        max_iterations=max_iterations,
        relative_tolerance=handoff_relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        save_steps=sorted(save_steps_set),
        progress_callback=phase_callback("pdhg", 0),
        **pdhg_options,
    )
    pdhg_wall_seconds = time.perf_counter() - pdhg_start
    fitted, gram, connectivity, pdhg_info = pdhg_result
    pdhg_iterations = int(pdhg_info["iterations"])
    handoff_residual = float(
        pdhg_info["relative_eliminated_kkt_residual"]
    )
    handoff = {
        "reached": handoff_residual <= handoff_relative_tolerance,
        "relative_tolerance": float(handoff_relative_tolerance),
        "relative_eliminated_kkt_residual": handoff_residual,
        "objective": float(pdhg_info["objective"]),
        "iterations": pdhg_iterations,
        "returned_iteration": int(pdhg_info["returned_iteration"]),
        "direct_gram": True,
        "scalar_calibration": False,
    }
    remaining_iterations = max_iterations - pdhg_iterations
    if not handoff["reached"] or remaining_iterations <= 0:
        status = (
            "pdhg_handoff_not_reached"
            if not handoff["reached"]
            else "iteration_budget_exhausted_at_handoff"
        )
        info = dict(pdhg_info)
        info.update(
            {
                "algorithm": "hybrid_whitened_newton",
                "converged": False,
                "status": status,
                "message": (
                    "PDHG did not reach the Newton handoff tolerance"
                    if not handoff["reached"]
                    else "PDHG exhausted the total iteration budget at handoff"
                ),
                "handoff": handoff,
                "phase_iterations": {
                    "pdhg": pdhg_iterations,
                    "newton": 0,
                },
                "phase_wall_seconds": {
                    "pdhg": pdhg_wall_seconds,
                    "newton": 0.0,
                },
                "history": _combine_hybrid_histories(
                    pdhg_info["history"],
                    {"iteration": np.asarray([], dtype=np.int64)},
                ),
                "wall_seconds": time.perf_counter() - start_time,
            }
        )
        warnings.warn(
            "hybrid solver stopped before whitened Newton refinement "
            f"(status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
        return fitted, gram, connectivity, info

    newton_save_steps = sorted(
        step - pdhg_iterations
        for step in save_steps_set
        if step > pdhg_iterations
    )
    newton_start = time.perf_counter()
    newton_result = fit_gaussian_noise_covariance_newton_whitened(
        squared_distances,
        noise_variance,
        relative_noise_std=relative_noise_std,
        initial_gram=gram,
        calibrate_initial_scale=False,
        expected_initial_objective=pdhg_info["objective"],
        expected_initial_relative_kkt=handoff_residual,
        use_gpu=use_gpu,
        max_iterations=remaining_iterations,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        save_steps=newton_save_steps,
        progress_callback=phase_callback("newton", pdhg_iterations),
        **newton_options,
    )
    newton_wall_seconds = time.perf_counter() - newton_start
    fitted, gram, connectivity, newton_info = newton_result
    newton_iterations = int(newton_info["iterations"])

    combined_history = _combine_hybrid_histories(
        pdhg_info["history"], newton_info["history"]
    )
    connectivity_at_steps = dict(
        pdhg_info["connectivity_matrix_at_steps"]
    )
    connectivity_at_steps.update(
        {
            pdhg_iterations + int(step): matrix
            for step, matrix in newton_info[
                "connectivity_matrix_at_steps"
            ].items()
        }
    )

    readiness_failed = newton_info["status"] == "newton_readiness_failed"
    if newton_info["converged"]:
        status = "optimality_tolerance"
        message = (
            "PDHG handoff and entropy-whitened Newton refinement converged"
        )
    elif readiness_failed:
        status = "newton_readiness_failed"
        message = (
            "PDHG reached the handoff threshold, but the capped whitened "
            "Newton readiness probe did not produce a KKT-improving step; "
            "the returned model remains the PDHG checkpoint"
        )
    else:
        status = f"newton_{newton_info['status']}"
        message = f"whitened Newton refinement stopped: {newton_info['message']}"

    info = dict(newton_info)
    info.update(
        {
            "converged": bool(newton_info["converged"]),
            "status": status,
            "message": message,
            "iterations": pdhg_iterations + newton_iterations,
            "returned_iteration": pdhg_iterations + newton_iterations,
            "algorithm": "hybrid_whitened_newton",
            "coordinate_parameterization": (
                "centered_pdhg_then_internal_entropy_whitened_svec"
            ),
            "handoff_relative_tolerance": float(
                handoff_relative_tolerance
            ),
            "handoff": handoff,
            "phase_iterations": {
                "pdhg": pdhg_iterations,
                "newton": newton_iterations,
            },
            "phase_wall_seconds": {
                "pdhg": pdhg_wall_seconds,
                "newton": newton_wall_seconds,
            },
            "pdhg": pdhg_info["pdhg"],
            "weighted_operator_norm": pdhg_info[
                "weighted_operator_norm"
            ],
            "initialization": pdhg_info["initialization"],
            "newton_initialization": newton_info["initialization"],
            "history": combined_history,
            "connectivity_matrix_at_steps": connectivity_at_steps,
            "wall_seconds": time.perf_counter() - start_time,
        }
    )
    if not info["converged"]:
        warnings.warn(
            "fit_gaussian_noise_covariance_hybrid_whitened_newton stopped "
            f"without satisfying the final KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted, gram, connectivity, info


fit_gaussian_noise_covariance_hybrid_scalable = (
    fit_gaussian_noise_covariance_hybrid_whitened_newton
)
