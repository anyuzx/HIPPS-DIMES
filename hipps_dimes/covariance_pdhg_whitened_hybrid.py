"""Hybrid wrapper using variance-whitened PDHG before Newton-CG."""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import covariance_pdhg as _base
from . import covariance_pdhg_whitened as _white
from . import numerics as _num


def fit_gaussian_noise_covariance_hybrid_whitened(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    handoff_relative_tolerance=1e-3,
    save_steps=None,
    progress_callback=None,
    pdhg_options=None,
    newton_options=None,
):
    """Use variance-whitened PDHG globally and existing Newton-CG locally."""
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    handoff_relative_tolerance = float(handoff_relative_tolerance)
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
            "handoff_relative_tolerance must not be smaller than the final "
            "relative_tolerance"
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
    reserved = {
        "relative_noise_std",
        "initial_connectivity",
        "use_gpu",
        "max_iterations",
        "relative_tolerance",
        "absolute_tolerance",
        "save_steps",
        "progress_callback",
    }
    conflict = reserved.intersection(pdhg_options)
    if conflict:
        raise ValueError(
            "pdhg_options contains reserved keys: "
            + ", ".join(sorted(conflict))
        )
    conflict = reserved.intersection(newton_options)
    if conflict:
        raise ValueError(
            "newton_options contains reserved keys: "
            + ", ".join(sorted(conflict))
        )

    def phase_callback(phase, offset):
        if progress_callback is None:
            return None

        def report(event):
            event = dict(event)
            event["optimizer"] = "hybrid_whitened"
            event["phase"] = phase
            if event.get("stage") == "covariance_optimization":
                event["phase_iteration"] = int(event["iteration"])
                event["phase_total"] = int(event["total"])
                event["iteration"] = offset + int(event["iteration"])
                event["total"] = int(max_iterations)
                event["method"] = "COV-hybrid-whitened"
            progress_callback(event)

        return report

    start_time = time.perf_counter()
    pdhg_start_time = time.perf_counter()
    fitted, gram, connectivity, pdhg_info = (
        _white.fit_gaussian_noise_covariance_pdhg_whitened(
            squared_distances,
            noise_variance,
            relative_noise_std=relative_noise_std,
            initial_connectivity=initial_connectivity,
            use_gpu=use_gpu,
            max_iterations=max_iterations,
            relative_tolerance=handoff_relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            save_steps=sorted(save_steps_set),
            progress_callback=phase_callback("pdhg_whitened", 0),
            **pdhg_options,
        )
    )
    pdhg_wall_seconds = time.perf_counter() - pdhg_start_time
    pdhg_iterations = int(pdhg_info["iterations"])
    handoff_residual = float(
        pdhg_info["relative_eliminated_kkt_residual"]
    )
    handoff = {
        "reached": handoff_residual <= handoff_relative_tolerance,
        "relative_tolerance": handoff_relative_tolerance,
        "relative_eliminated_kkt_residual": handoff_residual,
        "objective": float(pdhg_info["objective"]),
        "iterations": pdhg_iterations,
        "returned_iteration": int(pdhg_info["returned_iteration"]),
    }
    remaining_iterations = max_iterations - pdhg_iterations
    if not handoff["reached"] or remaining_iterations <= 0:
        status = (
            "pdhg_handoff_not_reached"
            if not handoff["reached"]
            else "iteration_budget_exhausted_at_handoff"
        )
        message = (
            "variance-whitened PDHG did not reach the hybrid handoff tolerance"
            if not handoff["reached"]
            else "PDHG reached handoff at the end of the iteration budget"
        )
        info = dict(pdhg_info)
        info.update(
            {
                "algorithm": "hybrid_whitened",
                "converged": False,
                "status": status,
                "message": message,
                "relative_tolerance": float(relative_tolerance),
                "handoff_relative_tolerance": handoff_relative_tolerance,
                "handoff": handoff,
                "phase_iterations": {
                    "pdhg_whitened": pdhg_iterations,
                    "newton": 0,
                },
                "phase_wall_seconds": {
                    "pdhg_whitened": pdhg_wall_seconds,
                    "newton": 0.0,
                },
                "history": _base._combine_hybrid_histories(
                    pdhg_info["history"],
                    {"iteration": np.asarray([], dtype=np.int64)},
                ),
                "wall_seconds": time.perf_counter() - start_time,
            }
        )
        warnings.warn(
            "fit_gaussian_noise_covariance_hybrid_whitened stopped before "
            f"Newton refinement (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
        return fitted, gram, connectivity, info

    newton_save_steps = sorted(
        step - pdhg_iterations
        for step in save_steps_set
        if step > pdhg_iterations
    )
    newton_start_time = time.perf_counter()
    fitted, gram, connectivity, newton_info = (
        _num.fit_gaussian_noise_covariance(
            squared_distances,
            noise_variance,
            relative_noise_std=relative_noise_std,
            initial_connectivity=connectivity,
            use_gpu=use_gpu,
            max_iterations=remaining_iterations,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            save_steps=newton_save_steps,
            progress_callback=phase_callback(
                "newton", pdhg_iterations
            ),
            **newton_options,
        )
    )
    newton_wall_seconds = time.perf_counter() - newton_start_time
    newton_iterations = int(newton_info["iterations"])

    (
        _,
        pair_i,
        pair_j,
        target_pairs,
        pair_variance,
        _,
        _,
    ) = _num._validate_gaussian_covariance_inputs(
        squared_distances, noise_variance, relative_noise_std
    )
    independent = _base._independent_eliminated_kkt_residuals(
        gram,
        target_pairs,
        pair_variance,
        pair_i,
        pair_j,
    )
    (
        independent_norm,
        independent_scale,
        independent_relative,
        _,
        independent_residual,
    ) = independent
    independent_converged = independent_norm <= (
        absolute_tolerance + relative_tolerance * independent_scale
    )
    internal_converged = bool(newton_info["converged"])
    converged = internal_converged and independent_converged
    if converged:
        status = "optimality_tolerance"
        message = (
            "whitened PDHG handoff, Newton refinement, and independent KKT "
            "certificate reached"
        )
    elif internal_converged:
        status = "independent_kkt_failed"
        message = (
            "Newton reached its internal tolerance, but the returned Gram "
            "matrix failed the independent eliminated KKT certificate"
        )
    else:
        status = f"newton_{newton_info['status']}"
        message = f"Newton refinement stopped: {newton_info['message']}"

    combined_history = _base._combine_hybrid_histories(
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
    info = dict(newton_info)
    info.update(
        {
            "converged": bool(converged),
            "status": status,
            "message": message,
            "iterations": pdhg_iterations + newton_iterations,
            "returned_iteration": pdhg_iterations + newton_iterations,
            "relative_gradient_norm": float(independent_relative),
            "relative_eliminated_kkt_residual": float(
                independent_relative
            ),
            "stationarity_residual_norm": float(independent_norm),
            "stationarity_residual_scale": float(independent_scale),
            "relative_stationarity_residual": float(independent_relative),
            "maximum_absolute_stationarity_residual": float(
                np.max(np.abs(independent_residual))
            ),
            "independent_kkt_converged": bool(independent_converged),
            "relative_tolerance": float(relative_tolerance),
            "absolute_tolerance": float(absolute_tolerance),
            "handoff_relative_tolerance": handoff_relative_tolerance,
            "handoff": handoff,
            "phase_iterations": {
                "pdhg_whitened": pdhg_iterations,
                "newton": newton_iterations,
            },
            "phase_wall_seconds": {
                "pdhg_whitened": pdhg_wall_seconds,
                "newton": newton_wall_seconds,
            },
            "initialization": pdhg_info["initialization"],
            "algorithm": "hybrid_whitened",
            "coordinate_parameterization": "centered_hybrid",
            "history": combined_history,
            "connectivity_matrix_at_steps": connectivity_at_steps,
            "wall_seconds": time.perf_counter() - start_time,
            "pdhg_whitened": pdhg_info["pdhg"],
            "weighted_operator_norm": pdhg_info[
                "weighted_operator_norm"
            ],
        }
    )
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_hybrid_whitened stopped without "
            f"the final KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted, gram, connectivity, info
