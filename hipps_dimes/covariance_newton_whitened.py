"""Entropy-whitened Newton-CG for the Gaussian noise-aware COV objective.

The solver minimizes

    -3/2 logdet(B) + 1/2 ||V^-1/2 (D(B) - d)||^2

on the positive-definite internal Gram matrix. Newton directions are written
as ``delta_B = L X L.T`` for ``B = L L.T``. In the variable ``X`` the
entropy Hessian is exactly ``3/2 I``; the remaining data curvature is applied
matrix-free. PCG operates on the orthonormal symmetric ``svec`` coordinates
and uses the exact diagonal of the whitened data Hessian.

This module is intentionally separate from the legacy Newton implementation so
large-system handoff behavior can be benchmarked without changing the existing
COV solver.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import numerics as _num


_ENTROPY_COEFFICIENT = 1.5
_SQRT_TWO = np.sqrt(2.0)


def _scalar(value) -> float:
    return float(value.item()) if hasattr(value, "item") else float(value)


def _norm(value, xp) -> float:
    return _scalar(xp.sqrt(xp.sum(value * value)))


def _symmetrize(value):
    return 0.5 * (value + value.T)


def _relative_gradient_state(state, xp):
    gradient = state[1]
    entropy_gradient = state[4][2]
    data_gradient = state[4][3]
    gradient_norm = _norm(gradient, xp)
    gradient_scale = max(
        1.0,
        _norm(entropy_gradient, xp),
        _norm(data_gradient, xp),
    )
    return gradient_norm, gradient_scale, gradient_norm / gradient_scale


def _forcing_tolerance(
    relative_kkt: float,
    minimum: float,
    maximum: float,
    coefficient: float,
) -> float:
    """Return an inexact-Newton forcing tolerance.

    ``minimum`` is the final linear-solve accuracy, while ``maximum`` permits
    inexpensive solves far from the solution. This intentionally does not use
    the legacy ``min(cg_relative_tolerance, adaptive_value)`` rule, which made
    every solve at least as strict as the final tolerance.
    """
    proposed = coefficient * np.sqrt(max(float(relative_kkt), 0.0))
    return float(min(maximum, max(minimum, proposed)))


def _svec_index_arrays(size: int, xp):
    diagonal = np.arange(size, dtype=np.int32)
    off_i, off_j = np.triu_indices(size, k=1)
    off_i = off_i.astype(np.int32, copy=False)
    off_j = off_j.astype(np.int32, copy=False)
    if xp is np:
        return diagonal, off_i, off_j
    return (
        xp.asarray(diagonal, dtype=xp.int32),
        xp.asarray(off_i, dtype=xp.int32),
        xp.asarray(off_j, dtype=xp.int32),
    )


def _svec(matrix, diagonal, off_i, off_j, xp):
    """Orthogonally vectorize a symmetric matrix."""
    return xp.concatenate(
        (
            matrix[diagonal, diagonal],
            _SQRT_TWO * matrix[off_i, off_j],
        )
    )


def _smat(vector, size, diagonal, off_i, off_j, xp):
    """Inverse of :func:`_svec`."""
    matrix = xp.zeros((size, size), dtype=xp.float64)
    matrix[diagonal, diagonal] = vector[:size]
    off_values = vector[size:] / _SQRT_TWO
    matrix[off_i, off_j] = off_values
    matrix[off_j, off_i] = off_values
    return matrix


def _whitened_svec_preconditioner(
    whitened_data_diagonal,
    diagonal,
    off_i,
    off_j,
    xp,
):
    """Return the exact Hessian diagonal in orthonormal svec coordinates."""
    diagonal_values = (
        _ENTROPY_COEFFICIENT
        + whitened_data_diagonal[diagonal, diagonal]
    )
    off_diagonal_values = (
        _ENTROPY_COEFFICIENT
        + 2.0 * whitened_data_diagonal[off_i, off_j]
    )
    result = xp.concatenate((diagonal_values, off_diagonal_values))
    return xp.maximum(result, np.finfo(np.float64).tiny)


def _pcg_svec(
    operator,
    right_hand_side,
    preconditioner_diagonal,
    relative_tolerance,
    max_iterations,
    *,
    xp,
    progress_callback=None,
    progress_interval=25,
    true_residual_interval=50,
    outer_iteration=1,
):
    """Preconditioned CG with explicit progress and a hard iteration cap."""
    solution = xp.zeros_like(right_hand_side)
    residual = right_hand_side.copy()
    right_norm = _norm(right_hand_side, xp)
    if right_norm == 0.0:
        return solution, {
            "iterations": 0,
            "relative_residual": 0.0,
            "best_relative_residual": 0.0,
            "converged": True,
            "termination_reason": "zero_right_hand_side",
            "hessian_actions": 0,
            "wall_seconds": 0.0,
        }

    tolerance = float(relative_tolerance) * right_norm
    preconditioned = residual / preconditioner_diagonal
    direction = preconditioned.copy()
    residual_product = _scalar(xp.sum(residual * preconditioned))
    best_relative = 1.0
    converged = False
    termination_reason = "maximum_iterations"
    hessian_actions = 0
    iteration = 0
    start_time = time.perf_counter()

    def emit(relative_residual, force=False):
        if progress_callback is None:
            return
        if not force and iteration % int(progress_interval) != 0:
            return
        progress_callback(
            {
                "stage": "covariance_newton_cg",
                "phase": "newton_cg",
                "outer_iteration": int(outer_iteration),
                "iteration": int(iteration),
                "total": int(max_iterations),
                "cg_relative_residual": float(relative_residual),
                "cg_best_relative_residual": float(best_relative),
                "cg_target": float(relative_tolerance),
                "hessian_actions": int(hessian_actions),
                "wall_seconds": time.perf_counter() - start_time,
                "method": "COV-Newton-whitened",
                "general_method": "covariance_optimization",
                "use_gpu": bool(xp is not np),
                "noisy": True,
            }
        )

    emit(1.0, force=True)
    for iteration in range(1, int(max_iterations) + 1):
        operator_direction = operator(direction)
        hessian_actions += 1
        curvature = _scalar(xp.sum(direction * operator_direction))
        if not np.isfinite(curvature) or curvature <= 0.0:
            termination_reason = "nonpositive_curvature"
            break
        step = residual_product / curvature
        solution += step * direction
        residual -= step * operator_direction

        replaced_residual = False
        if (
            true_residual_interval
            and iteration % int(true_residual_interval) == 0
        ):
            residual = right_hand_side - operator(solution)
            hessian_actions += 1
            replaced_residual = True

        residual_norm = _norm(residual, xp)
        relative_residual = residual_norm / right_norm
        best_relative = min(best_relative, relative_residual)
        emit(relative_residual)
        if residual_norm <= tolerance:
            converged = True
            termination_reason = "relative_tolerance"
            break

        preconditioned = residual / preconditioner_diagonal
        next_product = _scalar(xp.sum(residual * preconditioned))
        if not np.isfinite(next_product) or next_product <= 0.0:
            termination_reason = "invalid_preconditioned_residual"
            break
        if replaced_residual:
            direction = preconditioned.copy()
        else:
            direction = (
                preconditioned
                + (next_product / residual_product) * direction
            )
        residual_product = next_product

    residual = right_hand_side - operator(solution)
    hessian_actions += 1
    relative_residual = _norm(residual, xp) / right_norm
    best_relative = min(best_relative, relative_residual)
    emit(relative_residual, force=True)
    return solution, {
        "iterations": int(iteration),
        "relative_residual": float(relative_residual),
        "best_relative_residual": float(best_relative),
        "converged": bool(converged),
        "termination_reason": termination_reason,
        "hessian_actions": int(hessian_actions),
        "wall_seconds": time.perf_counter() - start_time,
    }


def _provided_gram_initialization(
    initial_gram,
    solver_basis,
    *,
    xp,
):
    """Project one physical centered Gram matrix into the internal basis."""
    n = solver_basis.shape[0]
    gram_cpu = np.asarray(initial_gram, dtype=np.float64)
    if gram_cpu.shape != (n, n) or not np.all(np.isfinite(gram_cpu)):
        raise ValueError("initial_gram must be a finite NxN matrix")
    if not np.allclose(gram_cpu, gram_cpu.T, rtol=1e-10, atol=1e-11):
        raise ValueError("initial_gram must be symmetric")
    gram_cpu = _symmetrize(gram_cpu)
    gram_scale = max(float(np.max(np.abs(gram_cpu))), 1.0)
    row_sum_error = float(np.max(np.abs(np.sum(gram_cpu, axis=1))))
    if row_sum_error > 1e-8 * gram_scale:
        raise ValueError("initial_gram must be centered (rows sum to zero)")

    gram = xp.asarray(gram_cpu, dtype=xp.float64)
    reduced = _symmetrize(solver_basis.T @ gram @ solver_basis)
    try:
        if xp is np:
            np.linalg.cholesky(reduced)
        else:
            with _num.cupyx.errstate(linalg="raise"):
                xp.linalg.cholesky(reduced)
    except np.linalg.LinAlgError as error:
        raise ValueError(
            "initial_gram must be positive definite on the internal subspace"
        ) from error

    projected = solver_basis @ reduced @ solver_basis.T
    projection_error = _norm(projected - gram, xp) / max(
        _norm(gram, xp), np.finfo(np.float64).tiny
    )
    return reduced, {
        "kind": "provided_gram",
        "scalar_calibration_applied": False,
        "relative_centered_projection_error": float(projection_error),
        "initial_row_sum_error": float(row_sum_error),
        "backend": "gpu" if xp is not np else "cpu",
        "wall_seconds": 0.0,
    }


def _build_whitened_preconditioner(
    solver_basis,
    cholesky_factor,
    solver_pair_i,
    solver_pair_j,
    solver_inverse_variance,
    diagonal,
    off_i,
    off_j,
    *,
    xp,
    use_gpu,
    progress_callback,
    outer_iteration,
):
    """Build the exact svec diagonal for the entropy-whitened Hessian."""
    whitened_basis = solver_basis @ cholesky_factor

    def report(event):
        if progress_callback is None:
            return
        event = dict(event)
        event["phase"] = "newton_preconditioner"
        event["outer_iteration"] = int(outer_iteration)
        event["method"] = "COV-Newton-whitened"
        progress_callback(event)

    data_diagonal, setup_seconds = (
        _num._gaussian_covariance_data_preconditioner_diagonal(
            whitened_basis,
            solver_pair_i,
            solver_pair_j,
            solver_inverse_variance,
            progress_callback=report if progress_callback is not None else None,
            array_module=xp,
            use_gpu=use_gpu,
        )
    )
    svec_diagonal = _whitened_svec_preconditioner(
        data_diagonal,
        diagonal,
        off_i,
        off_j,
        xp,
    )
    if progress_callback is not None:
        progress_callback(
            {
                "stage": "covariance_newton_preconditioner",
                "phase": "newton_preconditioner",
                "outer_iteration": int(outer_iteration),
                "iteration": 1,
                "total": 1,
                "completed": True,
                "wall_seconds": float(setup_seconds),
                "method": "COV-Newton-whitened",
                "general_method": "covariance_optimization",
                "use_gpu": bool(use_gpu),
                "noisy": True,
            }
        )
    return whitened_basis, svec_diagonal, float(setup_seconds)


def fit_gaussian_noise_covariance_newton_whitened(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_gram=None,
    initial_connectivity=None,
    calibrate_initial_scale=None,
    expected_initial_objective=None,
    expected_initial_relative_kkt=None,
    handoff_objective_relative_tolerance=1e-8,
    use_gpu=False,
    max_iterations=20,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    cg_forcing_min=1e-6,
    cg_forcing_max=0.1,
    cg_forcing_coefficient=0.1,
    cg_max_iterations=None,
    cg_progress_interval=25,
    cg_true_residual_interval=50,
    readiness_probe=True,
    readiness_cg_tolerance=0.1,
    readiness_cg_max_iterations=200,
    readiness_kkt_reduction=0.95,
    truncated_kkt_reduction=0.9,
    maximum_kkt_growth=1.0,
    line_search_max_iterations=40,
    save_steps=None,
    progress_callback=None,
):
    """Fit the COV objective with centered entropy-whitened Newton-PCG.

    When ``initial_gram`` is supplied it is handed to Newton directly and is
    *not* scalar calibrated by default. This is the intended PDHG-to-Newton
    path. The first Krylov solve can be capped as a readiness probe; failure to
    produce a KKT-improving step returns the untouched handoff point quickly.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "whitened Newton GPU optimization was requested, but CuPy and an "
            "accessible CUDA GPU are not available"
        )
    if initial_gram is not None and initial_connectivity is not None:
        raise ValueError(
            "provide at most one of initial_gram and initial_connectivity"
        )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    for name, value in (
        ("cg_forcing_min", cg_forcing_min),
        ("cg_forcing_max", cg_forcing_max),
        ("cg_forcing_coefficient", cg_forcing_coefficient),
        ("readiness_cg_tolerance", readiness_cg_tolerance),
        ("readiness_kkt_reduction", readiness_kkt_reduction),
        ("truncated_kkt_reduction", truncated_kkt_reduction),
        ("maximum_kkt_growth", maximum_kkt_growth),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if cg_forcing_min >= 1.0 or cg_forcing_max >= 1.0:
        raise ValueError("CG forcing tolerances must be below one")
    if cg_forcing_min > cg_forcing_max:
        raise ValueError("cg_forcing_min must not exceed cg_forcing_max")
    if readiness_cg_tolerance >= 1.0:
        raise ValueError("readiness_cg_tolerance must be below one")
    if readiness_kkt_reduction > 1.0 or truncated_kkt_reduction > 1.0:
        raise ValueError("KKT reduction factors must not exceed one")
    if maximum_kkt_growth < 1.0:
        raise ValueError("maximum_kkt_growth must be at least one")

    (
        observed,
        pair_i,
        pair_j,
        target_pairs,
        pair_variance,
        noise_model,
        noise_parameter,
    ) = _num._validate_gaussian_covariance_inputs(
        squared_distances, noise_variance, relative_noise_std
    )
    n = observed.shape[0]
    n_modes = n - 1
    if cg_max_iterations is None:
        cg_max_iterations = min(500, max(50, 2 * n_modes))
    if not isinstance(cg_max_iterations, (int, np.integer)) or cg_max_iterations <= 0:
        raise ValueError("cg_max_iterations must be a positive integer")
    if (
        not isinstance(readiness_cg_max_iterations, (int, np.integer))
        or readiness_cg_max_iterations <= 0
    ):
        raise ValueError("readiness_cg_max_iterations must be positive")
    if (
        not isinstance(line_search_max_iterations, (int, np.integer))
        or line_search_max_iterations <= 0
    ):
        raise ValueError("line_search_max_iterations must be positive")

    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer))
        or step < 1
        or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    basis = _num._centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = (
        _num._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    if use_gpu:
        xp = _num.cp
        solver_basis = xp.asarray(basis, dtype=xp.float64)
        solver_pair_i = xp.asarray(pair_i, dtype=xp.int32)
        solver_pair_j = xp.asarray(pair_j, dtype=xp.int32)
        solver_target_pairs = xp.asarray(target_pairs, dtype=xp.float64)
        solver_inverse_variance = xp.asarray(
            inverse_variance, dtype=xp.float64
        )
        solver_target_matrix = xp.asarray(target_matrix, dtype=xp.float64)
        solver_inverse_variance_matrix = xp.asarray(
            inverse_variance_matrix, dtype=xp.float64
        )
    else:
        xp = np
        solver_basis = basis
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target_pairs = target_pairs
        solver_inverse_variance = inverse_variance
        solver_target_matrix = target_matrix
        solver_inverse_variance_matrix = inverse_variance_matrix

    initialization_start = time.perf_counter()
    if initial_gram is not None:
        reduced_gram, initialization_info = _provided_gram_initialization(
            initial_gram,
            solver_basis,
            xp=xp,
        )
        if calibrate_initial_scale is None:
            calibrate_initial_scale = False
    else:
        reduced_cpu, initialization_info = (
            _num._initialize_gaussian_reduced_gram(
                observed,
                pair_i,
                pair_j,
                basis,
                initial_connectivity,
            )
        )
        reduced_gram = xp.asarray(reduced_cpu, dtype=xp.float64)
        if calibrate_initial_scale is None:
            calibrate_initial_scale = True

    if calibrate_initial_scale:
        reduced_gram, calibration = (
            _num._calibrate_gaussian_covariance_initial_scale(
                reduced_gram,
                solver_basis,
                solver_pair_i,
                solver_pair_j,
                solver_target_pairs,
                solver_inverse_variance,
                array_module=xp,
            )
        )
        initialization_info["scalar_calibration"] = calibration
        initialization_info["scalar_calibration_applied"] = True
    else:
        initialization_info["scalar_calibration"] = {
            "method": "disabled_for_direct_gram_handoff",
            "scale_factor": 1.0,
            "connectivity_scale_factor": 1.0,
            "objective_reduction": 0.0,
            "wall_seconds": 0.0,
        }
        initialization_info["scalar_calibration_applied"] = False
    initialization_info["wall_seconds"] = time.perf_counter() - initialization_start

    current_state = _num._gaussian_covariance_objective_gradient(
        reduced_gram,
        solver_basis,
        solver_target_matrix,
        solver_inverse_variance_matrix,
        solver_pair_i,
        solver_pair_j,
        array_module=xp,
    )
    if not np.isfinite(current_state[0]):
        raise RuntimeError("initial Gram matrix is outside the covariance cone")
    initial_objective = float(current_state[0])
    _, _, initial_relative_kkt = _relative_gradient_state(current_state, xp)
    initialization_info["objective"] = initial_objective
    initialization_info["relative_kkt"] = float(initial_relative_kkt)

    if expected_initial_objective is not None:
        objective_mismatch = abs(
            initial_objective - float(expected_initial_objective)
        ) / max(1.0, abs(float(expected_initial_objective)))
        initialization_info["expected_objective"] = float(
            expected_initial_objective
        )
        initialization_info["relative_objective_mismatch"] = float(
            objective_mismatch
        )
        if objective_mismatch > handoff_objective_relative_tolerance:
            raise RuntimeError(
                "Newton direct-Gram handoff does not reproduce the PDHG "
                "objective"
            )
    if expected_initial_relative_kkt is not None:
        initialization_info["expected_relative_kkt"] = float(
            expected_initial_relative_kkt
        )
        initialization_info["relative_kkt_difference"] = float(
            initial_relative_kkt - float(expected_initial_relative_kkt)
        )

    diagonal, off_i, off_j = _svec_index_arrays(n_modes, xp)
    history = {
        "iteration": [],
        "objective": [],
        "negative_entropy_objective": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "gradient_norm": [],
        "relative_gradient_norm": [],
        "relative_eliminated_kkt_residual": [],
        "newton_decrement": [],
        "relative_step": [],
        "step_size": [],
        "line_search_backtracks": [],
        "kkt_before": [],
        "kkt_after": [],
        "kkt_reduction_factor": [],
        "cg_iterations": [],
        "cg_relative_residual": [],
        "cg_best_relative_residual": [],
        "cg_target": [],
        "cg_converged": [],
        "cg_hessian_actions": [],
        "cg_termination_reason": [],
        "readiness_probe": [],
        "accepted_truncated_cg": [],
        "preconditioner_setup_seconds": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "gram_condition_number": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}
    status = "max_iterations"
    message = "maximum number of whitened Newton iterations reached"
    converged = False
    readiness_passed = not readiness_probe
    readiness_info = {
        "enabled": bool(readiness_probe),
        "passed": bool(readiness_passed),
        "cg_target": float(readiness_cg_tolerance),
        "cg_max_iterations": int(readiness_cg_max_iterations),
    }
    total_preconditioner_seconds = 0.0
    total_cg_seconds = 0.0
    start_time = time.perf_counter()
    accepted_iterations = 0

    for outer_iteration in range(1, max_iterations + 1):
        objective, gradient, _, _, _ = current_state
        gradient_norm, gradient_scale, relative_kkt = (
            _relative_gradient_state(current_state, xp)
        )
        if gradient_norm <= (
            absolute_tolerance + relative_tolerance * gradient_scale
        ):
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break

        if use_gpu:
            with _num.cupyx.errstate(linalg="raise"):
                cholesky_factor = xp.linalg.cholesky(reduced_gram)
        else:
            cholesky_factor = np.linalg.cholesky(reduced_gram)

        (
            whitened_basis,
            preconditioner_diagonal,
            preconditioner_seconds,
        ) = _build_whitened_preconditioner(
            solver_basis,
            cholesky_factor,
            solver_pair_i,
            solver_pair_j,
            solver_inverse_variance,
            diagonal,
            off_i,
            off_j,
            xp=xp,
            use_gpu=use_gpu,
            progress_callback=progress_callback,
            outer_iteration=outer_iteration,
        )
        total_preconditioner_seconds += preconditioner_seconds

        transformed_gradient = _symmetrize(
            cholesky_factor.T @ gradient @ cholesky_factor
        )
        right_hand_side = -_svec(
            transformed_gradient, diagonal, off_i, off_j, xp
        )

        def hessian_svec(vector):
            matrix_direction = _smat(
                vector, n_modes, diagonal, off_i, off_j, xp
            )
            data_action = _num._gaussian_covariance_data_hessian_action(
                matrix_direction,
                whitened_basis,
                solver_inverse_variance_matrix,
                array_module=xp,
            )
            result = _symmetrize(
                data_action + _ENTROPY_COEFFICIENT * matrix_direction
            )
            return _svec(result, diagonal, off_i, off_j, xp)

        probe_iteration = readiness_probe and not readiness_passed
        if probe_iteration:
            cg_target = max(
                readiness_cg_tolerance,
                _forcing_tolerance(
                    relative_kkt,
                    cg_forcing_min,
                    cg_forcing_max,
                    cg_forcing_coefficient,
                ),
            )
            cg_limit = min(
                int(cg_max_iterations), int(readiness_cg_max_iterations)
            )
        else:
            cg_target = _forcing_tolerance(
                relative_kkt,
                cg_forcing_min,
                cg_forcing_max,
                cg_forcing_coefficient,
            )
            cg_limit = int(cg_max_iterations)

        solution_vector, cg_info = _pcg_svec(
            hessian_svec,
            right_hand_side,
            preconditioner_diagonal,
            cg_target,
            cg_limit,
            xp=xp,
            progress_callback=progress_callback,
            progress_interval=cg_progress_interval,
            true_residual_interval=cg_true_residual_interval,
            outer_iteration=outer_iteration,
        )
        total_cg_seconds += cg_info["wall_seconds"]
        whitened_direction = _smat(
            solution_vector, n_modes, diagonal, off_i, off_j, xp
        )
        direction = _symmetrize(
            cholesky_factor
            @ whitened_direction
            @ cholesky_factor.T
        )
        directional_derivative = _scalar(xp.sum(gradient * direction))
        if (
            not np.isfinite(directional_derivative)
            or directional_derivative >= 0.0
        ):
            status = "non_descent_direction"
            message = "whitened PCG did not produce an objective descent direction"
            if probe_iteration:
                status = "newton_readiness_failed"
                message = "Newton readiness probe produced no descent direction"
            break

        if probe_iteration:
            required_kkt_factor = readiness_kkt_reduction
        elif cg_info["converged"]:
            required_kkt_factor = maximum_kkt_growth
        else:
            required_kkt_factor = truncated_kkt_reduction

        step_size = 1.0
        candidate = None
        candidate_state = None
        candidate_relative_kkt = np.inf
        backtracks = 0
        for backtracks in range(int(line_search_max_iterations)):
            trial = _symmetrize(reduced_gram + step_size * direction)
            trial_state = _num._gaussian_covariance_objective_gradient(
                trial,
                solver_basis,
                solver_target_matrix,
                solver_inverse_variance_matrix,
                solver_pair_i,
                solver_pair_j,
                array_module=xp,
            )
            if np.isfinite(trial_state[0]):
                _, _, trial_relative_kkt = _relative_gradient_state(
                    trial_state, xp
                )
                objective_ok = float(trial_state[0]) <= (
                    float(objective)
                    + 1e-4 * step_size * directional_derivative
                )
                kkt_limit = max(
                    relative_tolerance,
                    required_kkt_factor * relative_kkt,
                )
                kkt_ok = trial_relative_kkt <= kkt_limit
                if objective_ok and kkt_ok:
                    candidate = trial
                    candidate_state = trial_state
                    candidate_relative_kkt = trial_relative_kkt
                    break
            step_size *= 0.5

        if candidate is None:
            if probe_iteration:
                status = "newton_readiness_failed"
                message = (
                    "Newton readiness probe did not produce an objective- and "
                    "KKT-improving feasible step"
                )
                readiness_info.update(
                    {
                        "passed": False,
                        "cg_converged": bool(cg_info["converged"]),
                        "cg_relative_residual": float(
                            cg_info["relative_residual"]
                        ),
                    }
                )
            else:
                status = "kkt_safeguard_failed"
                message = (
                    "no feasible Armijo step satisfied the KKT safeguard"
                )
            break

        accepted_delta = candidate - reduced_gram
        relative_step = _norm(accepted_delta, xp) / max(
            _norm(reduced_gram, xp), 1.0
        )
        reduced_gram = candidate
        current_state = candidate_state
        accepted_iterations += 1
        if probe_iteration:
            readiness_passed = True
            readiness_info.update(
                {
                    "passed": True,
                    "cg_converged": bool(cg_info["converged"]),
                    "cg_relative_residual": float(
                        cg_info["relative_residual"]
                    ),
                    "kkt_before": float(relative_kkt),
                    "kkt_after": float(candidate_relative_kkt),
                }
            )

        (
            candidate_objective,
            candidate_gradient,
            candidate_pairs,
            candidate_inverse,
            candidate_components,
        ) = candidate_state
        (
            candidate_negative_entropy,
            candidate_data_objective,
            candidate_entropy_gradient,
            candidate_data_gradient,
        ) = candidate_components
        candidate_gradient_norm = _norm(candidate_gradient, xp)
        candidate_gradient_scale = max(
            1.0,
            _norm(candidate_entropy_gradient, xp),
            _norm(candidate_data_gradient, xp),
        )
        residual = candidate_pairs - solver_target_pairs
        residual_squared = xp.square(residual)
        relative_loss = _scalar(
            xp.sqrt(xp.mean(xp.square(residual / solver_target_pairs)))
        )
        weighted_rmse = _scalar(
            xp.sqrt(xp.mean(residual_squared * solver_inverse_variance))
        )
        distance_rmse = _scalar(xp.sqrt(xp.mean(residual_squared)))
        if use_gpu:
            with _num.cupyx.errstate(linalg="raise"):
                gram_eigenvalues = xp.linalg.eigvalsh(reduced_gram)
        else:
            gram_eigenvalues = np.linalg.eigvalsh(reduced_gram)
        logdet = _scalar(xp.sum(xp.log(gram_eigenvalues)))
        entropy = logdet - n_modes * np.log(3.0)
        reduced_connectivity = 3.0 * candidate_inverse
        basis_times_connectivity = solver_basis @ reduced_connectivity
        connectivity_diagonal = -xp.sum(
            basis_times_connectivity * solver_basis, axis=1
        )
        offdiagonal_norm_squared = 0.5 * max(
            _scalar(xp.sum(reduced_connectivity**2))
            - _scalar(xp.sum(connectivity_diagonal**2)),
            0.0,
        )
        connectivity_norm = float(np.sqrt(offdiagonal_norm_squared))

        history["iteration"].append(accepted_iterations)
        history["objective"].append(float(candidate_objective))
        history["negative_entropy_objective"].append(
            float(candidate_negative_entropy)
        )
        history["data_objective"].append(float(candidate_data_objective))
        history["loss"].append(relative_loss)
        history["entropy"].append(entropy)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["gradient_norm"].append(candidate_gradient_norm)
        history["relative_gradient_norm"].append(
            candidate_relative_kkt
        )
        history["relative_eliminated_kkt_residual"].append(
            candidate_relative_kkt
        )
        history["newton_decrement"].append(
            float(np.sqrt(max(-directional_derivative, 0.0)))
        )
        history["relative_step"].append(relative_step)
        history["step_size"].append(step_size)
        history["line_search_backtracks"].append(backtracks)
        history["kkt_before"].append(relative_kkt)
        history["kkt_after"].append(candidate_relative_kkt)
        history["kkt_reduction_factor"].append(
            candidate_relative_kkt / max(
                relative_kkt, np.finfo(np.float64).tiny
            )
        )
        history["cg_iterations"].append(cg_info["iterations"])
        history["cg_relative_residual"].append(
            cg_info["relative_residual"]
        )
        history["cg_best_relative_residual"].append(
            cg_info["best_relative_residual"]
        )
        history["cg_target"].append(cg_target)
        history["cg_converged"].append(cg_info["converged"])
        history["cg_hessian_actions"].append(
            cg_info["hessian_actions"]
        )
        history["cg_termination_reason"].append(
            cg_info["termination_reason"]
        )
        history["readiness_probe"].append(probe_iteration)
        history["accepted_truncated_cg"].append(
            not cg_info["converged"]
        )
        history["preconditioner_setup_seconds"].append(
            preconditioner_seconds
        )
        minimum_eigenvalue = _scalar(gram_eigenvalues[0])
        maximum_eigenvalue = _scalar(gram_eigenvalues[-1])
        history["minimum_internal_gram_eigenvalue"].append(
            minimum_eigenvalue
        )
        history["maximum_internal_gram_eigenvalue"].append(
            maximum_eigenvalue
        )
        history["gram_condition_number"].append(
            maximum_eigenvalue / minimum_eigenvalue
        )
        history["connectivity_offdiagonal_l2"].append(
            connectivity_norm
        )

        if accepted_iterations in save_steps_set:
            checkpoint_gram = (
                _num.cp.asnumpy(reduced_gram)
                if use_gpu
                else np.asarray(reduced_gram)
            )
            connectivity_at_steps[accepted_iterations] = (
                _num._connectivity_from_reduced_gram(
                    checkpoint_gram, basis
                )
            )

        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": "newton",
                    "iteration": accepted_iterations,
                    "total": max_iterations,
                    "objective": float(candidate_objective),
                    "loss": relative_loss,
                    "entropy": entropy,
                    "relative_gradient_norm": candidate_relative_kkt,
                    "relative_eliminated_kkt_residual": (
                        candidate_relative_kkt
                    ),
                    "step_size": float(step_size),
                    "cg_iterations": int(cg_info["iterations"]),
                    "cg_relative_residual": float(
                        cg_info["relative_residual"]
                    ),
                    "cg_target": float(cg_target),
                    "cg_converged": bool(cg_info["converged"]),
                    "readiness_probe": bool(probe_iteration),
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV-Newton-whitened",
                    "general_method": "covariance_optimization",
                }
            )

        if candidate_gradient_norm <= (
            absolute_tolerance
            + relative_tolerance * candidate_gradient_scale
        ):
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break

    for key, values in history.items():
        if key in {
            "iteration",
            "line_search_backtracks",
            "cg_iterations",
            "cg_hessian_actions",
        }:
            history[key] = np.asarray(values, dtype=np.int64)
        elif key in {
            "cg_converged",
            "readiness_probe",
            "accepted_truncated_cg",
        }:
            history[key] = np.asarray(values, dtype=bool)
        elif key == "cg_termination_reason":
            history[key] = np.asarray(values, dtype=object)
        else:
            history[key] = np.asarray(values, dtype=np.float64)

    if use_gpu:
        reduced_gram = _num.cp.asnumpy(reduced_gram)
    full_gram = _symmetrize(basis @ reduced_gram @ basis.T)
    fitted_squared_distances = _num._squared_distances_from_gram(full_gram)
    connectivity = _num._connectivity_from_reduced_gram(
        reduced_gram, basis
    )

    from .covariance_pdhg import _independent_eliminated_kkt_residuals

    independent = _independent_eliminated_kkt_residuals(
        full_gram,
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
    if independent_converged:
        converged = True
        status = "optimality_tolerance"
        message = "independent eliminated KKT tolerance reached"
    elif converged:
        converged = False
        status = "independent_kkt_failed"
        message = (
            "internal Newton tolerance was reached, but the returned Gram "
            "failed the independent KKT certificate"
        )

    fitted_pairs = fitted_squared_distances[pair_i, pair_j]
    pair_stationarity = (
        fitted_pairs
        - target_pairs
        - 0.5 * pair_variance * connectivity[pair_i, pair_j]
    )
    pair_stationarity_norm = float(np.linalg.norm(pair_stationarity))
    pair_stationarity_scale = max(
        1.0,
        float(np.linalg.norm(fitted_pairs - target_pairs)),
        float(
            np.linalg.norm(
                0.5 * pair_variance * connectivity[pair_i, pair_j]
            )
        ),
    )

    final_state = _num._gaussian_covariance_objective_gradient(
        np.asarray(reduced_gram),
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
        array_module=np,
    )
    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": int(accepted_iterations),
        "objective": float(final_state[0]),
        "relative_gradient_norm": float(independent_relative),
        "relative_eliminated_kkt_residual": float(independent_relative),
        "stationarity_residual_norm": float(independent_norm),
        "stationarity_residual_scale": float(independent_scale),
        "relative_stationarity_residual": float(independent_relative),
        "maximum_absolute_stationarity_residual": float(
            np.max(np.abs(independent_residual))
        ),
        "pair_stationarity_residual_norm": pair_stationarity_norm,
        "relative_pair_stationarity_residual": (
            pair_stationarity_norm / pair_stationarity_scale
        ),
        "maximum_absolute_pair_stationarity_residual": float(
            np.max(np.abs(pair_stationarity))
        ),
        "independent_kkt_converged": bool(independent_converged),
        "independent_kkt_recomputed_from_returned_gram": True,
        "observed_pair_count": len(target_pairs),
        "noise_model": noise_model,
        "noise_parameter": noise_parameter,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization_info,
        "algorithm": "entropy_whitened_newton_pcg",
        "coordinate_parameterization": "centered_internal_svec",
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "cupy_version": _num.cp.__version__ if use_gpu else None,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "readiness_probe": readiness_info,
        "newton": {
            "entropy_whitened": True,
            "svec": True,
            "exact_whitened_hessian_diagonal": True,
            "adaptive_forcing": True,
            "cg_forcing_min": float(cg_forcing_min),
            "cg_forcing_max": float(cg_forcing_max),
            "cg_forcing_coefficient": float(cg_forcing_coefficient),
            "cg_max_iterations": int(cg_max_iterations),
            "kkt_safeguard": True,
            "maximum_kkt_growth": float(maximum_kkt_growth),
            "truncated_kkt_reduction": float(
                truncated_kkt_reduction
            ),
            "direct_gram_handoff": bool(initial_gram is not None),
            "scalar_calibration_applied": bool(
                initialization_info["scalar_calibration_applied"]
            ),
        },
        "preconditioner_setup_seconds": float(
            total_preconditioner_seconds
        ),
        "cg_wall_seconds": float(total_cg_seconds),
        "wall_seconds": time.perf_counter() - start_time,
        "objective_definition": (
            "-1.5*logdet(B_internal) + 0.5*sum_unique_pairs("
            "(D_fit-D_obs)^2/noise_variance)"
        ),
        "stationarity_definition": (
            "D_adjoint((D(B)-D_obs)/noise_variance)-1.5*B_pseudoinverse"
        ),
        "pair_stationarity_definition": (
            "D_fit_ij-D_obs_ij-noise_variance_ij*A_ij/2"
        ),
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_newton_whitened stopped without "
            f"satisfying the KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, full_gram, connectivity, info
