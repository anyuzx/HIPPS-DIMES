"""Globalized covariance-cone solver for noisy HIPPS-DIMES constraints.

This module installs an improved replacement for
``numerics.fit_gaussian_noise_covariance``. The final objective is unchanged;
only the path used to reach its minimizer is modified. The implementation adds
variance continuation, monotone log-determinant proximal warm-up steps,
inexact Newton forcing, whitened Newton-CG coordinates, and a feasible scalar
line search. It intentionally reuses the validated objective, Hessian, and
initialization kernels in :mod:`hipps_dimes.numerics`.
"""

from __future__ import annotations

import time
import warnings
from collections.abc import Iterable

import numpy as np

from . import numerics as _num


_DEFAULT_CONTINUATION_FACTORS = (0.1, 1.0)


def _scalar(value) -> float:
    """Convert a NumPy/CuPy scalar to a Python float."""
    return float(value.item()) if hasattr(value, "item") else float(value)


def _symmetrize(matrix):
    return 0.5 * (matrix + matrix.T)


def _frobenius_norm(matrix, array_module) -> float:
    return _scalar(array_module.sqrt(array_module.sum(matrix * matrix)))


def _validate_continuation_factors(
    factors: Iterable[float] | None,
) -> tuple[float, ...]:
    """Return an increasing continuation schedule ending exactly at one."""
    if factors is None:
        factors = _DEFAULT_CONTINUATION_FACTORS
    try:
        schedule = tuple(float(value) for value in factors)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "continuation_factors must be an iterable of numbers"
        ) from error
    if not schedule:
        schedule = (1.0,)
    if any(
        not np.isfinite(value) or value <= 0.0 or value > 1.0
        for value in schedule
    ):
        raise ValueError("continuation_factors must lie in (0, 1]")
    schedule = tuple(sorted(set(schedule)))
    if schedule[-1] != 1.0:
        schedule = schedule + (1.0,)
    return schedule


def _inexact_newton_forcing_tolerance(
    relative_gradient_norm: float,
    minimum_tolerance: float,
    maximum_tolerance: float,
) -> float:
    """Eisenstat-Walker-like forcing sequence used by the inner CG solve."""
    candidate = 0.1 * np.sqrt(max(float(relative_gradient_norm), 0.0))
    return float(
        min(maximum_tolerance, max(minimum_tolerance, candidate))
    )


def _logdet_proximal_step(matrix, step_size: float, array_module):
    """Apply ``prox_{step*(-3/2 logdet)}`` to one symmetric matrix."""
    matrix = _symmetrize(matrix)
    if array_module is np:
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    else:
        with _num.cupyx.errstate(linalg="raise"):
            eigenvalues, eigenvectors = array_module.linalg.eigh(matrix)
    updated = 0.5 * (
        eigenvalues
        + array_module.sqrt(
            array_module.square(eigenvalues) + 6.0 * step_size
        )
    )
    result = (eigenvectors * updated) @ eigenvectors.T
    return _symmetrize(result)


def _stage_objective_gradient(
    reduced_gram,
    basis,
    target_matrix,
    base_inverse_variance_matrix,
    pair_i,
    pair_j,
    continuation_factor,
    *,
    array_module,
):
    """Evaluate one continuation-stage objective and gradient."""
    return _num._gaussian_covariance_objective_gradient(
        reduced_gram,
        basis,
        target_matrix,
        continuation_factor * base_inverse_variance_matrix,
        pair_i,
        pair_j,
        array_module=array_module,
    )


def _stage_data_hessian_action(
    direction,
    basis,
    base_inverse_variance_matrix,
    continuation_factor,
    *,
    array_module,
):
    return continuation_factor * _num._gaussian_covariance_data_hessian_action(
        direction,
        basis,
        base_inverse_variance_matrix,
        array_module=array_module,
    )


def _relative_gradient_from_state(
    state, array_module
) -> tuple[float, float, float]:
    """Return gradient norm, scale, and relative norm from an objective state."""
    gradient = state[1]
    components = state[4]
    entropy_gradient = components[2]
    data_gradient = components[3]
    gradient_norm = _frobenius_norm(gradient, array_module)
    gradient_scale = max(
        1.0,
        _frobenius_norm(entropy_gradient, array_module),
        _frobenius_norm(data_gradient, array_module),
    )
    return gradient_norm, gradient_scale, gradient_norm / gradient_scale


def _proximal_warmup(
    reduced_gram,
    state,
    basis,
    target_matrix,
    base_inverse_variance_matrix,
    pair_i,
    pair_j,
    continuation_factor,
    max_steps,
    switch_relative_gradient,
    initial_lipschitz,
    *,
    array_module,
):
    """Run monotone proximal-gradient steps for the exact stage objective."""
    records = []
    lipschitz = max(float(initial_lipschitz), np.finfo(np.float64).tiny)
    for _ in range(int(max_steps)):
        _, _, relative_gradient = _relative_gradient_from_state(
            state, array_module
        )
        if relative_gradient <= switch_relative_gradient:
            break
        objective = float(state[0])
        data_objective = float(state[4][1])
        data_gradient = state[4][3]
        accepted = False
        backtracks = 0
        for backtracks in range(60):
            step_size = 1.0 / lipschitz
            trial_center = reduced_gram - step_size * data_gradient
            candidate = _logdet_proximal_step(
                trial_center, step_size, array_module
            )
            candidate_state = _stage_objective_gradient(
                candidate,
                basis,
                target_matrix,
                base_inverse_variance_matrix,
                pair_i,
                pair_j,
                continuation_factor,
                array_module=array_module,
            )
            candidate_data_objective = float(candidate_state[4][1])
            delta = candidate - reduced_gram
            delta_norm_squared = _scalar(
                array_module.sum(delta * delta)
            )
            majorizer = (
                data_objective
                + _scalar(array_module.sum(data_gradient * delta))
                + 0.5 * lipschitz * delta_norm_squared
            )
            scale = max(
                1.0, abs(objective), abs(float(candidate_state[0]))
            )
            slack = 1e-12 * scale
            if (
                candidate_data_objective <= majorizer + slack
                and float(candidate_state[0]) <= objective + slack
            ):
                accepted = True
                break
            lipschitz *= 2.0
        if not accepted:
            break
        relative_step = _frobenius_norm(delta, array_module) / max(
            _frobenius_norm(reduced_gram, array_module), 1.0
        )
        reduced_gram = candidate
        state = candidate_state
        records.append(
            {
                "gram": reduced_gram.copy(),
                "state": state,
                "relative_step": relative_step,
                "step_size": step_size,
                "line_search_backtracks": int(backtracks),
                "full_step_accepted": True,
                "line_search_method": "logdet_proximal_backtracking",
            }
        )
        lipschitz = max(0.8 * lipschitz, np.finfo(np.float64).tiny)
    return reduced_gram, state, records, lipschitz


def _feasible_exact_scalar_step(
    reduced_gram,
    direction,
    whitened_direction,
    data_gradient,
    data_hessian_direction,
    *,
    array_module,
    maximum_step=1.0,
):
    """Minimize the objective exactly along a feasible Newton direction."""
    del reduced_gram
    if array_module is np:
        mode_eigenvalues = np.linalg.eigvalsh(
            _symmetrize(whitened_direction)
        )
    else:
        with _num.cupyx.errstate(linalg="raise"):
            mode_eigenvalues = array_module.linalg.eigvalsh(
                _symmetrize(whitened_direction)
            )
    mode_eigenvalues = np.asarray(
        mode_eigenvalues
        if array_module is np
        else _num.cp.asnumpy(mode_eigenvalues),
        dtype=np.float64,
    )
    minimum_mode = float(np.min(mode_eigenvalues))
    if minimum_mode < 0.0:
        feasible_limit = -1.0 / minimum_mode
        upper = min(float(maximum_step), np.nextafter(feasible_limit, 0.0))
    else:
        feasible_limit = np.inf
        upper = float(maximum_step)
    if not np.isfinite(upper) or upper <= 0.0:
        return 0.0, feasible_limit

    q1 = _scalar(array_module.sum(data_gradient * direction))
    q2 = _scalar(array_module.sum(direction * data_hessian_direction))
    q2 = max(q2, 0.0)

    def derivative(alpha: float) -> float:
        denominator = 1.0 + alpha * mode_eigenvalues
        return float(
            q1
            + alpha * q2
            - 1.5 * np.sum(mode_eigenvalues / denominator)
        )

    derivative_zero = derivative(0.0)
    if not np.isfinite(derivative_zero) or derivative_zero >= 0.0:
        return 0.0, feasible_limit
    derivative_upper = derivative(upper)
    if np.isfinite(derivative_upper) and derivative_upper <= 0.0:
        return upper, feasible_limit

    lower = 0.0
    alpha = min(1.0, 0.5 * upper)
    for _ in range(80):
        value = derivative(alpha)
        if abs(value) <= 1e-12 * max(1.0, abs(derivative_zero)):
            break
        if value > 0.0:
            upper = alpha
        else:
            lower = alpha
        denominator = 1.0 + alpha * mode_eigenvalues
        second = q2 + 1.5 * np.sum(
            np.square(mode_eigenvalues) / np.square(denominator)
        )
        proposal = alpha - value / second if second > 0.0 else np.nan
        if (
            not np.isfinite(proposal)
            or proposal <= lower
            or proposal >= upper
        ):
            proposal = 0.5 * (lower + upper)
        alpha = proposal
    return float(alpha), feasible_limit


def _armijo_fallback(
    reduced_gram,
    direction,
    objective,
    directional_derivative,
    evaluate_state,
    initial_step,
):
    """Fallback SPD/Armijo line search used if exact minimization fails."""
    step_size = float(initial_step)
    for backtracks in range(60):
        trial = _symmetrize(reduced_gram + step_size * direction)
        state = evaluate_state(trial)
        if np.isfinite(state[0]) and float(state[0]) <= (
            objective + 1e-4 * step_size * directional_derivative
        ):
            return trial, state, step_size, backtracks
        step_size *= 0.5
    return None, None, 0.0, 60


def _history_container():
    return {
        "iteration": [],
        "objective": [],
        "stage_objective": [],
        "negative_entropy_objective": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "gradient_norm": [],
        "relative_gradient_norm": [],
        "newton_decrement": [],
        "relative_step": [],
        "step_size": [],
        "cg_iterations": [],
        "cg_relative_residual": [],
        "cg_forcing_tolerance": [],
        "cg_converged": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "connectivity_offdiagonal_l2": [],
        "continuation_factor": [],
        "phase": [],
        "line_search_backtracks": [],
        "full_step_accepted": [],
        "line_search_method": [],
        "feasible_step_bound": [],
        "fallback_direction": [],
    }


def fit_gaussian_noise_covariance(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initialization="rouse",
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=100,
    relative_tolerance=1e-8,
    absolute_tolerance=1e-10,
    cg_relative_tolerance=1e-4,
    cg_max_iterations=None,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
    continuation_factors=_DEFAULT_CONTINUATION_FACTORS,
    continuation_intermediate_tolerance=1e-3,
    continuation_intermediate_newton_iterations=5,
    proximal_warmup_iterations=8,
    proximal_switch_relative_gradient=0.25,
    continuation_activation_relative_gradient=1.0,
    cg_forcing_max=0.5,
    exact_line_search=True,
    line_search_max_step=1.0,
    use_whitened_newton=True,
):
    """Fit Gaussian-noisy squared distances with a globalized COV solver.

    The final objective is identical to the legacy COV objective. Intermediate
    continuation stages multiply the inverse variances by a factor in ``(0, 1]``
    and are used only to obtain a robust path to the final ``factor=1`` solution.
    ``max_iterations`` counts all accepted proximal and Newton updates across all
    stages.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "COV GPU optimization was requested, but CuPy and an accessible "
            "CUDA GPU are not available"
        )
    (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        target_pairs,
        pair_variance,
        noise_model,
        noise_parameter,
    ) = _num._validate_gaussian_covariance_inputs(
        squared_distances, noise_variance, relative_noise_std
    )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(cg_relative_tolerance) or not (
        0.0 < cg_relative_tolerance < 1.0
    ):
        raise ValueError(
            "cg_relative_tolerance must lie strictly between zero and one"
        )
    if not np.isfinite(cg_forcing_max) or not (
        cg_relative_tolerance <= cg_forcing_max < 1.0
    ):
        raise ValueError(
            "cg_forcing_max must be at least cg_relative_tolerance and below one"
        )
    if (
        not np.isfinite(initial_gram_floor_relative)
        or initial_gram_floor_relative <= 0.0
    ):
        raise ValueError(
            "initial_gram_floor_relative must be positive and finite"
        )
    if (
        not isinstance(proximal_warmup_iterations, (int, np.integer))
        or proximal_warmup_iterations < 0
    ):
        raise ValueError(
            "proximal_warmup_iterations must be a nonnegative integer"
        )
    if (
        not np.isfinite(proximal_switch_relative_gradient)
        or proximal_switch_relative_gradient <= 0.0
    ):
        raise ValueError(
            "proximal_switch_relative_gradient must be positive"
        )
    if (
        not np.isfinite(continuation_activation_relative_gradient)
        or continuation_activation_relative_gradient <= 0.0
    ):
        raise ValueError(
            "continuation_activation_relative_gradient must be positive"
        )
    if not np.isfinite(line_search_max_step) or line_search_max_step <= 0.0:
        raise ValueError("line_search_max_step must be positive")
    schedule = _validate_continuation_factors(continuation_factors)
    if max_iterations < 10:
        schedule = (1.0,)

    n = observed.shape[0]
    n_modes = n - 1
    if cg_max_iterations is None:
        cg_max_iterations = 2 * n_modes
    if (
        not isinstance(cg_max_iterations, (int, np.integer))
        or cg_max_iterations <= 0
    ):
        raise ValueError("cg_max_iterations must be a positive integer")
    save_steps_set = set()
    if save_steps is not None:
        for step in save_steps:
            if not isinstance(step, (int, np.integer)) or not (
                1 <= step <= max_iterations
            ):
                raise ValueError(
                    "save_steps must contain integers between 1 and max_iterations"
                )
            save_steps_set.add(int(step))

    basis = _num._centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = (
        _num._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    reduced_gram, initialization_info = (
        _num._initialize_gaussian_reduced_gram(
            observed,
            pair_mask,
            pair_i,
            pair_j,
            inverse_variance,
            basis,
            initialization,
            initial_connectivity,
            initial_gram_floor_relative,
            use_gpu=use_gpu,
            progress_callback=progress_callback,
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
        reduced_gram = xp.asarray(reduced_gram, dtype=xp.float64)
        cg_solver = _num._preconditioned_conjugate_gradient_gpu
    else:
        xp = np
        solver_basis = basis
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target_pairs = target_pairs
        solver_inverse_variance = inverse_variance
        solver_target_matrix = target_matrix
        solver_inverse_variance_matrix = inverse_variance_matrix
        cg_solver = _num._preconditioned_conjugate_gradient

    reduced_gram, scalar_calibration = (
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
    initialization_info["scalar_calibration"] = scalar_calibration
    initialization_info["wall_seconds"] = (
        initialization_info.get("wall_seconds", 0.0)
        + scalar_calibration["wall_seconds"]
    )
    if initialization_info["kind"] == "rouse":
        initialization_info["effective_spring_constant"] = (
            initialization_info["spring_constant"]
            * scalar_calibration["connectivity_scale_factor"]
        )

    base_data_hessian_diagonal, preconditioner_setup_seconds = (
        _num._gaussian_covariance_data_preconditioner_diagonal(
            solver_basis,
            solver_pair_i,
            solver_pair_j,
            solver_inverse_variance,
            progress_callback=progress_callback,
            array_module=xp,
            use_gpu=use_gpu,
        )
    )
    base_lipschitz = 4.0 * _scalar(
        xp.max(xp.sum(solver_inverse_variance_matrix, axis=1))
    )
    base_lipschitz = max(base_lipschitz, np.finfo(np.float64).tiny)

    initial_final_state = _stage_objective_gradient(
        reduced_gram,
        solver_basis,
        solver_target_matrix,
        solver_inverse_variance_matrix,
        solver_pair_i,
        solver_pair_j,
        1.0,
        array_module=xp,
    )
    _, _, initial_final_relative_gradient = _relative_gradient_from_state(
        initial_final_state, xp
    )
    if (
        initial_final_relative_gradient
        <= continuation_activation_relative_gradient
    ):
        schedule = (1.0,)

    history = _history_container()
    connectivity_at_steps = {}
    continuation_info = []
    global_iteration = 0
    converged = False
    status = "max_iterations"
    message = "maximum number of iterations reached"
    current_state = None
    final_stage_reached = False
    start_time = time.perf_counter()

    def record_step(
        state,
        phase,
        factor,
        relative_step,
        step_size,
        cg_iterations=0,
        cg_residual=0.0,
        forcing_tolerance=np.nan,
        cg_converged=True,
        directional_derivative=0.0,
        backtracks=0,
        full_step=True,
        line_search_method="none",
        feasible_step_bound=np.inf,
        fallback_direction=False,
    ):
        nonlocal global_iteration
        global_iteration += 1
        objective, gradient, fitted_pairs, inverse_gram, components = state
        (
            negative_entropy,
            data_objective,
            entropy_gradient,
            data_gradient,
        ) = components
        gradient_norm, _, relative_gradient = _relative_gradient_from_state(
            state, xp
        )
        residual = fitted_pairs - solver_target_pairs
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
        reduced_connectivity = 3.0 * inverse_gram
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
        cg_relative_residual = float(cg_residual) / max(
            gradient_norm, np.finfo(np.float64).tiny
        )
        final_state = _stage_objective_gradient(
            reduced_gram,
            solver_basis,
            solver_target_matrix,
            solver_inverse_variance_matrix,
            solver_pair_i,
            solver_pair_j,
            1.0,
            array_module=xp,
        )
        history["iteration"].append(global_iteration)
        history["objective"].append(float(final_state[0]))
        history["stage_objective"].append(float(objective))
        history["negative_entropy_objective"].append(float(negative_entropy))
        history["data_objective"].append(float(data_objective))
        history["loss"].append(relative_loss)
        history["entropy"].append(entropy)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["gradient_norm"].append(gradient_norm)
        history["relative_gradient_norm"].append(relative_gradient)
        history["newton_decrement"].append(
            float(np.sqrt(max(-directional_derivative, 0.0)))
        )
        history["relative_step"].append(float(relative_step))
        history["step_size"].append(float(step_size))
        history["cg_iterations"].append(int(cg_iterations))
        history["cg_relative_residual"].append(cg_relative_residual)
        history["cg_forcing_tolerance"].append(float(forcing_tolerance))
        history["cg_converged"].append(bool(cg_converged))
        history["minimum_internal_gram_eigenvalue"].append(
            _scalar(gram_eigenvalues[0])
        )
        history["maximum_internal_gram_eigenvalue"].append(
            _scalar(gram_eigenvalues[-1])
        )
        history["connectivity_offdiagonal_l2"].append(connectivity_norm)
        history["continuation_factor"].append(float(factor))
        history["phase"].append(str(phase))
        history["line_search_backtracks"].append(int(backtracks))
        history["full_step_accepted"].append(bool(full_step))
        history["line_search_method"].append(str(line_search_method))
        history["feasible_step_bound"].append(float(feasible_step_bound))
        history["fallback_direction"].append(bool(fallback_direction))
        if global_iteration in save_steps_set:
            checkpoint = _num.cp.asnumpy(reduced_gram) if use_gpu else reduced_gram
            connectivity_at_steps[global_iteration] = (
                _num._connectivity_from_reduced_gram(checkpoint, basis)
            )
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": str(phase),
                    "continuation_factor": float(factor),
                    "iteration": global_iteration,
                    "total": max_iterations,
                    "objective": float(final_state[0]),
                    "stage_objective": float(objective),
                    "loss": relative_loss,
                    "entropy": entropy,
                    "relative_gradient_norm": relative_gradient,
                    "step_size": float(step_size),
                    "line_search_backtracks": int(backtracks),
                    "full_step_accepted": bool(full_step),
                    "cg_iterations": int(cg_iterations),
                    "cg_converged": bool(cg_converged),
                    "cg_forcing_tolerance": float(forcing_tolerance),
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV",
                    "general_method": "covariance_optimization",
                }
            )

    minimum_final_budget = min(10, max_iterations)
    for factor in schedule:
        if global_iteration >= max_iterations:
            break
        final_stage = factor == 1.0
        final_stage_reached = final_stage_reached or final_stage
        stage_start_iteration = global_iteration
        stage_inverse_variance = factor * solver_inverse_variance
        reduced_gram, stage_calibration = (
            _num._calibrate_gaussian_covariance_initial_scale(
                reduced_gram,
                solver_basis,
                solver_pair_i,
                solver_pair_j,
                solver_target_pairs,
                stage_inverse_variance,
                array_module=xp,
            )
        )
        current_state = _stage_objective_gradient(
            reduced_gram,
            solver_basis,
            solver_target_matrix,
            solver_inverse_variance_matrix,
            solver_pair_i,
            solver_pair_j,
            factor,
            array_module=xp,
        )
        stage_tolerance = (
            relative_tolerance
            if final_stage
            else max(relative_tolerance, continuation_intermediate_tolerance)
        )
        if final_stage:
            stage_budget = max_iterations - global_iteration
        else:
            available = max_iterations - global_iteration - minimum_final_budget
            stage_budget = max(
                0,
                min(
                    int(proximal_warmup_iterations)
                    + int(continuation_intermediate_newton_iterations),
                    available,
                ),
            )
        stage_record = {
            "factor": float(factor),
            "scalar_calibration": stage_calibration,
            "start_iteration": global_iteration + 1,
            "warmup_steps": 0,
            "newton_steps": 0,
            "converged": False,
        }
        if stage_budget <= 0 and not final_stage:
            continuation_info.append(stage_record)
            continue

        warmup_limit = min(int(proximal_warmup_iterations), stage_budget)
        reduced_gram, current_state, warmup_records, _ = _proximal_warmup(
            reduced_gram,
            current_state,
            solver_basis,
            solver_target_matrix,
            solver_inverse_variance_matrix,
            solver_pair_i,
            solver_pair_j,
            factor,
            warmup_limit,
            proximal_switch_relative_gradient,
            factor * base_lipschitz,
            array_module=xp,
        )
        for warmup_record in warmup_records:
            if global_iteration >= max_iterations:
                break
            reduced_gram = warmup_record["gram"]
            current_state = warmup_record["state"]
            record_step(
                current_state,
                "proximal_warmup",
                factor,
                warmup_record["relative_step"],
                warmup_record["step_size"],
                backtracks=warmup_record["line_search_backtracks"],
                full_step=warmup_record["full_step_accepted"],
                line_search_method=warmup_record["line_search_method"],
            )
            stage_record["warmup_steps"] += 1
        stage_budget -= stage_record["warmup_steps"]

        while stage_budget > 0 and global_iteration < max_iterations:
            objective, gradient, _, inverse_gram, components = current_state
            gradient_norm, gradient_scale, relative_gradient = (
                _relative_gradient_from_state(current_state, xp)
            )
            if gradient_norm <= (
                absolute_tolerance + stage_tolerance * gradient_scale
            ):
                stage_record["converged"] = True
                break
            forcing_tolerance = _inexact_newton_forcing_tolerance(
                relative_gradient,
                cg_relative_tolerance,
                cg_forcing_max,
            )
            fallback_direction = False
            if use_whitened_newton:
                if use_gpu:
                    with _num.cupyx.errstate(linalg="raise"):
                        cholesky_factor = xp.linalg.cholesky(reduced_gram)
                else:
                    cholesky_factor = np.linalg.cholesky(reduced_gram)
                transformed_gradient = (
                    cholesky_factor.T @ gradient @ cholesky_factor
                )

                def whitened_hessian_operator(matrix_direction):
                    matrix_direction = _symmetrize(matrix_direction)
                    gram_direction = (
                        cholesky_factor
                        @ matrix_direction
                        @ cholesky_factor.T
                    )
                    data_action = _stage_data_hessian_action(
                        gram_direction,
                        solver_basis,
                        solver_inverse_variance_matrix,
                        factor,
                        array_module=xp,
                    )
                    result = (
                        cholesky_factor.T
                        @ data_action
                        @ cholesky_factor
                        + 1.5 * matrix_direction
                    )
                    return _symmetrize(result)

                diagonal_scale = xp.square(xp.diag(cholesky_factor))
                approximate_data_diagonal = (
                    factor
                    * base_data_hessian_diagonal
                    * xp.outer(diagonal_scale, diagonal_scale)
                )
                preconditioner_diagonal = xp.maximum(
                    1.5 + approximate_data_diagonal,
                    np.finfo(np.float64).tiny,
                )
                (
                    whitened_direction,
                    cg_iterations,
                    cg_residual,
                    cg_converged,
                ) = cg_solver(
                    whitened_hessian_operator,
                    -transformed_gradient,
                    preconditioner_diagonal,
                    forcing_tolerance,
                    int(cg_max_iterations),
                )
                whitened_direction = _symmetrize(whitened_direction)
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
                    fallback_direction = True
                    whitened_direction = (
                        -transformed_gradient / preconditioner_diagonal
                    )
                    whitened_direction = _symmetrize(whitened_direction)
                    direction = _symmetrize(
                        cholesky_factor
                        @ whitened_direction
                        @ cholesky_factor.T
                    )
                    directional_derivative = _scalar(
                        xp.sum(gradient * direction)
                    )
                    cg_converged = False
            else:

                def hessian_operator(matrix_direction):
                    matrix_direction = _symmetrize(matrix_direction)
                    data_action = _stage_data_hessian_action(
                        matrix_direction,
                        solver_basis,
                        solver_inverse_variance_matrix,
                        factor,
                        array_module=xp,
                    )
                    entropy_action = (
                        1.5
                        * inverse_gram
                        @ matrix_direction
                        @ inverse_gram
                    )
                    return _symmetrize(data_action + entropy_action)

                inverse_diagonal = xp.diag(inverse_gram)
                preconditioner_diagonal = xp.maximum(
                    factor * base_data_hessian_diagonal
                    + 1.5 * xp.outer(inverse_diagonal, inverse_diagonal),
                    np.finfo(np.float64).tiny,
                )
                direction, cg_iterations, cg_residual, cg_converged = cg_solver(
                    hessian_operator,
                    -gradient,
                    preconditioner_diagonal,
                    forcing_tolerance,
                    int(cg_max_iterations),
                )
                direction = _symmetrize(direction)
                directional_derivative = _scalar(xp.sum(gradient * direction))
                if (
                    not np.isfinite(directional_derivative)
                    or directional_derivative >= 0.0
                ):
                    fallback_direction = True
                    direction = _symmetrize(-gradient / preconditioner_diagonal)
                    directional_derivative = _scalar(
                        xp.sum(gradient * direction)
                    )
                    cg_converged = False
                if use_gpu:
                    identity = xp.eye(n_modes, dtype=xp.float64)
                    inverse_cholesky = (
                        _num.cupyx_scipy_linalg.solve_triangular(
                            xp.linalg.cholesky(reduced_gram),
                            identity,
                            lower=True,
                            check_finite=False,
                        )
                    )
                    whitened_direction = (
                        inverse_cholesky
                        @ direction
                        @ inverse_cholesky.T
                    )
                else:
                    cholesky_factor = np.linalg.cholesky(reduced_gram)
                    whitened_direction = np.linalg.solve(
                        cholesky_factor,
                        np.linalg.solve(cholesky_factor, direction).T,
                    ).T
                    whitened_direction = _symmetrize(whitened_direction)

            def evaluate_stage(trial):
                return _stage_objective_gradient(
                    trial,
                    solver_basis,
                    solver_target_matrix,
                    solver_inverse_variance_matrix,
                    solver_pair_i,
                    solver_pair_j,
                    factor,
                    array_module=xp,
                )

            backtracks = 0
            feasible_bound = np.inf
            line_search_method = "exact_feasible_scalar"
            if exact_line_search and directional_derivative < 0.0:
                data_hessian_direction = _stage_data_hessian_action(
                    direction,
                    solver_basis,
                    solver_inverse_variance_matrix,
                    factor,
                    array_module=xp,
                )
                step_size, feasible_bound = _feasible_exact_scalar_step(
                    reduced_gram,
                    direction,
                    whitened_direction,
                    components[3],
                    data_hessian_direction,
                    array_module=xp,
                    maximum_step=line_search_max_step,
                )
                candidate = None
                candidate_state = None
                if step_size > 0.0:
                    candidate = _symmetrize(
                        reduced_gram + step_size * direction
                    )
                    candidate_state = evaluate_stage(candidate)
                    scale = max(
                        1.0,
                        abs(float(objective)),
                        abs(float(candidate_state[0])),
                    )
                    if (
                        not np.isfinite(candidate_state[0])
                        or float(candidate_state[0])
                        > float(objective) + 1e-11 * scale
                    ):
                        candidate = None
                        candidate_state = None
                if candidate is None:
                    line_search_method = "armijo_fallback"
                    initial_step = min(
                        1.0,
                        0.99 * feasible_bound
                        if np.isfinite(feasible_bound)
                        else 1.0,
                    )
                    (
                        candidate,
                        candidate_state,
                        step_size,
                        backtracks,
                    ) = _armijo_fallback(
                        reduced_gram,
                        direction,
                        float(objective),
                        directional_derivative,
                        evaluate_stage,
                        initial_step,
                    )
            else:
                line_search_method = "armijo"
                candidate, candidate_state, step_size, backtracks = (
                    _armijo_fallback(
                        reduced_gram,
                        direction,
                        float(objective),
                        directional_derivative,
                        evaluate_stage,
                        1.0,
                    )
                )
            if candidate is None:
                status = "line_search_failed"
                message = "covariance-cone line search failed"
                stage_budget = 0
                break
            delta = candidate - reduced_gram
            relative_step = _frobenius_norm(delta, xp) / max(
                _frobenius_norm(reduced_gram, xp), 1.0
            )
            reduced_gram = candidate
            current_state = candidate_state
            record_step(
                current_state,
                "newton",
                factor,
                relative_step,
                step_size,
                cg_iterations=cg_iterations,
                cg_residual=cg_residual,
                forcing_tolerance=forcing_tolerance,
                cg_converged=cg_converged,
                directional_derivative=directional_derivative,
                backtracks=backtracks,
                full_step=abs(step_size - 1.0) <= 1e-12,
                line_search_method=line_search_method,
                feasible_step_bound=feasible_bound,
                fallback_direction=fallback_direction,
            )
            stage_record["newton_steps"] += 1
            stage_budget -= 1

        if current_state is not None:
            gradient_norm, gradient_scale, relative_gradient = (
                _relative_gradient_from_state(current_state, xp)
            )
            stage_record["final_relative_gradient_norm"] = relative_gradient
            if gradient_norm <= (
                absolute_tolerance + stage_tolerance * gradient_scale
            ):
                stage_record["converged"] = True
        stage_record["end_iteration"] = global_iteration
        stage_record["accepted_steps"] = (
            global_iteration - stage_start_iteration
        )
        continuation_info.append(stage_record)
        if final_stage and stage_record["converged"]:
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break
        if status == "line_search_failed":
            break

    current_state = _stage_objective_gradient(
        reduced_gram,
        solver_basis,
        solver_target_matrix,
        solver_inverse_variance_matrix,
        solver_pair_i,
        solver_pair_j,
        1.0,
        array_module=xp,
    )
    final_gradient_norm, final_gradient_scale, final_relative_gradient_norm = (
        _relative_gradient_from_state(current_state, xp)
    )
    if final_stage_reached and final_gradient_norm <= (
        absolute_tolerance + relative_tolerance * final_gradient_scale
    ):
        converged = True
        status = "optimality_tolerance"
        message = "first-order optimality tolerance reached"

    for key, values in history.items():
        if key in {"iteration", "cg_iterations", "line_search_backtracks"}:
            history[key] = np.asarray(values, dtype=np.int64)
        elif key in {"cg_converged", "full_step_accepted", "fallback_direction"}:
            history[key] = np.asarray(values, dtype=bool)
        elif key in {"phase", "line_search_method"}:
            history[key] = np.asarray(values, dtype=object)
        else:
            history[key] = np.asarray(values, dtype=np.float64)

    if use_gpu:
        reduced_gram = _num.cp.asnumpy(reduced_gram)
    gram = basis @ reduced_gram @ basis.T
    gram = _symmetrize(gram)
    fitted_squared_distances = _num._squared_distances_from_gram(gram)
    connectivity = _num._connectivity_from_reduced_gram(reduced_gram, basis)
    fitted_pairs = fitted_squared_distances[pair_i, pair_j]
    stationarity = (
        fitted_pairs
        - target_pairs
        - 0.5 * pair_variance * connectivity[pair_i, pair_j]
    )
    stationarity_norm = float(np.linalg.norm(stationarity))
    stationarity_scale = max(
        1.0,
        float(np.linalg.norm(fitted_pairs - target_pairs)),
        float(
            np.linalg.norm(
                0.5 * pair_variance * connectivity[pair_i, pair_j]
            )
        ),
    )

    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": int(global_iteration),
        "objective": float(current_state[0]),
        "gradient_norm": final_gradient_norm,
        "relative_gradient_norm": final_relative_gradient_norm,
        "stationarity_residual_norm": stationarity_norm,
        "relative_stationarity_residual": stationarity_norm / stationarity_scale,
        "maximum_absolute_stationarity_residual": float(
            np.max(np.abs(stationarity))
        ),
        "observed_pair_count": len(target_pairs),
        "noise_model": noise_model,
        "noise_parameter": noise_parameter,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization_info,
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "cupy_version": _num.cp.__version__ if use_gpu else None,
        "preconditioner_pair_block_size": (
            _num._COVARIANCE_PRECONDITIONER_PAIR_BLOCK_SIZE
        ),
        "preconditioner_setup_seconds": float(preconditioner_setup_seconds),
        "preconditioner_data_setup_count": 1,
        "preconditioner_entropy_diagonal_updated_each_iteration": True,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "continuation": {
            "schedule": schedule,
            "stages": continuation_info,
            "final_stage_reached": bool(final_stage_reached),
        },
        "globalization": {
            "proximal_warmup_iterations": int(proximal_warmup_iterations),
            "proximal_switch_relative_gradient": float(
                proximal_switch_relative_gradient
            ),
            "continuation_activation_relative_gradient": float(
                continuation_activation_relative_gradient
            ),
            "initial_final_relative_gradient": float(
                initial_final_relative_gradient
            ),
            "cg_minimum_tolerance": float(cg_relative_tolerance),
            "cg_maximum_forcing_tolerance": float(cg_forcing_max),
            "use_whitened_newton": bool(use_whitened_newton),
            "exact_line_search": bool(exact_line_search),
            "line_search_max_step": float(line_search_max_step),
        },
        "wall_seconds": time.perf_counter() - start_time,
        "objective_definition": (
            "-1.5*logdet(B_internal) + 0.5*sum_unique_pairs("
            "(D_fit-D_obs)^2/noise_variance)"
        ),
        "stationarity_definition": (
            "D_fit_ij-D_obs_ij-noise_variance_ij*A_ij/2"
        ),
        "logged_metric_timing": "post-update accepted covariance iterate",
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance stopped without satisfying the "
            f"first-order optimality tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, gram, connectivity, info


def install_covariance_globalization() -> None:
    """Install this solver into ``hipps_dimes.numerics`` before API import."""
    if not hasattr(_num, "_legacy_fit_gaussian_noise_covariance"):
        _num._legacy_fit_gaussian_noise_covariance = (
            _num.fit_gaussian_noise_covariance
        )
    _num.fit_gaussian_noise_covariance = fit_gaussian_noise_covariance
