"""Production globalized solver for the noise-aware covariance objective."""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import covariance_geometry as _geo
from . import covariance_globalization as _base
from . import numerics as _num


_DEFAULT_CONTINUATION_FACTORS = (0.1, 1.0)


def _scalar(value):
    return float(value.item()) if hasattr(value, "item") else float(value)


def _state(
    gram,
    basis,
    target_matrix,
    inverse_variance_matrix,
    pair_i,
    pair_j,
    factor,
    parameterization,
    xp,
):
    """Evaluate one continuation-stage objective and gradient."""
    gram = xp.asarray(gram, dtype=xp.float64)
    if xp is np:
        try:
            factorization = np.linalg.cholesky(gram)
        except np.linalg.LinAlgError:
            return np.inf, None, None, None, None
        inverse = _num.scipy.linalg.cho_solve(
            (factorization, True), np.eye(gram.shape[0]), check_finite=False
        )
    else:
        try:
            with _num.cupyx.errstate(linalg="raise"):
                factorization = xp.linalg.cholesky(gram)
        except np.linalg.LinAlgError:
            return np.inf, None, None, None, None
        identity = xp.eye(gram.shape[0], dtype=xp.float64)
        inverse = _num.cupyx_scipy_linalg.solve_triangular(
            factorization, identity, lower=True, check_finite=False
        )
        inverse = _num.cupyx_scipy_linalg.solve_triangular(
            factorization.T, inverse, lower=False, check_finite=False
        )
    inverse = _geo.symmetrize(inverse)
    stage_weights = factor * inverse_variance_matrix
    if parameterization == "anchored":
        data_objective, data_gradient, fitted = (
            _geo.anchored_data_objective_gradient(
                gram,
                target_matrix,
                stage_weights,
                pair_i,
                pair_j,
                array_module=xp,
            )
        )
        gauge_constant = 1.5 * np.log(target_matrix.shape[0])
    else:
        data_objective, data_gradient, fitted = (
            _num._gaussian_covariance_data_objective_gradient(
                gram,
                basis,
                target_matrix,
                stage_weights,
                pair_i,
                pair_j,
                array_module=xp,
            )
        )
        gauge_constant = 0.0
    logdet = 2.0 * xp.sum(xp.log(xp.diag(factorization)))
    negative_entropy = -1.5 * logdet + gauge_constant
    entropy_gradient = -1.5 * inverse
    gradient = _geo.symmetrize(data_gradient + entropy_gradient)
    return (
        _scalar(negative_entropy + data_objective),
        gradient,
        fitted,
        inverse,
        (
            _scalar(negative_entropy),
            _scalar(data_objective),
            entropy_gradient,
            data_gradient,
        ),
    )


def _data_hessian(direction, basis, weights, factor, parameterization, xp):
    if parameterization == "anchored":
        return factor * _geo.anchored_data_hessian_action(
            direction, weights, array_module=xp
        )
    return factor * _num._gaussian_covariance_data_hessian_action(
        direction, basis, weights, array_module=xp
    )


def _adjust_calibration(calibration, parameterization, n):
    if parameterization != "anchored":
        return calibration
    calibration = dict(calibration)
    constant = 1.5 * np.log(n)
    calibration["objective_before"] += constant
    calibration["objective_after"] += constant
    calibration["objective_basis_constant"] = float(constant)
    return calibration


def _proximal_warmup(
    gram,
    state,
    basis,
    target_matrix,
    weights,
    pair_i,
    pair_j,
    factor,
    parameterization,
    max_steps,
    switch_tolerance,
    lipschitz,
    xp,
):
    """Take monotone proximal-gradient steps before Newton."""
    records = []
    lipschitz = max(float(lipschitz), np.finfo(np.float64).tiny)
    for _ in range(int(max_steps)):
        _, _, relative_gradient = _base._relative_gradient_from_state(state, xp)
        if relative_gradient <= switch_tolerance:
            break
        objective = float(state[0])
        data_objective = float(state[4][1])
        data_gradient = state[4][3]
        for backtracks in range(60):
            step = 1.0 / lipschitz
            candidate = _base._logdet_proximal_step(
                gram - step * data_gradient, step, xp
            )
            candidate_state = _state(
                candidate,
                basis,
                target_matrix,
                weights,
                pair_i,
                pair_j,
                factor,
                parameterization,
                xp,
            )
            delta = candidate - gram
            majorizer = (
                data_objective
                + _scalar(xp.sum(data_gradient * delta))
                + 0.5 * lipschitz * _scalar(xp.sum(delta * delta))
            )
            slack = 1e-12 * max(
                1.0, abs(objective), abs(float(candidate_state[0]))
            )
            if (
                float(candidate_state[4][1]) <= majorizer + slack
                and float(candidate_state[0]) <= objective + slack
            ):
                break
            lipschitz *= 2.0
        else:
            break
        relative_step = _base._frobenius_norm(delta, xp) / max(
            _base._frobenius_norm(gram, xp), 1.0
        )
        gram, state = candidate, candidate_state
        records.append((gram.copy(), state, relative_step, step, backtracks))
        lipschitz = max(0.8 * lipschitz, np.finfo(np.float64).tiny)
    return gram, state, records


def _history():
    keys = (
        "iteration objective stage_objective negative_entropy_objective "
        "data_objective stage_data_objective loss entropy weighted_rmse "
        "distance_rmse gradient_norm relative_gradient_norm newton_decrement "
        "relative_step step_size cg_iterations cg_relative_residual "
        "cg_forcing_tolerance cg_converged minimum_internal_gram_eigenvalue "
        "maximum_internal_gram_eigenvalue connectivity_offdiagonal_l2 "
        "continuation_factor phase line_search_backtracks full_step_accepted "
        "line_search_method feasible_step_bound fallback_direction"
    ).split()
    return {key: [] for key in keys}


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
    coordinate_parameterization="anchored",
):
    """Minimize the calibrated Gaussian covariance objective globally.

    The final objective is unchanged. The default anchored parameterization
    makes distance-gradient, Hessian-action, and data-preconditioner setup
    O(N^2), while continuation and proximal warm-up globalize Newton-CG.
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
    if max_iterations <= 0:
        raise ValueError("max_iterations must be positive")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not (0.0 < cg_relative_tolerance < 1.0):
        raise ValueError("cg_relative_tolerance must lie in (0, 1)")
    if not (cg_relative_tolerance <= cg_forcing_max < 1.0):
        raise ValueError("cg_forcing_max must lie in [cg_relative_tolerance, 1)")
    if coordinate_parameterization not in {"anchored", "centered"}:
        raise ValueError(
            "coordinate_parameterization must be 'anchored' or 'centered'"
        )
    schedule = _base._validate_continuation_factors(continuation_factors)
    if max_iterations < 10:
        schedule = (1.0,)

    n = observed.shape[0]
    m = n - 1
    if cg_max_iterations is None:
        cg_max_iterations = 2 * m
    save_steps_set = set(save_steps or ())
    if any(step < 1 or step > max_iterations for step in save_steps_set):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    orthonormal_basis = _num._centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = (
        _num._gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    orthonormal_gram, initialization_info = (
        _num._initialize_gaussian_reduced_gram(
            observed,
            pair_mask,
            pair_i,
            pair_j,
            inverse_variance,
            orthonormal_basis,
            initialization,
            initial_connectivity,
            initial_gram_floor_relative,
            use_gpu=use_gpu,
            progress_callback=progress_callback,
        )
    )
    if coordinate_parameterization == "anchored":
        centered = orthonormal_basis @ orthonormal_gram @ orthonormal_basis.T
        gram = _geo.centered_to_anchored_gram(centered)
        basis = _geo.anchored_basis(n)
    else:
        gram, basis = orthonormal_gram, orthonormal_basis
    initialization_info["coordinate_parameterization"] = coordinate_parameterization

    if use_gpu:
        xp = _num.cp
        solver_basis = xp.asarray(basis, dtype=xp.float64)
        solver_pair_i = xp.asarray(pair_i, dtype=xp.int32)
        solver_pair_j = xp.asarray(pair_j, dtype=xp.int32)
        solver_targets = xp.asarray(target_pairs, dtype=xp.float64)
        solver_weights = xp.asarray(inverse_variance, dtype=xp.float64)
        solver_target_matrix = xp.asarray(target_matrix, dtype=xp.float64)
        solver_weight_matrix = xp.asarray(
            inverse_variance_matrix, dtype=xp.float64
        )
        gram = xp.asarray(gram, dtype=xp.float64)
        cg_solver = _num._preconditioned_conjugate_gradient_gpu
    else:
        xp = np
        solver_basis = basis
        solver_pair_i, solver_pair_j = pair_i, pair_j
        solver_targets, solver_weights = target_pairs, inverse_variance
        solver_target_matrix = target_matrix
        solver_weight_matrix = inverse_variance_matrix
        cg_solver = _num._preconditioned_conjugate_gradient

    gram, calibration = _num._calibrate_gaussian_covariance_initial_scale(
        gram,
        solver_basis,
        solver_pair_i,
        solver_pair_j,
        solver_targets,
        solver_weights,
        array_module=xp,
    )
    calibration = _adjust_calibration(calibration, coordinate_parameterization, n)
    initialization_info["scalar_calibration"] = calibration
    initialization_info["wall_seconds"] = (
        initialization_info.get("wall_seconds", 0.0)
        + calibration["wall_seconds"]
    )
    if initialization_info["kind"] == "rouse":
        initialization_info["effective_spring_constant"] = (
            initialization_info["spring_constant"]
            * calibration["connectivity_scale_factor"]
        )

    setup_start = time.perf_counter()
    if coordinate_parameterization == "anchored":
        data_diagonal = _geo.anchored_data_svec_diagonal(
            solver_weight_matrix, xp
        )
        if use_gpu:
            xp.cuda.get_current_stream().synchronize()
        preconditioner_setup_seconds = time.perf_counter() - setup_start
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_preconditioner",
                    "iteration": 1,
                    "total": 1,
                    "pairs_completed": int(pair_i.size),
                    "pair_count": int(pair_i.size),
                    "method": "COV",
                    "general_method": "covariance_optimization",
                    "use_gpu": bool(use_gpu),
                    "noisy": True,
                    "coordinate_parameterization": "anchored",
                }
            )
    else:
        full_diagonal, preconditioner_setup_seconds = (
            _num._gaussian_covariance_data_preconditioner_diagonal(
                solver_basis,
                solver_pair_i,
                solver_pair_j,
                solver_weights,
                progress_callback=progress_callback,
                array_module=xp,
                use_gpu=use_gpu,
            )
        )
        data_diagonal = _geo.full_diagonal_to_svec(full_diagonal, xp)

    lipschitz = max(
        4.0 * _scalar(xp.max(xp.sum(solver_weight_matrix, axis=1))),
        np.finfo(np.float64).tiny,
    )
    initial_state = _state(
        gram,
        solver_basis,
        solver_target_matrix,
        solver_weight_matrix,
        solver_pair_i,
        solver_pair_j,
        1.0,
        coordinate_parameterization,
        xp,
    )
    _, _, initial_relative_gradient = _base._relative_gradient_from_state(
        initial_state, xp
    )
    if initial_relative_gradient <= continuation_activation_relative_gradient:
        schedule = (1.0,)

    history = _history()
    checkpoints = {}
    stage_info = []
    iteration = 0
    converged = False
    status = "max_iterations"
    message = "maximum number of iterations reached"
    state = None
    final_stage_reached = False
    start_time = time.perf_counter()

    def record(
        phase,
        factor,
        relative_step,
        step_size,
        cg_iterations=0,
        cg_residual=0.0,
        forcing=np.nan,
        cg_converged=True,
        directional_derivative=0.0,
        backtracks=0,
        line_method="none",
        feasible_bound=np.inf,
        fallback=False,
    ):
        nonlocal iteration
        iteration += 1
        objective, gradient, fitted, inverse, components = state
        negative_entropy, stage_data, entropy_gradient, data_gradient = components
        gradient_norm, _, relative_gradient = _base._relative_gradient_from_state(
            state, xp
        )
        residual = fitted - solver_targets
        residual_squared = xp.square(residual)
        final_data = float(stage_data) / float(factor)
        final_objective = float(negative_entropy) + final_data
        if use_gpu:
            with _num.cupyx.errstate(linalg="raise"):
                eigenvalues = xp.linalg.eigvalsh(gram)
        else:
            eigenvalues = np.linalg.eigvalsh(gram)
        entropy = _scalar(xp.sum(xp.log(eigenvalues))) - m * np.log(3.0)
        if coordinate_parameterization == "anchored":
            entropy -= np.log(n)
            full_connectivity = _geo.anchored_connectivity_from_inverse(
                inverse, xp
            )
            diagonal = xp.diag(full_connectivity)
            offdiag_sq = 0.5 * max(
                _scalar(xp.sum(full_connectivity**2))
                - _scalar(xp.sum(diagonal**2)),
                0.0,
            )
        else:
            reduced_precision = 3.0 * inverse
            basis_precision = solver_basis @ reduced_precision
            diagonal = -xp.sum(basis_precision * solver_basis, axis=1)
            offdiag_sq = 0.5 * max(
                _scalar(xp.sum(reduced_precision**2))
                - _scalar(xp.sum(diagonal**2)),
                0.0,
            )
        values = {
            "iteration": iteration,
            "objective": final_objective,
            "stage_objective": float(objective),
            "negative_entropy_objective": float(negative_entropy),
            "data_objective": final_data,
            "stage_data_objective": float(stage_data),
            "loss": _scalar(
                xp.sqrt(xp.mean(xp.square(residual / solver_targets)))
            ),
            "entropy": entropy,
            "weighted_rmse": _scalar(
                xp.sqrt(xp.mean(residual_squared * solver_weights))
            ),
            "distance_rmse": _scalar(xp.sqrt(xp.mean(residual_squared))),
            "gradient_norm": gradient_norm,
            "relative_gradient_norm": relative_gradient,
            "newton_decrement": float(
                np.sqrt(max(-directional_derivative, 0.0))
            ),
            "relative_step": float(relative_step),
            "step_size": float(step_size),
            "cg_iterations": int(cg_iterations),
            "cg_relative_residual": float(cg_residual)
            / max(gradient_norm, np.finfo(np.float64).tiny),
            "cg_forcing_tolerance": float(forcing),
            "cg_converged": bool(cg_converged),
            "minimum_internal_gram_eigenvalue": _scalar(eigenvalues[0]),
            "maximum_internal_gram_eigenvalue": _scalar(eigenvalues[-1]),
            "connectivity_offdiagonal_l2": float(np.sqrt(offdiag_sq)),
            "continuation_factor": float(factor),
            "phase": phase,
            "line_search_backtracks": int(backtracks),
            "full_step_accepted": bool(abs(step_size - 1.0) <= 1e-12),
            "line_search_method": line_method,
            "feasible_step_bound": float(feasible_bound),
            "fallback_direction": bool(fallback),
        }
        for key, value in values.items():
            history[key].append(value)
        if iteration in save_steps_set:
            checkpoint = _num.cp.asnumpy(gram) if use_gpu else gram
            if coordinate_parameterization == "anchored":
                checkpoints[iteration] = _geo.anchored_connectivity_from_gram(
                    checkpoint
                )
            else:
                checkpoints[iteration] = _num._connectivity_from_reduced_gram(
                    checkpoint, basis
                )
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": phase,
                    "continuation_factor": float(factor),
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": final_objective,
                    "stage_objective": float(objective),
                    "loss": values["loss"],
                    "entropy": entropy,
                    "relative_gradient_norm": relative_gradient,
                    "step_size": float(step_size),
                    "line_search_backtracks": int(backtracks),
                    "full_step_accepted": values["full_step_accepted"],
                    "cg_iterations": int(cg_iterations),
                    "cg_converged": bool(cg_converged),
                    "cg_forcing_tolerance": float(forcing),
                    "coordinate_parameterization": coordinate_parameterization,
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV",
                    "general_method": "covariance_optimization",
                }
            )

    reserve = min(10, max_iterations)
    for factor in schedule:
        if iteration >= max_iterations:
            break
        final_stage = factor == 1.0
        final_stage_reached |= final_stage
        stage_start = iteration
        stage_weights = factor * solver_weights
        gram, stage_calibration = (
            _num._calibrate_gaussian_covariance_initial_scale(
                gram,
                solver_basis,
                solver_pair_i,
                solver_pair_j,
                solver_targets,
                stage_weights,
                array_module=xp,
            )
        )
        stage_calibration = _adjust_calibration(
            stage_calibration, coordinate_parameterization, n
        )
        state = _state(
            gram,
            solver_basis,
            solver_target_matrix,
            solver_weight_matrix,
            solver_pair_i,
            solver_pair_j,
            factor,
            coordinate_parameterization,
            xp,
        )
        tolerance = (
            relative_tolerance
            if final_stage
            else max(relative_tolerance, continuation_intermediate_tolerance)
        )
        if final_stage:
            budget = max_iterations - iteration
        else:
            budget = max(
                0,
                min(
                    proximal_warmup_iterations
                    + continuation_intermediate_newton_iterations,
                    max_iterations - iteration - reserve,
                ),
            )
        stage_record = {
            "factor": float(factor),
            "scalar_calibration": stage_calibration,
            "start_iteration": iteration + 1,
            "warmup_steps": 0,
            "newton_steps": 0,
            "converged": False,
        }
        gram, state, warmup = _proximal_warmup(
            gram,
            state,
            solver_basis,
            solver_target_matrix,
            solver_weight_matrix,
            solver_pair_i,
            solver_pair_j,
            factor,
            coordinate_parameterization,
            min(proximal_warmup_iterations, budget),
            proximal_switch_relative_gradient,
            factor * lipschitz,
            xp,
        )
        for warm_gram, warm_state, relative_step, step, backtracks in warmup:
            gram, state = warm_gram, warm_state
            record(
                "proximal_warmup",
                factor,
                relative_step,
                step,
                backtracks=backtracks,
                line_method="logdet_proximal_backtracking",
            )
            stage_record["warmup_steps"] += 1
        budget -= stage_record["warmup_steps"]

        while budget > 0 and iteration < max_iterations:
            objective, gradient, _, inverse, components = state
            gradient_norm, gradient_scale, relative_gradient = (
                _base._relative_gradient_from_state(state, xp)
            )
            if gradient_norm <= absolute_tolerance + tolerance * gradient_scale:
                stage_record["converged"] = True
                break
            forcing = _base._inexact_newton_forcing_tolerance(
                relative_gradient, cg_relative_tolerance, cg_forcing_max
            )
            fallback = False
            if use_whitened_newton:
                if xp is np:
                    chol = np.linalg.cholesky(gram)
                else:
                    with _num.cupyx.errstate(linalg="raise"):
                        chol = xp.linalg.cholesky(gram)
                transformed_gradient = chol.T @ gradient @ chol

                def operator(matrix_direction):
                    matrix_direction = _geo.symmetrize(matrix_direction)
                    gram_direction = chol @ matrix_direction @ chol.T
                    data_action = _data_hessian(
                        gram_direction,
                        solver_basis,
                        solver_weight_matrix,
                        factor,
                        coordinate_parameterization,
                        xp,
                    )
                    return _geo.symmetrize(
                        chol.T @ data_action @ chol + 1.5 * matrix_direction
                    )

                scale = xp.square(xp.diag(chol))
                preconditioner = xp.maximum(
                    1.5
                    + factor
                    * data_diagonal
                    * xp.outer(scale, scale),
                    np.finfo(np.float64).tiny,
                )
                white_direction, cg_iterations, cg_residual, cg_ok = cg_solver(
                    operator,
                    -transformed_gradient,
                    preconditioner,
                    forcing,
                    int(cg_max_iterations),
                )
                white_direction = _geo.symmetrize(white_direction)
                direction = _geo.symmetrize(chol @ white_direction @ chol.T)
                derivative = _scalar(xp.sum(gradient * direction))
                if not np.isfinite(derivative) or derivative >= 0.0:
                    fallback = True
                    white_direction = _geo.symmetrize(
                        -transformed_gradient / preconditioner
                    )
                    direction = _geo.symmetrize(
                        chol @ white_direction @ chol.T
                    )
                    derivative = _scalar(xp.sum(gradient * direction))
                    cg_ok = False
            else:

                def operator(matrix_direction):
                    data_action = _data_hessian(
                        matrix_direction,
                        solver_basis,
                        solver_weight_matrix,
                        factor,
                        coordinate_parameterization,
                        xp,
                    )
                    return _geo.symmetrize(
                        data_action + 1.5 * inverse @ matrix_direction @ inverse
                    )

                preconditioner = xp.maximum(
                    factor * data_diagonal
                    + _geo.entropy_svec_diagonal(inverse, xp),
                    np.finfo(np.float64).tiny,
                )
                direction, cg_iterations, cg_residual, cg_ok = cg_solver(
                    operator,
                    -gradient,
                    preconditioner,
                    forcing,
                    int(cg_max_iterations),
                )
                direction = _geo.symmetrize(direction)
                derivative = _scalar(xp.sum(gradient * direction))
                if not np.isfinite(derivative) or derivative >= 0.0:
                    fallback = True
                    direction = _geo.symmetrize(-gradient / preconditioner)
                    derivative = _scalar(xp.sum(gradient * direction))
                    cg_ok = False
                if xp is np:
                    chol = np.linalg.cholesky(gram)
                    white_direction = np.linalg.solve(
                        chol, np.linalg.solve(chol, direction).T
                    ).T
                else:
                    identity = xp.eye(m, dtype=xp.float64)
                    inverse_chol = _num.cupyx_scipy_linalg.solve_triangular(
                        xp.linalg.cholesky(gram),
                        identity,
                        lower=True,
                        check_finite=False,
                    )
                    white_direction = (
                        inverse_chol @ direction @ inverse_chol.T
                    )
                white_direction = _geo.symmetrize(white_direction)

            def evaluate(trial):
                return _state(
                    trial,
                    solver_basis,
                    solver_target_matrix,
                    solver_weight_matrix,
                    solver_pair_i,
                    solver_pair_j,
                    factor,
                    coordinate_parameterization,
                    xp,
                )

            line_method = "exact_feasible_scalar"
            backtracks = 0
            feasible_bound = np.inf
            candidate = candidate_state = None
            if exact_line_search and derivative < 0.0:
                data_direction = _data_hessian(
                    direction,
                    solver_basis,
                    solver_weight_matrix,
                    factor,
                    coordinate_parameterization,
                    xp,
                )
                step, feasible_bound = _base._feasible_exact_scalar_step(
                    gram,
                    direction,
                    white_direction,
                    components[3],
                    data_direction,
                    array_module=xp,
                    maximum_step=line_search_max_step,
                )
                if step > 0.0:
                    candidate = _geo.symmetrize(gram + step * direction)
                    candidate_state = evaluate(candidate)
                    slack = 1e-11 * max(
                        1.0, abs(float(objective)), abs(float(candidate_state[0]))
                    )
                    if (
                        not np.isfinite(candidate_state[0])
                        or float(candidate_state[0]) > float(objective) + slack
                    ):
                        candidate = candidate_state = None
            if candidate is None:
                line_method = "armijo_fallback"
                initial_step = min(
                    1.0,
                    0.99 * feasible_bound
                    if np.isfinite(feasible_bound)
                    else 1.0,
                )
                candidate, candidate_state, step, backtracks = (
                    _base._armijo_fallback(
                        gram,
                        direction,
                        float(objective),
                        derivative,
                        evaluate,
                        initial_step,
                    )
                )
            if candidate is None:
                status = "line_search_failed"
                message = "covariance-cone line search failed"
                budget = 0
                break
            delta = candidate - gram
            relative_step = _base._frobenius_norm(delta, xp) / max(
                _base._frobenius_norm(gram, xp), 1.0
            )
            gram, state = candidate, candidate_state
            record(
                "newton",
                factor,
                relative_step,
                step,
                cg_iterations,
                cg_residual,
                forcing,
                cg_ok,
                derivative,
                backtracks,
                line_method,
                feasible_bound,
                fallback,
            )
            stage_record["newton_steps"] += 1
            budget -= 1

        gradient_norm, gradient_scale, relative_gradient = (
            _base._relative_gradient_from_state(state, xp)
        )
        stage_record["final_relative_gradient_norm"] = relative_gradient
        stage_record["converged"] = gradient_norm <= (
            absolute_tolerance + tolerance * gradient_scale
        )
        stage_record["end_iteration"] = iteration
        stage_record["accepted_steps"] = iteration - stage_start
        stage_info.append(stage_record)
        if final_stage and stage_record["converged"]:
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break
        if status == "line_search_failed":
            break

    state = _state(
        gram,
        solver_basis,
        solver_target_matrix,
        solver_weight_matrix,
        solver_pair_i,
        solver_pair_j,
        1.0,
        coordinate_parameterization,
        xp,
    )
    final_gradient_norm, final_gradient_scale, final_relative_gradient = (
        _base._relative_gradient_from_state(state, xp)
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
        gram = _num.cp.asnumpy(gram)
    if coordinate_parameterization == "anchored":
        full_gram = _geo.anchored_to_centered_gram(gram)
        fitted = _geo.anchored_squared_distances(gram, np)
        connectivity = _geo.anchored_connectivity_from_gram(gram)
    else:
        full_gram = _geo.symmetrize(basis @ gram @ basis.T)
        fitted = _num._squared_distances_from_gram(full_gram)
        connectivity = _num._connectivity_from_reduced_gram(gram, basis)
    fitted_pairs = fitted[pair_i, pair_j]
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
        "iterations": int(iteration),
        "objective": float(state[0]),
        "gradient_norm": final_gradient_norm,
        "relative_gradient_norm": final_relative_gradient,
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
        "preconditioner_kind": "svec_diagonal",
        "coordinate_parameterization": coordinate_parameterization,
        "history": history,
        "connectivity_matrix_at_steps": checkpoints,
        "continuation": {
            "schedule": schedule,
            "stages": stage_info,
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
            "initial_final_relative_gradient": float(initial_relative_gradient),
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
    return fitted, full_gram, connectivity, info


def install_covariance_solver():
    """Install this implementation while preserving the legacy reference."""
    if not hasattr(_num, "_legacy_fit_gaussian_noise_covariance"):
        _num._legacy_fit_gaussian_noise_covariance = (
            _num.fit_gaussian_noise_covariance
        )
    _num.fit_gaussian_noise_covariance = fit_gaussian_noise_covariance
