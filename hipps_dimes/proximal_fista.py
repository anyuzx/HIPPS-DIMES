"""Exact proximal-gradient and monotone FISTA for noisy HIPPS-DIMES.

This module optimizes the Gaussian-noise maximum-entropy dual in the natural
pair parameters ``theta_ij``.  The HIPPS-DIMES connectivity convention is

    theta_ij = -A_ij / 2,  i < j,

where positive ``A_ij`` is an attractive effective spring.  For independent
Gaussian uncertainty with variance ``v_ij`` in the observed mean-squared
pair distances, the objective is

    L(theta) = log Z(theta) - sum_ij theta_ij Dobs_ij
               + 1/2 sum_ij v_ij theta_ij**2.

The exact smooth gradient is

    grad f(theta)_ij = Dfit_ij(theta) - Dobs_ij,

and the quadratic proximal map is componentwise

    prox_{eta g}(x)_ij = x_ij / (1 + eta v_ij).

The implementation uses backtracking, spectral validity checks, monotone
FISTA restart, and an independently evaluated KKT residual.  Unlike the
legacy iterative-scaling-plus-shrinkage update, every accepted update is a
true proximal-gradient step for the stated objective.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import numerics as _num


def _scalar(value):
    """Convert a NumPy/CuPy scalar to a Python float."""
    return float(value.item()) if hasattr(value, "item") else float(value)


def _norm(value, xp):
    return _scalar(xp.sqrt(xp.sum(value * value)))


def _symmetrize(matrix):
    return 0.5 * (matrix + matrix.T)


def _centered_orthonormal_basis(n):
    """Return a Helmert basis for the center-of-mass-free subspace."""
    basis = np.zeros((n, n - 1), dtype=np.float64)
    for column in range(n - 1):
        denominator = np.sqrt((column + 1.0) * (column + 2.0))
        basis[: column + 1, column] = 1.0 / denominator
        basis[column + 1, column] = -(column + 1.0) / denominator
    return basis


def _validate_inputs(
    squared_distances,
    noise_variance,
    relative_noise_std,
):
    observed = np.asarray(squared_distances, dtype=np.float64)
    if observed.ndim != 2 or observed.shape[0] != observed.shape[1]:
        raise ValueError("squared_distances must be a square matrix")
    if observed.shape[0] < 2:
        raise ValueError("at least two loci are required")

    n = observed.shape[0]
    pair_i, pair_j = np.triu_indices(n, k=1)
    finite_forward = np.isfinite(observed[pair_i, pair_j])
    finite_reverse = np.isfinite(observed[pair_j, pair_i])
    pair_mask = finite_forward & finite_reverse
    pair_i = pair_i[pair_mask]
    pair_j = pair_j[pair_mask]
    if pair_i.size == 0:
        raise ValueError("no finite off-diagonal distance constraints were found")

    forward = observed[pair_i, pair_j]
    reverse = observed[pair_j, pair_i]
    symmetry_scale = max(float(np.max(np.abs(np.concatenate((forward, reverse))))), 1.0)
    if not np.allclose(forward, reverse, rtol=1e-10, atol=1e-12 * symmetry_scale):
        raise ValueError("squared_distances must be symmetric on observed pairs")
    target = 0.5 * (forward + reverse)

    if relative_noise_std is not None and noise_variance is not None:
        raise ValueError(
            "supply either noise_variance or relative_noise_std, not both"
        )
    if relative_noise_std is not None:
        relative_noise_std = float(relative_noise_std)
        if not np.isfinite(relative_noise_std) or relative_noise_std <= 0.0:
            raise ValueError("relative_noise_std must be positive and finite")
        if np.any(target == 0.0):
            raise ValueError(
                "relative_noise_std is undefined for zero observed distances"
            )
        pair_variance = np.square(relative_noise_std * target)
        noise_model = "heteroskedastic_relative_std"
        noise_parameter = relative_noise_std
    else:
        if noise_variance is None:
            raise ValueError("noise_variance or relative_noise_std is required")
        variance_array = np.asarray(noise_variance, dtype=np.float64)
        if variance_array.ndim == 0:
            scalar_variance = float(variance_array)
            if not np.isfinite(scalar_variance) or scalar_variance <= 0.0:
                raise ValueError("noise_variance must be positive and finite")
            pair_variance = np.full(target.shape, scalar_variance, dtype=np.float64)
            noise_model = "homoskedastic"
            noise_parameter = scalar_variance
        else:
            if variance_array.shape != observed.shape:
                raise ValueError(
                    "matrix noise_variance must have the same shape as squared_distances"
                )
            forward_variance = variance_array[pair_i, pair_j]
            reverse_variance = variance_array[pair_j, pair_i]
            if not np.allclose(
                forward_variance,
                reverse_variance,
                rtol=1e-10,
                atol=1e-14,
            ):
                raise ValueError("matrix noise_variance must be symmetric")
            pair_variance = 0.5 * (forward_variance + reverse_variance)
            if np.any(~np.isfinite(pair_variance)) or np.any(pair_variance <= 0.0):
                raise ValueError(
                    "noise variances must be positive and finite on observed pairs"
                )
            noise_model = "heteroskedastic_variance_matrix"
            noise_parameter = None

    return (
        observed,
        pair_i.astype(np.int64),
        pair_j.astype(np.int64),
        target.astype(np.float64),
        pair_variance.astype(np.float64),
        noise_model,
        noise_parameter,
    )


def _connectivity_from_theta(theta, fixed_offdiagonal, pair_i, pair_j, xp):
    """Build the full HIPPS-DIMES connectivity matrix from natural parameters."""
    connectivity = fixed_offdiagonal.copy()
    edge_values = -2.0 * theta
    connectivity[pair_i, pair_j] = edge_values
    connectivity[pair_j, pair_i] = edge_values
    connectivity[xp.diag_indices(connectivity.shape[0])] = 0.0
    diagonal = -xp.sum(connectivity, axis=1)
    connectivity[xp.diag_indices(connectivity.shape[0])] = diagonal
    return _symmetrize(connectivity)


def _rouse_initial_connectivity(n, pair_i, pair_j, target):
    positive = target > 0.0
    if not np.any(positive):
        raise ValueError(
            "Rouse initialization requires at least one positive observed distance"
        )
    mean_target = float(np.mean(target[positive]))
    mean_separation = float(np.mean((pair_j - pair_i)[positive]))
    spring_constant = 3.0 * mean_separation / mean_target
    connectivity = np.zeros((n, n), dtype=np.float64)
    indices = np.arange(n - 1)
    connectivity[indices, indices + 1] = spring_constant
    connectivity[indices + 1, indices] = spring_constant
    np.fill_diagonal(connectivity, -np.sum(connectivity, axis=1))
    return connectivity, spring_constant


def _prepare_initial_connectivity(
    n,
    pair_i,
    pair_j,
    target,
    initial_connectivity,
):
    if initial_connectivity is None:
        connectivity, spring_constant = _rouse_initial_connectivity(
            n, pair_i, pair_j, target
        )
        initialization = {
            "kind": "rouse",
            "spring_constant": spring_constant,
        }
    else:
        connectivity = np.asarray(initial_connectivity, dtype=np.float64)
        if connectivity.shape != (n, n):
            raise ValueError(
                "initial_connectivity must match squared_distances shape"
            )
        if not np.all(np.isfinite(connectivity)):
            raise ValueError("initial_connectivity must be finite")
        connectivity = _symmetrize(connectivity.copy())
        np.fill_diagonal(connectivity, 0.0)
        np.fill_diagonal(connectivity, -np.sum(connectivity, axis=1))
        initialization = {"kind": "provided_connectivity"}
    return connectivity, initialization


def _evaluate_state(
    theta,
    fixed_offdiagonal,
    basis,
    pair_i,
    pair_j,
    target,
    variance,
    xp,
):
    """Evaluate the exact dual objective and gradient at one valid parameter vector."""
    connectivity = _connectivity_from_theta(
        theta, fixed_offdiagonal, pair_i, pair_j, xp
    )
    precision = _symmetrize(-basis.T @ connectivity @ basis)
    try:
        factor = xp.linalg.cholesky(precision)
    except Exception:
        return None

    diagonal_factor = xp.diag(factor)
    if _scalar(xp.min(diagonal_factor)) <= 0.0:
        return None
    logdet_precision = 2.0 * _scalar(xp.sum(xp.log(diagonal_factor)))
    identity = xp.eye(precision.shape[0], dtype=xp.float64)
    try:
        inverse_precision = xp.linalg.solve(
            factor.T,
            xp.linalg.solve(factor, identity),
        )
    except Exception:
        return None
    inverse_precision = _symmetrize(inverse_precision)

    gram = _symmetrize(3.0 * basis @ inverse_precision @ basis.T)
    gram_diagonal = xp.diag(gram)
    fitted = (
        gram_diagonal[pair_i]
        + gram_diagonal[pair_j]
        - 2.0 * gram[pair_i, pair_j]
    )

    log_partition = -1.5 * logdet_precision
    smooth_objective = log_partition - _scalar(xp.sum(theta * target))
    regularizer = 0.5 * _scalar(xp.sum(variance * theta * theta))
    objective = smooth_objective + regularizer
    smooth_gradient = fitted - target
    total_gradient = smooth_gradient + variance * theta

    return {
        "objective": float(objective),
        "smooth_objective": float(smooth_objective),
        "regularizer": float(regularizer),
        "log_partition": float(log_partition),
        "smooth_gradient": smooth_gradient,
        "total_gradient": total_gradient,
        "fitted": fitted,
        "connectivity": connectivity,
        "precision": precision,
        "inverse_precision": inverse_precision,
        "gram": gram,
    }


def _quadratic_proximal_map(point, step_size, variance, enforce_nonnegative, xp):
    result = point / (1.0 + step_size * variance)
    if enforce_nonnegative:
        # Positive HIPPS-DIMES A_ij corresponds to theta_ij <= 0.
        result = xp.minimum(result, 0.0)
    return result


def _relative_kkt(state, theta, variance, xp):
    residual = state["total_gradient"]
    residual_norm = _norm(residual, xp)
    data_part = state["smooth_gradient"]
    regularization_part = variance * theta
    scale = max(
        _norm(data_part, xp),
        _norm(regularization_part, xp),
        1.0,
    )
    return residual_norm, scale, residual_norm / scale


def _backtracking_step(
    extrapolated,
    extrapolated_state,
    current_objective,
    step_size,
    variance,
    evaluate,
    enforce_nonnegative,
    require_monotone,
    backtracking_factor,
    maximum_backtracks,
    xp,
):
    """Return one feasible exact proximal-gradient step with majorization."""
    smooth_gradient = extrapolated_state["smooth_gradient"]
    smooth_objective = extrapolated_state["smooth_objective"]
    step = float(step_size)
    for backtracks in range(maximum_backtracks + 1):
        forward = extrapolated - step * smooth_gradient
        candidate = _quadratic_proximal_map(
            forward, step, variance, enforce_nonnegative, xp
        )
        candidate_state = evaluate(candidate)
        if candidate_state is not None:
            delta = candidate - extrapolated
            quadratic_majorizer = (
                smooth_objective
                + _scalar(xp.sum(smooth_gradient * delta))
                + 0.5 * _scalar(xp.sum(delta * delta)) / step
            )
            scale = max(
                1.0,
                abs(smooth_objective),
                abs(candidate_state["smooth_objective"]),
                abs(current_objective),
            )
            slack = 1e-12 * scale
            majorization_ok = (
                candidate_state["smooth_objective"]
                <= quadratic_majorizer + slack
            )
            monotone_ok = (
                not require_monotone
                or candidate_state["objective"] <= current_objective + slack
            )
            if majorization_ok and monotone_ok:
                return candidate, candidate_state, step, backtracks
        step *= backtracking_factor
        if step <= np.finfo(np.float64).tiny:
            break
    return None, None, step, maximum_backtracks + 1


def fit_gaussian_noise_dual_fista(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    initial_step_size=None,
    maximum_step_size=None,
    step_growth=1.1,
    backtracking_factor=0.5,
    maximum_backtracks=80,
    accelerated=True,
    monotone=True,
    enforce_nonnegative_connectivity=False,
    save_steps=None,
    progress_callback=None,
):
    """Minimize the exact Gaussian-noise dual with monotone FISTA.

    Parameters
    ----------
    squared_distances : array_like, shape (N, N)
        Observed mean-squared distance map. Finite off-diagonal pairs are used.
        The map may be an invalid EDM.
    noise_variance : float or array_like
        Positive variance of each observed squared-distance constraint. A scalar
        gives homoskedastic noise; a symmetric matrix gives pair-specific noise.
    relative_noise_std : float, optional
        Alternative to ``noise_variance``. Uses
        ``variance_ij = (relative_noise_std * Dobs_ij)**2``.
    initial_connectivity : array_like, optional
        Valid HIPPS-DIMES connectivity used as initialization. If omitted, an
        observed-pair-calibrated Rouse chain is used.
    use_gpu : bool, default=False
        Use CuPy float64 calculations when an accessible CUDA GPU is available.
    max_iterations : int, default=1000
        Maximum accepted proximal-gradient updates.
    relative_tolerance : float, default=1e-5
        Relative exact KKT tolerance.
    initial_step_size : float, optional
        Proximal-gradient step in natural-parameter units. If omitted, a
        dimensionally scaled value based on the median target MSD is used.
    accelerated : bool, default=True
        Use FISTA extrapolation. ``False`` gives exact proximal gradient.
    monotone : bool, default=True
        Restart extrapolation whenever the proposed objective exceeds the last
        accepted objective.

    Returns
    -------
    fitted_squared_distances : ndarray
        Full fitted mean-squared-distance map.
    gram : ndarray
        Centered three-dimensional Gram covariance.
    connectivity : ndarray
        Final HIPPS-DIMES connectivity matrix.
    info : dict
        Convergence certificate, history, initialization, and checkpoints.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "FISTA GPU optimization was requested, but CuPy and an accessible "
            "CUDA GPU are not available"
        )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(step_growth) or step_growth < 1.0:
        raise ValueError("step_growth must be finite and at least one")
    if (
        not np.isfinite(backtracking_factor)
        or not (0.0 < backtracking_factor < 1.0)
    ):
        raise ValueError("backtracking_factor must lie strictly between zero and one")
    if not isinstance(maximum_backtracks, (int, np.integer)) or maximum_backtracks < 0:
        raise ValueError("maximum_backtracks must be a nonnegative integer")

    (
        observed,
        pair_i,
        pair_j,
        target,
        pair_variance,
        noise_model,
        noise_parameter,
    ) = _validate_inputs(
        squared_distances,
        noise_variance,
        relative_noise_std,
    )
    n = observed.shape[0]
    basis_cpu = _centered_orthonormal_basis(n)
    initial_cpu, initialization = _prepare_initial_connectivity(
        n,
        pair_i,
        pair_j,
        target,
        initial_connectivity,
    )

    fixed_offdiagonal_cpu = _symmetrize(initial_cpu.copy())
    np.fill_diagonal(fixed_offdiagonal_cpu, 0.0)
    initial_theta_cpu = -0.5 * initial_cpu[pair_i, pair_j]

    if use_gpu:
        xp = _num.cp
        basis = xp.asarray(basis_cpu, dtype=xp.float64)
        solver_pair_i = xp.asarray(pair_i, dtype=xp.int64)
        solver_pair_j = xp.asarray(pair_j, dtype=xp.int64)
        solver_target = xp.asarray(target, dtype=xp.float64)
        solver_variance = xp.asarray(pair_variance, dtype=xp.float64)
        fixed_offdiagonal = xp.asarray(
            fixed_offdiagonal_cpu, dtype=xp.float64
        )
        theta = xp.asarray(initial_theta_cpu, dtype=xp.float64)
    else:
        xp = np
        basis = basis_cpu
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target = target
        solver_variance = pair_variance
        fixed_offdiagonal = fixed_offdiagonal_cpu
        theta = initial_theta_cpu.copy()

    typical_distance = float(np.median(np.abs(target)))
    if not np.isfinite(typical_distance) or typical_distance <= 0.0:
        typical_distance = float(np.mean(np.abs(target)))
    if not np.isfinite(typical_distance) or typical_distance <= 0.0:
        raise ValueError("cannot infer a positive distance scale for step selection")
    if initial_step_size is None:
        initial_step_size = 1.0 / (typical_distance * typical_distance)
    initial_step_size = float(initial_step_size)
    if not np.isfinite(initial_step_size) or initial_step_size <= 0.0:
        raise ValueError("initial_step_size must be positive and finite")
    if maximum_step_size is None:
        maximum_step_size = 1e6 * initial_step_size
    maximum_step_size = float(maximum_step_size)
    if not np.isfinite(maximum_step_size) or maximum_step_size < initial_step_size:
        raise ValueError(
            "maximum_step_size must be finite and at least initial_step_size"
        )

    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer))
        or step < 1
        or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must contain integers in [1, max_iterations]")

    def evaluate(value):
        return _evaluate_state(
            value,
            fixed_offdiagonal,
            basis,
            solver_pair_i,
            solver_pair_j,
            solver_target,
            solver_variance,
            xp,
        )

    state = evaluate(theta)
    if state is None:
        raise ValueError(
            "initial_connectivity is not normalizable: its internal precision "
            "is not positive definite"
        )

    previous_theta = theta.copy()
    fista_parameter = 1.0
    step_size = initial_step_size
    connectivity_at_steps = {}
    history = {
        "iteration": [],
        "objective": [],
        "smooth_objective": [],
        "regularizer": [],
        "log_partition": [],
        "loss": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "kkt_residual_norm": [],
        "relative_kkt_residual": [],
        "proximal_gradient_mapping_norm": [],
        "step_size": [],
        "momentum_coefficient": [],
        "restarted": [],
        "restart_reason": [],
        "line_search_backtracks": [],
        "minimum_internal_precision_eigenvalue": [],
        "maximum_internal_precision_eigenvalue": [],
        "precision_condition_number": [],
        "connectivity_offdiagonal_l2": [],
    }
    converged = False
    status = "max_iterations"
    message = "maximum number of FISTA iterations reached"
    start_time = time.perf_counter()

    for iteration in range(1, max_iterations + 1):
        current_theta = theta
        current_state = state
        restarted = False
        restart_reason = "none"

        if accelerated and iteration > 1:
            next_fista_parameter = 0.5 * (
                1.0 + np.sqrt(1.0 + 4.0 * fista_parameter * fista_parameter)
            )
            momentum_coefficient = (
                (fista_parameter - 1.0) / next_fista_parameter
            )
            extrapolated = current_theta + momentum_coefficient * (
                current_theta - previous_theta
            )
            extrapolated_state = evaluate(extrapolated)
            if extrapolated_state is None:
                restarted = True
                restart_reason = "invalid_extrapolation"
                momentum_coefficient = 0.0
                next_fista_parameter = 1.0
                extrapolated = current_theta
                extrapolated_state = current_state
        else:
            next_fista_parameter = 1.0
            momentum_coefficient = 0.0
            extrapolated = current_theta
            extrapolated_state = current_state

        trial_step = min(step_size * step_growth, maximum_step_size)
        candidate, candidate_state, accepted_step, backtracks = (
            _backtracking_step(
                extrapolated,
                extrapolated_state,
                current_state["objective"],
                trial_step,
                solver_variance,
                evaluate,
                enforce_nonnegative_connectivity,
                False,
                backtracking_factor,
                maximum_backtracks,
                xp,
            )
        )

        if candidate is None:
            restarted = True
            restart_reason = "line_search_failure"
        elif monotone:
            scale = max(
                1.0,
                abs(current_state["objective"]),
                abs(candidate_state["objective"]),
            )
            if candidate_state["objective"] > current_state["objective"] + 1e-12 * scale:
                restarted = True
                restart_reason = "objective_increase"

        if restarted:
            momentum_coefficient = 0.0
            next_fista_parameter = 1.0
            candidate, candidate_state, accepted_step, restart_backtracks = (
                _backtracking_step(
                    current_theta,
                    current_state,
                    current_state["objective"],
                    min(trial_step, step_size),
                    solver_variance,
                    evaluate,
                    enforce_nonnegative_connectivity,
                    True,
                    backtracking_factor,
                    maximum_backtracks,
                    xp,
                )
            )
            backtracks += restart_backtracks

        if candidate is None or candidate_state is None:
            status = "line_search_failed"
            message = "no feasible monotone proximal-gradient step was found"
            break

        previous_theta = current_theta.copy()
        theta = candidate
        state = candidate_state
        fista_parameter = next_fista_parameter
        step_size = accepted_step

        kkt_norm, kkt_scale, relative_kkt = _relative_kkt(
            state, theta, solver_variance, xp
        )
        mapping_point = _quadratic_proximal_map(
            theta - accepted_step * state["smooth_gradient"],
            accepted_step,
            solver_variance,
            enforce_nonnegative_connectivity,
            xp,
        )
        mapping_norm = _norm(theta - mapping_point, xp) / accepted_step
        residual = state["fitted"] - solver_target
        weighted_rmse = _scalar(
            xp.sqrt(xp.mean(xp.square(residual) / solver_variance))
        )
        distance_rmse = _scalar(xp.sqrt(xp.mean(xp.square(residual))))
        denominator = xp.where(
            xp.abs(solver_target) > np.finfo(np.float64).tiny,
            xp.abs(solver_target),
            1.0,
        )
        relative_loss = _scalar(
            xp.sqrt(xp.mean(xp.square(residual / denominator)))
        )
        precision_eigenvalues = xp.linalg.eigvalsh(state["precision"])
        minimum_precision = _scalar(precision_eigenvalues[0])
        maximum_precision = _scalar(precision_eigenvalues[-1])
        condition_number = maximum_precision / minimum_precision
        offdiagonal = state["connectivity"][xp.triu_indices(n, k=1)]
        connectivity_norm = _norm(offdiagonal, xp)

        history["iteration"].append(iteration)
        history["objective"].append(state["objective"])
        history["smooth_objective"].append(state["smooth_objective"])
        history["regularizer"].append(state["regularizer"])
        history["log_partition"].append(state["log_partition"])
        history["loss"].append(relative_loss)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["kkt_residual_norm"].append(kkt_norm)
        history["relative_kkt_residual"].append(relative_kkt)
        history["proximal_gradient_mapping_norm"].append(mapping_norm)
        history["step_size"].append(accepted_step)
        history["momentum_coefficient"].append(momentum_coefficient)
        history["restarted"].append(restarted)
        history["restart_reason"].append(restart_reason)
        history["line_search_backtracks"].append(backtracks)
        history["minimum_internal_precision_eigenvalue"].append(
            minimum_precision
        )
        history["maximum_internal_precision_eigenvalue"].append(
            maximum_precision
        )
        history["precision_condition_number"].append(condition_number)
        history["connectivity_offdiagonal_l2"].append(connectivity_norm)

        if iteration in save_steps_set:
            checkpoint = state["connectivity"]
            if use_gpu:
                checkpoint = _num.cp.asnumpy(checkpoint)
            connectivity_at_steps[iteration] = np.asarray(checkpoint).copy()

        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "optimization",
                    "phase": "monotone_fista" if accelerated else "proximal_gradient",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": state["objective"],
                    "loss": relative_loss,
                    "entropy": state["log_partition"],
                    "relative_gradient_norm": relative_kkt,
                    "relative_kkt_residual": relative_kkt,
                    "step_size": accepted_step,
                    "momentum_coefficient": momentum_coefficient,
                    "restarted": restarted,
                    "restart_reason": restart_reason,
                    "line_search_backtracks": backtracks,
                    "method": "FISTA" if accelerated else "PG",
                    "general_method": "optimization",
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                }
            )

        if kkt_norm <= absolute_tolerance + relative_tolerance * kkt_scale:
            converged = True
            status = "optimality_tolerance"
            message = "exact noise-aware KKT tolerance reached"
            break

    if use_gpu:
        final_gram = _num.cp.asnumpy(state["gram"])
        final_connectivity = _num.cp.asnumpy(state["connectivity"])
        fitted_pairs = _num.cp.asnumpy(state["fitted"])
        final_theta = _num.cp.asnumpy(theta)
    else:
        final_gram = np.asarray(state["gram"])
        final_connectivity = np.asarray(state["connectivity"])
        fitted_pairs = np.asarray(state["fitted"])
        final_theta = np.asarray(theta)

    fitted_squared_distances = np.zeros((n, n), dtype=np.float64)
    fitted_squared_distances[pair_i, pair_j] = fitted_pairs
    fitted_squared_distances[pair_j, pair_i] = fitted_pairs
    gram_diagonal = np.diag(final_gram)
    fitted_squared_distances = (
        gram_diagonal[:, None]
        + gram_diagonal[None, :]
        - 2.0 * final_gram
    )
    fitted_squared_distances = _symmetrize(fitted_squared_distances)
    np.fill_diagonal(fitted_squared_distances, 0.0)

    for key, values in history.items():
        if key in {"iteration", "line_search_backtracks"}:
            history[key] = np.asarray(values, dtype=np.int64)
        elif key == "restarted":
            history[key] = np.asarray(values, dtype=bool)
        elif key == "restart_reason":
            history[key] = np.asarray(values, dtype=object)
        else:
            history[key] = np.asarray(values, dtype=np.float64)

    final_kkt_norm, final_kkt_scale, final_relative_kkt = _relative_kkt(
        state, theta, solver_variance, xp
    )
    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": int(len(history["iteration"])),
        "objective": float(state["objective"]),
        "kkt_residual_norm": float(final_kkt_norm),
        "relative_kkt_residual": float(final_relative_kkt),
        "maximum_absolute_kkt_residual": float(
            np.max(
                np.abs(
                    fitted_pairs
                    - target
                    + pair_variance * final_theta
                )
            )
        ),
        "algorithm": "monotone_fista" if accelerated else "proximal_gradient",
        "parameterization": "theta_ij=-A_ij/2",
        "noise_model": noise_model,
        "noise_parameter": noise_parameter,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization,
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "wall_seconds": time.perf_counter() - start_time,
        "objective_definition": (
            "logZ(theta)-sum(theta_ij*Dobs_ij)+0.5*sum(variance_ij*theta_ij^2)"
        ),
        "stationarity_definition": (
            "Dfit_ij-Dobs_ij+variance_ij*theta_ij=0; theta_ij=-A_ij/2"
        ),
        "monotone": bool(monotone),
        "accelerated": bool(accelerated),
        "initial_step_size": float(initial_step_size),
        "final_step_size": float(step_size),
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_dual_fista stopped without satisfying the "
            f"KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, final_gram, final_connectivity, info


# Descriptive alias for users who think in connectivity rather than natural parameters.
fit_gaussian_noise_connectivity_fista = fit_gaussian_noise_dual_fista
