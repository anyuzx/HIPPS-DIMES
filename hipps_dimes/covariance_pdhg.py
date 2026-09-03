"""Variance-whitened first-order solvers for the noise-aware COV objective.

The solver minimizes

    -3/2 logdet'(B) + 1/2 sum_a (D_a(B) - d_a)^2 / v_a

on centered Gram matrices that are positive definite on the internal subspace.
Writing ``A = V^-1/2 D`` and ``b = V^-1/2 d`` gives a dual problem with unit
quadratic curvature for every observed pair. Runtime KKT diagnostics follow
from proximal optimality without reconstructing ``B^-1`` at every iteration.
Each PDHG iteration has no inner linear solve and cannot leave the valid
covariance cone.

The module also contains the low-level direct-Gram monotone FISTA refinement
used by the public hybrid solver. It deliberately remains below the public
``HippsDimes`` and CLI surfaces: users select ``hybrid`` rather than invoking a
standalone FISTA optimizer.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import numerics as _num

_ENTROPY_COEFFICIENT = 1.5


def _scalar(value) -> float:
    return float(value.item()) if hasattr(value, "item") else float(value)


def _norm(value, xp) -> float:
    return _scalar(xp.sqrt(xp.sum(value * value)))


def _symmetrize(value):
    return 0.5 * (value + value.T)


def _center(value, xp):
    """Apply J value J without dense matrix multiplications."""
    row_mean = xp.mean(value, axis=1, keepdims=True)
    centered = value - row_mean - row_mean.T + xp.mean(value)
    return _symmetrize(centered)


def _householder_vector(n: int, xp):
    """Return w such that H=I-2ww.T maps e_last to 1/sqrt(n)."""
    unit_last = xp.zeros(n, dtype=xp.float64)
    unit_last[-1] = 1.0
    translation = xp.ones(n, dtype=xp.float64) / np.sqrt(n)
    vector = unit_last - translation
    return vector / xp.sqrt(xp.sum(vector * vector))


def _householder_congruence(matrix, vector, xp):
    """Return H matrix H for H=I-2vv.T in O(N^2) work."""
    matrix_vector = matrix @ vector
    scalar = xp.sum(vector * matrix_vector)
    return _symmetrize(
        matrix
        - 2.0 * xp.outer(vector, matrix_vector)
        - 2.0 * xp.outer(matrix_vector, vector)
        + 4.0 * scalar * xp.outer(vector, vector)
    )


def _internal_block(centered_matrix, householder_vector, xp):
    transformed = _householder_congruence(centered_matrix, householder_vector, xp)
    return _symmetrize(transformed[:-1, :-1])


def _full_centered_from_internal(internal_matrix, householder_vector, xp):
    n = internal_matrix.shape[0] + 1
    transformed = xp.zeros((n, n), dtype=xp.float64)
    transformed[:-1, :-1] = internal_matrix
    return _center(_householder_congruence(transformed, householder_vector, xp), xp)


def _internal_eigh(centered_matrix, householder_vector, xp):
    internal = _internal_block(centered_matrix, householder_vector, xp)
    if xp is np:
        eigenvalues, eigenvectors = np.linalg.eigh(internal)
    else:
        with _num.cupyx.errstate(linalg="raise"):
            eigenvalues, eigenvectors = xp.linalg.eigh(internal)
    return eigenvalues, eigenvectors


def _internal_inverse_and_logdet(centered_matrix, householder_vector, xp):
    eigenvalues, eigenvectors = _internal_eigh(centered_matrix, householder_vector, xp)
    minimum = _scalar(eigenvalues[0])
    if not np.isfinite(minimum) or minimum <= 0.0:
        raise RuntimeError(
            "PDHG requires a Gram matrix that is positive definite on the "
            "internal subspace"
        )
    inverse_internal = (eigenvectors * (1.0 / eigenvalues)) @ eigenvectors.T
    inverse_full = _full_centered_from_internal(
        inverse_internal, householder_vector, xp
    )
    logdet = _scalar(xp.sum(xp.log(eigenvalues)))
    return inverse_full, logdet, eigenvalues


def _distance_operator(gram, pair_i, pair_j, xp):
    diagonal = xp.diag(gram)
    return diagonal[pair_i] + diagonal[pair_j] - 2.0 * gram[pair_i, pair_j]


def _distance_adjoint(values, n, pair_i, pair_j, xp):
    """Adjoint of the unique-pair squared-distance operator."""
    pair_matrix = xp.zeros((n, n), dtype=xp.float64)
    pair_matrix[pair_i, pair_j] = values
    pair_matrix[pair_j, pair_i] = values
    row_sum = xp.sum(pair_matrix, axis=1)
    result = -pair_matrix
    result[xp.diag_indices(n)] = row_sum
    return _symmetrize(result)


def _eliminated_kkt_residuals(
    gram,
    inverse_gram,
    target,
    variance,
    n,
    pair_i,
    pair_j,
    xp,
    *,
    fitted=None,
):
    """Evaluate primal stationarity after analytically eliminating the dual."""
    if fitted is None:
        fitted = _distance_operator(gram, pair_i, pair_j, xp)
    eliminated_dual = (fitted - target) / variance
    data_gradient = _distance_adjoint(eliminated_dual, n, pair_i, pair_j, xp)
    entropy_gradient = _ENTROPY_COEFFICIENT * inverse_gram
    residual = data_gradient - entropy_gradient
    residual_norm = _norm(residual, xp)
    residual_scale = max(
        _norm(data_gradient, xp),
        _norm(entropy_gradient, xp),
        np.finfo(np.float64).tiny,
    )
    return (
        residual_norm,
        residual_scale,
        residual_norm / residual_scale,
        fitted,
        residual,
    )


def _independent_eliminated_kkt_residuals(
    gram,
    target,
    variance,
    pair_i,
    pair_j,
):
    """Recompute eliminated KKT stationarity from a returned physical Gram."""
    gram = _symmetrize(np.asarray(gram, dtype=np.float64))
    target = np.asarray(target, dtype=np.float64)
    variance = np.asarray(variance, dtype=np.float64)
    pair_i = np.asarray(pair_i, dtype=np.int64)
    pair_j = np.asarray(pair_j, dtype=np.int64)
    n = gram.shape[0]
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    zero_index = int(np.argmin(np.abs(eigenvalues)))
    internal_mask = np.ones(n, dtype=bool)
    internal_mask[zero_index] = False
    internal_eigenvalues = eigenvalues[internal_mask]
    minimum = float(np.min(internal_eigenvalues))
    if not np.isfinite(minimum) or minimum <= 0.0:
        raise RuntimeError(
            "Cannot certify a returned Gram matrix that is not positive "
            "definite on the internal subspace"
        )
    internal_vectors = eigenvectors[:, internal_mask]
    inverse_gram = (
        internal_vectors * (1.0 / internal_eigenvalues)
    ) @ internal_vectors.T
    return _eliminated_kkt_residuals(
        gram,
        inverse_gram,
        target,
        variance,
        n,
        pair_i,
        pair_j,
        np,
    )


def _connectivity_from_scaled_inverse(inverse_gram, distance_scale, xp):
    connectivity = -(3.0 / distance_scale) * inverse_gram
    connectivity = _symmetrize(connectivity)
    return _num.a2a(connectivity)


def _weighted_distance_operator(gram, inverse_standard_deviation, pair_i, pair_j, xp):
    return _distance_operator(gram, pair_i, pair_j, xp) * inverse_standard_deviation


def _weighted_distance_adjoint(
    values, inverse_standard_deviation, n, pair_i, pair_j, xp
):
    return _distance_adjoint(
        values * inverse_standard_deviation,
        n,
        pair_i,
        pair_j,
        xp,
    )


def _prox_centered_negative_logdet_inverse_free(
    matrix, step_size: float, householder_vector, xp
):
    """Apply the log-determinant proximal map without reconstructing B^-1."""
    matrix = _center(matrix, xp)
    eigenvalues, eigenvectors = _internal_eigh(matrix, householder_vector, xp)
    root = xp.hypot(
        eigenvalues,
        np.sqrt(4.0 * _ENTROPY_COEFFICIENT * step_size),
    )
    nonnegative = xp.maximum(eigenvalues, 0.0)
    nonpositive = xp.minimum(eigenvalues, 0.0)
    updated = xp.where(
        eigenvalues >= 0.0,
        0.5 * (nonnegative + root),
        (2.0 * _ENTROPY_COEFFICIENT * step_size) / (root - nonpositive),
    )
    internal = (eigenvectors * updated) @ eigenvectors.T
    gram = _full_centered_from_internal(internal, householder_vector, xp)
    logdet = _scalar(xp.sum(xp.log(updated)))
    return gram, logdet, updated, eigenvectors


def _inverse_from_internal_spectrum(eigenvalues, eigenvectors, householder_vector, xp):
    inverse_internal = (eigenvectors * (1.0 / eigenvalues)) @ eigenvectors.T
    return _full_centered_from_internal(inverse_internal, householder_vector, xp)


def _certify_weighted_operator_norm_squared(
    inverse_standard_deviation,
    n,
    pair_i,
    pair_j,
    xp,
    *,
    max_iterations,
    relative_tolerance,
    safety_factor,
    phase,
    progress_callback=None,
):
    """Bound ``||V^-1/2 D||^2`` by edge-space Collatz iteration.

    The edge-space normal operator ``G = A A*`` is entrywise nonnegative.
    Starting from a strictly positive vector makes the maximum component of
    ``Gx / x`` a certified upper bound on its Perron eigenvalue.  This avoids
    the dominant-eigenspace failure possible for an arbitrary Gram-space
    power-iteration start while retaining matrix-free O(E + N^2) updates.
    """
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("operator_norm_power_iterations must be positive")
    if not np.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise ValueError("operator_norm_tolerance must be positive and finite")
    if not np.isfinite(safety_factor) or safety_factor <= 1.0:
        raise ValueError("operator_norm_safety must exceed one")

    vector = xp.ones(len(pair_i), dtype=xp.float64)
    vector = vector / _norm(vector, xp)
    lower_bound = 0.0
    upper_bound = np.inf
    first_row_sum_bound = None
    relative_bound_gap = np.inf
    bound_tolerance_reached = False
    used_row_sum_fallback = False
    iteration = 0
    start_time = time.perf_counter()
    for iteration in range(1, int(max_iterations) + 1):
        adjoint = _weighted_distance_adjoint(
            vector,
            inverse_standard_deviation,
            n,
            pair_i,
            pair_j,
            xp,
        )
        next_vector = _weighted_distance_operator(
            adjoint,
            inverse_standard_deviation,
            pair_i,
            pair_j,
            xp,
        )
        ratios = next_vector / vector
        valid_ratios = bool(_scalar(xp.all(xp.isfinite(ratios) & (ratios > 0.0))))
        if not valid_ratios:
            if first_row_sum_bound is None:
                raise RuntimeError(
                    "weighted distance operator produced invalid Collatz ratios"
                )
            upper_bound = first_row_sum_bound
            lower_bound = 0.0
            relative_bound_gap = 1.0
            used_row_sum_fallback = True
            break

        lower_bound = _scalar(xp.min(ratios))
        upper_bound = _scalar(xp.max(ratios))
        if first_row_sum_bound is None:
            first_row_sum_bound = upper_bound
        relative_bound_gap = (upper_bound - lower_bound) / max(
            upper_bound, np.finfo(np.float64).tiny
        )
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_operator_norm",
                    "phase": phase,
                    "iteration": iteration,
                    "total": int(max_iterations),
                    "operator_norm_lower_bound": float(lower_bound),
                    "operator_norm_upper_bound": float(upper_bound),
                    "operator_norm_relative_bound_gap": float(relative_bound_gap),
                    "method": "COV-FISTA" if phase == "fista" else "COV-PDHG",
                    "general_method": "covariance_optimization",
                    "use_gpu": bool(xp is not np),
                    "noisy": True,
                }
            )
        if relative_bound_gap <= relative_tolerance:
            bound_tolerance_reached = True
            break
        next_norm = _norm(next_vector, xp)
        if not np.isfinite(next_norm) or next_norm <= 0.0:
            upper_bound = first_row_sum_bound
            lower_bound = 0.0
            relative_bound_gap = 1.0
            used_row_sum_fallback = True
            break
        vector = next_vector / next_norm

    safe_bound = safety_factor * upper_bound
    return safe_bound, {
        "method": "edge_space_collatz",
        "iterations": int(iteration),
        "bound_tolerance_reached": bool(bound_tolerance_reached),
        "lower_bound": float(lower_bound),
        "certified_upper_bound": float(upper_bound),
        "relative_bound_gap": float(relative_bound_gap),
        "first_row_sum_upper_bound": float(first_row_sum_bound),
        "used_row_sum_fallback": bool(used_row_sum_fallback),
        "safety_factor": float(safety_factor),
        "safe_operator_norm_squared": float(safe_bound),
        "wall_seconds": time.perf_counter() - start_time,
    }


def _choose_initial_whitened_dual(
    gram,
    inverse_gram,
    whitened_target,
    inverse_standard_deviation,
    standard_deviation,
    pair_i,
    pair_j,
    xp,
    mode,
):
    fitted = _weighted_distance_operator(
        gram,
        inverse_standard_deviation,
        pair_i,
        pair_j,
        xp,
    )
    if mode == "zero":
        dual = xp.zeros_like(whitened_target)
    elif mode == "residual":
        dual = fitted - whitened_target
    elif mode == "connectivity":
        dual = standard_deviation * (
            -_ENTROPY_COEFFICIENT * inverse_gram[pair_i, pair_j]
        )
    else:
        raise ValueError(
            "dual_initialization must be 'zero', 'residual', or 'connectivity'"
        )
    return dual.copy(), mode


def _inverse_free_residuals(
    previous_gram,
    gram,
    weighted_dual_adjoint,
    dual,
    whitened_fitted,
    whitened_target,
    inverse_standard_deviation,
    primal_step,
    n,
    pair_i,
    pair_j,
    xp,
):
    """Compute split and eliminated KKT residuals without forming B^-1."""
    gram_rate = (gram - previous_gram) / primal_step
    # Proximal optimality gives 1.5 B^-1 = A* u + gram_rate.
    entropy_gradient = weighted_dual_adjoint + gram_rate
    primal_residual = -gram_rate
    primal_norm = _norm(primal_residual, xp)
    primal_scale = max(
        _norm(weighted_dual_adjoint, xp),
        _norm(entropy_gradient, xp),
        np.finfo(np.float64).tiny,
    )

    dual_residual = whitened_fitted - whitened_target - dual
    dual_norm = _norm(dual_residual, xp)
    fitted_residual = whitened_fitted - whitened_target
    dual_scale = max(
        _norm(fitted_residual, xp),
        _norm(dual, xp),
        np.finfo(np.float64).tiny,
    )
    dual_component = _weighted_distance_adjoint(
        dual_residual,
        inverse_standard_deviation,
        n,
        pair_i,
        pair_j,
        xp,
    )
    data_gradient = weighted_dual_adjoint + dual_component
    eliminated_residual = dual_component - gram_rate
    eliminated_norm = _norm(eliminated_residual, xp)
    eliminated_scale = max(
        _norm(data_gradient, xp),
        _norm(entropy_gradient, xp),
        np.finfo(np.float64).tiny,
    )
    return {
        "primal_residual": primal_residual,
        "primal_norm": primal_norm,
        "primal_scale": primal_scale,
        "primal_relative": primal_norm / primal_scale,
        "dual_residual": dual_residual,
        "dual_norm": dual_norm,
        "dual_scale": dual_scale,
        "dual_relative": dual_norm / dual_scale,
        "dual_component": dual_component,
        "dual_component_norm": _norm(dual_component, xp),
        "eliminated_residual": eliminated_residual,
        "eliminated_norm": eliminated_norm,
        "eliminated_scale": eliminated_scale,
        "eliminated_relative": eliminated_norm / eliminated_scale,
    }


def _connectivity_checkpoint(
    eigenvalues,
    eigenvectors,
    householder,
    distance_scale,
    xp,
):
    inverse_gram = _inverse_from_internal_spectrum(
        eigenvalues, eigenvectors, householder, xp
    )
    connectivity = _connectivity_from_scaled_inverse(inverse_gram, distance_scale, xp)
    if xp is np:
        return np.asarray(connectivity).copy(), inverse_gram
    return _num.cp.asnumpy(connectivity), inverse_gram


def fit_gaussian_noise_covariance_fista(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_gram,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    backtracking_factor=0.5,
    max_backtracking_steps=80,
    distance_scale=None,
    operator_norm_power_iterations=25,
    operator_norm_tolerance=1e-4,
    operator_norm_safety=1.10,
    _precomputed_operator_norm=None,
    save_steps=None,
    progress_callback=None,
):
    """Refine a physical centered Gram matrix with monotone FISTA.

    This low-level solver minimizes the same Gaussian COV
    objective as :func:`fit_gaussian_noise_covariance_pdhg`. ``initial_gram``
    is used directly, without a connectivity round trip or scalar calibration.
    The function is intentionally not re-exported through the package API or
    CLI; it is the refinement implementation behind the public hybrid solver.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "FISTA GPU optimization was requested, but CuPy and an "
            "accessible CUDA GPU are not available"
        )
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
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(backtracking_factor) or not (0.0 < backtracking_factor < 1.0):
        raise ValueError("backtracking_factor must lie strictly between zero and one")
    if (
        not isinstance(max_backtracking_steps, (int, np.integer))
        or max_backtracking_steps <= 0
    ):
        raise ValueError("max_backtracking_steps must be a positive integer")

    n = observed.shape[0]
    n_modes = n - 1
    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer)) or step < 1 or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")
    physical_initial_gram = np.asarray(initial_gram, dtype=np.float64)
    if physical_initial_gram.shape != (n, n):
        raise ValueError("initial_gram must have the same square shape as the target")
    if not np.all(np.isfinite(physical_initial_gram)):
        raise ValueError("initial_gram must contain only finite values")
    gram_norm = max(
        float(np.linalg.norm(physical_initial_gram, ord="fro")),
        np.finfo(np.float64).tiny,
    )
    symmetry_residual = float(
        np.linalg.norm(physical_initial_gram - physical_initial_gram.T, ord="fro")
        / gram_norm
    )
    centering_residual = float(
        np.linalg.norm(np.sum(physical_initial_gram, axis=1)) / gram_norm
    )
    if symmetry_residual > 1e-10:
        raise ValueError("initial_gram must be symmetric to numerical precision")
    if centering_residual > 1e-10:
        raise ValueError("initial_gram must be centered to numerical precision")
    physical_initial_gram = _center(_symmetrize(physical_initial_gram), np)

    if distance_scale is None:
        distance_scale = float(np.median(target_pairs))
    if not np.isfinite(distance_scale) or distance_scale <= 0.0:
        raise ValueError("distance_scale must be positive and finite")

    solver_target_cpu = target_pairs / distance_scale
    solver_variance_cpu = pair_variance / (distance_scale * distance_scale)
    inverse_standard_deviation_cpu = 1.0 / np.sqrt(solver_variance_cpu)
    whitened_target_cpu = solver_target_cpu * inverse_standard_deviation_cpu
    gram_cpu = physical_initial_gram / distance_scale
    if use_gpu:
        xp = _num.cp
        solver_pair_i = xp.asarray(pair_i, dtype=xp.int32)
        solver_pair_j = xp.asarray(pair_j, dtype=xp.int32)
        solver_target = xp.asarray(solver_target_cpu, dtype=xp.float64)
        solver_variance = xp.asarray(solver_variance_cpu, dtype=xp.float64)
        inverse_standard_deviation = xp.asarray(
            inverse_standard_deviation_cpu, dtype=xp.float64
        )
        whitened_target = xp.asarray(whitened_target_cpu, dtype=xp.float64)
        gram = xp.asarray(gram_cpu, dtype=xp.float64)
    else:
        xp = np
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target = solver_target_cpu
        solver_variance = solver_variance_cpu
        inverse_standard_deviation = inverse_standard_deviation_cpu
        whitened_target = whitened_target_cpu
        gram = gram_cpu

    householder = _householder_vector(n, xp)
    inverse_gram, logdet, eigenvalues = _internal_inverse_and_logdet(
        gram, householder, xp
    )
    if _precomputed_operator_norm is None:
        operator_norm_squared, operator_norm_info = (
            _certify_weighted_operator_norm_squared(
                inverse_standard_deviation,
                n,
                solver_pair_i,
                solver_pair_j,
                xp,
                max_iterations=operator_norm_power_iterations,
                relative_tolerance=operator_norm_tolerance,
                safety_factor=operator_norm_safety,
                phase="fista",
                progress_callback=progress_callback,
            )
        )
    else:
        operator_norm_squared, operator_norm_info = _precomputed_operator_norm
        operator_norm_squared = float(operator_norm_squared)
        if not np.isfinite(operator_norm_squared) or operator_norm_squared <= 0.0:
            raise ValueError("precomputed operator norm must be positive and finite")
        operator_norm_info = dict(operator_norm_info)
    step_size = 1.0 / operator_norm_squared
    initial_step_size = step_size
    logdet_offset = n_modes * np.log(distance_scale)

    def evaluate(
        state_gram,
        state_inverse,
        state_logdet,
        state_eigenvalues,
        whitened_fitted=None,
    ):
        if whitened_fitted is None:
            whitened_fitted = _weighted_distance_operator(
                state_gram,
                inverse_standard_deviation,
                solver_pair_i,
                solver_pair_j,
                xp,
            )
        standardized_residual = whitened_fitted - whitened_target
        data_objective = 0.5 * _scalar(
            xp.sum(standardized_residual * standardized_residual)
        )
        fitted = whitened_fitted / inverse_standard_deviation
        data_residual = fitted - solver_target
        relative_loss = _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual / solver_target)))
        )
        weighted_rmse = _scalar(xp.sqrt(xp.mean(xp.square(standardized_residual))))
        distance_rmse = distance_scale * _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual)))
        )
        kkt = _eliminated_kkt_residuals(
            state_gram,
            state_inverse,
            solver_target,
            solver_variance,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
            fitted=fitted,
        )
        minimum_eigenvalue = _scalar(state_eigenvalues[0])
        maximum_eigenvalue = _scalar(state_eigenvalues[-1])
        negative_entropy = -_ENTROPY_COEFFICIENT * (state_logdet + logdet_offset)
        entropy = state_logdet + logdet_offset - n_modes * np.log(3.0)
        return {
            "objective": negative_entropy + data_objective,
            "negative_entropy_objective": negative_entropy,
            "data_objective": data_objective,
            "loss": relative_loss,
            "entropy": entropy,
            "weighted_rmse": weighted_rmse,
            "distance_rmse": distance_rmse,
            "relative_eliminated_kkt_residual": kkt[2],
            "stationarity_residual_norm": kkt[0],
            "stationarity_residual_scale": kkt[1],
            "minimum_internal_gram_eigenvalue": (distance_scale * minimum_eigenvalue),
            "maximum_internal_gram_eigenvalue": (distance_scale * maximum_eigenvalue),
            "gram_condition_number": maximum_eigenvalue / minimum_eigenvalue,
        }

    current = evaluate(gram, inverse_gram, logdet, eigenvalues)
    initial_objective = float(current["objective"])
    initial_relative_kkt = float(current["relative_eliminated_kkt_residual"])
    history = {
        "iteration": [],
        "objective": [],
        "negative_entropy_objective": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "relative_eliminated_kkt_residual": [],
        "relative_gradient_norm": [],
        "relative_stationarity_residual": [],
        "stationarity_residual_norm": [],
        "stationarity_residual_scale": [],
        "relative_step": [],
        "step_size": [],
        "backtracking_steps": [],
        "momentum_parameter": [],
        "momentum_coefficient": [],
        "restarted": [],
        "objective_decrease": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "gram_condition_number": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}
    search_gram = gram.copy()
    search_is_extrapolated = False
    momentum_parameter = 1.0
    momentum_coefficient = 0.0
    restart_count = 0
    backtracking_reductions = 0
    converged = current["stationarity_residual_norm"] <= (
        absolute_tolerance + relative_tolerance * current["stationarity_residual_scale"]
    )
    status = "optimality_tolerance" if converged else "max_iterations"
    message = (
        "KKT tolerance reached at the provided Gram"
        if converged
        else "maximum number of monotone FISTA iterations reached"
    )
    start_time = time.perf_counter()

    for iteration in range(1, max_iterations + 1):
        if converged:
            break
        accepted = False
        restarted = False
        iteration_backtracking = 0
        used_momentum_parameter = momentum_parameter
        used_momentum_coefficient = momentum_coefficient

        for _ in range(max_backtracking_steps):
            search_fitted = _weighted_distance_operator(
                search_gram,
                inverse_standard_deviation,
                solver_pair_i,
                solver_pair_j,
                xp,
            )
            search_residual = search_fitted - whitened_target
            search_data_objective = 0.5 * _scalar(
                xp.sum(search_residual * search_residual)
            )
            search_data_gradient = _weighted_distance_adjoint(
                search_residual,
                inverse_standard_deviation,
                n,
                solver_pair_i,
                solver_pair_j,
                xp,
            )
            (
                candidate,
                candidate_logdet,
                candidate_eigenvalues,
                candidate_eigenvectors,
            ) = _prox_centered_negative_logdet_inverse_free(
                search_gram - step_size * search_data_gradient,
                step_size,
                householder,
                xp,
            )
            candidate_fitted = _weighted_distance_operator(
                candidate,
                inverse_standard_deviation,
                solver_pair_i,
                solver_pair_j,
                xp,
            )
            candidate_residual = candidate_fitted - whitened_target
            candidate_data_objective = 0.5 * _scalar(
                xp.sum(candidate_residual * candidate_residual)
            )
            delta = candidate - search_gram
            quadratic_upper_bound = (
                search_data_objective
                + _scalar(xp.sum(search_data_gradient * delta))
                + 0.5 * _scalar(xp.sum(delta * delta)) / step_size
            )
            majorization_tolerance = 1e-12 * max(
                1.0,
                abs(search_data_objective),
                abs(candidate_data_objective),
            )
            if candidate_data_objective > (
                quadratic_upper_bound + majorization_tolerance
            ):
                step_size *= backtracking_factor
                iteration_backtracking += 1
                backtracking_reductions += 1
                continue

            candidate_objective = (
                -_ENTROPY_COEFFICIENT * (candidate_logdet + logdet_offset)
                + candidate_data_objective
            )
            monotonicity_tolerance = 1e-12 * max(
                1.0, abs(current["objective"]), abs(candidate_objective)
            )
            if candidate_objective > (current["objective"] + monotonicity_tolerance):
                if search_is_extrapolated:
                    restart_count += 1
                    restarted = True
                    search_gram = gram.copy()
                    search_is_extrapolated = False
                    momentum_parameter = 1.0
                    momentum_coefficient = 0.0
                    used_momentum_parameter = 1.0
                    used_momentum_coefficient = 0.0
                else:
                    step_size *= backtracking_factor
                    iteration_backtracking += 1
                    backtracking_reductions += 1
                continue
            accepted = True
            break

        if not accepted:
            status = "backtracking_failed"
            message = "monotone FISTA backtracking failed"
            break

        candidate_inverse = _inverse_from_internal_spectrum(
            candidate_eigenvalues,
            candidate_eigenvectors,
            householder,
            xp,
        )
        candidate_metric = evaluate(
            candidate,
            candidate_inverse,
            candidate_logdet,
            candidate_eigenvalues,
            candidate_fitted,
        )
        previous_gram = gram
        objective_decrease = current["objective"] - candidate_metric["objective"]
        relative_step = _norm(candidate - previous_gram, xp) / max(
            _norm(previous_gram, xp), np.finfo(np.float64).tiny
        )
        gram = candidate
        inverse_gram = candidate_inverse
        current = candidate_metric

        history["iteration"].append(iteration)
        for key in (
            "objective",
            "negative_entropy_objective",
            "data_objective",
            "loss",
            "entropy",
            "weighted_rmse",
            "distance_rmse",
            "relative_eliminated_kkt_residual",
            "stationarity_residual_norm",
            "stationarity_residual_scale",
            "minimum_internal_gram_eigenvalue",
            "maximum_internal_gram_eigenvalue",
            "gram_condition_number",
        ):
            history[key].append(current[key])
        history["relative_gradient_norm"].append(
            current["relative_eliminated_kkt_residual"]
        )
        history["relative_stationarity_residual"].append(
            current["relative_eliminated_kkt_residual"]
        )
        history["relative_step"].append(relative_step)
        history["step_size"].append(step_size)
        history["backtracking_steps"].append(iteration_backtracking)
        history["momentum_parameter"].append(used_momentum_parameter)
        history["momentum_coefficient"].append(used_momentum_coefficient)
        history["restarted"].append(restarted)
        history["objective_decrease"].append(max(objective_decrease, 0.0))
        history["connectivity_offdiagonal_l2"].append(np.nan)

        if iteration in save_steps_set:
            checkpoint = _connectivity_from_scaled_inverse(
                candidate_inverse, distance_scale, xp
            )
            if xp is np:
                checkpoint = np.asarray(checkpoint).copy()
            else:
                checkpoint = _num.cp.asnumpy(checkpoint)
            connectivity_at_steps[iteration] = checkpoint
            offdiagonal = checkpoint[np.triu_indices(n, k=1)]
            history["connectivity_offdiagonal_l2"][-1] = float(
                np.linalg.norm(offdiagonal)
            )

        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": "fista",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": current["objective"],
                    "relative_gradient_norm": current[
                        "relative_eliminated_kkt_residual"
                    ],
                    "relative_eliminated_kkt_residual": current[
                        "relative_eliminated_kkt_residual"
                    ],
                    "relative_stationarity_residual": current[
                        "relative_eliminated_kkt_residual"
                    ],
                    "step_size": step_size,
                    "backtracking_steps": iteration_backtracking,
                    "restarted": restarted,
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV-FISTA",
                    "general_method": "covariance_optimization",
                }
            )

        converged = current["stationarity_residual_norm"] <= (
            absolute_tolerance
            + relative_tolerance * current["stationarity_residual_scale"]
        )
        if converged:
            status = "optimality_tolerance"
            message = "monotone FISTA KKT tolerance reached"
            break

        next_momentum_parameter = 0.5 * (
            1.0 + np.sqrt(1.0 + 4.0 * momentum_parameter**2)
        )
        next_momentum_coefficient = (momentum_parameter - 1.0) / next_momentum_parameter
        search_gram = _center(
            gram + next_momentum_coefficient * (gram - previous_gram), xp
        )
        search_is_extrapolated = next_momentum_coefficient > 0.0
        momentum_parameter = next_momentum_parameter
        momentum_coefficient = next_momentum_coefficient

    optimization_wall_seconds = time.perf_counter() - start_time
    full_gram = distance_scale * gram
    full_fitted = _num._squared_distances_from_gram(full_gram, array_module=xp)
    connectivity = _connectivity_from_scaled_inverse(inverse_gram, distance_scale, xp)
    if use_gpu:
        full_fitted = _num.cp.asnumpy(full_fitted)
        full_gram = _num.cp.asnumpy(full_gram)
        connectivity = _num.cp.asnumpy(connectivity)
    else:
        full_fitted = np.asarray(full_fitted)
        full_gram = np.asarray(full_gram)
        connectivity = np.asarray(connectivity)

    independent = _independent_eliminated_kkt_residuals(
        full_gram, target_pairs, pair_variance, pair_i, pair_j
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
    runtime_converged = current["stationarity_residual_norm"] <= (
        absolute_tolerance + relative_tolerance * current["stationarity_residual_scale"]
    )
    converged = bool(runtime_converged and independent_converged)
    if converged:
        status = "optimality_tolerance"
        message = "monotone FISTA and independent KKT tolerance reached"
    elif status == "optimality_tolerance":
        status = "independent_kkt_failed"
        message = (
            "runtime FISTA KKT tolerance was reached, but the returned Gram "
            "failed the independent certificate"
        )

    fitted_pairs = full_fitted[pair_i, pair_j]
    pair_stationarity = (
        fitted_pairs - target_pairs - 0.5 * pair_variance * connectivity[pair_i, pair_j]
    )
    pair_stationarity_norm = float(np.linalg.norm(pair_stationarity))
    pair_stationarity_scale = max(
        1.0,
        float(np.linalg.norm(fitted_pairs - target_pairs)),
        float(np.linalg.norm(0.5 * pair_variance * connectivity[pair_i, pair_j])),
    )
    for key, values in history.items():
        if key in {"iteration", "backtracking_steps"}:
            history[key] = np.asarray(values, dtype=np.int64)
        elif key == "restarted":
            history[key] = np.asarray(values, dtype=bool)
        else:
            history[key] = np.asarray(values, dtype=np.float64)
    iterations = int(history["iteration"][-1]) if len(history["iteration"]) else 0
    history["is_returned_iterate"] = np.zeros(iterations, dtype=bool)
    if iterations:
        history["is_returned_iterate"][-1] = True
    info = {
        "converged": converged,
        "status": status,
        "message": message,
        "iterations": iterations,
        "returned_iteration": iterations,
        "objective": float(current["objective"]),
        "loss": float(current["loss"]),
        "entropy": float(current["entropy"]),
        "initial_objective": initial_objective,
        "initial_relative_eliminated_kkt_residual": initial_relative_kkt,
        "relative_gradient_norm": float(independent_relative),
        "relative_eliminated_kkt_residual": float(independent_relative),
        "runtime_relative_eliminated_kkt_residual": float(
            current["relative_eliminated_kkt_residual"]
        ),
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
            np.max(np.abs(pair_stationarity), initial=0.0)
        ),
        "termination_internal_kkt_converged": bool(runtime_converged),
        "independent_kkt_converged": bool(independent_converged),
        "independent_kkt_recomputed_from_returned_gram": True,
        "relative_tolerance": float(relative_tolerance),
        "absolute_tolerance": float(absolute_tolerance),
        "observed_pair_count": len(target_pairs),
        "noise_model": noise_model,
        "noise_parameter": noise_parameter,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": {
            "kind": "provided_gram",
            "scalar_calibration": None,
            "physical_gram_used_directly": True,
            "input_symmetry_relative_residual": symmetry_residual,
            "input_centering_relative_residual": centering_residual,
        },
        "algorithm": "direct_b_monotone_fista",
        "coordinate_parameterization": "centered_householder",
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "cupy_version": _num.cp.__version__ if use_gpu else None,
        "distance_scale": float(distance_scale),
        "weighted_operator_norm": operator_norm_info,
        "initial_step_size": float(initial_step_size),
        "final_step_size": float(step_size),
        "backtracking_factor": float(backtracking_factor),
        "maximum_backtracking_steps": int(max_backtracking_steps),
        "backtracking_reductions": int(backtracking_reductions),
        "restart_count": int(restart_count),
        "accelerated": True,
        "monotone_restart": True,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        # Initial evaluation, one reconstruction per accepted iterate, and
        # the independent certificate recomputed from the returned Gram.
        "inverse_reconstruction_count": iterations + 2,
        "inverse_reconstructed_each_iteration": True,
        "wall_seconds": optimization_wall_seconds,
        "objective_definition": (
            "-1.5*logdet(B_internal) + " "0.5*||V^(-1/2)*(D(B)-D_obs)||^2"
        ),
        "stationarity_definition": (
            "D_adjoint((D(B)-D_obs)/noise_variance)-1.5*B_pseudoinverse"
        ),
        "pair_stationarity_definition": ("D_fit_ij-D_obs_ij-noise_variance_ij*A_ij/2"),
        "logged_metric_timing": "post-update accepted monotone FISTA iterate",
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_fista stopped without satisfying "
            f"the independent KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return full_fitted, full_gram, connectivity, info


def fit_gaussian_noise_covariance_pdhg(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    save_steps=None,
    progress_callback=None,
    theta=1.0,
    step_safety=0.95,
    step_ratio=None,
    initial_dual_step=1.0,
    adaptive_steps=True,
    adaptation_interval=20,
    adaptation_threshold=3.0,
    adaptation_factor=np.sqrt(2.0),
    minimum_step_ratio=1e-16,
    maximum_step_ratio=1e16,
    dual_initialization="residual",
    distance_scale=None,
    return_best=True,
    operator_norm_power_iterations=25,
    operator_norm_tolerance=1e-4,
    operator_norm_safety=1.10,
):
    """Fit the noise-aware COV objective with variance-whitened PDHG.

    The distance likelihood is standardized with
    ``A = V^-1/2 D`` and ``b = V^-1/2 d``. The dual update therefore has
    unit curvature for every observed pair. During the iterations the exact
    eliminated KKT residual is obtained from proximal optimality, so no Gram
    inverse or connectivity matrix is reconstructed unless a checkpoint or the
    final result is requested.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "PDHG GPU optimization was requested, but CuPy and an "
            "accessible CUDA GPU are not available"
        )
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
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(theta) or not (0.0 <= theta <= 1.0):
        raise ValueError("theta must lie in [0, 1]")
    if not np.isfinite(step_safety) or not (0.0 < step_safety < 1.0):
        raise ValueError("step_safety must lie strictly between zero and one")
    if not np.isfinite(initial_dual_step) or initial_dual_step <= 0.0:
        raise ValueError("initial_dual_step must be positive and finite")
    if (
        not isinstance(adaptation_interval, (int, np.integer))
        or adaptation_interval <= 0
    ):
        raise ValueError("adaptation_interval must be a positive integer")
    if not np.isfinite(adaptation_threshold) or adaptation_threshold <= 1.0:
        raise ValueError("adaptation_threshold must exceed one")
    if not np.isfinite(adaptation_factor) or adaptation_factor <= 1.0:
        raise ValueError("adaptation_factor must exceed one")
    if (
        not np.isfinite(minimum_step_ratio)
        or not np.isfinite(maximum_step_ratio)
        or minimum_step_ratio <= 0.0
        or maximum_step_ratio < minimum_step_ratio
    ):
        raise ValueError("invalid PDHG step-ratio bounds")

    n = observed.shape[0]
    n_modes = n - 1
    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer)) or step < 1 or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    basis = _num._centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    reduced_gram, initialization_info = _num._initialize_gaussian_reduced_gram(
        observed,
        pair_i,
        pair_j,
        basis,
        initial_connectivity,
    )

    if use_gpu:
        xp = _num.cp
        solver_basis = xp.asarray(basis, dtype=xp.float64)
        solver_pair_i = xp.asarray(pair_i, dtype=xp.int32)
        solver_pair_j = xp.asarray(pair_j, dtype=xp.int32)
        solver_target_pairs = xp.asarray(target_pairs, dtype=xp.float64)
        solver_inverse_variance = xp.asarray(inverse_variance, dtype=xp.float64)
        reduced_gram = xp.asarray(reduced_gram, dtype=xp.float64)
    else:
        xp = np
        solver_basis = basis
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target_pairs = target_pairs
        solver_inverse_variance = inverse_variance

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

    full_gram = _center(solver_basis @ reduced_gram @ solver_basis.T, xp)
    if distance_scale is None:
        distance_scale = float(np.median(target_pairs))
    if not np.isfinite(distance_scale) or distance_scale <= 0.0:
        raise ValueError("distance_scale must be positive and finite")

    gram = full_gram / distance_scale
    target = solver_target_pairs / distance_scale
    variance = xp.asarray(
        pair_variance / (distance_scale * distance_scale),
        dtype=xp.float64,
    )
    standard_deviation = xp.sqrt(variance)
    inverse_standard_deviation = 1.0 / standard_deviation
    whitened_target = target * inverse_standard_deviation
    householder = _householder_vector(n, xp)

    inverse_reconstruction_count = 0
    inverse_gram, logdet, internal_eigenvalues = _internal_inverse_and_logdet(
        gram, householder, xp
    )
    inverse_reconstruction_count += 1
    dual, chosen_dual_initialization = _choose_initial_whitened_dual(
        gram,
        inverse_gram,
        whitened_target,
        inverse_standard_deviation,
        standard_deviation,
        solver_pair_i,
        solver_pair_j,
        xp,
        dual_initialization,
    )

    operator_norm_squared, operator_norm_info = _certify_weighted_operator_norm_squared(
        inverse_standard_deviation,
        n,
        solver_pair_i,
        solver_pair_j,
        xp,
        max_iterations=operator_norm_power_iterations,
        relative_tolerance=operator_norm_tolerance,
        safety_factor=operator_norm_safety,
        phase="pdhg",
        progress_callback=progress_callback,
    )
    step_product = (step_safety * step_safety) / operator_norm_squared
    if step_ratio is None:
        dual_step = float(initial_dual_step)
        primal_step = step_product / dual_step
        step_ratio = primal_step / dual_step
        if step_ratio < minimum_step_ratio or step_ratio > maximum_step_ratio:
            step_ratio = float(
                np.clip(step_ratio, minimum_step_ratio, maximum_step_ratio)
            )
            primal_step = np.sqrt(step_product * step_ratio)
            dual_step = np.sqrt(step_product / step_ratio)
    else:
        step_ratio = float(step_ratio)
        if (
            not np.isfinite(step_ratio)
            or step_ratio < minimum_step_ratio
            or step_ratio > maximum_step_ratio
        ):
            raise ValueError("step_ratio lies outside its configured bounds")
        primal_step = np.sqrt(step_product * step_ratio)
        dual_step = np.sqrt(step_product / step_ratio)

    extrapolated_gram = gram.copy()
    history = {
        "iteration": [],
        "objective": [],
        "negative_entropy_objective": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "relative_primal_kkt_residual": [],
        "relative_dual_kkt_residual": [],
        "relative_eliminated_kkt_residual": [],
        "relative_gradient_norm": [],
        "relative_stationarity_residual": [],
        "primal_component_norm": [],
        "dual_component_primal_norm": [],
        "component_balance_ratio": [],
        "primal_step": [],
        "dual_step": [],
        "step_ratio": [],
        "step_adapted": [],
        "relative_primal_change": [],
        "relative_dual_change": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "gram_condition_number": [],
        "internal_precision_frobenius_norm": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}
    best = None
    stopping_criterion_reached = False
    status = "max_iterations"
    message = "maximum number of variance-whitened PDHG iterations reached"
    start_time = time.perf_counter()

    for iteration in range(1, max_iterations + 1):
        previous_gram = gram
        extrapolated_distances = _weighted_distance_operator(
            extrapolated_gram,
            inverse_standard_deviation,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        next_dual = (dual + dual_step * (extrapolated_distances - whitened_target)) / (
            1.0 + dual_step
        )
        weighted_dual_adjoint = _weighted_distance_adjoint(
            next_dual,
            inverse_standard_deviation,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        (
            next_gram,
            next_logdet,
            next_eigenvalues,
            next_eigenvectors,
        ) = _prox_centered_negative_logdet_inverse_free(
            gram - primal_step * weighted_dual_adjoint,
            primal_step,
            householder,
            xp,
        )
        next_extrapolated = next_gram + theta * (next_gram - gram)

        primal_change = _norm(next_gram - gram, xp) / max(
            _norm(gram, xp), np.finfo(np.float64).tiny
        )
        dual_change = _norm(next_dual - dual, xp) / max(_norm(dual, xp), 1.0)
        whitened_fitted = _weighted_distance_operator(
            next_gram,
            inverse_standard_deviation,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        residuals = _inverse_free_residuals(
            previous_gram,
            next_gram,
            weighted_dual_adjoint,
            next_dual,
            whitened_fitted,
            whitened_target,
            inverse_standard_deviation,
            primal_step,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        stationarity_score = residuals["eliminated_relative"]
        standardized_data_residual = whitened_fitted - whitened_target
        data_objective = 0.5 * _scalar(xp.sum(standardized_data_residual**2))
        negative_entropy = -_ENTROPY_COEFFICIENT * (
            next_logdet + n_modes * np.log(distance_scale)
        )
        objective = negative_entropy + data_objective
        fitted_normalized = whitened_fitted * standard_deviation
        data_residual = fitted_normalized - target
        relative_loss = _scalar(xp.sqrt(xp.mean(xp.square(data_residual / target))))
        weighted_rmse = _scalar(xp.sqrt(xp.mean(standardized_data_residual**2)))
        distance_rmse = distance_scale * _scalar(xp.sqrt(xp.mean(data_residual**2)))
        entropy = next_logdet + n_modes * np.log(distance_scale) - n_modes * np.log(3.0)
        minimum_eigenvalue = _scalar(next_eigenvalues[0])
        maximum_eigenvalue = _scalar(next_eigenvalues[-1])
        condition_number = maximum_eigenvalue / minimum_eigenvalue
        precision_frobenius = (
            3.0 / distance_scale * _scalar(xp.sqrt(xp.sum(1.0 / next_eigenvalues**2)))
        )

        used_primal_step = primal_step
        used_dual_step = dual_step
        used_step_ratio = step_ratio
        primal_component_norm = residuals["primal_norm"]
        dual_component_norm = residuals["dual_component_norm"]
        component_balance_ratio = dual_component_norm / max(
            primal_component_norm, np.finfo(np.float64).tiny
        )
        adapted = False

        gram = next_gram
        logdet = next_logdet
        internal_eigenvalues = next_eigenvalues
        dual = next_dual
        extrapolated_gram = next_extrapolated

        history["iteration"].append(iteration)
        history["objective"].append(objective)
        history["negative_entropy_objective"].append(negative_entropy)
        history["data_objective"].append(data_objective)
        history["loss"].append(relative_loss)
        history["entropy"].append(entropy)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["relative_primal_kkt_residual"].append(residuals["primal_relative"])
        history["relative_dual_kkt_residual"].append(residuals["dual_relative"])
        history["relative_eliminated_kkt_residual"].append(stationarity_score)
        history["relative_gradient_norm"].append(stationarity_score)
        history["relative_stationarity_residual"].append(stationarity_score)
        history["primal_component_norm"].append(primal_component_norm)
        history["dual_component_primal_norm"].append(dual_component_norm)
        history["component_balance_ratio"].append(component_balance_ratio)
        history["primal_step"].append(used_primal_step)
        history["dual_step"].append(used_dual_step)
        history["step_ratio"].append(used_step_ratio)
        history["step_adapted"].append(adapted)
        history["relative_primal_change"].append(primal_change)
        history["relative_dual_change"].append(dual_change)
        history["minimum_internal_gram_eigenvalue"].append(
            distance_scale * minimum_eigenvalue
        )
        history["maximum_internal_gram_eigenvalue"].append(
            distance_scale * maximum_eigenvalue
        )
        history["gram_condition_number"].append(condition_number)
        history["internal_precision_frobenius_norm"].append(precision_frobenius)
        history["connectivity_offdiagonal_l2"].append(np.nan)

        if adaptive_steps and iteration % adaptation_interval == 0:
            proposed_ratio = step_ratio
            factor_squared = adaptation_factor * adaptation_factor
            if dual_component_norm > adaptation_threshold * primal_component_norm:
                proposed_ratio = step_ratio / factor_squared
            elif primal_component_norm > adaptation_threshold * dual_component_norm:
                proposed_ratio = step_ratio * factor_squared
            proposed_ratio = float(
                np.clip(
                    proposed_ratio,
                    minimum_step_ratio,
                    maximum_step_ratio,
                )
            )
            adapted = proposed_ratio != step_ratio
            if adapted:
                step_ratio = proposed_ratio
                primal_step = np.sqrt(step_product * step_ratio)
                dual_step = np.sqrt(step_product / step_ratio)
                extrapolated_gram = gram.copy()
        history["step_adapted"][-1] = adapted

        if best is None or stationarity_score < best[0]:
            best = (
                stationarity_score,
                gram.copy(),
                dual.copy(),
                logdet,
                internal_eigenvalues.copy(),
                iteration,
            )

        if iteration in save_steps_set:
            checkpoint, _ = _connectivity_checkpoint(
                next_eigenvalues,
                next_eigenvectors,
                householder,
                distance_scale,
                xp,
            )
            inverse_reconstruction_count += 1
            connectivity_at_steps[iteration] = checkpoint
            offdiagonal = checkpoint[np.triu_indices(n, k=1)]
            history["connectivity_offdiagonal_l2"][-1] = float(
                np.linalg.norm(offdiagonal)
            )

        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": "pdhg",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": objective,
                    "loss": relative_loss,
                    "entropy": entropy,
                    "relative_gradient_norm": stationarity_score,
                    "relative_primal_kkt_residual": residuals["primal_relative"],
                    "relative_dual_kkt_residual": residuals["dual_relative"],
                    "relative_eliminated_kkt_residual": (stationarity_score),
                    "primal_component_norm": primal_component_norm,
                    "dual_component_primal_norm": dual_component_norm,
                    "component_balance_ratio": component_balance_ratio,
                    "primal_step": used_primal_step,
                    "dual_step": used_dual_step,
                    "step_ratio": used_step_ratio,
                    "next_step_ratio": step_ratio,
                    "step_adapted": adapted,
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV-PDHG",
                    "general_method": "covariance_optimization",
                }
            )

        eliminated_converged = residuals["eliminated_norm"] <= (
            absolute_tolerance + relative_tolerance * residuals["eliminated_scale"]
        )
        if eliminated_converged:
            stopping_criterion_reached = True
            status = "eliminated_kkt_tolerance"
            message = "variance-whitened PDHG eliminated KKT tolerance reached"
            break

    if return_best and best is not None:
        (
            _,
            gram,
            dual,
            logdet,
            internal_eigenvalues,
            best_iteration,
        ) = best
    else:
        best_iteration = iteration

    final_eigenvalues, final_eigenvectors = _internal_eigh(gram, householder, xp)
    inverse_gram = _inverse_from_internal_spectrum(
        final_eigenvalues,
        final_eigenvectors,
        householder,
        xp,
    )
    inverse_reconstruction_count += 1
    full_gram = distance_scale * gram
    full_fitted = _num._squared_distances_from_gram(full_gram, array_module=xp)
    connectivity = _connectivity_from_scaled_inverse(inverse_gram, distance_scale, xp)

    if use_gpu:
        full_fitted = _num.cp.asnumpy(full_fitted)
        full_gram = _num.cp.asnumpy(full_gram)
        connectivity = _num.cp.asnumpy(connectivity)
        dual = _num.cp.asnumpy(dual)

    fitted_pairs_cpu = full_fitted[pair_i, pair_j]
    pair_stationarity = (
        fitted_pairs_cpu
        - target_pairs
        - 0.5 * pair_variance * connectivity[pair_i, pair_j]
    )
    pair_stationarity_norm = float(np.linalg.norm(pair_stationarity))
    pair_stationarity_scale = max(
        1.0,
        float(np.linalg.norm(fitted_pairs_cpu - target_pairs)),
        float(np.linalg.norm(0.5 * pair_variance * connectivity[pair_i, pair_j])),
    )

    for key, values in history.items():
        if key == "iteration":
            history[key] = np.asarray(values, dtype=np.int64)
        elif key == "step_adapted":
            history[key] = np.asarray(values, dtype=bool)
        else:
            history[key] = np.asarray(values, dtype=np.float64)
    history["is_returned_iterate"] = history["iteration"] == int(best_iteration)

    independent_residuals = _independent_eliminated_kkt_residuals(
        full_gram,
        target_pairs,
        pair_variance,
        pair_i,
        pair_j,
    )
    inverse_reconstruction_count += 1
    (
        independent_norm,
        independent_scale,
        independent_relative,
        _,
        independent_residual,
    ) = independent_residuals
    runtime_relative = float(
        history["relative_eliminated_kkt_residual"][best_iteration - 1]
    )
    independent_kkt_converged = independent_norm <= (
        absolute_tolerance + relative_tolerance * independent_scale
    )
    converged = bool(independent_kkt_converged)
    if converged:
        status = "optimality_tolerance"
        message = "independent eliminated KKT tolerance reached"
    elif stopping_criterion_reached:
        status = "independent_kkt_failed"
        message = (
            "runtime eliminated KKT tolerance was reached, but the returned "
            "Gram matrix failed the independent certificate"
        )

    final_offdiagonal = connectivity[np.triu_indices(n, k=1)]
    if len(history["connectivity_offdiagonal_l2"]) > 0:
        history["connectivity_offdiagonal_l2"][best_iteration - 1] = float(
            np.linalg.norm(final_offdiagonal)
        )

    info = {
        "converged": converged,
        "status": status,
        "message": message,
        "iterations": int(iteration),
        "returned_iteration": int(best_iteration),
        "objective": float(history["objective"][best_iteration - 1]),
        "loss": float(history["loss"][best_iteration - 1]),
        "entropy": float(history["entropy"][best_iteration - 1]),
        "relative_gradient_norm": float(independent_relative),
        "relative_primal_kkt_residual": float(
            history["relative_primal_kkt_residual"][best_iteration - 1]
        ),
        "relative_dual_kkt_residual": float(
            history["relative_dual_kkt_residual"][best_iteration - 1]
        ),
        "relative_eliminated_kkt_residual": float(independent_relative),
        "runtime_relative_eliminated_kkt_residual": runtime_relative,
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
        "termination_internal_kkt_converged": bool(stopping_criterion_reached),
        "independent_kkt_converged": bool(independent_kkt_converged),
        "independent_kkt_recomputed_from_returned_gram": True,
        "relative_tolerance": float(relative_tolerance),
        "absolute_tolerance": float(absolute_tolerance),
        "observed_pair_count": len(target_pairs),
        "noise_model": noise_model,
        "noise_parameter": noise_parameter,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization_info,
        "dual_initialization": chosen_dual_initialization,
        "algorithm": "pdhg",
        "coordinate_parameterization": "centered_householder",
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "cupy_version": _num.cp.__version__ if use_gpu else None,
        "distance_scale": float(distance_scale),
        "weighted_operator_norm": operator_norm_info,
        "inverse_reconstruction_count": int(inverse_reconstruction_count),
        "inverse_reconstructed_each_iteration": False,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "pdhg": {
            "variance_whitened": True,
            "inverse_free_runtime_kkt": True,
            "weighted_residual_balancing": True,
            "theta": float(theta),
            "step_safety": float(step_safety),
            "step_product": float(step_product),
            "initial_step_ratio": float(history["step_ratio"][0]),
            "final_step_ratio": float(step_ratio),
            "initial_dual_step": float(initial_dual_step),
            "adaptive_steps": bool(adaptive_steps),
            "adaptation_interval": int(adaptation_interval),
            "adaptation_threshold": float(adaptation_threshold),
            "adaptation_factor": float(adaptation_factor),
            "return_best": bool(return_best),
            "stopping_criterion": "eliminated_kkt",
        },
        "wall_seconds": time.perf_counter() - start_time,
        "objective_definition": (
            "-1.5*logdet(B_internal) + 0.5*||V^(-1/2)*(D(B)-D_obs)||^2"
        ),
        "stationarity_definition": (
            "D_adjoint((D(B)-D_obs)/noise_variance)-1.5*B_pseudoinverse"
        ),
        "pair_stationarity_definition": ("D_fit_ij-D_obs_ij-noise_variance_ij*A_ij/2"),
        "logged_metric_timing": (
            "post-update whitened PDHG iterate; runtime KKT obtained from "
            "proximal optimality without B inverse"
        ),
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_pdhg stopped without "
            f"satisfying the KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return full_fitted, full_gram, connectivity, info


def _combine_hybrid_histories(pdhg_history, fista_history):
    """Combine phase histories on one global accepted-update axis."""
    pdhg_count = len(pdhg_history["iteration"])
    fista_count = len(fista_history["iteration"])
    combined = {
        "iteration": np.arange(1, pdhg_count + fista_count + 1, dtype=np.int64),
        "phase": np.asarray(["pdhg"] * pdhg_count + ["fista"] * fista_count),
        "phase_iteration": np.concatenate(
            (
                np.asarray(pdhg_history["iteration"], dtype=np.int64),
                np.asarray(fista_history["iteration"], dtype=np.int64),
            )
        ),
    }
    history_keys = (set(pdhg_history) | set(fista_history)) - {"iteration"}
    for key in sorted(history_keys):
        pdhg_values = pdhg_history.get(key)
        fista_values = fista_history.get(key)
        if pdhg_values is None:
            pdhg_values = np.full(pdhg_count, np.nan)
        if fista_values is None:
            fista_values = np.full(fista_count, np.nan)
        combined[key] = np.concatenate(
            (np.asarray(pdhg_values), np.asarray(fista_values))
        )
    return combined


def fit_gaussian_noise_covariance_hybrid(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    handoff_relative_tolerance=1e-2,
    save_steps=None,
    progress_callback=None,
):
    """Fit COV globally with PDHG and finish with direct-Gram FISTA.

    PDHG starts from a Rouse chain or an explicitly supplied connectivity
    matrix and runs until its independently certified KKT residual reaches
    ``handoff_relative_tolerance``. The returned physical Gram matrix is handed
    directly to monotone FISTA without a connectivity round trip or scalar
    recalibration. ``max_iterations`` is a single total accepted-update budget
    shared by both phases. The final returned Gram matrix is independently
    certified against ``relative_tolerance``.
    """
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if isinstance(handoff_relative_tolerance, (bool, np.bool_)) or not np.isscalar(
        handoff_relative_tolerance
    ):
        raise ValueError("handoff_relative_tolerance must be a positive finite scalar")
    handoff_relative_tolerance = float(handoff_relative_tolerance)
    if not np.isfinite(handoff_relative_tolerance) or handoff_relative_tolerance <= 0.0:
        raise ValueError("handoff_relative_tolerance must be a positive finite scalar")
    if relative_tolerance > 0.0 and handoff_relative_tolerance < relative_tolerance:
        raise ValueError(
            "handoff_relative_tolerance must not be smaller than the final "
            "relative_tolerance"
        )

    save_steps_set = set(save_steps or ())
    if any(
        not isinstance(step, (int, np.integer)) or step < 1 or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    def phase_progress_callback(phase, offset):
        if progress_callback is None:
            return None

        def report(event):
            event = dict(event)
            event["optimizer"] = "hybrid"
            event["phase"] = phase
            if event.get("stage") == "covariance_optimization":
                event["phase_iteration"] = int(event["iteration"])
                event["phase_total"] = int(event["total"])
                event["iteration"] = offset + int(event["iteration"])
                event["total"] = int(max_iterations)
                event["method"] = "COV-hybrid"
            progress_callback(event)

        return report

    start_time = time.perf_counter()
    pdhg_start_time = time.perf_counter()
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
        progress_callback=phase_progress_callback("pdhg", 0),
    )
    pdhg_wall_seconds = time.perf_counter() - pdhg_start_time
    fitted, gram, connectivity, pdhg_info = pdhg_result
    pdhg_iterations = int(pdhg_info["iterations"])
    handoff_residual = float(pdhg_info["relative_eliminated_kkt_residual"])
    residual_norm = float(pdhg_info["stationarity_residual_norm"])
    residual_scale = float(pdhg_info["stationarity_residual_scale"])
    handoff_reached = residual_norm <= (
        absolute_tolerance + handoff_relative_tolerance * residual_scale
    )
    final_kkt_converged = residual_norm <= (
        absolute_tolerance + relative_tolerance * residual_scale
    )
    handoff = {
        "reached": bool(handoff_reached),
        "relative_tolerance": handoff_relative_tolerance,
        "absolute_tolerance": float(absolute_tolerance),
        "relative_eliminated_kkt_residual": handoff_residual,
        "objective": float(pdhg_info["objective"]),
        "iterations": pdhg_iterations,
        "returned_iteration": int(pdhg_info["returned_iteration"]),
    }
    remaining_iterations = max_iterations - pdhg_iterations
    if final_kkt_converged or not handoff_reached or remaining_iterations <= 0:
        if final_kkt_converged:
            status = "optimality_tolerance"
            message = "PDHG returned Gram passed the final independent KKT tolerance"
        elif not handoff_reached:
            status = "pdhg_handoff_not_reached"
            message = "PDHG did not reach the hybrid handoff tolerance"
        else:
            status = "iteration_budget_exhausted_at_handoff"
            message = "PDHG reached the handoff at the end of the iteration budget"
        combined_history = _combine_hybrid_histories(
            pdhg_info["history"],
            {"iteration": np.asarray([], dtype=np.int64)},
        )
        combined_history["is_returned_iterate"] = combined_history["iteration"] == int(
            pdhg_info["returned_iteration"]
        )
        info = dict(pdhg_info)
        info.update(
            {
                "algorithm": "hybrid",
                "coordinate_parameterization": "centered_hybrid",
                "converged": bool(final_kkt_converged),
                "status": status,
                "message": message,
                "relative_tolerance": float(relative_tolerance),
                "absolute_tolerance": float(absolute_tolerance),
                "termination_internal_kkt_converged": bool(final_kkt_converged),
                "independent_kkt_converged": bool(final_kkt_converged),
                "handoff_relative_tolerance": (handoff_relative_tolerance),
                "handoff": handoff,
                "phase_iterations": {
                    "pdhg": pdhg_iterations,
                    "fista": 0,
                },
                "phase_wall_seconds": {
                    "pdhg": pdhg_wall_seconds,
                    "fista": 0.0,
                },
                "history": combined_history,
                "wall_seconds": time.perf_counter() - start_time,
            }
        )
        if not final_kkt_converged:
            warnings.warn(
                "fit_gaussian_noise_covariance_hybrid stopped before FISTA "
                f"refinement (status={status})",
                RuntimeWarning,
                stacklevel=2,
            )
        return fitted, gram, connectivity, info

    fista_save_steps = sorted(
        step - pdhg_iterations for step in save_steps_set if step > pdhg_iterations
    )
    fista_start_time = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(
                "fit_gaussian_noise_covariance_fista stopped without " "satisfying.*"
            ),
            category=RuntimeWarning,
        )
        fista_result = fit_gaussian_noise_covariance_fista(
            squared_distances,
            noise_variance,
            relative_noise_std=relative_noise_std,
            initial_gram=gram,
            use_gpu=use_gpu,
            max_iterations=remaining_iterations,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            distance_scale=pdhg_info["distance_scale"],
            _precomputed_operator_norm=(
                pdhg_info["weighted_operator_norm"]["safe_operator_norm_squared"],
                pdhg_info["weighted_operator_norm"],
            ),
            save_steps=fista_save_steps,
            progress_callback=phase_progress_callback("fista", pdhg_iterations),
        )
    fista_wall_seconds = time.perf_counter() - fista_start_time
    fitted, gram, connectivity, fista_info = fista_result
    fista_iterations = int(fista_info["iterations"])
    independent_norm = float(fista_info["stationarity_residual_norm"])
    independent_scale = float(fista_info["stationarity_residual_scale"])
    independent_relative = float(fista_info["relative_eliminated_kkt_residual"])
    independent_kkt_converged = bool(fista_info["independent_kkt_converged"])
    internal_kkt_converged = bool(fista_info["termination_internal_kkt_converged"])
    converged = internal_kkt_converged and independent_kkt_converged
    if converged:
        status = "optimality_tolerance"
        message = (
            "PDHG handoff, FISTA refinement, and independent final KKT "
            "certificate reached"
        )
    elif internal_kkt_converged:
        status = "independent_kkt_failed"
        message = (
            "FISTA reached its internal tolerance, but the returned Gram "
            "matrix failed the independent eliminated KKT certificate"
        )
    else:
        status = f"fista_{fista_info['status']}"
        message = f"FISTA refinement stopped: {fista_info['message']}"

    combined_history = _combine_hybrid_histories(
        pdhg_info["history"], fista_info["history"]
    )
    returned_iteration = pdhg_iterations + fista_iterations
    combined_history["is_returned_iterate"] = (
        combined_history["iteration"] == returned_iteration
    )
    connectivity_at_steps = dict(pdhg_info["connectivity_matrix_at_steps"])
    connectivity_at_steps.update(
        {
            pdhg_iterations + int(step): matrix
            for step, matrix in fista_info["connectivity_matrix_at_steps"].items()
        }
    )
    info = dict(fista_info)
    info.update(
        {
            "converged": bool(converged),
            "status": status,
            "message": message,
            "iterations": pdhg_iterations + fista_iterations,
            "returned_iteration": returned_iteration,
            "loss": float(fista_info["loss"]),
            "entropy": float(fista_info["entropy"]),
            "relative_gradient_norm": float(independent_relative),
            "relative_eliminated_kkt_residual": float(independent_relative),
            "runtime_relative_eliminated_kkt_residual": float(
                fista_info["runtime_relative_eliminated_kkt_residual"]
            ),
            "stationarity_residual_norm": float(independent_norm),
            "stationarity_residual_scale": float(independent_scale),
            "relative_stationarity_residual": float(independent_relative),
            "maximum_absolute_stationarity_residual": float(
                fista_info["maximum_absolute_stationarity_residual"]
            ),
            "pair_stationarity_residual_norm": float(
                fista_info["pair_stationarity_residual_norm"]
            ),
            "relative_pair_stationarity_residual": float(
                fista_info["relative_pair_stationarity_residual"]
            ),
            "maximum_absolute_pair_stationarity_residual": float(
                fista_info["maximum_absolute_pair_stationarity_residual"]
            ),
            "termination_internal_kkt_converged": internal_kkt_converged,
            "independent_kkt_converged": bool(independent_kkt_converged),
            "independent_kkt_recomputed_from_returned_gram": True,
            "relative_tolerance": float(relative_tolerance),
            "absolute_tolerance": float(absolute_tolerance),
            "handoff_relative_tolerance": handoff_relative_tolerance,
            "handoff": {
                **handoff,
                "physical_gram_used_directly": True,
                "scalar_recalibration": None,
            },
            "phase_iterations": {
                "pdhg": pdhg_iterations,
                "fista": fista_iterations,
            },
            "phase_wall_seconds": {
                "pdhg": pdhg_wall_seconds,
                "fista": fista_wall_seconds,
            },
            "initialization": pdhg_info["initialization"],
            "algorithm": "hybrid",
            "coordinate_parameterization": "centered_hybrid",
            "history": combined_history,
            "connectivity_matrix_at_steps": connectivity_at_steps,
            "wall_seconds": time.perf_counter() - start_time,
            "pdhg": pdhg_info["pdhg"],
            "weighted_operator_norm": pdhg_info["weighted_operator_norm"],
            "inverse_reconstruction_count": (
                pdhg_info["inverse_reconstruction_count"]
                + fista_info["inverse_reconstruction_count"]
            ),
            "inverse_reconstructed_each_iteration": False,
            "stationarity_definition": pdhg_info["stationarity_definition"],
            "pair_stationarity_definition": pdhg_info["pair_stationarity_definition"],
            "logged_metric_timing": (
                "post-update PDHG then FISTA iterates on one global axis; "
                "final certificate recomputed from returned Gram matrix"
            ),
        }
    )
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_hybrid stopped without satisfying "
            f"the final KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted, gram, connectivity, info
