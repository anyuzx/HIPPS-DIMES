"""Variance-whitened PDHG for the Gaussian noise-aware COV objective.

This module implements three large-system improvements over the scalar-step
PDHG solver in :mod:`hipps_dimes.covariance_pdhg`:

1. variance whitening of the distance residuals,
2. inverse-free eliminated-KKT diagnostics during the iterations, and
3. residual balancing in the common primal matrix space.

The statistical objective is unchanged::

    -3/2 logdet'(B) + 1/2 sum_a (D_a(B) - d_a)^2 / v_a.

Writing A = V^{-1/2} D and b = V^{-1/2} d gives the equivalent form::

    -3/2 logdet'(B) + 1/2 ||A(B) - b||^2.

The whitened dual has unit quadratic curvature, so every dual coordinate has
the same fixed-primal contraction factor. The primal log-determinant proximal
map preserves positive definiteness on the center-of-mass-free subspace.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from . import covariance_pdhg as _base
from . import numerics as _num


_ENTROPY_COEFFICIENT = 1.5


def _scalar(value) -> float:
    return float(value.item()) if hasattr(value, "item") else float(value)


def _norm(value, xp) -> float:
    return _scalar(xp.sqrt(xp.sum(value * value)))


def _weighted_distance_operator(
    gram, inverse_standard_deviation, pair_i, pair_j, xp
):
    return (
        _base._distance_operator(gram, pair_i, pair_j, xp)
        * inverse_standard_deviation
    )


def _weighted_distance_adjoint(
    values, inverse_standard_deviation, n, pair_i, pair_j, xp
):
    return _base._distance_adjoint(
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
    matrix = _base._center(matrix, xp)
    eigenvalues, eigenvectors = _base._internal_eigh(
        matrix, householder_vector, xp
    )
    root = xp.hypot(
        eigenvalues,
        np.sqrt(4.0 * _ENTROPY_COEFFICIENT * step_size),
    )
    nonnegative = xp.maximum(eigenvalues, 0.0)
    nonpositive = xp.minimum(eigenvalues, 0.0)
    updated = xp.where(
        eigenvalues >= 0.0,
        0.5 * (nonnegative + root),
        (2.0 * _ENTROPY_COEFFICIENT * step_size)
        / (root - nonpositive),
    )
    internal = (eigenvectors * updated) @ eigenvectors.T
    gram = _base._full_centered_from_internal(
        internal, householder_vector, xp
    )
    logdet = _scalar(xp.sum(xp.log(updated)))
    return gram, logdet, updated, eigenvectors


def _inverse_from_internal_spectrum(
    eigenvalues, eigenvectors, householder_vector, xp
):
    inverse_internal = (
        eigenvectors * (1.0 / eigenvalues)
    ) @ eigenvectors.T
    return _base._full_centered_from_internal(
        inverse_internal, householder_vector, xp
    )


def _estimate_weighted_operator_norm_squared(
    initial_gram,
    inverse_standard_deviation,
    n,
    pair_i,
    pair_j,
    xp,
    *,
    max_iterations,
    relative_tolerance,
    safety_factor,
    progress_callback=None,
):
    """Estimate ``||V^-1/2 D||^2`` by matrix-free power iteration."""
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("operator_norm_power_iterations must be positive")
    if not np.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise ValueError("operator_norm_tolerance must be positive and finite")
    if not np.isfinite(safety_factor) or safety_factor <= 1.0:
        raise ValueError("operator_norm_safety must exceed one")

    vector = _base._center(initial_gram.copy(), xp)
    vector_norm = _norm(vector, xp)
    if not np.isfinite(vector_norm) or vector_norm <= 0.0:
        vector = xp.eye(n, dtype=xp.float64)
        vector = _base._center(vector, xp)
        vector_norm = _norm(vector, xp)
    vector = vector / vector_norm

    rayleigh = 0.0
    residual_norm = np.inf
    relative_residual = np.inf
    iteration = 0
    start_time = time.perf_counter()
    for iteration in range(1, int(max_iterations) + 1):
        transformed = _weighted_distance_operator(
            vector,
            inverse_standard_deviation,
            pair_i,
            pair_j,
            xp,
        )
        normal_vector = _weighted_distance_adjoint(
            transformed,
            inverse_standard_deviation,
            n,
            pair_i,
            pair_j,
            xp,
        )
        normal_vector = _base._center(normal_vector, xp)
        rayleigh = _scalar(xp.sum(transformed * transformed))
        residual = normal_vector - rayleigh * vector
        residual_norm = _norm(residual, xp)
        relative_residual = residual_norm / max(
            abs(rayleigh), np.finfo(np.float64).tiny
        )
        normal_norm = _norm(normal_vector, xp)
        if not np.isfinite(normal_norm) or normal_norm <= 0.0:
            raise RuntimeError(
                "weighted distance operator has a nonpositive power-iteration norm"
            )
        vector = normal_vector / normal_norm
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "phase": "weighted_operator_norm",
                    "iteration": iteration,
                    "total": int(max_iterations),
                    "operator_norm_squared_rayleigh": float(rayleigh),
                    "operator_norm_relative_residual": float(
                        relative_residual
                    ),
                    "method": "COV-PDHG-whitened",
                }
            )
        if iteration >= 3 and relative_residual <= relative_tolerance:
            break

    transformed = _weighted_distance_operator(
        vector,
        inverse_standard_deviation,
        pair_i,
        pair_j,
        xp,
    )
    normal_vector = _weighted_distance_adjoint(
        transformed,
        inverse_standard_deviation,
        n,
        pair_i,
        pair_j,
        xp,
    )
    normal_vector = _base._center(normal_vector, xp)
    rayleigh = _scalar(xp.sum(transformed * transformed))
    residual_norm = _norm(normal_vector - rayleigh * vector, xp)
    relative_residual = residual_norm / max(
        abs(rayleigh), np.finfo(np.float64).tiny
    )
    estimated_upper = max(
        rayleigh + residual_norm,
        rayleigh,
        np.finfo(np.float64).tiny,
    )
    safe_estimate = safety_factor * estimated_upper
    conservative_bound = 2.0 * n * _scalar(
        xp.max(inverse_standard_deviation**2)
    )
    return safe_estimate, {
        "iterations": int(iteration),
        "rayleigh_quotient": float(rayleigh),
        "residual_norm": float(residual_norm),
        "relative_residual": float(relative_residual),
        "estimated_upper_before_safety": float(estimated_upper),
        "safety_factor": float(safety_factor),
        "safe_operator_norm_squared": float(safe_estimate),
        "conservative_operator_norm_squared_bound": float(
            conservative_bound
        ),
        "wall_seconds": time.perf_counter() - start_time,
    }


def _choose_initial_whitened_dual(
    gram,
    inverse_gram,
    whitened_target,
    inverse_standard_deviation,
    standard_deviation,
    n,
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
    candidates = {
        "zero": xp.zeros_like(whitened_target),
        "residual": fitted - whitened_target,
        "connectivity": (
            standard_deviation
            * (
                -_ENTROPY_COEFFICIENT
                * inverse_gram[pair_i, pair_j]
            )
        ),
    }
    if mode != "auto":
        if mode not in candidates:
            raise ValueError(
                "dual_initialization must be 'auto', 'zero', 'residual', or "
                "'connectivity'"
            )
        chosen = mode
    else:
        scores = {}
        entropy_gradient = _ENTROPY_COEFFICIENT * inverse_gram
        for name, candidate in candidates.items():
            dual_adjoint = _weighted_distance_adjoint(
                candidate,
                inverse_standard_deviation,
                n,
                pair_i,
                pair_j,
                xp,
            )
            primal_residual = dual_adjoint - entropy_gradient
            dual_residual = fitted - whitened_target - candidate
            dual_residual_adjoint = _weighted_distance_adjoint(
                dual_residual,
                inverse_standard_deviation,
                n,
                pair_i,
                pair_j,
                xp,
            )
            scale = max(
                _norm(dual_adjoint, xp),
                _norm(entropy_gradient, xp),
                np.finfo(np.float64).tiny,
            )
            eliminated = primal_residual + dual_residual_adjoint
            scores[name] = _norm(eliminated, xp) / scale
        chosen = min(scores, key=scores.get)
    return candidates[chosen].copy(), chosen


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
    connectivity = _base._connectivity_from_scaled_inverse(
        inverse_gram, distance_scale, xp
    )
    if xp is np:
        return np.asarray(connectivity).copy(), inverse_gram
    return _num.cp.asnumpy(connectivity), inverse_gram


def fit_gaussian_noise_covariance_pdhg_whitened(
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
    dual_initialization="auto",
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
            "whitened PDHG GPU optimization was requested, but CuPy and an "
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
        not isinstance(step, (int, np.integer))
        or step < 1
        or step > max_iterations
        for step in save_steps_set
    ):
        raise ValueError("save_steps must lie between 1 and max_iterations")

    basis = _num._centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    reduced_gram, initialization_info = (
        _num._initialize_gaussian_reduced_gram(
            observed,
            pair_i,
            pair_j,
            basis,
            initial_connectivity,
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

    full_gram = _base._center(
        solver_basis @ reduced_gram @ solver_basis.T, xp
    )
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
    householder = _base._householder_vector(n, xp)

    inverse_reconstruction_count = 0
    inverse_gram, logdet, internal_eigenvalues = (
        _base._internal_inverse_and_logdet(gram, householder, xp)
    )
    inverse_reconstruction_count += 1
    dual, chosen_dual_initialization = _choose_initial_whitened_dual(
        gram,
        inverse_gram,
        whitened_target,
        inverse_standard_deviation,
        standard_deviation,
        n,
        solver_pair_i,
        solver_pair_j,
        xp,
        dual_initialization,
    )

    operator_norm_squared, operator_norm_info = (
        _estimate_weighted_operator_norm_squared(
            gram,
            inverse_standard_deviation,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
            max_iterations=operator_norm_power_iterations,
            relative_tolerance=operator_norm_tolerance,
            safety_factor=operator_norm_safety,
            progress_callback=progress_callback,
        )
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
        next_dual = (
            dual
            + dual_step * (extrapolated_distances - whitened_target)
        ) / (1.0 + dual_step)
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
        dual_change = _norm(next_dual - dual, xp) / max(
            _norm(dual, xp), 1.0
        )
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
        data_objective = 0.5 * _scalar(
            xp.sum(standardized_data_residual**2)
        )
        negative_entropy = (
            -_ENTROPY_COEFFICIENT
            * (next_logdet + n_modes * np.log(distance_scale))
        )
        objective = negative_entropy + data_objective
        fitted_normalized = whitened_fitted * standard_deviation
        data_residual = fitted_normalized - target
        relative_loss = _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual / target)))
        )
        weighted_rmse = _scalar(
            xp.sqrt(xp.mean(standardized_data_residual**2))
        )
        distance_rmse = distance_scale * _scalar(
            xp.sqrt(xp.mean(data_residual**2))
        )
        entropy = (
            next_logdet
            + n_modes * np.log(distance_scale)
            - n_modes * np.log(3.0)
        )
        minimum_eigenvalue = _scalar(next_eigenvalues[0])
        maximum_eigenvalue = _scalar(next_eigenvalues[-1])
        condition_number = maximum_eigenvalue / minimum_eigenvalue
        precision_frobenius = (
            3.0
            / distance_scale
            * _scalar(xp.sqrt(xp.sum(1.0 / next_eigenvalues**2)))
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
        history["relative_primal_kkt_residual"].append(
            residuals["primal_relative"]
        )
        history["relative_dual_kkt_residual"].append(
            residuals["dual_relative"]
        )
        history["relative_eliminated_kkt_residual"].append(
            stationarity_score
        )
        history["relative_gradient_norm"].append(stationarity_score)
        history["relative_stationarity_residual"].append(
            stationarity_score
        )
        history["primal_component_norm"].append(primal_component_norm)
        history["dual_component_primal_norm"].append(
            dual_component_norm
        )
        history["component_balance_ratio"].append(
            component_balance_ratio
        )
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
        history["internal_precision_frobenius_norm"].append(
            precision_frobenius
        )
        history["connectivity_offdiagonal_l2"].append(np.nan)

        if adaptive_steps and iteration % adaptation_interval == 0:
            proposed_ratio = step_ratio
            factor_squared = adaptation_factor * adaptation_factor
            if (
                dual_component_norm
                > adaptation_threshold * primal_component_norm
            ):
                proposed_ratio = step_ratio / factor_squared
            elif (
                primal_component_norm
                > adaptation_threshold * dual_component_norm
            ):
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
                    "phase": "pdhg_whitened",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": objective,
                    "loss": relative_loss,
                    "entropy": entropy,
                    "relative_gradient_norm": stationarity_score,
                    "relative_primal_kkt_residual": residuals[
                        "primal_relative"
                    ],
                    "relative_dual_kkt_residual": residuals[
                        "dual_relative"
                    ],
                    "relative_eliminated_kkt_residual": (
                        stationarity_score
                    ),
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
                    "method": "COV-PDHG-whitened",
                    "general_method": "covariance_optimization",
                }
            )

        eliminated_converged = residuals["eliminated_norm"] <= (
            absolute_tolerance
            + relative_tolerance * residuals["eliminated_scale"]
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

    final_eigenvalues, final_eigenvectors = _base._internal_eigh(
        gram, householder, xp
    )
    inverse_gram = _inverse_from_internal_spectrum(
        final_eigenvalues,
        final_eigenvectors,
        householder,
        xp,
    )
    inverse_reconstruction_count += 1
    full_gram = distance_scale * gram
    full_fitted = _num._squared_distances_from_gram(
        full_gram, array_module=xp
    )
    connectivity = _base._connectivity_from_scaled_inverse(
        inverse_gram, distance_scale, xp
    )

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
        float(
            np.linalg.norm(
                0.5 * pair_variance * connectivity[pair_i, pair_j]
            )
        ),
    )

    for key, values in history.items():
        if key == "iteration":
            history[key] = np.asarray(values, dtype=np.int64)
        elif key == "step_adapted":
            history[key] = np.asarray(values, dtype=bool)
        else:
            history[key] = np.asarray(values, dtype=np.float64)

    independent_residuals = _base._independent_eliminated_kkt_residuals(
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
        "algorithm": "pdhg_whitened",
        "coordinate_parameterization": "centered_householder",
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": _num.get_gpu_name() if use_gpu else None,
        "cupy_version": _num.cp.__version__ if use_gpu else None,
        "distance_scale": float(distance_scale),
        "weighted_operator_norm": operator_norm_info,
        "inverse_reconstruction_count": int(
            inverse_reconstruction_count
        ),
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
        "pair_stationarity_definition": (
            "D_fit_ij-D_obs_ij-noise_variance_ij*A_ij/2"
        ),
        "logged_metric_timing": (
            "post-update whitened PDHG iterate; runtime KKT obtained from "
            "proximal optimality without B inverse"
        ),
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_pdhg_whitened stopped without "
            f"satisfying the KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return full_fitted, full_gram, connectivity, info


fit_gaussian_noise_covariance_preconditioned_pdhg = (
    fit_gaussian_noise_covariance_pdhg_whitened
)
