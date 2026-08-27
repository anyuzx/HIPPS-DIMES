"""PDHG and PDHG-to-Newton solvers for the noise-aware COV objective.

The solver minimizes

    -3/2 logdet'(B) + 1/2 sum_a (D_a(B) - d_a)^2 / v_a

on centered Gram matrices that are positive definite on the internal subspace.
It uses the Chambolle--Pock/PDHG splitting with analytic proximal operators for
both the Gaussian data likelihood and the log-determinant entropy barrier.
Unlike Newton--CG, each iteration has no inner linear solve and cannot leave the
valid covariance cone.
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
    transformed = _householder_congruence(
        centered_matrix, householder_vector, xp
    )
    return _symmetrize(transformed[:-1, :-1])


def _full_centered_from_internal(internal_matrix, householder_vector, xp):
    n = internal_matrix.shape[0] + 1
    transformed = xp.zeros((n, n), dtype=xp.float64)
    transformed[:-1, :-1] = internal_matrix
    return _center(
        _householder_congruence(transformed, householder_vector, xp), xp
    )


def _internal_eigh(centered_matrix, householder_vector, xp):
    internal = _internal_block(centered_matrix, householder_vector, xp)
    if xp is np:
        eigenvalues, eigenvectors = np.linalg.eigh(internal)
    else:
        with _num.cupyx.errstate(linalg="raise"):
            eigenvalues, eigenvectors = xp.linalg.eigh(internal)
    return eigenvalues, eigenvectors


def _internal_inverse_and_logdet(
    centered_matrix, householder_vector, xp
):
    eigenvalues, eigenvectors = _internal_eigh(
        centered_matrix, householder_vector, xp
    )
    minimum = _scalar(eigenvalues[0])
    if not np.isfinite(minimum) or minimum <= 0.0:
        raise RuntimeError(
            "PDHG requires a Gram matrix that is positive definite on the "
            "internal subspace"
        )
    inverse_internal = (
        eigenvectors * (1.0 / eigenvalues)
    ) @ eigenvectors.T
    inverse_full = _full_centered_from_internal(
        inverse_internal, householder_vector, xp
    )
    logdet = _scalar(xp.sum(xp.log(eigenvalues)))
    return inverse_full, logdet, eigenvalues


def _prox_centered_negative_logdet(
    matrix, step_size: float, householder_vector, xp
):
    """Apply prox of ``-3/2 logdet'`` while preserving the COM zero mode."""
    matrix = _center(matrix, xp)
    eigenvalues, eigenvectors = _internal_eigh(
        matrix, householder_vector, xp
    )
    updated = 0.5 * (
        eigenvalues
        + xp.sqrt(
            xp.square(eigenvalues)
            + 4.0 * _ENTROPY_COEFFICIENT * step_size
        )
    )
    internal = (eigenvectors * updated) @ eigenvectors.T
    inverse_internal = (
        eigenvectors * (1.0 / updated)
    ) @ eigenvectors.T
    gram = _full_centered_from_internal(
        internal, householder_vector, xp
    )
    inverse = _full_centered_from_internal(
        inverse_internal, householder_vector, xp
    )
    logdet = _scalar(xp.sum(xp.log(updated)))
    return gram, inverse, logdet, updated


def _distance_operator(gram, pair_i, pair_j, xp):
    diagonal = xp.diag(gram)
    return (
        diagonal[pair_i]
        + diagonal[pair_j]
        - 2.0 * gram[pair_i, pair_j]
    )


def _distance_adjoint(values, n, pair_i, pair_j, xp):
    """Adjoint of the unique-pair squared-distance operator."""
    pair_matrix = xp.zeros((n, n), dtype=xp.float64)
    pair_matrix[pair_i, pair_j] = values
    pair_matrix[pair_j, pair_i] = values
    row_sum = xp.sum(pair_matrix, axis=1)
    result = -pair_matrix
    result[xp.diag_indices(n)] = row_sum
    return _symmetrize(result)


def _relative_kkt_residuals(
    gram,
    inverse_gram,
    dual,
    target,
    variance,
    n,
    pair_i,
    pair_j,
    xp,
):
    dual_adjoint = _distance_adjoint(
        dual, n, pair_i, pair_j, xp
    )
    entropy_precision = _ENTROPY_COEFFICIENT * inverse_gram
    primal_residual = dual_adjoint - entropy_precision
    primal_norm = _norm(primal_residual, xp)
    primal_scale = max(
        _norm(dual_adjoint, xp),
        _norm(entropy_precision, xp),
        np.finfo(np.float64).tiny,
    )

    fitted = _distance_operator(gram, pair_i, pair_j, xp)
    distance_residual = fitted - target
    variance_dual = variance * dual
    dual_residual = distance_residual - variance_dual
    dual_norm = _norm(dual_residual, xp)
    dual_scale = max(
        _norm(distance_residual, xp),
        _norm(variance_dual, xp),
        np.finfo(np.float64).tiny,
    )
    return (
        primal_norm,
        primal_scale,
        primal_norm / primal_scale,
        dual_norm,
        dual_scale,
        dual_norm / dual_scale,
        fitted,
    )


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
    data_gradient = _distance_adjoint(
        eliminated_dual, n, pair_i, pair_j, xp
    )
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


def _choose_initial_dual(
    gram,
    inverse_gram,
    target,
    variance,
    n,
    pair_i,
    pair_j,
    xp,
    mode,
):
    fitted = _distance_operator(gram, pair_i, pair_j, xp)
    candidates = {
        "zero": xp.zeros_like(target),
        "residual": (fitted - target) / variance,
        "connectivity": (
            -_ENTROPY_COEFFICIENT
            * inverse_gram[pair_i, pair_j]
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
        for name, candidate in candidates.items():
            residuals = _relative_kkt_residuals(
                gram,
                inverse_gram,
                candidate,
                target,
                variance,
                n,
                pair_i,
                pair_j,
                xp,
            )
            scores[name] = max(residuals[2], residuals[5])
        chosen = min(scores, key=scores.get)
    return candidates[chosen].copy(), chosen


def _initial_step_ratio(
    variance,
    primal_relative_residual,
    dual_relative_residual,
    operator_product,
    minimum,
    maximum,
    xp,
):
    """Choose a dimensionless primal/dual step ratio from local scales."""
    typical_variance = _scalar(xp.median(variance))
    ratio = operator_product * typical_variance
    residual_balance = np.sqrt(
        (float(primal_relative_residual) + 1e-16)
        / (float(dual_relative_residual) + 1e-16)
    )
    ratio *= float(np.clip(residual_balance, 0.1, 10.0))
    return float(np.clip(ratio, minimum, maximum))


def _connectivity_from_scaled_inverse(inverse_gram, distance_scale, xp):
    connectivity = -(3.0 / distance_scale) * inverse_gram
    connectivity = _symmetrize(connectivity)
    return _num.a2a(connectivity)


def fit_gaussian_noise_covariance_pdhg(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initialization="rouse",
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
    theta=1.0,
    step_safety=0.99,
    step_ratio=None,
    adaptive_steps=True,
    adaptation_interval=20,
    adaptation_threshold=5.0,
    adaptation_factor=2.0,
    minimum_step_ratio=1e-8,
    maximum_step_ratio=1e8,
    dual_initialization="auto",
    distance_scale=None,
    return_best=True,
):
    """Fit the Gaussian noise-aware COV objective with PDHG.

    Parameters shared with ``fit_gaussian_noise_covariance`` have the same
    meaning. ``step_ratio`` is ``tau / sigma``. When omitted, a scale-aware
    value is selected automatically and, by default, balanced using the KKT
    residuals while the product ``tau * sigma`` remains fixed. The returned
    Gram matrix is accepted only after a fresh, dual-eliminated KKT check. The
    bound
    ``||D||^2 <= 2N`` for the complete centered distance operator guarantees
    ``tau * sigma * ||D||^2 < 1`` for every observed-pair subset.
    """
    if use_gpu and not _num.is_gpu_available():
        raise RuntimeError(
            "PDHG GPU optimization was requested, but CuPy and an accessible "
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
    if not np.isfinite(theta) or not (0.0 <= theta <= 1.0):
        raise ValueError("theta must lie in [0, 1]")
    if not np.isfinite(step_safety) or not (0.0 < step_safety < 1.0):
        raise ValueError("step_safety must lie strictly between zero and one")
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

    full_gram = _center(
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
    householder = _householder_vector(n, xp)
    inverse_gram, logdet, internal_eigenvalues = (
        _internal_inverse_and_logdet(gram, householder, xp)
    )
    dual, chosen_dual_initialization = _choose_initial_dual(
        gram,
        inverse_gram,
        target,
        variance,
        n,
        solver_pair_i,
        solver_pair_j,
        xp,
        dual_initialization,
    )
    initial_residuals = _relative_kkt_residuals(
        gram,
        inverse_gram,
        dual,
        target,
        variance,
        n,
        solver_pair_i,
        solver_pair_j,
        xp,
    )

    operator_norm_squared_bound = 2.0 * n
    base_step = step_safety / np.sqrt(operator_norm_squared_bound)
    operator_product = base_step * base_step
    if step_ratio is None:
        step_ratio = _initial_step_ratio(
            variance,
            initial_residuals[2],
            initial_residuals[5],
            operator_product,
            minimum_step_ratio,
            maximum_step_ratio,
            xp,
        )
    else:
        step_ratio = float(step_ratio)
        if (
            not np.isfinite(step_ratio)
            or step_ratio < minimum_step_ratio
            or step_ratio > maximum_step_ratio
        ):
            raise ValueError("step_ratio lies outside its configured bounds")

    def update_steps(ratio):
        root = np.sqrt(ratio)
        return base_step * root, base_step / root

    primal_step, dual_step = update_steps(step_ratio)
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
        "primal_step": [],
        "dual_step": [],
        "step_ratio": [],
        "step_adapted": [],
        "relative_primal_change": [],
        "relative_dual_change": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "gram_condition_number": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}
    best = None
    converged = False
    status = "max_iterations"
    message = "maximum number of PDHG iterations reached"
    start_time = time.perf_counter()

    for iteration in range(1, max_iterations + 1):
        extrapolated_distances = _distance_operator(
            extrapolated_gram,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        next_dual = (
            dual
            + dual_step * (extrapolated_distances - target)
        ) / (1.0 + dual_step * variance)
        dual_adjoint = _distance_adjoint(
            next_dual,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        next_gram, next_inverse, next_logdet, next_eigenvalues = (
            _prox_centered_negative_logdet(
                gram - primal_step * dual_adjoint,
                primal_step,
                householder,
                xp,
            )
        )
        next_extrapolated = next_gram + theta * (next_gram - gram)

        primal_change = _norm(next_gram - gram, xp) / max(
            _norm(gram, xp), np.finfo(np.float64).tiny
        )
        dual_change = _norm(next_dual - dual, xp) / max(
            _norm(dual, xp), 1.0
        )
        gram = next_gram
        inverse_gram = next_inverse
        logdet = next_logdet
        internal_eigenvalues = next_eigenvalues
        dual = next_dual
        extrapolated_gram = next_extrapolated

        residuals = _relative_kkt_residuals(
            gram,
            inverse_gram,
            dual,
            target,
            variance,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
        )
        (
            primal_norm,
            primal_scale,
            primal_relative,
            dual_norm,
            dual_scale,
            dual_relative,
            fitted_normalized,
        ) = residuals
        eliminated_residuals = _eliminated_kkt_residuals(
            gram,
            inverse_gram,
            target,
            variance,
            n,
            solver_pair_i,
            solver_pair_j,
            xp,
            fitted=fitted_normalized,
        )
        (
            eliminated_norm,
            eliminated_scale,
            eliminated_relative,
            _,
            _,
        ) = eliminated_residuals
        stationarity_score = eliminated_relative
        data_residual = fitted_normalized - target
        data_objective = 0.5 * _scalar(
            xp.sum(xp.square(data_residual) / variance)
        )
        negative_entropy = (
            -_ENTROPY_COEFFICIENT
            * (logdet + n_modes * np.log(distance_scale))
        )
        objective = negative_entropy + data_objective
        relative_loss = _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual / target)))
        )
        weighted_rmse = _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual) / variance))
        )
        distance_rmse = distance_scale * _scalar(
            xp.sqrt(xp.mean(xp.square(data_residual)))
        )
        entropy = (
            logdet
            + n_modes * np.log(distance_scale)
            - n_modes * np.log(3.0)
        )
        minimum_eigenvalue = _scalar(internal_eigenvalues[0])
        maximum_eigenvalue = _scalar(internal_eigenvalues[-1])
        condition_number = maximum_eigenvalue / minimum_eigenvalue
        connectivity = _connectivity_from_scaled_inverse(
            inverse_gram, distance_scale, xp
        )
        connectivity_diagonal = xp.diag(connectivity)
        offdiagonal_norm_squared = 0.5 * max(
            _scalar(xp.sum(connectivity * connectivity))
            - _scalar(xp.sum(connectivity_diagonal * connectivity_diagonal)),
            0.0,
        )
        connectivity_norm = float(np.sqrt(offdiagonal_norm_squared))

        used_primal_step = primal_step
        used_dual_step = dual_step
        used_step_ratio = step_ratio
        adapted = False

        history["iteration"].append(iteration)
        history["objective"].append(objective)
        history["negative_entropy_objective"].append(negative_entropy)
        history["data_objective"].append(data_objective)
        history["loss"].append(relative_loss)
        history["entropy"].append(entropy)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["relative_primal_kkt_residual"].append(primal_relative)
        history["relative_dual_kkt_residual"].append(dual_relative)
        history["relative_eliminated_kkt_residual"].append(
            eliminated_relative
        )
        history["relative_gradient_norm"].append(stationarity_score)
        history["relative_stationarity_residual"].append(
            stationarity_score
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
        history["connectivity_offdiagonal_l2"].append(
            connectivity_norm
        )

        if adaptive_steps and iteration % adaptation_interval == 0:
            if primal_relative > adaptation_threshold * dual_relative:
                step_ratio = min(
                    step_ratio * adaptation_factor, maximum_step_ratio
                )
                adapted = step_ratio != used_step_ratio
            elif dual_relative > adaptation_threshold * primal_relative:
                step_ratio = max(
                    step_ratio / adaptation_factor, minimum_step_ratio
                )
                adapted = step_ratio != used_step_ratio
            if adapted:
                primal_step, dual_step = update_steps(step_ratio)
                extrapolated_gram = gram.copy()
        history["step_adapted"][-1] = adapted

        if best is None or stationarity_score < best[0]:
            best = (
                stationarity_score,
                gram.copy(),
                inverse_gram.copy(),
                dual.copy(),
                logdet,
                internal_eigenvalues.copy(),
                iteration,
            )

        if iteration in save_steps_set:
            checkpoint = (
                _num.cp.asnumpy(connectivity)
                if use_gpu
                else np.asarray(connectivity).copy()
            )
            connectivity_at_steps[iteration] = checkpoint

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
                    "relative_primal_kkt_residual": primal_relative,
                    "relative_dual_kkt_residual": dual_relative,
                    "relative_eliminated_kkt_residual": (
                        eliminated_relative
                    ),
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

        primal_converged = primal_norm <= (
            absolute_tolerance + relative_tolerance * primal_scale
        )
        dual_converged = dual_norm <= (
            absolute_tolerance + relative_tolerance * dual_scale
        )
        eliminated_converged = eliminated_norm <= (
            absolute_tolerance + relative_tolerance * eliminated_scale
        )
        if primal_converged and dual_converged and eliminated_converged:
            converged = True
            status = "optimality_tolerance"
            message = "PDHG primal, dual, and eliminated KKT tolerances reached"
            break

    if return_best and best is not None:
        (
            _,
            gram,
            inverse_gram,
            dual,
            logdet,
            internal_eigenvalues,
            best_iteration,
        ) = best
    else:
        best_iteration = iteration

    full_fitted = _num._squared_distances_from_gram(
        distance_scale * gram, array_module=xp
    )
    connectivity = _connectivity_from_scaled_inverse(
        inverse_gram, distance_scale, xp
    )
    full_gram = distance_scale * gram

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

    if use_gpu:
        final_primal_relative = float(
            history["relative_primal_kkt_residual"][best_iteration - 1]
        )
        final_dual_relative = float(
            history["relative_dual_kkt_residual"][best_iteration - 1]
        )
    else:
        final_residuals = _relative_kkt_residuals(
            gram,
            inverse_gram,
            dual,
            target_pairs / distance_scale,
            pair_variance / (distance_scale**2),
            n,
            pair_i,
            pair_j,
            np,
        )
        final_primal_relative = float(final_residuals[2])
        final_dual_relative = float(final_residuals[5])

    independent_residuals = _independent_eliminated_kkt_residuals(
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
    ) = independent_residuals
    runtime_eliminated_relative = float(
        history["relative_eliminated_kkt_residual"][best_iteration - 1]
    )
    internal_kkt_converged = bool(converged)
    if relative_tolerance > 0.0:
        independent_kkt_converged = (
            independent_relative <= relative_tolerance
        )
    else:
        independent_kkt_converged = independent_norm <= absolute_tolerance
    converged = internal_kkt_converged and independent_kkt_converged
    if internal_kkt_converged and not independent_kkt_converged:
        status = "independent_kkt_failed"
        message = (
            "PDHG internal tolerances were reached, but the returned Gram "
            "matrix failed the independent eliminated KKT certificate"
        )

    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": int(iteration),
        "returned_iteration": int(best_iteration),
        "objective": float(history["objective"][best_iteration - 1]),
        "relative_gradient_norm": float(independent_relative),
        "relative_primal_kkt_residual": final_primal_relative,
        "relative_dual_kkt_residual": final_dual_relative,
        "relative_eliminated_kkt_residual": float(independent_relative),
        "runtime_relative_eliminated_kkt_residual": (
            runtime_eliminated_relative
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
            np.max(np.abs(pair_stationarity))
        ),
        "termination_internal_kkt_converged": internal_kkt_converged,
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
        "operator_norm_squared_bound": float(
            operator_norm_squared_bound
        ),
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "pdhg": {
            "theta": float(theta),
            "step_safety": float(step_safety),
            "initial_step_ratio": float(history["step_ratio"][0]),
            "final_step_ratio": float(step_ratio),
            "adaptive_steps": bool(adaptive_steps),
            "adaptation_interval": int(adaptation_interval),
            "adaptation_threshold": float(adaptation_threshold),
            "adaptation_factor": float(adaptation_factor),
            "return_best": bool(return_best),
        },
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
        "logged_metric_timing": (
            "post-update PDHG iterate; final certificate recomputed from "
            "returned Gram matrix"
        ),
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_covariance_pdhg stopped without satisfying "
            f"the KKT tolerance (status={status})",
            RuntimeWarning,
            stacklevel=2,
        )
    return full_fitted, full_gram, connectivity, info


def _combine_hybrid_histories(pdhg_history, newton_history):
    """Combine phase histories on one global accepted-update axis."""
    pdhg_count = len(pdhg_history["iteration"])
    newton_count = len(newton_history["iteration"])
    combined = {
        "iteration": np.arange(
            1, pdhg_count + newton_count + 1, dtype=np.int64
        ),
        "phase": np.asarray(
            ["pdhg"] * pdhg_count + ["newton"] * newton_count
        ),
        "phase_iteration": np.concatenate(
            (
                np.asarray(pdhg_history["iteration"], dtype=np.int64),
                np.asarray(newton_history["iteration"], dtype=np.int64),
            )
        ),
    }
    history_keys = (
        set(pdhg_history) | set(newton_history)
    ) - {"iteration"}
    for key in sorted(history_keys):
        pdhg_values = pdhg_history.get(key)
        newton_values = newton_history.get(key)
        if key == "relative_eliminated_kkt_residual" and newton_values is None:
            newton_values = newton_history.get("relative_gradient_norm")
        if pdhg_values is None:
            pdhg_values = np.full(pdhg_count, np.nan)
        if newton_values is None:
            newton_values = np.full(newton_count, np.nan)
        combined[key] = np.concatenate(
            (np.asarray(pdhg_values), np.asarray(newton_values))
        )
    return combined


def fit_gaussian_noise_covariance_hybrid(
    squared_distances,
    noise_variance=None,
    *,
    relative_noise_std=None,
    initialization="rouse",
    initial_connectivity=None,
    use_gpu=False,
    max_iterations=1000,
    relative_tolerance=1e-5,
    absolute_tolerance=1e-10,
    handoff_relative_tolerance=1e-3,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
):
    """Fit COV globally with PDHG and finish locally with Newton-CG.

    PDHG starts from the requested physical initialization and runs until its
    independently certified KKT residual reaches
    ``handoff_relative_tolerance``. The returned PDHG connectivity then
    initializes the existing centered Newton-CG solver. ``max_iterations`` is
    a single total update budget shared by both phases. The final returned
    Gram matrix is independently certified against ``relative_tolerance``.
    """
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if (
        isinstance(handoff_relative_tolerance, (bool, np.bool_))
        or not np.isscalar(handoff_relative_tolerance)
    ):
        raise ValueError(
            "handoff_relative_tolerance must be a positive finite scalar"
        )
    handoff_relative_tolerance = float(handoff_relative_tolerance)
    if (
        not np.isfinite(handoff_relative_tolerance)
        or handoff_relative_tolerance <= 0.0
    ):
        raise ValueError(
            "handoff_relative_tolerance must be a positive finite scalar"
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
        initialization=initialization,
        initial_connectivity=initial_connectivity,
        use_gpu=use_gpu,
        max_iterations=max_iterations,
        relative_tolerance=handoff_relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        initial_gram_floor_relative=initial_gram_floor_relative,
        save_steps=sorted(save_steps_set),
        progress_callback=phase_progress_callback("pdhg", 0),
    )
    pdhg_wall_seconds = time.perf_counter() - pdhg_start_time
    fitted, gram, connectivity, pdhg_info = pdhg_result
    pdhg_iterations = int(pdhg_info["iterations"])
    handoff = {
        "reached": bool(pdhg_info["converged"]),
        "relative_tolerance": handoff_relative_tolerance,
        "relative_eliminated_kkt_residual": float(
            pdhg_info["relative_eliminated_kkt_residual"]
        ),
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
            "PDHG did not reach the hybrid handoff tolerance"
            if not handoff["reached"]
            else "PDHG reached the handoff at the end of the iteration budget"
        )
        info = dict(pdhg_info)
        info.update(
            {
                "algorithm": "hybrid",
                "coordinate_parameterization": "centered_hybrid",
                "converged": False,
                "status": status,
                "message": message,
                "relative_tolerance": float(relative_tolerance),
                "handoff_relative_tolerance": (
                    handoff_relative_tolerance
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
            "fit_gaussian_noise_covariance_hybrid stopped before Newton "
            f"refinement (status={status})",
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
    newton_result = _num.fit_gaussian_noise_covariance(
        squared_distances,
        noise_variance,
        relative_noise_std=relative_noise_std,
        initialization="rouse",
        initial_connectivity=connectivity,
        use_gpu=use_gpu,
        max_iterations=remaining_iterations,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        initial_gram_floor_relative=initial_gram_floor_relative,
        save_steps=newton_save_steps,
        progress_callback=phase_progress_callback(
            "newton", pdhg_iterations
        ),
    )
    newton_wall_seconds = time.perf_counter() - newton_start_time
    fitted, gram, connectivity, newton_info = newton_result
    newton_iterations = int(newton_info["iterations"])

    (
        _,
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
    independent_residuals = _independent_eliminated_kkt_residuals(
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
    ) = independent_residuals
    if relative_tolerance > 0.0:
        independent_kkt_converged = (
            independent_relative <= relative_tolerance
        )
    else:
        independent_kkt_converged = (
            independent_norm <= absolute_tolerance
        )
    internal_kkt_converged = bool(newton_info["converged"])
    converged = internal_kkt_converged and independent_kkt_converged
    if converged:
        status = "optimality_tolerance"
        message = (
            "PDHG handoff, Newton refinement, and independent final KKT "
            "certificate reached"
        )
    elif internal_kkt_converged:
        status = "independent_kkt_failed"
        message = (
            "Newton reached its internal tolerance, but the returned Gram "
            "matrix failed the independent eliminated KKT certificate"
        )
    else:
        status = f"newton_{newton_info['status']}"
        message = f"Newton refinement stopped: {newton_info['message']}"

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
            "runtime_relative_eliminated_kkt_residual": float(
                newton_info["relative_gradient_norm"]
            ),
            "stationarity_residual_norm": float(independent_norm),
            "stationarity_residual_scale": float(independent_scale),
            "relative_stationarity_residual": float(independent_relative),
            "maximum_absolute_stationarity_residual": float(
                np.max(np.abs(independent_residual))
            ),
            "pair_stationarity_residual_norm": float(
                newton_info["stationarity_residual_norm"]
            ),
            "relative_pair_stationarity_residual": float(
                newton_info["relative_stationarity_residual"]
            ),
            "maximum_absolute_pair_stationarity_residual": float(
                newton_info["maximum_absolute_stationarity_residual"]
            ),
            "termination_internal_kkt_converged": internal_kkt_converged,
            "independent_kkt_converged": bool(independent_kkt_converged),
            "independent_kkt_recomputed_from_returned_gram": True,
            "relative_tolerance": float(relative_tolerance),
            "absolute_tolerance": float(absolute_tolerance),
            "handoff_relative_tolerance": handoff_relative_tolerance,
            "handoff": handoff,
            "phase_iterations": {
                "pdhg": pdhg_iterations,
                "newton": newton_iterations,
            },
            "phase_wall_seconds": {
                "pdhg": pdhg_wall_seconds,
                "newton": newton_wall_seconds,
            },
            "initialization": pdhg_info["initialization"],
            "algorithm": "hybrid",
            "coordinate_parameterization": "centered_hybrid",
            "history": combined_history,
            "connectivity_matrix_at_steps": connectivity_at_steps,
            "wall_seconds": time.perf_counter() - start_time,
            "stationarity_definition": pdhg_info[
                "stationarity_definition"
            ],
            "pair_stationarity_definition": pdhg_info[
                "pair_stationarity_definition"
            ],
            "logged_metric_timing": (
                "post-update PDHG then Newton iterates on one global axis; "
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
