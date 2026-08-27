"""Efficient covariance-coordinate operators for noise-aware HIPPS-DIMES."""

from __future__ import annotations

import numpy as np


def symmetrize(matrix):
    return 0.5 * (matrix + matrix.T)


def anchored_basis(n: int) -> np.ndarray:
    """Return the gauge basis in which the final locus is fixed at the origin."""
    basis = np.zeros((n, n - 1), dtype=np.float64)
    basis[: n - 1] = np.eye(n - 1, dtype=np.float64)
    return basis


def centered_to_anchored_gram(centered_gram: np.ndarray) -> np.ndarray:
    """Convert a centered full Gram matrix to final-locus relative coordinates."""
    centered_gram = np.asarray(centered_gram, dtype=np.float64)
    anchor = centered_gram.shape[0] - 1
    anchored = (
        centered_gram[:anchor, :anchor]
        - centered_gram[:anchor, anchor, None]
        - centered_gram[anchor, :anchor][None, :]
        + centered_gram[anchor, anchor]
    )
    return symmetrize(anchored)


def anchored_to_centered_gram(anchored_gram, array_module=np):
    """Convert a final-locus anchored Gram matrix to a centered full Gram."""
    m = anchored_gram.shape[0]
    full = array_module.zeros((m + 1, m + 1), dtype=anchored_gram.dtype)
    full[:m, :m] = anchored_gram
    row_mean = array_module.mean(full, axis=1, keepdims=True)
    return symmetrize(full - row_mean - row_mean.T + array_module.mean(full))


def anchored_squared_distances(anchored_gram, array_module=np):
    """Return all squared distances induced by an anchored Gram matrix."""
    m = anchored_gram.shape[0]
    diagonal = array_module.diag(anchored_gram)
    distances = array_module.zeros((m + 1, m + 1), dtype=anchored_gram.dtype)
    distances[:m, :m] = (
        diagonal[:, None] + diagonal[None, :] - 2.0 * anchored_gram
    )
    distances[:m, m] = diagonal
    distances[m, :m] = diagonal
    array_module.fill_diagonal(distances, 0.0)
    return symmetrize(distances)


def anchored_data_objective_gradient(
    anchored_gram,
    target_matrix,
    inverse_variance_matrix,
    pair_i,
    pair_j,
    *,
    array_module=np,
):
    """Evaluate weighted distance least squares in anchor-relative coordinates."""
    fitted = anchored_squared_distances(anchored_gram, array_module)
    residual = fitted - target_matrix
    weighted = inverse_variance_matrix * residual
    objective = 0.25 * array_module.sum(weighted * residual)
    m = anchored_gram.shape[0]
    gradient = (
        array_module.diag(array_module.sum(weighted[:m, :], axis=1))
        - weighted[:m, :m]
    )
    return objective, symmetrize(gradient), fitted[pair_i, pair_j]


def anchored_data_hessian_action(
    direction, inverse_variance_matrix, *, array_module=np
):
    """Apply the constant weighted distance Hessian in O(N^2) work."""
    distance_direction = anchored_squared_distances(
        symmetrize(direction), array_module
    )
    weighted = inverse_variance_matrix * distance_direction
    m = direction.shape[0]
    action = (
        array_module.diag(array_module.sum(weighted[:m, :], axis=1))
        - weighted[:m, :m]
    )
    return symmetrize(action)


def anchored_data_svec_diagonal(inverse_variance_matrix, array_module=np):
    """Exact data-Hessian diagonal in the orthonormal symmetric svec basis."""
    m = inverse_variance_matrix.shape[0] - 1
    diagonal = 2.0 * inverse_variance_matrix[:m, :m].copy()
    diagonal[array_module.diag_indices(m)] = array_module.sum(
        inverse_variance_matrix[:m, :], axis=1
    )
    return symmetrize(diagonal)


def entropy_svec_diagonal(inverse_gram, array_module=np):
    """Exact svec diagonal of H -> 3/2 B^{-1} H B^{-1}."""
    diagonal = array_module.diag(inverse_gram)
    result = 1.5 * (
        array_module.outer(diagonal, diagonal)
        + array_module.square(inverse_gram)
    )
    result[array_module.diag_indices(inverse_gram.shape[0])] = (
        1.5 * array_module.square(diagonal)
    )
    return symmetrize(result)


def full_diagonal_to_svec(full_coordinate_diagonal, array_module=np):
    """Convert a full-matrix Hessian diagonal to normalized svec scaling."""
    result = 2.0 * full_coordinate_diagonal.copy()
    indices = array_module.diag_indices(result.shape[0])
    result[indices] = full_coordinate_diagonal[indices]
    return symmetrize(result)


def anchored_connectivity_from_inverse(inverse_gram, array_module=np):
    """Return the full row-sum-zero connectivity from inverse anchored Gram."""
    m = inverse_gram.shape[0]
    precision = 3.0 * inverse_gram
    row_sum = array_module.sum(precision, axis=1)
    connectivity = array_module.empty((m + 1, m + 1), dtype=inverse_gram.dtype)
    connectivity[:m, :m] = -precision
    connectivity[:m, m] = row_sum
    connectivity[m, :m] = row_sum
    connectivity[m, m] = -array_module.sum(row_sum)
    return symmetrize(connectivity)


def anchored_connectivity_from_gram(anchored_gram: np.ndarray) -> np.ndarray:
    """Recover the full connectivity from an anchored Gram matrix."""
    anchored_gram = np.asarray(anchored_gram, dtype=np.float64)
    inverse = np.linalg.solve(
        anchored_gram, np.eye(anchored_gram.shape[0], dtype=np.float64)
    )
    connectivity = anchored_connectivity_from_inverse(inverse, np)
    np.fill_diagonal(connectivity, 0.0)
    np.fill_diagonal(connectivity, -np.sum(connectivity, axis=1))
    return symmetrize(connectivity)
