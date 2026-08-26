"""
Reconstruction of 3D genome organization using the Maximum Entropy Principle

Reference:
1. Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051.
https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051
2. Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants." Nature Communications 14.1 (2023): 1150.
3. Shi, Guang, Sucheol, Shin, and D. Thirumalai. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advance 11 (31), eadx1763
"""

import os
import sys
import time
import warnings
import itertools

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import numpy as np
import scipy
import scipy.linalg
import scipy.interpolate
import scipy.optimize
import scipy.spatial.distance
import pandas as pd

try:
    import rich_click as click
except ImportError:
    import click  # type: ignore[no-redef]
try:
    import cooler  # for cooler format of HiC data
except Exception:
    cooler = None
try:
    import hicstraw  # for .hic format of HiC data
except Exception:
    hicstraw = None
from rich import print
from rich.panel import Panel
from rich.text import Text
from rich.console import Console
from rich.table import Table

# from tqdm.rich import trange, tqdm
from tqdm import trange, tqdm


console = Console()

# ------------------------------------------------------------------#
# GPU Support (CuPy)
# ------------------------------------------------------------------#
# Check if CuPy is available for GPU acceleration
_CUPY_AVAILABLE = False
_CUPY_GPU_NAME = None
try:
    import cupy as cp
    import cupyx
    import cupyx.scipy.linalg as cupyx_scipy_linalg
    import cupyx.scipy.sparse.linalg as cupyx_sparse_linalg

    # Test if GPU is actually accessible
    cp.cuda.runtime.getDeviceCount()
    _CUPY_AVAILABLE = True
    try:
        gpu_props = cp.cuda.runtime.getDeviceProperties(0)
        _CUPY_GPU_NAME = (
            gpu_props["name"].decode()
            if isinstance(gpu_props["name"], bytes)
            else gpu_props["name"]
        )
    except:
        _CUPY_GPU_NAME = "Unknown GPU"
except ImportError:
    pass
except Exception:
    pass


def is_gpu_available():
    """Check if GPU acceleration is available via CuPy.

    Returns
    -------
    bool
        True if CuPy is installed and a GPU is accessible.
    """
    return _CUPY_AVAILABLE


def get_gpu_name():
    """Get the name of the available GPU.

    Returns
    -------
    str or None
        GPU name if available, None otherwise.
    """
    return _CUPY_GPU_NAME


def _a2dmap_theory_gpu(A_gpu, force_positive_definite=False, return_eigenvalues=False):
    """
    GPU version of a2dmap_theory using CuPy.

    Parameters
    ----------
    A_gpu : cupy.ndarray
        Connectivity matrix on GPU
    force_positive_definite : bool
        If True, replace negative temp values with zero
    return_eigenvalues : bool
        If True, also return the eigenvalues

    Returns
    -------
    dmap : cupy.ndarray
        Mean distance map on GPU
    eigvalue : cupy.ndarray, optional
        Eigenvalues (only if return_eigenvalues=True)
    """
    TOL = 10**8
    eigvalue, eigvector = cp.linalg.eigh(A_gpu)

    temp = -1.0 / eigvalue
    temp = cp.where(cp.isinf(temp), 0.0, temp)
    temp = cp.where(temp >= TOL, 0.0, temp)
    temp = cp.where(temp <= -TOL, 0.0, temp)

    if force_positive_definite:
        temp = cp.where(temp < 0.0, 0.0, temp)

    # Avoid materializing a dense diagonal matrix:
    # V @ diag(temp) @ V.T == (V * temp) @ V.T  (broadcast scales columns of V)
    Omega = (eigvector * temp) @ eigvector.T
    Omega_diag = cp.diag(Omega)
    sigma = cp.sqrt(Omega_diag[:, cp.newaxis] + Omega_diag - 2.0 * Omega)
    dmap = 2.0 * cp.sqrt(2.0 / cp.pi) * sigma

    if return_eigenvalues:
        return dmap, eigvalue
    return dmap


# ------------------------------------------------------------------#
# Thread control for eigenvalue (eigh) and other BLAS/LAPACK operations
def set_eigh_num_threads(n):
    """
    Set the number of threads used by BLAS/LAPACK for eigenvalue and linear algebra.
    Must be called before any np.linalg.eigh (or similar) to take effect.
    Sets: OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, OMP_NUM_THREADS, VECLIB_MAXIMUM_THREADS.
    """
    s = str(int(n))
    os.environ["OPENBLAS_NUM_THREADS"] = s
    os.environ["MKL_NUM_THREADS"] = s
    os.environ["OMP_NUM_THREADS"] = s
    os.environ["VECLIB_MAXIMUM_THREADS"] = s


# ------------------------------------------------------------------#
# Helper functions
def restore_matrix_with_nans(small, removed_idx, original_size):
    """Restore a reduced square matrix to its original size, filling removed loci with NaN.

    Parameters
    ----------
    small : array_like
        Square matrix after removing fully missing loci.
    removed_idx : array_like
        Zero-based indices removed from the original matrix.
    original_size : int
        Size of the original square matrix.

    Returns
    -------
    np.ndarray
        ``(original_size, original_size)`` float array with the reduced matrix
        inserted at the kept loci positions and NaN-filled rows/columns for the
        removed loci.
    """
    small = np.asarray(small)
    if small.ndim != 2 or small.shape[0] != small.shape[1]:
        raise ValueError(f"small must be a square 2D array, got shape {small.shape}")

    original_size = int(original_size)
    if original_size < 0:
        raise ValueError("original_size must be non-negative")

    removed_idx = np.asarray(removed_idx, dtype=int)
    if removed_idx.ndim != 1:
        raise ValueError("removed_idx must be a 1D array-like of indices")
    if removed_idx.size and np.unique(removed_idx).size != removed_idx.size:
        raise ValueError("removed_idx must not contain duplicate indices")
    if np.any(removed_idx < 0) or np.any(removed_idx >= original_size):
        raise ValueError(
            f"removed_idx entries must be in [0, {max(original_size - 1, 0)}]"
        )

    expected_size = original_size - removed_idx.size
    if small.shape != (expected_size, expected_size):
        raise ValueError(
            "small has incompatible shape for the provided removed_idx/original_size: "
            f"expected {(expected_size, expected_size)}, got {small.shape}"
        )

    keep_mask = np.ones(original_size, dtype=bool)
    keep_mask[removed_idx] = False

    out_dtype = np.result_type(small.dtype, float)
    restored = np.full((original_size, original_size), np.nan, dtype=out_dtype)
    restored[np.ix_(keep_mask, keep_mask)] = small
    return restored


def compute_acf_general_theory(i, j, t, a, zeta=1.0):
    """
    Numerically compute the autocorrelation function (ACF) for monomers i, j
    using the connectivity matrix `a`. Returns both the time-dependent
    ACF and the corresponding MSD for each time point.

    Parameters
    ----------
    i, j : int
        Indices of the monomers for which to calculate the correlation function.
    t : array_like
        A 1D array of time points (lag times).
    a : np.ndarray
        The connectivity (or "Laplacian") matrix for the polymer/chain.
    zeta : float
        Friction coefficient (if part of the model). Default is 1.0.

    Returns
    -------
    two_point_acf : np.ndarray
        2D array. First column is time `t`, second column is the ACF values at each time.
    two_point_msd : np.ndarray
        2D array. First column is time `t`, second column is the two-point MSD at each time
    """
    eigvalue, eigvector = np.linalg.eigh(a)
    eigvalue_inv = 1.0 / eigvalue

    # difference in eigenvector components for monomers i and j
    vpi_vpj = eigvector[i, :] - eigvector[j, :]

    # normal_modes_square_mean = -(1 / eigenvalue) but filter out any inf
    normal_modes_square_mean = -np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)

    # Expand time dimension for broadcast
    t_reshaped = np.expand_dims(t, axis=-1)
    # Effective relaxation times
    tau_p = -zeta / eigvalue
    decay_factor = np.exp(-t_reshaped / tau_p)

    # ACF(t)
    res = 3.0 * np.sum(vpi_vpj**2 * decay_factor * normal_modes_square_mean, axis=-1)
    # Equilibrium part
    res_eq = 3.0 * np.sum(vpi_vpj**2 * normal_modes_square_mean, axis=-1)

    # Combine results: first column is time, second is res
    two_point_acf = np.column_stack((t, res))

    # The mean-square displacement from the correlation function
    two_point_msd = 2.0 * (res_eq - two_point_acf[:, 1])
    two_point_msd = np.column_stack((t, two_point_msd))

    return two_point_acf, two_point_msd


def compute_m1_i(i, t, a, zeta=1.0):
    """
    Compute the single-monomer mean-square displacement (MSD) for monomer i,
    given the connectivity matrix `a`.

    Parameters
    ----------
    i : int
        Index of the monomer.
    t : array_like
        A 1D array of time points (lag times).
    a : np.ndarray
        The connectivity (or "Laplacian") matrix for the polymer/chain.
    zeta : float
        Friction coefficient. Default is 1.0.

    Returns
    -------
    msd : np.ndarray
        2D array. First column is time `t`, second column is the MSD for
        monomer i at those times.
    """
    t = np.asarray(t)
    n = a.shape[0]

    lam, V = np.linalg.eigh(a)

    # Identify and remove center-of-mass (near-zero) mode.
    p0 = np.argmin(np.abs(lam))
    p_nz = np.arange(n) != p0
    lam_nz = lam[p_nz]
    V_nz = V[:, p_nz]

    sigma2_nz = -1.0 / lam_nz
    tau_nz = -zeta / lam_nz

    vpi2 = V_nz[i, :] ** 2
    t_col = t[:, None]
    decay = np.exp(-t_col / tau_nz[None, :])

    # Internal mode contribution (excluding CM mode).
    res = 3.0 * np.sum(vpi2[None, :] * decay * sigma2_nz[None, :], axis=1)
    r2_eq = 3.0 * np.sum(vpi2 * sigma2_nz)

    # Center-of-mass diffusion contribution.
    D_cm = 1.0 / (zeta * n)
    msd_cm = 6.0 * D_cm * t

    # Total MSD = CM diffusion + internal motion
    msd_data = msd_cm + 2.0 * (r2_eq - res)

    # Combine time with MSD
    msd = np.column_stack((t, msd_data))
    return msd


def compute_m1_all(a, t, zeta=1.0, tol=1e-12):
    """
    Compute single-monomer MSD for all monomers, including center-of-mass motion.
    The zero (CM) mode is handled analytically rather than dropped.

    Parameters
    ----------
    a : np.ndarray, shape (n, n)
        Connectivity (Laplacian) matrix, negative-semidefinite.
    t : array_like, shape (len(t),)
        A 1D array of lag times.
    zeta : float, optional
        Friction coefficient. Default is 1.0.
    tol : float, optional
        Tolerance to identify zero eigenvalue. Default is 1e-12.

    Returns
    -------
    msd_all : np.ndarray, shape (n, len(t), 2)
        For each monomer m=0..n-1 and each lag time t[k], returns [t[k], MSD(m, t[k])].
        MSD includes both internal motion and center-of-mass diffusion.

    Notes
    -----
    The MSD calculation includes:
    1. Internal motion from non-zero modes
    2. Center-of-mass diffusion (6D_cm*t where D_cm = 1/(zeta*n))
    """
    # Input validation
    if not isinstance(a, np.ndarray) or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")
    if not isinstance(zeta, (int, float)) or zeta <= 0:
        raise ValueError("zeta must be a positive number")
    if not isinstance(tol, (int, float)) or tol <= 0:
        raise ValueError("tol must be a positive number")

    t = np.asarray(t)
    if t.ndim != 1:
        raise ValueError("t must be a 1D array")

    n = a.shape[0]

    # Eigendecomposition
    lam, V = np.linalg.eigh(a)  # lam: eigenvalues, V: eigenvectors

    # Find zero (CM) mode
    p0 = np.argmin(np.abs(lam))
    if abs(lam[p0]) > tol:
        raise ValueError(
            "No near-zero eigenvalue found; CM mode missing or 'a' not Rouse-type"
        )

    # Get non-zero modes
    p_nz = np.arange(n) != p0  # boolean mask for non-zero modes
    lam_nz = lam[p_nz]  # non-zero eigenvalues
    V_nz = V[:, p_nz]  # non-zero eigenvectors

    # Compute equilibrium variances for non-zero modes
    sigma2_nz = -1.0 / lam_nz  # shape (n-1,)

    # Compute relaxation times
    tau = np.full(n, np.inf)  # initialize all to infinity
    tau[p_nz] = -zeta / lam_nz  # set non-zero mode times

    # Compute decay factors for non-zero modes
    t_col = t[:, None]  # shape (len(t), 1)
    decay = np.exp(-t_col / tau[p_nz])  # shape (len(t), n-1)

    # Compute squared eigenvector components
    V2 = V_nz**2  # shape (n, n-1)

    # Compute time-dependent part from non-zero modes
    res = 3.0 * np.sum(
        V2[:, None, :] * decay[None, :, :] * sigma2_nz[None, None, :], axis=2
    )  # shape (n, len(t))

    # Compute equilibrium variance for each monomer
    r2_eq = 3.0 * np.sum(V2 * sigma2_nz[None, :], axis=1)  # shape (n,)

    # Compute center-of-mass diffusion
    D_cm = 1.0 / (zeta * n)
    msd_cm = 6.0 * D_cm * t  # shape (len(t),)

    # Total MSD = CM diffusion + internal motion
    msd = msd_cm[None, :] + 2.0 * (r2_eq[:, None] - res)  # shape (n, len(t))

    # Stack time with MSD
    msd_all = np.stack([np.tile(t, (n, 1)), msd], axis=-1)

    return msd_all


def compute_mcom_segment(i, j, t, a, zeta=1.0):
    """
    Compute the center-of-mass mean-square displacement (MSD) for the subsegment
    running from monomer i to monomer j (inclusive), given the connectivity matrix `a`.

    Parameters
    ----------
    i : int
        Index of the first monomer in the segment.
    j : int
        Index of the last monomer in the segment (inclusive).
    t : array_like
        A 1D array of time points (lag times).
    a : np.ndarray
        The connectivity (or "Laplacian") matrix for the polymer/chain.
    zeta : float
        Friction coefficient. Default is 1.0.

    Returns
    -------
    msd : np.ndarray
        2D array. First column is time `t`, second column is the MSD for the
        segment COM at those times.
    """
    # Input validation
    if not isinstance(i, int) or not isinstance(j, int):
        raise TypeError("i and j must be integers")
    if i < 0 or j < i or j >= a.shape[0]:
        raise ValueError("Invalid indices i and j")
    if not isinstance(zeta, (int, float)) or zeta <= 0:
        raise ValueError("zeta must be a positive number")
    if not isinstance(t, np.ndarray) or t.ndim != 1:
        raise ValueError("t must be a 1D array")
    if not isinstance(a, np.ndarray) or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")

    # Eigen-decomposition of the connectivity matrix
    eigvalue, eigvector = np.linalg.eigh(a)

    # Discard zero mode
    eigvalue = eigvalue[1:]
    eigvector = eigvector[:, 1:]

    # Compute weights W_p = (1/(j-i+1)) * sum_{m=i}^j V_{p,m}
    segment_length = j - i + 1
    W = np.sum(eigvector[i : j + 1, :], axis=0) / segment_length

    # 1/eigvalue, filtering infinities (zero-mode)
    eigvalue_inv = 1.0 / eigvalue
    normal_modes_square_mean = -np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)

    # Expand time dimension for broadcasting
    t_reshaped = np.expand_dims(t, axis=-1)  # shape (len(t), 1)

    # Relaxation times tau_p
    tau_p = -zeta / eigvalue
    tau_p = np.nan_to_num(
        tau_p, posinf=0.0, neginf=0.0
    )  # Handle potential division by zero

    # Decay factor for each mode at each time
    decay_factor = np.exp(-t_reshaped / tau_p)  # shape (len(t), n)

    # Compute time-dependent part: 3 * sum(W^2 * decay_factor * normal_modes_square_mean, over modes)
    # Ensure proper broadcasting by reshaping W
    W_reshaped = W[None, :]  # shape (1, n)
    res = 3.0 * np.sum(
        (W_reshaped**2) * decay_factor * normal_modes_square_mean, axis=-1
    )

    # Equilibrium variance of COM: 3 * sum(W^2 * normal_modes_square_mean)
    r2_eq = 3.0 * np.sum((W**2) * normal_modes_square_mean)

    # MSD = 2 * (r2_eq - res)
    msd_data = 2.0 * (r2_eq - res)

    # Combine time with MSD into a 2D array
    msd = np.column_stack((t, msd_data))
    return msd


def compute_all_tau_1e(K, num_t=200, factor=10.0, tol=1e-8):
    """
    Vectorized computation of 1/e relaxation times tau_{ij} for all pairs (i,j)
    by sampling the normalized correlator on a time grid and linearly interpolating
    to find when G_norm(t) crosses 1/e.

    Parameters
    ----------
    K : (n×n) ndarray
        Connectivity matrix (symmetric).
    num_t : int
        Number of time points in the log-spaced grid.
    factor : float
        Multiplier for raising the upper time bound relative to the slowest-mode time.
    tol : float
        Small epsilon to avoid division by zero.

    Returns
    -------
    tau_mat : (n×n) ndarray
        Approximate 1/e times tau_{ij} for each pair (i,j).
    """
    n = K.shape[0]
    # 1. Diagonalize -K
    eigvals, eigvecs = np.linalg.eigh(-K)
    lam = eigvals[1:]  # drop zero mode
    vec = eigvecs[:, 1:]  # shape (n, m)

    # 2. Precompute E_p(i,j) and W0 = E_p/λ_p
    diff = vec[:, None, :] - vec[None, :, :]  # shape (n,n,m)
    E = diff * diff  # (n,n,m)
    W0 = E / lam[None, None, :]  # (n,n,m)

    # 3. Time grid
    t_min = 0.0
    tau_slow = 1.0 / np.min(lam)
    t_max = factor * tau_slow
    t_vals = np.logspace(np.log10(t_min + 1e-12), np.log10(t_max), num_t)

    # 4. Compute normalized G for each t
    denom = np.sum(W0, axis=2)  # (n,n)
    G_norm = np.empty((n, n, num_t))
    for idx, t in enumerate(tqdm(t_vals)):
        num = np.sum(W0 * np.exp(-lam[None, None, :] * t), axis=2)
        G_norm[..., idx] = num / (denom + tol)

    # 5. Find first crossing index where G_norm <= 1/e
    mask = G_norm <= (1.0 / np.e)
    # argmax returns first True along axis=2, but if no True, gives 0.
    cross_idx = mask.argmax(axis=2)

    # 6. Interpolate tau
    tau_mat = np.zeros((n, n))
    # handle crossing at first bin or no crossing
    tau_mat[mask[..., 0] == False] = t_vals[0]  # if first bin is true or none
    for i in range(n):
        for j in range(n):
            k = cross_idx[i, j]
            if k == 0:
                tau_mat[i, j] = t_vals[0]
            else:
                t0, t1 = t_vals[k - 1], t_vals[k]
                G0, G1 = G_norm[i, j, k - 1], G_norm[i, j, k]
                # linear interpolation
                tau_mat[i, j] = t0 + (1 / np.e - G0) * (t1 - t0) / (G1 - G0 + tol)

    np.fill_diagonal(tau_mat, 0.0)

    return tau_mat


def compute_all_tau_1e_memory_efficient(
    K, num_t=200, factor=10.0, tol=1e-8, block_size=128
):
    """
    Memory-efficient computation of 1/e relaxation times tau_{ij} for all loci pairs.

    This function avoids allocating any (n, n, n) tensor by processing the matrix in
    symmetric blocks and streaming over time points.

    Parameters
    ----------
    K : (n×n) ndarray
        Connectivity matrix (symmetric).
    num_t : int, optional
        Number of time points in the log-spaced grid. Default is 200.
    factor : float, optional
        Multiplier for the upper time bound relative to the slowest-mode time.
        Default is 10.0.
    tol : float, optional
        Small epsilon for numerical stability. Default is 1e-8.
    block_size : int, optional
        Block size for pairwise computation. Smaller values use less memory but can
        be slower. Default is 128.

    Returns
    -------
    tau_mat : (n×n) ndarray
        Approximate 1/e times tau_{ij} for each pair (i,j).
    """
    K = np.asarray(K)
    if K.ndim != 2 or K.shape[0] != K.shape[1]:
        raise ValueError("K must be a square matrix")
    if num_t < 1:
        raise ValueError("num_t must be >= 1")
    if block_size < 1:
        raise ValueError("block_size must be >= 1")

    n = K.shape[0]
    if n == 0:
        return np.empty((0, 0), dtype=float)
    if n == 1:
        return np.zeros((1, 1), dtype=float)

    eigvals, eigvecs = np.linalg.eigh(-K)
    lam = eigvals[1:]  # drop zero mode
    vec = eigvecs[:, 1:]  # shape (n, n-1)

    if lam.size == 0:
        return np.zeros((n, n), dtype=float)

    # Build time grid.
    tau_slow = 1.0 / np.min(lam)
    t_vals = np.logspace(np.log10(1e-12), np.log10(factor * tau_slow), num_t)
    if num_t == 1:
        tau_single = np.full((n, n), t_vals[0], dtype=float)
        np.fill_diagonal(tau_single, 0.0)
        return tau_single

    inv_lam = 1.0 / lam
    target = 1.0 / np.e
    # Precompute mode weights for each sampled time: w_p(t) = exp(-lam_p t) / lam_p.
    weights_t = np.exp(-np.outer(t_vals, lam)) * inv_lam[None, :]

    tau_mat = np.full((n, n), t_vals[-1], dtype=float)

    for i0 in range(0, n, block_size):
        i1 = min(i0 + block_size, n)
        Vi = vec[i0:i1, :]  # (bi, m)
        Vi_sq = Vi * Vi
        Vi_inv = Vi * inv_lam[None, :]
        s0_i = np.sum(Vi_sq * inv_lam[None, :], axis=1)  # (bi,)

        for j0 in range(i0, n, block_size):
            j1 = min(j0 + block_size, n)
            Vj = vec[j0:j1, :]  # (bj, m)
            Vj_sq = Vj * Vj
            s0_j = np.sum(Vj_sq * inv_lam[None, :], axis=1)  # (bj,)

            # Denominator G_ij(0): sum_p ((V_ip - V_jp)^2 / lam_p)
            M0 = Vi_inv @ Vj.T
            denom = s0_i[:, None] + s0_j[None, :] - 2.0 * M0
            denom = np.maximum(denom, tol)

            block_tau = np.full((i1 - i0, j1 - j0), t_vals[-1], dtype=float)
            done = np.zeros_like(block_tau, dtype=bool)

            # t = t_vals[0]
            w_prev = weights_t[0]
            si_prev = np.sum(Vi_sq * w_prev[None, :], axis=1)
            sj_prev = np.sum(Vj_sq * w_prev[None, :], axis=1)
            M_prev = (Vi * w_prev[None, :]) @ Vj.T
            g_prev = (si_prev[:, None] + sj_prev[None, :] - 2.0 * M_prev) / denom

            crossed = g_prev <= target
            block_tau[crossed] = t_vals[0]
            done[crossed] = True

            for t_idx in range(1, num_t):
                if done.all():
                    break

                w_curr = weights_t[t_idx]
                si_curr = np.sum(Vi_sq * w_curr[None, :], axis=1)
                sj_curr = np.sum(Vj_sq * w_curr[None, :], axis=1)
                M_curr = (Vi * w_curr[None, :]) @ Vj.T
                g_curr = (si_curr[:, None] + sj_curr[None, :] - 2.0 * M_curr) / denom

                new_cross = (~done) & (g_curr <= target)
                if np.any(new_cross):
                    g0 = g_prev[new_cross]
                    g1 = g_curr[new_cross]
                    dg = g1 - g0
                    safe_dg = np.where(np.abs(dg) < tol, -tol, dg)
                    frac = np.clip((target - g0) / safe_dg, 0.0, 1.0)
                    t0 = t_vals[t_idx - 1]
                    dt = t_vals[t_idx] - t0
                    block_tau[new_cross] = t0 + frac * dt
                    done[new_cross] = True

                g_prev = g_curr

            tau_mat[i0:i1, j0:j1] = block_tau
            if j0 != i0:
                tau_mat[j0:j1, i0:i1] = block_tau.T

    np.fill_diagonal(tau_mat, 0.0)
    return tau_mat


def compute_all_tau_integral(a):
    """
    Compute first‐moment (tau1) and second‐moment (tau2) relaxation times
    for all pairs of loci (i, j) in a single vectorized pass.

    Parameters
    ----------
    K : (n×n) ndarray
        Connectivity matrix.

    Returns
    -------
    tau1 : (n×n) ndarray
        First‐moment relaxation time for each pair (i,j):
          tau1[i,j] = (∑_p E_p(i,j)/λ_p^2) / (∑_p E_p(i,j)/λ_p)
    tau2 : (n×n) ndarray
        Second‐moment relaxation time ratio for each pair (i,j):
          tau2[i,j] = (∑_p E_p(i,j)/λ_p^3) / (∑_p E_p(i,j)/λ_p^2)
    """
    # 1) Eigendecompose -a
    eigvals, eigvecs = np.linalg.eigh(-a)  # eigvals ascending

    # 2) Discard the zero‐mode
    lam = eigvals[1:]  # shape (n-1,)
    vec = eigvecs[:, 1:]  # shape (n,   n-1)

    # 3) Precompute inverse powers
    inv_lam = 1.0 / lam  # 1/λ_p
    inv_lam2 = inv_lam**2  # 1/λ_p^2
    inv_lam3 = inv_lam**3  # 1/λ_p^3

    # 4) Compute squared differences for all (i,j,p)
    #    E[i,j,p] = (V[p,i] - V[p,j])^2
    #    vec shape -> (n, m), we want (n, n, m)
    diff = vec[:, None, :] - vec[None, :, :]  # shape (n, n, m)
    E = diff * diff  # squared differences

    # 5) Weighted sums via tensordot over the mode axis
    sum0 = np.tensordot(E, inv_lam, axes=([2], [0]))  # ∑ E/λ_p
    sum1 = np.tensordot(E, inv_lam2, axes=([2], [0]))  # ∑ E/λ_p^2
    sum2 = np.tensordot(E, inv_lam3, axes=([2], [0]))  # ∑ E/λ_p^3

    # 6) Form relaxation times
    tau1 = sum1 / sum0  # first moment
    tau2 = sum2 / sum1  # second moment ratio

    np.fill_diagonal(tau1, 0.0)
    np.fill_diagonal(tau2, 0.0)

    return tau1, tau2


def compute_all_tau_integral_memory_efficient(a, block_size=256, tol=1e-14):
    """
    Memory-efficient computation of first/second-moment relaxation times
    for all loci pairs (i, j), without materializing an (n, n, n) tensor.

    Parameters
    ----------
    a : (n×n) ndarray
        Connectivity matrix.
    block_size : int, optional
        Block size used for pairwise computations. Smaller values reduce peak
        memory but can be slower. Default is 256.
    tol : float, optional
        Small threshold for safe division. Default is 1e-14.

    Returns
    -------
    tau1 : (n×n) ndarray
        First-moment relaxation time for each pair (i,j):
          tau1[i,j] = (∑_p E_p(i,j)/λ_p^2) / (∑_p E_p(i,j)/λ_p)
    tau2 : (n×n) ndarray
        Second-moment relaxation time ratio for each pair (i,j):
          tau2[i,j] = (∑_p E_p(i,j)/λ_p^3) / (∑_p E_p(i,j)/λ_p^2)
    """
    a = np.asarray(a)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")
    if block_size < 1:
        raise ValueError("block_size must be >= 1")

    n = a.shape[0]
    if n == 0:
        empty = np.empty((0, 0), dtype=float)
        return empty, empty
    if n == 1:
        zero = np.zeros((1, 1), dtype=float)
        return zero, zero

    # 1) Eigendecompose -a
    eigvals, eigvecs = np.linalg.eigh(-a)

    # 2) Discard the zero mode
    lam = eigvals[1:]  # shape (n-1,)
    vec = eigvecs[:, 1:]  # shape (n, n-1)
    if lam.size == 0:
        zero = np.zeros((n, n), dtype=float)
        return zero, zero

    # 3) Precompute inverse powers and per-row weighted norms.
    inv_lam = 1.0 / lam
    inv_lam2 = inv_lam * inv_lam
    inv_lam3 = inv_lam2 * inv_lam

    vec_sq = vec * vec
    s0 = np.sum(vec_sq * inv_lam[None, :], axis=1)
    s1 = np.sum(vec_sq * inv_lam2[None, :], axis=1)
    s2 = np.sum(vec_sq * inv_lam3[None, :], axis=1)

    # Weighted eigenvector matrices used in block GEMMs.
    vec_w0 = vec * inv_lam[None, :]
    vec_w1 = vec * inv_lam2[None, :]
    vec_w2 = vec * inv_lam3[None, :]

    tau1 = np.zeros((n, n), dtype=float)
    tau2 = np.zeros((n, n), dtype=float)

    for i0 in range(0, n, block_size):
        i1 = min(i0 + block_size, n)
        v0_i = vec_w0[i0:i1, :]
        v1_i = vec_w1[i0:i1, :]
        v2_i = vec_w2[i0:i1, :]
        s0_i = s0[i0:i1]
        s1_i = s1[i0:i1]
        s2_i = s2[i0:i1]

        for j0 in range(i0, n, block_size):
            j1 = min(j0 + block_size, n)
            v_j = vec[j0:j1, :]
            s0_j = s0[j0:j1]
            s1_j = s1[j0:j1]
            s2_j = s2[j0:j1]

            # sumk = sum_p (v_ip - v_jp)^2 / lam_p^(k+1), with k=0,1,2
            m0 = v0_i @ v_j.T
            m1 = v1_i @ v_j.T
            m2 = v2_i @ v_j.T

            sum0 = s0_i[:, None] + s0_j[None, :] - 2.0 * m0
            sum1 = s1_i[:, None] + s1_j[None, :] - 2.0 * m1
            sum2 = s2_i[:, None] + s2_j[None, :] - 2.0 * m2

            with np.errstate(divide="ignore", invalid="ignore"):
                t1_blk = np.divide(
                    sum1, sum0, out=np.zeros_like(sum1), where=np.abs(sum0) > tol
                )
                t2_blk = np.divide(
                    sum2, sum1, out=np.zeros_like(sum2), where=np.abs(sum1) > tol
                )

            tau1[i0:i1, j0:j1] = t1_blk
            tau2[i0:i1, j0:j1] = t2_blk
            if j0 != i0:
                tau1[j0:j1, i0:i1] = t1_blk.T
                tau2[j0:j1, i0:i1] = t2_blk.T

    np.fill_diagonal(tau1, 0.0)
    np.fill_diagonal(tau2, 0.0)
    return tau1, tau2


def compute_tau_ij_integral(a, i, j):
    # a is the connectivity matrix
    eigval, eigvec = np.linalg.eigh(-a)  # note negative sign

    lam = eigval[1:]  # ignore the p=0 mode
    vec = eigvec[:, 1:]  # remove the eigvector corresponding to p=0 mode

    v_pi = vec[i, :]  # V_{p,i}
    v_pj = vec[j, :]  # V_{p,j}
    E_ps = (v_pi - v_pj) ** 2

    sum0 = np.sum(E_ps / lam)
    sum1 = np.sum(E_ps / lam**2.0)
    sum2 = np.sum(E_ps / lam**3.0)

    tau1 = sum1 / sum0
    tau2 = sum2 / sum1
    return tau1, tau2


def compute_tau_ij_1e(K, i, j, t_min=0.0, t_max=None, factor=10.0, tol=1e-8):
    """
    Compute the 1/e relaxation time tau_{ij} for loci i and j,
    solving G2_ij(tau) / G2_ij(0) = 1/e.

    Parameters
    ----------
    K : (n×n) ndarray
        Connectivity (spring‐constant) matrix (symmetric, positive semidef.).
    i, j : int
        Indices of the two loci (0 <= i, j < n, i != j).
    t_min : float, optional
        Lower bracket for root‐finding (default: 0.0).
    t_max : float or None, optional
        Upper bracket for root‐finding. If None, set to `factor`×slowest‐mode time.
    factor : float, optional
        Multiplier for slowest‐mode timescale to set t_max (default: 10).
    tol : float, optional
        Root‐finder convergence tolerance.

    Returns
    -------
    tau_ij : float
        Relaxation time τ_{ij} where the normalized two‐point correlator equals 1/e.
    """
    # 1) Diagonalize -K so eigenvalues are positive
    eigvals, eigvecs = np.linalg.eigh(-K)

    # 2) Discard the zero (COM) mode
    lam = eigvals[1:]  # shape (n-1,)
    vec = eigvecs[:, 1:]  # shape (n,   n-1)

    # 3) Mode-difference squared for pair (i,j)
    v_i = vec[i, :]  # (n-1,)
    v_j = vec[j, :]  # (n-1,)
    E_p = (v_i - v_j) ** 2  # (n-1,)

    # 4) Compute G2_ij(0)
    G2_0 = np.sum(E_p / lam)

    # 5) Normalized correlator function
    def G_norm(t):
        return np.sum((E_p / lam) * np.exp(-lam * t)) / G2_0

    # 6) Determine t_max if needed
    if t_max is None:
        tau_slow = 1.0 / np.min(lam)
        t_max = factor * tau_slow

    # 7) Root‐find G_norm(t) = 1/e
    tau_ij = scipy.optimize.brentq(
        lambda t: G_norm(t) - 1 / np.e, t_min, t_max, xtol=tol
    )
    return tau_ij


# ------------------------------
# Stress relaxation G(t)
def compute_stress_relaxation(
    a: np.ndarray,
    t: np.ndarray,
    zeta: float = 1.0,
    stretched_exponent: float = 1.0,
):
    """
    Compute dimensionless stress relaxation function G(t) in HIPPS-DIMES model.
    G(t) here is defined as the sum over internal modes normalized by chain length

    G(t) = (1 / N) * sum_p exp(-2 * (t/tau_p)^stretched_exponent)

    summation over modes p ignore the p=0 mode (center of mass mode)
    N is the chain length, i.e. len(a)
    tau_p is defined as tau_p = - zeta / lambda_p
    """
    if not isinstance(a, np.ndarray) or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")
    if not isinstance(zeta, (int, float)) or zeta <= 0:
        raise ValueError("zeta must be a positive number")
    if not isinstance(stretched_exponent, (int, float)) or stretched_exponent <= 0:
        raise ValueError("stretched_exponent must be a positive number")

    t = np.asarray(t, dtype=float)
    if t.ndim != 1:
        raise ValueError("t must be a 1D array")

    eigvalue = np.linalg.eigvalsh(a)
    nonzero_modes = np.abs(eigvalue) > 1e-12
    if not np.any(nonzero_modes):
        raise ValueError("a must contain at least one non-zero internal mode")

    tau_p = -zeta / eigvalue[nonzero_modes]
    G_t = np.sum(
        np.exp(-2.0 * np.power(t[:, None] / tau_p[None, :], stretched_exponent)),
        axis=1,
    ) / len(a)
    return np.column_stack((t, G_t))


# Modulus
def compute_modulus(
    a: np.ndarray, freq: np.ndarray, zeta: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the storage and loss moduli for a polymer system.
    The zero eigenvalue (corresponding to center-of-mass motion) is excluded.

    Parameters
    ----------
    a : np.ndarray
        Connectivity matrix of the polymer system
    freq : np.ndarray
        Array of frequencies at which to compute the moduli
    zeta : float, optional
        Friction coefficient, by default 1.0

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Two 2D arrays containing:
        - First array: frequencies and storage modulus
        - Second array: frequencies and loss modulus

    Raises
    ------
    ValueError
        If input matrices have incorrect shapes or types
    """
    # Input validation
    if not isinstance(a, np.ndarray) or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")
    if not isinstance(freq, np.ndarray) or freq.ndim != 1:
        raise ValueError("freq must be a 1D array")
    if not isinstance(zeta, (int, float)) or zeta <= 0:
        raise ValueError("zeta must be a positive number")

    # Compute eigenvalues and eigenvectors
    eigvalue, eigvector = np.linalg.eigh(a)

    # Exclude the last eigenvalue (zero eigenvalue) and corresponding eigenvector
    eigvalue = eigvalue[:-1]  # Shape: (n_modes,)
    eigvector = eigvector[:, :-1]  # Shape: (n_monomers, n_modes)

    # Compute normal modes and relaxation times
    eigvalue_inv = 1.0 / eigvalue
    normal_modes_square_mean = -eigvalue_inv
    tau_p = -zeta / eigvalue

    # Reshape frequency for broadcasting
    freq_reshaped = np.expand_dims(freq, axis=-1)

    # Stress correlations decay as exp(-2t/tau_p), so each mode has an
    # effective Maxwell relaxation time tau_p / 2.
    omega_tau = freq_reshaped * tau_p / 2.0

    # Compute moduli
    storage_modulus = np.sum(
        omega_tau**2.0 / (1 + omega_tau**2), axis=-1
    )
    loss_modulus = np.sum(omega_tau / (1 + omega_tau**2), axis=-1)

    return (
        np.column_stack((freq, storage_modulus)),
        np.column_stack((freq, loss_modulus)),
    )


def compute_monomer_mechanical_susceptibility(
    a: np.ndarray, freq: np.ndarray, zeta: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the complex mechanical susceptibility of each monomer.
    The zero eigenvalue (corresponding to center-of-mass motion) is excluded.

    chi_i'(omega) = (1 / zeta) * sum_{p>0}
        v_{pi}^2 * tau_p / (1 + (omega * tau_p)^2)

    chi_i''(omega) = (1 / zeta) * sum_{p>0}
        v_{pi}^2 * omega * tau_p^2 / (1 + (omega * tau_p)^2)

    where tau_p = -zeta / lambda_p. Unlike the system-level modulus,
    monomer susceptibility uses tau_p directly, without a factor of two.

    Parameters
    ----------
    a : np.ndarray
        Connectivity matrix of the polymer system
    freq : np.ndarray
        Array of angular frequencies at which to compute the susceptibility
    zeta : float, optional
        Monomer friction coefficient, by default 1.0

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Three arrays containing:
        - First array: frequencies
        - Second array: real susceptibility for each monomer
          (shape: n_freqs × n_monomers)
        - Third array: imaginary susceptibility for each monomer
          (shape: n_freqs × n_monomers)

    Raises
    ------
    ValueError
        If input matrices have incorrect shapes or types
    """
    # Input validation
    if not isinstance(a, np.ndarray) or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a must be a square matrix")
    if not isinstance(freq, np.ndarray) or freq.ndim != 1:
        raise ValueError("freq must be a 1D array")
    if not isinstance(zeta, (int, float)) or zeta <= 0:
        raise ValueError("zeta must be a positive number")

    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(a)

    # Exclude the last eigenvalue (zero eigenvalue) and corresponding eigenvector
    eigvals = eigvals[:-1]  # Shape: (n_modes,)
    eigvecs = eigvecs[:, :-1]  # Shape: (n_monomers, n_modes)

    # Relaxation times τ_p
    tau_p = -zeta / eigvals  # Shape: (n_modes,)

    # Compute ωτ_p for all frequencies and modes
    omega_tau_p = (
        freq[:, np.newaxis] * tau_p[np.newaxis, :]
    )  # Shape: (n_freqs, n_modes)

    # Modal contributions to chi_i'(omega) and chi_i''(omega)
    denominator = 1 + omega_tau_p**2
    chi_prime_p = tau_p[np.newaxis, :] / (zeta * denominator)
    chi_double_prime_p = omega_tau_p * chi_prime_p

    # Square of eigenvector components (V_{pi}^2)
    eigvecs_squared = eigvecs**2  # Shape: (n_monomers, n_modes)

    # Weight each mode by v_{pi}^2 for every monomer i
    chi_prime_i = np.einsum(
        "mp,fp->fm", eigvecs_squared, chi_prime_p
    )  # Shape: (n_freqs, n_monomers)
    chi_double_prime_i = np.einsum(
        "mp,fp->fm", eigvecs_squared, chi_double_prime_p
    )  # Shape: (n_freqs, n_monomers)

    return freq, chi_prime_i, chi_double_prime_i


# ------------------------------


def Ornstein_Uhlenbeck_update(
    x,
    dt,
    k,
    zeta,
    beta,
    force_projection=0.0,
    method="euler-maruyama",
    update_zero_modes=True,
    zero_mode_tol=1e-10,
):
    """
    Update variable x for a Ornstein Uhlenbeck process

    Parameters
    ----------
    x : Array for value of x of each degree of freedom
    dt : time step
    k : Array for spring constant for each degree of freedom (eigenvalues)
    zeta : friction coefficient
    beta : inverse temperature
    force_projection : V^T * B (force projected onto eigenmodes)
    method : 'euler-maruyama' or 'exact'
    update_zero_modes : if True, zero eigenvalue modes (COM) will diffuse; if False, COM is fixed
    zero_mode_tol : tolerance for identifying zero eigenvalue modes

    Notes
    -----
    For OU process: dX = θ(μ - X)dt + σ dW
    where θ = k/ζ, σ = sqrt(2/(ζβ))

    With external force, the long-term mean μ satisfies:
    θμ = force_projection/ζ
    => μ = force_projection/(ζθ) = force_projection/k

    For Euler-Maruyama: directly use force_projection/ζ
    For exact method: calculate μ = force_projection/k (handling division by zero for zero modes)

    For standard connectivity matrices (e.g. Rouse), eigenvalues are non-positive, so
    k = -eigvalue ≥ 0 and all modes are stable or zero; both methods are then appropriate.
    """
    if isinstance(x, np.ndarray):
        rand_noise = np.random.randn(*x.shape)
    else:
        rand_noise = np.random.randn()

    theta = k[:, np.newaxis] / zeta

    # Identify zero eigenvalue modes
    zero_modes_mask = np.abs(k) < zero_mode_tol

    if method == "euler-maruyama":
        # Euler-Maruyama: dX = -θX dt + (force_projection/ζ) dt + σ dW
        # No need to divide by eigenvalue here
        dx = (
            -theta * x * dt
            + force_projection * dt / zeta
            + np.sqrt(2.0 * dt / (zeta * beta)) * rand_noise
        )
        x_new = x + dx
    elif method == "exact":
        sigma = (2.0 / (zeta * beta)) ** 0.5
        mu = np.exp(-theta * dt)

        # Calculate force drift term: μ(1 - e^(-θt))
        # For non-zero modes: μ = force_projection/k, so drift = (force_projection/k)(1 - e^(-θt))
        # For zero modes (k ≈ 0): drift reduces to (force_projection/ζ) * t
        if isinstance(force_projection, np.ndarray):
            force_drift = np.zeros_like(force_projection)
            non_zero_mask = ~zero_modes_mask

            # Non-zero modes: (force_projection/k)(1 - e^(-θt))
            force_drift[non_zero_mask] = (
                force_projection[non_zero_mask] / k[non_zero_mask, np.newaxis]
            ) * (1 - mu[non_zero_mask])

            # Zero modes: (force_projection/ζ) * t
            force_drift[zero_modes_mask] = (
                force_projection[zero_modes_mask] / zeta
            ) * dt
        else:
            # force_projection is a scalar (0.0)
            force_drift = 0.0

        # For non-zero modes: use standard Gillespie formula
        noise_term = np.sqrt((sigma**2.0 / (2.0 * theta)) * (1.0 - mu**2.0))

        # For zero modes (theta ≈ 0): analytically reduce to sqrt(sigma^2 * t)
        # As theta -> 0: sqrt(sigma^2/(2*theta) * (1 - exp(-2*theta*t))) -> sqrt(sigma^2 * t)
        noise_term[zero_modes_mask] = sigma * np.sqrt(dt)

        # Gillespie exact solution: X(t) = X(0)e^(-θt) + μ(1-e^(-θt)) + noise
        x_new = x * mu + force_drift + noise_term * rand_noise

    # If user wants to fix COM (zero modes), set them to their initial values
    if not update_zero_modes:
        x_new[zero_modes_mask] = x[zero_modes_mask]

    return x_new


def construct_connectivity_matrix_rouse(n, k):
    """
    Function to construct a ideal chain connectivity matrix given the number of monomers and the spring constant
    """
    A = np.diag(np.full(n - 1, k), 1)
    A += A.T
    A[np.diag_indices(n)] = -2 * k
    A[0, 0] = -k
    A[n - 1, n - 1] = -k
    return A


def construct_connectivity_matrix_random(n, m, k):
    """
    Generate random connected rouse chain
    n: the length of chain
    k: the spring constant
    m: number of non-consecutive bonds
    """
    A = construct_connectivity_matrix_rouse(n, k)
    pairs = list(itertools.combinations(np.arange(n), 2))

    for pair in pairs:
        if pair[1] - pair[0] == 1:
            pairs.remove(pair)
    pairs_indices = np.random.choice(len(pairs), m, replace=False)

    for idx in pairs_indices:
        pair = pairs[idx]
        A[pair[0], pair[1]] = k
        A[pair[1], pair[0]] = k

    for i in range(A.shape[0]):
        A[i, i] = -np.sum(np.delete(A[:, i], i))

    return A


def sigma2omega(sigma_mtx):
    """
    Return Omega matrix given the sigma matrix
    """
    n = sigma_mtx.shape[0]
    sigma_mtx_square = np.power(sigma_mtx, 2.0)
    sigma_row_sum = np.sum(sigma_mtx_square, axis=1)
    sigma_sum = np.sum(sigma_mtx_square)
    return (sigma_row_sum[:, np.newaxis] + sigma_row_sum - sigma_sum / n) / (
        2 * n
    ) - sigma_mtx_square / 2.0


def dmap2a_direct(dmap):
    """
    Return connectivity matrix A given the mean distance map directly through matrix peusudo inversion
    """
    sigma_mtx = 0.5 * np.sqrt(np.pi / 2.0) * dmap
    Omega = sigma2omega(sigma_mtx)
    a_direct = nearestNSD(-scipy.linalg.pinvh(Omega), 0.0)

    return a_direct


def ddmap2a_direct(ddmap):
    sigma_mtx = np.sqrt(ddmap / 3.0)
    Omega = sigma2omega(sigma_mtx)
    a_direct = nearestNSD(-scipy.linalg.pinvh(Omega), 0.0)

    return a_direct


def a2dmap_theory(A, force_positive_definite=False, return_eigenvalues=False):
    """
    Return mean distance map given the connectivity matrix A theoretically

    Parameters
    ----------
    A : np.ndarray
        Connectivity matrix (should be symmetric, negative semi-definite)
    force_positive_definite : bool, optional
        If True, replace negative temp values with zero. Default is False.
    return_eigenvalues : bool, optional
        If True, also return the eigenvalues of A. This is useful for avoiding
        redundant eigendecomposition when eigenvalues are needed elsewhere
        (e.g., for entropy computation). Default is False.

    Returns
    -------
    dmap : np.ndarray
        Mean distance map
    eigvalue : np.ndarray, optional
        Eigenvalues of A (only returned if return_eigenvalues=True)
    """
    TOL = 10**8
    # check_finite=False avoids extra scans of the matrix and is safe here
    # (A is generated/updated numerically and we already nan_to_num() in optimization)
    eigvalue, eigvector = np.linalg.eigh(A)

    temp = -1.0 / eigvalue

    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[temp >= TOL] = 0.0
    temp[temp <= -TOL] = 0.0
    # temp[np.abs(temp) <= 10**-7] = 0.0

    # replace all positive element to be zero
    if force_positive_definite:
        temp[temp < 0.0] = 0.0

    # Avoid materializing a dense diagonal matrix:
    # V @ diag(temp) @ V.T == (V * temp) @ V.T  (broadcast scales columns of V)
    Omega = (eigvector * temp) @ eigvector.T
    Omega_diag = np.diag(Omega)
    sigma = np.sqrt(Omega_diag[:, np.newaxis] + Omega_diag - 2.0 * Omega)

    dmap = 2.0 * np.sqrt(2.0 / np.pi) * sigma

    if return_eigenvalues:
        return dmap, eigvalue
    return dmap


def a2dmap_theory_with_force_applied(A, force):
    """
    Return mean distance map given the connectivity matrix A with external force applied.

    This function computes the theoretical mean distance map when the polymer is under
    constant external force. The equilibrium positions shift due to the force, and the
    distance map reflects both this shift and thermal fluctuations.

    Theory:
    -------
    Without force: <(R_i - R_j)²> = Ω_ii + Ω_jj - 2Ω_ij, where Ω = -V * Λ^(-1) * V^T

    With force: Equilibrium positions shift by R_eq = A^(-1) * F = V * (V^T * F) / Λ
    The mean squared distance becomes:
    <(R_i - R_j)²> = (R_eq,i - R_eq,j)² + (Ω_ii + Ω_jj - 2Ω_ij)
                   = equilibrium shift² + thermal fluctuations

    Parameters
    ----------
    A : (N, N) array_like
        Connectivity matrix (Laplacian-like, should have one zero eigenvalue)
    force : dict
        Dictionary specifying the force with keys:
        - 'loci': list of int, indices where force is applied
        - 'amplitude': float, magnitude of force
        - 'direction': (3,) array_like, direction of force (will be normalized)

    Returns
    -------
    dmap : (N, N) ndarray
        Mean distance map with force applied. dmap[i,j] is the mean distance
        between beads i and j under the applied force.

    Examples
    --------
    >>> A = construct_connectivity_matrix_rouse(10, 1.0)
    >>> force = {'loci': [0], 'amplitude': 5.0, 'direction': [1, 0, 0]}
    >>> dmap = a2dmap_theory_with_force_applied(A, force)
    >>> print(dmap.shape)  # (10, 10)

    Notes
    -----
    - For zero eigenvalues (COM mode), we set 1/λ = 0
    - The function returns the total mean distance including both equilibrium
      displacement and thermal fluctuations
    - The mean distance is computed as: <d> = 2*sqrt(2/π) * σ
      where σ² is the mean squared distance
    """
    TOL = 1e8

    # Validate force input
    if not isinstance(force, dict):
        raise TypeError(
            "force must be a dictionary with keys: 'loci', 'amplitude', 'direction'"
        )
    required_keys = ["loci", "amplitude", "direction"]
    for key in required_keys:
        if key not in force:
            raise ValueError(f"force dictionary must contain key: {key}")

    # Extract force parameters
    force_loci = force["loci"]
    force_amplitude = force["amplitude"]
    force_direction = np.array(force["direction"])

    # Normalize force direction
    force_direction = force_direction / np.linalg.norm(force_direction)

    # Create force vector B in real space
    N = A.shape[0]
    B = np.zeros((N, 3))
    for locus in force_loci:
        if locus < 0 or locus >= N:
            raise ValueError(f"force locus {locus} is out of range [0, {N - 1}]")
        B[locus] = force_amplitude * force_direction

    # Eigendecomposition
    eigvalue, eigvector = np.linalg.eigh(A)

    # Compute -1/eigenvalue for thermal fluctuations (Ω matrix)
    temp = -1.0 / eigvalue

    # Handle infinities and large values (zero eigenvalue handling)
    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[np.abs(temp) >= TOL] = 0.0

    # Compute Ω matrix for thermal fluctuations: Ω = V * diag(-1/λ) * V^T
    # Avoid materializing a dense diagonal matrix:
    Omega = (eigvector * temp) @ eigvector.T
    Omega_diag = np.diag(Omega)

    # Thermal fluctuation contribution (variance in each dimension)
    # σ² (per dimension) = Ω_ii + Ω_jj - 2Ω_ij
    # This represents <(R_i - R_j)²> / 3 for thermal fluctuations
    sigma_thermal = np.sqrt(Omega_diag[:, np.newaxis] + Omega_diag - 2.0 * Omega)

    # Compute equilibrium displacement due to force
    # R_eq = A^(-1) * B = V * (V^T * B) / λ
    # For zero eigenvalues, set 1/λ = 0
    temp_force = 1.0 / eigvalue
    temp_force[np.abs(temp_force) >= TOL] = 0.0
    temp_force[np.isinf(temp_force)] = 0.0

    # Project force to eigenmode space: V^T * B
    force_projection = eigvector.T @ B  # Shape: (n_modes, 3)

    # Compute equilibrium positions: R_eq = V * (V^T * B / λ)
    X_eq = force_projection * temp_force[:, np.newaxis]  # Shape: (n_modes, 3)
    R_eq = eigvector @ X_eq  # Shape: (N, 3)

    # For the distance calculation, we need to consider that:
    # - Thermal fluctuations contribute equally in all 3 dimensions: 3*σ²_thermal
    # - Equilibrium shift is a fixed vector displacement
    #
    # The mean distance for a 3D Gaussian with mean μ and isotropic variance σ² in each dim:
    # If there's also a displacement d, the distribution becomes non-central chi
    # For simplicity, we compute the squared distance from equilibrium positions
    # and add thermal fluctuation variance

    # Compute squared equilibrium separations for all pairs.
    delta_R_eq_sq = scipy.spatial.distance.cdist(
        R_eq, R_eq, metric="sqeuclidean"
    )  # Shape: (N, N)

    # Combined variance: thermal variance (3*σ²) + equilibrium shift²
    # σ²_total (in each dimension) remains sigma_thermal²
    # But equilibrium shift adds a constant offset to the distance

    # For mean distance with both effects:
    # We use: <|r|> = 2*sqrt(2/π) * sqrt(σ²_thermal + (Δr_eq/sqrt(3))²)
    # where the equilibrium shift is distributed across 3 dimensions

    sigma_combined = np.sqrt(sigma_thermal**2 + delta_R_eq_sq / 3.0)

    # Mean distance: <d> = 2*sqrt(2/π) * σ_combined
    dmap = 2.0 * np.sqrt(2.0 / np.pi) * sigma_combined

    return dmap


def dmap2cmap(dmap, rc):
    """
    Return contact map given the mean distance map and the contact threshold
    """
    sigma_mtx = 0.5 * np.sqrt(np.pi / 2.0) * dmap
    cmap = (
        scipy.special.erf(rc / (np.sqrt(2) * sigma_mtx))
        - np.sqrt(2.0 / np.pi)
        * np.exp(-0.5 * rc**2.0 / np.power(sigma_mtx, 2.0))
        * rc
        / sigma_mtx
    )
    np.fill_diagonal(cmap, 1.0)
    return cmap


def a2cmap_theory(A, rc):
    """
    Return contact map given the connectivity matrix and contact threshold, theoretically
    """
    dmap = a2dmap_theory(A)
    cmap = dmap2cmap(dmap, rc)
    return cmap


def a2a(a, fill_negative=False):
    """Return a Laplacian-form connectivity matrix on the input array backend."""
    xp = cp.get_array_module(a) if _CUPY_AVAILABLE else np
    temp = xp.copy(a)
    if fill_negative:
        temp[temp < 0.0] = 0.0
    # Zero out diagonal first, then set diagonal to negative sum of off-diagonal elements
    xp.fill_diagonal(temp, 0.0)
    xp.fill_diagonal(temp, -xp.sum(temp, axis=1))
    return temp


def optimal_rotate(P, Q, return_rotation=False, allow_reflection=False):
    """
    Return aligned matrix referred to Q
    Can return rotation matrix if return_rotation is set True
    """
    # P and Q are two sets of vectors
    P = np.matrix(P)
    Q = np.matrix(Q)

    assert P.shape == Q.shape

    Qc = np.mean(Q, axis=0)

    P = P - np.mean(P, axis=0)
    Q = Q - np.mean(Q, axis=0)

    # calculate covariance matrix A = (P^T)Q
    A = P.T * Q

    # SVD for matrix A
    V, S, Wt = np.linalg.svd(A)

    # correct rotation matrix to ensure a right-handed system if necessary
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0

    if not allow_reflection:
        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

    # calculate the final rotation matrix U
    # U = V * Wt
    U = np.dot(V, Wt)

    if not return_rotation:
        return np.array(P * U + Qc)
    else:
        return np.array(P * U + Qc), U


def write2xyz(fout, xyzs, alignment=True, allow_reflection=True):
    natoms = xyzs.shape[1]  # number of atoms
    xyz0 = xyzs[0]

    with open(fout, "w") as f:
        for snapshot in xyzs:
            if alignment:
                xyz = optimal_rotate(snapshot, xyz0, allow_reflection=allow_reflection)
            else:
                xyz = snapshot
            f.write("{}\n\n".format(natoms))
            for idx, item in enumerate(xyz):
                f.write("{} {} {} {}\n".format("C", item[0], item[1], item[2]))


def write2xyz_traj(fout, traj, recenter=False):
    natoms = traj.shape[1]

    with open(fout, "w") as f:
        for snapshot in traj:
            if recenter:
                snapshot = snapshot - snapshot.mean(axis=0)
            f.write(f"{natoms}\n\n")
            for i, atom in enumerate(snapshot):
                f.write(f"C {atom[0]} {atom[1]} {atom[2]}\n")


def a2xyz_sample(A, ensemble=1, force_positive_definite=False):
    """
    Function to generate an ensemble of configurations given the connectivity matrix
    """
    TOL = 10**8.0
    eigvalue, eigvector = np.linalg.eigh(A)
    temp = 1.0 / eigvalue[:, np.newaxis]

    # replace close zero eigvenvalue with zero
    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[temp >= TOL] = 0.0
    temp[temp <= -TOL] = 0.0

    # replace all positive element to be zero
    if force_positive_definite:
        temp[temp > 0.0] = 0.0

    # get positions
    positions = []
    for _ in range(ensemble):
        position = eigvector @ (np.sqrt(-temp) * np.random.randn(len(eigvalue), 3))
        positions.append(position)

    return np.array(positions)


def a2xyz_sample_with_force_applied(A, force, ensemble=1):
    """
    Generate ensemble of 3D structures with external force applied.

    This function generates equilibrium structures under constant external force.
    The equilibrium position is determined by: R_eq = A^(-1) * F
    where A is the connectivity matrix and F is the force vector.

    In eigenmode space: X_eq = (V^T * F) / λ, where λ are eigenvalues.
    For zero eigenvalues (COM mode), we set 1/λ = 0, which means the COM
    drifts freely and we sample it from the thermal distribution.

    Parameters
    ----------
    A : (N, N) array_like
        Connectivity matrix (Laplacian-like, should have one zero eigenvalue)
    force : dict
        Dictionary specifying the force with keys:
        - 'loci': list of int, indices where force is applied
        - 'amplitude': float, magnitude of force
        - 'direction': (3,) array_like, direction of force (will be normalized)
    ensemble : int, optional
        Number of independent samples to generate. Default is 1.

    Returns
    -------
    positions : (ensemble, N, 3) ndarray
        Sampled bead coordinates with force applied

    Examples
    --------
    >>> A = construct_connectivity_matrix_rouse(10)
    >>> force = {'loci': [0], 'amplitude': 1.0, 'direction': [1, 0, 0]}
    >>> positions = a2xyz_sample_with_force_applied(A, force, ensemble=100)
    >>> print(positions.shape)  # (100, 10, 3)

    Notes
    -----
    - The connectivity matrix A should have exactly one zero eigenvalue (COM mode)
    - Zero eigenvalues are identified and their inverse set to 0
    - The equilibrium structure represents the balance between elastic forces
      and external forces
    - Thermal fluctuations around equilibrium are sampled using Gaussian noise
    """
    TOL = 1e8

    # Validate force input
    if not isinstance(force, dict):
        raise TypeError(
            "force must be a dictionary with keys: 'loci', 'amplitude', 'direction'"
        )
    required_keys = ["loci", "amplitude", "direction"]
    for key in required_keys:
        if key not in force:
            raise ValueError(f"force dictionary must contain key: {key}")

    # Extract force parameters
    force_loci = force["loci"]
    force_amplitude = force["amplitude"]
    force_direction = np.array(force["direction"])

    # Normalize force direction
    force_direction = force_direction / np.linalg.norm(force_direction)

    # Create force vector B in real space
    N = A.shape[0]
    B = np.zeros((N, 3))
    for locus in force_loci:
        if locus < 0 or locus >= N:
            raise ValueError(f"force locus {locus} is out of range [0, {N - 1}]")
        B[locus] = force_amplitude * force_direction

    # Eigendecomposition
    eigvalue, eigvector = np.linalg.eigh(A)

    # Compute 1/eigenvalue, handling zero eigenvalues
    temp = 1.0 / eigvalue[:, np.newaxis]

    # Replace infinities and large values with zero (zero eigenvalue handling)
    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[np.abs(temp) >= TOL] = 0.0

    # Project force into eigenmode space: V^T * B
    force_projection = eigvector.T @ B  # Shape: (n_modes, 3)

    # Compute equilibrium eigenmode values: X_eq = (V^T * B) / λ
    # For zero eigenvalues, this will be zero (COM drifts freely)
    X_eq = force_projection * temp  # Broadcasting: (n_modes, 3) * (n_modes, 1)

    # Sample thermal fluctuations around equilibrium
    positions = []
    for _ in range(ensemble):
        # Random thermal fluctuations in eigenmode space
        # The variance is proportional to sqrt(-1/λ) for each mode
        thermal_fluctuations = np.sqrt(-temp) * np.random.randn(len(eigvalue), 3)

        # Total eigenmode values = equilibrium + fluctuations
        X_total = X_eq + thermal_fluctuations

        # Transform back to real space: R = V * X
        position = eigvector @ X_total
        positions.append(position)

    return np.array(positions)


def a2xyz_sample_fixed_end(
    A, xyz_start, xyz_end, ensemble=1, force_positive_definite=False
):
    """
    Generate `ensemble` random polymer configurations from connectivity matrix A,
    *with* bead 0 fixed at xyz_start and bead n−1 fixed at xyz_end.

    Parameters
    ----------
    A : (n,n) array_like
        Connectivity matrix.
    xyz_start : length‐3 array_like
        3D coordinates for bead 0.
    xyz_end : length‐3 array_like
        3D coordinates for bead n−1.
    ensemble : int, optional
        Number of independent samples to draw.  Default is 1.
    force_positive_definite : bool, optional

    Returns
    -------
    positions : (ensemble, n, 3) ndarray
        Sampled bead coordinates.
    """
    # 1) Make a copy of A and add a large negative diag‐entry at the two ends
    A_copy = np.array(A, dtype=float)
    w = -1e5
    A_copy[0, 0] += w
    A_copy[-1, -1] += w

    # 2) Eigendecompose
    evals, evecs = np.linalg.eigh(A_copy)

    # 3) Build the "temp" = 1/λ matrix, clean up infinities and huge values
    temp = 1.0 / evals[:, None]
    temp[np.isinf(temp)] = 0.0
    TOL = 1e8
    temp[np.abs(temp) >= TOL] = 0.0
    if force_positive_definite:
        temp[temp > 0.0] = 0.0

    # 4) Compute end‐to‐end distance and linear "b" term
    xyz_start = np.array(xyz_start, float)
    xyz_end = np.array(xyz_end, float)
    L = np.linalg.norm(xyz_end - xyz_start)

    n = len(evals)
    b = np.zeros((n, 3))
    b[-1, 2] = -w * L  # ensures the last bead sits at z=L in the internal frame

    # 5) Sample `ensemble` positions
    out = []
    for _ in range(ensemble):
        # random + shift in eigenspace
        coeff = np.sqrt(-temp) * np.random.randn(n, 3)
        shift = (evecs.T @ b) * (-1.0 / evals)[:, None]
        xyz = evecs @ (coeff + shift)

        # 6) Hard‐set the two ends along the z‐axis. This is not strictly necessary if w is large enough.
        xyz[0] = [0.0, 0.0, 0.0]
        xyz[-1] = [0.0, 0.0, L]

        # 7) xyz now is oriented in such way that the first monomer is at [0,0,0] and last is at [0,0,L]
        #    Hence we need to rotate+translate into the user‐specified endpoints
        ref1 = np.array([[0, 0, 0], [0, 0, L]])
        ref2 = np.vstack((xyz_start, xyz_end))
        _, R = optimal_rotate(ref1, ref2, return_rotation=True)

        c1 = ref1.mean(axis=0)
        c2 = ref2.mean(axis=0)
        xyz = (xyz - c1) @ R + c2

        out.append(xyz)
    out = np.array(out)

    return out


def interpolate_missing(matrix):
    matrix_copy = np.copy(matrix)
    x = np.arange(0, matrix_copy.shape[1])
    y = np.arange(0, matrix_copy.shape[0])
    # mask invalid values
    matrix_copy = np.ma.masked_invalid(matrix_copy)
    xx, yy = np.meshgrid(x, y)
    # get only the valid values
    x1 = xx[~matrix_copy.mask]
    y1 = yy[~matrix_copy.mask]
    newarr = matrix_copy[~matrix_copy.mask]

    GD1 = scipy.interpolate.griddata(
        (x1, y1), newarr.ravel(), (xx, yy), method="nearest"
    )
    return GD1


def objective_func(rc, A_mtx, cmap_exp):
    x = a2cmap_theory(A_mtx, rc)
    y = cmap_exp / np.nanmax(cmap_exp)
    logx = interpolate_missing(np.log(x))
    logy = interpolate_missing(np.log(y))
    res = (
        np.power(
            logx[np.triu_indices_from(logx, k=1)]
            - logy[np.triu_indices_from(logy, k=1)],
            2.0,
        ).mean()
        ** 0.5
    )
    return res


# FUNCTION TO CONVERT CMAP TO DMAP
def cmap2dmap_core(cmap_exp, rc, alpha, not_normalize, norm_max=1.0, mode="log"):
    # rc is the prefactor
    # norm_max is the maximum contact probability
    if mode == "raw":
        if not_normalize:
            log10_pmap = np.log10(cmap_exp)
        else:
            log10_pmap = (
                np.log10(cmap_exp) + np.log10(norm_max) - np.log10(np.nanmax(cmap_exp))
            )
    elif mode == "log":
        if not_normalize:
            log10_pmap = np.copy(cmap_exp)
        else:
            log10_pmap = cmap_exp + np.log10(norm_max) - np.nanmax(cmap_exp)

    return rc * 10 ** (-1.0 / alpha * log10_pmap)


def cmap2dmap(cmap, alpha, not_normalize):
    # cmap is the raw data
    # we take log on contact map
    # and then interpolate the missing data. Any zero contact pair will be interpolated
    cmap_log = interpolate_missing(np.log10(cmap))
    cmap_log = np.array((cmap_log + cmap_log.T) / 2.0)
    # lastly, convert to distance map using value of alpha
    dmap = cmap2dmap_core(cmap_log, 1.0, alpha, not_normalize)
    return dmap


def cmap2dmap_missing_data(cmap, alpha, not_normalize):
    # cmap is the raw data
    # we take log on contact map
    # unlike cmap2dmap(), this function does not interpolate the missing data. Just leave the missing data as is
    cmap_log = np.log10(cmap)
    cmap_log = np.array((cmap_log + cmap_log.T) / 2.0)
    # convert to distance map using value of alpha
    dmap = cmap2dmap_core(cmap_log, 1.0, alpha, not_normalize)
    return dmap


def nearestNSD(X, delta):
    v, w = np.linalg.eigh(X)
    v_new = np.minimum(v, delta)
    return w @ np.diag(v_new) @ w.T


def _squared_distances_from_gram(B, array_module=np):
    """Return the squared Euclidean distance matrix induced by Gram matrix ``B``."""
    diagonal = array_module.diag(B)
    ddmap = diagonal[:, array_module.newaxis] + diagonal - 2.0 * B
    ddmap = 0.5 * (ddmap + ddmap.T)
    array_module.fill_diagonal(ddmap, 0.0)
    return ddmap


def _project_centered_psd(X, gram_eigenvalue_floor=0.0, array_module=np):
    """Project onto ``B @ 1 = 0`` and ``B >= floor * J``."""
    X = 0.5 * (X + X.T)
    row_mean = array_module.mean(X, axis=1, keepdims=True)
    centered = X - row_mean - row_mean.T + array_module.mean(X)

    # The floored feasible set is a translation of the centered PSD cone:
    # {B : B >= floor * J, B @ 1 = 0} = floor * J + {C : C >= 0, C @ 1 = 0}.
    # Express J without matrix multiplication so the translational mode remains
    # exactly excluded from the positive spectral floor.
    n = centered.shape[0]
    shifted = centered.copy()
    if gram_eigenvalue_floor != 0.0:
        shifted += gram_eigenvalue_floor / n
        shifted[array_module.diag_indices(n)] -= gram_eigenvalue_floor

    if array_module is np:
        eigenvalues, eigenvectors = np.linalg.eigh(shifted)
    else:
        with cupyx.errstate(linalg="raise"):
            eigenvalues, eigenvectors = array_module.linalg.eigh(shifted)
    projected = (
        eigenvectors * array_module.maximum(eigenvalues, 0.0)
    ) @ eigenvectors.T

    # Re-center by the congruence J @ projected @ J, expressed with means to
    # avoid two dense matrix multiplications. This preserves PSD analytically.
    row_mean = array_module.mean(projected, axis=1, keepdims=True)
    projected = projected - row_mean - row_mean.T + array_module.mean(projected)
    if gram_eigenvalue_floor != 0.0:
        projected -= gram_eigenvalue_floor / n
        projected[array_module.diag_indices(n)] += gram_eigenvalue_floor
    return 0.5 * (projected + projected.T)


def _nearest_edm_objective_gradient(
    B, squared_distances, weights, array_module=np
):
    """Evaluate the unique-pair weighted EDM objective and its gradient."""
    fitted = _squared_distances_from_gram(B, array_module=array_module)
    residual = fitted - squared_distances
    weighted_residual = weights * residual

    # weights and residual are symmetric. Multiplication by 1/4 therefore
    # evaluates 1/2 * sum_{i < j} w_ij * residual_ij**2.
    objective = 0.25 * array_module.sum(weighted_residual * residual)
    gradient = (
        array_module.diag(array_module.sum(weighted_residual, axis=1))
        - weighted_residual
    )
    gradient = 0.5 * (gradient + gradient.T)
    objective = (
        float(objective.item()) if hasattr(objective, "item") else float(objective)
    )
    return objective, gradient, fitted


def nearest_edm(
    squared_distances,
    weights=None,
    *,
    gram_eigenvalue_floor=0.0,
    max_iterations=1000,
    relative_tolerance=1e-8,
    absolute_tolerance=1e-10,
    use_gpu=False,
    progress_callback=None,
):
    """Find a weighted least-squares Euclidean distance matrix.

    This function solves

    ``min_B 0.5 * sum_{i < j} w_ij * (D(B)_ij - D_obs_ij)**2``

    subject to ``B @ 1 = 0`` and ``B >= gram_eigenvalue_floor * J``, where
    ``J = I - 11.T / N`` and
    ``D(B)_ij = B_ii + B_jj - 2 * B_ij``. The input and returned distance
    matrices contain *squared* distances. A positive floor gives exactly
    ``N - 1`` positive internal Gram modes while preserving the translational
    zero mode. The default zero floor imposes no embedding-rank constraint.

    Parameters
    ----------
    squared_distances : (N, N) array_like
        Symmetric observed squared distances. NaNs may mark unobserved pairs.
        Infinite values are rejected.
    weights : (N, N) array_like, optional
        Symmetric, finite, nonnegative pair weights. A zero weight excludes a
        pair. When omitted, every finite off-diagonal pair has unit weight.
    gram_eigenvalue_floor : float, optional
        Nonnegative lower bound for every Gram eigenvalue in the centered
        subspace. It has the same units as the squared distances. The
        translational eigenvalue remains zero. Default is 0.
    max_iterations : int, optional
        Maximum number of projected-gradient iterations.
    relative_tolerance, absolute_tolerance : float, optional
        Stop when the projected-gradient norm is no larger than
        ``absolute_tolerance + relative_tolerance * ||gradient||_F``.
    use_gpu : bool, optional
        If True, run the float64 projected-gradient solver and dense PSD
        projections with CuPy. An accessible CUDA GPU is required; there is no
        silent CPU fallback. Returned matrices remain NumPy arrays.
    progress_callback : callable, optional
        Called after each accepted iteration with scalar diagnostics.

    Returns
    -------
    fitted_squared_distances : (N, N) ndarray
        The fitted symmetric, hollow Euclidean distance matrix.
    gram_matrix : (N, N) ndarray
        The fitted centered positive-semidefinite Gram matrix.
    info : dict
        Convergence status and scalar iteration history.

    Notes
    -----
    The monotone accelerated projected-gradient iteration uses a dense
    eigendecomposition for each PSD projection, so the runtime is O(N^3) per
    iteration and memory use is O(N^2) on either backend.
    """
    if use_gpu and not is_gpu_available():
        raise RuntimeError(
            "nearest EDM GPU optimization was requested, but CuPy and an "
            "accessible CUDA GPU are not available"
        )

    observed = np.asarray(squared_distances, dtype=np.float64)
    if observed.ndim != 2 or observed.shape[0] != observed.shape[1]:
        raise ValueError("squared_distances must be a square matrix")
    if observed.shape[0] < 2:
        raise ValueError("squared_distances must contain at least two loci")
    if np.any(np.isinf(observed)):
        raise ValueError("squared_distances must not contain infinite values")

    n = observed.shape[0]
    off_diagonal = ~np.eye(n, dtype=bool)

    if weights is None:
        finite = np.isfinite(observed)
        if not np.array_equal(finite, finite.T):
            raise ValueError("the observed-pair pattern must be symmetric")
        pair_mask = finite & off_diagonal
        pair_weights = np.zeros_like(observed)
        pair_weights[pair_mask] = 1.0
    else:
        pair_weights = np.asarray(weights, dtype=np.float64)
        if pair_weights.shape != observed.shape:
            raise ValueError("weights must have the same shape as squared_distances")
        if not np.all(np.isfinite(pair_weights)):
            raise ValueError("weights must be finite")
        if np.any(pair_weights < 0.0):
            raise ValueError("weights must be nonnegative")
        if not np.allclose(pair_weights, pair_weights.T, rtol=1e-12, atol=1e-14):
            raise ValueError("weights must be symmetric")
        pair_weights = 0.5 * (pair_weights + pair_weights.T)
        np.fill_diagonal(pair_weights, 0.0)
        pair_mask = (pair_weights > 0.0) & off_diagonal
        if np.any(pair_mask & ~np.isfinite(observed)):
            raise ValueError("every positive-weight pair must have a finite distance")

    if not np.any(np.triu(pair_mask, k=1)):
        raise ValueError("at least one off-diagonal pair must be observed")
    if not np.array_equal(pair_mask, pair_mask.T):
        raise ValueError("the observed-pair pattern must be symmetric")
    if not np.allclose(
        observed[pair_mask], observed.T[pair_mask], rtol=1e-10, atol=1e-12
    ):
        raise ValueError("observed squared distances must be symmetric")

    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if isinstance(gram_eigenvalue_floor, (bool, np.bool_)) or not np.isscalar(
        gram_eigenvalue_floor
    ):
        raise ValueError("gram_eigenvalue_floor must be a finite nonnegative scalar")
    try:
        gram_eigenvalue_floor = float(gram_eigenvalue_floor)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "gram_eigenvalue_floor must be a finite nonnegative scalar"
        ) from error
    if not np.isfinite(gram_eigenvalue_floor) or gram_eigenvalue_floor < 0.0:
        raise ValueError("gram_eigenvalue_floor must be a finite nonnegative scalar")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")

    target = np.zeros_like(observed)
    target[pair_mask] = observed[pair_mask]
    target = 0.5 * (target + target.T)

    weight_sum = float(0.5 * np.sum(pair_weights))
    array_module = cp if use_gpu else np
    gpu_memory_pool = cp.get_default_memory_pool() if use_gpu else None
    gpu_memory_pool_baseline_used_bytes = (
        int(gpu_memory_pool.used_bytes()) if use_gpu else None
    )
    gpu_memory_pool_baseline_total_bytes = (
        int(gpu_memory_pool.total_bytes()) if use_gpu else None
    )
    gpu_memory_pool_maximum_used_bytes = gpu_memory_pool_baseline_used_bytes
    gpu_memory_pool_maximum_total_bytes = gpu_memory_pool_baseline_total_bytes

    if use_gpu:
        cp.cuda.get_current_stream().synchronize()
    start_time = time.perf_counter()
    solver_target = array_module.asarray(target, dtype=array_module.float64)
    solver_weights = array_module.asarray(pair_weights, dtype=array_module.float64)

    def record_gpu_memory():
        nonlocal gpu_memory_pool_maximum_total_bytes
        nonlocal gpu_memory_pool_maximum_used_bytes
        if use_gpu:
            gpu_memory_pool_maximum_used_bytes = max(
                gpu_memory_pool_maximum_used_bytes,
                int(gpu_memory_pool.used_bytes()),
            )
            gpu_memory_pool_maximum_total_bytes = max(
                gpu_memory_pool_maximum_total_bytes,
                int(gpu_memory_pool.total_bytes()),
            )

    def solver_sum(value):
        result = array_module.sum(value)
        return float(result.item()) if hasattr(result, "item") else float(result)

    def solver_norm(value):
        result = array_module.linalg.norm(value, ord="fro")
        return float(result.item()) if hasattr(result, "item") else float(result)

    projection_count = 1
    line_search_projection_count = 0
    certificate_projection_count = 0
    backtracking_rejection_count = 0
    monotone_restart_count = 0
    B = _project_centered_psd(
        array_module.zeros_like(solver_target),
        gram_eigenvalue_floor,
        array_module=array_module,
    )
    record_gpu_memory()
    extrapolated = B.copy()
    momentum = 1.0
    current_objective, _, _ = _nearest_edm_objective_gradient(
        B,
        solver_target,
        solver_weights,
        array_module=array_module,
    )
    record_gpu_memory()

    # For H -> D(H), 4 * max_i sum_j w_ij is a safe Hessian-norm bound.
    # Backtracking protects the iteration against roundoff and permits the
    # bound to shrink toward a less conservative local value.
    lipschitz = 4.0 * np.max(np.sum(pair_weights, axis=1))
    lipschitz = max(float(lipschitz), np.finfo(np.float64).tiny)

    history = {
        "iteration": [],
        "objective": [],
        "weighted_rmse": [],
        "projected_gradient_norm": [],
        "relative_projected_gradient_norm": [],
        "relative_step": [],
        "step_size": [],
        "restarted": [],
    }
    converged = False
    status = "max_iterations"

    for iteration in range(1, max_iterations + 1):
        restarted = False
        while True:
            objective, gradient, _ = _nearest_edm_objective_gradient(
                extrapolated,
                solver_target,
                solver_weights,
                array_module=array_module,
            )

            for _ in range(60):
                step_size = 1.0 / lipschitz
                candidate = _project_centered_psd(
                    extrapolated - step_size * gradient,
                    gram_eigenvalue_floor,
                    array_module=array_module,
                )
                projection_count += 1
                line_search_projection_count += 1
                candidate_objective, candidate_gradient, _ = (
                    _nearest_edm_objective_gradient(
                        candidate,
                        solver_target,
                        solver_weights,
                        array_module=array_module,
                    )
                )
                record_gpu_memory()
                proximal_delta = candidate - extrapolated
                proximal_delta_norm_squared = solver_sum(
                    proximal_delta * proximal_delta
                )
                majorizer = (
                    objective
                    + solver_sum(gradient * proximal_delta)
                    + 0.5 * lipschitz * proximal_delta_norm_squared
                )
                slack = 1e-12 * max(1.0, abs(objective), abs(candidate_objective))
                if candidate_objective <= majorizer + slack:
                    break
                lipschitz *= 2.0
                backtracking_rejection_count += 1
            else:
                raise RuntimeError("nearest_edm backtracking line search failed")

            monotone_slack = 1e-12 * max(
                1.0, abs(current_objective), abs(candidate_objective)
            )
            if candidate_objective <= current_objective + monotone_slack:
                break
            if restarted:
                raise RuntimeError("nearest_edm monotone restart failed")

            extrapolated = B
            momentum = 1.0
            restarted = True
            monotone_restart_count += 1

        certificate = _project_centered_psd(
            candidate - step_size * candidate_gradient,
            gram_eigenvalue_floor,
            array_module=array_module,
        )
        projection_count += 1
        certificate_projection_count += 1
        record_gpu_memory()
        certificate_delta = candidate - certificate
        projected_gradient_norm = lipschitz * solver_norm(certificate_delta)
        gradient_norm = solver_norm(candidate_gradient)
        relative_projected_gradient_norm = projected_gradient_norm / max(
            gradient_norm, np.finfo(np.float64).tiny
        )
        accepted_delta = candidate - B
        relative_step = solver_norm(accepted_delta) / max(solver_norm(B), 1.0)
        weighted_rmse = np.sqrt(2.0 * candidate_objective / weight_sum)

        history["iteration"].append(iteration)
        history["objective"].append(candidate_objective)
        history["weighted_rmse"].append(weighted_rmse)
        history["projected_gradient_norm"].append(projected_gradient_norm)
        history["relative_projected_gradient_norm"].append(
            relative_projected_gradient_norm
        )
        history["relative_step"].append(relative_step)
        history["step_size"].append(step_size)
        history["restarted"].append(restarted)

        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "nearest_edm_initialization",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": candidate_objective,
                    "weighted_rmse": weighted_rmse,
                    "relative_projected_gradient_norm": (
                        relative_projected_gradient_norm
                    ),
                    "relative_step": relative_step,
                    "step_size": step_size,
                    "restarted": restarted,
                    "method": "nearest_edm",
                    "general_method": "covariance_initialization",
                    "use_gpu": bool(use_gpu),
                }
            )

        previous_B = B
        B = candidate
        current_objective = candidate_objective
        if projected_gradient_norm <= (
            absolute_tolerance + relative_tolerance * gradient_norm
        ):
            converged = True
            status = "optimality_tolerance"
            break

        next_momentum = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * momentum**2))
        extrapolated = B + ((momentum - 1.0) / next_momentum) * (B - previous_B)
        momentum = next_momentum
        lipschitz *= 0.9

    for key, values in history.items():
        history[key] = np.asarray(values)

    fitted_squared_distances = _squared_distances_from_gram(
        B, array_module=array_module
    )
    if use_gpu:
        fitted_squared_distances = cp.asnumpy(fitted_squared_distances)
        B = cp.asnumpy(B)
        cp.cuda.get_current_stream().synchronize()
    wall_seconds = time.perf_counter() - start_time
    final_objective = float(history["objective"][-1])
    final_weighted_rmse = float(history["weighted_rmse"][-1])
    final_projected_gradient_norm = float(history["projected_gradient_norm"][-1])
    final_relative_projected_gradient_norm = float(
        history["relative_projected_gradient_norm"][-1]
    )
    info = {
        "converged": converged,
        "status": status,
        "message": (
            "projected-gradient optimality tolerance reached"
            if converged
            else "maximum number of iterations reached"
        ),
        "iterations": int(history["iteration"][-1]),
        "objective": final_objective,
        "weighted_rmse": final_weighted_rmse,
        "projected_gradient_norm": final_projected_gradient_norm,
        "relative_projected_gradient_norm": final_relative_projected_gradient_norm,
        "gram_eigenvalue_floor": gram_eigenvalue_floor,
        "weight_sum": float(weight_sum),
        "backend": "gpu" if use_gpu else "cpu",
        "dtype": "float64",
        "gpu_device": get_gpu_name() if use_gpu else None,
        "cupy_version": cp.__version__ if use_gpu else None,
        "wall_seconds": float(wall_seconds),
        "projection_count": projection_count,
        "line_search_projection_count": line_search_projection_count,
        "certificate_projection_count": certificate_projection_count,
        "backtracking_rejection_count": backtracking_rejection_count,
        "monotone_restart_count": monotone_restart_count,
        "gpu_memory_pool_baseline_used_bytes": (
            gpu_memory_pool_baseline_used_bytes
        ),
        "gpu_memory_pool_maximum_used_bytes": (
            gpu_memory_pool_maximum_used_bytes
        ),
        "gpu_memory_pool_baseline_total_bytes": (
            gpu_memory_pool_baseline_total_bytes
        ),
        "gpu_memory_pool_maximum_total_bytes": (
            gpu_memory_pool_maximum_total_bytes
        ),
        "history": history,
    }
    if not converged:
        warnings.warn(
            "nearest_edm reached max_iterations before satisfying the "
            "projected-gradient optimality tolerance",
            RuntimeWarning,
            stacklevel=2,
        )

    return fitted_squared_distances, B, info


_COVARIANCE_PRECONDITIONER_PAIR_BLOCK_SIZE = 4096


def _centered_orthonormal_basis(n):
    """Return a deterministic orthonormal basis for ``1.T @ x = 0``."""
    basis = np.zeros((n, n - 1), dtype=np.float64)
    for column in range(n - 1):
        denominator = np.sqrt((column + 1) * (column + 2))
        basis[: column + 1, column] = 1.0 / denominator
        basis[column + 1, column] = -(column + 1) / denominator
    return basis


def _positive_finite_scalar(value, name):
    """Validate and return one strictly positive finite scalar."""
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a positive finite scalar")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a positive finite scalar") from error
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return value


def _validate_gaussian_covariance_inputs(
    squared_distances, noise_variance, relative_noise_std
):
    """Validate observed pairs and construct their Gaussian variances."""
    observed = np.asarray(squared_distances, dtype=np.float64)
    if observed.ndim != 2 or observed.shape[0] != observed.shape[1]:
        raise ValueError("squared_distances must be a square matrix")
    if observed.shape[0] < 2:
        raise ValueError("squared_distances must contain at least two loci")
    if np.any(np.isinf(observed)):
        raise ValueError("squared_distances must not contain infinite values")

    n = observed.shape[0]
    off_diagonal = ~np.eye(n, dtype=bool)
    finite = np.isfinite(observed)
    if not np.array_equal(finite & off_diagonal, finite.T & off_diagonal):
        raise ValueError("the observed-pair pattern must be symmetric")
    pair_mask = finite & off_diagonal
    upper_mask = np.triu(pair_mask, k=1)
    if not np.any(upper_mask):
        raise ValueError("at least one off-diagonal pair must be observed")
    if not np.allclose(
        observed[pair_mask], observed.T[pair_mask], rtol=1e-10, atol=1e-12
    ):
        raise ValueError("observed squared distances must be symmetric")
    if np.any(observed[pair_mask] <= 0.0):
        raise ValueError("observed off-diagonal squared distances must be positive")

    has_absolute = noise_variance is not None
    has_relative = relative_noise_std is not None
    if has_absolute == has_relative:
        raise ValueError(
            "provide exactly one of noise_variance or relative_noise_std"
        )

    pair_i, pair_j = np.where(upper_mask)
    target_pairs = observed[pair_i, pair_j]
    if has_absolute:
        noise_parameter = _positive_finite_scalar(
            noise_variance, "noise_variance"
        )
        pair_variance = np.full(
            target_pairs.shape, noise_parameter, dtype=np.float64
        )
        noise_model = "homoskedastic_absolute_variance"
    else:
        noise_parameter = _positive_finite_scalar(
            relative_noise_std, "relative_noise_std"
        )
        pair_variance = np.square(noise_parameter * target_pairs)
        noise_model = "heteroskedastic_relative_std"

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        inverse_variance = 1.0 / pair_variance
    if not np.all(np.isfinite(pair_variance)) or not np.all(
        np.isfinite(inverse_variance)
    ):
        raise ValueError(
            "the noise specification must produce positive finite pair variances "
            "and inverse variances"
        )

    return (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        target_pairs,
        pair_variance,
        noise_model,
        noise_parameter,
    )


def _gaussian_covariance_pair_matrices(
    n, pair_i, pair_j, target_pairs, inverse_variance
):
    """Expand validated unique-pair Gaussian data into symmetric matrices."""
    target_matrix = np.zeros((n, n), dtype=np.float64)
    target_matrix[pair_i, pair_j] = target_pairs
    target_matrix[pair_j, pair_i] = target_pairs

    inverse_variance_matrix = np.zeros((n, n), dtype=np.float64)
    inverse_variance_matrix[pair_i, pair_j] = inverse_variance
    inverse_variance_matrix[pair_j, pair_i] = inverse_variance
    return target_matrix, inverse_variance_matrix


def _gaussian_covariance_data_objective_gradient(
    reduced_gram,
    basis,
    target_matrix,
    inverse_variance_matrix,
    pair_i,
    pair_j,
    *,
    array_module=np,
):
    """Evaluate the Gaussian data term through the structured EDM operator."""
    reduced_gram = array_module.asarray(reduced_gram, dtype=array_module.float64)
    gram = basis @ reduced_gram @ basis.T
    fitted_matrix = _squared_distances_from_gram(gram, array_module)
    residual_matrix = fitted_matrix - target_matrix
    weighted_residual = inverse_variance_matrix * residual_matrix
    data_objective = 0.25 * array_module.sum(
        weighted_residual * residual_matrix
    )
    full_gradient = (
        array_module.diag(array_module.sum(weighted_residual, axis=1))
        - weighted_residual
    )
    data_gradient = basis.T @ full_gradient @ basis
    data_gradient = 0.5 * (data_gradient + data_gradient.T)
    return data_objective, data_gradient, fitted_matrix[pair_i, pair_j]


def _gaussian_covariance_data_hessian_action(
    direction, basis, inverse_variance_matrix, *, array_module=np
):
    """Apply the Gaussian data Hessian without materializing pair tensors."""
    direction = array_module.asarray(direction, dtype=array_module.float64)
    direction = 0.5 * (direction + direction.T)
    full_direction = basis @ direction @ basis.T
    distance_direction = _squared_distances_from_gram(
        full_direction, array_module
    )
    weighted_direction = inverse_variance_matrix * distance_direction
    full_action = (
        array_module.diag(array_module.sum(weighted_direction, axis=1))
        - weighted_direction
    )
    reduced_action = basis.T @ full_action @ basis
    return 0.5 * (reduced_action + reduced_action.T)


def _gaussian_covariance_objective_gradient(
    reduced_gram,
    basis,
    target_matrix,
    inverse_variance_matrix,
    pair_i,
    pair_j,
    *,
    array_module=np,
):
    """Evaluate the Gaussian soft-constraint objective in the covariance cone."""
    reduced_gram = array_module.asarray(reduced_gram, dtype=array_module.float64)
    if _CUPY_AVAILABLE and array_module is cp:
        try:
            with cupyx.errstate(linalg="raise"):
                cholesky_factor = cp.linalg.cholesky(reduced_gram)
        except np.linalg.LinAlgError:
            return np.inf, None, None, None, None
        identity = cp.eye(reduced_gram.shape[0], dtype=cp.float64)
        inverse_gram = cupyx_scipy_linalg.solve_triangular(
            cholesky_factor,
            identity,
            lower=True,
            check_finite=False,
        )
        inverse_gram = cupyx_scipy_linalg.solve_triangular(
            cholesky_factor.T,
            inverse_gram,
            lower=False,
            check_finite=False,
        )
    else:
        try:
            cholesky_factor = np.linalg.cholesky(reduced_gram)
        except np.linalg.LinAlgError:
            return np.inf, None, None, None, None
        inverse_gram = scipy.linalg.cho_solve(
            (cholesky_factor, True),
            np.eye(reduced_gram.shape[0], dtype=np.float64),
            check_finite=False,
        )

    logdet = 2.0 * array_module.sum(
        array_module.log(array_module.diag(cholesky_factor))
    )
    inverse_gram = 0.5 * (inverse_gram + inverse_gram.T)
    data_objective, data_gradient, fitted_pairs = (
        _gaussian_covariance_data_objective_gradient(
            reduced_gram,
            basis,
            target_matrix,
            inverse_variance_matrix,
            pair_i,
            pair_j,
            array_module=array_module,
        )
    )
    negative_entropy = -1.5 * logdet
    entropy_gradient = -1.5 * inverse_gram
    gradient = data_gradient + entropy_gradient
    gradient = 0.5 * (gradient + gradient.T)
    return (
        float((negative_entropy + data_objective).item()),
        gradient,
        fitted_pairs,
        inverse_gram,
        (
            float(negative_entropy.item()),
            float(data_objective.item()),
            entropy_gradient,
            data_gradient,
        ),
    )


def _gaussian_covariance_data_preconditioner_diagonal(
    basis,
    pair_i,
    pair_j,
    inverse_variance,
    *,
    block_size=_COVARIANCE_PRECONDITIONER_PAIR_BLOCK_SIZE,
    progress_callback=None,
    array_module=np,
    use_gpu=False,
):
    """Build the exact data-Hessian diagonal with bounded pair-block memory."""
    n_modes = basis.shape[1]
    diagonal = array_module.zeros(
        (n_modes, n_modes), dtype=array_module.float64
    )
    pair_count = int(pair_i.size)
    block_count = (pair_count + block_size - 1) // block_size
    if use_gpu:
        cp.cuda.get_current_stream().synchronize()
    start_time = time.perf_counter()
    for block_index, start in enumerate(range(0, pair_count, block_size), start=1):
        stop = min(start + block_size, pair_count)
        pair_vectors = basis[pair_i[start:stop]] - basis[pair_j[start:stop]]
        pair_vectors_squared = array_module.square(pair_vectors)
        diagonal += pair_vectors_squared.T @ (
            inverse_variance[start:stop, array_module.newaxis]
            * pair_vectors_squared
        )
        if progress_callback is not None:
            if use_gpu:
                cp.cuda.get_current_stream().synchronize()
            progress_callback(
                {
                    "stage": "covariance_preconditioner",
                    "iteration": block_index,
                    "total": block_count,
                    "pairs_completed": stop,
                    "pair_count": pair_count,
                    "method": "COV",
                    "general_method": "covariance_optimization",
                    "use_gpu": bool(use_gpu),
                    "noisy": True,
                }
            )
    if use_gpu:
        cp.cuda.get_current_stream().synchronize()
    return diagonal, time.perf_counter() - start_time


def _preconditioned_conjugate_gradient(
    operator,
    right_hand_side,
    preconditioner_diagonal,
    relative_tolerance,
    max_iterations,
):
    """Solve a symmetric positive-definite matrix equation without forming it."""
    right_hand_side = np.asarray(right_hand_side, dtype=np.float64)
    diagonal = np.asarray(preconditioner_diagonal, dtype=np.float64)
    solution = np.zeros_like(right_hand_side)
    residual = right_hand_side.copy()
    scale = float(np.linalg.norm(right_hand_side, ord="fro"))
    if scale == 0.0:
        return solution, 0, 0.0, True
    tolerance = float(relative_tolerance) * scale

    preconditioned = residual / diagonal
    direction = preconditioned.copy()
    residual_product = float(np.sum(residual * preconditioned))
    converged = False
    iteration = 0
    for iteration in range(1, max_iterations + 1):
        operator_direction = operator(direction)
        curvature = float(np.sum(direction * operator_direction))
        if not np.isfinite(curvature) or curvature <= 0.0:
            break
        step = residual_product / curvature
        solution += step * direction
        residual -= step * operator_direction
        residual_norm = float(np.linalg.norm(residual, ord="fro"))
        if residual_norm <= tolerance:
            converged = True
            break
        preconditioned = residual / diagonal
        next_product = float(np.sum(residual * preconditioned))
        if not np.isfinite(next_product) or next_product <= 0.0:
            break
        direction = preconditioned + (next_product / residual_product) * direction
        residual_product = next_product

    solution = 0.5 * (solution + solution.T)
    residual_norm = float(np.linalg.norm(residual, ord="fro"))
    return solution, iteration, residual_norm, converged


def _preconditioned_conjugate_gradient_gpu(
    operator,
    right_hand_side,
    preconditioner_diagonal,
    relative_tolerance,
    max_iterations,
):
    """Solve the COV Newton equation with CuPy's matrix-free CG."""
    matrix_shape = right_hand_side.shape
    vector_size = int(right_hand_side.size)
    right_hand_side_vector = right_hand_side.reshape(vector_size)
    diagonal_vector = preconditioner_diagonal.reshape(vector_size)
    scale = float(cp.linalg.norm(right_hand_side_vector).item())
    if scale == 0.0:
        return cp.zeros_like(right_hand_side), 0, 0.0, True

    def matrix_vector_product(vector):
        return operator(vector.reshape(matrix_shape)).reshape(vector_size)

    linear_operator = cupyx_sparse_linalg.LinearOperator(
        (vector_size, vector_size),
        matvec=matrix_vector_product,
        dtype=cp.float64,
    )
    preconditioner = cupyx_sparse_linalg.LinearOperator(
        (vector_size, vector_size),
        matvec=lambda vector: vector / diagonal_vector,
        dtype=cp.float64,
    )
    iteration_count = [0]

    def count_iteration(_):
        iteration_count[0] += 1

    solution_vector, solver_status = cupyx_sparse_linalg.cg(
        linear_operator,
        right_hand_side_vector,
        tol=float(relative_tolerance),
        atol=0.0,
        maxiter=int(max_iterations),
        M=preconditioner,
        callback=count_iteration,
    )
    solution = solution_vector.reshape(matrix_shape)
    solution = 0.5 * (solution + solution.T)
    residual = right_hand_side - operator(solution)
    residual_norm = float(cp.linalg.norm(residual).item())
    converged = int(solver_status) == 0
    return solution, iteration_count[0], residual_norm, converged


def _rouse_initial_connectivity(squared_distances, pair_i, pair_j):
    """Return a Rouse chain matched to the observed-pair distance mean."""
    squared_distances = np.asarray(squared_distances, dtype=np.float64)
    observed_pair_mean_squared_distance = _positive_finite_scalar(
        np.mean(squared_distances[pair_i, pair_j]),
        "observed-pair mean squared distance",
    )
    unit_spring_rouse_pair_mean_squared_distance = 3.0 * float(
        np.mean(pair_j - pair_i)
    )
    spring_constant = (
        unit_spring_rouse_pair_mean_squared_distance
        / observed_pair_mean_squared_distance
    )
    connectivity = construct_connectivity_matrix_rouse(
        squared_distances.shape[0], spring_constant
    )
    return (
        connectivity,
        observed_pair_mean_squared_distance,
        unit_spring_rouse_pair_mean_squared_distance,
        spring_constant,
    )


def _reduced_gram_from_connectivity(connectivity, basis):
    """Validate a centered stable connectivity and return its internal Gram."""
    n = basis.shape[0]
    connectivity = np.asarray(connectivity, dtype=np.float64)
    if connectivity.shape != (n, n) or not np.all(np.isfinite(connectivity)):
        raise ValueError("initial_connectivity must be a finite NxN matrix")
    if not np.allclose(connectivity, connectivity.T, rtol=1e-10, atol=1e-12):
        raise ValueError("initial_connectivity must be symmetric")
    row_scale = max(float(np.max(np.abs(connectivity))), 1.0)
    row_sum_error = float(np.max(np.abs(np.sum(connectivity, axis=1))))
    if row_sum_error > 1e-10 * row_scale:
        raise ValueError("initial_connectivity rows must sum to zero")
    reduced_connectivity = basis.T @ (-connectivity) @ basis
    eigenvalues = np.linalg.eigvalsh(reduced_connectivity)
    if float(eigenvalues[0]) <= 0.0:
        raise ValueError("initial_connectivity must be stable in every internal mode")
    reduced_gram = 3.0 * np.linalg.solve(
        reduced_connectivity,
        np.eye(n - 1, dtype=np.float64),
    )
    return reduced_gram, eigenvalues


def _connectivity_from_reduced_gram(reduced_gram, basis):
    """Return the physical connectivity associated with an internal Gram matrix."""
    reduced_connectivity = 3.0 * np.linalg.solve(
        reduced_gram, np.eye(reduced_gram.shape[0], dtype=np.float64)
    )
    connectivity = -(basis @ reduced_connectivity @ basis.T)
    connectivity = 0.5 * (connectivity + connectivity.T)
    return a2a(connectivity)


def _calibrate_gaussian_covariance_initial_scale(
    reduced_gram,
    basis,
    pair_i,
    pair_j,
    target_pairs,
    inverse_variance,
    *,
    array_module=np,
):
    """Exactly minimize the COV objective over one positive Gram scale."""
    use_gpu = _CUPY_AVAILABLE and array_module is cp
    if use_gpu:
        cp.cuda.get_current_stream().synchronize()
    start_time = time.perf_counter()

    gram = basis @ reduced_gram @ basis.T
    fitted_pairs = _squared_distances_from_gram(
        gram, array_module=array_module
    )[pair_i, pair_j]
    distance_square_coefficient = float(
        array_module.sum(
            inverse_variance * array_module.square(fitted_pairs)
        ).item()
    )
    distance_target_coefficient = float(
        array_module.sum(inverse_variance * fitted_pairs * target_pairs).item()
    )
    if (
        not np.isfinite(distance_square_coefficient)
        or distance_square_coefficient <= 0.0
        or not np.isfinite(distance_target_coefficient)
        or distance_target_coefficient <= 0.0
    ):
        raise RuntimeError(
            "initial scalar calibration requires positive finite COV coefficients"
        )

    n_modes = reduced_gram.shape[0]
    entropy_scale_coefficient = 1.5 * n_modes
    square_root_term = 2.0 * np.sqrt(distance_square_coefficient) * np.sqrt(
        entropy_scale_coefficient
    )
    discriminant_root = np.hypot(
        distance_target_coefficient, square_root_term
    )
    scale_factor = 0.5 * (
        distance_target_coefficient / distance_square_coefficient
        + discriminant_root / distance_square_coefficient
    )
    if not np.isfinite(scale_factor) or scale_factor <= 0.0:
        raise RuntimeError(
            "initial scalar calibration produced a nonpositive or nonfinite scale"
        )

    determinant_sign, logdet = array_module.linalg.slogdet(reduced_gram)
    determinant_sign = float(determinant_sign.item())
    logdet = float(logdet.item())
    if determinant_sign <= 0.0 or not np.isfinite(logdet):
        raise RuntimeError(
            "initial scalar calibration requires a positive-definite Gram matrix"
        )
    residual_before = fitted_pairs - target_pairs
    residual_after = scale_factor * fitted_pairs - target_pairs
    data_objective_before = 0.5 * float(
        array_module.sum(
            inverse_variance * array_module.square(residual_before)
        ).item()
    )
    data_objective_after = 0.5 * float(
        array_module.sum(
            inverse_variance * array_module.square(residual_after)
        ).item()
    )
    objective_before = -1.5 * logdet + data_objective_before
    objective_after = (
        -1.5 * (logdet + n_modes * np.log(scale_factor))
        + data_objective_after
    )
    derivative_residual = (
        distance_square_coefficient * scale_factor
        - distance_target_coefficient
        - entropy_scale_coefficient / scale_factor
    )
    derivative_scale = max(
        abs(distance_square_coefficient * scale_factor),
        abs(distance_target_coefficient),
        abs(entropy_scale_coefficient / scale_factor),
        1.0,
    )
    reduced_gram = scale_factor * reduced_gram
    if use_gpu:
        cp.cuda.get_current_stream().synchronize()

    return reduced_gram, {
        "method": "exact_covariance_objective_ray_minimum",
        "scale_factor": float(scale_factor),
        "connectivity_scale_factor": float(1.0 / scale_factor),
        "objective_before": float(objective_before),
        "objective_after": float(objective_after),
        "objective_reduction": float(objective_before - objective_after),
        "relative_derivative_residual": float(
            abs(derivative_residual) / derivative_scale
        ),
        "backend": "gpu" if use_gpu else "cpu",
        "wall_seconds": time.perf_counter() - start_time,
    }


def _initialize_gaussian_reduced_gram(
    observed,
    pair_mask,
    pair_i,
    pair_j,
    inverse_variance,
    basis,
    initialization,
    initial_connectivity,
    initial_gram_floor_relative,
    use_gpu=False,
    progress_callback=None,
):
    """Return the requested physical initialization for COV."""
    if initialization not in {"rouse", "nearest_edm"}:
        raise ValueError("initialization must be 'rouse' or 'nearest_edm'")
    if initial_connectivity is not None and initialization == "nearest_edm":
        raise ValueError(
            "initial_connectivity cannot be combined with nearest_edm initialization"
        )

    start_time = time.perf_counter()
    n = observed.shape[0]
    n_modes = n - 1
    if initial_connectivity is not None:
        reduced_gram, eigenvalues = _reduced_gram_from_connectivity(
            initial_connectivity, basis
        )
        initialization_info = {
            "kind": "provided_connectivity",
            "minimum_internal_stiffness": float(eigenvalues[0]),
            "maximum_internal_stiffness": float(eigenvalues[-1]),
        }
    elif initialization == "rouse":
        (
            connectivity,
            observed_pair_mean_squared_distance,
            unit_spring_rouse_pair_mean_squared_distance,
            spring_constant,
        ) = _rouse_initial_connectivity(observed, pair_i, pair_j)
        reduced_gram, eigenvalues = _reduced_gram_from_connectivity(
            connectivity, basis
        )
        initialization_info = {
            "kind": "rouse",
            "observed_pair_mean_squared_distance": (
                observed_pair_mean_squared_distance
            ),
            "unit_spring_rouse_pair_mean_squared_distance": (
                unit_spring_rouse_pair_mean_squared_distance
            ),
            "spring_constant": spring_constant,
            "minimum_internal_stiffness": float(eigenvalues[0]),
            "maximum_internal_stiffness": float(eigenvalues[-1]),
        }
    else:
        normalized_weights = np.zeros_like(observed)
        normalized_pair_weights = inverse_variance / float(
            np.mean(inverse_variance)
        )
        normalized_weights[pair_i, pair_j] = normalized_pair_weights
        normalized_weights[pair_j, pair_i] = normalized_pair_weights
        projection_target = observed.copy()
        projection_target[~pair_mask] = np.nan
        np.fill_diagonal(projection_target, 0.0)
        _, closest_gram, initializer = nearest_edm(
            projection_target,
            normalized_weights,
            max_iterations=2000,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
            use_gpu=use_gpu,
            progress_callback=progress_callback,
        )
        gram_scale = float(np.trace(closest_gram) / n_modes)
        if not np.isfinite(gram_scale) or gram_scale <= 0.0:
            raise RuntimeError("nearest-EDM initialization has a nonpositive scale")
        initial_floor = initial_gram_floor_relative * gram_scale
        closest_gram = _project_centered_psd(closest_gram, initial_floor)
        reduced_gram = basis.T @ closest_gram @ basis
        initialization_info = {
            "kind": "weighted_nearest_edm",
            "nearest_edm_converged": bool(initializer["converged"]),
            "nearest_edm_status": initializer["status"],
            "nearest_edm_iterations": int(initializer["iterations"]),
            "nearest_edm_weighted_rmse": float(initializer["weighted_rmse"]),
            "nearest_edm_wall_seconds": float(initializer["wall_seconds"]),
            "nearest_edm_projection_count": int(initializer["projection_count"]),
            "dtype": initializer["dtype"],
            "gpu_device": initializer["gpu_device"],
            "cupy_version": initializer["cupy_version"],
            "gpu_memory_pool_baseline_used_bytes": initializer[
                "gpu_memory_pool_baseline_used_bytes"
            ],
            "gpu_memory_pool_maximum_used_bytes": initializer[
                "gpu_memory_pool_maximum_used_bytes"
            ],
            "gpu_memory_pool_baseline_total_bytes": initializer[
                "gpu_memory_pool_baseline_total_bytes"
            ],
            "gpu_memory_pool_maximum_total_bytes": initializer[
                "gpu_memory_pool_maximum_total_bytes"
            ],
            "gram_scale": gram_scale,
            "gram_eigenvalue_floor": initial_floor,
            "backend": initializer["backend"],
        }

    reduced_gram = 0.5 * (reduced_gram + reduced_gram.T)
    try:
        np.linalg.cholesky(reduced_gram)
    except np.linalg.LinAlgError as error:
        raise RuntimeError(
            "initial covariance is not strictly positive definite"
        ) from error
    initialization_info.setdefault("backend", "cpu")
    initialization_info["wall_seconds"] = time.perf_counter() - start_time
    return reduced_gram, initialization_info


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
):
    """Fit Gaussian-noisy squared distances inside the covariance cone.

    The strictly positive internal Gram matrix ``B`` minimizes

    ``-3/2 logdet(B) + 1/2 sum_{i<j} (D(B)_ij-Dobs_ij)^2 / v_ij``.

    Supply either one absolute variance or one relative standard deviation
    ``sigma_ij / Dobs_ij`` shared by all observed pairs. The latter produces
    ``v_ij = (relative_noise_std * Dobs_ij)**2``. With ``use_gpu=True``,
    Newton-CG and the once-per-fit exact blockwise data-Hessian setup run in
    float64 on the GPU. The default Rouse initialization remains on the CPU;
    the optional nearest-EDM initialization uses the requested backend. Every
    initial Gram shape is exactly rescaled along its positive ray for the COV
    objective before Newton starts.
    """
    if use_gpu and not is_gpu_available():
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
    ) = _validate_gaussian_covariance_inputs(
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
    if not np.isfinite(initial_gram_floor_relative) or (
        initial_gram_floor_relative <= 0.0
    ):
        raise ValueError("initial_gram_floor_relative must be positive and finite")

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

    basis = _centered_orthonormal_basis(n)
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = _gaussian_covariance_pair_matrices(
        n, pair_i, pair_j, target_pairs, inverse_variance
    )
    reduced_gram, initialization_info = _initialize_gaussian_reduced_gram(
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

    if use_gpu:
        solver_basis = cp.asarray(basis, dtype=cp.float64)
        solver_pair_i = cp.asarray(pair_i, dtype=cp.int32)
        solver_pair_j = cp.asarray(pair_j, dtype=cp.int32)
        solver_target_pairs = cp.asarray(target_pairs, dtype=cp.float64)
        solver_inverse_variance = cp.asarray(inverse_variance, dtype=cp.float64)
        solver_target_matrix = cp.asarray(target_matrix, dtype=cp.float64)
        solver_inverse_variance_matrix = cp.asarray(
            inverse_variance_matrix, dtype=cp.float64
        )
        reduced_gram = cp.asarray(reduced_gram, dtype=cp.float64)
        conjugate_gradient_function = (
            _preconditioned_conjugate_gradient_gpu
        )
        array_module = cp
    else:
        solver_basis = basis
        solver_pair_i = pair_i
        solver_pair_j = pair_j
        solver_target_pairs = target_pairs
        solver_inverse_variance = inverse_variance
        solver_target_matrix = target_matrix
        solver_inverse_variance_matrix = inverse_variance_matrix
        conjugate_gradient_function = _preconditioned_conjugate_gradient
        array_module = np

    reduced_gram, scalar_calibration = (
        _calibrate_gaussian_covariance_initial_scale(
            reduced_gram,
            solver_basis,
            solver_pair_i,
            solver_pair_j,
            solver_target_pairs,
            solver_inverse_variance,
            array_module=array_module,
        )
    )
    initialization_info["scalar_calibration"] = scalar_calibration
    initialization_info["wall_seconds"] += scalar_calibration["wall_seconds"]
    if initialization_info["kind"] == "rouse":
        initialization_info["effective_spring_constant"] = (
            initialization_info["spring_constant"]
            * scalar_calibration["connectivity_scale_factor"]
        )

    data_hessian_diagonal, preconditioner_setup_seconds = (
        _gaussian_covariance_data_preconditioner_diagonal(
            solver_basis,
            solver_pair_i,
            solver_pair_j,
            solver_inverse_variance,
            progress_callback=progress_callback,
            array_module=array_module,
            use_gpu=use_gpu,
        )
    )

    def solver_scalar(value):
        return float(value.item()) if hasattr(value, "item") else float(value)

    def solver_norm(value):
        return solver_scalar(array_module.sqrt(array_module.sum(value * value)))

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
        "newton_decrement": [],
        "relative_step": [],
        "step_size": [],
        "cg_iterations": [],
        "cg_relative_residual": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}
    converged = False
    status = "max_iterations"
    message = "maximum number of iterations reached"
    current_state = _gaussian_covariance_objective_gradient(
        reduced_gram,
        solver_basis,
        solver_target_matrix,
        solver_inverse_variance_matrix,
        solver_pair_i,
        solver_pair_j,
        array_module=array_module,
    )

    for iteration in range(1, max_iterations + 1):
        objective, gradient, _, inverse_gram, components = current_state
        negative_entropy, data_objective, entropy_gradient, data_gradient = components
        gradient_norm = solver_norm(gradient)
        gradient_scale = max(
            1.0,
            solver_norm(entropy_gradient),
            solver_norm(data_gradient),
        )
        relative_gradient_norm = gradient_norm / gradient_scale
        if gradient_norm <= absolute_tolerance + relative_tolerance * gradient_scale:
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break

        def hessian_operator(direction):
            direction = 0.5 * (direction + direction.T)
            data_term = _gaussian_covariance_data_hessian_action(
                direction,
                solver_basis,
                solver_inverse_variance_matrix,
                array_module=array_module,
            )
            entropy_term = 1.5 * inverse_gram @ direction @ inverse_gram
            result = data_term + entropy_term
            return 0.5 * (result + result.T)

        inverse_diagonal = array_module.diag(inverse_gram)
        preconditioner_diagonal = data_hessian_diagonal + 1.5 * array_module.outer(
            inverse_diagonal, inverse_diagonal
        )
        preconditioner_diagonal = array_module.maximum(
            preconditioner_diagonal, np.finfo(np.float64).tiny
        )
        forcing_tolerance = min(
            cg_relative_tolerance,
            max(1e-8, 0.1 * np.sqrt(max(relative_gradient_norm, 0.0))),
        )
        direction, cg_iterations, cg_residual, cg_converged = (
            conjugate_gradient_function(
                hessian_operator,
                -gradient,
                preconditioner_diagonal,
                forcing_tolerance,
                int(cg_max_iterations),
            )
        )
        directional_derivative = solver_scalar(
            array_module.sum(gradient * direction)
        )
        if (
            not np.isfinite(directional_derivative)
            or directional_derivative >= 0.0
        ):
            direction = -gradient / preconditioner_diagonal
            direction = 0.5 * (direction + direction.T)
            directional_derivative = solver_scalar(
                array_module.sum(gradient * direction)
            )
            cg_converged = False

        step_size = 1.0
        candidate = None
        candidate_state = None
        for _ in range(60):
            trial = reduced_gram + step_size * direction
            trial = 0.5 * (trial + trial.T)
            trial_state = _gaussian_covariance_objective_gradient(
                trial,
                solver_basis,
                solver_target_matrix,
                solver_inverse_variance_matrix,
                solver_pair_i,
                solver_pair_j,
                array_module=array_module,
            )
            trial_objective = trial_state[0]
            if not np.isfinite(trial_objective):
                step_size *= 0.5
                continue
            if trial_objective <= (
                objective + 1e-4 * step_size * directional_derivative
            ):
                candidate = trial
                candidate_state = trial_state
                break
            step_size *= 0.5
        if candidate is None:
            status = "line_search_failed"
            message = "covariance-cone Armijo line search failed"
            break

        accepted_delta = candidate - reduced_gram
        relative_step = solver_norm(accepted_delta) / max(
            solver_norm(reduced_gram), 1.0
        )
        reduced_gram = candidate
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
        current_state = candidate_state
        candidate_gradient_norm = solver_norm(candidate_gradient)
        candidate_gradient_scale = max(
            1.0,
            solver_norm(candidate_entropy_gradient),
            solver_norm(candidate_data_gradient),
        )
        candidate_relative_gradient = (
            candidate_gradient_norm / candidate_gradient_scale
        )
        residual = candidate_pairs - solver_target_pairs
        residual_squared = array_module.square(residual)
        relative_loss = solver_scalar(
            array_module.sqrt(
                array_module.mean(
                    array_module.square(residual / solver_target_pairs)
                )
            )
        )
        weighted_rmse = solver_scalar(
            array_module.sqrt(
                array_module.mean(residual_squared * solver_inverse_variance)
            )
        )
        distance_rmse = solver_scalar(
            array_module.sqrt(array_module.mean(residual_squared))
        )
        if use_gpu:
            with cupyx.errstate(linalg="raise"):
                gram_eigenvalues = cp.linalg.eigvalsh(reduced_gram)
        else:
            gram_eigenvalues = np.linalg.eigvalsh(reduced_gram)
        logdet = solver_scalar(array_module.sum(array_module.log(gram_eigenvalues)))
        entropy = logdet - n_modes * np.log(3.0)
        reduced_connectivity = 3.0 * candidate_inverse
        basis_times_connectivity = solver_basis @ reduced_connectivity
        connectivity_diagonal = -array_module.sum(
            basis_times_connectivity * solver_basis, axis=1
        )
        offdiagonal_norm_squared = 0.5 * max(
            solver_scalar(array_module.sum(reduced_connectivity**2))
            - solver_scalar(array_module.sum(connectivity_diagonal**2)),
            0.0,
        )
        connectivity_norm = float(np.sqrt(offdiagonal_norm_squared))
        cg_relative_residual = cg_residual / max(
            gradient_norm, np.finfo(np.float64).tiny
        )

        history["iteration"].append(iteration)
        history["objective"].append(candidate_objective)
        history["negative_entropy_objective"].append(candidate_negative_entropy)
        history["data_objective"].append(candidate_data_objective)
        history["loss"].append(relative_loss)
        history["entropy"].append(entropy)
        history["weighted_rmse"].append(weighted_rmse)
        history["distance_rmse"].append(distance_rmse)
        history["gradient_norm"].append(candidate_gradient_norm)
        history["relative_gradient_norm"].append(candidate_relative_gradient)
        history["newton_decrement"].append(
            float(np.sqrt(max(-directional_derivative, 0.0)))
        )
        history["relative_step"].append(relative_step)
        history["step_size"].append(step_size)
        history["cg_iterations"].append(cg_iterations)
        history["cg_relative_residual"].append(cg_relative_residual)
        history["minimum_internal_gram_eigenvalue"].append(
            solver_scalar(gram_eigenvalues[0])
        )
        history["maximum_internal_gram_eigenvalue"].append(
            solver_scalar(gram_eigenvalues[-1])
        )
        history["connectivity_offdiagonal_l2"].append(connectivity_norm)

        if iteration in save_steps_set:
            checkpoint_reduced_gram = (
                cp.asnumpy(reduced_gram) if use_gpu else reduced_gram
            )
            connectivity_at_steps[iteration] = _connectivity_from_reduced_gram(
                checkpoint_reduced_gram, basis
            )
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_optimization",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": float(candidate_objective),
                    "loss": relative_loss,
                    "entropy": entropy,
                    "relative_gradient_norm": candidate_relative_gradient,
                    "step_size": step_size,
                    "cg_iterations": cg_iterations,
                    "cg_converged": bool(cg_converged),
                    "noisy": True,
                    "use_gpu": bool(use_gpu),
                    "method": "COV",
                    "general_method": "covariance_optimization",
                }
            )

        if candidate_gradient_norm <= (
            absolute_tolerance + relative_tolerance * candidate_gradient_scale
        ):
            converged = True
            status = "optimality_tolerance"
            message = "first-order optimality tolerance reached"
            break

    for key, values in history.items():
        dtype = np.int64 if key in {"iteration", "cg_iterations"} else np.float64
        history[key] = np.asarray(values, dtype=dtype)

    if use_gpu:
        reduced_gram = cp.asnumpy(reduced_gram)
    gram = basis @ reduced_gram @ basis.T
    gram = 0.5 * (gram + gram.T)
    fitted_squared_distances = _squared_distances_from_gram(gram)
    connectivity = _connectivity_from_reduced_gram(reduced_gram, basis)
    if history["iteration"].size == 0:
        final_state = current_state
        final_objective = float(final_state[0])
        final_gradient_norm = solver_norm(final_state[1])
        final_relative_gradient_norm = final_gradient_norm / max(
            1.0,
            solver_norm(final_state[4][2]),
            solver_norm(final_state[4][3]),
        )
    else:
        final_objective = float(history["objective"][-1])
        final_gradient_norm = float(history["gradient_norm"][-1])
        final_relative_gradient_norm = float(history["relative_gradient_norm"][-1])

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
        "iterations": (
            int(history["iteration"][-1]) if history["iteration"].size else 0
        ),
        "objective": final_objective,
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
        "gpu_device": get_gpu_name() if use_gpu else None,
        "cupy_version": cp.__version__ if use_gpu else None,
        "preconditioner_pair_block_size": (
            _COVARIANCE_PRECONDITIONER_PAIR_BLOCK_SIZE
        ),
        "preconditioner_setup_seconds": float(preconditioner_setup_seconds),
        "preconditioner_data_setup_count": 1,
        "preconditioner_entropy_diagonal_updated_each_iteration": True,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
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


def ddmap2cov(ddmap):
    # convert a squared distance map to covariance matrix
    n = ddmap.shape[0]
    omega2_mtx = ddmap / 3.0
    omega2_row_sum = np.sum(omega2_mtx, axis=1)
    omega2_sum = np.sum(omega2_mtx)
    return (omega2_row_sum[:, np.newaxis] + omega2_row_sum - omega2_sum / n) / (
        2 * n
    ) - omega2_mtx / 2.0


def dmap2cov(dmap):
    # convert a distance map to covariance matrix
    ddmap = (3.0 * np.pi / 8.0) * np.power(dmap, 2.0)
    return ddmap2cov(ddmap)


def checkEMD(ddmap, neg_tol=1e-10):
    # check whether a squared distance map is a valid EDM (avg of EDMs is EDM in some dim; tiny neg evals = numerical noise)
    cov = ddmap2cov(ddmap)
    eigvalue, eigvector = np.linalg.eigh(cov)
    min_eig = eigvalue.min()
    if min_eig < -0.1:  # unacceptably non-PSD
        return False
    if min_eig < -neg_tol:  # meaningfully negative but within tolerance → warn
        print(
            "[red]Warning: The smallest eigenvalue of the covariance matrix is negative. \
Direct inversion method [italic]may[/italic] not work. Check the final results. \
If the results are not good enough, please try iterative scaling or gradient descent method."
        )
    return True


def subnetwork_schur(A, keep, tol=1e-12):
    """
    Compute the marginal (Schur‐complement) connectivity matrix on a subset of nodes.

    Parameters
    ----------
    A : (N, N) array_like
        Full connectivity (precision) matrix.
    keep : 1D array_like of ints
        Indices of the nodes to keep (subset S).
    tol : float, optional
        Threshold below which A_CC is considered singular; uses pinv accordingly.

    Returns
    -------
    A_eff : (len(keep), len(keep)) ndarray
        Effective connectivity on the subset S,
        A_eff = A_SS - A_SC @ inv(A_CC) @ A_CS
    """
    A = np.asarray(A)
    keep = np.asarray(keep, dtype=int)
    N = A.shape[0]
    # build complement indices
    all_idx = np.arange(N)
    remove = np.setdiff1d(all_idx, keep)

    # partition
    A_SS = A[np.ix_(keep, keep)]
    A_SC = A[np.ix_(keep, remove)]
    A_CC = A[np.ix_(remove, remove)]
    A_CS = A_SC.T

    # invert (or pseudo-invert) A_CC
    # if it's near-singular, pinv is safer
    if np.linalg.cond(A_CC) > 1 / tol:
        A_CC_inv = np.linalg.pinv(A_CC, rcond=tol)
    else:
        A_CC_inv = np.linalg.inv(A_CC)

    # Schur complement
    A_eff = A_SS - A_SC @ A_CC_inv @ A_CS
    return A_eff


def neighbor_balance_symmetric(
    C, *, not_normalize=False, circular=False, epsilon=1e-12, return_scales=False
):
    """
    Symmetric neighbor balancing:

        s_i = 0.5 * (C[i, i-1] + C[i, i+1])         (with edge handling)
        C_bal[i, j] = C[i, j] / sqrt(s_i * s_j)

    Parameters
    ----------
    C : (N, N) array_like
        Contact matrix (typically symmetric, nonnegative).
    circular : bool, default False
        If True, treat indices modulo N (wrap-around neighbors).
        If False, use the single available neighbor at the ends:
            s_0   = C[0, 1]
            s_{N-1} = C[N-1, N-2]
    epsilon : float, default 1e-12
        Stabilizer to avoid division by zero and negative scales.
    return_scales : bool, default False
        If True, also return the vector s of per-locus neighbor averages.

    Returns
    -------
    C_bal : (N, N) ndarray
        Neighbor-balanced matrix using your symmetric formula.
    s : (N,) ndarray, optional
        The per-locus neighbor averages used (only if return_scales=True).
    """
    X = np.asarray(C, dtype=float)
    if X.ndim != 2 or X.shape[0] != X.shape[1]:
        raise ValueError("C must be a square 2D array")

    n = X.shape[0]
    idx = np.arange(n)

    if circular:
        left_idx = (idx - 1) % n
        right_idx = (idx + 1) % n
        s = 0.5 * (X[idx, left_idx] + X[idx, right_idx])
    else:
        s = np.empty(n, dtype=float)
        # interior i=1..n-2 have both neighbors
        if n >= 3:
            s[1:-1] = 0.5 * (np.diag(X, k=-1)[:-1] + np.diag(X, k=1)[1:])
        elif n == 2:
            # degenerate case: each has only one neighbor
            s[:] = np.array([X[0, 1], X[1, 0]])
        else:
            # n == 1: no neighbors; avoid divide by zero by using epsilon
            s[:] = epsilon

        # edges: use the single available neighbor
        if n >= 2:
            s[0] = X[0, 1]
            s[-1] = X[-1, -2]

    # Guard against zeros/negatives
    s = np.maximum(s, epsilon)
    scale = np.sqrt(s)

    C_bal = X / (scale[:, None] * scale[None, :])

    # multiply by the mean value of diagonal contacts if not_normalize is set to be True
    if not_normalize:
        C_bal = np.nanmean(np.diag(C)) * C_bal

    if return_scales:
        return C_bal, s
    return C_bal


# ------------------------------------------------------------------#


def compute_entropy_from_A(A, zero_tol=1e-12, eigvals=None):
    """Compute the entropy of the maximum entropy model from connectivity matrix A.

    The model is a multivariate Gaussian distribution. The entropy is:
        H = constant + log(det(K+))
    where K = -A is the positive semi-definite matrix, and K+ is its pseudo-inverse.

    For numerical stability, we compute:
        log(det(K+)) = sum_i log(1/λ_i) = -sum_i log(λ_i)
    for all non-zero eigenvalues λ_i of K.

    Parameters
    ----------
    A : (n,n) array
        Connectivity matrix (negative semi-definite, Laplacian-like).
    zero_tol : float, default=1e-12
        Tolerance for identifying zero eigenvalues.
    eigvals : array_like, optional
        Precomputed eigenvalues of K = -A. If provided, avoids recomputing eigenvalues.

    Returns
    -------
    entropy : float
        log(det(K+)) where K = -A and K+ is pseudo-inverse.
    """
    K = -A

    if eigvals is not None:
        eigvals = np.asarray(eigvals)
    else:
        eigvals = np.linalg.eigvalsh(K)

    positive_mask = eigvals > zero_tol
    positive_eigvals = eigvals[positive_mask]

    if len(positive_eigvals) == 0:
        max_eigval = np.max(eigvals) if len(eigvals) > 0 else 0.0
        if max_eigval < zero_tol:
            return -np.inf
        return np.nan

    with np.errstate(divide="ignore", invalid="ignore"):
        log_terms = -np.log(positive_eigvals)

    if np.any(~np.isfinite(log_terms)):
        return np.nan

    entropy = np.sum(log_terms)
    if not np.isfinite(entropy):
        return np.nan
    return entropy


# --------
# Sample HIPPS-DIMES structures conditioned on single pairwise distance
def _conditioned_pair_gain(K, pair):
    """Return the linear conditioning gain g for a selected locus pair."""
    i, j = pair
    N = K.shape[0]

    Sigma = scipy.linalg.pinvh(-K)

    delta = np.zeros(N)
    delta[i] = 1.0
    delta[j] = -1.0

    denom = delta @ Sigma @ delta
    if denom <= 0:
        raise ValueError("Invalid pair variance denominator.")

    return (Sigma @ delta) / denom


def sample_conditioned_pair_vector(r_eq, K, pair, b):
    """
    Exact conditional equilibrium sample for Gaussian HIPPS-DIMES,
    conditioned on r_ij = b, where

        r_ij = r_i - r_j = (delta^T r)^T.

    Parameters
    ----------
    r_eq : ndarray, shape (N, 3)
        Unconstrained equilibrium sample, i.e. r_eq.
    K : ndarray, shape (N, N)
        Connectivity matrix.
    pair : tuple of int
        (i, j)
    b : ndarray, shape (3,)
        Desired pair vector.

    Returns
    -------
    r_cond : ndarray, shape (N, 3)
        Exact conditional sample satisfying
            r_cond[i] - r_cond[j] = b.
    """
    i, j = pair
    g = _conditioned_pair_gain(K, pair)

    # current pair vector: r_ij = (delta^T r_eq)^T = r_eq[i] - r_eq[j]
    r_ij_eq = r_eq[i] - r_eq[j]  # shape (3,)

    # r_cond = r_eq + g (b - r_ij_eq)^T
    r_cond = r_eq + g[:, None] * (b - r_ij_eq)[None, :]

    return r_cond


def random_unit_vector():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)


def random_unit_vectors(count):
    """Return ``count`` random unit vectors in R^3."""
    v = np.random.normal(size=(count, 3))
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    norm = np.maximum(norm, np.finfo(float).tiny)
    return v / norm


def sample_conditioned_pair_distance(X_eq, K, pair, b_scalar):
    """
    Exact conditional equilibrium sample for standard isotropic HIPPS-DIMES,
    conditioned on |r_i - r_j| = R0.
    """
    n = random_unit_vector()
    b = b_scalar * n
    return sample_conditioned_pair_vector(X_eq, K, pair, b)


def sample_conditioned_pair_distance_batch(X_eq, K, pair, b_scalar):
    """
    Vectorized conditional sampling for an ensemble constrained on one pair distance.

    Parameters
    ----------
    X_eq : ndarray, shape (ensemble, N, 3)
        Unconstrained equilibrium samples.
    K : ndarray, shape (N, N)
        Connectivity matrix.
    pair : tuple of int
        Selected locus pair ``(i, j)``.
    b_scalar : float
        Desired pair distance.

    Returns
    -------
    ndarray, shape (ensemble, N, 3)
        Conditional samples with ``|r_i - r_j| = b_scalar`` for each sample.
    """
    X_eq = np.asarray(X_eq)
    if X_eq.ndim != 3 or X_eq.shape[2] != 3:
        raise ValueError(f"X_eq must have shape (ensemble, N, 3), got {X_eq.shape}")

    i, j = pair
    g = _conditioned_pair_gain(K, pair)
    b = b_scalar * random_unit_vectors(X_eq.shape[0])
    pair_vectors = X_eq[:, i, :] - X_eq[:, j, :]
    return X_eq + g[np.newaxis, :, np.newaxis] * (b - pair_vectors)[:, np.newaxis, :]


def a2xyz_sample_conditioned_pair_distance(K, pair, b_scalar, ensemble=1):
    xyzs_uncond = a2xyz_sample(K, ensemble=ensemble)
    return sample_conditioned_pair_distance_batch(xyzs_uncond, K, pair, b_scalar)
