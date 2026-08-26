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


def _squared_distances_from_gram(B):
    """Return the squared Euclidean distance matrix induced by Gram matrix ``B``."""
    diagonal = np.diag(B)
    ddmap = diagonal[:, np.newaxis] + diagonal - 2.0 * B
    ddmap = 0.5 * (ddmap + ddmap.T)
    np.fill_diagonal(ddmap, 0.0)
    return ddmap


def _project_centered_psd(X, gram_eigenvalue_floor=0.0):
    """Project onto ``B @ 1 = 0`` and ``B >= floor * J``."""
    X = 0.5 * (X + X.T)
    row_mean = np.mean(X, axis=1, keepdims=True)
    centered = X - row_mean - row_mean.T + np.mean(X)

    # The floored feasible set is a translation of the centered PSD cone:
    # {B : B >= floor * J, B @ 1 = 0} = floor * J + {C : C >= 0, C @ 1 = 0}.
    # Express J without matrix multiplication so the translational mode remains
    # exactly excluded from the positive spectral floor.
    n = centered.shape[0]
    shifted = centered.copy()
    if gram_eigenvalue_floor != 0.0:
        shifted += gram_eigenvalue_floor / n
        shifted[np.diag_indices(n)] -= gram_eigenvalue_floor

    eigenvalues, eigenvectors = np.linalg.eigh(shifted)
    projected = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T

    # Re-center by the congruence J @ projected @ J, expressed with means to
    # avoid two dense matrix multiplications. This preserves PSD analytically.
    row_mean = np.mean(projected, axis=1, keepdims=True)
    projected = projected - row_mean - row_mean.T + np.mean(projected)
    if gram_eigenvalue_floor != 0.0:
        projected -= gram_eigenvalue_floor / n
        projected[np.diag_indices(n)] += gram_eigenvalue_floor
    return 0.5 * (projected + projected.T)


def _nearest_edm_objective_gradient(B, squared_distances, weights):
    """Evaluate the unique-pair weighted EDM objective and its gradient."""
    fitted = _squared_distances_from_gram(B)
    residual = fitted - squared_distances
    weighted_residual = weights * residual

    # weights and residual are symmetric. Multiplication by 1/4 therefore
    # evaluates 1/2 * sum_{i < j} w_ij * residual_ij**2.
    objective = 0.25 * np.sum(weighted_residual * residual)
    gradient = np.diag(np.sum(weighted_residual, axis=1)) - weighted_residual
    gradient = 0.5 * (gradient + gradient.T)
    return float(objective), gradient, fitted


def _centered_orthonormal_basis(n):
    """Return a deterministic orthonormal basis for ``1.T @ x = 0``."""
    basis = np.zeros((n, n - 1), dtype=np.float64)
    for column in range(n - 1):
        denominator = np.sqrt((column + 1) * (column + 2))
        basis[: column + 1, column] = 1.0 / denominator
        basis[column + 1, column] = -(column + 1) / denominator
    return basis


def _rouse_reference_gram(n, spring_constant, basis=None):
    """Return the Rouse Gram matrix and its reduced-space representation."""
    if basis is None:
        basis = _centered_orthonormal_basis(n)
    connectivity = -construct_connectivity_matrix_rouse(n, spring_constant)
    reduced_connectivity = basis.T @ connectivity @ basis
    reduced_gram_inverse = reduced_connectivity / 3.0
    reduced_gram = np.linalg.inv(reduced_gram_inverse)
    gram = basis @ reduced_gram @ basis.T
    return 0.5 * (gram + gram.T), reduced_gram, reduced_gram_inverse


def _mean_rouse_kl(reduced_gram, reference_inverse, reference_logdet):
    """Return Gaussian KL divergence per centered internal mode."""
    sign, logdet = np.linalg.slogdet(reduced_gram)
    if sign <= 0.0:
        return np.inf
    n_modes = reduced_gram.shape[0]
    divergence = 0.5 * (
        np.sum(reference_inverse * reduced_gram.T) - logdet + reference_logdet - n_modes
    )
    roundoff = 1e-12 * max(1.0, float(n_modes))
    if divergence < 0.0 and divergence >= -roundoff:
        divergence = 0.0
    return float(divergence / n_modes)


def _rouse_kl_prox(trial, step_size, prior_coefficient, reference_inverse):
    """Apply the reduced-space proximal map for the Rouse KL penalty."""
    shifted = trial - step_size * prior_coefficient * reference_inverse
    shifted = 0.5 * (shifted + shifted.T)
    eigenvalues, eigenvectors = np.linalg.eigh(shifted)
    product = step_size * prior_coefficient
    discriminant = np.sqrt(eigenvalues * eigenvalues + 4.0 * product)
    proximal_eigenvalues = np.empty_like(eigenvalues)
    nonnegative = eigenvalues >= 0.0
    proximal_eigenvalues[nonnegative] = 0.5 * (
        eigenvalues[nonnegative] + discriminant[nonnegative]
    )
    proximal_eigenvalues[~nonnegative] = (
        2.0 * product / (discriminant[~nonnegative] - eigenvalues[~nonnegative])
    )
    proximal = (eigenvectors * proximal_eigenvalues) @ eigenvectors.T
    return 0.5 * (proximal + proximal.T)


def _nearest_edm_with_rouse_prior(
    target,
    pair_weights,
    pair_count,
    spring_constant,
    prior_weight,
    max_iterations,
    relative_tolerance,
    absolute_tolerance,
):
    """Solve the positive-definite closest-EDM problem with a Rouse KL prior."""
    n = target.shape[0]
    basis = _centered_orthonormal_basis(n)
    _, reference_reduced, reference_inverse = _rouse_reference_gram(
        n, spring_constant, basis
    )
    reference_sign, reference_logdet = np.linalg.slogdet(reference_reduced)
    if reference_sign <= 0.0:
        raise RuntimeError("the reduced Rouse reference Gram matrix is not positive")

    n_modes = n - 1
    prior_coefficient = prior_weight / (2.0 * n_modes)

    def data_objective_gradient(reduced_gram):
        gram = basis @ reduced_gram @ basis.T
        objective, gradient, fitted = _nearest_edm_objective_gradient(
            gram, target, pair_weights
        )
        return (
            objective / pair_count,
            (basis.T @ gradient @ basis) / pair_count,
            objective,
            fitted,
        )

    def total_objective(reduced_gram, normalized_data_objective):
        mean_kl = _mean_rouse_kl(reduced_gram, reference_inverse, reference_logdet)
        return (
            normalized_data_objective + prior_weight * mean_kl,
            mean_kl,
        )

    reduced_gram = reference_reduced.copy()
    extrapolated = reduced_gram.copy()
    momentum = 1.0
    current_data_objective, _, _, _ = data_objective_gradient(reduced_gram)
    current_total_objective, _ = total_objective(reduced_gram, current_data_objective)

    # The unnormalized data-term Hessian bound is divided by the observed
    # unordered-pair count to match the mean data objective.
    lipschitz = 4.0 * np.max(np.sum(pair_weights, axis=1)) / pair_count
    lipschitz = max(float(lipschitz), np.finfo(np.float64).tiny)
    weight_sum = 0.5 * np.sum(pair_weights)

    history = {
        "iteration": [],
        "objective": [],
        "normalized_data_objective": [],
        "weighted_rmse": [],
        "mean_rouse_kl": [],
        "rouse_prior_penalty": [],
        "total_objective": [],
        "proximal_gradient_norm": [],
        "relative_proximal_gradient_norm": [],
        "relative_step": [],
        "step_size": [],
        "restarted": [],
    }
    converged = False
    status = "max_iterations"

    for iteration in range(1, max_iterations + 1):
        restarted = False
        while True:
            data_objective, data_gradient, _, _ = data_objective_gradient(extrapolated)

            for _ in range(60):
                step_size = 1.0 / lipschitz
                candidate = _rouse_kl_prox(
                    extrapolated - step_size * data_gradient,
                    step_size,
                    prior_coefficient,
                    reference_inverse,
                )
                (
                    candidate_data_objective,
                    candidate_data_gradient,
                    candidate_raw_objective,
                    _,
                ) = data_objective_gradient(candidate)
                proximal_delta = candidate - extrapolated
                majorizer = (
                    data_objective
                    + np.sum(data_gradient * proximal_delta)
                    + 0.5 * lipschitz * np.sum(proximal_delta * proximal_delta)
                )
                slack = 1e-12 * max(
                    1.0, abs(data_objective), abs(candidate_data_objective)
                )
                if candidate_data_objective <= majorizer + slack:
                    break
                lipschitz *= 2.0
            else:
                raise RuntimeError("nearest_edm Rouse-prior backtracking failed")

            candidate_total_objective, candidate_mean_kl = total_objective(
                candidate, candidate_data_objective
            )
            monotone_slack = 1e-12 * max(
                1.0,
                abs(current_total_objective),
                abs(candidate_total_objective),
            )
            if candidate_total_objective <= current_total_objective + monotone_slack:
                break
            if restarted:
                raise RuntimeError("nearest_edm Rouse-prior monotone restart failed")
            extrapolated = reduced_gram
            momentum = 1.0
            restarted = True

        certificate = _rouse_kl_prox(
            candidate - step_size * candidate_data_gradient,
            step_size,
            prior_coefficient,
            reference_inverse,
        )
        proximal_gradient_norm = lipschitz * np.linalg.norm(
            candidate - certificate, ord="fro"
        )
        candidate_inverse = np.linalg.inv(candidate)
        prior_gradient = prior_coefficient * (reference_inverse - candidate_inverse)
        optimality_scale = max(
            np.linalg.norm(candidate_data_gradient, ord="fro"),
            np.linalg.norm(prior_gradient, ord="fro"),
            np.finfo(np.float64).tiny,
        )
        relative_proximal_gradient_norm = proximal_gradient_norm / optimality_scale
        relative_step = np.linalg.norm(candidate - reduced_gram, ord="fro") / max(
            np.linalg.norm(reduced_gram, ord="fro"), 1.0
        )
        weighted_rmse = np.sqrt(2.0 * candidate_raw_objective / weight_sum)

        history["iteration"].append(iteration)
        history["objective"].append(candidate_raw_objective)
        history["normalized_data_objective"].append(candidate_data_objective)
        history["weighted_rmse"].append(weighted_rmse)
        history["mean_rouse_kl"].append(candidate_mean_kl)
        history["rouse_prior_penalty"].append(prior_weight * candidate_mean_kl)
        history["total_objective"].append(candidate_total_objective)
        history["proximal_gradient_norm"].append(proximal_gradient_norm)
        history["relative_proximal_gradient_norm"].append(
            relative_proximal_gradient_norm
        )
        history["relative_step"].append(relative_step)
        history["step_size"].append(step_size)
        history["restarted"].append(restarted)

        previous_reduced_gram = reduced_gram
        reduced_gram = candidate
        current_total_objective = candidate_total_objective
        if proximal_gradient_norm <= (
            absolute_tolerance + relative_tolerance * optimality_scale
        ):
            converged = True
            status = "optimality_tolerance"
            break

        next_momentum = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * momentum**2))
        extrapolated = reduced_gram + ((momentum - 1.0) / next_momentum) * (
            reduced_gram - previous_reduced_gram
        )
        momentum = next_momentum
        lipschitz *= 0.9

    for key, values in history.items():
        history[key] = np.asarray(values)

    gram = basis @ reduced_gram @ basis.T
    gram = 0.5 * (gram + gram.T)
    fitted_squared_distances = _squared_distances_from_gram(gram)
    final = {
        key: float(values[-1])
        for key, values in history.items()
        if key not in {"iteration", "restarted"}
    }
    info = {
        "converged": converged,
        "status": status,
        "message": (
            "proximal-gradient optimality tolerance reached"
            if converged
            else "maximum number of iterations reached"
        ),
        "iterations": int(history["iteration"][-1]),
        "objective": final["objective"],
        "normalized_data_objective": final["normalized_data_objective"],
        "weighted_rmse": final["weighted_rmse"],
        "mean_rouse_kl": final["mean_rouse_kl"],
        "rouse_prior_penalty": final["rouse_prior_penalty"],
        "total_objective": final["total_objective"],
        "proximal_gradient_norm": final["proximal_gradient_norm"],
        "relative_proximal_gradient_norm": final["relative_proximal_gradient_norm"],
        "gram_eigenvalue_floor": 0.0,
        "rouse_prior_weight": prior_weight,
        "rouse_spring_constant": spring_constant,
        "observed_pair_count": int(pair_count),
        "weight_sum": float(weight_sum),
        "history": history,
    }
    if not converged:
        warnings.warn(
            "nearest_edm reached max_iterations before satisfying the "
            "proximal-gradient optimality tolerance",
            RuntimeWarning,
            stacklevel=3,
        )
    return fitted_squared_distances, gram, info


def nearest_edm(
    squared_distances,
    weights=None,
    *,
    gram_eigenvalue_floor=0.0,
    rouse_prior_weight=0.0,
    rouse_spring_constant=None,
    max_iterations=1000,
    relative_tolerance=1e-8,
    absolute_tolerance=1e-10,
):
    """Find a weighted least-squares Euclidean distance matrix.

    This function solves

    ``min_B 0.5 * sum_{i < j} w_ij * (D(B)_ij - D_obs_ij)**2``

    subject to ``B @ 1 = 0`` and ``B >= gram_eigenvalue_floor * J``, where
    ``J = I - 11.T / N`` and
    ``D(B)_ij = B_ii + B_jj - 2 * B_ij``. The input and returned distance
    matrices contain *squared* distances. A positive floor gives exactly
    ``N - 1`` positive internal Gram modes while preserving the translational
    zero mode. Alternatively, a positive ``rouse_prior_weight`` adds a
    normalized Gaussian KL penalty relative to a full-rank free Rouse chain.
    The hard floor and soft prior cannot be enabled together.

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
    rouse_prior_weight : float, optional
        Nonnegative dimensionless weight for the mean per-mode Gaussian KL
        divergence from a free Rouse-chain Gram matrix. With a positive value,
        the data objective is divided by the number of observed unordered
        pairs. Default is 0, which retains the projected-gradient algorithm.
    rouse_spring_constant : float, optional
        Positive Rouse spring constant in reduced units with ``kBT = 1``.
        When the prior is active and this is omitted, use ``3 / median(Delta_i,
        i+1)`` over finite, positive observed adjacent squared distances.
    max_iterations : int, optional
        Maximum number of projected-gradient iterations.
    relative_tolerance, absolute_tolerance : float, optional
        Stop when the applicable projected- or proximal-gradient norm is no
        larger than an absolute plus relative optimality scale.

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
    iteration and memory use is O(N^2).
    """
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
    if isinstance(rouse_prior_weight, (bool, np.bool_)) or not np.isscalar(
        rouse_prior_weight
    ):
        raise ValueError("rouse_prior_weight must be a finite nonnegative scalar")
    try:
        rouse_prior_weight = float(rouse_prior_weight)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "rouse_prior_weight must be a finite nonnegative scalar"
        ) from error
    if not np.isfinite(rouse_prior_weight) or rouse_prior_weight < 0.0:
        raise ValueError("rouse_prior_weight must be a finite nonnegative scalar")
    if rouse_spring_constant is not None:
        if isinstance(rouse_spring_constant, (bool, np.bool_)) or not np.isscalar(
            rouse_spring_constant
        ):
            raise ValueError("rouse_spring_constant must be a finite positive scalar")
        try:
            rouse_spring_constant = float(rouse_spring_constant)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "rouse_spring_constant must be a finite positive scalar"
            ) from error
        if not np.isfinite(rouse_spring_constant) or rouse_spring_constant <= 0.0:
            raise ValueError("rouse_spring_constant must be a finite positive scalar")
    if rouse_prior_weight > 0.0 and gram_eigenvalue_floor > 0.0:
        raise ValueError(
            "gram_eigenvalue_floor and rouse_prior_weight cannot both be positive"
        )
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")

    target = np.zeros_like(observed)
    target[pair_mask] = observed[pair_mask]
    target = 0.5 * (target + target.T)

    weight_sum = 0.5 * np.sum(pair_weights)
    pair_count = int(np.count_nonzero(np.triu(pair_mask, k=1)))
    if rouse_prior_weight > 0.0:
        if rouse_spring_constant is None:
            adjacent = target[np.arange(n - 1), np.arange(1, n)]
            adjacent_observed = pair_mask[np.arange(n - 1), np.arange(1, n)]
            adjacent = adjacent[adjacent_observed & (adjacent > 0.0)]
            if adjacent.size == 0:
                raise ValueError(
                    "rouse_spring_constant cannot be inferred without a finite "
                    "positive observed adjacent squared distance"
                )
            rouse_spring_constant = 3.0 / float(np.median(adjacent))
        return _nearest_edm_with_rouse_prior(
            target,
            pair_weights,
            pair_count,
            rouse_spring_constant,
            rouse_prior_weight,
            max_iterations,
            relative_tolerance,
            absolute_tolerance,
        )

    B = _project_centered_psd(np.zeros_like(observed), gram_eigenvalue_floor)
    extrapolated = B.copy()
    momentum = 1.0
    current_objective, _, _ = _nearest_edm_objective_gradient(B, target, pair_weights)

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
                extrapolated, target, pair_weights
            )

            for _ in range(60):
                step_size = 1.0 / lipschitz
                candidate = _project_centered_psd(
                    extrapolated - step_size * gradient,
                    gram_eigenvalue_floor,
                )
                candidate_objective, candidate_gradient, _ = (
                    _nearest_edm_objective_gradient(candidate, target, pair_weights)
                )
                proximal_delta = candidate - extrapolated
                proximal_delta_norm_squared = np.sum(proximal_delta * proximal_delta)
                majorizer = (
                    objective
                    + np.sum(gradient * proximal_delta)
                    + 0.5 * lipschitz * proximal_delta_norm_squared
                )
                slack = 1e-12 * max(1.0, abs(objective), abs(candidate_objective))
                if candidate_objective <= majorizer + slack:
                    break
                lipschitz *= 2.0
            else:
                raise RuntimeError("nearest_edm backtracking line search failed")

            monotone_slack = 1e-12 * max(
                1.0, abs(current_objective), abs(candidate_objective)
            )
            if candidate_objective <= current_objective + monotone_slack:
                break
            if restarted:
                raise RuntimeError("nearest_edm monotone restart failed")

            # Restart acceleration from the last accepted feasible iterate.
            extrapolated = B
            momentum = 1.0
            restarted = True

        # Evaluate the projected-gradient mapping at the returned feasible
        # candidate, rather than using relative iterate change as a proxy for
        # first-order optimality.
        certificate = _project_centered_psd(
            candidate - step_size * candidate_gradient,
            gram_eigenvalue_floor,
        )
        certificate_delta = candidate - certificate
        projected_gradient_norm = lipschitz * np.linalg.norm(
            certificate_delta, ord="fro"
        )
        gradient_norm = np.linalg.norm(candidate_gradient, ord="fro")
        relative_projected_gradient_norm = projected_gradient_norm / max(
            gradient_norm, np.finfo(np.float64).tiny
        )
        accepted_delta = candidate - B
        relative_step = np.linalg.norm(accepted_delta, ord="fro") / max(
            np.linalg.norm(B, ord="fro"), 1.0
        )
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

        # Try a modestly larger step on the next iteration. Backtracking will
        # restore the majorization bound when this is too aggressive.
        lipschitz *= 0.9

    for key, values in history.items():
        history[key] = np.asarray(values)

    fitted_squared_distances = _squared_distances_from_gram(B)
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
        "relative_projected_gradient_norm": (final_relative_projected_gradient_norm),
        "gram_eigenvalue_floor": gram_eigenvalue_floor,
        "weight_sum": float(weight_sum),
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


def _validate_gaussian_noise_variance(squared_distances, noise_variance):
    """Return observed-pair indices and variances for Gaussian soft constraints."""
    observed = np.asarray(squared_distances, dtype=np.float64)
    if observed.ndim != 2 or observed.shape[0] != observed.shape[1]:
        raise ValueError("squared_distances must be a square matrix")
    if observed.shape[0] < 2:
        raise ValueError("squared_distances must contain at least two loci")
    if np.any(np.isinf(observed)):
        raise ValueError("squared_distances must not contain infinite values")

    n = observed.shape[0]
    upper = np.triu_indices(n, k=1)
    finite = np.isfinite(observed)
    if not np.array_equal(finite, finite.T):
        raise ValueError("the observed-pair pattern must be symmetric")
    pair_mask = finite & ~np.eye(n, dtype=bool)
    if not np.any(np.triu(pair_mask, k=1)):
        raise ValueError("at least one off-diagonal pair must be observed")
    if not np.allclose(
        observed[pair_mask], observed.T[pair_mask], rtol=1e-10, atol=1e-12
    ):
        raise ValueError("observed squared distances must be symmetric")
    if np.any(observed[pair_mask] <= 0.0):
        raise ValueError("observed off-diagonal squared distances must be positive")

    if isinstance(noise_variance, (bool, np.bool_)):
        raise ValueError("noise_variance must be positive and finite")
    if np.isscalar(noise_variance):
        try:
            scalar_variance = float(noise_variance)
        except (TypeError, ValueError) as error:
            raise ValueError("noise_variance must be positive and finite") from error
        if not np.isfinite(scalar_variance) or scalar_variance <= 0.0:
            raise ValueError("noise_variance must be positive and finite")
        pair_variance = np.full(
            int(np.count_nonzero(np.triu(pair_mask, k=1))),
            scalar_variance,
            dtype=np.float64,
        )
        variance_kind = "scalar"
    else:
        variance_matrix = np.asarray(noise_variance, dtype=np.float64)
        if variance_matrix.shape != observed.shape:
            raise ValueError(
                "noise_variance matrix must have the same shape as squared_distances"
            )
        if not np.all(np.isfinite(variance_matrix)):
            raise ValueError("noise_variance matrix must be finite")
        if np.any(variance_matrix < 0.0):
            raise ValueError("noise_variance matrix must be nonnegative")
        if not np.allclose(
            variance_matrix, variance_matrix.T, rtol=1e-12, atol=1e-14
        ):
            raise ValueError("noise_variance matrix must be symmetric")
        variance_matrix = 0.5 * (variance_matrix + variance_matrix.T)
        pair_variance = variance_matrix[upper][np.triu(pair_mask, k=1)[upper]]
        if np.any(pair_variance <= 0.0):
            raise ValueError(
                "noise_variance must be positive for every observed pair"
            )
        variance_kind = "matrix"

    pair_selector = np.triu(pair_mask, k=1)[upper]
    pair_i = upper[0][pair_selector]
    pair_j = upper[1][pair_selector]
    return observed, pair_mask, pair_i, pair_j, pair_variance, variance_kind


def _validate_connectivity_l1(connectivity_l1):
    """Return a finite nonnegative scalar connectivity L1 coefficient."""
    if isinstance(connectivity_l1, (bool, np.bool_)) or not np.isscalar(
        connectivity_l1
    ):
        raise ValueError("connectivity_l1 must be a nonnegative finite scalar")
    try:
        coefficient = float(connectivity_l1)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "connectivity_l1 must be a nonnegative finite scalar"
        ) from error
    if not np.isfinite(coefficient) or coefficient < 0.0:
        raise ValueError("connectivity_l1 must be a nonnegative finite scalar")
    return coefficient


def _soft_threshold(values, threshold):
    """Apply the elementwise signed soft-threshold operator."""
    values = np.asarray(values, dtype=np.float64)
    return np.sign(values) * np.maximum(np.abs(values) - threshold, 0.0)


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
    connectivity_l1=0.0,
):
    """Evaluate the Gaussian/L1 data term through the structured EDM operator."""
    reduced_gram = np.asarray(reduced_gram, dtype=np.float64)
    gram = basis @ reduced_gram @ basis.T
    fitted_matrix = _squared_distances_from_gram(gram)
    residual_matrix = fitted_matrix - target_matrix
    if connectivity_l1 == 0.0:
        effective_residual = residual_matrix
    else:
        effective_residual = _soft_threshold(
            residual_matrix, 2.0 * connectivity_l1
        )
    weighted_residual = inverse_variance_matrix * effective_residual
    # Each observed unordered pair appears twice in the symmetric matrices.
    data_objective = 0.25 * float(
        np.sum(weighted_residual * effective_residual)
    )
    full_gradient = np.diag(np.sum(weighted_residual, axis=1)) - weighted_residual
    data_gradient = basis.T @ full_gradient @ basis
    data_gradient = 0.5 * (data_gradient + data_gradient.T)
    return (
        data_objective,
        data_gradient,
        fitted_matrix[pair_i, pair_j],
        residual_matrix[pair_i, pair_j],
    )


def _gaussian_covariance_data_hessian_action(
    direction, basis, inverse_variance_matrix
):
    """Apply the Gaussian data Hessian without materializing pair vectors."""
    direction = np.asarray(direction, dtype=np.float64)
    direction = 0.5 * (direction + direction.T)
    full_direction = basis @ direction @ basis.T
    distance_direction = _squared_distances_from_gram(full_direction)
    weighted_direction = inverse_variance_matrix * distance_direction
    full_action = np.diag(np.sum(weighted_direction, axis=1)) - weighted_direction
    reduced_action = basis.T @ full_action @ basis
    return 0.5 * (reduced_action + reduced_action.T)


def _gaussian_covariance_objective_gradient(
    reduced_gram,
    basis,
    target_matrix,
    inverse_variance_matrix,
    pair_i,
    pair_j,
    connectivity_l1=0.0,
):
    """Evaluate the Gaussian soft-constraint objective in the covariance cone."""
    reduced_gram = np.asarray(reduced_gram, dtype=np.float64)
    try:
        cholesky_factor = np.linalg.cholesky(reduced_gram)
    except np.linalg.LinAlgError:
        return np.inf, None, None, None, None

    logdet = 2.0 * float(np.sum(np.log(np.diag(cholesky_factor))))
    inverse_gram = scipy.linalg.cho_solve(
        (cholesky_factor, True),
        np.eye(reduced_gram.shape[0], dtype=np.float64),
        check_finite=False,
    )
    inverse_gram = 0.5 * (inverse_gram + inverse_gram.T)
    data_objective, data_gradient, fitted_pairs, _ = (
        _gaussian_covariance_data_objective_gradient(
            reduced_gram,
            basis,
            target_matrix,
            inverse_variance_matrix,
            pair_i,
            pair_j,
            connectivity_l1,
        )
    )
    negative_entropy = -1.5 * float(logdet)
    entropy_gradient = -1.5 * inverse_gram
    gradient = data_gradient + entropy_gradient
    gradient = 0.5 * (gradient + gradient.T)
    return (
        negative_entropy + data_objective,
        gradient,
        fitted_pairs,
        inverse_gram,
        (negative_entropy, data_objective, entropy_gradient, data_gradient),
    )


def _proximal_gaussian_negative_entropy(matrix, step_size):
    """Apply the analytic proximal map of ``-3/2 logdet`` to a symmetric matrix."""
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("matrix must be finite")
    if not np.isfinite(step_size) or step_size <= 0.0:
        raise ValueError("step_size must be positive and finite")

    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    radius = np.hypot(
        eigenvalues, np.sqrt(6.0) * np.sqrt(step_size)
    )
    updated_eigenvalues = np.empty_like(eigenvalues)
    nonnegative = eigenvalues >= 0.0
    updated_eigenvalues[nonnegative] = (
        0.5 * eigenvalues[nonnegative] + 0.5 * radius[nonnegative]
    )
    updated_eigenvalues[~nonnegative] = 3.0 * (
        step_size / radius[~nonnegative]
    ) / (
        1.0 - eigenvalues[~nonnegative] / radius[~nonnegative]
    )
    if (
        not np.all(np.isfinite(updated_eigenvalues))
        or np.any(updated_eigenvalues <= 0.0)
    ):
        raise RuntimeError("log-determinant proximal map lost positive definiteness")
    proximal = (eigenvectors * updated_eigenvalues) @ eigenvectors.T
    proximal = 0.5 * (proximal + proximal.T)
    return proximal, eigenvalues, updated_eigenvalues


def _preconditioned_conjugate_gradient(
    operator, right_hand_side, preconditioner_diagonal, relative_tolerance, max_iterations
):
    """Solve a symmetric positive-definite matrix equation without materializing it."""
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


def _connectivity_from_reduced_gram(reduced_gram, basis):
    """Return the physical connectivity associated with an internal EDM Gram matrix."""
    reduced_connectivity = 3.0 * np.linalg.solve(
        reduced_gram, np.eye(reduced_gram.shape[0], dtype=np.float64)
    )
    connectivity = -(basis @ reduced_connectivity @ basis.T)
    connectivity = 0.5 * (connectivity + connectivity.T)
    # Rebuild the dependent diagonal so rows sum to zero to machine precision.
    return a2a(connectivity)


def _initialize_gaussian_reduced_gram(
    observed,
    pair_mask,
    pair_i,
    pair_j,
    inverse_variance,
    basis,
    initial_connectivity,
    initial_gram_floor_relative,
):
    """Return a shared physical initialization for Gaussian-noise solvers."""
    n = observed.shape[0]
    n_modes = n - 1
    if initial_connectivity is None:
        normalized_weights = np.zeros_like(observed)
        normalized_pair_weights = inverse_variance / float(
            np.mean(inverse_variance)
        )
        normalized_weights[pair_i, pair_j] = normalized_pair_weights
        normalized_weights[pair_j, pair_i] = normalized_pair_weights
        projection_target = np.asarray(observed, dtype=np.float64).copy()
        projection_target[~pair_mask] = np.nan
        np.fill_diagonal(projection_target, 0.0)
        _, closest_gram, initializer = nearest_edm(
            projection_target,
            normalized_weights,
            max_iterations=2000,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
        )
        gram_scale = float(np.trace(closest_gram) / n_modes)
        if not np.isfinite(gram_scale) or gram_scale <= 0.0:
            raise RuntimeError("closest-EDM initialization has a nonpositive scale")
        initial_floor = initial_gram_floor_relative * gram_scale
        closest_gram = _project_centered_psd(closest_gram, initial_floor)
        reduced_gram = basis.T @ closest_gram @ basis
        initialization = {
            "kind": "weighted_closest_edm",
            "closest_edm_converged": bool(initializer["converged"]),
            "closest_edm_status": initializer["status"],
            "closest_edm_iterations": int(initializer["iterations"]),
            "closest_edm_weighted_rmse": float(initializer["weighted_rmse"]),
            "gram_scale": gram_scale,
            "gram_eigenvalue_floor": initial_floor,
        }
    else:
        connectivity = np.asarray(initial_connectivity, dtype=np.float64)
        if connectivity.shape != (n, n) or not np.all(np.isfinite(connectivity)):
            raise ValueError("initial_connectivity must be a finite NxN matrix")
        if not np.allclose(
            connectivity, connectivity.T, rtol=1e-10, atol=1e-12
        ):
            raise ValueError("initial_connectivity must be symmetric")
        row_scale = max(float(np.max(np.abs(connectivity))), 1.0)
        row_sum_error = float(np.max(np.abs(np.sum(connectivity, axis=1))))
        if row_sum_error > 1e-10 * row_scale:
            raise ValueError("initial_connectivity rows must sum to zero")
        reduced_connectivity = basis.T @ (-connectivity) @ basis
        eigenvalues = np.linalg.eigvalsh(reduced_connectivity)
        if float(eigenvalues[0]) <= 0.0:
            raise ValueError(
                "initial_connectivity must be stable in every internal mode"
            )
        reduced_gram = 3.0 * np.linalg.solve(
            reduced_connectivity, np.eye(n_modes, dtype=np.float64)
        )
        initialization = {
            "kind": "provided_connectivity",
            "minimum_internal_stiffness": float(eigenvalues[0]),
            "maximum_internal_stiffness": float(eigenvalues[-1]),
        }

    reduced_gram = 0.5 * (reduced_gram + reduced_gram.T)
    try:
        np.linalg.cholesky(reduced_gram)
    except np.linalg.LinAlgError as error:
        raise RuntimeError(
            "initial covariance is not strictly positive definite"
        ) from error
    return reduced_gram, initialization


def fit_gaussian_noise_covariance(
    squared_distances,
    noise_variance,
    *,
    initial_connectivity=None,
    max_iterations=100,
    relative_tolerance=1e-8,
    absolute_tolerance=1e-10,
    cg_relative_tolerance=1e-4,
    cg_max_iterations=None,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
):
    """Fit Gaussian-noisy squared-distance constraints inside the covariance cone.

    The fitted centered EDM Gram matrix ``B`` minimizes

    ``-3/2 logdet(B) + 1/2 sum_{i<j} (D(B)_ij-Dobs_ij)^2 / variance_ij``

    over matrices that are strictly positive definite on the centered subspace.
    This is the convex soft-constraint maximum-entropy objective for independent
    Gaussian errors on the observed mean-squared distances. A damped Newton-CG
    iteration with an Armijo line search keeps every accepted iterate in the
    covariance cone. The default initialization is the existing weighted
    closest-EDM solution with a small internal spectral floor.

    Parameters
    ----------
    squared_distances : (N, N) array_like
        Symmetric mean-squared distance observations. NaNs may mark missing pairs.
    noise_variance : float or (N, N) array_like
        Positive Gaussian variance. A matrix specifies heteroscedastic pair
        variances; its diagonal is ignored.
    initial_connectivity : (N, N) array_like, optional
        Stable centered connectivity used instead of closest-EDM initialization.
    max_iterations : int, default=100
        Maximum accepted damped-Newton iterations.
    relative_tolerance, absolute_tolerance : float
        First-order optimality tolerances.
    cg_relative_tolerance : float, default=1e-4
        Maximum relative residual requested from each Newton-CG solve.
    cg_max_iterations : int, optional
        Maximum preconditioned-CG iterations. Defaults to twice the number of
        internal modes.
    initial_gram_floor_relative : float, default=1e-8
        Spectral floor relative to the closest-EDM mean internal eigenvalue.
    save_steps : iterable of int, optional
        Accepted iterations whose connectivity matrices should be returned.
    progress_callback : callable, optional
        Receives a dictionary after each accepted iteration.

    Returns
    -------
    fitted_squared_distances, gram_matrix, connectivity_matrix, info
        ``info`` contains the convergence certificate, scalar iteration history,
        initialization record, and any requested connectivity checkpoints.
    """
    (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        pair_variance,
        variance_kind,
    ) = _validate_gaussian_noise_variance(squared_distances, noise_variance)

    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(cg_relative_tolerance) or not (
        0.0 < cg_relative_tolerance < 1.0
    ):
        raise ValueError("cg_relative_tolerance must lie strictly between zero and one")
    if not np.isfinite(initial_gram_floor_relative) or (
        initial_gram_floor_relative <= 0.0
    ):
        raise ValueError("initial_gram_floor_relative must be positive and finite")

    n = observed.shape[0]
    n_modes = n - 1
    if cg_max_iterations is None:
        cg_max_iterations = 2 * n_modes
    if not isinstance(cg_max_iterations, (int, np.integer)) or cg_max_iterations <= 0:
        raise ValueError("cg_max_iterations must be a positive integer")

    save_steps_set = set()
    if save_steps is not None:
        for step in save_steps:
            if not isinstance(step, (int, np.integer)) or not (1 <= step <= max_iterations):
                raise ValueError("save_steps must contain integers between 1 and max_iterations")
            save_steps_set.add(int(step))

    basis = _centered_orthonormal_basis(n)
    pair_vectors = basis[pair_i] - basis[pair_j]
    target_pairs = observed[pair_i, pair_j]
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = (
        _gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    pair_vectors_squared = np.square(pair_vectors)
    data_hessian_diagonal = pair_vectors_squared.T @ (
        inverse_variance[:, np.newaxis] * pair_vectors_squared
    )
    del pair_vectors, pair_vectors_squared

    reduced_gram, initialization = _initialize_gaussian_reduced_gram(
        observed,
        pair_mask,
        pair_i,
        pair_j,
        inverse_variance,
        basis,
        initial_connectivity,
        initial_gram_floor_relative,
    )

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
    pair_count = len(target_pairs)
    current_state = _gaussian_covariance_objective_gradient(
        reduced_gram,
        basis,
        target_matrix,
        inverse_variance_matrix,
        pair_i,
        pair_j,
    )

    for iteration in range(1, max_iterations + 1):
        (
            objective,
            gradient,
            fitted_pairs,
            inverse_gram,
            components,
        ) = current_state
        (
            negative_entropy,
            data_objective,
            entropy_gradient,
            data_gradient,
        ) = components
        gradient_norm = float(np.linalg.norm(gradient, ord="fro"))
        gradient_scale = max(
            1.0,
            float(np.linalg.norm(entropy_gradient, ord="fro")),
            float(np.linalg.norm(data_gradient, ord="fro")),
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
                direction, basis, inverse_variance_matrix
            )
            entropy_term = 1.5 * inverse_gram @ direction @ inverse_gram
            result = data_term + entropy_term
            return 0.5 * (result + result.T)

        preconditioner_diagonal = data_hessian_diagonal + 1.5 * np.outer(
            np.diag(inverse_gram), np.diag(inverse_gram)
        )
        preconditioner_diagonal = np.maximum(
            preconditioner_diagonal, np.finfo(np.float64).tiny
        )
        forcing_tolerance = min(
            cg_relative_tolerance,
            max(1e-8, 0.1 * np.sqrt(max(relative_gradient_norm, 0.0))),
        )
        direction, cg_iterations, cg_residual, cg_converged = (
            _preconditioned_conjugate_gradient(
                hessian_operator,
                -gradient,
                preconditioner_diagonal,
                forcing_tolerance,
                int(cg_max_iterations),
            )
        )
        directional_derivative = float(np.sum(gradient * direction))
        if not np.isfinite(directional_derivative) or directional_derivative >= 0.0:
            direction = -gradient / preconditioner_diagonal
            direction = 0.5 * (direction + direction.T)
            directional_derivative = float(np.sum(gradient * direction))
            cg_converged = False

        step_size = 1.0
        candidate = None
        candidate_state = None
        for _ in range(60):
            trial = reduced_gram + step_size * direction
            trial = 0.5 * (trial + trial.T)
            trial_state = _gaussian_covariance_objective_gradient(
                trial,
                basis,
                target_matrix,
                inverse_variance_matrix,
                pair_i,
                pair_j,
            )
            trial_objective = trial_state[0]
            if not np.isfinite(trial_objective):
                step_size *= 0.5
                continue
            if trial_objective <= objective + 1e-4 * step_size * directional_derivative:
                candidate = trial
                candidate_state = trial_state
                break
            step_size *= 0.5
        if candidate is None:
            status = "line_search_failed"
            message = "covariance-cone Armijo line search failed"
            break

        accepted_delta = candidate - reduced_gram
        relative_step = float(np.linalg.norm(accepted_delta, ord="fro")) / max(
            float(np.linalg.norm(reduced_gram, ord="fro")), 1.0
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
        candidate_gradient_norm = float(np.linalg.norm(candidate_gradient, ord="fro"))
        candidate_gradient_scale = max(
            1.0,
            float(np.linalg.norm(candidate_entropy_gradient, ord="fro")),
            float(np.linalg.norm(candidate_data_gradient, ord="fro")),
        )
        candidate_relative_gradient = candidate_gradient_norm / candidate_gradient_scale
        residual = candidate_pairs - target_pairs
        relative_loss = float(
            np.sqrt(np.mean(np.square(residual / target_pairs)))
        )
        weighted_rmse = float(
            np.sqrt(np.mean(np.square(residual) * inverse_variance))
        )
        distance_rmse = float(np.sqrt(np.mean(np.square(residual))))
        gram_eigenvalues = np.linalg.eigvalsh(reduced_gram)
        logdet = float(np.sum(np.log(gram_eigenvalues)))
        entropy = logdet - n_modes * np.log(3.0)
        reduced_connectivity = 3.0 * candidate_inverse
        full_connectivity = -(basis @ reduced_connectivity @ basis.T)
        offdiagonal = full_connectivity[np.triu_indices(n, k=1)]
        connectivity_norm = float(np.linalg.norm(offdiagonal))
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
            float(gram_eigenvalues[0])
        )
        history["maximum_internal_gram_eigenvalue"].append(
            float(gram_eigenvalues[-1])
        )
        history["connectivity_offdiagonal_l2"].append(connectivity_norm)

        if iteration in save_steps_set:
            connectivity_at_steps[iteration] = _connectivity_from_reduced_gram(
                reduced_gram, basis
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
                    "use_gpu": False,
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

    gram = basis @ reduced_gram @ basis.T
    gram = 0.5 * (gram + gram.T)
    fitted_squared_distances = _squared_distances_from_gram(gram)
    connectivity = _connectivity_from_reduced_gram(reduced_gram, basis)
    if history["iteration"].size == 0:
        final_state = current_state
        final_objective = float(final_state[0])
        final_gradient_norm = float(np.linalg.norm(final_state[1], ord="fro"))
        final_relative_gradient_norm = final_gradient_norm / max(
            1.0,
            float(np.linalg.norm(final_state[4][2], ord="fro")),
            float(np.linalg.norm(final_state[4][3], ord="fro")),
        )
    else:
        final_objective = float(history["objective"][-1])
        final_gradient_norm = float(history["gradient_norm"][-1])
        final_relative_gradient_norm = float(history["relative_gradient_norm"][-1])

    info = {
        "converged": converged,
        "status": status,
        "message": message,
        "iterations": int(history["iteration"][-1]) if history["iteration"].size else 0,
        "objective": final_objective,
        "gradient_norm": final_gradient_norm,
        "relative_gradient_norm": final_relative_gradient_norm,
        "observed_pair_count": pair_count,
        "noise_variance_kind": variance_kind,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "objective_definition": (
            "-1.5*logdet(B_internal) + 0.5*sum_unique_pairs("
            "(D_fit-D_obs)^2/noise_variance)"
        ),
        "logged_metric_timing": "post-update accepted covariance iterate",
    }
    if not converged and status == "max_iterations":
        warnings.warn(
            "fit_gaussian_noise_covariance reached max_iterations before "
            "satisfying the first-order optimality tolerance",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, gram, connectivity, info


def fit_gaussian_noise_covariance_fista(
    squared_distances,
    noise_variance,
    *,
    connectivity_l1=0.0,
    initial_connectivity=None,
    max_iterations=1000,
    relative_tolerance=1e-8,
    absolute_tolerance=1e-10,
    initial_step_size=1.0,
    backtracking_factor=0.5,
    max_backtracking_steps=80,
    accelerated=True,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
):
    """Fit Gaussian-noisy distances by covariance-space proximal gradient.

    With ``connectivity_l1=0``, the solver minimizes the same strictly convex
    objective as :func:`fit_gaussian_noise_covariance`,

    ``0.5*sum((D(B)-Dobs)**2/variance) - 1.5*logdet(B)``,

    exactly preserving the original Gaussian solver. A positive
    ``connectivity_l1=lambda`` adds ``lambda*sum(abs(k_ij))`` to the equivalent
    signed-connectivity dual. In covariance space this replaces each quadratic
    residual by the differentiable epsilon-insensitive loss

    ``soft_threshold(D(B)_ij-Dobs_ij, 2*lambda)**2 / (2*variance_ij)``.

    For either coefficient, the solver splits the covariance objective into a
    differentiable data term and the negative Gaussian entropy. Backtracking
    controls a forward step on the data term. The analytic proximal map of
    ``-1.5*logdet`` transforms every eigenvalue ``y`` of that step to

    ``(y + sqrt(y**2 + 6*eta))/2 > 0``.

    Therefore every accepted covariance is strictly positive definite without
    an SPD projection or a final eigenvalue cutoff. By default, monotone FISTA
    acceleration is used: an extrapolated step that raises the full objective is
    restarted from the last accepted covariance. Set ``accelerated=False`` for
    plain proximal gradient. The current implementation requires a complete
    pairwise squared-distance map so the quadratic measurement operator is
    coercive on the internal symmetric-matrix space.

    Returns
    -------
    fitted_squared_distances, gram_matrix, connectivity_matrix, info
        Fitted total 3D squared distances, centered total-distance Gram matrix,
        signed connectivity, and a convergence/backtracking certificate.
    """
    (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        pair_variance,
        variance_kind,
    ) = _validate_gaussian_noise_variance(squared_distances, noise_variance)
    connectivity_l1 = _validate_connectivity_l1(connectivity_l1)

    n = observed.shape[0]
    n_modes = n - 1
    pair_count = n * (n - 1) // 2
    if len(pair_i) != pair_count:
        raise ValueError(
            "fit_gaussian_noise_covariance_fista currently requires every "
            "off-diagonal squared distance to be observed"
        )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("convergence tolerances must be nonnegative")
    if relative_tolerance == 0.0 and absolute_tolerance == 0.0:
        raise ValueError("at least one convergence tolerance must be positive")
    if not np.isfinite(initial_step_size) or initial_step_size <= 0.0:
        raise ValueError("initial_step_size must be positive and finite")
    if not np.isfinite(backtracking_factor) or not (
        0.0 < backtracking_factor < 1.0
    ):
        raise ValueError("backtracking_factor must lie strictly between zero and one")
    if (
        not isinstance(max_backtracking_steps, (int, np.integer))
        or max_backtracking_steps <= 0
    ):
        raise ValueError("max_backtracking_steps must be a positive integer")
    if not isinstance(accelerated, (bool, np.bool_)):
        raise ValueError("accelerated must be boolean")
    if not np.isfinite(initial_gram_floor_relative) or (
        initial_gram_floor_relative <= 0.0
    ):
        raise ValueError("initial_gram_floor_relative must be positive and finite")

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
    target_pairs = observed[pair_i, pair_j]
    inverse_variance = 1.0 / pair_variance
    target_matrix, inverse_variance_matrix = (
        _gaussian_covariance_pair_matrices(
            n, pair_i, pair_j, target_pairs, inverse_variance
        )
    )
    reduced_gram, initialization = _initialize_gaussian_reduced_gram(
        observed,
        pair_mask,
        pair_i,
        pair_j,
        inverse_variance,
        basis,
        initial_connectivity,
        initial_gram_floor_relative,
    )

    def diagnostics(gram, state=None):
        if state is None:
            state = _gaussian_covariance_objective_gradient(
                gram,
                basis,
                target_matrix,
                inverse_variance_matrix,
                pair_i,
                pair_j,
                connectivity_l1,
            )
        objective, gradient, fitted_pairs, inverse_gram, components = state
        if gradient is None:
            raise RuntimeError("FISTA produced a non-positive-definite covariance")
        negative_entropy, data_objective, entropy_gradient, data_gradient = (
            components
        )
        gradient_norm = float(np.linalg.norm(gradient, ord="fro"))
        gradient_scale = max(
            1.0,
            float(np.linalg.norm(entropy_gradient, ord="fro")),
            float(np.linalg.norm(data_gradient, ord="fro")),
        )
        residual = fitted_pairs - target_pairs
        gram_eigenvalues = np.linalg.eigvalsh(gram)
        if float(gram_eigenvalues[0]) <= 0.0:
            raise RuntimeError("FISTA produced a non-positive covariance eigenvalue")
        entropy = float(
            np.sum(np.log(gram_eigenvalues)) - n_modes * np.log(3.0)
        )
        reduced_connectivity = 3.0 * inverse_gram
        connectivity = -(basis @ reduced_connectivity @ basis.T)
        connectivity = a2a(0.5 * (connectivity + connectivity.T))
        connectivity_pairs = connectivity[np.triu_indices(n, k=1)]
        effective_residual = (
            residual
            if connectivity_l1 == 0.0
            else _soft_threshold(residual, 2.0 * connectivity_l1)
        )
        stationarity = (
            effective_residual - 0.5 * pair_variance * connectivity_pairs
        )
        stationarity_scale = max(float(np.linalg.norm(target_pairs)), 1.0)
        grounded_precision = -connectivity[:n_modes, :n_modes]
        precision_sign, precision_logdet = np.linalg.slogdet(
            grounded_precision
        )
        if precision_sign <= 0.0:
            raise RuntimeError("FISTA produced an invalid grounded precision")
        connectivity_l1_norm = float(np.sum(np.abs(connectivity_pairs)))
        connectivity_l1_penalty = connectivity_l1 * connectivity_l1_norm
        connectivity_dual_objective = (
            -1.5 * float(precision_logdet)
            + 0.5 * float(np.dot(connectivity_pairs, target_pairs))
            + 0.125
            * float(np.dot(pair_variance, np.square(connectivity_pairs)))
            + connectivity_l1_penalty
        )
        metric = {
            "objective": float(objective),
            "negative_entropy_objective": float(negative_entropy),
            "data_objective": float(data_objective),
            "loss": float(
                np.sqrt(np.mean(np.square(residual / target_pairs)))
            ),
            "entropy": entropy,
            "weighted_rmse": float(
                np.sqrt(np.mean(np.square(residual) * inverse_variance))
            ),
            "distance_rmse": float(np.sqrt(np.mean(np.square(residual)))),
            "gradient_norm": gradient_norm,
            "relative_gradient_norm": gradient_norm / gradient_scale,
            "gradient_scale": gradient_scale,
            "minimum_internal_gram_eigenvalue": float(gram_eigenvalues[0]),
            "maximum_internal_gram_eigenvalue": float(gram_eigenvalues[-1]),
            "connectivity_offdiagonal_l2": float(
                np.linalg.norm(connectivity_pairs)
            ),
            "connectivity_offdiagonal_l1": connectivity_l1_norm,
            "connectivity_l1_penalty": connectivity_l1_penalty,
            "connectivity_dual_objective": connectivity_dual_objective,
            "stationarity_residual_norm": float(
                np.linalg.norm(stationarity)
            ),
            "relative_stationarity_residual": float(
                np.linalg.norm(stationarity) / stationarity_scale
            ),
            "maximum_absolute_stationarity_residual": float(
                np.max(np.abs(stationarity))
            ),
        }
        return metric, state, connectivity

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
        "relative_step": [],
        "step_size": [],
        "backtracking_steps": [],
        "proximal_gradient_mapping_norm": [],
        "relative_proximal_gradient_mapping": [],
        "momentum_parameter": [],
        "momentum_coefficient": [],
        "restarted": [],
        "objective_decrease": [],
        "minimum_internal_gram_eigenvalue": [],
        "maximum_internal_gram_eigenvalue": [],
        "connectivity_offdiagonal_l2": [],
        "connectivity_offdiagonal_l1": [],
        "connectivity_l1_penalty": [],
        "connectivity_dual_objective": [],
        "stationarity_residual_norm": [],
        "relative_stationarity_residual": [],
        "maximum_absolute_stationarity_residual": [],
    }
    connectivity_at_steps = {}
    current_metric, current_state, current_connectivity = diagnostics(reduced_gram)
    current_objective = current_metric["objective"]
    search_gram = reduced_gram.copy()
    search_is_extrapolated = False
    momentum_parameter = 1.0
    momentum_coefficient = 0.0
    step_size = float(initial_step_size)
    restart_count = 0
    backtracking_reductions = 0
    stationarity_scale = max(float(np.linalg.norm(target_pairs)), 1.0)

    def convergence_reached(metric):
        if connectivity_l1 > 0.0:
            return metric["stationarity_residual_norm"] <= (
                absolute_tolerance + relative_tolerance * stationarity_scale
            )
        return metric["gradient_norm"] <= (
            absolute_tolerance
            + relative_tolerance * metric["gradient_scale"]
        )

    convergence_status = (
        "stationarity_tolerance"
        if connectivity_l1 > 0.0
        else "optimality_tolerance"
    )
    converged = convergence_reached(current_metric)
    status = convergence_status if converged else "max_iterations"
    message = (
        "first-order optimality tolerance reached at initialization"
        if converged
        else "maximum number of proximal-gradient iterations reached"
    )

    for iteration in range(1, max_iterations + 1):
        if converged:
            break
        restarted = False
        accepted = False
        iteration_backtracking = 0

        while not accepted:
            search_data_objective, search_data_gradient, _, _ = (
                _gaussian_covariance_data_objective_gradient(
                    search_gram,
                    basis,
                    target_matrix,
                    inverse_variance_matrix,
                    pair_i,
                    pair_j,
                    connectivity_l1,
                )
            )
            restart_requested = False
            for _ in range(max_backtracking_steps):
                forward_step = search_gram - step_size * search_data_gradient
                candidate, _, _ = _proximal_gaussian_negative_entropy(
                    forward_step, step_size
                )
                candidate_data_objective = (
                    _gaussian_covariance_data_objective_gradient(
                        candidate,
                        basis,
                        target_matrix,
                        inverse_variance_matrix,
                        pair_i,
                        pair_j,
                        connectivity_l1,
                    )[0]
                )
                search_delta = candidate - search_gram
                quadratic_upper_bound = (
                    search_data_objective
                    + float(np.sum(search_data_gradient * search_delta))
                    + 0.5
                    * float(np.sum(np.square(search_delta)))
                    / step_size
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

                candidate_state = _gaussian_covariance_objective_gradient(
                    candidate,
                    basis,
                    target_matrix,
                    inverse_variance_matrix,
                    pair_i,
                    pair_j,
                    connectivity_l1,
                )
                candidate_objective = float(candidate_state[0])
                monotonicity_tolerance = 1e-12 * max(
                    1.0, abs(current_objective), abs(candidate_objective)
                )
                if candidate_objective > (
                    current_objective + monotonicity_tolerance
                ):
                    if accelerated and search_is_extrapolated:
                        restart_requested = True
                        break
                    step_size *= backtracking_factor
                    iteration_backtracking += 1
                    backtracking_reductions += 1
                    continue
                accepted = True
                break

            if accepted:
                break
            if restart_requested:
                restart_count += 1
                restarted = True
                momentum_parameter = 1.0
                momentum_coefficient = 0.0
                search_gram = reduced_gram.copy()
                search_is_extrapolated = False
                continue
            status = "backtracking_failed"
            message = "proximal-gradient backtracking failed"
            break

        if not accepted:
            break

        previous_gram = reduced_gram
        accepted_delta = candidate - previous_gram
        search_delta = candidate - search_gram
        relative_step = float(np.linalg.norm(accepted_delta, ord="fro")) / max(
            float(np.linalg.norm(previous_gram, ord="fro")), 1.0
        )
        proximal_mapping_norm = float(
            np.linalg.norm(search_delta, ord="fro") / step_size
        )
        reduced_gram = candidate
        candidate_metric, current_state, current_connectivity = diagnostics(
            reduced_gram, candidate_state
        )
        objective_decrease = current_objective - candidate_metric["objective"]
        current_objective = candidate_metric["objective"]
        candidate_metric.update(
            {
                "relative_step": relative_step,
                "step_size": step_size,
                "backtracking_steps": iteration_backtracking,
                "proximal_gradient_mapping_norm": proximal_mapping_norm,
                "relative_proximal_gradient_mapping": proximal_mapping_norm
                / max(candidate_metric["gradient_scale"], 1.0),
                "momentum_parameter": momentum_parameter,
                "momentum_coefficient": momentum_coefficient,
                "restarted": float(restarted),
                "objective_decrease": max(objective_decrease, 0.0),
            }
        )

        history["iteration"].append(iteration)
        for key in history:
            if key != "iteration":
                history[key].append(candidate_metric[key])
        if iteration in save_steps_set:
            connectivity_at_steps[iteration] = current_connectivity.copy()
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "covariance_proximal_optimization",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": current_objective,
                    "loss": candidate_metric["loss"],
                    "entropy": candidate_metric["entropy"],
                    "relative_gradient_norm": candidate_metric[
                        "relative_gradient_norm"
                    ],
                    "relative_stationarity_residual": candidate_metric[
                        "relative_stationarity_residual"
                    ],
                    "step_size": step_size,
                    "backtracking_steps": iteration_backtracking,
                    "restarted": restarted,
                    "noisy": True,
                    "use_gpu": False,
                    "method": "FISTA",
                    "general_method": "covariance_proximal_optimization",
                }
            )

        converged = convergence_reached(candidate_metric)
        current_metric = candidate_metric
        if converged:
            status = convergence_status
            message = "first-order optimality tolerance reached"
            break

        if accelerated:
            next_momentum_parameter = 0.5 * (
                1.0 + np.sqrt(1.0 + 4.0 * momentum_parameter**2)
            )
            next_momentum_coefficient = (
                (momentum_parameter - 1.0) / next_momentum_parameter
            )
            search_gram = reduced_gram + next_momentum_coefficient * (
                reduced_gram - previous_gram
            )
            search_gram = 0.5 * (search_gram + search_gram.T)
            search_is_extrapolated = next_momentum_coefficient > 0.0
            momentum_parameter = next_momentum_parameter
            momentum_coefficient = next_momentum_coefficient
        else:
            search_gram = reduced_gram.copy()
            search_is_extrapolated = False
            momentum_parameter = 1.0
            momentum_coefficient = 0.0

    for key, values in history.items():
        if key in {"iteration", "backtracking_steps"}:
            dtype = np.int64
        elif key == "restarted":
            dtype = np.bool_
        else:
            dtype = np.float64
        history[key] = np.asarray(values, dtype=dtype)

    gram = basis @ reduced_gram @ basis.T
    gram = 0.5 * (gram + gram.T)
    fitted_squared_distances = _squared_distances_from_gram(gram)
    final_metric, _, connectivity = diagnostics(reduced_gram, current_state)
    initialization = dict(initialization)
    initialization["proximal_initial_step_size"] = float(initial_step_size)
    iterations = int(history["iteration"][-1]) if history["iteration"].size else 0
    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": iterations,
        "objective": final_metric["objective"],
        "connectivity_dual_objective": final_metric[
            "connectivity_dual_objective"
        ],
        "gradient_norm": final_metric["gradient_norm"],
        "relative_gradient_norm": final_metric["relative_gradient_norm"],
        "stationarity_residual_norm": final_metric[
            "stationarity_residual_norm"
        ],
        "relative_stationarity_residual": final_metric[
            "relative_stationarity_residual"
        ],
        "maximum_absolute_stationarity_residual": final_metric[
            "maximum_absolute_stationarity_residual"
        ],
        "connectivity_l1": float(connectivity_l1),
        "connectivity_l1_penalty": final_metric[
            "connectivity_l1_penalty"
        ],
        "relative_tolerance": float(relative_tolerance),
        "absolute_tolerance": float(absolute_tolerance),
        "accelerated": bool(accelerated),
        "monotone_restart": bool(accelerated),
        "restart_count": int(restart_count),
        "backtracking_reductions": int(backtracking_reductions),
        "initial_step_size": float(initial_step_size),
        "final_step_size": float(step_size),
        "backtracking_factor": float(backtracking_factor),
        "observed_pair_count": pair_count,
        "noise_variance_kind": variance_kind,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "objective_definition": (
            "-1.5*logdet(B_internal) + 0.5*sum_unique_pairs("
            "soft_threshold(D_fit-D_obs, 2*connectivity_l1)^2/"
            "noise_variance)"
            if connectivity_l1 > 0.0
            else "-1.5*logdet(B_internal) + 0.5*sum_unique_pairs("
            "(D_fit-D_obs)^2/noise_variance)"
        ),
        "connectivity_dual_objective_definition": (
            "-1.5*logdet(P_grounded) + 0.5*sum_unique_pairs("
            "k_ij*D_obs_ij) + 0.125*sum_unique_pairs("
            "noise_variance_ij*k_ij**2) + connectivity_l1*"
            "sum_unique_pairs(abs(k_ij))"
        ),
        "stationarity_definition": (
            "soft_threshold(D_fit_ij-D_obs_ij, 2*connectivity_l1)-"
            "noise_variance_ij*k_ij/2"
        ),
        "proximal_eigenvalue_update": (
            "x=(y+sqrt(y^2+6*step_size))/2, evaluated with a "
            "cancellation-safe negative-y branch"
        ),
        "convergence_definition": (
            "physical L1 KKT stationarity residual norm"
            if connectivity_l1 > 0.0
            else "Frobenius norm of the full covariance-space objective gradient"
        ),
        "allows_signed_offdiagonal_connectivity": True,
        "requires_complete_pairwise_observations": True,
        "logged_metric_timing": "post-update accepted proximal iterate",
    }
    if not converged and status == "max_iterations":
        warnings.warn(
            "fit_gaussian_noise_covariance_fista reached max_iterations before "
            "satisfying the first-order optimality tolerance",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, gram, connectivity, info


def _connectivity_from_grounded_precision(grounded_precision):
    """Reconstruct the full row-sum-zero connectivity from one grounded minor."""
    grounded_precision = np.asarray(grounded_precision, dtype=np.float64)
    n_modes = grounded_precision.shape[0]
    connectivity = np.zeros((n_modes + 1, n_modes + 1), dtype=np.float64)
    row_sums = np.sum(grounded_precision, axis=1)
    connectivity[:n_modes, :n_modes] = -grounded_precision
    connectivity[:n_modes, n_modes] = row_sums
    connectivity[n_modes, :n_modes] = row_sums
    connectivity[n_modes, n_modes] = -float(np.sum(row_sums))
    return a2a(0.5 * (connectivity + connectivity.T))


def _squared_distances_from_grounded_covariance(covariance):
    """Return total 3D squared distances from a grounded one-coordinate covariance."""
    covariance = np.asarray(covariance, dtype=np.float64)
    n_modes = covariance.shape[0]
    covariance = 0.5 * (covariance + covariance.T)
    diagonal = np.diag(covariance)
    squared_distances = np.zeros((n_modes + 1, n_modes + 1), dtype=np.float64)
    squared_distances[:n_modes, :n_modes] = 3.0 * (
        diagonal[:, np.newaxis] + diagonal - 2.0 * covariance
    )
    squared_distances[:n_modes, n_modes] = 3.0 * diagonal
    squared_distances[n_modes, :n_modes] = 3.0 * diagonal
    squared_distances = 0.5 * (squared_distances + squared_distances.T)
    np.fill_diagonal(squared_distances, 0.0)
    return squared_distances


def _squared_distances_from_grounded_precision(grounded_precision):
    """Return total 3D squared distances from a grounded one-coordinate precision."""
    grounded_precision = np.asarray(grounded_precision, dtype=np.float64)
    n_modes = grounded_precision.shape[0]
    covariance = np.linalg.solve(
        grounded_precision, np.eye(n_modes, dtype=np.float64)
    )
    return _squared_distances_from_grounded_covariance(covariance)


def _centered_gram_from_squared_distances(squared_distances):
    """Return the centered total-distance Gram matrix induced by an EDM."""
    squared_distances = np.asarray(squared_distances, dtype=np.float64)
    row_mean = np.mean(squared_distances, axis=1, keepdims=True)
    gram = -0.5 * (
        squared_distances
        - row_mean
        - row_mean.T
        + np.mean(squared_distances)
    )
    return 0.5 * (gram + gram.T)


def _gaussian_connectivity_physical_diagnostics(
    grounded_precision,
    observed,
    variance_matrix,
    objective_components,
    *,
    parameter_gradient=None,
    parameter_scale=1.0,
    fitted_squared_distances=None,
    connectivity_l1=0.0,
):
    """Evaluate shared physical diagnostics for exact connectivity solvers."""
    grounded_precision = np.asarray(grounded_precision, dtype=np.float64)
    observed = np.asarray(observed, dtype=np.float64)
    variance_matrix = np.asarray(variance_matrix, dtype=np.float64)
    n = observed.shape[0]
    upper = np.triu_indices(n, k=1)
    observed_pairs = observed[upper]
    variance_pairs = variance_matrix[upper]
    connectivity = _connectivity_from_grounded_precision(grounded_precision)
    if fitted_squared_distances is None:
        fitted = _squared_distances_from_grounded_precision(grounded_precision)
    else:
        fitted = np.asarray(fitted_squared_distances, dtype=np.float64)
    residual = fitted[upper] - observed_pairs
    connectivity_pairs = connectivity[upper]
    effective_residual = (
        residual
        if connectivity_l1 == 0.0
        else _soft_threshold(residual, 2.0 * connectivity_l1)
    )
    stationarity = (
        effective_residual - 0.5 * variance_pairs * connectivity_pairs
    )

    if parameter_gradient is None:
        parameter_gradient = -0.5 * stationarity
    parameter_gradient = np.asarray(parameter_gradient, dtype=np.float64)
    gradient_norm = float(np.linalg.norm(parameter_gradient))
    gradient_infinity_norm = float(np.max(np.abs(parameter_gradient)))
    precision_eigenvalues = np.linalg.eigvalsh(grounded_precision)
    sign, precision_logdet = np.linalg.slogdet(grounded_precision)
    if sign <= 0.0:
        raise RuntimeError("exact connectivity solver produced an invalid precision")

    negative_entropy = float(
        objective_components["negative_entropy_objective"]
    )
    distance_linear = float(objective_components["distance_linear_objective"])
    gaussian_penalty = float(
        objective_components["gaussian_connectivity_penalty"]
    )
    l1_penalty = float(
        objective_components.get("connectivity_l1_penalty", 0.0)
    )
    objective = (
        negative_entropy + distance_linear + gaussian_penalty + l1_penalty
    )
    metric = {
        "objective": objective,
        "negative_entropy_objective": negative_entropy,
        "distance_linear_objective": distance_linear,
        "gaussian_connectivity_penalty": gaussian_penalty,
        "connectivity_l1_penalty": l1_penalty,
        "data_objective": distance_linear + gaussian_penalty + l1_penalty,
        "loss": float(np.sqrt(np.mean(np.square(residual / observed_pairs)))),
        "entropy": float(-precision_logdet - np.log(n)),
        "weighted_rmse": float(
            np.sqrt(np.mean(np.square(residual) / variance_pairs))
        ),
        "distance_rmse": float(np.sqrt(np.mean(np.square(residual)))),
        "gradient_norm": gradient_norm,
        "gradient_infinity_norm": gradient_infinity_norm,
        "relative_gradient_norm": gradient_norm
        / max(float(parameter_scale), 1.0),
        "stationarity_residual_norm": float(np.linalg.norm(stationarity)),
        "relative_stationarity_residual": float(
            np.linalg.norm(stationarity)
            / max(float(np.linalg.norm(observed_pairs)), 1.0)
        ),
        "maximum_absolute_stationarity_residual": float(
            np.max(np.abs(stationarity))
        ),
        "minimum_grounded_precision_eigenvalue": float(
            precision_eigenvalues[0]
        ),
        "maximum_grounded_precision_eigenvalue": float(
            precision_eigenvalues[-1]
        ),
        "connectivity_offdiagonal_l2": float(np.linalg.norm(connectivity_pairs)),
        "connectivity_offdiagonal_l1": float(
            np.sum(np.abs(connectivity_pairs))
        ),
    }
    return metric, fitted, connectivity, stationarity


def _pack_log_cholesky(cholesky_factor, lower_indices, diagonal_mask):
    """Pack a lower Cholesky factor, using logs for its positive diagonal."""
    parameters = np.asarray(cholesky_factor[lower_indices], dtype=np.float64).copy()
    parameters[diagonal_mask] = np.log(parameters[diagonal_mask])
    return parameters


def _unpack_log_cholesky(parameters, dimension, lower_indices, diagonal_mask):
    """Unpack unconstrained parameters into a lower factor with positive diagonal."""
    values = np.asarray(parameters, dtype=np.float64).copy()
    values[diagonal_mask] = np.exp(values[diagonal_mask])
    cholesky_factor = np.zeros((dimension, dimension), dtype=np.float64)
    cholesky_factor[lower_indices] = values
    return cholesky_factor


def _gaussian_connectivity_cholesky_objective_gradient(
    parameters,
    lower_indices,
    diagonal_mask,
    inverse_sqrt_reference_covariance,
    whitened_observed_covariance,
    internal_variance,
    pinned_variance,
):
    """Evaluate the calibrated Gaussian dual objective and its exact gradient."""
    inverse_sqrt_reference_covariance = np.asarray(
        inverse_sqrt_reference_covariance, dtype=np.float64
    )
    dimension = inverse_sqrt_reference_covariance.shape[0]
    cholesky_factor = _unpack_log_cholesky(
        parameters, dimension, lower_indices, diagonal_mask
    )
    whitened_precision = cholesky_factor @ cholesky_factor.T
    grounded_precision = (
        inverse_sqrt_reference_covariance
        @ whitened_precision
        @ inverse_sqrt_reference_covariance
    )
    grounded_precision = 0.5 * (grounded_precision + grounded_precision.T)

    diagonal = np.diag(cholesky_factor)
    negative_entropy = -3.0 * float(np.sum(np.log(diagonal)))
    distance_linear = 1.5 * float(
        np.sum(whitened_precision * whitened_observed_covariance)
    )

    internal_upper = np.triu_indices(dimension, k=1)
    internal_values = grounded_precision[internal_upper]
    internal_penalty = 0.125 * float(
        np.dot(internal_variance[internal_upper], np.square(internal_values))
    )
    ones = np.ones(dimension, dtype=np.float64)
    pinned_connectivity = grounded_precision @ ones
    pinned_penalty = 0.125 * float(
        np.dot(pinned_variance, np.square(pinned_connectivity))
    )
    gaussian_penalty = internal_penalty + pinned_penalty
    objective = negative_entropy + distance_linear + gaussian_penalty

    # The full-matrix gradient is defined with symmetric perturbations. For an
    # internal edge, k_ij=-P_ij; for an edge to the pinned locus, k_iN=(P1)_i.
    penalty_gradient = np.zeros_like(grounded_precision)
    internal_gradient = (
        0.125 * internal_variance[internal_upper] * internal_values
    )
    penalty_gradient[internal_upper] = internal_gradient
    penalty_gradient[(internal_upper[1], internal_upper[0])] = internal_gradient
    weighted_pinned = pinned_variance * pinned_connectivity
    penalty_gradient += 0.125 * (
        np.outer(weighted_pinned, ones) + np.outer(ones, weighted_pinned)
    )
    whitened_penalty_gradient = (
        inverse_sqrt_reference_covariance
        @ penalty_gradient
        @ inverse_sqrt_reference_covariance
    )

    factor_gradient = (
        3.0 * whitened_observed_covariance @ cholesky_factor
        + 2.0 * whitened_penalty_gradient @ cholesky_factor
    )
    factor_gradient[np.diag_indices(dimension)] -= 3.0 / diagonal
    gradient = np.asarray(factor_gradient[lower_indices], dtype=np.float64)
    gradient[diagonal_mask] *= diagonal

    state = {
        "cholesky_factor": cholesky_factor,
        "whitened_precision": whitened_precision,
        "grounded_precision": grounded_precision,
        "negative_entropy_objective": negative_entropy,
        "distance_linear_objective": distance_linear,
        "gaussian_connectivity_penalty": gaussian_penalty,
    }
    return float(objective), gradient, state


def fit_gaussian_noise_connectivity_cholesky(
    squared_distances,
    noise_variance,
    *,
    initial_connectivity=None,
    max_iterations=500,
    gradient_tolerance=1e-8,
    function_tolerance=1e-12,
    stationarity_tolerance=1e-7,
    history_size=20,
    initial_gram_floor_relative=1e-8,
    save_steps=None,
    progress_callback=None,
):
    """Fit Gaussian-noisy distances through a signed connectivity Cholesky model.

    One locus is pinned and the resulting grounded precision is parameterized as
    ``P = C0**(-1/2) L L.T C0**(-1/2)``, with an exponential diagonal for
    ``L``. This spans every stable row-sum-zero connectivity matrix and places no
    sign constraint on its off-diagonal couplings. The standard (not doubled)
    dual objective is

    ``-3/2 logdet(P) + 1/2 sum k_ij Dobs_ij + 1/8 sum v_ij k_ij**2``.

    The evaluated objective omits only the fixed log-determinant constant from
    whitening by ``C0``; its minimizer and derivatives are unchanged.

    The ``1/8`` coefficient is the connectivity-space form of the calibrated
    Gaussian error model. Its stationarity condition is
    ``Dfit_ij-Dobs_ij = v_ij*k_ij/2`` under the HIPPS-DIMES Hamiltonian
    convention. Optimization uses analytic gradients and L-BFGS-B. Unlike
    :func:`fit_gaussian_noise_covariance`, this implementation currently
    requires a complete squared-distance matrix.

    Returns
    -------
    fitted_squared_distances, gram_matrix, connectivity_matrix, info
        The fitted total 3D squared distances, centered total-distance Gram
        matrix, signed connectivity, and convergence/provenance diagnostics.
    """
    (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        pair_variance,
        variance_kind,
    ) = _validate_gaussian_noise_variance(squared_distances, noise_variance)

    n = observed.shape[0]
    n_modes = n - 1
    expected_pair_count = n * (n - 1) // 2
    if len(pair_i) != expected_pair_count:
        raise ValueError(
            "fit_gaussian_noise_connectivity_cholesky currently requires "
            "every off-diagonal squared distance to be observed"
        )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    for value, name in (
        (gradient_tolerance, "gradient_tolerance"),
        (function_tolerance, "function_tolerance"),
        (stationarity_tolerance, "stationarity_tolerance"),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not isinstance(history_size, (int, np.integer)) or history_size <= 0:
        raise ValueError("history_size must be a positive integer")
    if not np.isfinite(initial_gram_floor_relative) or (
        initial_gram_floor_relative <= 0.0
    ):
        raise ValueError("initial_gram_floor_relative must be positive and finite")

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
    initial_reduced_gram, initialization = _initialize_gaussian_reduced_gram(
        observed,
        pair_mask,
        pair_i,
        pair_j,
        inverse_variance,
        basis,
        initial_connectivity,
        initial_gram_floor_relative,
    )
    initial_full_connectivity = _connectivity_from_reduced_gram(
        initial_reduced_gram, basis
    )
    initial_grounded_precision = -initial_full_connectivity[:n_modes, :n_modes]
    reference_covariance = np.linalg.solve(
        initial_grounded_precision, np.eye(n_modes, dtype=np.float64)
    )
    reference_covariance = 0.5 * (reference_covariance + reference_covariance.T)
    reference_eigenvalues, reference_eigenvectors = np.linalg.eigh(
        reference_covariance
    )
    if float(reference_eigenvalues[0]) <= 0.0:
        raise RuntimeError("initial grounded covariance is not positive definite")
    inverse_sqrt_reference_covariance = (
        reference_eigenvectors
        * (1.0 / np.sqrt(reference_eigenvalues))
    ) @ reference_eigenvectors.T

    # Total 3D squared distances become one-coordinate variances after division
    # by three; polarization contributes the additional factor of one half.
    # This /6 conversion is where the familiar variance /9 scaling enters.
    pinned_distances = observed[:n_modes, n_modes]
    observed_covariance = (
        pinned_distances[:, np.newaxis]
        + pinned_distances
        - observed[:n_modes, :n_modes]
    ) / 6.0
    observed_covariance = 0.5 * (observed_covariance + observed_covariance.T)
    whitened_observed_covariance = (
        inverse_sqrt_reference_covariance
        @ observed_covariance
        @ inverse_sqrt_reference_covariance
    )
    whitened_observed_covariance = 0.5 * (
        whitened_observed_covariance + whitened_observed_covariance.T
    )

    variance_matrix = np.zeros_like(observed)
    variance_matrix[pair_i, pair_j] = pair_variance
    variance_matrix[pair_j, pair_i] = pair_variance
    internal_variance = variance_matrix[:n_modes, :n_modes]
    pinned_variance = variance_matrix[:n_modes, n_modes]

    lower_indices = np.tril_indices(n_modes)
    diagonal_mask = lower_indices[0] == lower_indices[1]
    initial_factor = np.eye(n_modes, dtype=np.float64)
    initial_parameters = _pack_log_cholesky(
        initial_factor, lower_indices, diagonal_mask
    )

    history = {
        "iteration": [],
        "objective": [],
        "negative_entropy_objective": [],
        "distance_linear_objective": [],
        "gaussian_connectivity_penalty": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "gradient_norm": [],
        "gradient_infinity_norm": [],
        "relative_gradient_norm": [],
        "stationarity_residual_norm": [],
        "relative_stationarity_residual": [],
        "maximum_absolute_stationarity_residual": [],
        "minimum_grounded_precision_eigenvalue": [],
        "maximum_grounded_precision_eigenvalue": [],
        "connectivity_offdiagonal_l2": [],
    }
    connectivity_at_steps = {}

    objective_arguments = (
        lower_indices,
        diagonal_mask,
        inverse_sqrt_reference_covariance,
        whitened_observed_covariance,
        internal_variance,
        pinned_variance,
    )

    def evaluate(parameters):
        return _gaussian_connectivity_cholesky_objective_gradient(
            parameters, *objective_arguments
        )

    def diagnostics(parameters):
        objective, gradient, state = evaluate(parameters)
        objective_components = {
            "negative_entropy_objective": state["negative_entropy_objective"],
            "distance_linear_objective": state["distance_linear_objective"],
            "gaussian_connectivity_penalty": state[
                "gaussian_connectivity_penalty"
            ],
        }
        metric, fitted, connectivity, stationarity = (
            _gaussian_connectivity_physical_diagnostics(
                state["grounded_precision"],
                observed,
                variance_matrix,
                objective_components,
                parameter_gradient=gradient,
                parameter_scale=np.linalg.norm(parameters),
            )
        )
        # The shared helper reconstructs the same objective from its components.
        if not np.isclose(metric["objective"], objective, rtol=1e-12, atol=1e-12):
            raise RuntimeError("inconsistent Cholesky objective diagnostics")
        return metric, fitted, connectivity, stationarity

    def callback(parameters):
        metric, _, connectivity, _ = diagnostics(parameters)
        iteration = len(history["iteration"]) + 1
        history["iteration"].append(iteration)
        for key in history:
            if key != "iteration":
                history[key].append(metric[key])
        if iteration in save_steps_set:
            connectivity_at_steps[iteration] = connectivity.copy()
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "connectivity_cholesky_optimization",
                    "iteration": iteration,
                    "total": max_iterations,
                    "objective": metric["objective"],
                    "loss": metric["loss"],
                    "entropy": metric["entropy"],
                    "relative_gradient_norm": metric["relative_gradient_norm"],
                    "relative_stationarity_residual": metric[
                        "relative_stationarity_residual"
                    ],
                    "noisy": True,
                    "use_gpu": False,
                    "method": "CHOL",
                    "general_method": "connectivity_cholesky_optimization",
                }
            )

    max_line_search_steps = 50
    max_function_evaluations = (
        int(max_iterations) * (max_line_search_steps + 1) + 1
    )
    optimization_result = scipy.optimize.minimize(
        evaluate,
        initial_parameters,
        method="L-BFGS-B",
        jac=True,
        callback=callback,
        options={
            "maxiter": int(max_iterations),
            "maxcor": int(history_size),
            "gtol": float(gradient_tolerance),
            "ftol": float(function_tolerance),
            "maxls": max_line_search_steps,
            # SciPy's default maxfun=15000 can terminate a requested long run
            # before maxiter. The line-search-derived ceiling leaves maxiter as
            # the effective user budget while retaining a finite safeguard.
            "maxfun": max_function_evaluations,
        },
    )

    final_metric, fitted_squared_distances, connectivity, _ = diagnostics(
        optimization_result.x
    )
    gram = _centered_gram_from_squared_distances(fitted_squared_distances)

    for key, values in history.items():
        dtype = np.int64 if key == "iteration" else np.float64
        history[key] = np.asarray(values, dtype=dtype)

    converged = (
        final_metric["relative_stationarity_residual"] <= stationarity_tolerance
    )
    if converged:
        status = "stationarity_tolerance"
        message = "physical first-order stationarity tolerance reached"
    elif optimization_result.success:
        status = "optimizer_stopped_without_stationarity"
        message = (
            "L-BFGS-B stopped successfully, but the physical stationarity "
            "residual remains above tolerance"
        )
    elif int(optimization_result.nit) >= max_iterations:
        status = "max_iterations"
        message = "maximum number of L-BFGS-B iterations reached"
    else:
        status = "optimizer_failed"
        message = str(optimization_result.message)

    initialization = dict(initialization)
    initialization.update(
        {
            "pinned_locus": n - 1,
            "whitening_minimum_covariance_eigenvalue": float(
                reference_eigenvalues[0]
            ),
            "whitening_maximum_covariance_eigenvalue": float(
                reference_eigenvalues[-1]
            ),
        }
    )
    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": int(optimization_result.nit),
        "optimizer_success": bool(optimization_result.success),
        "optimizer_status": int(optimization_result.status),
        "optimizer_message": str(optimization_result.message),
        "optimizer_function_evaluations": int(optimization_result.nfev),
        "optimizer_gradient_evaluations": int(optimization_result.njev),
        "maximum_function_evaluations": max_function_evaluations,
        "objective": final_metric["objective"],
        "gradient_norm": final_metric["gradient_norm"],
        "gradient_infinity_norm": final_metric["gradient_infinity_norm"],
        "relative_gradient_norm": final_metric["relative_gradient_norm"],
        "stationarity_residual_norm": final_metric[
            "stationarity_residual_norm"
        ],
        "relative_stationarity_residual": final_metric[
            "relative_stationarity_residual"
        ],
        "maximum_absolute_stationarity_residual": final_metric[
            "maximum_absolute_stationarity_residual"
        ],
        "stationarity_tolerance": float(stationarity_tolerance),
        "observed_pair_count": expected_pair_count,
        "noise_variance_kind": variance_kind,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "objective_definition": (
            "-1.5*logdet(P_grounded) + 0.5*sum_unique_pairs("
            "k_ij*D_obs_ij) + 0.125*sum_unique_pairs("
            "noise_variance_ij*k_ij**2), up to the fixed whitening constant"
        ),
        "stationarity_definition": (
            "D_fit_ij-D_obs_ij-noise_variance_ij*k_ij/2"
        ),
        "parameterization": (
            "P=C0^(-1/2)*L*L.T*C0^(-1/2), with exp(log_diagonal(L))"
        ),
        "allows_signed_offdiagonal_connectivity": True,
        "logged_metric_timing": "post-update accepted L-BFGS-B iterate",
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_connectivity_cholesky stopped before satisfying "
            "the physical stationarity tolerance",
            RuntimeWarning,
            stacklevel=2,
        )
    return fitted_squared_distances, gram, connectivity, info


def _gaussian_connectivity_coordinate_minimum(
    current_connectivity,
    current_squared_distance,
    observed_squared_distance,
    noise_variance,
    connectivity_l1=0.0,
):
    """Return the exact stable one-edge Gaussian/L1 dual minimizer."""
    rho = current_squared_distance / 3.0
    if not np.isfinite(rho) or rho <= 0.0:
        raise ValueError("current_squared_distance must be positive and finite")
    if not np.isfinite(noise_variance) or noise_variance <= 0.0:
        raise ValueError("noise_variance must be positive and finite")

    def smooth_candidate(effective_observed):
        linear_coefficient = (
            0.5 * effective_observed
            + 0.25 * noise_variance * current_connectivity
            - 0.25 * noise_variance / rho
        )
        discriminant = np.hypot(
            linear_coefficient, np.sqrt(1.5 * noise_variance)
        )
        if linear_coefficient >= 0.0:
            determinant_ratio = 3.0 * rho / (
                linear_coefficient + discriminant
            )
        else:
            determinant_ratio = (
                (-linear_coefficient + discriminant)
                / (noise_variance / (2.0 * rho))
            )
        coordinate_delta = (determinant_ratio - 1.0) / rho
        updated_connectivity = current_connectivity + coordinate_delta
        return updated_connectivity, coordinate_delta, determinant_ratio

    def full_objective_change(updated_connectivity, determinant_ratio):
        coordinate_delta = updated_connectivity - current_connectivity
        return (
            -1.5 * np.log(determinant_ratio)
            + 0.5 * observed_squared_distance * coordinate_delta
            + 0.125
            * noise_variance
            * (
                updated_connectivity * updated_connectivity
                - current_connectivity * current_connectivity
            )
            + connectivity_l1
            * (abs(updated_connectivity) - abs(current_connectivity))
        )

    if connectivity_l1 == 0.0:
        updated_connectivity, coordinate_delta, determinant_ratio = (
            smooth_candidate(observed_squared_distance)
        )
        objective_change = full_objective_change(
            updated_connectivity, determinant_ratio
        )
        return (
            float(coordinate_delta),
            float(determinant_ratio),
            float(objective_change),
        )

    candidates = []
    for sign, effective_observed in (
        (1.0, observed_squared_distance + 2.0 * connectivity_l1),
        (-1.0, observed_squared_distance - 2.0 * connectivity_l1),
    ):
        updated_connectivity, _, determinant_ratio = smooth_candidate(
            effective_observed
        )
        if sign * updated_connectivity > 0.0 and determinant_ratio > 0.0:
            candidates.append(
                (
                    full_objective_change(
                        updated_connectivity, determinant_ratio
                    ),
                    updated_connectivity,
                    determinant_ratio,
                )
            )

    zero_determinant_ratio = 1.0 - rho * current_connectivity
    if zero_determinant_ratio > 0.0:
        candidates.append(
            (
                full_objective_change(0.0, zero_determinant_ratio),
                0.0,
                zero_determinant_ratio,
            )
        )
    if not candidates:
        raise RuntimeError("L1 coordinate minimization found no stable candidate")

    objective_change, updated_connectivity, determinant_ratio = min(
        candidates, key=lambda candidate: candidate[0]
    )
    coordinate_delta = updated_connectivity - current_connectivity
    return (
        float(coordinate_delta),
        float(determinant_ratio),
        float(objective_change),
    )


def fit_gaussian_noise_connectivity_coordinate_descent(
    squared_distances,
    noise_variance,
    *,
    connectivity_l1=0.0,
    initial_connectivity=None,
    max_iterations=500,
    stationarity_tolerance=1e-7,
    initial_gram_floor_relative=1e-2,
    save_steps=None,
    progress_callback=None,
):
    """Fit the calibrated Gaussian model by exact cyclic coordinate minimization.

    One iteration is one deterministic sweep through all unique pair couplings.
    For edge ``e`` with incidence vector ``b_e``, the grounded precision update
    is ``P_new=P+delta*b_e*b_e.T``. The matrix determinant lemma and
    Sherman-Morrison identity reduce the exact dual line minimization to a
    closed-form scalar update. Every coordinate step obeys
    ``1+delta*b_e.T@inv(P)@b_e > 0``, so the grounded precision stays strictly
    positive definite while off-diagonal connectivity values remain signed.

    The minimized standard dual objective is

    ``-3/2 logdet(P) + 1/2 sum k_ij Dobs_ij + 1/8 sum v_ij k_ij**2``

    ``+ connectivity_l1 * sum abs(k_ij)``.

    Positive variance for every pair makes this objective strictly convex in
    the independent connectivity coordinates. Each L1 coordinate update checks
    the negative branch, the exact-zero kink, and the positive branch. Exact
    cyclic minimization therefore targets the same unique solution as
    :func:`fit_gaussian_noise_covariance_fista` with the same L1 coefficient.
    The nearest-EDM start uses a moderate spectral floor because an invalid EDM
    otherwise projects onto the covariance-cone boundary and can give a severely
    ill-conditioned initial precision. This implementation currently requires a
    complete squared-distance matrix.

    Returns
    -------
    fitted_squared_distances, gram_matrix, connectivity_matrix, info
        Fitted total 3D squared distances, centered total-distance Gram matrix,
        signed connectivity, and sweep-level convergence diagnostics.
    """
    (
        observed,
        pair_mask,
        pair_i,
        pair_j,
        pair_variance,
        variance_kind,
    ) = _validate_gaussian_noise_variance(squared_distances, noise_variance)
    connectivity_l1 = _validate_connectivity_l1(connectivity_l1)

    n = observed.shape[0]
    n_modes = n - 1
    pair_count = n * (n - 1) // 2
    if len(pair_i) != pair_count:
        raise ValueError(
            "fit_gaussian_noise_connectivity_coordinate_descent currently "
            "requires every off-diagonal squared distance to be observed"
        )
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if not np.isfinite(stationarity_tolerance) or stationarity_tolerance <= 0.0:
        raise ValueError("stationarity_tolerance must be positive and finite")
    if not np.isfinite(initial_gram_floor_relative) or (
        initial_gram_floor_relative <= 0.0
    ):
        raise ValueError("initial_gram_floor_relative must be positive and finite")

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
    initial_reduced_gram, initialization = _initialize_gaussian_reduced_gram(
        observed,
        pair_mask,
        pair_i,
        pair_j,
        inverse_variance,
        basis,
        initial_connectivity,
        initial_gram_floor_relative,
    )
    connectivity = _connectivity_from_reduced_gram(initial_reduced_gram, basis)
    upper = (pair_i, pair_j)
    connectivity_values = np.asarray(connectivity[upper], dtype=np.float64).copy()
    observed_pairs = observed[upper]
    variance_matrix = np.zeros_like(observed)
    variance_matrix[pair_i, pair_j] = pair_variance
    variance_matrix[pair_j, pair_i] = pair_variance

    grounded_precision = -connectivity[:n_modes, :n_modes]
    grounded_precision = 0.5 * (grounded_precision + grounded_precision.T)
    try:
        precision_factor = scipy.linalg.cho_factor(
            grounded_precision, lower=True, check_finite=False
        )
    except np.linalg.LinAlgError as error:
        raise RuntimeError(
            "initial grounded precision is not strictly positive definite"
        ) from error
    grounded_covariance = scipy.linalg.cho_solve(
        precision_factor,
        np.eye(n_modes, dtype=np.float64),
        check_finite=False,
    )
    grounded_covariance = 0.5 * (
        grounded_covariance + grounded_covariance.T
    )

    history = {
        "iteration": [],
        "objective": [],
        "negative_entropy_objective": [],
        "distance_linear_objective": [],
        "gaussian_connectivity_penalty": [],
        "connectivity_l1_penalty": [],
        "data_objective": [],
        "loss": [],
        "entropy": [],
        "weighted_rmse": [],
        "distance_rmse": [],
        "gradient_norm": [],
        "gradient_infinity_norm": [],
        "relative_gradient_norm": [],
        "stationarity_residual_norm": [],
        "relative_stationarity_residual": [],
        "maximum_absolute_stationarity_residual": [],
        "minimum_grounded_precision_eigenvalue": [],
        "maximum_grounded_precision_eigenvalue": [],
        "connectivity_offdiagonal_l2": [],
        "connectivity_offdiagonal_l1": [],
        "objective_decrease": [],
        "relative_objective_decrease": [],
        "maximum_absolute_coordinate_update": [],
        "minimum_coordinate_determinant_ratio": [],
        "sweep_relaxation": [],
    }
    connectivity_at_steps = {}

    def diagnostics():
        sign, precision_logdet = np.linalg.slogdet(grounded_precision)
        if sign <= 0.0:
            raise RuntimeError("coordinate descent produced an invalid precision")
        objective_components = {
            "negative_entropy_objective": -1.5 * float(precision_logdet),
            "distance_linear_objective": 0.5
            * float(np.dot(connectivity_values, observed_pairs)),
            "gaussian_connectivity_penalty": 0.125
            * float(
                np.dot(pair_variance, np.square(connectivity_values))
            ),
            "connectivity_l1_penalty": connectivity_l1
            * float(np.sum(np.abs(connectivity_values))),
        }
        fitted = _squared_distances_from_grounded_covariance(
            grounded_covariance
        )
        return _gaussian_connectivity_physical_diagnostics(
            grounded_precision,
            observed,
            variance_matrix,
            objective_components,
            parameter_scale=np.linalg.norm(connectivity_values),
            fitted_squared_distances=fitted,
            connectivity_l1=connectivity_l1,
        )

    previous_metric, fitted_squared_distances, connectivity, _ = diagnostics()
    converged = (
        previous_metric["relative_stationarity_residual"]
        <= stationarity_tolerance
    )
    status = "stationarity_tolerance" if converged else "max_iterations"
    message = (
        "physical first-order stationarity tolerance reached at initialization"
        if converged
        else "maximum number of coordinate sweeps reached"
    )
    coordinate_updates = 0
    sweep_retries = 0
    damped_sweeps = 0
    minimum_accepted_sweep_relaxation = 1.0

    for sweep in range(1, max_iterations + 1):
        if converged:
            break
        sweep_start_connectivity_values = connectivity_values.copy()
        sweep_start_precision = grounded_precision.copy()
        sweep_start_covariance = grounded_covariance.copy()
        sweep_relaxation = 1.0

        while True:
            connectivity_values = sweep_start_connectivity_values.copy()
            grounded_precision = sweep_start_precision.copy()
            grounded_covariance = sweep_start_covariance.copy()
            maximum_absolute_update = 0.0
            minimum_determinant_ratio = np.inf
            attempt_failed = False

            for edge_index, (i, j) in enumerate(zip(pair_i, pair_j)):
                if j == n_modes:
                    covariance_times_incidence = (
                        grounded_covariance[:, i].copy()
                    )
                    rho = float(grounded_covariance[i, i])
                else:
                    covariance_times_incidence = (
                        grounded_covariance[:, i]
                        - grounded_covariance[:, j]
                    )
                    rho = float(
                        grounded_covariance[i, i]
                        + grounded_covariance[j, j]
                        - 2.0 * grounded_covariance[i, j]
                    )
                if not np.isfinite(rho) or rho <= 0.0:
                    attempt_failed = True
                    break

                current_squared_distance = 3.0 * rho
                exact_delta, _, _ = _gaussian_connectivity_coordinate_minimum(
                    connectivity_values[edge_index],
                    current_squared_distance,
                    observed_pairs[edge_index],
                    pair_variance[edge_index],
                    connectivity_l1,
                )
                coordinate_delta = sweep_relaxation * exact_delta
                determinant_ratio = 1.0 + rho * coordinate_delta
                if (
                    not np.isfinite(determinant_ratio)
                    or determinant_ratio <= 0.0
                ):
                    attempt_failed = True
                    break

                connectivity_values[edge_index] += coordinate_delta
                if j == n_modes:
                    grounded_precision[i, i] += coordinate_delta
                else:
                    grounded_precision[i, i] += coordinate_delta
                    grounded_precision[j, j] += coordinate_delta
                    grounded_precision[i, j] -= coordinate_delta
                    grounded_precision[j, i] -= coordinate_delta
                grounded_covariance -= (
                    coordinate_delta / determinant_ratio
                ) * np.outer(
                    covariance_times_incidence, covariance_times_incidence
                )
                grounded_covariance = 0.5 * (
                    grounded_covariance + grounded_covariance.T
                )
                maximum_absolute_update = max(
                    maximum_absolute_update, abs(coordinate_delta)
                )
                minimum_determinant_ratio = min(
                    minimum_determinant_ratio, determinant_ratio
                )

            if not attempt_failed:
                # Reconstruct and refactor once per sweep to remove accumulated
                # Sherman-Morrison roundoff while preserving the edge values.
                connectivity = _connectivity_from_grounded_precision(
                    grounded_precision
                )
                connectivity_values = np.asarray(
                    connectivity[upper], dtype=np.float64
                ).copy()
                grounded_precision = -connectivity[:n_modes, :n_modes]
                grounded_precision = 0.5 * (
                    grounded_precision + grounded_precision.T
                )
                try:
                    precision_factor = scipy.linalg.cho_factor(
                        grounded_precision, lower=True, check_finite=False
                    )
                except np.linalg.LinAlgError:
                    attempt_failed = True

            if not attempt_failed:
                grounded_covariance = scipy.linalg.cho_solve(
                    precision_factor,
                    np.eye(n_modes, dtype=np.float64),
                    check_finite=False,
                )
                grounded_covariance = 0.5 * (
                    grounded_covariance + grounded_covariance.T
                )
                metric, fitted_squared_distances, connectivity, _ = diagnostics()
                objective_decrease = (
                    previous_metric["objective"] - metric["objective"]
                )
                objective_scale = max(abs(previous_metric["objective"]), 1.0)
                monotonicity_tolerance = 1e-11 * objective_scale
                if objective_decrease < -monotonicity_tolerance:
                    attempt_failed = True

            if not attempt_failed:
                break
            sweep_relaxation *= 0.5
            sweep_retries += 1
            if sweep_relaxation < 2.0 ** -20:
                raise RuntimeError(
                    "coordinate sweep could not maintain a numerically stable "
                    "grounded precision after safeguarded under-relaxation"
                )

        coordinate_updates += pair_count
        if sweep_relaxation < 1.0:
            damped_sweeps += 1
            minimum_accepted_sweep_relaxation = min(
                minimum_accepted_sweep_relaxation, sweep_relaxation
            )
        metric["objective_decrease"] = max(objective_decrease, 0.0)
        metric["relative_objective_decrease"] = max(
            objective_decrease, 0.0
        ) / objective_scale
        metric["maximum_absolute_coordinate_update"] = (
            maximum_absolute_update
        )
        metric["minimum_coordinate_determinant_ratio"] = (
            minimum_determinant_ratio
        )
        metric["sweep_relaxation"] = sweep_relaxation

        history["iteration"].append(sweep)
        for key in history:
            if key != "iteration":
                history[key].append(metric[key])
        if sweep in save_steps_set:
            connectivity_at_steps[sweep] = connectivity.copy()
        if progress_callback is not None:
            progress_callback(
                {
                    "stage": "connectivity_coordinate_optimization",
                    "iteration": sweep,
                    "total": max_iterations,
                    "objective": metric["objective"],
                    "loss": metric["loss"],
                    "entropy": metric["entropy"],
                    "relative_gradient_norm": metric[
                        "relative_gradient_norm"
                    ],
                    "relative_stationarity_residual": metric[
                        "relative_stationarity_residual"
                    ],
                    "sweep_relaxation": sweep_relaxation,
                    "noisy": True,
                    "use_gpu": False,
                    "method": "CIS",
                    "general_method": "connectivity_coordinate_optimization",
                }
            )

        previous_metric = metric
        if (
            metric["relative_stationarity_residual"]
            <= stationarity_tolerance
        ):
            converged = True
            status = "stationarity_tolerance"
            message = "physical first-order stationarity tolerance reached"

    final_metric, fitted_squared_distances, connectivity, _ = diagnostics()
    gram = _centered_gram_from_squared_distances(fitted_squared_distances)
    for key, values in history.items():
        dtype = np.int64 if key == "iteration" else np.float64
        history[key] = np.asarray(values, dtype=dtype)

    initialization = dict(initialization)
    initialization["pinned_locus"] = n - 1
    iterations = int(history["iteration"][-1]) if history["iteration"].size else 0
    info = {
        "converged": bool(converged),
        "status": status,
        "message": message,
        "iterations": iterations,
        "coordinate_updates": int(coordinate_updates),
        "pairs_per_sweep": int(pair_count),
        "sweep_retries": int(sweep_retries),
        "damped_sweeps": int(damped_sweeps),
        "minimum_accepted_sweep_relaxation": float(
            minimum_accepted_sweep_relaxation
        ),
        "objective": final_metric["objective"],
        "connectivity_dual_objective": final_metric["objective"],
        "gradient_norm": final_metric["gradient_norm"],
        "gradient_infinity_norm": final_metric["gradient_infinity_norm"],
        "relative_gradient_norm": final_metric["relative_gradient_norm"],
        "stationarity_residual_norm": final_metric[
            "stationarity_residual_norm"
        ],
        "relative_stationarity_residual": final_metric[
            "relative_stationarity_residual"
        ],
        "maximum_absolute_stationarity_residual": final_metric[
            "maximum_absolute_stationarity_residual"
        ],
        "connectivity_l1": float(connectivity_l1),
        "connectivity_l1_penalty": final_metric[
            "connectivity_l1_penalty"
        ],
        "stationarity_tolerance": float(stationarity_tolerance),
        "observed_pair_count": pair_count,
        "noise_variance_kind": variance_kind,
        "noise_variance_minimum": float(np.min(pair_variance)),
        "noise_variance_median": float(np.median(pair_variance)),
        "noise_variance_maximum": float(np.max(pair_variance)),
        "initialization": initialization,
        "history": history,
        "connectivity_matrix_at_steps": connectivity_at_steps,
        "objective_definition": (
            "-1.5*logdet(P_grounded) + 0.5*sum_unique_pairs("
            "k_ij*D_obs_ij) + 0.125*sum_unique_pairs("
            "noise_variance_ij*k_ij**2) + connectivity_l1*"
            "sum_unique_pairs(abs(k_ij))"
        ),
        "stationarity_definition": (
            "soft_threshold(D_fit_ij-D_obs_ij, 2*connectivity_l1)-"
            "noise_variance_ij*k_ij/2"
        ),
        "coordinate_update_definition": (
            "exact cyclic minimization over the negative branch, exact-zero "
            "L1 kink, and positive branch using the matrix determinant lemma "
            "and Sherman-Morrison covariance update; a failed numerical SPD "
            "check retries the sweep under-relaxed toward each exact minimizer"
        ),
        "sweep_definition": "one deterministic pass over all unique pairs",
        "allows_signed_offdiagonal_connectivity": True,
        "logged_metric_timing": "after each complete coordinate sweep",
    }
    if not converged:
        warnings.warn(
            "fit_gaussian_noise_connectivity_coordinate_descent reached "
            "max_iterations before satisfying the physical stationarity "
            "tolerance",
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

    if eigvals.shape != (K.shape[0],):
        raise ValueError("eigvals must contain one stiffness eigenvalue per locus")

    # A centered connectivity matrix always has one translational (COM) mode.
    # Its numerically computed value can exceed a fixed absolute tolerance when
    # the stiffness spectrum becomes very large, so remove it by identity rather
    # than relying on the cutoff. Additional zero modes remain cutoff based.
    internal_eigvals = np.delete(eigvals, int(np.argmin(np.abs(eigvals))))
    positive_mask = internal_eigvals > zero_tol
    positive_eigvals = internal_eigvals[positive_mask]

    if len(positive_eigvals) == 0:
        max_eigval = np.max(internal_eigvals) if len(internal_eigvals) > 0 else 0.0
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
