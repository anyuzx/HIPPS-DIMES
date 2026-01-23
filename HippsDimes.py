"""
Reconstruction of 3D genome organization using the Maximum Entropy Principle

Reference:
1. Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051.
https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051
2. Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants." Nature Communications 14.1 (2023): 1150.
3. Shi, Guang, Sucheol, Shin, and D. Thirumalai. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advance 11 (31), eadx1763
"""

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
import click
import cooler # for cooler format of HiC data
import hicstraw # for .hic format of HiC data
from rich import print
from rich.panel import Panel
from rich.text import Text
from rich.console import Console
from rich.table import Table
#from tqdm.rich import trange, tqdm
from tqdm import trange, tqdm


console = Console()

#------------------------------------------------------------------#
# GPU Support (CuPy)
#------------------------------------------------------------------#
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
        _CUPY_GPU_NAME = gpu_props["name"].decode() if isinstance(gpu_props["name"], bytes) else gpu_props["name"]
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
    
    Omega = eigvector @ cp.diag(temp) @ eigvector.T
    Omega_diag = cp.diag(Omega)
    sigma = cp.sqrt(Omega_diag[:, cp.newaxis] + Omega_diag - 2.0 * Omega)
    dmap = 2.0 * cp.sqrt(2.0 / cp.pi) * sigma
    
    if return_eigenvalues:
        return dmap, eigvalue
    return dmap


#------------------------------------------------------------------#
# Helper functions
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
    eigvalue, eigvector = scipy.linalg.eigh(a)
    eigvalue_inv = 1.0 / eigvalue

    # difference in eigenvector components for monomers i and j
    vpi_vpj = eigvector[i, :] - eigvector[j, :]
    
    # normal_modes_square_mean = -(1 / eigenvalue) but filter out any inf
    normal_modes_square_mean = - np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)
    
    # Expand time dimension for broadcast
    t_reshaped = np.expand_dims(t, axis=-1)
    # Effective relaxation times
    tau_p = - zeta / eigvalue
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
    eigvalue, eigvector = scipy.linalg.eigh(a)
    eigvalue_inv = 1.0 / eigvalue
    vpi = eigvector[i, :]

    # Filter out infinities
    normal_modes_square_mean = - np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)

    # Expand time dimension for broadcast
    t_reshaped = np.expand_dims(t, axis=-1)
    tau_p = - zeta / eigvalue
    decay_factor = np.exp(-t_reshaped / tau_p)

    # The time-dependent part
    res = 3.0 * np.sum(vpi**2 * decay_factor * normal_modes_square_mean, axis=-1)
    # Equilibrium radius
    r2_eq = 3.0 * np.sum(vpi**2 * normal_modes_square_mean, axis=-1)
    # MSD
    msd_data = 2.0 * (r2_eq - res)

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
    lam, V = scipy.linalg.eigh(a)  # lam: eigenvalues, V: eigenvectors
    
    # Find zero (CM) mode
    p0 = np.argmin(np.abs(lam))
    if abs(lam[p0]) > tol:
        raise ValueError("No near-zero eigenvalue found; CM mode missing or 'a' not Rouse-type")
    
    # Get non-zero modes
    p_nz = np.arange(n) != p0  # boolean mask for non-zero modes
    lam_nz = lam[p_nz]        # non-zero eigenvalues
    V_nz = V[:, p_nz]         # non-zero eigenvectors
    
    # Compute equilibrium variances for non-zero modes
    sigma2_nz = -1.0 / lam_nz  # shape (n-1,)
    
    # Compute relaxation times
    tau = np.full(n, np.inf)   # initialize all to infinity
    tau[p_nz] = -zeta / lam_nz # set non-zero mode times
    
    # Compute decay factors for non-zero modes
    t_col = t[:, None]         # shape (len(t), 1)
    decay = np.exp(-t_col / tau[p_nz])  # shape (len(t), n-1)
    
    # Compute squared eigenvector components
    V2 = V_nz**2               # shape (n, n-1)
    
    # Compute time-dependent part from non-zero modes
    res = 3.0 * np.sum(
        V2[:, None, :] * decay[None, :, :] * sigma2_nz[None, None, :],
        axis=2
    )  # shape (n, len(t))
    
    # Compute equilibrium variance for each monomer
    r2_eq = 3.0 * np.sum(V2 * sigma2_nz[None, :], axis=1)  # shape (n,)
    
    # Compute center-of-mass diffusion
    D_cm = 1.0 / (zeta * n)
    msd_cm = 6.0 * D_cm * t    # shape (len(t),)
    
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
    eigvalue, eigvector = scipy.linalg.eigh(a)
    
    # Discard zero mode
    eigvalue = eigvalue[1:]
    eigvector = eigvector[:, 1:]
    
    # Compute weights W_p = (1/(j-i+1)) * sum_{m=i}^j V_{p,m}
    segment_length = j - i + 1
    W = np.sum(eigvector[i:j+1, :], axis=0) / segment_length

    # 1/eigvalue, filtering infinities (zero-mode)
    eigvalue_inv = 1.0 / eigvalue
    normal_modes_square_mean = -np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)

    # Expand time dimension for broadcasting
    t_reshaped = np.expand_dims(t, axis=-1)  # shape (len(t), 1)

    # Relaxation times tau_p
    tau_p = -zeta / eigvalue
    tau_p = np.nan_to_num(tau_p, posinf=0.0, neginf=0.0)  # Handle potential division by zero

    # Decay factor for each mode at each time
    decay_factor = np.exp(-t_reshaped / tau_p)  # shape (len(t), n)

    # Compute time-dependent part: 3 * sum(W^2 * decay_factor * normal_modes_square_mean, over modes)
    # Ensure proper broadcasting by reshaping W
    W_reshaped = W[None, :]  # shape (1, n)
    res = 3.0 * np.sum((W_reshaped**2) * decay_factor * normal_modes_square_mean, axis=-1)

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
    lam = eigvals[1:]          # drop zero mode
    vec = eigvecs[:, 1:]       # shape (n, m)

    # 2. Precompute E_p(i,j) and W0 = E_p/λ_p
    diff = vec[:, None, :] - vec[None, :, :]  # shape (n,n,m)
    E = diff * diff                           # (n,n,m)
    W0 = E / lam[None, None, :]               # (n,n,m)

    # 3. Time grid
    t_min = 0.0
    tau_slow = 1.0 / np.min(lam)
    t_max = factor * tau_slow
    t_vals = np.logspace(np.log10(t_min + 1e-12), np.log10(t_max), num_t)

    # 4. Compute normalized G for each t
    denom = np.sum(W0, axis=2)               # (n,n)
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
                t0, t1 = t_vals[k-1], t_vals[k]
                G0, G1 = G_norm[i, j, k-1], G_norm[i, j, k]
                # linear interpolation
                tau_mat[i, j] = t0 + (1/np.e - G0) * (t1 - t0) / (G1 - G0 + tol)
                
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
    lam = eigvals[1:]             # shape (n-1,)
    vec = eigvecs[:, 1:]          # shape (n,   n-1)

    # 3) Precompute inverse powers
    inv_lam  = 1.0 / lam           # 1/λ_p
    inv_lam2 = inv_lam ** 2        # 1/λ_p^2
    inv_lam3 = inv_lam ** 3        # 1/λ_p^3

    # 4) Compute squared differences for all (i,j,p)
    #    E[i,j,p] = (V[p,i] - V[p,j])^2
    #    vec shape -> (n, m), we want (n, n, m)
    diff = vec[:, None, :] - vec[None, :, :]  # shape (n, n, m)
    E = diff * diff                          # squared differences

    # 5) Weighted sums via tensordot over the mode axis
    sum0 = np.tensordot(E, inv_lam,  axes=([2], [0]))  # ∑ E/λ_p
    sum1 = np.tensordot(E, inv_lam2, axes=([2], [0]))  # ∑ E/λ_p^2
    sum2 = np.tensordot(E, inv_lam3, axes=([2], [0]))  # ∑ E/λ_p^3

    # 6) Form relaxation times
    tau1 = sum1 / sum0  # first moment
    tau2 = sum2 / sum1  # second moment ratio

    np.fill_diagonal(tau1, 0.0)
    np.fill_diagonal(tau2, 0.0)

    return tau1, tau2

def compute_tau_ij_integral(a, i, j):
    # a is the connectivity matrix
    eigval, eigvec = np.linalg.eigh(-a) # note negative sign

    lam = eigval[1:] # ignore the p=0 mode
    vec = eigvec[:, 1:] # remove the eigvector corresponding to p=0 mode

    v_pi = vec[i, :] # V_{p,i}
    v_pj = vec[j, :] # V_{p,j}
    E_ps = (v_pi - v_pj) ** 2

    sum0 = np.sum(E_ps / lam)
    sum1 = np.sum(E_ps / lam ** 2.)
    sum2 = np.sum(E_ps / lam ** 3.)

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
    lam = eigvals[1:]                # shape (n-1,)
    vec = eigvecs[:, 1:]             # shape (n,   n-1)

    # 3) Mode-difference squared for pair (i,j)
    v_i = vec[i, :]                  # (n-1,)
    v_j = vec[j, :]                  # (n-1,)
    E_p = (v_i - v_j)**2             # (n-1,)

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
    tau_ij = scipy.optimize.brentq(lambda t: G_norm(t) - 1/np.e, t_min, t_max, xtol=tol)
    return tau_ij

# ------------------------------
# Modulus
def compute_modulus(a: np.ndarray, freq: np.ndarray, zeta: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
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
    eigvalue, eigvector = scipy.linalg.eigh(a)
    
    # Exclude the last eigenvalue (zero eigenvalue) and corresponding eigenvector
    eigvalue = eigvalue[:-1]             # Shape: (n_modes,)
    eigvector = eigvector[:, :-1]        # Shape: (n_monomers, n_modes)
    
    # Compute normal modes and relaxation times
    eigvalue_inv = 1.0 / eigvalue
    normal_modes_square_mean = -eigvalue_inv
    tau_p = -zeta / eigvalue
    
    # Reshape frequency for broadcasting
    freq_reshaped = np.expand_dims(freq, axis=-1)
    
    # Compute moduli
    storage_modulus = np.sum((freq_reshaped * tau_p) ** 2. / 
                           (1 + (freq_reshaped * tau_p) ** 2), axis=-1)
    loss_modulus = np.sum((freq_reshaped * tau_p) / 
                         (1 + (freq_reshaped * tau_p) ** 2), axis=-1)
    
    return (np.column_stack((freq, storage_modulus)), 
            np.column_stack((freq, loss_modulus)))

def compute_monomer_modulus(a: np.ndarray, freq: np.ndarray, zeta: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the storage and loss moduli for individual monomers in a polymer system.
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
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Three arrays containing:
        - First array: frequencies
        - Second array: storage modulus for each monomer (shape: n_freqs × n_monomers)
        - Third array: loss modulus for each monomer (shape: n_freqs × n_monomers)
        
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
    eigvals, eigvecs = scipy.linalg.eigh(a)
    
    # Exclude the last eigenvalue (zero eigenvalue) and corresponding eigenvector
    eigvals = eigvals[:-1]             # Shape: (n_modes,)
    eigvecs = eigvecs[:, :-1]          # Shape: (n_monomers, n_modes)
    
    # Relaxation times τ_p
    tau_p = -zeta / eigvals           # Shape: (n_modes,)
    
    # Compute ωτ_p for all frequencies and modes
    omega_tau_p = freq[:, np.newaxis] * tau_p[np.newaxis, :]  # Shape: (n_freqs, n_modes)
    
    # Calculate f_p(ω) and g_p(ω)
    f_p = (omega_tau_p ** 2) / (1 + omega_tau_p ** 2)         # For storage modulus
    g_p = omega_tau_p / (1 + omega_tau_p ** 2)                # For loss modulus
    
    # Square of eigenvector components (V_{pi}^2)
    eigvecs_squared = eigvecs ** 2                            # Shape: (n_monomers, n_modes)
    
    # Compute G'_i(ω) and G''_i(ω) using Einstein summation
    G_prime_i = np.einsum('mp,fp->fm', eigvecs_squared, f_p)        # Shape: (n_freqs, n_monomers)
    G_double_prime_i = np.einsum('mp,fp->fm', eigvecs_squared, g_p) # Shape: (n_freqs, n_monomers)
    
    return freq, G_prime_i, G_double_prime_i

# ------------------------------


def Ornstein_Uhlenbeck_update(x, dt, k, zeta, beta, force_projection = 0.0, method='euler-maruyama', update_zero_modes=True, zero_mode_tol=1e-10):
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
    
    **Unstable Modes (k < 0):**
    When the connectivity matrix has positive eigenvalues (e.g., from broken bonds), k < 0.
    This causes exponential growth in those modes. The 'euler-maruyama' method handles this
    with controlled growth over small time steps. The 'exact' method should NOT be used with
    unstable modes as the noise term becomes undefined.
    """
    if isinstance(x, np.ndarray):
        rand_noise = np.random.randn(*x.shape)
    else:
        rand_noise = np.random.randn()
    
    theta = k[:, np.newaxis] / zeta
    
    # Identify zero eigenvalue modes and unstable modes (k < 0)
    zero_modes_mask = np.abs(k) < zero_mode_tol
    unstable_modes_mask = k < -zero_mode_tol  # Negative spring constants = positive eigenvalues
    
    # Warn if using exact method with unstable modes
    if method == 'exact' and np.any(unstable_modes_mask):
        import warnings
        warnings.warn(
            f"Detected {np.sum(unstable_modes_mask)} unstable modes (positive eigenvalues) with 'exact' method. "
            "This may cause numerical issues. Consider using 'euler-maruyama' method for systems with broken bonds.",
            RuntimeWarning
        )

    if method == 'euler-maruyama':
        # Euler-Maruyama: dX = -θX dt + (force_projection/ζ) dt + σ dW
        # No need to divide by eigenvalue here
        dx = - theta * x * dt + force_projection * dt / zeta + np.sqrt(2.0 * dt / (zeta * beta)) * rand_noise
        x_new = x + dx
    elif method == 'exact':
        sigma = (2. / (zeta * beta)) ** .5
        mu = np.exp(- theta * dt)
        
        # Calculate force drift term: μ(1 - e^(-θt))
        # For non-zero modes: μ = force_projection/k, so drift = (force_projection/k)(1 - e^(-θt))
        # For zero modes (k ≈ 0): drift reduces to (force_projection/ζ) * t
        if isinstance(force_projection, np.ndarray):
            force_drift = np.zeros_like(force_projection)
            non_zero_mask = ~zero_modes_mask
            
            # Non-zero modes: (force_projection/k)(1 - e^(-θt))
            force_drift[non_zero_mask] = (force_projection[non_zero_mask] / k[non_zero_mask, np.newaxis]) * (1 - mu[non_zero_mask])
            
            # Zero modes: (force_projection/ζ) * t
            force_drift[zero_modes_mask] = (force_projection[zero_modes_mask] / zeta) * dt
        else:
            # force_projection is a scalar (0.0)
            force_drift = 0.0
        
        # For non-zero modes: use standard Gillespie formula
        noise_term = np.sqrt((sigma ** 2. / (2. * theta)) * (1. - mu ** 2.))
        
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
    A = np.diag(np.full(n-1, k), 1)
    A += A.T
    A[np.diag_indices(n)] = -2*k
    A[0, 0] = -k
    A[n-1, n-1] = -k
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
    return (sigma_row_sum[:, np.newaxis] + sigma_row_sum - sigma_sum / n) / (2 * n) - sigma_mtx_square / 2.0


def dmap2a_direct(dmap):
    """
    Return connectivity matrix A given the mean distance map directly through matrix peusudo inversion
    """
    sigma_mtx = 0.5 * np.sqrt(np.pi / 2.0) * dmap
    Omega = sigma2omega(sigma_mtx)
    a_direct = nearestNSD(- scipy.linalg.pinvh(Omega), 0.0)

    return a_direct


def ddmap2a_direct(ddmap):
    sigma_mtx = np.sqrt(ddmap / 3.0)
    Omega = sigma2omega(sigma_mtx)
    a_direct = nearestNSD(- scipy.linalg.pinvh(Omega), 0.0)

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
    eigvalue, eigvector = scipy.linalg.eigh(A)

    temp = -1.0 / eigvalue

    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[temp >= TOL] = 0.0
    temp[temp <= -TOL] = 0.0
    #temp[np.abs(temp) <= 10**-7] = 0.0

    # replace all positive element to be zero
    if force_positive_definite:
        temp[temp < 0.0] = 0.0

    Omega = eigvector @ np.diag(temp) @ eigvector.T
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
        raise TypeError("force must be a dictionary with keys: 'loci', 'amplitude', 'direction'")
    required_keys = ['loci', 'amplitude', 'direction']
    for key in required_keys:
        if key not in force:
            raise ValueError(f"force dictionary must contain key: {key}")
    
    # Extract force parameters
    force_loci = force['loci']
    force_amplitude = force['amplitude']
    force_direction = np.array(force['direction'])
    
    # Normalize force direction
    force_direction = force_direction / np.linalg.norm(force_direction)
    
    # Create force vector B in real space
    N = A.shape[0]
    B = np.zeros((N, 3))
    for locus in force_loci:
        if locus < 0 or locus >= N:
            raise ValueError(f"force locus {locus} is out of range [0, {N-1}]")
        B[locus] = force_amplitude * force_direction
    
    # Eigendecomposition
    eigvalue, eigvector = scipy.linalg.eigh(A)
    
    # Compute -1/eigenvalue for thermal fluctuations (Ω matrix)
    temp = -1.0 / eigvalue
    
    # Handle infinities and large values (zero eigenvalue handling)
    temp[temp == -np.inf] = 0.0
    temp[temp == np.inf] = 0.0
    temp[np.abs(temp) >= TOL] = 0.0
    
    # Compute Ω matrix for thermal fluctuations: Ω = V * diag(-1/λ) * V^T
    Omega = eigvector @ np.diag(temp) @ eigvector.T
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
    
    # Equilibrium displacement contribution: (R_eq,i - R_eq,j)
    # This is the shift in mean position due to force
    delta_R_eq = R_eq[:, np.newaxis, :] - R_eq[np.newaxis, :, :]  # Shape: (N, N, 3)
    
    # For the distance calculation, we need to consider that:
    # - Thermal fluctuations contribute equally in all 3 dimensions: 3*σ²_thermal
    # - Equilibrium shift is a fixed vector displacement
    # 
    # The mean distance for a 3D Gaussian with mean μ and isotropic variance σ² in each dim:
    # If there's also a displacement d, the distribution becomes non-central chi
    # For simplicity, we compute the squared distance from equilibrium positions
    # and add thermal fluctuation variance
    
    # Compute squared equilibrium separation
    delta_R_eq_sq = np.sum(delta_R_eq**2, axis=2)  # Shape: (N, N)
    
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
    cmap = scipy.special.erf(rc/(np.sqrt(2) * sigma_mtx)) - \
        np.sqrt(2.0/np.pi) * np.exp(-0.5 * rc**2.0 /
                                    np.power(sigma_mtx, 2.0)) * rc / sigma_mtx
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
    """
    Correct the connectivity matrx. Make it Laplacian, and non negative (options)
    """
    temp = np.copy(a)
    if fill_negative:
        temp[temp < 0.0] = 0.0
    # Zero out diagonal first, then set diagonal to negative sum of off-diagonal elements
    np.fill_diagonal(temp, 0.0)
    np.fill_diagonal(temp, -np.sum(temp, axis=1))
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
    #U = V * Wt
    U = np.dot(V, Wt)

    if not return_rotation:
        return np.array(P * U + Qc)
    else:
        return np.array(P * U + Qc), U


def write2xyz(fout, xyzs):
    natoms = xyzs.shape[1]

    xyz0 = xyzs[0]

    with open(fout, 'w') as f:
        for snapshot in xyzs:
            xyz = optimal_rotate(snapshot, xyz0, allow_reflection=True)
            f.write('{}\n\n'.format(natoms))
            for idx, item in enumerate(xyz):
                f.write('{} {} {} {}\n'.format('C', item[0], item[1], item[2]))

def write2xyz_traj(fout, traj):
    natoms = traj.shape[1]

    with open(fout, 'w') as f:
        for snapshot in traj:
            f.write(f"{natoms}\n\n")
            for i, atom in enumerate(snapshot):
                f.write(f"C {atom[0]} {atom[1]} {atom[2]}\n")


def a2xyz_sample(A, ensemble=1, force_positive_definite=False):
    """
    Function to generate an ensemble of configurations given the connectivity matrix
    """
    TOL = 10**8.0
    eigvalue, eigvector = scipy.linalg.eigh(A)
    temp = 1.0/eigvalue[:, np.newaxis]

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
        position = eigvector @ (np.sqrt(-temp) *
                                np.random.randn(len(eigvalue), 3))
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
        raise TypeError("force must be a dictionary with keys: 'loci', 'amplitude', 'direction'")
    required_keys = ['loci', 'amplitude', 'direction']
    for key in required_keys:
        if key not in force:
            raise ValueError(f"force dictionary must contain key: {key}")
    
    # Extract force parameters
    force_loci = force['loci']
    force_amplitude = force['amplitude']
    force_direction = np.array(force['direction'])
    
    # Normalize force direction
    force_direction = force_direction / np.linalg.norm(force_direction)
    
    # Create force vector B in real space
    N = A.shape[0]
    B = np.zeros((N, 3))
    for locus in force_loci:
        if locus < 0 or locus >= N:
            raise ValueError(f"force locus {locus} is out of range [0, {N-1}]")
        B[locus] = force_amplitude * force_direction
    
    # Eigendecomposition
    eigvalue, eigvector = scipy.linalg.eigh(A)
    
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

def a2xyz_sample_fixed_end(A,
                           xyz_start,
                           xyz_end,
                           ensemble=1,
                           force_positive_definite=False):
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
    A_copy[0,   0]   += w
    A_copy[-1, -1]   += w

    # 2) Eigendecompose
    evals, evecs = scipy.linalg.eigh(A_copy)

    # 3) Build the "temp" = 1/λ matrix, clean up infinities and huge values
    temp = 1.0 / evals[:, None]
    temp[np.isinf(temp)] = 0.0
    TOL = 1e8
    temp[np.abs(temp) >= TOL] = 0.0
    if force_positive_definite:
        temp[temp > 0.0] = 0.0

    # 4) Compute end‐to‐end distance and linear "b" term
    xyz_start = np.array(xyz_start, float)
    xyz_end   = np.array(xyz_end,   float)
    L = np.linalg.norm(xyz_end - xyz_start)

    n = len(evals)
    b = np.zeros((n, 3))
    b[-1, 2] = -w * L   # ensures the last bead sits at z=L in the internal frame

    # 5) Sample `ensemble` positions
    out = []
    for _ in range(ensemble):
        # random + shift in eigenspace
        coeff = np.sqrt(-temp) * np.random.randn(n, 3)
        shift = (evecs.T @ b) * (-1.0 / evals)[:, None]
        xyz   = evecs @ (coeff + shift)

        # 6) Hard‐set the two ends along the z‐axis. This is not strictly necessary if w is large enough.
        xyz[0]    = [0.0, 0.0, 0.0]
        xyz[-1]   = [0.0, 0.0, L]

        # 7) xyz now is oriented in such way that the first monomer is at [0,0,0] and last is at [0,0,L]
        #    Hence we need to rotate+translate into the user‐specified endpoints
        ref1 = np.array([[0,0,0], [0,0,L]])
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
        (x1, y1), newarr.ravel(), (xx, yy), method='nearest')
    return GD1


def objective_func(rc, A_mtx, cmap_exp):
    x = a2cmap_theory(A_mtx, rc)
    y = cmap_exp / np.nanmax(cmap_exp)
    logx = interpolate_missing(np.log(x))
    logy = interpolate_missing(np.log(y))
    res = np.power(logx[np.triu_indices_from(logx, k=1)] -
                   logy[np.triu_indices_from(logy, k=1)], 2.).mean()**0.5
    return res

# FUNCTION TO CONVERT CMAP TO DMAP
def cmap2dmap_core(cmap_exp, rc, alpha, not_normalize, norm_max=1.0, mode='log'):
    # rc is the prefactor
    # norm_max is the maximum contact probability
    if mode == 'raw':
        if not_normalize:
            log10_pmap = np.log10(cmap_exp)
        else:
            log10_pmap = np.log10(
                cmap_exp) + np.log10(norm_max) - np.log10(np.nanmax(cmap_exp))
    elif mode == 'log':
        if not_normalize:
            log10_pmap = np.copy(cmap_exp)
        else:
            log10_pmap = cmap_exp + np.log10(norm_max) - np.nanmax(cmap_exp)

    return rc * 10 ** (-1.0/alpha * log10_pmap)


def cmap2dmap(cmap, alpha, not_normalize):
    # cmap is the raw data
    # we take log on contact map
    # and then interpolate the missing data. Any zero contact pair will be interpolated
    cmap_log = interpolate_missing(np.log10(cmap))
    cmap_log = np.array((cmap_log + cmap_log.T) / 2.)
    # lastly, convert to distance map using value of alpha
    dmap = cmap2dmap_core(cmap_log, 1.0, alpha, not_normalize)
    return dmap


def cmap2dmap_missing_data(cmap, alpha, not_normalize):
    # cmap is the raw data
    # we take log on contact map
    # unlike cmap2dmap(), this function does not interpolate the missing data. Just leave the missing data as is
    cmap_log = np.log10(cmap)
    cmap_log = np.array((cmap_log + cmap_log.T) / 2.)
    # convert to distance map using value of alpha
    dmap = cmap2dmap_core(cmap_log, 1.0, alpha, not_normalize)
    return dmap


def nearestNSD(X, delta):
    v, w = scipy.linalg.eigh(X)
    v_new = np.minimum(v, delta)
    return w @ np.diag(v_new) @ w.T


def ddmap2cov(ddmap):
    # convert a squared distance map to covariance matrix
    n = ddmap.shape[0]
    omega2_mtx = ddmap / 3.
    omega2_row_sum = np.sum(omega2_mtx, axis=1)
    omega2_sum = np.sum(omega2_mtx)
    return (omega2_row_sum[:, np.newaxis] + omega2_row_sum - omega2_sum / n) / (2 * n) - omega2_mtx / 2.0


def dmap2cov(dmap):
    # convert a distance map to covariance matrix
    ddmap = (3. * np.pi / 8.) * np.power(dmap, 2.)
    return ddmap2cov(ddmap)


def checkEMD(ddmap):
    # check whether a squared distance map is a valid Euclidean matrix
    cov = ddmap2cov(ddmap)
    eigvalue, eigvector = scipy.linalg.eigh(cov)
    # print(eigvalue)
    if np.all(eigvalue >= -0.1):
        if eigvalue.min() < 0.0:
            console.print("[red]Warning: The smallest eigenvalue of the covariance matrix is negative. \
Direct inversion method [italic]may[/italic] not work. Check the final results. \
If the results are not good enough, please try iterative scaling or gradient descent method.")
        return True
    else:
        return False
    
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
    if np.linalg.cond(A_CC) > 1/tol:
        A_CC_inv = np.linalg.pinv(A_CC, rcond=tol)
    else:
        A_CC_inv = np.linalg.inv(A_CC)

    # Schur complement
    A_eff = A_SS - A_SC @ A_CC_inv @ A_CS
    return A_eff

def neighbor_balance_symmetric(C, *, not_normalize=False, circular=False, epsilon=1e-12, return_scales=False):
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
        left_idx  = (idx - 1) % n
        right_idx = (idx + 1) % n
        s = 0.5 * (X[idx, left_idx] + X[idx, right_idx])
    else:
        s = np.empty(n, dtype=float)
        # interior i=1..n-2 have both neighbors
        if n >= 3:
            s[1:-1] = 0.5 * (np.diag(X, k=-1)[:-1] + np.diag(X, k=1)[1:])
        elif n == 2:
            # degenerate case: each has only one neighbor
            s[:] = np.array([X[0,1], X[1,0]])
        else:
            # n == 1: no neighbors; avoid divide by zero by using epsilon
            s[:] = epsilon

        # edges: use the single available neighbor
        if n >= 2:
            s[0]   = X[0, 1]
            s[-1]  = X[-1, -2]

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

#------------------------------------------------------------------#


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
        eigvals = scipy.linalg.eigh(K, eigvals_only=True)
    
    positive_mask = eigvals > zero_tol
    positive_eigvals = eigvals[positive_mask]
    
    if len(positive_eigvals) == 0:
        max_eigval = np.max(eigvals) if len(eigvals) > 0 else 0.0
        if max_eigval < zero_tol:
            return -np.inf
        return np.nan
    
    with np.errstate(divide='ignore', invalid='ignore'):
        log_terms = -np.log(positive_eigvals)
    
    if np.any(~np.isfinite(log_terms)):
        return np.nan
    
    entropy = np.sum(log_terms)
    if not np.isfinite(entropy):
        return np.nan
    return entropy


class Optimize:
    """
    Optimizer for finding connectivity matrix A that matches target distance constraints.
    
    Uses Maximum Entropy principle to find the Gaussian polymer model that best
    matches the target mean squared distance map.
    
    Optimization Methods
    --------------------
    - 'IS' (Iterative Scaling): Default method. Multiplicative updates based on
      log-ratio of current vs target distances. Generally stable and reliable.
    - 'GD' (Gradient Descent): Additive updates with Nesterov momentum.
      Requires smaller learning rates (e.g., 1e-8).
    - 'DI' (Direct Inversion): One-shot pseudo-inverse of covariance matrix.
      Fast but may not satisfy all constraints; best for valid Euclidean distance matrices.
    
    Convergence Acceleration
    ------------------------
    **Momentum with Nesterov** (RECOMMENDED for fastest convergence):
    - Nesterov uses "look-ahead" correction to prevent overshooting
    - Enables higher momentum values (0.95) that would diverge with standard momentum
    - BEST setting: momentum=0.95, nesterov=True (~50% faster than momentum=0.9)
    
    **Standard Momentum** (fallback option):
    - Uses Polyak's heavy ball method
    - Safe with momentum=0.9, provides ~2-4x speedup over baseline
    
    Examples
    --------
    >>> # Standard optimization (no acceleration)
    >>> opt = Optimize(ddmap_target)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS')
    
    >>> # With Nesterov momentum (RECOMMENDED - fastest)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS', momentum=0.95, nesterov=True)
    
    >>> # With standard momentum (fallback)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS', momentum=0.9)
    """
    
    def __init__(self, ddmap_target, connectivity_matrix=None, use_gpu=False):
        # ddmap_target is the targeted matrix we would like to match
        # note that ddmap_taret is the mean SQUARED distance matrix, not mean distance matrix
        self.ddmap_target = ddmap_target

        # get the size of system
        self.n = ddmap_target.shape[0]

        if connectivity_matrix is None:
            # initialize the connectivity matrix
            # here the connectivity matrix is initialized as a simple rouse chain whose spring constant is determined such\
            # that its radius of gyration is close to the target
            # we need to filter out both NaN and inf entries
            rg2 = .5 * np.nanmean(self.ddmap_target[~np.isinf(self.ddmap_target)])
            k = self.n / (4. * rg2)
            self.A = construct_connectivity_matrix_rouse(self.n, k)
        else:
            self.A = connectivity_matrix

        # initialize the loss
        self.loss = None


        # Optional off-diagonal mask: True means keep/update A[i,j]; False means freeze at 0
        # Diagonal entries are always recomputed by a2a()
        self.edge_mask = None
        
        # GPU support
        self.use_gpu = use_gpu and is_gpu_available()
        if use_gpu and not is_gpu_available():
            console.print("[yellow]Warning: use_gpu=True but CuPy is not available. Falling back to CPU.[/yellow]")
        
        if self.use_gpu:
            # Move data to GPU
            self._A_gpu = cp.asarray(self.A)
            self._ddmap_target_gpu = cp.asarray(self.ddmap_target)
            self._velocity_gpu = None
            cp.cuda.Stream.null.synchronize()

    def set_edge_mask(self, edge_mask):
        """Set/update the off-diagonal mask for optimization.

        Parameters
        ----------
        edge_mask : (n, n) array_like of bool or None
            True => keep/allow this off-diagonal entry A[i,j] to be optimized (i!=j).
            False => freeze A[i,j]=0 (i!=j). The diagonal is ignored and recomputed by a2a().
            If None, no masking is applied (all off-diagonals are free).
        """
        if edge_mask is None:
            self.edge_mask = None
            return

        m = np.array(edge_mask, dtype=bool, copy=True)
        if m.shape != (self.n, self.n):
            raise ValueError(f"edge_mask must have shape ({self.n},{self.n}), got {m.shape}")

        # Symmetrize and ignore diagonal
        m = np.logical_and(m, m.T)
        np.fill_diagonal(m, False)
        self.edge_mask = m

    def _freeze_masked_edges(self):
        """Force masked (disallowed) off-diagonals to be exactly zero."""
        if self.edge_mask is None:
            return
        freeze = ~self.edge_mask
        np.fill_diagonal(freeze, False)
        self.A[freeze] = 0.0

    #def run_masked(self, epoch, edge_mask, ddmap_target_masked=None, general_method='optimization', **kwargs):
        """Run optimization with a fixed off-diagonal mask (useful for greedy pruning).

        Parameters
        ----------
        epoch : int
            Number of iterations.
        edge_mask : (n,n) bool
            Off-diagonal mask; False => freeze A_ij = 0 for i!=j.
        ddmap_target_masked : (n,n) float or None
            If provided, overrides self.ddmap_target. Use NaN on unconstrained pairs.
        general_method : str
            Only 'optimization' is supported for masked runs. ('direct' requires full EMD).
        """
    #    if ddmap_target_masked is not None:
    #        if ddmap_target_masked.shape != (self.n, self.n):
    #            raise ValueError("ddmap_target_masked must match the current system size")
    #        self.ddmap_target = ddmap_target_masked

    #    self.set_edge_mask(edge_mask)
    #    return self.run(epoch, general_method=general_method, **kwargs)

    def run_masked(self, epoch, edge_mask, ddmap_target_masked=None, general_method='optimization', **kwargs):
        """Run optimization with a fixed off-diagonal mask, optimized for *sparse* constraints.

        Replaces the previous wrapper that called `run()` (which recomputed the full
        distance map via an eigendecomposition at every iteration).

        When the number of constrained pairs is small (typical in minimax / greedy growing),
        we compute only the needed mean-squared distances using grounded Laplacian solves:

            ddmap(i,j) = <||x_i-x_j||^2> = 3 * (e_i-e_j)^T (-A)^+ (e_i-e_j)

        Parameters
        ----------
        epoch : int
            Number of iterations.
        edge_mask : (n,n) bool
            Off-diagonal mask; True => allow/update A_ij for i!=j, False => freeze A_ij=0.
        ddmap_target_masked : (n,n) float or None
            If provided, overrides self.ddmap_target. Use NaN/inf on unconstrained pairs.
        general_method : str
            Only 'optimization' is supported here.
        kwargs :
            learning_rate : float
            lamd : float
            reg : {'L1','L2'}
            method : {'IS'}  (other methods fall back to the full `run()` implementation)
            enforce_nonnegative_connectivity_matrix : bool
            pin : int (default 0)  -- grounded node for Laplacian solve
            jitter : float (default 1e-10) -- diagonal jitter for numerical stability
        """
        if ddmap_target_masked is not None:
            if ddmap_target_masked.shape != (self.n, self.n):
                raise ValueError("ddmap_target_masked must match the current system size")
            self.ddmap_target = ddmap_target_masked

        if general_method != 'optimization':
            self.set_edge_mask(edge_mask)
            loss_array, entropy_array, dmap_maxent, A_final = self.run(epoch, general_method=general_method, **kwargs)
            return loss_array, entropy_array, dmap_maxent, A_final

        method = kwargs.get('method', 'IS')
        if method != 'IS':
            self.set_edge_mask(edge_mask)
            loss_array, entropy_array, dmap_maxent, A_final = self.run(epoch, general_method=general_method, **kwargs)
            return loss_array, entropy_array, dmap_maxent, A_final

        learning_rate = float(kwargs.get('learning_rate', 1.0))
        lamd = float(kwargs.get('lamd', 0.0))
        reg = str(kwargs.get('reg', 'L2')).upper()
        enforce_nonnegative = bool(kwargs.get('enforce_nonnegative_connectivity_matrix', False))
        pin = int(kwargs.get('pin', 0))
        jitter = float(kwargs.get('jitter', 1e-10))

        self.set_edge_mask(edge_mask)

        # -------- Determine constrained pairs (i<j) --------
        tgt = np.array(self.ddmap_target, dtype=float, copy=False)
        if self.edge_mask is None:
            return self.run(epoch, general_method=general_method, **kwargs)

        finite = np.isfinite(tgt) & (tgt > 0.0)
        allow = self.edge_mask & finite
        ii, jj = np.where(np.triu(allow, 1))
        m = ii.size
        if m == 0:
            self._freeze_masked_edges()
            self.A = a2a(0.5 * (self.A + self.A.T), fill_negative=enforce_nonnegative)
            dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)
            # Compute entropy
            K = -self.A
            eigvals_K = scipy.linalg.eigh(K, eigvals_only=True)
            entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)
            return [], [entropy], dmap_maxent, self.A

        dd_target_edges = tgt[ii, jj].astype(float)

        # -------- Precompute grounded incidence matrix B for constrained pairs --------
        if not (0 <= pin < self.n):
            raise ValueError(f"pin must be in [0,{self.n-1}]")

        keep = np.ones(self.n, dtype=bool)
        keep[pin] = False

        pos = -np.ones(self.n, dtype=int)
        pos[np.where(keep)[0]] = np.arange(self.n - 1)

        B = np.zeros((self.n - 1, m), dtype=float)
        cols = np.arange(m)

        gi = pos[ii]
        gj = pos[jj]
        valid_i = gi >= 0
        valid_j = gj >= 0
        if np.any(valid_i):
            B[gi[valid_i], cols[valid_i]] += 1.0
        if np.any(valid_j):
            B[gj[valid_j], cols[valid_j]] -= 1.0

        # -------- Main loop --------
        loss_array = []
        entropy_array = []
        with trange(epoch, desc="Performing masked optimization (sparse)", unit="iteration") as pbar:
            for t in pbar:
                # Grounded Laplacian Lg = (-A) with pinned node removed
                A_sym = 0.5 * (self.A + self.A.T)
                L = -A_sym
                Lg = L[np.ix_(keep, keep)].copy()
                if jitter > 0.0:
                    Lg.flat[::Lg.shape[0] + 1] += jitter

                # Solve Lg X = B for all constrained edges at once
                try:
                    c, lower = scipy.linalg.cho_factor(Lg, lower=True, check_finite=False, overwrite_a=False)
                    X = scipy.linalg.cho_solve((c, lower), B, check_finite=False)
                except np.linalg.LinAlgError:
                    # Fallback: eigen-based pseudo-inverse on the grounded system
                    w, V = scipy.linalg.eigh(Lg, check_finite=False)
                    w_inv = np.zeros_like(w)
                    good = w > jitter
                    w_inv[good] = 1.0 / w[good]
                    X = (V * w_inv) @ (V.T @ B)

                # ddmap on constrained edges: dd = 3 * b^T Lg^{-1} b
                rho = np.sum(B * X, axis=0)
                dd_cur_edges = 3.0 * np.maximum(rho, 0.0)

                # Iterative-scaling gradient on constrained edges
                with np.errstate(divide='ignore', invalid='ignore'):
                    log_ratio = np.log(dd_cur_edges / dd_target_edges)
                log_ratio = np.nan_to_num(log_ratio, posinf=0.0, neginf=0.0)

                fhash = 0.5 * float(np.sum(dd_cur_edges))
                if fhash <= 0.0:
                    fhash = 1.0

                g = log_ratio

                # Regularization (matches __update_parameter() on allowed edges)
                if lamd > 0.0:
                    aij = A_sym[ii, jj]
                    if reg == 'L2':
                        g = g - 2.0 * lamd * aij
                    elif reg == 'L1':
                        g = g + lamd * np.sign(-aij)

                g = g / fhash

                # Update only allowed off-diagonals
                self.A[ii, jj] += learning_rate * g
                self.A[jj, ii] += learning_rate * g

                # Cleanup + enforce mask + rebuild diagonal
                self.A = np.nan_to_num(self.A)
                self.A = 0.5 * (self.A + self.A.T)
                self._freeze_masked_edges()
                self.A = a2a(self.A, fill_negative=enforce_nonnegative)

                # Loss on constrained edges only
                with np.errstate(divide='ignore', invalid='ignore'):
                    rel = (dd_cur_edges - dd_target_edges) / dd_target_edges
                    loss = float(np.nanmean(rel * rel) ** 0.5)
                self.loss = loss
                
                # Compute entropy using eigenvalues of -A (K = -A)
                K = -self.A
                eigvals_K = scipy.linalg.eigh(K, eigvals_only=True)
                entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)
                self.entropy = entropy
                
                pbar.set_postfix(loss=self.loss, entropy=self.entropy if self.entropy is not None else np.nan)
                loss_array.append(self.loss)
                entropy_array.append(self.entropy if self.entropy is not None else np.nan)

        # Compute full mean distance map once at the end (compat)
        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)
        return loss_array, entropy_array, dmap_maxent, self.A
    def __compute_loss(self, ddmap_t):
        with np.errstate(divide='ignore', invalid='ignore'):
            loss = np.nanmean(
                np.power((ddmap_t - self.ddmap_target)/self.ddmap_target, 2.)) ** .5
        return loss

    def __update_parameter(self, t, learning_rate, lamd=0.0, reg='l2', method='IS', enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        # updating using Iterative Scaling or Gradient Descent
        
        if self.use_gpu and method == 'IS':
            # GPU path for IS using CuPy
            return self.__update_parameter_gpu(t, learning_rate, lamd, reg, enforce_nonnegative_connectivity_matrix, momentum, nesterov)
        elif self.use_gpu and method == 'GD':
            # GPU path for GD using CuPy
            return self.__update_parameter_gpu_gd(t, learning_rate, lamd, reg, enforce_nonnegative_connectivity_matrix)

        # CPU path (original implementation)
        # compute the mean squared distance matrix at current iteration step
        # Also get eigenvalues of A to reuse for entropy computation (avoiding double eigendecomposition)
        dmap_t, eigvals_A = a2dmap_theory(self.A, force_positive_definite=True, return_eigenvalues=True)
        ddmap_t = ((3. * np.pi) / 8.) * np.power(dmap_t, 2.)
        # compute the ratio between the current value and the target
        compare_ratio = ddmap_t / self.ddmap_target
        # compute the prefactor for iterative scaling
        fhash = np.nansum(ddmap_t) / 2.

        if method == 'IS':
            # compute the gradient
            if lamd > 0.0:
                if reg == 'L2':
                    gradient_t = (np.nan_to_num(
                        np.log(compare_ratio), posinf=0., neginf=0.) - 2. * lamd * self.A) / fhash
                elif reg == 'L1':
                    gradient_t = (np.nan_to_num(
                        np.log(compare_ratio), posinf=0., neginf=0.) + lamd * np.sign(- self.A)) / fhash
            elif lamd == 0.0:
                gradient_t = np.nan_to_num(
                    np.log(compare_ratio), posinf=0., neginf=0.) / fhash

            # enforce symmetry and apply optional off-diagonal mask
            gradient_t = 0.5 * (gradient_t + gradient_t.T)
            if self.edge_mask is not None:
                gradient_t *= self.edge_mask

            # Apply update with optional momentum (standard or Nesterov)
            if momentum > 0.0:
                # Initialize velocity on first iteration
                if t == 0:
                    self.velocity = np.zeros_like(self.A)
                # Update velocity
                self.velocity = momentum * self.velocity + gradient_t
                
                if nesterov:
                    # Nesterov Accelerated Gradient: use look-ahead correction
                    # Update = gradient + momentum * velocity (after velocity update)
                    # This is equivalent to computing gradient at A + momentum * velocity
                    self.A = self.A + learning_rate * (gradient_t + momentum * self.velocity)
                else:
                    # Standard momentum (Polyak's heavy ball)
                    self.A = self.A + learning_rate * self.velocity
            else:
                # Standard update without momentum
                self.A = self.A + learning_rate * gradient_t
        elif method == 'GD':
            if t == 0:
                self.theta = np.copy(self.A)

            # compute the gradient
            if lamd > 0.0:
                if reg == 'L2':
                    gradient_t = (ddmap_t - self.ddmap_target -
                                  2. * lamd * self.A)
                elif reg == 'L1':
                    gradient_t = (ddmap_t - self.ddmap_target +
                                  lamd * np.sign(- self.A))
            elif lamd == 0.0:
                gradient_t = (ddmap_t - self.ddmap_target)

            # perform Nesterov update rule
            # gradient descent state
            theta_previous = np.copy(self.theta)

            #self.theta = self.A + np.maximum(np.minimum(learning_rate * gradient_t, step_cap),-step_cap)
            self.theta = self.A + learning_rate * gradient_t

            # if momentum_rate == None:
            #    momentum_rate = t/(t+3)

            # update the connectivity matrix
            self.A = self.theta + (t/(t+3)) * (self.theta - theta_previous)

        # convert all nan to zero
        self.A = np.nan_to_num(self.A)

        self.A = a2a(self.A, fill_negative=enforce_nonnegative_connectivity_matrix)
        # project to be negative semidefinite
        #self.A = nearestNSD(self.A, 0.0)

        # compute the loss
        self.loss = self.__compute_loss(ddmap_t)
        
        # Compute entropy reusing eigenvalues from a2dmap_theory
        # Note: eigenvalues of K = -A are simply -eigvals_A (no need for second eigendecomposition)
        eigvals_K = -eigvals_A
        self.entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)

    def __update_parameter_gpu(self, t, learning_rate, lamd=0.0, reg='l2', enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        """GPU-accelerated version of __update_parameter for IS method using CuPy."""
        # Compute mean squared distance matrix on GPU
        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        ddmap_t_gpu = ((3. * cp.pi) / 8.) * cp.power(dmap_t_gpu, 2.)
        
        # Compute ratio and gradient on GPU
        compare_ratio_gpu = ddmap_t_gpu / self._ddmap_target_gpu
        fhash = float(cp.nansum(ddmap_t_gpu) / 2.)
        
        # Compute gradient
        if lamd > 0.0:
            if reg == 'L2':
                gradient_t_gpu = (cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) - 2. * lamd * self._A_gpu) / fhash
            elif reg == 'L1':
                gradient_t_gpu = (cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) + lamd * cp.sign(-self._A_gpu)) / fhash
        else:
            gradient_t_gpu = cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) / fhash
        
        # Enforce symmetry
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        
        # Apply optional edge mask (if set, need to transfer to GPU)
        if self.edge_mask is not None:
            edge_mask_gpu = cp.asarray(self.edge_mask)
            gradient_t_gpu *= edge_mask_gpu
        
        # Apply update with optional momentum (standard or Nesterov)
        if momentum > 0.0:
            if t == 0:
                self._velocity_gpu = cp.zeros_like(self._A_gpu)
            self._velocity_gpu = momentum * self._velocity_gpu + gradient_t_gpu
            
            if nesterov:
                self._A_gpu = self._A_gpu + learning_rate * (gradient_t_gpu + momentum * self._velocity_gpu)
            else:
                self._A_gpu = self._A_gpu + learning_rate * self._velocity_gpu
        else:
            self._A_gpu = self._A_gpu + learning_rate * gradient_t_gpu
        
        # Clean up NaN values
        self._A_gpu = cp.nan_to_num(self._A_gpu)
        
        # Keep symmetry
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        
        # Apply edge mask freeze (if set)
        if self.edge_mask is not None:
            freeze_gpu = ~cp.asarray(self.edge_mask)
            cp.fill_diagonal(freeze_gpu, False)
            self._A_gpu = cp.where(freeze_gpu, 0.0, self._A_gpu)
        
        # Rebuild diagonal (a2a operation on GPU)
        # a2a sets diagonal to negative sum of off-diagonal elements in each row
        # and optionally clamps negative off-diagonal values
        if enforce_nonnegative_connectivity_matrix:
            # Clamp off-diagonal to be non-positive (spring constants >= 0)
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)
        
        # Rebuild diagonal: A_ii = -sum(A_ij for j != i)
        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)
        
        # Compute loss on GPU (avoid CPU transfer for ddmap)
        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self.loss = float(loss_gpu)
        
        # Compute entropy on GPU using eigenvalues already on GPU
        # K = -A, so eigenvalues of K = -eigenvalues of A
        # entropy = -sum(log(λ_i)) for positive eigenvalues λ_i of K
        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)
        self.entropy = float(cp.sum(log_terms))
        
        # Sync A back to CPU (needed for saving/display)
        # Note: asnumpy() already synchronizes, no need for explicit sync
        self.A = cp.asnumpy(self._A_gpu)

    def __update_parameter_gpu_gd(self, t, learning_rate, lamd=0.0, reg='l2', enforce_nonnegative_connectivity_matrix=False):
        """GPU-accelerated version of __update_parameter for GD method using CuPy."""
        # Initialize theta on first iteration
        if t == 0:
            self._theta_gpu = self._A_gpu.copy()
        
        # Compute mean squared distance matrix on GPU
        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        ddmap_t_gpu = ((3. * cp.pi) / 8.) * cp.power(dmap_t_gpu, 2.)
        
        # Compute gradient for GD
        if lamd > 0.0:
            if reg == 'L2':
                gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu - 2. * lamd * self._A_gpu
            elif reg == 'L1':
                gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu + lamd * cp.sign(-self._A_gpu)
        else:
            gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu
        
        # Enforce symmetry
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        
        # Apply optional edge mask
        if self.edge_mask is not None:
            edge_mask_gpu = cp.asarray(self.edge_mask)
            gradient_t_gpu *= edge_mask_gpu
        
        # Nesterov-like update (built into GD method)
        theta_previous_gpu = self._theta_gpu.copy()
        self._theta_gpu = self._A_gpu + learning_rate * gradient_t_gpu
        
        # Momentum update: A = theta + (t/(t+3)) * (theta - theta_previous)
        momentum_rate = t / (t + 3)
        self._A_gpu = self._theta_gpu + momentum_rate * (self._theta_gpu - theta_previous_gpu)
        
        # Clean up NaN values
        self._A_gpu = cp.nan_to_num(self._A_gpu)
        
        # Keep symmetry
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        
        # Apply edge mask freeze (if set)
        if self.edge_mask is not None:
            freeze_gpu = ~cp.asarray(self.edge_mask)
            cp.fill_diagonal(freeze_gpu, False)
            self._A_gpu = cp.where(freeze_gpu, 0.0, self._A_gpu)
        
        # Rebuild diagonal (a2a operation on GPU)
        if enforce_nonnegative_connectivity_matrix:
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)
        
        # Rebuild diagonal: A_ii = -sum(A_ij for j != i)
        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)
        
        # Compute loss on GPU (avoid CPU transfer for ddmap)
        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self.loss = float(loss_gpu)
        
        # Compute entropy on GPU using eigenvalues already on GPU
        # K = -A, so eigenvalues of K = -eigenvalues of A
        # entropy = -sum(log(λ_i)) for positive eigenvalues λ_i of K
        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)
        self.entropy = float(cp.sum(log_terms))
        
        # Sync A back to CPU (needed for saving/display)
        # Note: asnumpy() already synchronizes, no need for explicit sync
        self.A = cp.asnumpy(self._A_gpu)

    def run(self, epoch, general_method='optimization', save_steps=None, output_prefix=None, **kwargs):
        """
        Main function to run the optimization
        
        Parameters
        ----------
        epoch : int
            Number of iterations
        general_method : str
            'optimization' or 'direct'
        **kwargs
            Additional arguments passed to __update_parameter:
            - learning_rate : float
                Learning rate for optimization
            - lamd : float
                Regularization weight
            - reg : str
                Regularization type ('L1' or 'L2')
            - method : str
                Optimization method ('IS' or 'GD')
            - enforce_nonnegative_connectivity_matrix : bool
                Enforce non-negative spring constants
            - momentum : float, optional
                Momentum coefficient for IS method (default: 0.0, i.e., no momentum).
                RECOMMENDED: Use 0.95 with nesterov=True for fastest convergence.
                Use 0.9 for conservative settings. Only applies when method='IS'.
            - nesterov : bool, optional
                If True and momentum > 0, use Nesterov Accelerated Gradient (NAG).
                NAG enables higher momentum (0.95) without divergence.
                RECOMMENDED: Use with momentum=0.95 for ~50% faster convergence.
        """

        console = Console()

        loss_array = []

        if general_method == 'optimization':
            with trange(epoch, desc="Performing optimization", unit="iteration") as pbar:
                for t in pbar:
                    self.__update_parameter(t, **kwargs)
                    # display loss at each iterations
                    pbar.set_postfix(loss=self.loss)
                    loss_array.append(self.loss)
        elif general_method == 'direct':
            if not checkEMD(self.ddmap_target):
                raise ValueError(
                    'The distance matrix is a not valid Euclidean distance matrix. Direct inversion method is not applicable. Please use optimization method such as Iterative scaling or gradient descent')
            self.A = ddmap2a_direct(self.ddmap_target)
            # Compute loss for direct inversion
            ddmap_t = ((3. * np.pi) / 8.) * np.power(a2dmap_theory(self.A, force_positive_definite=True), 2.)
            loss_array.append(self.__compute_loss(ddmap_t))

        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)

        return loss_array, dmap_maxent, self.A

class Dynamics:
    def __init__(self, input, M=None, k=None, model=None):
        if isinstance(input, int) and M is None and k is not None:
            if not isinstance(k, int) and not isinstance(k, float):
                sys.stdout.write('Spring constant should be a number')
                sys.exit(0)
            if model != 'rouse' and model is not None:
                sys.stdout.write("Please specify model to be 'rouse'")
                sys.exit(0)

            self.A = construct_connectivity_matrix_rouse(input, k)
            self.eigvalue, self.eigvector = scipy.linalg.eigh(self.A)
            self.N = input
        elif isinstance(input, int) and M is not None and k is not None:
            if not isinstance(k, int) and not isinstance(k, float):
                sys.stdout.write('Spring constant should be a number')
                sys.exit(0)
            if not isinstance(M, int):
                sys.stdout.write('Number of random cross links need to be an integer')
                sys.exit(0)
            if model != 'random':
                sys.stdout.write("Please specify model to be 'random'")
                sys.exit(0)
            self.A = construct_connectivity_matrix_random(input, M, k)
            self.eigvalue, self.eigvector = scipy.linalg.eigh(self.A)
            self.N = input
        elif isinstance(input, np.ndarray) and M is None and k is None:
            if len(input.shape) !=2 or input.shape[0] != input.shape[1]:
                sys.stdout.write('The connectivity matrix should be a square matrix')
                sys.exit(0)
            if not np.allclose(input, input.T):
                sys.stdout.write('The connectivity matrix should be a symmetrix real matrix')
                sys.exit(0)

            self.A = input
            self.eigvalue, self.eigvector = scipy.linalg.eigh(self.A)
            self.N = input.shape[0]

    def generateXYZ(self, force_positive_definite = False):
        self.xyz = a2xyz_sample(self.A, force_positive_definite = force_positive_definite)[0]
        self.modes = self.eigvector.T @ self.xyz

    def initialize(self, dt, zeta, beta):
        if not isinstance(dt, int) and not isinstance(dt, float):
            sys.stdout.write('Time step should be a number')
            sys.exit(0)
        if not isinstance(zeta, int) and not isinstance(zeta, float):
            sys.stdout.write('Friction coefficient step should be a number')
            sys.exit(0)
        if not isinstance(beta, int) and not isinstance(beta, float):
            sys.stdout.write('Temperature step should be a number')
            sys.exit(0)
        elif beta <= 0.0:
            sys.stdout.write('Temperature should be positive')
            sys.exit(0)

        self.zeta = zeta
        self.beta = beta
        self.dt = dt

    def updateModes(self, method='euler-maruyama', force_projection=None, update_zero_modes=True):
        try:
            self.zeta
            self.beta
            self.dt
        except AttributeError:
            sys.stdout.write('Please run initialize() first')
            sys.exit(0)

        # Pass force_projection directly to Ornstein_Uhlenbeck_update
        # The function will handle the conversion to long-term mean internally
        if force_projection is None:
            force_projection = 0.0
            
        self.modes = Ornstein_Uhlenbeck_update(self.modes, self.dt, - self.eigvalue, self.zeta, self.beta, 
                                             force_projection=force_projection, method=method, 
                                             update_zero_modes=update_zero_modes)

    def updateXYZ(self):
        self.xyz = self.eigvector @ self.modes

    def run(self, T, update=1, every=1, initial_conformation=None, method='euler-maruyama', update_zero_modes=True):
        """
        T: number of timesteps
        update: update x,y,z positions every this many timesteps
        every: save x,y,z positions to the trajectory every this many timesteps
        initial_conformation: initial conformation of the simulation
        update_zero_modes: if True, COM will diffuse; if False, COM is fixed at initial position
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)

        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation

        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz)
                #sys.stdout.write('\rTimestep {}'.format(t+1))
                #sys.stdout.flush()
            self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def run_with_force(self, T, force_loci, force_amplitude, force_direction, force_duration=None, 
                      update=1, every=1, initial_conformation=None, method='euler-maruyama', update_zero_modes=True):
        """
        Run dynamics simulation with applied forces following the correct derivation.
        
        The equation of motion is: ξ dR/dt = A*R + B + f
        In eigenmodes: ξ dX/dt = Λ*X + V^T*B + f̃
        
        Parameters
        ----------
        T : int
            Total number of timesteps
        force_loci : list of int
            Indices of loci where force is applied
        force_amplitude : float
            Magnitude of the force
        force_direction : (3,) array_like
            Direction vector of the force (will be normalized)
        force_duration : int, optional
            Number of timesteps to apply force. If None, force is applied for entire simulation.
        update : int, optional
            Update x,y,z positions every this many timesteps
        every : int, optional
            Save x,y,z positions to the trajectory every this many timesteps
        initial_conformation : (N, 3) array_like, optional
            Initial conformation of the simulation
        method : str, optional
            Integration method ('euler-maruyama' or 'exact')
        update_zero_modes : bool, optional
            If True, COM will diffuse; if False, COM is fixed at initial position
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)

        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation

        # Normalize force direction
        force_direction = np.array(force_direction)
        force_direction = force_direction / np.linalg.norm(force_direction)
        
        # Create force vector B (following correct derivation)
        B = np.zeros((self.N, 3))
        for locus in force_loci:
            B[locus] = force_amplitude * force_direction

        # Project force onto eigenmodes: V^T * B
        # V is self.eigvector, so V^T is self.eigvector.T
        # This transforms forces from real space to eigenmode space
        force_projection = self.eigvector.T @ B  # Shape: (n_modes, 3)

        # Allow COM to move naturally due to force
        # No COM constraint applied
        
        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz.copy())
            
            # Apply force only during specified duration
            if force_duration is None or t < force_duration:
                self.updateModes(method=method, force_projection=force_projection, update_zero_modes=update_zero_modes)
            else:
                self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def reset(self):
        self.generateXYZ()

    def run_breakable_bond(self, T, cutoff_distance, update=1, every=1, initial_conformation=None, 
                          method='euler-maruyama', update_zero_modes=True):
        """
        Run dynamics simulation with breakable bonds based on distance cutoff.
        
        IMPORTANT: Nearest-neighbor bonds (i, i+1) representing the polymer backbone are 
        NEVER broken, regardless of distance. Only non-local interactions (i, j where j > i+1) 
        are subject to the distance cutoff.
        
        At each time step:
        1. Check distances between all pairs of loci
        2. For nearest-neighbor bonds (i, i+1): always keep original connectivity value
        3. For non-local bonds (i, j where j > i+1):
           - If distance < cutoff: bond is present (use original connectivity value)
           - If distance >= cutoff: bond is broken (set connectivity to zero)
        4. Update connectivity matrix and recompute eigenvalues/eigenvectors
        5. Update normal modes using the new connectivity matrix
        
        Parameters
        ----------
        T : int
            Total number of timesteps
        cutoff_distance : float
            Distance threshold for breaking non-local bonds. Only applies to bonds 
            between non-adjacent loci (i, j where j > i+1).
        update : int, optional
            Update x,y,z positions every this many timesteps
        every : int, optional
            Save x,y,z positions to the trajectory every this many timesteps
        initial_conformation : (N, 3) array_like, optional
            Initial conformation of the simulation
        method : str, optional
            Integration method ('euler-maruyama' or 'exact'). Default is 'euler-maruyama'.
            NOTE: 'euler-maruyama' is STRONGLY RECOMMENDED for breakable bond simulations,
            as broken bonds can create positive eigenvalues (unstable modes). The euler-maruyama
            method handles these gracefully with controlled growth over small timesteps, while
            the 'exact' method may produce numerical issues.
        update_zero_modes : bool, optional
            If True, COM will diffuse; if False, COM is fixed at initial position
        
        Notes
        -----
        The input connectivity matrix is preserved and not modified. All calculations 
        use an internal copy (A_internal). Nearest-neighbor bonds are never broken 
        to maintain polymer chain connectivity.
        
        **Handling Positive Eigenvalues:**
        When bonds break, the connectivity matrix may develop positive eigenvalues (unstable modes).
        This causes exponential growth in those modes. However, the feedback mechanism (bonds
        break as distances increase) combined with small timesteps prevents catastrophic blowup.
        The system remains stable because: (1) growth is incremental with small dt, and (2) as
        distances grow, more bonds break, reducing the instability.
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)
        if not isinstance(cutoff_distance, (int, float)) or cutoff_distance <= 0:
            sys.stdout.write('Cutoff distance should be a positive number')
            sys.exit(0)

        # Store original connectivity matrix for reference (keep input unchanged)
        self.A_original = np.copy(self.A)
        
        # Create internal connectivity matrix for breakable bond calculations
        self.A_internal = np.copy(self.A)
        
        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation
            self.modes = self.eigvector.T @ self.xyz

        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz.copy())
            
            # Update internal connectivity matrix based on current distances
            self._update_connectivity_from_distances(cutoff_distance)
            
            # Update normal modes (eigvalue and eigvector already updated in _update_connectivity_from_distances)
            self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def _update_connectivity_from_distances(self, cutoff_distance):
        """
        Update internal connectivity matrix based on current pairwise distances.
        
        Nearest-neighbor bonds (i, i+1) are never broken - they represent the polymer backbone.
        Only non-local interactions are subject to the distance cutoff.
        
        Parameters
        ----------
        cutoff_distance : float
            Distance threshold for bond breaking (only applies to non-nearest-neighbor bonds)
        """
        # Calculate pairwise distances
        distances = self._calculate_pairwise_distances()
        
        # Update internal connectivity matrix
        # Only loop over non-local bonds (j > i+1) since nearest-neighbor bonds never break
        for i in range(self.N):
            for j in range(i+2, self.N):
                # Non-local bonds: apply distance cutoff
                if distances[i, j] < cutoff_distance:
                    # Bond is present - use original connectivity value
                    self.A_internal[i, j] = self.A_original[i, j]
                    self.A_internal[j, i] = self.A_original[j, i]
                else:
                    # Bond is broken - set to zero
                    self.A_internal[i, j] = 0.0
                    self.A_internal[j, i] = 0.0
        
        # Update diagonal elements to maintain Laplacian property using a2a function
        self.A_internal = a2a(self.A_internal, fill_negative=False)
        
        # Recompute eigenvalues and eigenvectors using internal matrix
        self.eigvalue, self.eigvector = scipy.linalg.eigh(self.A_internal)
        
        # Update modes to maintain consistency
        self.modes = self.eigvector.T @ self.xyz

    def _calculate_pairwise_distances(self):
        """
        Calculate pairwise distances between all loci using scipy.
        
        Returns
        -------
        distances : (N, N) ndarray
            Matrix of pairwise distances. distances[i,j] is the distance between loci i and j.
        """
        # Use pdist+squareform for better performance on larger systems
        # pdist computes condensed distance matrix, squareform converts to full symmetric matrix
        distances = scipy.spatial.distance.squareform(
            scipy.spatial.distance.pdist(self.xyz, metric='euclidean')
        )
        return distances

    def restore_original_connectivity(self):
        """
        Restore the original connectivity matrix and recompute eigenvalues/eigenvectors.
        This is useful after running breakable bond simulations to return to the original state.
        """
        if hasattr(self, 'A_original'):
            self.A = np.copy(self.A_original)
            self.eigvalue, self.eigvector = scipy.linalg.eigh(self.A)
            # Update modes to maintain consistency
            self.modes = self.eigvector.T @ self.xyz
        else:
            sys.stdout.write('No original connectivity matrix found. Cannot restore.')


def run_optimization(input_path=None,
                     output_prefix=None,
                     input_matrix=None,
                     connectivity_matrix=None,
                     ensemble=1000,
                     alpha=4.0,
                     selection=None,
                     method='IS',
                     lamd=0.0,
                     reg='L2',
                     iteration=10000,
                     learning_rate=10.0,
                     momentum=0.0,
                     nesterov=False,
                     use_gpu=False,
                     input_type='cmap',
                     input_format='text',
                     binsize=25000,
                     hic_norm='KR',
                     hic_unit='BP',
                     log=False,
                     no_xyzs=False,
                     ignore_missing_data=False,
                     balance=False,
                     not_normalize=False,
                     neighbor_balance=False,
                     enforce_nonnegative_connectivity_matrix=False,
                     verbose=True):
    """
    Core function to run HIPPS/DIMES optimization that can be called programmatically or from CLI.
    
    Parameters
    ----------
    input_path : str, optional
        Path to the input file (required if input_matrix is not provided)
    output_prefix : str, optional
        Prefix for output files (if None, results are only returned, not saved)
    input_matrix : np.ndarray, optional
        Input matrix (contact map or distance map). If provided, input_path is ignored
    connectivity_matrix : np.ndarray or str, optional
        Initial connectivity matrix or path to file containing it
    ensemble : int, default=1000
        Number of conformations to generate
    alpha : float, default=4.0
        Exponent for cmap-to-dmap conversion
    selection : str, optional
        Region selection for cooler/hic files
    method : str, default='IS'
        Optimization method: 'IS' (Iterative Scaling), 'GD' (Gradient Descent), or 'DI' (Direct Inversion)
    lamd : float, default=0.0
        Regularization weight
    reg : str, default='L2'
        Regularization type: 'L1' or 'L2'
    iteration : int, default=10000
        Number of optimization iterations
    learning_rate : float, default=10.0
        Learning rate for optimization
    momentum : float, default=0.0
        Momentum coefficient for IS method.
        RECOMMENDED: Use 0.95 with nesterov=True for fastest convergence (~50% faster).
        Use 0.9 for more conservative settings. Only applies when method='IS'.
    nesterov : bool, default=False
        If True and momentum > 0, use Nesterov Accelerated Gradient (NAG).
        NAG's look-ahead correction enables higher momentum (0.95) without divergence.
        RECOMMENDED: Use with momentum=0.95 for best performance.
    use_gpu : bool, default=False
        If True and CuPy is installed, use GPU acceleration for eigendecomposition.
        Provides 40-180x speedup for matrices with n >= 200.
        Requires: conda install -c conda-forge cupy
    input_type : str, default='cmap'
        Type of input: 'cmap' (contact map) or 'dmap' (distance map)
    input_format : str, default='text'
        Format of input file: 'text', 'cooler', or 'hic'
    binsize : int, default=25000
        Bin size for .hic format in bp
    hic_norm : str, default='KR'
        Normalization for .hic: 'KR', 'VC', 'NONE'
    hic_unit : str, default='BP'
        Unit for .hic: 'BP' or 'FRAG'
    log : bool, default=False
        Whether to write log file
    no_xyzs : bool, default=False
        If True, skip writing conformations to file
    ignore_missing_data : bool, default=False
        Whether to ignore missing elements in contact/distance map
    balance : bool, default=False
        Apply matrix balancing for contact map (cooler format)
    not_normalize : bool, default=False
        Turn off auto normalization of contact map
    neighbor_balance : bool, default=False
        Apply neighbor balancing for contact map
    enforce_nonnegative_connectivity_matrix : bool, default=False
        Enforce non-negative spring constants
    verbose : bool, default=True
        Whether to print status messages
        
    Returns
    -------
    results : dict
        Dictionary containing:
        - 'loss': Loss values during optimization (pandas DataFrame or list)
        - 'dmap_final': Final distance map (numpy array)
        - 'connectivity_matrix': Final connectivity matrix (numpy array)
        - 'cmap_final': Final contact map (numpy array, only if input_type=='cmap')
        - 'xyzs': Generated conformations (numpy array, only if no_xyzs==False)
        - 'rc_optimal': Optimal contact threshold (float, only if input_type=='cmap')
    
    Notes
    -----
    **Convergence Acceleration**:
    
    **Nesterov Momentum** (RECOMMENDED for fastest convergence):
    - Use momentum=0.95 with nesterov=True for best performance
    - ~50% faster than standard momentum at 0.9
    - Nesterov's look-ahead correction enables higher momentum without divergence
    - Example: momentum=0.95, nesterov=True
    
    **Standard Momentum** (fallback):
    - Use momentum=0.9 if you prefer more conservative settings
    - Example: momentum=0.9
    
    **GPU Acceleration** (for large matrices):
    - Use use_gpu=True when CuPy is installed
    - Provides 40-180x speedup for matrices with n >= 200
    - For n < 200, CPU may be faster due to GPU transfer overhead
    - Install CuPy: conda install -c conda-forge cupy
    
    Examples
    --------
    >>> # Basic usage with numpy array
    >>> cmap = np.loadtxt('contact_map.txt')
    >>> results = run_optimization(input_matrix=cmap, input_type='cmap', 
    ...                            method='IS', iteration=5000)
    >>> connectivity_matrix = results['connectivity_matrix']
    >>> xyzs = results['xyzs']
    
    >>> # With momentum for faster convergence (recommended)
    >>> results = run_optimization(input_matrix=cmap, input_type='cmap',
    ...                            method='IS', iteration=5000,
    ...                            learning_rate=10.0, momentum=0.9)
    
    >>> # Use as a library with a cooler file
    >>> results = run_optimization(input_path='data.cool', 
    ...                            input_type='cmap',
    ...                            input_format='cooler',
    ...                            selection='chr21',
    ...                            output_prefix='output/chr21',
    ...                            momentum=0.9)
    """
    # Validate inputs
    if input_matrix is None and input_path is None:
        raise ValueError("Either input_matrix or input_path must be provided")
    
    # Initialize console for output
    if verbose:
        console = Console()
        title = Text.assemble(("HIPPS-DIMES", "bold yellow"),
                              ": Maximum Entropy Based HI-C/Distance Map - Polymer Physics - Structures Reconstruction - Dynamics Prediction\n",
                              "1. Shi, Guang, and D. Thirumalai. From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes. Physical Review X 11.1 (2021): 011051.\n",
                              "2. Shi, Guang, and D. Thirumalai. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nature Communications 14.1 (2023): 1150.\n",
                              "3. Shi, Guang, Shin, Sucheol, and D. Thirumalai. Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes. Science Advances 11.31 (2025): eadx1763.")
        console.print(Panel(title))
    else:
        console = None

    # Load or use input matrix
    if verbose and console:
        status = console.status("[bold green]System initialization...")
        status.start()
    
    try:
        if input_type == 'dmap':
            if verbose and console:
                console.print("Reading distance matrix")
            if input_matrix is not None:
                # Use provided matrix
                dmap_target = input_matrix
                if dmap_target.ndim != 2 or dmap_target.shape[0] != dmap_target.shape[1]:
                    raise ValueError("Distance map must be a square 2D array")
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
            elif input_format == 'text':
                dmap_target = np.loadtxt(input_path)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
            elif input_format == 'cooler':
                raise ValueError('input-type=dmap only supports text format file or numpy array')
        elif input_type == 'cmap':
            if verbose and console:
                console.print("Reading contact map")
            if input_matrix is not None:
                # Use provided matrix
                cmap = input_matrix
                if cmap.ndim != 2 or cmap.shape[0] != cmap.shape[1]:
                    raise ValueError("Contact map must be a square 2D array")
            elif input_format == 'text':
                cmap = np.loadtxt(input_path)
            elif input_format == 'cooler':
                cmap_data = cooler.Cooler(input_path)
                if verbose and console:
                    console.print("Cooler file read completed")
                cmap = cmap_data.matrix(balance=balance).fetch(selection)
                if verbose and console:
                    console.print("Cooler file selection completed")
                if len(cmap) >= 5000:
                    warning_msg = "The matrix size is {}x{}. It is too large. Please use smaller matrix".format(
                        len(cmap), len(cmap))
                    if verbose and console:
                        console.print(warning_msg)
                    raise ValueError(warning_msg)
            elif input_format == 'hic':
                # parse selection for .hic
                if not selection or ',' not in selection:
                    raise ValueError("For .hic input, --selection must be 'chr1:start1-end1,chr2:start2-end2'")
                hic = hicstraw.HiCFile(input_path)
                if verbose and console:
                    console.print(".hic format file read completed")
                reg1, reg2 = selection.split(',')
                # strip optional 'chr' prefix for hicstraw
                raw_chrom1, r1 = reg1.split(':')
                raw_chrom2, r2 = reg2.split(':')
                chrom1 = raw_chrom1[3:] if raw_chrom1.lower().startswith('chr') else raw_chrom1
                chrom2 = raw_chrom2[3:] if raw_chrom2.lower().startswith('chr') else raw_chrom2
                start1, end1 = map(int, r1.split('-'))
                start2, end2 = map(int, r2.split('-'))
                
                # try efficient random-access API
                matrix_obj = hic.getMatrixZoomData(chrom1, chrom2, 'observed', hic_norm, hic_unit, binsize)
                if verbose and console:
                    console.print("Fetched hic matrix zoom data")
                try:
                    cmap = matrix_obj.getRecordsAsMatrix(start1, end1, start2, end2)
                except Exception:
                    # fallback manual assembly via straw
                    if verbose and console:
                        console.print("Falling back to manual assembly via hicstraw.straw()...")
                    region1 = f"{chrom1}:{start1}:{end1}"
                    region2 = f"{chrom2}:{start2}:{end2}"
                    result = hicstraw.straw('observed', hic_norm, input_path,
                                             region1, region2,
                                             hic_unit, binsize)
                    # compute dimensions
                    dim1 = (end1 - start1) // binsize + 1
                    dim2 = (end2 - start2) // binsize + 1
                    cmap = np.zeros((dim1, dim2))
                    # build map
                    desc = "Building hic contact map" if verbose else None
                    for pt in tqdm(result, desc=desc, disable=not verbose):
                        i = int((pt.binX - start1) / binsize)
                        j = int((pt.binY - start2) / binsize)
                        cmap[i, j] = pt.counts
                    cmap = cmap + cmap.T

                if verbose and console:
                    console.print(".hic contact map extracted")
            
            # Apply neighbor balancing if requested
            if neighbor_balance:
                if verbose and console:
                    console.print("Applying neighbor balancing to contact map (see Paggi, Zhang 2025)")
                cmap = neighbor_balance_symmetric(cmap, not_normalize=not_normalize)
            
            if ignore_missing_data:
                dmap_target = cmap2dmap_missing_data(cmap, alpha, not_normalize)
            else:
                dmap_target = cmap2dmap(cmap, alpha, not_normalize)
            dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        else:
            raise ValueError("Invalid input_type. Must be 'cmap' or 'dmap'")

        # Load connectivity matrix if provided
        if connectivity_matrix is not None:
            if isinstance(connectivity_matrix, str):
                connectivity_matrix = np.loadtxt(connectivity_matrix)
            if verbose and console:
                console.print("Loaded the provided connectivity matrix and will use it as initialization.")
        
        if verbose and console:
            console.print("Initialization completed")
            if hasattr(status, 'stop'):
                status.stop()
    except Exception as e:
        if verbose and console and hasattr(status, 'stop'):
            status.stop()
        raise e

    # Print parameter table if verbose
    if verbose and console:
        table = Table(title="Some Basic Parameters")
        table.add_column("Input Source", no_wrap=False)
        table.add_column("Input Type", no_wrap=False)
        table.add_column("Input Format", no_wrap=False)
        table.add_column("Optimization method", no_wrap=False)
        table.add_column("Matrix Size", no_wrap=False)
        table.add_column("Number of Iterations", no_wrap=False)
        table.add_column("Regularization", no_wrap=False)
        table.add_column("Ignore Missing Data", no_wrap=False)
        table.add_column("Matrix Balancing", no_wrap=False)
        table.add_column("Neighbor Balancing", no_wrap=False)
        table.add_column("Matrix Normalization", no_wrap=False)
        
        input_source = input_path if input_path else "NumPy array"
        table.add_row(input_source,
                      "{}".format("Contact Map" if input_type ==
                                  'cmap' else "Distance Map" if input_type == 'dmap' else "Unknown"),
                      "{}".format("Text" if input_format ==
                                  'text' else "Cooler File" if input_format == 'cooler' else ".hic file" if input_format == 'hic' else "Unknown"),
                      "{}".format("Iterative Scaling" if method == 'IS' else "Gradient Descent" if method ==
                                  'GD' else "Direct Inversion" if method == 'DI' else "Unknown"),
                      "{}".format("{}x{}".format(
                          dmap_target.shape[0], dmap_target.shape[1])),
                      "{}".format(iteration),
                      "{}".format(reg if lamd > 0.0 else "No" if lamd ==
                                  0.0 else "Unknown"),
                      "{}".format("Yes" if ignore_missing_data else "No"),
                      "{}".format("Yes" if balance else "No" if (
                          balance is False and input_format == 'cooler') else "N/A"),
                      "{}".format("Yes" if neighbor_balance else "No" if (
                          neighbor_balance is False and input_type == 'cmap') else "N/A"),
                      "{}".format("No" if (not_normalize is True and input_type == 'cmap') else "Yes" if (
                          not_normalize is False and input_type == 'cmap') else "N/A")
                      )
        console.print(table)

    # GPU availability reminder
    if verbose:
        if is_gpu_available() and not use_gpu:
            gpu_name = get_gpu_name()
            console.print(f"[cyan]Tip: GPU detected ({gpu_name}). Use --use-gpu (CLI) or use_gpu=True for 40-180x speedup on large matrices (n >= 200).[/cyan]")
        elif not is_gpu_available() and dmap_target.shape[0] >= 200:
            console.print("[cyan]Tip: For large matrices, GPU acceleration can provide 40-180x speedup. Install CuPy to enable: conda install -c conda-forge cupy[/cyan]")
    
    # Run optimization
    model = Optimize(dmap_target, connectivity_matrix=connectivity_matrix, use_gpu=use_gpu)
    
    if use_gpu and model.use_gpu and verbose:
        console.print(f"[green]GPU acceleration enabled ({get_gpu_name()})[/green]")
    
    keyword_arguments = {'learning_rate': learning_rate, 'lamd': lamd, 'reg': reg, 'method': method,
                         'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix,
                         'momentum': momentum, 'nesterov': nesterov}

    if method == 'IS' or method == 'GD':
        general_method = 'optimization'
    elif method == 'DI':
        general_method = 'direct'

    loss, dmap_maxent, final_connectivity_matrix = model.run(
        iteration, general_method=general_method, **keyword_arguments)
    
    # Format loss data
    try:
        loss_df = pd.DataFrame(
            np.dstack((np.arange(1, len(loss)+1), loss))[0], columns=['iteration', 'loss'])
    except IndexError:
        loss_df = None

    # Print regularization norms if requested
    if verbose:
        if reg == 'L2':
            print('L2 norm of the connectivity matrix:', np.linalg.norm(
                final_connectivity_matrix[np.triu_indices_from(final_connectivity_matrix, k=1)]))
        elif reg == 'L1':
            print('L1 norm of the connectivity matrix:', np.abs(
                final_connectivity_matrix[np.triu_indices_from(final_connectivity_matrix, k=1)]).sum())

        if loss_df is not None and console:
            console.print("Final loss: {}".format(loss_df['loss'].values[-1]))

    # Finalize results
    results = {
        'loss': loss_df if loss_df is not None else loss,
        'dmap_final': dmap_maxent,
        'connectivity_matrix': final_connectivity_matrix
    }
    
    # Compute contact map if input was contact map
    cmap_maxent = None
    rc_optimal = None
    if input_type == 'cmap':
        cmap_rc_minimize_res = scipy.optimize.minimize_scalar(
            objective_func, args=(final_connectivity_matrix, cmap))
        rc_optimal = cmap_rc_minimize_res.x
        if verbose and console:
            console.print('Optimized contact threshold distance: {}\n'.format(rc_optimal))
        cmap_maxent = a2cmap_theory(final_connectivity_matrix, rc_optimal)
        results['cmap_final'] = cmap_maxent
        results['rc_optimal'] = rc_optimal

    # Generate conformations if requested
    xyzs = None
    if not no_xyzs:
        xyzs = a2xyz_sample(final_connectivity_matrix, ensemble=ensemble)
        results['xyzs'] = xyzs

    # Save output files if output_prefix is provided
    if output_prefix is not None:
        if verbose and console:
            status = console.status("[bold green]Saving results...")
            status.start()
        
        try:
            if log and loss_df is not None:
                loss_df.to_csv('{}_loss_function_iteration.csv'.format(output_prefix))
                if verbose and console:
                    console.print(
                        "Loss function saved to file: [bold magenta]{}_loss_function_iteration.csv[/bold magenta]".format(output_prefix))

            np.savetxt('{}_dmap_final.txt'.format(output_prefix), dmap_maxent)
            if verbose and console:
                console.print(
                    "Final distance map saved to file: [bold magenta]{}_dmap_final.txt[/bold magenta]".format(output_prefix))
            
            if input_type == 'cmap':
                np.savetxt('{}_dmap_target.txt'.format(output_prefix), np.sqrt((8./(3.*np.pi))* dmap_target))
                if verbose and console:
                    console.print(
                        "Target distance map saved to file: [bold magenta]{}_dmap_target.txt[/bold magenta]".format(output_prefix))
                np.savetxt('{}_cmap_final.txt'.format(output_prefix), cmap_maxent)
                if verbose and console:
                    console.print(
                        "Final contact map saved to file: [bold magenta]{}_cmap_final.txt[/bold magenta]".format(output_prefix))
            
            np.savetxt('{}_connectivity_matrix.txt'.format(output_prefix), final_connectivity_matrix)
            if verbose and console:
                console.print(
                    'Connectivity matrix saved to file: [bold magenta]{}_connectivity_matrix.txt[/bold magenta]'.format(output_prefix))

            if not no_xyzs and xyzs is not None:
                write2xyz('{}.xyz'.format(output_prefix), xyzs)
                if verbose and console:
                    console.print(
                        "Ensemble of structures saved to file: [bold magenta]{}.xyz[/bold magenta]".format(output_prefix))
            
            if verbose and console:
                status.stop()
        except Exception as e:
            if verbose and console and hasattr(status, 'stop'):
                status.stop()
            raise e
    
    return results


@click.command()
@click.argument('input', nargs=1)
@click.argument('output-prefix', nargs=1)
@click.option('-k', '--connectivity-matrix', type=str, required=False, help='Use provided connectivity matrix as initialization. Useful when restart from previous run')
@click.option('-e', '--ensemble', type=int, default=1000, show_default=True, help='specify the number of conformations generated')
@click.option('-a', '--alpha', type=float, default=4.0, show_default=True, help='specify the value of cmap-to-dmap conversion exponent')
@click.option('-s', '--selection', type=str, required=False, help='For cooler: any valid selector for cooler.Cooler.matrix().fetch(), e.g. "chr1" or "chr1::start-end". For .hic: use "chr1:start1-end1,chr2:start2-end2"')
@click.option('-m', '--method', type=click.Choice(['IS', 'GD', 'DI'], case_sensitive=True), default='IS', show_default=True, help='specify the method. IS: Iterative Scaling. GD: Gradient Descent. DI: Direct Inversion. When using\
    Direct Inversion, no iterations are performed. The connectivity matrix is obtained by direct Moore–Penrose inverse of the covariance matrix. Note that the resulting connectivity matrix using Direct Inversion can be very different from the results obtained by GD or IS method.')
@click.option('-l', '--lamd', type=click.FloatRange(0, max=None), default=0.0, show_default=True, help='Specify the weight for the regularization.')
@click.option('-r', '--reg', type=click.Choice(['L1', 'L2'], case_sensitive=True), default='L2', show_default=True, required=False, help='specify the type of regularization. Currently support L1 and L2 regularization. Note that this option should be used together with option -l')
@click.option('-i', '--iteration', type=int, default=10000, show_default=True, help='Number of iterations')
@click.option('-r', '--learning-rate', type=float, default=10.0, show_default=True, help='Learning rate. This hyperparameter controls the speed of convergence. \
    If its value is too small, then convergence is very slow. If its value is too large, the program may never converge. Typically, learning rate can be set to be 1-30 if use Iterative scaling method. \
        It should be a very small value (such as 1e-8) when using gradient descent optimization')
@click.option('--momentum', type=click.FloatRange(0, 1), default=0.0, show_default=True, help='Momentum coefficient for IS method. \
    RECOMMENDED: Use 0.95 with --nesterov for fastest convergence (~50%% faster). Use 0.9 for conservative settings. Only applies when method=IS.')
@click.option('--nesterov', is_flag=True, default=False, show_default=True, help='Use Nesterov Accelerated Gradient (NAG). \
    Enables higher momentum (0.95) without divergence. RECOMMENDED: Use with --momentum 0.95 for fastest convergence.')
@click.option('--use-gpu', is_flag=True, default=False, show_default=True, help='Use GPU acceleration via CuPy. \
    Provides 40-180x speedup for large matrices (n >= 200). Requires CuPy: conda install -c conda-forge cupy')
@click.option('--input-type', required=True, type=click.Choice(['cmap', 'dmap'], case_sensitive=False), help='Specify the type of the input. cmap: contact map or dmap: distance map')
@click.option('--input-format', required=True, type=click.Choice(['text', 'cooler', 'hic'], case_sensitive=False), help='Format of input: text, cooler, or hic')
@click.option('--binsize', type=int, default=25000, show_default=True, help='Bin size (resolution) for .hic format in bp')
@click.option('--norm', 'hic_norm', type=str, default='KR', show_default=True, help='Normalization for .hic: KR, VC, NONE')
@click.option('--unit', 'hic_unit', type=click.Choice(['BP', 'FRAG'], case_sensitive=False), default='BP', show_default=True, help='Unit for .hic: BP or FRAG')
@click.option('--log', is_flag=True, default=False, show_default=True, help='Write a log file')
@click.option('--no-xyzs', is_flag=True, default=False, show_default=True, help='Turn off writing conformations to .xyz file')
@click.option('--ignore-missing-data', is_flag=True, default=False, show_default=True, help='Turn on this argument will let the program ignore the missing elementsin the contact map or distance map')
@click.option('--balance', is_flag=True, default=False, show_default=True, help='Turn on the matrix balance for contact map. Only effective when input_type == cmap and input_format == cooler')
@click.option('--neighbor-balance', is_flag=True, default=False, show_default=True, help='Turn on neighbor balancing for contact map. Only effective when input_type == cmap. Normalizes contact between i and j by dividing it by the geometric mean of neighbor contact for i and j. see Paggi, Zhang 2025 for method details')
@click.option('--not-normalize', is_flag=True, default=False, show_default=True, help='Turn off auto normalization of contact map. Only effective when the input is contact map')
@click.option('--enforce-nonnegative-connectivity-matrix', is_flag=True, default=False, show_default=True, help='Enforcing that the "spring constants" in the connectivity matrix can only be nonnegative')
@click.option('--save-steps', type=str, default=None, help='Comma-separated list of iteration steps at which to save the connectivity matrix. Example: --save-steps 1000,5000,10000,50000')
def main(input, output_prefix, connectivity_matrix, ensemble, alpha, selection, method, lamd, reg, iteration, learning_rate, momentum, nesterov, use_gpu, input_type, \
    input_format, binsize, hic_norm, hic_unit, log, no_xyzs, ignore_missing_data, balance, not_normalize, neighbor_balance, enforce_nonnegative_connectivity_matrix, save_steps):
    """
    Command-line interface for HIPPS/DIMES to generate ensemble of genome structures from either contact map or mean distance map.
    
    INPUT: Specify the path to the input file
    
    OUTPUT_PREFIX: Specify the prefix for output files
    
    If you use this program in your publication, please cite this paper: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051
    """
    # Call the core function
    run_optimization(
        input_path=input,
        output_prefix=output_prefix,
        input_matrix=None,
        connectivity_matrix=connectivity_matrix,
        ensemble=ensemble,
        alpha=alpha,
        selection=selection,
        method=method,
        lamd=lamd,
        reg=reg,
        iteration=iteration,
        learning_rate=learning_rate,
        momentum=momentum,
        nesterov=nesterov,
        use_gpu=use_gpu,
        input_type=input_type,
        input_format=input_format,
        binsize=binsize,
        hic_norm=hic_norm,
        hic_unit=hic_unit,
        log=log,
        no_xyzs=no_xyzs,
        ignore_missing_data=ignore_missing_data,
        balance=balance,
        not_normalize=not_normalize,
        neighbor_balance=neighbor_balance,
        enforce_nonnegative_connectivity_matrix=enforce_nonnegative_connectivity_matrix,
        verbose=True
    )


if __name__ == '__main__':
    main()
