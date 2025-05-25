"""
Reconstruction of 3D genome organization using the Maximum Entropy Principle

Reference:
1. Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051.
https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051
2. Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants." Nature Communications 14.1 (2023): 1150.
"""

import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import numpy as np
import scipy
import scipy.linalg
import scipy.interpolate
import scipy.optimize
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

def compute_m1_all(a, t, zeta=1.0):
    """
    Compute the single-monomer mean-square displacement (MSD) for all monomers, 
    given the connectivity matrix `a`. This is a vectorized version that computes
    MSD for all monomers simultaneously.

    Parameters
    ----------
    a : np.ndarray
        The connectivity (or "Laplacian") matrix for the polymer/chain.
    t : array_like
        A 1D array of time points (lag times).
    zeta : float
        Friction coefficient. Default is 1.0.

    Returns
    -------
    msd_all : np.ndarray
        3D array of shape (n, len(t), 2) where:
        - First dimension is the monomer index
        - Second dimension is the time points
        - Third dimension contains [time, msd] pairs
    """
    n = a.shape[0]
    eigvalue, eigvector = scipy.linalg.eigh(a)
    eigvalue_inv = 1.0 / eigvalue

    # Filter out infinities
    normal_modes_square_mean = - np.nan_to_num(eigvalue_inv, posinf=0.0, neginf=0.0)

    # Expand time dimension for broadcast
    t_reshaped = np.expand_dims(t, axis=-1)
    tau_p = - zeta / eigvalue
    decay_factor = np.exp(-t_reshaped / tau_p)

    # Compute for all monomers simultaneously
    # Shape: (n, len(t), n-1) where n-1 is number of modes
    vpi_squared = np.power(eigvector, 2)  # Shape: (n, n-1)
    
    # Time-dependent part for all monomers
    # Shape: (n, len(t))
    res = 3.0 * np.sum(vpi_squared[:, None, :] * decay_factor[None, :, :] * normal_modes_square_mean[None, None, :], axis=-1)
    
    # Equilibrium radius for all monomers
    # Shape: (n,)
    r2_eq = 3.0 * np.sum(vpi_squared * normal_modes_square_mean[None, :], axis=-1)
    
    # MSD for all monomers
    # Shape: (n, len(t))
    msd_data = 2.0 * (r2_eq[:, None] - res)

    # Combine time with MSD for all monomers
    # Shape: (n, len(t), 2)
    msd_all = np.stack([np.tile(t, (n, 1)), msd_data], axis=-1)
    
    return msd_all

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


def Ornstein_Uhlenbeck_update(x, dt, k, zeta, beta, b = 0.0, method='euler-maruyama'):
    """
    Update variable x for a Ornstein Uhlenbeck process
    x: Array for value of x of each degree of freedom
    k: Array for spring constant for each degree of freedom
    zeta: one value
    beta: one value
    """
    if isinstance(x, np.ndarray):
        rand_noise = np.random.randn(*x.shape)
    else:
        rand_noise = np.random.randn()
    
    if method == 'euler-maruyama':
        dx = - k[:, np.newaxis] * x * dt / zeta + b * dt / zeta + np.sqrt(2.0 * dt / (zeta * beta)) * rand_noise
        x_new = x + dx
    elif method == 'exact':
        theta = k[:, np.newaxis] / zeta
        sigma = (2. / (zeta * beta)) ** .5
        mu = np.exp(- theta * dt)
        x_new = x * mu + np.nan_to_num(np.sqrt((sigma ** 2. / (2. * theta)) * (1. - mu ** 2.))) * rand_noise
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


def a2dmap_theory(A, force_positive_definite=False):
    """
    Return mean distance map given the connectivity matrix A theoretically
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
    np.fill_diagonal(temp, - np.sum(a, axis=1) + a.diagonal())
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

def a2xyz_sample_fixed_end(A,
                           xyz_start,
                           xyz_end,
                           ensemble=1,
                           force_positive_definite=False):
    """
    Generate `ensemble` random polymer configurations from connectivity matrix A,
    *with* bead 0 fixed at xyz_start and bead n−1 fixed at xyz_end.

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

#------------------------------------------------------------------#


class Optimize:
    def __init__(self, ddmap_target, connectivity_matrix=None):
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


    def __compute_loss(self, ddmap_t):
        with np.errstate(divide='ignore', invalid='ignore'):
            loss = np.nanmean(
                np.power((ddmap_t - self.ddmap_target)/self.ddmap_target, 2.)) ** .5
        return loss

    def __update_parameter(self, t, learning_rate, lamd=0.0, reg='l2', method='IS', enforce_nonnegative_connectivity_matrix=False):
        # updating using Iterative Scaling

        # compute the mean squared distance matrix at current iteration step
        ddmap_t = ((3. * np.pi) / 8.) * \
            np.power(a2dmap_theory(self.A, force_positive_definite=True), 2.)
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

            # update the connectivity matrix
            self.A += learning_rate * gradient_t
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

    def run(self, epoch, general_method='optimization', **kwargs):
        """
        Main function to run the optimization
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
            loss_array.append(self.__compute_loss())

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

    def updateModes(self, method='euler-maruyama'):
        try:
            self.zeta
            self.beta
            self.dt
        except AttributeError:
            sys.stdout.write('Please run initialize() first')
            sys.exit(0)

        self.modes = Ornstein_Uhlenbeck_update(self.modes, self.dt, - self.eigvalue, self.zeta, self.beta, method=method)
        # self.modes = OU.OU(self.modes, self.dt, - self.eigvalue, self.zeta, self.beta)

    def updateXYZ(self):
        self.xyz = self.eigvector @ self.modes

    def run(self, T, update=1, every=1, initial_conformation=None, method='euler-maruyama'):
        """
        T: number of timesteps
        update: update x,y,z positions every this many timesteps
        every: save x,y,z positions to the trajectory every this many timesteps
        initial_conformation: initial conformation of the simulation
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
            self.updateModes(method=method)

        self.traj = np.array(self.traj)

    def reset(self):
        self.generateXYZ()


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
@click.option('--input-type', required=True, type=click.Choice(['cmap', 'dmap'], case_sensitive=False), help='Specify the type of the input. cmap: contact map or dmap: distance map')
@click.option('--input-format', required=True, type=click.Choice(['text', 'cooler', 'hic'], case_sensitive=False), help='Format of input: text, cooler, or hic')
@click.option('--binsize', type=int, default=25000, show_default=True, help='Bin size (resolution) for .hic format in bp')
@click.option('--norm', 'hic_norm', type=str, default='KR', show_default=True, help='Normalization for .hic: KR, VC, NONE')
@click.option('--unit', 'hic_unit', type=click.Choice(['BP', 'FRAG'], case_sensitive=False), default='BP', show_default=True, help='Unit for .hic: BP or FRAG')
@click.option('--log', is_flag=True, default=False, show_default=True, help='Write a log file')
@click.option('--no-xyzs', is_flag=True, default=False, show_default=True, help='Turn off writing conformations to .xyz file')
@click.option('--ignore-missing-data', is_flag=True, default=False, show_default=True, help='Turn on this argument will let the program ignore the missing elementsin the contact map or distance map')
@click.option('--balance', is_flag=True, default=False, show_default=True, help='Turn on the matrix balance for contact map. Only effective when input_type == cmap and input_format == cooler')
@click.option('--not-normalize', is_flag=True, default=False, show_default=True, help='Turn off auto normalization of contact map. Only effective when the input is contact map')
@click.option('--enforce-nonnegative-connectivity-matrix', is_flag=True, default=False, show_default=True, help='Enforcing that the "spring constants" in the connectivity matrix can only be nonnegative')
def main(input, output_prefix, connectivity_matrix, ensemble, alpha, selection, method, lamd, reg, iteration, learning_rate, input_type, \
    input_format, binsize, hic_norm, hic_unit, log, no_xyzs, ignore_missing_data, balance, not_normalize, enforce_nonnegative_connectivity_matrix):
    """
    Script to run HIPPS/DIMES to generate ensemble of genome structures from either contact map or mean distance map\n
    INPUT: Specify the path to the input file\n
    OUTPUT_PREFIX: Specify the prefix for output files\n\n
    If you use this program in your publication, please cite this paper: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051\n
    """
    console = Console()

    title = Text.assemble(("HIPPS-DIMES", "bold yellow"),
                          ": Maximum Entropy Based HI-C/Distance Map - Polymer Physics - Structures Reconstruction\n",
                          "1. Shi, Guang, and D. Thirumalai. From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes. Physical Review X 11.1 (2021): 011051.\n",
                          "2. Shi, Guang, and D. Thirumalai. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nature Communications 14.1 (2023): 1150.")
    console.print(Panel(title))

    with console.status("[bold green]System initialization...") as status:
        if input_type == 'dmap':
            console.print("Reading distance matrix from file")
            if input_format == 'text':
                dmap_target = np.loadtxt(input)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
            elif input_format == 'cooler':
                click.echo('input-type=dmap only support text format file')
        elif input_type == 'cmap':
            console.print("Reading contact map from file")
            if input_format == 'text':
                cmap = np.loadtxt(input)
            elif input_format == 'cooler':
                cmap = cooler.Cooler(input)
                console.print("Cooler file read completed")
                cmap = cmap.matrix(balance=balance).fetch(selection)
                console.print("Cooler file selection completed")
                if len(cmap) >= 5000:
                    console.print("The matrix size is {}x{}. It is too large. Please use smaller matrix".format(
                        len(cmap), len(cmap)))
                    exit(0)
            elif input_format == 'hic':
                # parse selection for .hic
                if not selection or ',' not in selection:
                    console.print("[red]For .hic input, --selection must be 'chr1:start1-end1,chr2:start2-end2'")
                    sys.exit(1)
                hic = hicstraw.HiCFile(input)
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
                console.print("Fetched hic matrix zoom data")
                try:
                    cmap = matrix_obj.getRecordsAsMatrix(start1, end1, start2, end2)
                except Exception:
                    # fallback manual assembly via straw
                    console.print("Falling back to manual assembly via hicstraw.straw()...")
                    region1 = f"{chrom1}:{start1}:{end1}"
                    region2 = f"{chrom2}:{start2}:{end2}"
                    result = hicstraw.straw('observed', hic_norm, input,
                                             region1, region2,
                                             hic_unit, binsize)
                    # compute dimensions
                    dim1 = (end1 - start1) // binsize + 1
                    dim2 = (end2 - start2) // binsize + 1
                    cmap = np.zeros((dim1, dim2))
                    # build map
                    for pt in tqdm(result, desc="Building hic contact map"):
                        i = int((pt.binX - start1) / binsize)
                        j = int((pt.binY - start2) / binsize)
                        cmap[i, j] = pt.counts
                    cmap = cmap + cmap.T

                console.print(".hic contact map extracted")
            
            if ignore_missing_data:
                dmap_target = cmap2dmap_missing_data(cmap, alpha, not_normalize)
            else:
                dmap_target = cmap2dmap(cmap, alpha, not_normalize)
            dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        else:
            console.print("[red]Invalid input_type")
            sys.exit(1)

        if connectivity_matrix is not None:
            connectivity_matrix = np.loadtxt(connectivity_matrix)
            console.print("Load the provided connectivity matrix and will use it as initialization.")
        console.print("Initialization completed")


    table = Table(title="Some Basic Parameters")
    table.add_column("Input File", no_wrap=False)
    table.add_column("Input Type", no_wrap=False)
    table.add_column("Input Format", no_wrap=False)
    table.add_column("Optimization method", no_wrap=False)
    table.add_column("Matrix Size", no_wrap=False)
    table.add_column("Number of Iterations", no_wrap=False)
    table.add_column("Regularization", no_wrap=False)
    table.add_column("Ignore Missing Data", no_wrap=False)
    table.add_column("Matrix Balancing", no_wrap=False)
    table.add_column("Matrix Normalization", no_wrap=False)
    table.add_row(input,
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
                  "{}".format("No" if (not_normalize is True and input_type == 'cmap') else "Yes" if (
                      not_normalize is False and input_type == 'cmap') else "N/A")
                  )
    console.print(table)

    model = Optimize(dmap_target, connectivity_matrix=connectivity_matrix)
    keyword_arguments = {'learning_rate': learning_rate, 'lamd': lamd, 'reg': reg, 'method': method,
                         'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix}

    if method == 'IS' or method == 'GD':
        general_method = 'optimization'
    elif method == 'DI':
        general_method = 'direct'

    loss, dmap_maxent, connectivity_matrix = model.run(
        iteration, general_method=general_method, **keyword_arguments)
    try:
        loss = pd.DataFrame(
            np.dstack((np.arange(1, len(loss)+1), loss))[0], columns=['iteration', 'loss'])
    except IndexError:
        pass

    if reg == 'L2':
        print('L2 norm of the connectivity matrix:', np.linalg.norm(
            connectivity_matrix[np.triu_indices_from(connectivity_matrix, k=1)]))
    elif reg == 'L1':
        print('L1 norm of the connectivity matrix:', np.abs(
            connectivity_matrix[np.triu_indices_from(connectivity_matrix, k=1)]).sum())

    console.print("Final loss: {}".format(loss['loss'].values[-1]))

    with console.status("[bold green]System finalizing...") as status:
        if input_type == 'cmap':
            cmap_rc_minimize_res = scipy.optimize.minimize_scalar(
                objective_func, args=(connectivity_matrix, cmap))
            console.print('Optimized contact threshold distance: {}\n'.format(
                cmap_rc_minimize_res.x))
            cmap_maxent = a2cmap_theory(
                connectivity_matrix, cmap_rc_minimize_res.x)

        if log:
            loss.to_csv('{}_loss_function_iteration.csv'.format(output_prefix))
            console.print(
                "Loss function saved to file: [bold magenta]{}_loss_function_iteration.csv[/bold magenta]".format(output_prefix))

        np.savetxt('{}_dmap_final.txt'.format(output_prefix), dmap_maxent)
        console.print(
            "Final distance map saved to file: [bold magenta]{}_dmap_final.txt[/bold magenta]".format(output_prefix))
        if input_type == 'cmap':
            np.savetxt('{}_dmap_target.txt'.format(output_prefix), np.sqrt((8./(3.*np.pi))* dmap_target))
            console.print(
                "Target distance map saved to file: [bold magenta]{}_dmap_target.txt[/bold magenta]".format(output_prefix))
            np.savetxt('{}_cmap_final.txt'.format(output_prefix), cmap_maxent)
            console.print(
                "Final contact map saved to file: [bold magenta]{}_cmap_final.txt[/bold magenta]".format(output_prefix))
        np.savetxt('{}_connectivity_matrix.txt'.format(
            output_prefix), connectivity_matrix)
        console.print(
            'Connectivity matrix saved to file: [bold magenta]{}_connectivity_matrix.txt[/bold magenta]'.format(output_prefix))

        if not no_xyzs:
            xyzs = a2xyz_sample(connectivity_matrix, ensemble=ensemble)
            write2xyz('{}.xyz'.format(output_prefix), xyzs)
            console.print(
                "Ensemble of structures saved to file: [bold magenta]{}.xyz[/bold magenta]".format(output_prefix))


if __name__ == '__main__':
    main()
