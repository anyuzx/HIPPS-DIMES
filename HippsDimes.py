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
import cooler
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
    output : np.ndarray
        2D array. First column is time `t`, second column is the ACF values at each time.
    msd : np.ndarray
        1D array of the same length as `t`, giving the mean-square displacement 
        inferred from the ACF. This is returned as a second column alongside `t`.
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

def compute_m1_general_theory(i, t, a, zeta=1.0):
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
        self.xyz = a2xyz(self.A, force_positive_definite = force_positive_definite)
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
@click.option('-s', '--selection', type=str, help='specify which chromosome or region to run the model on if the input file is Hi-C data in cooler format. Accept any valid options for [fetch] method in cooler.Cooler.matrix() selector')
@click.option('-m', '--method', type=click.Choice(['IS', 'GD', 'DI'], case_sensitive=True), default='IS', show_default=True, help='specify the method. IS: Iterative Scaling. GD: Gradient Descent. DI: Direct Inversion. When using\
    Direct Inversion, no iterations are performed. The connectivity matrix is obtained by direct Moore–Penrose inverse of the covariance matrix. Note that the resulting connectivity matrix using Direct Inversion can be very different from the results obtained by GD or IS method.')
@click.option('-l', '--lamd', type=click.FloatRange(0, max=None), default=0.0, show_default=True, help='Specify the weight for the regularization.')
@click.option('-r', '--reg', type=click.Choice(['L1', 'L2'], case_sensitive=True), default='L2', show_default=True, required=False, help='specify the type of regularization. Currently support L1 and L2 regularization. Note that this option should be used together with option -l')
@click.option('-i', '--iteration', type=int, default=10000, show_default=True, help='Number of iterations')
@click.option('-r', '--learning-rate', type=float, default=10.0, show_default=True, help='Learning rate. This hyperparameter controls the speed of convergence. \
    If its value is too small, then convergence is very slow. If its value is too large, the program may never converge. Typically, learning rate can be set to be 1-30 if use Iterative scaling method. \
        It should be a very small value (such as 1e-8) when using gradient descent optimization')
@click.option('--input-type', required=True, type=click.Choice(['cmap', 'dmap'], case_sensitive=False), help='Specify the type of the input. cmap: contact map or dmap: distance map')
@click.option('--input-format', required=True, type=click.Choice(['text', 'cooler'], case_sensitive=False), help='Specify the format of the input. Support pure text format or cooler Hi-C contact map')
@click.option('--log', is_flag=True, default=False, show_default=True, help='Write a log file')
@click.option('--no-xyzs', is_flag=True, default=False, show_default=True, help='Turn off writing conformations to .xyz file')
@click.option('--ignore-missing-data', is_flag=True, default=False, show_default=True, help='Turn on this argument will let the program ignore the missing elementsin the contact map or distance map')
@click.option('--balance', is_flag=True, default=False, show_default=True, help='Turn on the matrix balance for contact map. Only effective when input_type == cmap and input_format == cooler')
@click.option('--not-normalize', is_flag=True, default=False, show_default=True, help='Turn off auto normalization of contact map. Only effective when the input is contact map')
@click.option('--enforce-nonnegative-connectivity-matrix', is_flag=True, default=False, show_default=True, help='Enforcing that the "spring constants" in the connectivity matrix can only be nonnegative')
def main(input, output_prefix, connectivity_matrix, ensemble, alpha, selection, method, lamd, reg, iteration, learning_rate, input_type, \
    input_format, log, no_xyzs, ignore_missing_data, balance, not_normalize, enforce_nonnegative_connectivity_matrix):
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
                if ignore_missing_data:
                    dmap_target = cmap2dmap_missing_data(
                        cmap, alpha, not_normalize)
                else:
                    dmap_target = cmap2dmap(cmap, alpha, not_normalize)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
            elif input_format == 'cooler':
                cmap = cooler.Cooler(input)
                console.print("Cooler file read completed")
                cmap = cmap.matrix(balance=balance).fetch(selection)
                console.print("Cooler file selection completed")
                if len(cmap) >= 5000:
                    console.print("The matrix size is {}x{}. It is too large. Please use smaller matrix".format(
                        len(cmap), len(cmap)))
                    exit(0)
                if ignore_missing_data:
                    dmap_target = cmap2dmap_missing_data(
                        cmap, alpha, not_normalize)
                else:
                    dmap_target = cmap2dmap(cmap, alpha, not_normalize)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        
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
                              'text' else "Cooler File" if input_format == 'cooler' else "Unknown"),
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
