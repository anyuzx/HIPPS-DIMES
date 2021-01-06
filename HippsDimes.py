"""
Reconstruction of 3D genome organization using the Maximum Entropy Principle

Reference: https://www.biorxiv.org/content/10.1101/2020.05.21.109421v1

Usage:
    input 
"""

import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import numpy as np
import scipy
import scipy.linalg
import pandas as pd
import click
import cooler
from tqdm import tqdm

#------------------------------------------------------------------#
# Helper functions
def construct_connectivity_matrix_rouse(n,k):
    """
    Function to construct a ideal chain connectivity matrix given the number of monomers and the spring constant
    """
    A = np.diag(np.full(n-1,k),1)
    A += A.T
    A[np.diag_indices(n)] = -2*k
    A[0,0]=-k
    A[n-1,n-1]=-k
    return A

def a2dmap_theory(A, force_positive_definite = False):
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

def a2a(a, fill_negative=False):
    """
    Correct the connectivity matrx. Make it Laplacian, and non negative (options)
    """
    temp = a
    if fill_negative:
        temp[temp < 0.0] = 0.0
    np.fill_diagonal(temp, - np.sum(a, axis=1) + a.diagonal())
    return temp

def optimal_rotate(P, Q, return_rotation = False, allow_reflection=False):
    """
    Return aligned matrix referred to Q
    Can return rotation matrix if return_rotation is set True
    """
    # P and Q are two sets of vectors
    P = np.matrix(P)
    Q = np.matrix(Q)

    assert P.shape == Q.shape

    Qc = np.mean(Q,axis=0)

    P = P - np.mean(P,axis=0)
    Q = Q - np.mean(Q,axis=0)

    # calculate covariance matrix A = (P^T)Q
    A = P.T * Q

    # SVD for matrix A
    V, S, Wt = np.linalg.svd(A)

    # correct rotation matrix to ensure a right-handed system if necessary
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0

    if not allow_reflection:
        if d:
            S[-1] = -S[-1]
            V[:,-1] = -V[:,-1]

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
                f.write('{} {} {} {}\n'.format(int(idx+1), item[0], item[1], item[2]))

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
        position = eigvector @ (np.sqrt(-temp) * np.random.randn(len(eigvalue), 3))
        positions.append(position)

    return np.array(positions)

def interpolate_miss(matrix):
    matrix_copy = np.copy(matrix)
    x = np.arange(0, matrix_copy.shape[1])
    y = np.arange(0, matrix_copy.shape[0])
    #mask invalid values
    matrix_copy = np.ma.masked_invalid(matrix_copy)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~matrix_copy.mask]
    y1 = yy[~matrix_copy.mask]
    newarr = matrix_copy[~matrix_copy.mask]

    GD1 = scipy.interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy), method='nearest')
    return GD1


# FUNCTION TO CONVERT CMAP TO DMAP
def cmap2dmap_core(cmap_exp, rc, alpha, norm_max=1.0, mode='log'):
    # rc is the prefactor
    # norm_max is the maximum contact probability
    if mode == 'raw':
        log10_pmap = np.log10(cmap_exp) + np.log10(norm_max) - np.log10(np.max(cmap_exp))
    elif mode == 'log':
        log10_pmap = cmap_exp + np.log10(norm_max) - np.max(cmap_exp)

    return rc * 10 ** (-1.0/alpha * log10_pmap)

def cmap2dmap(cmap, alpha):
    # cmap is the raw data
    # we take log on contact map
    # and then interpolate the missing data. Any zero contact pair will be interpolated
    cmap_log = interpolate_missing(np.log10(cmap))
    cmap_log = np.array((cmap_log + cmap_log.T) / 2.)
    # lastly, convert to distance map using value of alpha
    dmap = cmap2dmap_core(cmap_log, 1.0, alpha)
    return dmap
#------------------------------------------------------------------#

class optimize:
    def __init__(self, dmap_target, mode='second moment'):
        # dmap_target is the targeted matrix we would like to match
        # the provided dmap_target should be consistent with the mode argument
        # For instance, if the mode is the second moment, then the provided dmap_target should be averaged squared distance matrix
        # if the mode is the first moment, then the provided dmap_target should be the averaged distance matrix

        self.dmap_target = dmap_target
        self.n = dmap_target.shape[0]
        self.mode = mode

        self.A = np.zeros((self.n, self.n))
        self.A_tot = construct_connectivity_matrix_rouse(self.n, self.n / self.dmap_target.max())

    def __update_parameter(self, learning_rate):
        rij = ((3. * np.pi) / 8.)  * np.power(a2dmap_theory(a2a(self.A_tot)), 2.)
        compare_ratio = rij / self.dmap_target
        fhash = np.sum(rij) / 2.
        self.A_tot += learning_rate * np.log(compare_ratio) / fhash
        np.fill_diagonal(self.A_tot, 0.0)
        self.A_tot = a2a(self.A_tot)
        self.cost = np.sqrt(np.mean(np.power(rij - self.dmap_target, 2.)) / np.mean(np.power(self.dmap_target, 2.)))

        return a2dmap_theory(self.A_tot)

    def run(self, epoch, learning_rate):
        """
        Main function to run the optimization
        """

        cost_array = []
        for _ in tqdm(range(epoch)):
            dmap_maxent = self.__update_parameter(learning_rate)
            cost_array.append(self.cost)

        return cost_array, dmap_maxent, self.A_tot

@click.command()
@click.argument('input', nargs=1)
@click.argument('output-prefix', nargs=1)
@click.option('-e', '--ensemble', type=int, default=1000, show_default=True, help='specify the number of conformations generated')
@click.option('-a', '--alpha', type=float, default=4.0, show_default=True, help='specify the value of cmap-to-dmap conversion exponent')
@click.option('-s', '--selection', type=str, help='specify which chromosome or region to run the model on if the input file is Hi-C data in cooler format. Accept any valid options for [fetch] method in cooler.Cooler.matrix() selector')
@click.option('-i', '--iteration', type=int, default=10000, show_default=True)
@click.option('-r', '--learning-rate', type=float, default=10.0, show_default=True)
@click.option('--input-type', required=True, type=click.Choice(['cmap', 'dmap'], case_sensitive=False))
@click.option('--input-format', required=True, type=click.Choice(['text', 'cooler'], case_sensitive=False))
@click.option('--log', is_flag=True, default=False, show_default=True, help='write a log file')
@click.option('--no-xyzs', is_flag=True, default=False, show_default=True, help='turn off writing conformations to .xyz file')
def main(input, output_prefix, ensemble, alpha, selection, iteration, learning_rate, input_type, input_format, log, no_xyzs):
    """
    Script to run HIPPS/DIMES to generate ensemble of genome structures from either contact map or mean distance map\n
    INPUT: Specify the path to the input file\n
    OUTPUT_PREFIX: Specify the prefix for output files\n\n
    Reference: https://www.biorxiv.org/content/10.1101/2020.05.21.109421v1\n
    """
    if input_type == 'dmap':
        if input_format == 'text':
            dmap_target = np.loadtxt(input)
            dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        elif input_format == 'cooler':
            click.echo('input-type=dmap only support text format file')
    elif input_type == 'cmap':
        if input_format == 'text':
            cmap = np.loadtxt(input)
            dmap_target = cmap2dmap(cmap, alpha)
            dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        elif input_format == 'cooler':
            cmap = cooler.Cooler(input)
            cmap = cmap.matrix(balance=False).fetch(selection)
            dmap_target = cmap2dmap(cmap, alpha)
            dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)

    model = optimize(dmap_target)
    cost, dmap_maxent, connectivity_matrix = model.run(iteration, learning_rate)
    cost = pd.DataFrame(np.dstack((np.arange(1, iteration+1), cost))[0], columns=['iteration', 'cost'])

    if log:
        cost.to_csv('cost_function_iteration.csv')
    pass
    
    np.savetxt('{}_dmap_final.txt'.format(output_prefix), dmap_maxent)
    np.savetxt('{}_connectivity_matrix.txt'.format(output_prefix), connectivity_matrix)

    if not no_xyzs:
        xyzs = a2xyz_sample(connectivity_matrix, ensemble = ensemble)
        write2xyz('{}.xyz'.format(output_prefix), xyzs)
    

if __name__ == '__main__':
    main()
