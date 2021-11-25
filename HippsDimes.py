"""
Reconstruction of 3D genome organization using the Maximum Entropy Principle

Reference:
Shi, Guang, and Dave Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051.
https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051
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
    a_direct = nearestNSD(- scipy.linalg.pinv(Omega), 0.0)

    return a_direct

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

def dmap2cmap(dmap, rc):
    """
    Return contact map given the mean distance map and the contact threshold
    """
    sigma_mtx = 0.5 * np.sqrt(np.pi / 2.0) * dmap
    cmap = scipy.special.erf(rc/(np.sqrt(2) * sigma_mtx)) - \
            np.sqrt(2.0/np.pi) * np.exp(-0.5 * rc**2.0/np.power(sigma_mtx, 2.0)) * rc / sigma_mtx
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
        position = eigvector @ (np.sqrt(-temp) * np.random.randn(len(eigvalue), 3))
        positions.append(position)

    return np.array(positions)

def interpolate_missing(matrix):
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

def objective_func(rc, A_mtx, cmap_exp):
    x = a2cmap_theory(A_mtx, rc)
    y = cmap_exp / np.nanmax(cmap_exp)
    logx = interpolate_missing(np.log(x))
    logy = interpolate_missing(np.log(y))
    res = np.power(logx[np.triu_indices_from(logx, k=1)] - logy[np.triu_indices_from(logy, k=1)], 2.).mean()**0.5
    return res

# FUNCTION TO CONVERT CMAP TO DMAP
def cmap2dmap_core(cmap_exp, rc, alpha, not_normalize, norm_max=1.0, mode='log'):
    # rc is the prefactor
    # norm_max is the maximum contact probability
    if mode == 'raw':
        if not_normalize:
            log10_pmap = np.log10(cmap_exp)
        else:
            log10_pmap = np.log10(cmap_exp) + np.log10(norm_max) - np.log10(np.nanmax(cmap_exp))
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
#------------------------------------------------------------------#

class Optimize:
    def __init__(self, ddmap_target):
        # ddmap_target is the targeted matrix we would like to match
        # note that ddmap_taret is the mean SQUARED distance matrix, not mean distance matrix
        self.ddmap_target = ddmap_target

        # get the size of system
        self.n = ddmap_target.shape[0]

        # initialize the connectivity matrix
        # here the connectivity matrix is initialized as a simple rouse chain whose spring constant is determined such\
        # that its radius of gyration is close to the target
        rg2 = .5 * np.nanmean(self.ddmap_target)
        k = self.n / (4. * rg2)
        self.A = construct_connectivity_matrix_rouse(self.n, k)

        # initialize the loss
        self.loss = None

    def __compute_loss(self):
        ddmap_t = ((3. * np.pi) / 8.)  * np.power(a2dmap_theory(self.A, force_positive_definite=True), 2.)
        loss = np.nanmean(np.power((ddmap_t - self.ddmap_target)/self.ddmap_target, 2.)) ** .5
        return loss

    def __update_parameter(self, t, learning_rate, lamd=0.0, reg='l2', method='IS'):
        # updating using Iterative Scaling

        # compute the mean squared distance matrix at current iteration step
        ddmap_t = ((3. * np.pi) / 8.)  * np.power(a2dmap_theory(self.A, force_positive_definite=True), 2.)
        # compute the ratio between the current value and the target
        compare_ratio = ddmap_t / self.ddmap_target
        # compute the prefactor for iterative scaling
        fhash = np.nansum(ddmap_t) / 2.

        if method == 'IS':
            # compute the gradient
            if lamd > 0.0:
                if reg == 'L2':
                    gradient_t = (np.nan_to_num(np.log(compare_ratio), posinf=0., neginf=0.) - 2. * lamd * self.A) / fhash
                elif reg == 'L1':
                    gradient_t = (np.nan_to_num(np.log(compare_ratio), posinf=0., neginf=0.) + lamd * np.sign(- self.A)) / fhash
            elif lamd == 0.0:
                gradient_t = np.nan_to_num(np.log(compare_ratio), posinf=0., neginf=0.) / fhash

            # update the connectivity matrix
            self.A += learning_rate * gradient_t
        elif method == 'GD':
            if t == 0:
                self.theta = np.copy(self.A)

            # compute the gradient
            if lamd > 0.0:
                if reg == 'L2':
                    gradient_t = (ddmap_t - self.ddmap_target - 2. * lamd * self.A)
                elif reg == 'L1':
                    gradient_t = (ddmap_t - self.ddmap_target + lamd * np.sign(- self.A))
            elif lamd == 0.0:
                gradient_t = (ddmap_t - self.ddmap_target)

            # perform Nesterov update rule
            # gradient descent state
            theta_previous = np.copy(self.theta)

            #self.theta = self.A + np.maximum(np.minimum(learning_rate * gradient_t, step_cap),-step_cap)
            self.theta = self.A + learning_rate * gradient_t

            #if momentum_rate == None:
            #    momentum_rate = t/(t+3)

            # update the connectivity matrix
            self.A = self.theta + (t/(t+3)) * (self.theta - theta_previous)


        self.A = a2a(self.A)
        # project to be negative semidefinite
        #self.A = nearestNSD(self.A, 0.0)

        # compute the loss
        self.loss = self.__compute_loss()

    def run(self, epoch, **kwargs):
        """
        Main function to run the optimization
        """

        console = Console()

        loss_array = []

        #for t in trange(epoch):
        with trange(epoch, desc="Performing optimization", unit="iteration") as pbar:
            for t in pbar:
                self.__update_parameter(t, **kwargs)
                # display loss at each iterations
                pbar.set_postfix(loss=self.loss)
                loss_array.append(self.loss)

        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)

        return loss_array, dmap_maxent, self.A

@click.command()
@click.argument('input', nargs=1)
@click.argument('output-prefix', nargs=1)
@click.option('-e', '--ensemble', type=int, default=1000, show_default=True, help='specify the number of conformations generated')
@click.option('-a', '--alpha', type=float, default=4.0, show_default=True, help='specify the value of cmap-to-dmap conversion exponent')
@click.option('-s', '--selection', type=str, help='specify which chromosome or region to run the model on if the input file is Hi-C data in cooler format. Accept any valid options for [fetch] method in cooler.Cooler.matrix() selector')
@click.option('-m', '--method', type=click.Choice(['IS','GD'],case_sensitive=True), default='IS', show_default=True,help='specify the method for optimization. IS: Iterative Scaling. GD: Gradient Descent')
@click.option('-l', '--lamd', type=click.FloatRange(0,max=None), default=0.0,show_default=True, help='Specify the weight for the regularization.')
@click.option('-r', '--reg', type=click.Choice(['L1','L2'], case_sensitive=True), default='L2', show_default=True, required=False, help='specify the type of regularization. Currently support L1 and L2 regularization. Note that this option should be used together with option -l')
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
def main(input, output_prefix, ensemble, alpha, selection, method, lamd, reg, iteration, learning_rate, input_type, input_format, log, no_xyzs, ignore_missing_data, balance, not_normalize):
    """
    Script to run HIPPS/DIMES to generate ensemble of genome structures from either contact map or mean distance map\n
    INPUT: Specify the path to the input file\n
    OUTPUT_PREFIX: Specify the prefix for output files\n\n
    If you use this program in your publication, please cite this paper: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051\n
    """
    console = Console()

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
                    dmap_target = cmap2dmap_missing_data(cmap, alpha, not_normalize)
                else:
                    dmap_target = cmap2dmap(cmap, alpha, not_normalize)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
            elif input_format == 'cooler':
                cmap = cooler.Cooler(input)
                console.print("Cooler file read completed")
                cmap = cmap.matrix(balance=balance).fetch(selection)
                console.print("Cooler file selection completed")
                if len(cmap) >= 5000:
                    console.print("The matrix size is {}x{}. It is too large. Please use smaller matrix".format(len(cmap), len(cmap)))
                    exit(0)
                if ignore_missing_data:
                    dmap_target = cmap2dmap_missing_data(cmap, alpha, not_normalize)
                else:
                    dmap_target = cmap2dmap(cmap, alpha, not_normalize)
                dmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
        console.print("Initialization completed")

    title = Text.assemble(("HIPPS-DIMES", "bold yellow"), \
        ": Maximum Entropy Based HI-C/Distance Map - Polymer Physics - Structures Reconstruction\n",\
        "Shi, Guang, and Dave Thirumalai. From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes. Physical Review X 11.1 (2021): 011051.")
    console.print(Panel(title))

    table = Table(title="Information")
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
    table.add_row(input, \
        "{}".format("Contact Map" if input_type=='cmap' else "Distance Map" if input_type=='dmap' else "Unknown"), \
        "{}".format("Text" if input_format=='text' else "Cooler File" if input_format=='cooler' else "Unknown"),\
        "{}".format("Iterative Scaling" if method == 'IS' else "Gradient Descent" if method == 'GD' else "Unknown"),\
        "{}".format("{}x{}".format(dmap_target.shape[0], dmap_target.shape[1])),\
        "{}".format(iteration),\
        "{}".format(reg if lamd > 0.0 else "No" if lamd == 0.0 else "Unknown"),\
        "{}".format("Yes" if ignore_missing_data else "No"),\
        "{}".format("Yes" if balance else "No" if (balance is False and input_format == 'cooler') else "N/A"),\
        "{}".format("No" if (not_normalize is True and input_type == 'cmap') else "Yes" if (not_normalize is False and input_type == 'cmap') else "N/A")
        )
    console.print(table)

    model = Optimize(dmap_target)
    keyword_arguments = {'learning_rate': learning_rate, 'lamd': lamd, 'reg': reg, 'method': method}

    loss, dmap_maxent, connectivity_matrix = model.run(iteration, **keyword_arguments)
    loss = pd.DataFrame(np.dstack((np.arange(1, len(loss)+1), loss))[0], columns=['iteration', 'loss'])

    if reg == 'L2':
        print('L2 norm of the connectivity matrix:', np.linalg.norm(connectivity_matrix[np.triu_indices_from(connectivity_matrix, k=1)]))
    elif reg == 'L1':
        print('L1 norm of the connectivity matrix:', np.abs(connectivity_matrix[np.triu_indices_from(connectivity_matrix, k=1)]).sum())

    console.print("Final loss: {}".format(loss['loss'].values[-1]))

    if input_type == 'cmap':
        cmap_rc_minimize_res = scipy.optimize.minimize_scalar(objective_func, args=(connectivity_matrix, cmap))
        console.print('Optimized contact threshold distance: {}\n'.format(cmap_rc_minimize_res.x))
        cmap_maxent = a2cmap_theory(connectivity_matrix, cmap_rc_minimize_res.x)

    if log:
        loss.to_csv('{}_loss_function_iteration.csv'.format(output_prefix))
        console.print("Loss function saved to file: [bold magenta]{}_loss_function_iteration.csv[/bold magenta]".format(output_prefix))

    np.savetxt('{}_dmap_final.txt'.format(output_prefix), dmap_maxent)
    console.print("Final distance map saved to file: [bold magenta]{}_dmap_final.txt[/bold magenta]".format(output_prefix))
    if input_type == 'cmap':
        np.savetxt('{}_cmap_final.txt'.format(output_prefix), cmap_maxent)
        console.print("Final contact map saved to file: [bold magenta]{}_cmap_final.txt[/bold magenta]".format(output_prefix))
    np.savetxt('{}_connectivity_matrix.txt'.format(output_prefix), connectivity_matrix)
    console.print('Connectivity matrix saved to file: [bold magenta]{}_connectivity_matrix.txt[/bold magenta]'.format(output_prefix))

    if not no_xyzs:
        xyzs = a2xyz_sample(connectivity_matrix, ensemble = ensemble)
        write2xyz('{}.xyz'.format(output_prefix), xyzs)
        console.print("Ensemble of structures saved to file: [bold magenta]{}.xyz[/bold magenta]".format(output_prefix))

if __name__ == '__main__':
    main()
