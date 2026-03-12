"""Command-line interface for HIPPS-DIMES."""

try:
    import rich_click as click
except ImportError:
    import click

from .core import run_optimization, set_eigh_num_threads


def _parse_save_steps(save_steps_str):
    """Parse comma-separated save_steps string into list of integers."""
    if not save_steps_str:
        return None

    try:
        steps = [int(s.strip()) for s in save_steps_str.split(',')]
        return sorted(list(set(steps)))
    except ValueError as e:
        raise ValueError(
            f"Invalid save_steps format: '{save_steps_str}'. "
            "Expected comma-separated integers, e.g., '1000,5000,10000'"
        ) from e


@click.command()
@click.argument('input', nargs=1)
@click.argument('output-prefix', nargs=1)
@click.option('--input-type', required=True, type=click.Choice(['cmap', 'dmap', 'ddmap'], case_sensitive=False), help='Specify the type of the input. cmap: contact map, dmap: mean distance map, or ddmap: mean squared distance map')
@click.option('--input-format', required=True, type=click.Choice(['text', 'npy', 'cooler', 'hic'], case_sensitive=False), help='Format of input: text, npy, cooler, or hic')
@click.option('-i', '--iteration', type=int, default=10000, show_default=True, help='Number of iterations')
@click.option('-k', '--connectivity-matrix', type=str, required=False, help='Use provided connectivity matrix as initialization. Useful when restart from previous run')
@click.option('-e', '--ensemble', type=int, default=1000, show_default=True, help='specify the number of conformations generated')
@click.option('-a', '--alpha', type=float, default=4.0, show_default=True, help='specify the value of cmap-to-dmap conversion exponent')
@click.option('-s', '--selection', type=str, required=False, help='For cooler: any valid selector for cooler.Cooler.matrix().fetch(), e.g. "chr1" or "chr1::start-end". For .hic: use "chr1:start1-end1,chr2:start2-end2"')
@click.option('-m', '--method', type=click.Choice(['IS', 'GD', 'DI'], case_sensitive=True), default='IS', show_default=True, help='specify the method. IS: Iterative Scaling. GD: Gradient Descent. DI: Direct Inversion. When using Direct Inversion, no iterations are performed. The connectivity matrix is obtained by direct Moore–Penrose inverse of the covariance matrix. Note that the resulting connectivity matrix using Direct Inversion can be very different from the results obtained by GD or IS method.')
@click.option('-l', '--lamd', type=click.FloatRange(0, max=None), default=0.0, show_default=True, help='Specify the weight for the regularization.')
@click.option('-r', '--reg', type=click.Choice(['L1', 'L2'], case_sensitive=True), default='L2', show_default=True, required=False, help='specify the type of regularization. Currently support L1 and L2 regularization. Note that this option should be used together with option -l')
@click.option('--gaussian-noise-variance', type=click.FloatRange(0, max=None), default=0.0, show_default=True, help='Noise variance for independent Gaussian noise on constraints (IS/GD only). Cannot be combined with --lamd.')
@click.option('--learning-rate', type=float, default=10.0, show_default=True, help='Learning rate. This hyperparameter controls the speed of convergence. If its value is too small, then convergence is very slow. If its value is too large, the program may never converge. Typically, learning rate can be set to be 1-30 if use Iterative scaling method. It should be a very small value (such as 1e-8) when using gradient descent optimization')
@click.option('--momentum', type=click.FloatRange(0, 1), default=0.0, show_default=True, help='Momentum coefficient for IS method. RECOMMENDED: Use 0.95 with --nesterov for fastest convergence (~50%% faster). Use 0.9 for conservative settings. Only applies when method=IS.')
@click.option('--nesterov', is_flag=True, default=False, show_default=True, help='Use Nesterov Accelerated Gradient (NAG). Enables higher momentum (0.95) without divergence. RECOMMENDED: Use with --momentum 0.95 for fastest convergence.')
@click.option('--use-gpu', is_flag=True, default=False, show_default=True, help='Use GPU acceleration via CuPy. Provides 2-4x speedup for large matrices (n >= 200). Requires CuPy: conda install -c conda-forge cupy')
@click.option('--gpu-float32', is_flag=True, default=False, show_default=True, help='When using --use-gpu, run GPU math/eigendecomposition in float32 (often faster, slightly different numerics).')
@click.option('--binsize', type=int, default=25000, show_default=True, help='Bin size (resolution) for .hic format in bp')
@click.option('--norm', 'hic_norm', type=str, default='KR', show_default=True, help='Normalization for .hic: KR, VC, NONE')
@click.option('--unit', 'hic_unit', type=click.Choice(['BP', 'FRAG'], case_sensitive=False), default='BP', show_default=True, help='Unit for .hic: BP or FRAG')
@click.option('--no-log', is_flag=True, default=False, show_default=True, help='Disable writing run-parameter and iteration-series log files')
@click.option('--no-xyzs', is_flag=True, default=False, show_default=True, help='Turn off writing conformations to .xyz file')
@click.option('--ignore-missing-data', is_flag=True, default=False, show_default=True, help='Turn on this argument will let the program ignore the missing elementsin the contact map or distance map')
@click.option('--balance', is_flag=True, default=False, show_default=True, help='Turn on the matrix balance for contact map. Only effective when input_type == cmap and input_format == cooler')
@click.option('--neighbor-balance', is_flag=True, default=False, show_default=True, help='Turn on neighbor balancing for contact map. Only effective when input_type == cmap. Normalizes contact between i and j by dividing it by the geometric mean of neighbor contact for i and j. see Paggi, Zhang 2025 for method details')
@click.option('--not-normalize', is_flag=True, default=False, show_default=True, help='Turn off auto normalization of contact map. Only effective when the input is contact map')
@click.option('--enforce-nonnegative-connectivity-matrix', is_flag=True, default=False, show_default=True, help='Enforcing that the "spring constants" in the connectivity matrix can only be nonnegative')
@click.option('--save-steps', type=str, default=None, help='Comma-separated list of iteration steps at which to save the connectivity matrix. Example: --save-steps 1000,5000,10000,50000')
@click.option('--save-pickle', is_flag=True, default=False, show_default=True, help='Save the returned results dictionary to {output_prefix}_HIPPS_DIMES_results.pkl and suppress the default text/CSV/XYZ file outputs')
@click.option('--eigh-threads', type=int, default=None, help='Number of threads for eigenvalue (eigh) and BLAS/LAPACK. If not set, backend default is used. Set to 1 for single-threaded.')
@click.option('--quiet', '-q', is_flag=True, default=False, show_default=True, help='Quiet mode: disable fancy tables display, keep only the progress bar.')
def main(input, output_prefix, connectivity_matrix, ensemble, alpha, selection, method, lamd, reg, gaussian_noise_variance, iteration, learning_rate, momentum, nesterov, use_gpu, input_type, gpu_float32, input_format, binsize, hic_norm, hic_unit, no_log, no_xyzs, ignore_missing_data, balance, not_normalize, neighbor_balance, enforce_nonnegative_connectivity_matrix, save_steps, save_pickle, eigh_threads, quiet):
    """CLI for HIPPS-DIMES.

    INPUT: Path to the input file.
    OUTPUT_PREFIX: Prefix for output files.
    """
    if eigh_threads is not None:
        set_eigh_num_threads(eigh_threads)

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
        gaussian_noise_variance=gaussian_noise_variance,
        iteration=iteration,
        learning_rate=learning_rate,
        momentum=momentum,
        nesterov=nesterov,
        use_gpu=use_gpu,
        gpu_float32=gpu_float32,
        input_type=input_type,
        input_format=input_format,
        binsize=binsize,
        hic_norm=hic_norm,
        hic_unit=hic_unit,
        no_log=no_log,
        no_xyzs=no_xyzs,
        ignore_missing_data=ignore_missing_data,
        balance=balance,
        not_normalize=not_normalize,
        neighbor_balance=neighbor_balance,
        enforce_nonnegative_connectivity_matrix=enforce_nonnegative_connectivity_matrix,
        save_steps=_parse_save_steps(save_steps) if save_steps else None,
        save_pickle=save_pickle,
        eigh_threads=eigh_threads,
        verbose=not quiet,
    )


if __name__ == '__main__':
    main()
