"""High-level programmatic API for HIPPS-DIMES."""

import json
import time

import numpy as np
import pandas as pd
import scipy
from rich import print

from .models import Optimize
from .numerics import *  # noqa: F401,F403


def _build_iteration_series_frame(loss, entropy, extra_series=None):
    """Build a DataFrame for per-iteration scalar outputs."""
    series_data = {
        'iteration': np.arange(1, len(loss) + 1),
        'loss': loss,
        'entropy': entropy,
    }
    if extra_series:
        expected_length = len(loss)
        for column_name, values in extra_series.items():
            if len(values) != expected_length:
                raise ValueError(
                    f"Iteration series '{column_name}' must have length {expected_length}, "
                    f"got {len(values)}"
                )
            series_data[column_name] = values
    return pd.DataFrame(series_data)


def _serialize_run_parameter(value):
    """Serialize parameter values for storage in a simple key/value CSV."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value)
    return value


def _build_run_parameters_frame(parameters):
    """Build a DataFrame for run parameter logging."""
    return pd.DataFrame({
        'parameter': list(parameters.keys()),
        'value': [_serialize_run_parameter(v) for v in parameters.values()],
    })


def _ensure_square_matrix(matrix, matrix_name):
    """Return a NumPy array and validate that it is a square 2D matrix."""
    array = np.asarray(matrix)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(
            f"{matrix_name} must be a square 2D array, got shape {array.shape}"
        )
    return array


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
                     gaussian_noise_variance=0.0,
                     iteration=10000,
                     learning_rate=10.0,
                     momentum=0.0,
                     nesterov=False,
                     use_gpu=False,
                     gpu_float32=False,
                     input_type='cmap',
                     input_format='text',
                     binsize=25000,
                     hic_norm='KR',
                     hic_unit='BP',
                     no_log=False,
                     no_xyzs=False,
                     ignore_missing_data=False,
                     balance=False,
                     not_normalize=False,
                     neighbor_balance=False,
                     enforce_nonnegative_connectivity_matrix=False,
                     save_steps=None,
                     eigh_threads=None,
                     verbose=True,
                     log=None,
                     progress_callback=None,
                     show_progress=True):
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
    gaussian_noise_variance : float, default=0.0
        Noise variance for independent Gaussian noise on constraints.
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
        Provides 2-4x speedup for matrices with n >= 200.
        Requires: conda install -c conda-forge cupy
    input_type : str, default='cmap'
        Type of input:
        - 'cmap': contact map
        - 'dmap': mean distance map (converted internally to mean squared distance map)
        - 'ddmap': mean squared distance map (used directly)
    input_format : str, default='text'
        Format of input file: 'text', 'npy', 'cooler', or 'hic'
    binsize : int, default=25000
        Bin size for .hic format in bp
    hic_norm : str, default='KR'
        Normalization for .hic: 'KR', 'VC', 'NONE'
    hic_unit : str, default='BP'
        Unit for .hic: 'BP' or 'FRAG'
    no_log : bool, default=False
        If True, skip writing both log files when output_prefix is provided
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
    save_steps : list of int, optional
        Iteration steps at which to capture the connectivity matrix.
        When set, results include ``connectivity_matrix_at_steps`` (dict: step -> matrix)
        for library use. If ``output_prefix`` is also set, files are saved as
        ``{output_prefix}_connectivity_matrix_iter{step}.txt``.
    eigh_threads : int, optional
        Number of threads for eigenvalue (eigh) and BLAS/LAPACK.
        If None, backend default is used. Set to 1 for single-threaded.
        Must be set before any eigh call (e.g. via set_eigh_num_threads) for effect.
    verbose : bool, default=True
        Whether to print status messages
    progress_callback : callable, optional
        Callback receiving structured per-iteration progress dictionaries.
        Payload keys include iteration, total, loss, entropy,
        iterations_per_sec, method, general_method, stage, noisy, and use_gpu.
    show_progress : bool, default=True
        Whether to render solver progress bars. Set to False when consuming
        progress_callback programmatically.
        
    Returns
    -------
    results : dict
        Dictionary containing:
        - 'iteration_series': Optimization iteration-series data as a DataFrame
        - 'log': Alias for 'iteration_series' (backward compatibility)
        - 'run_parameters': Run parameters as a DataFrame with columns ['parameter', 'value']
        - 'dmap_final': Final distance map (numpy array)
        - 'connectivity_matrix': Final connectivity matrix (numpy array)
        - 'connectivity_matrix_at_steps': Dict of step -> connectivity matrix (only if save_steps was set)
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
    - In practice, this provides 2-4x speedup for matrices with n >= 200
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
    start_time = time.perf_counter()

    # Validate inputs
    valid_input_types = {'cmap', 'dmap', 'ddmap'}
    valid_input_formats = {'text', 'npy', 'cooler', 'hic'}
    if input_matrix is None and input_path is None:
        raise ValueError("Either input_matrix or input_path must be provided")
    if input_type not in valid_input_types:
        raise ValueError(
            f"Invalid input_type '{input_type}'. Must be one of {sorted(valid_input_types)}"
        )
    if input_format not in valid_input_formats:
        # Common user mistake: swapping input_type and input_format when using numpy arrays.
        if input_matrix is not None and input_format in valid_input_types:
            raise ValueError(
                f"input_format='{input_format}' is invalid. Did you mean input_type='{input_format}'? "
                "For numpy arrays, keep input_format as 'text' (default)."
            )
        raise ValueError(
            f"Invalid input_format '{input_format}'. Must be one of {sorted(valid_input_formats)}"
        )
    if gaussian_noise_variance < 0.0:
        raise ValueError("gaussian_noise_variance must be non-negative")
    if gaussian_noise_variance > 0.0 and method == 'DI':
        raise ValueError("gaussian_noise_variance is only supported for optimization methods (IS/GD), not DI")
    if gaussian_noise_variance > 0.0 and lamd > 0.0:
        raise ValueError("gaussian_noise_variance (noise variance) cannot be combined with lamd regularization")
    if log is not None:
        no_log = not log
    if eigh_threads is not None:
        set_eigh_num_threads(eigh_threads)

    input_source = input_path if input_path else "NumPy array"
    connectivity_matrix_source = "default initialization"
    if connectivity_matrix is not None:
        connectivity_matrix_source = connectivity_matrix if isinstance(connectivity_matrix, str) else "provided matrix"
    save_target_cmap = input_type == 'cmap' and input_format in {'cooler', 'hic'}
    
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
                console.print("Reading mean distance matrix")
            if input_matrix is not None:
                # Use provided matrix
                ddmap_target = _ensure_square_matrix(input_matrix, "Distance map")
                ddmap_target = ((3. * np.pi) / 8.) * np.power(ddmap_target, 2.)
            elif input_format == 'text':
                ddmap_target = _ensure_square_matrix(
                    np.loadtxt(input_path), "Distance map"
                )
                ddmap_target = ((3. * np.pi) / 8.) * np.power(ddmap_target, 2.)
            elif input_format == 'npy':
                ddmap_target = _ensure_square_matrix(
                    np.load(input_path), "Distance map"
                )
                ddmap_target = ((3. * np.pi) / 8.) * np.power(ddmap_target, 2.)
            elif input_format in {'cooler', 'hic'}:
                raise ValueError("input_type='dmap' only supports input_format='text' or 'npy' (or provide input_matrix)")
            else:
                raise ValueError(
                    f"Invalid input_format '{input_format}' for input_type='dmap'. "
                    "Supported: 'text' or 'npy' (or provide input_matrix)."
                )
        elif input_type == 'ddmap':
            if verbose and console:
                console.print("Reading mean squared distance matrix")
            if input_matrix is not None:
                # Use provided matrix directly as ddmap constraints
                ddmap_target = _ensure_square_matrix(
                    input_matrix, "Squared distance map"
                )
            elif input_format == 'text':
                ddmap_target = _ensure_square_matrix(
                    np.loadtxt(input_path), "Squared distance map"
                )
            elif input_format == 'npy':
                ddmap_target = _ensure_square_matrix(
                    np.load(input_path), "Squared distance map"
                )
            elif input_format in {'cooler', 'hic'}:
                raise ValueError("input_type='ddmap' only supports input_format='text' or 'npy' (or provide input_matrix)")
            else:
                raise ValueError(
                    f"Invalid input_format '{input_format}' for input_type='ddmap'. "
                    "Supported: 'text' or 'npy' (or provide input_matrix)."
                )
        elif input_type == 'cmap':
            if verbose and console:
                console.print("Reading contact map")
            if input_matrix is not None:
                # Use provided matrix
                cmap = _ensure_square_matrix(input_matrix, "Contact map")
            elif input_format == 'text':
                cmap = _ensure_square_matrix(np.loadtxt(input_path), "Contact map")
            elif input_format == 'npy':
                cmap = _ensure_square_matrix(np.load(input_path), "Contact map")
            elif input_format == 'cooler':
                if cooler is None:
                    raise ImportError(
                        "cooler support is unavailable in this environment. "
                        "Install/import `cooler` and required system resources, "
                        "or use input_format='text'."
                    )
                cmap_data = cooler.Cooler(input_path)
                if verbose and console:
                    console.print("Cooler file read completed")
                cmap = cmap_data.matrix(balance=balance).fetch(selection)
                cmap = _ensure_square_matrix(cmap, "Contact map")
                if verbose and console:
                    console.print("Cooler file selection completed")
                if len(cmap) >= 5000:
                    warning_msg = "The matrix size is {}x{}. It is too large. Please use smaller matrix".format(
                        len(cmap), len(cmap))
                    if verbose and console:
                        console.print(warning_msg)
                    raise ValueError(warning_msg)
            elif input_format == 'hic':
                if hicstraw is None:
                    raise ImportError(
                        ".hic support is unavailable in this environment. "
                        "Install/import `hicstraw`, or use input_format='text'/'cooler'."
                    )
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
                same_region = (
                    chrom1 == chrom2 and start1 == start2 and end1 == end2
                )
                
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
                    if same_region:
                        cmap = cmap + cmap.T

                cmap = _ensure_square_matrix(cmap, "Contact map")
                if verbose and console:
                    console.print(".hic contact map extracted")
            else:
                raise ValueError(
                    f"Invalid input_format '{input_format}' for input_type='cmap'. "
                    "Supported: 'text', 'npy', 'cooler', 'hic' (or provide input_matrix)."
                )
            
            # Apply neighbor balancing if requested
            if neighbor_balance:
                if verbose and console:
                    console.print("Applying neighbor balancing to contact map (see Paggi, Zhang 2025)")
                cmap = neighbor_balance_symmetric(cmap, not_normalize=not_normalize)
            
            if ignore_missing_data:
                ddmap_target = cmap2dmap_missing_data(cmap, alpha, not_normalize)
            else:
                ddmap_target = cmap2dmap(cmap, alpha, not_normalize)
            ddmap_target = ((3. * np.pi) / 8.) * np.power(ddmap_target, 2.)
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

    # Print parameter tables if verbose
    if verbose and console:
        # Table 1: Input & Data Settings
        input_table = Table(title="Input & Data Settings", show_header=True, header_style="bold cyan")
        input_table.add_column("Parameter", style="dim", width=20)
        input_table.add_column("Value", style="green")
        
        input_type_str = (
            "Contact Map" if input_type == 'cmap'
            else "Distance Map" if input_type == 'dmap'
            else "Squared Distance Map" if input_type == 'ddmap'
            else "Unknown"
        )
        input_format_str = (
            "Text" if input_format == 'text'
            else "NumPy (.npy)" if input_format == 'npy'
            else "Cooler File" if input_format == 'cooler'
            else ".hic file" if input_format == 'hic'
            else "Unknown"
        )
        
        input_table.add_row("Input Source", input_source)
        input_table.add_row("Input Type", input_type_str)
        input_table.add_row("Input Format", input_format_str)
        input_table.add_row("Matrix Size", f"{ddmap_target.shape[0]} × {ddmap_target.shape[1]}")
        if input_type == 'cmap':
            input_table.add_row("Alpha (cmap→dmap)", f"{alpha}")
            input_table.add_row("Matrix Balancing", "Yes" if balance else "No" if input_format == 'cooler' else "N/A")
            input_table.add_row("Neighbor Balancing", "Yes" if neighbor_balance else "No")
            input_table.add_row("Auto Normalization", "No" if not_normalize else "Yes")
        input_table.add_row("Ignore Missing Data", "Yes" if ignore_missing_data else "No")
        
        console.print(input_table)
        console.print()  # spacing
        
        # Table 2: Optimization Settings
        opt_table = Table(title="Optimization Settings", show_header=True, header_style="bold cyan")
        opt_table.add_column("Parameter", style="dim", width=20)
        opt_table.add_column("Value", style="green")
        
        method_str = "Iterative Scaling (IS)" if method == 'IS' else "Gradient Descent (GD)" if method == 'GD' else "Direct Inversion (DI)" if method == 'DI' else "Unknown"
        opt_table.add_row("Method", method_str)
        
        if method != 'DI':
            opt_table.add_row("Iterations", f"{iteration:,}")
            opt_table.add_row("Learning Rate", f"{learning_rate}")
            
            # Momentum and Nesterov (only for IS)
            if method == 'IS':
                if momentum > 0:
                    momentum_str = f"{momentum}"
                    if nesterov:
                        momentum_str += " [bold green](+ Nesterov)[/bold green]"
                    opt_table.add_row("Momentum", momentum_str)
                else:
                    opt_table.add_row("Momentum", "Disabled")
            
            # Regularization
            if lamd > 0.0:
                opt_table.add_row("Regularization", f"{reg} (λ = {lamd})")
            else:
                opt_table.add_row("Regularization", "None")
            if gaussian_noise_variance > 0.0:
                opt_table.add_row("Noise Variance", f"{gaussian_noise_variance}")
        
        # GPU
        gpu_status = "[green]Enabled[/green]" if use_gpu and is_gpu_available() else "[yellow]Disabled[/yellow]"
        if use_gpu and is_gpu_available():
            gpu_status += f" ({get_gpu_name()})"
        elif use_gpu and not is_gpu_available():
            gpu_status = "[red]Requested but CuPy not available[/red]"
        opt_table.add_row("GPU Acceleration", gpu_status)
        
        # Eigh/BLAS threads
        eigh_threads_str = str(eigh_threads) if eigh_threads is not None else "backend default"
        opt_table.add_row("Eigh/BLAS Threads", eigh_threads_str)
        
        # Constraints
        opt_table.add_row("Nonnegative Springs", "Yes" if enforce_nonnegative_connectivity_matrix else "No")
        
        console.print(opt_table)
        console.print()  # spacing
        
        # Table 3: Output Settings
        output_table = Table(title="Output Settings", show_header=True, header_style="bold cyan")
        output_table.add_column("Parameter", style="dim", width=20)
        output_table.add_column("Value", style="green")
        
        output_table.add_row("Ensemble Size", f"{ensemble:,} structures")
        output_table.add_row("Output Prefix", output_prefix if output_prefix else "[dim]None (results returned only)[/dim]")
        output_table.add_row("Write XYZ File", "No" if no_xyzs else "Yes")
        output_table.add_row("Write Log Files", "No" if no_log else "Yes")
        if save_steps:
            output_table.add_row("Save Steps", ", ".join(str(s) for s in save_steps))
        
        console.print(output_table)

    # GPU availability tip (only if not using GPU and it's available or could help)
    if verbose:
        if is_gpu_available() and not use_gpu:
            gpu_name = get_gpu_name()
            console.print(f"\n[cyan]💡 Tip: GPU detected ({gpu_name}). Use --use-gpu for 2-4x speedup on large matrices.[/cyan]")
        elif not is_gpu_available() and ddmap_target.shape[0] >= 200:
            console.print("\n[cyan]💡 Tip: For large matrices, GPU can provide 2-4x speedup. Install CuPy: conda install -c conda-forge cupy[/cyan]")
    
    # Run optimization
    model = Optimize(ddmap_target, connectivity_matrix=connectivity_matrix, use_gpu=use_gpu, gpu_float32=gpu_float32)
    
    if use_gpu and model.use_gpu and verbose:
        dtype_str = "float32" if getattr(model, "gpu_float32", False) else "float64"
        console.print(f"[green]GPU acceleration enabled ({get_gpu_name()}), dtype={dtype_str}[/green]")
    
    keyword_arguments = {'learning_rate': learning_rate, 'lamd': lamd, 'reg': reg, 'method': method,
                         'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix,
                         'momentum': momentum, 'nesterov': nesterov}

    if method == 'IS' or method == 'GD':
        general_method = 'optimization'
    elif method == 'DI':
        general_method = 'direct'

    if gaussian_noise_variance > 0.0:
        keyword_arguments_noisy = {
            'learning_rate': learning_rate,
            'method': method,
            'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix,
            'momentum': momentum,
            'nesterov': nesterov
        }
        loss, entropy, dmap_maxent, final_connectivity_matrix, connectivity_at_steps = model.run_noisy(
            iteration, gaussian_noise_variance=gaussian_noise_variance, general_method=general_method, save_steps=save_steps,
            output_prefix=output_prefix, progress_callback=progress_callback, show_progress=show_progress,
            **keyword_arguments_noisy)
    else:
        loss, entropy, dmap_maxent, final_connectivity_matrix, connectivity_at_steps = model.run(
            iteration, general_method=general_method, save_steps=save_steps,
            output_prefix=output_prefix, progress_callback=progress_callback, show_progress=show_progress,
            **keyword_arguments)
    
    # Format per-iteration scalar outputs.
    iteration_series_df = _build_iteration_series_frame(loss, entropy)

    # Print regularization norms if requested
    if verbose:
        if reg == 'L2':
            print('L2 norm of the connectivity matrix:', np.linalg.norm(
                final_connectivity_matrix[np.triu_indices_from(final_connectivity_matrix, k=1)]))
        elif reg == 'L1':
            print('L1 norm of the connectivity matrix:', np.abs(
                final_connectivity_matrix[np.triu_indices_from(final_connectivity_matrix, k=1)]).sum())

        if len(iteration_series_df) > 0 and console:
            console.print("Final loss: {}".format(iteration_series_df['loss'].values[-1]))
            console.print("Final entropy: {}".format(iteration_series_df['entropy'].values[-1]))

    run_parameters_df = _build_run_parameters_frame({
        'input_source': input_source,
        'output_prefix': output_prefix,
        'input_type': input_type,
        'input_format': input_format,
        'matrix_rows': ddmap_target.shape[0],
        'matrix_cols': ddmap_target.shape[1],
        'connectivity_matrix_source': connectivity_matrix_source,
        'ensemble': ensemble,
        'alpha': alpha,
        'selection': selection,
        'method': method,
        'lamd': lamd,
        'reg': reg,
        'gaussian_noise_variance': gaussian_noise_variance,
        'iteration': iteration,
        'learning_rate': learning_rate,
        'momentum': momentum,
        'nesterov': nesterov,
        'use_gpu_requested': use_gpu,
        'use_gpu_enabled': model.use_gpu,
        'gpu_float32': gpu_float32,
        'binsize': binsize,
        'hic_norm': hic_norm,
        'hic_unit': hic_unit,
        'no_log': no_log,
        'no_xyzs': no_xyzs,
        'ignore_missing_data': ignore_missing_data,
        'balance': balance,
        'not_normalize': not_normalize,
        'neighbor_balance': neighbor_balance,
        'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix,
        'save_target_cmap': save_target_cmap,
        'save_steps': save_steps or [],
        'eigh_threads': eigh_threads,
        'verbose': verbose,
        'show_progress': show_progress,
    })

    results = {
        'iteration_series': iteration_series_df,
        'log': iteration_series_df,
        'run_parameters': run_parameters_df,
        'dmap_final': dmap_maxent,
        'connectivity_matrix': final_connectivity_matrix
    }
    if connectivity_at_steps:
        results['connectivity_matrix_at_steps'] = connectivity_at_steps
    
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
            if not no_log:
                run_parameters_df.to_csv('{}_run_parameters.csv'.format(output_prefix), index=False)
                if verbose and console:
                    console.print(
                        "Run parameters saved to file: [bold magenta]{}_run_parameters.csv[/bold magenta]".format(output_prefix))
                iteration_series_df.to_csv('{}_iteration_series.csv'.format(output_prefix), index=False)
                if verbose and console:
                    console.print(
                        "Iteration-series data saved to file: [bold magenta]{}_iteration_series.csv[/bold magenta]".format(output_prefix))

            np.savetxt('{}_dmap_final.txt'.format(output_prefix), dmap_maxent)
            if verbose and console:
                console.print(
                    "Final distance map saved to file: [bold magenta]{}_dmap_final.txt[/bold magenta]".format(output_prefix))
            
            if save_target_cmap:
                np.savetxt('{}_cmap_target.txt'.format(output_prefix), cmap)
                if verbose and console:
                    console.print(
                        "Internal target contact map saved to file: [bold magenta]{}_cmap_target.txt[/bold magenta]".format(output_prefix))

            if input_type == 'cmap':
                np.savetxt('{}_cmap_final.txt'.format(output_prefix), cmap_maxent)
                if verbose and console:
                    console.print(
                        "Final contact map saved to file: [bold magenta]{}_cmap_final.txt[/bold magenta]".format(output_prefix))
            
            np.savetxt('{}_connectivity_matrix.txt'.format(output_prefix), final_connectivity_matrix)
            if verbose and console:
                console.print(
                    'Connectivity matrix saved to file: [bold magenta]{}_connectivity_matrix.txt[/bold magenta]'.format(output_prefix))

            if not no_xyzs and xyzs is not None:
                write2xyz('{}.xyz'.format(output_prefix), xyzs, alignment=True, allow_reflection=True)
                if verbose and console:
                    console.print(
                        "Ensemble of (aligned) structures saved to file: [bold magenta]{}.xyz[/bold magenta]".format(output_prefix))
            
            if verbose and console:
                status.stop()
        except Exception as e:
            if verbose and console and hasattr(status, 'stop'):
                status.stop()
            raise e

    elapsed_seconds = time.perf_counter() - start_time
    if verbose and console:
        console.print(f"Total walltime: {elapsed_seconds:.2f} seconds")
    else:
        print(f"Total walltime: {elapsed_seconds:.2f} seconds")
    
    return results
