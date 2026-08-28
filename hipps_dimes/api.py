"""High-level programmatic API for HIPPS-DIMES."""

import json
import pickle
import time

import numpy as np
import pandas as pd
from rich import print
from tqdm import tqdm

from .covariance_pdhg import (
    fit_gaussian_noise_covariance_hybrid,
    fit_gaussian_noise_covariance_pdhg,
)
from .models import Optimize
from .numerics import *  # noqa: F401,F403
from .numerics import _optimize_contact_threshold

_COVARIANCE_PROGRESS_STAGES = {
    'covariance_operator_norm': ('COV operator norm', 'iteration'),
}


def _make_covariance_progress_callback(progress_callback, show_progress):
    """Render COV solver stages while preserving structured callbacks."""
    if not show_progress:
        return progress_callback, lambda: None

    progress_bar = None
    current_stage = None

    def close_progress_bar():
        nonlocal progress_bar
        if progress_bar is not None:
            progress_bar.close()
            progress_bar = None

    def report_progress(event):
        nonlocal current_stage, progress_bar
        if progress_callback is not None:
            progress_callback(event)

        stage = event.get('stage')
        phase = event.get('phase')
        if stage == 'covariance_optimization':
            description = (
                'COV PDHG optimization'
                if phase == 'pdhg'
                else 'COV FISTA optimization'
            )
            unit = 'iteration'
        elif stage in _COVARIANCE_PROGRESS_STAGES:
            description, unit = _COVARIANCE_PROGRESS_STAGES[stage]
        else:
            return

        stage_key = (stage, phase)
        if stage_key != current_stage:
            close_progress_bar()
            progress_bar = tqdm(
                total=int(event['total']),
                desc=description,
                unit=unit,
            )
            current_stage = stage_key

        if stage == 'covariance_operator_norm':
            progress_bar.set_postfix(
                relative_residual=(
                    f"{event['operator_norm_relative_residual']:.3e}"
                ),
                refresh=False,
            )
        else:
            postfix = {
                'objective': f"{event['objective']:.3e}",
                'relative_kkt': f"{event['relative_gradient_norm']:.3e}",
            }
            progress_bar.set_postfix(**postfix, refresh=False)

        completed = int(event['iteration'])
        progress_bar.update(max(0, completed - progress_bar.n))

    return report_progress, close_progress_bar


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


def _compute_missing_mask(matrix, matrix_kind):
    """Compute an off-diagonal missing-data mask for the given matrix kind."""
    matrix = np.asarray(matrix)
    if matrix_kind == 'cmap':
        missing_mask = ~np.isfinite(matrix) | (matrix <= 0.0)
    elif matrix_kind in {'dmap', 'ddmap'}:
        missing_mask = ~np.isfinite(matrix)
    else:
        raise ValueError(f"Unsupported matrix_kind '{matrix_kind}'")
    np.fill_diagonal(missing_mask, False)
    return missing_mask


def _summarize_missing_data(matrix, matrix_kind):
    """Summarize missing-data statistics for a square target matrix."""
    matrix = np.asarray(matrix)
    missing_mask = _compute_missing_mask(matrix, matrix_kind)
    n = matrix.shape[0]
    total_pairs = n * (n - 1) // 2
    missing_pairs = int(np.count_nonzero(np.triu(missing_mask, k=1)))
    missing_fraction = 0.0 if total_pairs == 0 else missing_pairs / total_pairs
    fully_missing_loci = np.where(np.count_nonzero(missing_mask, axis=1) == max(n - 1, 0))[0].tolist()
    return {
        'missing_mask': missing_mask,
        'missing_pairs': missing_pairs,
        'total_pairs': total_pairs,
        'missing_fraction': missing_fraction,
        'fully_missing_loci': fully_missing_loci,
    }


def _format_index_list(indices, limit=12):
    """Format a list of locus indices for console display."""
    if not indices:
        return "None"
    if len(indices) <= limit:
        return ", ".join(str(i) for i in indices)
    preview = ", ".join(str(i) for i in indices[:limit])
    return f"{preview}, ... ({len(indices)} total)"


def _repair_fully_missing_loci_nearest_neighbors(matrix, matrix_kind, missing_mask, fully_missing_loci):
    """Fill nearest-neighbor entries for fully missing loci using nearest interpolation."""
    repaired = np.array(matrix, dtype=float, copy=True)
    n = repaired.shape[0]

    if matrix_kind == 'cmap':
        with np.errstate(divide='ignore', invalid='ignore'):
            transformed = np.log10(repaired)
        transformed[missing_mask] = np.nan
        interpolated = interpolate_missing(transformed)
    elif matrix_kind in {'dmap', 'ddmap'}:
        transformed = np.array(repaired, dtype=float, copy=True)
        transformed[missing_mask] = np.nan
        interpolated = interpolate_missing(transformed)
    else:
        raise ValueError(f"Unsupported matrix_kind '{matrix_kind}'")

    repaired_pairs = set()
    for i in fully_missing_loci:
        for j in (i - 1, i + 1):
            if not (0 <= j < n):
                continue
            if not missing_mask[i, j]:
                continue

            candidates = [interpolated[i, j], interpolated[j, i]]
            finite_candidates = []
            for value in candidates:
                if not np.isfinite(value):
                    continue
                if matrix_kind == 'cmap':
                    finite_candidates.append(10.0 ** value)
                else:
                    finite_candidates.append(float(value))

            if not finite_candidates:
                raise ValueError(
                    f"Nearest-neighbor interpolation failed for fully missing locus {i} at pair ({i}, {j})"
                )

            fill_value = float(np.mean(finite_candidates))
            if matrix_kind == 'cmap' and fill_value <= 0.0:
                raise ValueError(
                    f"Nearest-neighbor interpolation produced a nonpositive contact value for pair ({i}, {j})"
                )

            repaired[i, j] = fill_value
            repaired[j, i] = fill_value
            repaired_pairs.add((min(i, j), max(i, j)))

    repaired_summary = _summarize_missing_data(repaired, matrix_kind)
    failed_pairs = [
        pair for pair in sorted(repaired_pairs)
        if repaired_summary['missing_mask'][pair[0], pair[1]]
    ]
    if failed_pairs:
        raise ValueError(
            f"Nearest-neighbor repair verification failed for pairs: {failed_pairs}"
        )

    return repaired, {
        'nearest_neighbor_repaired_pairs': sorted(repaired_pairs),
        'nearest_neighbor_repair_count': len(repaired_pairs),
        'remaining_fully_missing_loci': repaired_summary['fully_missing_loci'],
        'missing_pairs_after_repair': repaired_summary['missing_pairs'],
        'missing_fraction_after_repair': repaired_summary['missing_fraction'],
    }


def _remove_fully_missing_loci(matrix, matrix_kind):
    """Remove fully missing loci before optimization and report original-index mapping."""
    reduced = np.array(matrix, dtype=float, copy=True)
    original_size = reduced.shape[0]
    kept_loci = np.arange(original_size, dtype=int)
    removed_loci = []

    while True:
        summary = _summarize_missing_data(reduced, matrix_kind)
        current_fully_missing = summary['fully_missing_loci']
        if not current_fully_missing:
            break
        if reduced.shape[0] - len(current_fully_missing) < 2:
            raise ValueError(
                "Removing fully missing loci leaves fewer than 2 loci; optimization cannot continue"
            )

        keep_mask = np.ones(reduced.shape[0], dtype=bool)
        keep_mask[current_fully_missing] = False
        removed_loci.extend(kept_loci[~keep_mask].tolist())
        kept_loci = kept_loci[keep_mask]
        reduced = reduced[np.ix_(keep_mask, keep_mask)]

    reduced_summary = _summarize_missing_data(reduced, matrix_kind)
    return reduced, {
        'removed_fully_missing_loci': removed_loci,
        'removed_fully_missing_loci_count': len(removed_loci),
        'kept_loci': kept_loci.tolist(),
        'remaining_fully_missing_loci': reduced_summary['fully_missing_loci'],
        'missing_pairs_after_removal': reduced_summary['missing_pairs'],
        'missing_fraction_after_removal': reduced_summary['missing_fraction'],
    }


def _subset_connectivity_matrix(connectivity_matrix, kept_loci, original_size):
    """Subset a provided connectivity matrix to the kept loci when needed."""
    connectivity_matrix = _ensure_square_matrix(connectivity_matrix, "Connectivity matrix")
    optimized_size = len(kept_loci)
    if connectivity_matrix.shape[0] == original_size:
        kept_loci = np.asarray(kept_loci, dtype=int)
        subset = connectivity_matrix[np.ix_(kept_loci, kept_loci)]
        return a2a(subset), True
    if connectivity_matrix.shape[0] == optimized_size:
        return connectivity_matrix, False
    raise ValueError(
        "Connectivity matrix must match either the original input size "
        f"({original_size}x{original_size}) or the reduced optimization size "
        f"({optimized_size}x{optimized_size}); got {connectivity_matrix.shape}"
    )


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
                     gaussian_noise_relative_std=None,
                     covariance_optimizer='hybrid',
                     covariance_relative_tolerance=1e-5,
                     covariance_absolute_tolerance=1e-10,
                     covariance_handoff_relative_tolerance=1e-2,
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
                     remove_fully_missing_loci=False,
                     balance=False,
                     not_normalize=False,
                     neighbor_balance=False,
                     enforce_nonnegative_connectivity_matrix=False,
                     save_steps=None,
                     save_pickle=False,
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
        Optimization method: 'IS', 'GD', 'DI', or calibrated Gaussian 'COV'
    lamd : float, default=0.0
        Regularization weight
    reg : str, default='L2'
        Regularization type: 'L1' or 'L2'
    gaussian_noise_variance : float, default=0.0
        Positive homoskedastic variance on squared-distance constraints. COV only.
    gaussian_noise_relative_std : float, optional
        Positive shared relative standard deviation ``sigma_ij / Dobs_ij``.
        COV converts it to ``variance_ij = (value * Dobs_ij)**2`` after input
        conversion and missing-data handling.
    covariance_optimizer : {'hybrid', 'pdhg'}, default='hybrid'
        Optimizer used for the Gaussian COV objective. PDHG is
        variance-whitened and uses inverse-free runtime KKT diagnostics. The
        hybrid default uses PDHG globally and direct-Gram monotone FISTA for
        refinement.
    covariance_relative_tolerance : float, default=1e-5
        Relative COV KKT tolerance. Hybrid and PDHG results must pass a
        freshly recomputed dual-eliminated KKT certificate at this tolerance
        before the run is reported as converged.
    covariance_absolute_tolerance : float, default=1e-10
        Absolute tolerance used by the selected COV optimizer's internal KKT
        checks.
    covariance_handoff_relative_tolerance : float, default=1e-2
        Relative KKT threshold at which the hybrid optimizer switches from
        PDHG to FISTA. Ignored by standalone PDHG.
    iteration : int, default=10000
        Maximum optimization iterations. Hybrid PDHG and FISTA updates share
        this single total budget.
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
        If True, use CuPy GPU acceleration. All COV optimizers use float64
        throughout and fail clearly if no CUDA GPU is accessible. Legacy
        IS/GD retain their existing behavior.
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
    remove_fully_missing_loci : bool, default=False
        If True, and ``ignore_missing_data`` is also True, remove loci whose
        entire off-diagonal row/column is missing before optimization.
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
    save_pickle : bool, default=False
        If True, save the returned results dictionary to
        ``{output_prefix}_HIPPS_DIMES_results.pkl`` and suppress the default
        text/CSV/XYZ file outputs. Requires ``output_prefix``.
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
        Whether to render solver progress bars. COV reports the selected
        optimizer phases and operator-norm estimation separately. Set to False
        when consuming progress_callback programmatically.
        
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
        - 'kept_loci': Original locus indices retained for optimization (only if loci were removed)
        - 'removed_fully_missing_loci': Original fully missing loci removed before optimization
    
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
    - All COV optimizers run GPU matrix operations in float64
    - Hybrid COV hands the physical centered PDHG Gram directly to FISTA
    - For small matrices, CPU may be faster due to GPU setup overhead
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
    valid_methods = {'IS', 'GD', 'DI', 'COV'}
    if method not in valid_methods:
        raise ValueError(
            f"Invalid method '{method}'. Must be one of {sorted(valid_methods)}"
        )
    if isinstance(gaussian_noise_variance, (bool, np.bool_)) or not np.isscalar(
        gaussian_noise_variance
    ):
        raise ValueError("gaussian_noise_variance must be a nonnegative finite scalar")
    try:
        gaussian_noise_variance = float(gaussian_noise_variance)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "gaussian_noise_variance must be a nonnegative finite scalar"
        ) from error
    if not np.isfinite(gaussian_noise_variance) or gaussian_noise_variance < 0.0:
        raise ValueError("gaussian_noise_variance must be a nonnegative finite scalar")
    if gaussian_noise_relative_std is not None:
        if isinstance(gaussian_noise_relative_std, (bool, np.bool_)) or not np.isscalar(
            gaussian_noise_relative_std
        ):
            raise ValueError(
                "gaussian_noise_relative_std must be a positive finite scalar"
            )
        try:
            gaussian_noise_relative_std = float(gaussian_noise_relative_std)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "gaussian_noise_relative_std must be a positive finite scalar"
            ) from error
        if (
            not np.isfinite(gaussian_noise_relative_std)
            or gaussian_noise_relative_std <= 0.0
        ):
            raise ValueError(
                "gaussian_noise_relative_std must be a positive finite scalar"
            )
    if covariance_optimizer not in {'hybrid', 'pdhg'}:
        raise ValueError(
            "covariance_optimizer must be 'hybrid' or 'pdhg'"
        )
    covariance_tolerances = {
        'covariance_relative_tolerance': covariance_relative_tolerance,
        'covariance_absolute_tolerance': covariance_absolute_tolerance,
    }
    for tolerance_name, tolerance_value in covariance_tolerances.items():
        if isinstance(tolerance_value, (bool, np.bool_)) or not np.isscalar(
            tolerance_value
        ):
            raise ValueError(f"{tolerance_name} must be a nonnegative finite scalar")
        try:
            covariance_tolerances[tolerance_name] = float(tolerance_value)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{tolerance_name} must be a nonnegative finite scalar"
            ) from error
        if (
            not np.isfinite(covariance_tolerances[tolerance_name])
            or covariance_tolerances[tolerance_name] < 0.0
        ):
            raise ValueError(
                f"{tolerance_name} must be a nonnegative finite scalar"
            )
    covariance_relative_tolerance = covariance_tolerances[
        'covariance_relative_tolerance'
    ]
    covariance_absolute_tolerance = covariance_tolerances[
        'covariance_absolute_tolerance'
    ]
    if (
        covariance_relative_tolerance == 0.0
        and covariance_absolute_tolerance == 0.0
    ):
        raise ValueError("at least one covariance convergence tolerance must be positive")
    if (
        isinstance(covariance_handoff_relative_tolerance, (bool, np.bool_))
        or not np.isscalar(covariance_handoff_relative_tolerance)
    ):
        raise ValueError(
            "covariance_handoff_relative_tolerance must be a positive finite scalar"
        )
    try:
        covariance_handoff_relative_tolerance = float(
            covariance_handoff_relative_tolerance
        )
    except (TypeError, ValueError) as error:
        raise ValueError(
            "covariance_handoff_relative_tolerance must be a positive finite scalar"
        ) from error
    if (
        not np.isfinite(covariance_handoff_relative_tolerance)
        or covariance_handoff_relative_tolerance <= 0.0
    ):
        raise ValueError(
            "covariance_handoff_relative_tolerance must be a positive finite scalar"
        )
    if (
        covariance_optimizer == 'hybrid'
        and covariance_relative_tolerance > 0.0
        and covariance_handoff_relative_tolerance
        < covariance_relative_tolerance
    ):
        raise ValueError(
            "covariance_handoff_relative_tolerance must not be smaller than "
            "covariance_relative_tolerance"
        )

    has_absolute_noise = gaussian_noise_variance > 0.0
    has_relative_noise = gaussian_noise_relative_std is not None
    if method == 'COV':
        if has_absolute_noise == has_relative_noise:
            raise ValueError(
                "COV requires exactly one of gaussian_noise_variance or "
                "gaussian_noise_relative_std"
            )
        if lamd > 0.0:
            raise ValueError("COV cannot be combined with lamd regularization")
        if enforce_nonnegative_connectivity_matrix:
            raise ValueError(
                "COV cannot be combined with nonnegative-connectivity enforcement"
            )
        if gpu_float32:
            raise ValueError("COV supports float64 only; gpu_float32 is not allowed")
    else:
        if has_absolute_noise or has_relative_noise:
            raise ValueError(
                "Gaussian-noise options are supported only with method='COV'; "
                "legacy noisy IS/GD does not optimize the calibrated Gaussian objective"
            )
        if covariance_optimizer != 'hybrid':
            raise ValueError(
                "covariance_optimizer is supported only with method='COV'"
            )
        if (
            covariance_relative_tolerance != 1e-5
            or covariance_absolute_tolerance != 1e-10
        ):
            raise ValueError(
                "covariance convergence tolerances are supported only with "
                "method='COV'"
            )
        if covariance_handoff_relative_tolerance != 1e-2:
            raise ValueError(
                "covariance_handoff_relative_tolerance is supported only "
                "with method='COV'"
            )
    if log is not None:
        no_log = not log
    if save_pickle and output_prefix is None:
        raise ValueError("output_prefix must be provided when save_pickle=True")
    if remove_fully_missing_loci and not ignore_missing_data:
        raise ValueError("remove_fully_missing_loci=True requires ignore_missing_data=True")
    if eigh_threads is not None:
        set_eigh_num_threads(eigh_threads)

    input_source = input_path if input_path else "NumPy array"
    connectivity_matrix_source = "default initialization"
    if connectivity_matrix is not None:
        connectivity_matrix_source = connectivity_matrix if isinstance(connectivity_matrix, str) else "provided matrix"
    save_target_cmap = input_type == 'cmap' and input_format in {'cooler', 'hic'}
    missing_analysis = None
    original_matrix_rows = None
    original_matrix_cols = None
    
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
                dmap_target = _ensure_square_matrix(input_matrix, "Distance map")
            elif input_format == 'text':
                dmap_target = _ensure_square_matrix(
                    np.loadtxt(input_path), "Distance map"
                )
            elif input_format == 'npy':
                dmap_target = _ensure_square_matrix(
                    np.load(input_path), "Distance map"
                )
            elif input_format in {'cooler', 'hic'}:
                raise ValueError("input_type='dmap' only supports input_format='text' or 'npy' (or provide input_matrix)")
            else:
                raise ValueError(
                    f"Invalid input_format '{input_format}' for input_type='dmap'. "
                    "Supported: 'text' or 'npy' (or provide input_matrix)."
                )
            original_matrix_rows, original_matrix_cols = dmap_target.shape
            missing_analysis = _summarize_missing_data(dmap_target, 'dmap')
            handling_info = {
                'nearest_neighbor_repaired_pairs': [],
                'nearest_neighbor_repair_count': 0,
                'removed_fully_missing_loci': [],
                'removed_fully_missing_loci_count': 0,
                'remaining_fully_missing_loci': missing_analysis['fully_missing_loci'],
                'missing_pairs_after_repair': missing_analysis['missing_pairs'],
                'missing_fraction_after_repair': missing_analysis['missing_fraction'],
                'missing_pairs_after_removal': missing_analysis['missing_pairs'],
                'missing_fraction_after_removal': missing_analysis['missing_fraction'],
            }
            if ignore_missing_data and remove_fully_missing_loci and missing_analysis['fully_missing_loci']:
                dmap_target, removal_info = _remove_fully_missing_loci(dmap_target, 'dmap')
                handling_info.update(removal_info)
            elif ignore_missing_data and missing_analysis['fully_missing_loci']:
                dmap_target, repair_update = _repair_fully_missing_loci_nearest_neighbors(
                    dmap_target,
                    'dmap',
                    missing_analysis['missing_mask'],
                    missing_analysis['fully_missing_loci'],
                )
                handling_info.update(repair_update)
            missing_analysis.update(handling_info)
            ddmap_target = ((3. * np.pi) / 8.) * np.power(dmap_target, 2.)
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
            original_matrix_rows, original_matrix_cols = ddmap_target.shape
            missing_analysis = _summarize_missing_data(ddmap_target, 'ddmap')
            handling_info = {
                'nearest_neighbor_repaired_pairs': [],
                'nearest_neighbor_repair_count': 0,
                'removed_fully_missing_loci': [],
                'removed_fully_missing_loci_count': 0,
                'remaining_fully_missing_loci': missing_analysis['fully_missing_loci'],
                'missing_pairs_after_repair': missing_analysis['missing_pairs'],
                'missing_fraction_after_repair': missing_analysis['missing_fraction'],
                'missing_pairs_after_removal': missing_analysis['missing_pairs'],
                'missing_fraction_after_removal': missing_analysis['missing_fraction'],
            }
            if ignore_missing_data and remove_fully_missing_loci and missing_analysis['fully_missing_loci']:
                ddmap_target, removal_info = _remove_fully_missing_loci(ddmap_target, 'ddmap')
                handling_info.update(removal_info)
            elif ignore_missing_data and missing_analysis['fully_missing_loci']:
                ddmap_target, repair_update = _repair_fully_missing_loci_nearest_neighbors(
                    ddmap_target,
                    'ddmap',
                    missing_analysis['missing_mask'],
                    missing_analysis['fully_missing_loci'],
                )
                handling_info.update(repair_update)
            missing_analysis.update(handling_info)
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
            original_matrix_rows, original_matrix_cols = cmap.shape
            missing_analysis = _summarize_missing_data(cmap, 'cmap')
            handling_info = {
                'nearest_neighbor_repaired_pairs': [],
                'nearest_neighbor_repair_count': 0,
                'removed_fully_missing_loci': [],
                'removed_fully_missing_loci_count': 0,
                'remaining_fully_missing_loci': missing_analysis['fully_missing_loci'],
                'missing_pairs_after_repair': missing_analysis['missing_pairs'],
                'missing_fraction_after_repair': missing_analysis['missing_fraction'],
                'missing_pairs_after_removal': missing_analysis['missing_pairs'],
                'missing_fraction_after_removal': missing_analysis['missing_fraction'],
            }
            if ignore_missing_data and remove_fully_missing_loci and missing_analysis['fully_missing_loci']:
                cmap, removal_info = _remove_fully_missing_loci(cmap, 'cmap')
                handling_info.update(removal_info)
            elif ignore_missing_data and missing_analysis['fully_missing_loci']:
                cmap, repair_update = _repair_fully_missing_loci_nearest_neighbors(
                    cmap,
                    'cmap',
                    missing_analysis['missing_mask'],
                    missing_analysis['fully_missing_loci'],
                )
                handling_info.update(repair_update)
            missing_analysis.update(handling_info)
            
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

        # Legacy missing-data paths use infinity as their sentinel, while COV
        # uses NaN so finite entries alone define the observed pair set.
        if (
            method == 'COV'
            and ignore_missing_data
            and np.any(np.isinf(ddmap_target))
        ):
            ddmap_target = np.asarray(ddmap_target, dtype=np.float64).copy()
            ddmap_target[np.isinf(ddmap_target)] = np.nan

        # Load connectivity matrix if provided
        if connectivity_matrix is not None:
            if isinstance(connectivity_matrix, str):
                connectivity_matrix = np.loadtxt(connectivity_matrix)
            if missing_analysis is not None and missing_analysis['removed_fully_missing_loci_count'] > 0:
                connectivity_matrix, subset_applied = _subset_connectivity_matrix(
                    connectivity_matrix,
                    missing_analysis['kept_loci'],
                    original_matrix_rows,
                )
                if subset_applied:
                    connectivity_matrix_source = (
                        f"{connectivity_matrix_source} (subset to kept loci)"
                    )
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
        if missing_analysis is not None and missing_analysis['removed_fully_missing_loci_count'] > 0:
            input_table.add_row("Original Matrix Size", f"{original_matrix_rows} × {original_matrix_cols}")
            input_table.add_row("Optimized Matrix Size", f"{ddmap_target.shape[0]} × {ddmap_target.shape[1]}")
        else:
            input_table.add_row("Matrix Size", f"{ddmap_target.shape[0]} × {ddmap_target.shape[1]}")
        if input_type == 'cmap':
            input_table.add_row("Alpha (cmap→dmap)", f"{alpha}")
            input_table.add_row("Matrix Balancing", "Yes" if balance else "No" if input_format == 'cooler' else "N/A")
            input_table.add_row("Neighbor Balancing", "Yes" if neighbor_balance else "No")
            input_table.add_row("Auto Normalization", "No" if not_normalize else "Yes")
        input_table.add_row("Ignore Missing Data", "Yes" if ignore_missing_data else "No")
        input_table.add_row("Remove Fully Missing Loci", "Yes" if remove_fully_missing_loci else "No")
        
        console.print(input_table)
        console.print()  # spacing

        if missing_analysis is not None:
            missing_table = Table(title="Missing Data Analysis", show_header=True, header_style="bold cyan")
            missing_table.add_column("Parameter", style="dim", width=24)
            missing_table.add_column("Value", style="green")
            missing_table.add_row(
                "Missing Pairs",
                f"{missing_analysis['missing_pairs']:,} / {missing_analysis['total_pairs']:,} ({100.0 * missing_analysis['missing_fraction']:.2f}%)",
            )
            missing_table.add_row(
                "Fully Missing Loci",
                _format_index_list(missing_analysis['fully_missing_loci']),
            )
            if ignore_missing_data and remove_fully_missing_loci and missing_analysis['removed_fully_missing_loci_count'] > 0:
                missing_table.add_row(
                    "Removed Fully Missing",
                    _format_index_list(missing_analysis['removed_fully_missing_loci']),
                )
                missing_table.add_row(
                    "Removed Loci Count",
                    str(missing_analysis['removed_fully_missing_loci_count']),
                )
                missing_table.add_row(
                    "Remaining Fully Missing",
                    _format_index_list(missing_analysis['remaining_fully_missing_loci']),
                )
                missing_table.add_row(
                    "Missing Pairs After Removal",
                    f"{missing_analysis['missing_pairs_after_removal']:,} / {ddmap_target.shape[0] * (ddmap_target.shape[0] - 1) // 2:,} ({100.0 * missing_analysis['missing_fraction_after_removal']:.2f}%)",
                )
            elif ignore_missing_data and missing_analysis['fully_missing_loci']:
                missing_table.add_row(
                    "NN Repaired Pairs",
                    _format_index_list(
                        [f"{i}-{j}" for i, j in missing_analysis['nearest_neighbor_repaired_pairs']]
                    ),
                )
                missing_table.add_row(
                    "NN Repair Count",
                    str(missing_analysis['nearest_neighbor_repair_count']),
                )
                missing_table.add_row(
                    "Remaining Fully Missing",
                    _format_index_list(missing_analysis['remaining_fully_missing_loci']),
                )
                missing_table.add_row(
                    "Missing Pairs After Repair",
                    f"{missing_analysis['missing_pairs_after_repair']:,} / {missing_analysis['total_pairs']:,} ({100.0 * missing_analysis['missing_fraction_after_repair']:.2f}%)",
                )

            console.print(missing_table)
            console.print()  # spacing
        
        # Table 2: Optimization Settings
        opt_table = Table(title="Optimization Settings", show_header=True, header_style="bold cyan")
        opt_table.add_column("Parameter", style="dim", width=20)
        opt_table.add_column("Value", style="green")
        
        method_str = (
            "Iterative Scaling (IS)" if method == 'IS'
            else "Gradient Descent (GD)" if method == 'GD'
            else "Direct Inversion (DI)" if method == 'DI'
            else "Covariance-cone Gaussian fit (COV)" if method == 'COV'
            else "Unknown"
        )
        opt_table.add_row("Method", method_str)
        
        if method != 'DI':
            opt_table.add_row("Iterations", f"{iteration:,}")
            if method in {'IS', 'GD'}:
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
            if method == 'COV':
                opt_table.add_row("Optimizer", covariance_optimizer)
                opt_table.add_row(
                    "Initialization",
                    "provided connectivity"
                    if connectivity_matrix is not None
                    else "Rouse",
                )
                if has_absolute_noise:
                    opt_table.add_row(
                        "Noise Model", f"absolute variance {gaussian_noise_variance}"
                    )
                else:
                    opt_table.add_row(
                        "Noise Model",
                        f"relative std {gaussian_noise_relative_std}",
                    )
                opt_table.add_row(
                    "Relative KKT Tolerance",
                    f"{covariance_relative_tolerance:.3e}",
                )
                opt_table.add_row(
                    "Absolute KKT Tolerance",
                    f"{covariance_absolute_tolerance:.3e}",
                )
                if covariance_optimizer == 'hybrid':
                    opt_table.add_row(
                        "PDHG-to-FISTA Handoff",
                        f"{covariance_handoff_relative_tolerance:.3e}",
                    )
        
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
        output_table.add_row("Save Results Pickle", "Yes" if save_pickle else "No")
        if save_pickle:
            output_table.add_row("Default Output Files", "Disabled")
            output_table.add_row("Write XYZ File", "Stored in pickle only" if not no_xyzs else "No")
            output_table.add_row("Write Log Files", "Stored in pickle only")
        else:
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
    solver_output_prefix = None if save_pickle else output_prefix
    covariance_optimization_info = None
    fitted_gram_matrix = None
    iteration_extra_series = None
    use_gpu_enabled = False

    if method == 'COV':
        covariance_solver = {
            'hybrid': fit_gaussian_noise_covariance_hybrid,
            'pdhg': fit_gaussian_noise_covariance_pdhg,
        }[covariance_optimizer]
        covariance_progress_callback, close_covariance_progress = (
            _make_covariance_progress_callback(progress_callback, show_progress)
        )
        try:
            covariance_solver_arguments = {
                'relative_noise_std': gaussian_noise_relative_std,
                'initial_connectivity': connectivity_matrix,
                'use_gpu': use_gpu,
                'max_iterations': iteration,
                'relative_tolerance': covariance_relative_tolerance,
                'absolute_tolerance': covariance_absolute_tolerance,
                'save_steps': save_steps,
                'progress_callback': covariance_progress_callback,
            }
            if covariance_optimizer == 'hybrid':
                covariance_solver_arguments[
                    'handoff_relative_tolerance'
                ] = covariance_handoff_relative_tolerance
            (
                fitted_ddmap,
                fitted_gram_matrix,
                final_connectivity_matrix,
                covariance_optimization_info,
            ) = covariance_solver(
                ddmap_target,
                gaussian_noise_variance if has_absolute_noise else None,
                **covariance_solver_arguments,
            )
        finally:
            close_covariance_progress()
        dmap_maxent = a2dmap_theory(
            final_connectivity_matrix, force_positive_definite=True
        )
        reconstructed_ddmap = (3.0 * np.pi / 8.0) * np.square(dmap_maxent)
        consistency_scale = max(float(np.max(np.abs(fitted_ddmap))), 1.0)
        if not np.allclose(
            reconstructed_ddmap,
            fitted_ddmap,
            rtol=1e-8,
            atol=1e-10 * consistency_scale,
        ):
            raise RuntimeError(
                "COV connectivity does not reproduce its fitted squared-distance map"
            )
        history = covariance_optimization_info['history']
        loss = history['loss'].tolist()
        entropy = history['entropy'].tolist()
        iteration_extra_series = {
            key: values
            for key, values in history.items()
            if key not in {'iteration', 'loss', 'entropy'}
        }
        connectivity_at_steps = covariance_optimization_info[
            'connectivity_matrix_at_steps'
        ]
        use_gpu_enabled = covariance_optimization_info['backend'] == 'gpu'
        if solver_output_prefix is not None:
            for step, matrix in connectivity_at_steps.items():
                np.savetxt(
                    f'{solver_output_prefix}_connectivity_matrix_iter{step}.txt',
                    matrix,
                )
    else:
        model = Optimize(
            ddmap_target,
            connectivity_matrix=connectivity_matrix,
            use_gpu=use_gpu,
            gpu_float32=gpu_float32,
        )
        use_gpu_enabled = model.use_gpu
        if use_gpu and model.use_gpu and verbose:
            dtype_str = (
                "float32" if getattr(model, "gpu_float32", False) else "float64"
            )
            console.print(
                f"[green]GPU acceleration enabled ({get_gpu_name()}), "
                f"dtype={dtype_str}[/green]"
            )
        keyword_arguments = {
            'learning_rate': learning_rate,
            'lamd': lamd,
            'reg': reg,
            'method': method,
            'enforce_nonnegative_connectivity_matrix': (
                enforce_nonnegative_connectivity_matrix
            ),
            'momentum': momentum,
            'nesterov': nesterov,
        }
        general_method = 'optimization' if method in {'IS', 'GD'} else 'direct'
        loss, entropy, dmap_maxent, final_connectivity_matrix, connectivity_at_steps = model.run(
            iteration, general_method=general_method, save_steps=save_steps,
            output_prefix=solver_output_prefix, progress_callback=progress_callback, show_progress=show_progress,
            **keyword_arguments)
    
    # Format per-iteration scalar outputs.
    iteration_series_df = _build_iteration_series_frame(
        loss, entropy, iteration_extra_series
    )

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
        'matrix_rows_original': original_matrix_rows if original_matrix_rows is not None else ddmap_target.shape[0],
        'matrix_cols_original': original_matrix_cols if original_matrix_cols is not None else ddmap_target.shape[1],
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
        'gaussian_noise_relative_std': gaussian_noise_relative_std,
        'gaussian_noise_model': (
            covariance_optimization_info['noise_model']
            if covariance_optimization_info is not None
            else None
        ),
        'gaussian_noise_variance_minimum': (
            covariance_optimization_info['noise_variance_minimum']
            if covariance_optimization_info is not None
            else None
        ),
        'gaussian_noise_variance_median': (
            covariance_optimization_info['noise_variance_median']
            if covariance_optimization_info is not None
            else None
        ),
        'gaussian_noise_variance_maximum': (
            covariance_optimization_info['noise_variance_maximum']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_optimizer_requested': (
            covariance_optimizer if method == 'COV' else None
        ),
        'covariance_optimizer_resolved': (
            covariance_optimization_info['algorithm']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_relative_tolerance': (
            covariance_relative_tolerance if method == 'COV' else None
        ),
        'covariance_absolute_tolerance': (
            covariance_absolute_tolerance if method == 'COV' else None
        ),
        'covariance_handoff_relative_tolerance': (
            covariance_handoff_relative_tolerance
            if method == 'COV' and covariance_optimizer == 'hybrid'
            else None
        ),
        'covariance_phase_iterations': (
            covariance_optimization_info.get('phase_iterations')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_phase_wall_seconds': (
            covariance_optimization_info.get('phase_wall_seconds')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_converged': (
            covariance_optimization_info['converged']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_status': (
            covariance_optimization_info['status']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_relative_eliminated_kkt_residual': (
            covariance_optimization_info.get(
                'relative_eliminated_kkt_residual'
            )
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_independent_kkt_converged': (
            covariance_optimization_info.get('independent_kkt_converged')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_initialization_resolved': (
            covariance_optimization_info['initialization']['kind']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_initialization_backend': (
            covariance_optimization_info['initialization']['backend']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_initialization_wall_seconds': (
            covariance_optimization_info['initialization']['wall_seconds']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_backend': (
            covariance_optimization_info['backend']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_dtype': (
            covariance_optimization_info['dtype']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_gpu_device': (
            covariance_optimization_info['gpu_device']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_cupy_version': (
            covariance_optimization_info['cupy_version']
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_fista_initial_step_size': (
            covariance_optimization_info.get('initial_step_size')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_fista_final_step_size': (
            covariance_optimization_info.get('final_step_size')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_fista_backtracking_reductions': (
            covariance_optimization_info.get('backtracking_reductions')
            if covariance_optimization_info is not None
            else None
        ),
        'covariance_fista_restart_count': (
            covariance_optimization_info.get('restart_count')
            if covariance_optimization_info is not None
            else None
        ),
        'iteration': iteration,
        'learning_rate': learning_rate,
        'momentum': momentum,
        'nesterov': nesterov,
        'use_gpu_requested': use_gpu,
        'use_gpu_enabled': use_gpu_enabled,
        'gpu_float32': gpu_float32,
        'binsize': binsize,
        'hic_norm': hic_norm,
        'hic_unit': hic_unit,
        'no_log': no_log,
        'no_xyzs': no_xyzs,
        'ignore_missing_data': ignore_missing_data,
        'remove_fully_missing_loci': remove_fully_missing_loci,
        'balance': balance,
        'not_normalize': not_normalize,
        'neighbor_balance': neighbor_balance,
        'enforce_nonnegative_connectivity_matrix': enforce_nonnegative_connectivity_matrix,
        'save_target_cmap': save_target_cmap,
        'save_steps': save_steps or [],
        'save_pickle': save_pickle,
        'missing_pairs': missing_analysis['missing_pairs'] if missing_analysis is not None else 0,
        'missing_pair_fraction': missing_analysis['missing_fraction'] if missing_analysis is not None else 0.0,
        'fully_missing_loci': missing_analysis['fully_missing_loci'] if missing_analysis is not None else [],
        'removed_fully_missing_loci': missing_analysis['removed_fully_missing_loci'] if missing_analysis is not None else [],
        'removed_fully_missing_loci_count': missing_analysis['removed_fully_missing_loci_count'] if missing_analysis is not None else 0,
        'nearest_neighbor_repair_count': missing_analysis['nearest_neighbor_repair_count'] if missing_analysis is not None else 0,
        'nearest_neighbor_repaired_pairs': missing_analysis['nearest_neighbor_repaired_pairs'] if missing_analysis is not None else [],
        'missing_pairs_after_repair': missing_analysis['missing_pairs_after_repair'] if missing_analysis is not None else 0,
        'missing_pair_fraction_after_repair': missing_analysis['missing_fraction_after_repair'] if missing_analysis is not None else 0.0,
        'missing_pairs_after_removal': missing_analysis['missing_pairs_after_removal'] if missing_analysis is not None else 0,
        'missing_pair_fraction_after_removal': missing_analysis['missing_fraction_after_removal'] if missing_analysis is not None else 0.0,
        'remaining_fully_missing_loci': missing_analysis['remaining_fully_missing_loci'] if missing_analysis is not None else [],
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
    if covariance_optimization_info is not None:
        results['gram_matrix'] = fitted_gram_matrix
        results['covariance_optimization'] = covariance_optimization_info
    if connectivity_at_steps:
        results['connectivity_matrix_at_steps'] = connectivity_at_steps
    if missing_analysis is not None and missing_analysis['removed_fully_missing_loci_count'] > 0:
        results['kept_loci'] = missing_analysis['kept_loci']
        results['removed_fully_missing_loci'] = missing_analysis['removed_fully_missing_loci']
    
    # Compute contact map if input was contact map
    cmap_maxent = None
    rc_optimal = None
    if input_type == 'cmap':
        rc_optimal = _optimize_contact_threshold(dmap_maxent, cmap)
        if verbose and console:
            console.print('Optimized contact threshold distance: {}\n'.format(rc_optimal))
        cmap_maxent = dmap2cmap(dmap_maxent, rc_optimal)
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
            if save_pickle:
                pickle_path = f"{output_prefix}_HIPPS_DIMES_results.pkl"
                with open(pickle_path, 'wb') as fout:
                    pickle.dump(results, fout, protocol=pickle.HIGHEST_PROTOCOL)
                if verbose and console:
                    console.print(
                        f"Results pickle saved to file: [bold magenta]{pickle_path}[/bold magenta]"
                    )
            else:
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
