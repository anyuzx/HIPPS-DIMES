![movie](https://raw.githubusercontent.com/anyuzx/HIPPS-DIMES/main/doc/source/movie.gif)

[**Installation**](#install) | [**Quick Start**](#quick-start) | [**API**](#api) | [**Chromatin Dynamics**](#dynamics-prediction-functionality) | [**Chromatin Mechanics**](#linear-mechanical-response)

HIPPS-DIMES is a Python implementation of a maximum-entropy method for
constructing ensembles of <ins>3D chromatin structures</ins> from experimentally
measured contact maps or pairwise distance constraints.[^1][^2][^3] It accepts
a mean distance map, a mean-squared-distance map, or a Hi-C contact map
(converted internally into distance constraints) and generates conformations
whose ensemble statistics fit those constraints under the selected optimization
method.

In addition to reconstructing static 3D chromatin structures, the model predicts
**chromatin-locus dynamics** and **chromatin-locus mechanics** using polymer
physics and the Ornstein–Uhlenbeck process. Available observables include
autocorrelation functions (ACFs), mean-square displacements (MSDs), system-level
viscoelastic moduli, and per-locus mechanical susceptibilities.

![schematic](https://raw.githubusercontent.com/anyuzx/HIPPS-DIMES/main/doc/source/schematic.jpg)

The theory and applications of this method are described in the following
publications:

- Shi, Guang, and D. Thirumalai. "Epigenetic state encodes locus-specific chromatin mechanics." bioRxiv 2025.12.27.696709 (2025). [link](https://www.biorxiv.org/content/10.64898/2025.12.27.696709v1.abstract)
- Shi, Guang, and D. Thirumalai. "From Hi-C contact map to three-dimensional organization of interphase human chromosomes." Physical Review X 11.1 (2021): 011051. [link](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051)
- Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants." Nature Communications 14.1 (2023): 1150. [link](https://www.nature.com/articles/s41467-023-36412-4)
- Shi, Guang, Sucheol Shin, and D. Thirumalai. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advances 11.31 (2025): eadx1763. [link](https://www.science.org/doi/full/10.1126/sciadv.adx1763)

Other applications of this method can be found in:

- Dey, Atreya, et al. "Structural changes in chromosomes driven by multiple condensin motors during mitosis." Cell Reports 42.4 (2023).
- Jeong, Davin, et al. "Structural basis for the preservation of a subset of topologically associating domains in interphase chromosomes upon cohesin depletion." eLife 12 (2024): RP88564.

## Release 4.0

HIPPS-DIMES 4.0 adds the calibrated, noise-aware `COV` method with CPU and
float64 GPU execution. Its default hybrid optimizer combines
variance-whitened PDHG with direct-Gram monotone FISTA and reports an
independently recomputed KKT convergence certificate. This release also makes
fully missing-locus repair explicit, validates that COV observations connect
all retained loci, and makes the returned COV iterate unambiguous in logs and
results. Missing-pair handling is now consistent across contact, distance, and
squared-distance inputs: pairs are either interpolated before optimization or
explicitly excluded.

# Documentation

## Requirements

- Python 3.9+

## Install

First, clone this repository:

```bash
git clone https://github.com/anyuzx/HIPPS-DIMES
cd HIPPS-DIMES
```

### Using pip

```bash
pip install -e .
```

### Using uv

```bash
uv pip install -e .
```

Either command installs the dependencies declared in `pyproject.toml`, installs
HIPPS-DIMES as a Python package, and provides the `HippsDimes` command-line
entry point.

### Dependencies

The package requires:

- Python 3.9+
- `Click` - Command-line interface
- `rich-click` - Rich-styled command-line interface
- `NumPy` - Numerical computing
- `SciPy` - Scientific computing
- `pandas` - Data manipulation
- `tqdm` - Progress bars
- `cooler` - `.mcool` and `.cool` data-format support
- `Rich` - Rich terminal output
- `hic-straw` - `.hic` data-format support

**Optional (for GPU acceleration):**

- `CuPy` - GPU-accelerated computing via CUDA

To install CuPy for GPU support:

```bash
conda install -c conda-forge cupy
# or
pip install cupy-cuda11x  # Replace with your CUDA version
```

## Quick start

The quickest way to get started is with one of the Google Colab notebooks:

- [Quick start](https://colab.research.google.com/drive/1w7cK6S3z2_D5Mzgq2-SZgl4FoDjT9EC5?usp=sharing)
- [Chromatin Dynamics and Mechanics](https://colab.research.google.com/drive/1DsTIJTkiKc1vRq6lBncEg4nvPh4hQ4hp?usp=sharing)

To view the available arguments and options, run:

```bash
HippsDimes --help
```

## API

### Input files

HIPPS-DIMES accepts contact maps (`cmap`), mean distance maps (`dmap`), and
mean-squared-distance maps (`ddmap`). Contact maps may use `.cool` or `.mcool`
format through [cooler](https://github.com/open2c/cooler), `.hic` format through
hicstraw, a plain-text matrix, or a NumPy `.npy` matrix. Mean distance and
mean-squared-distance maps may use plain text or NumPy `.npy` format.

In a plain-text matrix, each line contains one space-separated row. For
example, a mean distance map could look like:

```text
0  2  3
2  0  2
3  2  0
```

### Output files

This script will generate several files:

- A text file containing the final model mean distance map.
- For contact-map inputs, a text file containing the model contact map at the
  optimized contact threshold. This threshold maximizes agreement with the
  normalized experimental contact map.
- For `cooler` and `.hic` contact-map inputs, a text file containing the
  internal target contact map used to construct the distance constraints:
  `{output_prefix}_cmap_target.txt`.
- A text file containing the connectivity matrix.
- An `.xyz` file containing the generated ensemble (optional).
- A CSV file containing run parameters:
  `{output_prefix}_run_parameters.csv` (optional).
- A CSV file containing iteration-series scalar data:
  `{output_prefix}_iteration_series.csv` (optional).

### Explanation of the arguments and options

#### Arguments

- `INPUT`: Path to a contact map, mean distance map, or mean-squared-distance
  map.
- `OUTPUT_PREFIX`: Prefix for output files. For example, specifying `TEST`
  makes the output filenames start with `TEST_`.

#### Options

- `-k, --connectivity-matrix`: Path to a connectivity matrix to use as the
  initialization, for example when restarting from a previous result.
- `-e, --ensemble`: Number of individual conformations to generate from the
  fitted model. Default: 1000.
- `-a, --alpha`: Contact-to-distance conversion exponent. For contact-map
  input, HIPPS-DIMES uses $d_{ij} \propto c_{ij}^{-1/\alpha}$. The default
  value, $\alpha=4.0$, was estimated in
  [this work](https://doi.org/10.1126/science.aaf8084). Default: 4.0.
- `-s, --selection`: Specify chromosome or region. This option is required
  when the input file has `cooler` or `.hic` format. For cooler files, the value is passed to the `cooler.Cooler.matrix().fetch()` method. For .hic files, use format "chr1:start1-end1,chr2:start2-end2". For details on cooler selectors, please refer to their [documentation](https://cooler.readthedocs.io/en/latest/concepts.html#matrix-selector).
- `-m, --method`: Select IS (Iterative Scaling, default), GD (Gradient Descent), DI (Direct Inversion), or COV (calibrated Gaussian covariance-cone optimization).
- `-l, --lamd`: L1 or L2 regularization weight. A value of `0.0` disables
  regularization. Default: `0.0`.
- `-r, --reg`: Regularization type: `L1` or `L2` (default). Use this option
  together with `--lamd`.
- `--gaussian-noise-variance`: Positive scalar absolute variance on squared-distance constraints. COV only.
- `--gaussian-noise-relative-std`: Positive scalar relative standard deviation `sigma_ij / Dobs_ij`. COV converts it to pair variance `(value * Dobs_ij)^2` after preprocessing. Mutually exclusive with `--gaussian-noise-variance`.
- `--covariance-optimizer`: COV optimizer: `hybrid` (default) or `pdhg`.
  PDHG is variance-whitened and inverse-free during its runtime KKT checks.
  The hybrid runs PDHG globally, then refines the same physical Gram matrix
  with direct-Gram monotone FISTA. Both phases support CPU and float64 GPU
  execution.
- `--covariance-relative-tolerance`: Relative COV KKT tolerance. Default: `1e-5`.
- `--covariance-absolute-tolerance`: Absolute internal COV KKT tolerance. Default: `1e-10`.
- `--covariance-handoff-relative-tolerance`: Relative KKT threshold for the
  default PDHG-to-FISTA handoff. Default: `1e-2`.
- `-i, --iteration`: Maximum optimizer iterations. Default: 10000.
- `--learning-rate`: Learning rate for IS or GD. Typical IS values are 1–30;
  GD generally requires a much smaller value, such as `1e-8`. Default: `10.0`.
- `--momentum`: Momentum coefficient for IS, between `0.0` and `1.0`.
  **Recommended: use `0.95` with `--nesterov` for the fastest convergence
  observed in benchmarks.** Use `0.9` for a more conservative setting.
  Default: `0.0`.
- `--nesterov`: Use Nesterov Accelerated Gradient with IS. Recommended with
  `--momentum 0.95`.
- `--use-gpu`: Enable GPU acceleration through CuPy. All COV optimizers use
  `float64`. COV fails rather than silently falling back when CUDA is
  unavailable.
- `--gpu-float32`: Use `float32` for legacy GPU IS/GD. COV is `float64`-only.
- `--save-steps`: Comma-separated list of iteration steps at which to save the connectivity matrix. Example: `--save-steps 1000,5000,10000`. Files are saved as `{output_prefix}_connectivity_matrix_iter{step}.txt`. When used as a library (without `output_prefix`), connectivity matrices at these steps are still returned in `results['connectivity_matrix_at_steps']`.
- `--eigh-threads`: Number of eigendecomposition and BLAS/LAPACK threads. If
  unset, the backend default is used. Set to `1` for single-threaded runs.
- `--input-type`: Required input type: `cmap` (contact map), `dmap` (mean
  distance map), or `ddmap` (mean-squared-distance map).
- `--input-format`: Required input format: `text`, `npy`, `cooler`, or `hic`.
  Contact maps support all four formats; `dmap` and `ddmap` support `text` and
  `npy`.
- `--binsize`: Bin size for `.hic` input, in bp. Default: `25000`.
- `--norm`: `.hic` normalization: `KR`, `VC`, or `NONE`. Default: `KR`.
- `--unit`: `.hic` unit: `BP` or `FRAG`. Default: `BP`.
- `--no-log`: By default, the program writes two log files when `output-prefix` is provided: `{output_prefix}_run_parameters.csv` and `{output_prefix}_iteration_series.csv`. Use `--no-log` to disable writing both files.
- `--no-xyzs`: Do not write generated conformations to an `.xyz` file.
- `--ignore-missing-data`: Exclude remaining missing pair constraints. Without
  this flag, every remaining missing pair is interpolated before optimization:
  contact maps in log-contact space, distance maps in distance space, and
  squared-distance maps in squared-distance space. The completed map is the
  optimization target.
- `--repair-fully-missing-loci`: Impute only the genomic nearest-neighbor
  constraints `(i, i-1)` and `(i, i+1)` needed to reconnect a locus with no
  observed off-diagonal pairs. The repaired values become ordinary target
  constraints in the selected objective and, for COV, in its noise model.
- `--remove-fully-missing-loci`: Remove loci that have no observed
  off-diagonal pairs before applying the remaining missing-pair policy.
- `--balance`: Balance a cooler-format contact map before optimization.
- `--neighbor-balance`: Apply neighbor balancing to a contact map by dividing
  each pair value by the geometric mean of the corresponding neighbor-contact
  values. See Paggi and Zhang (2025) for details.
- `--not-normalize`: Disable automatic maximum-value normalization of a contact
  map.
- `--enforce-nonnegative-connectivity-matrix`: Constrain all off-diagonal
  spring constants to be nonnegative.
- `--save-pickle`: Save the returned results dictionary as
  `{output_prefix}_HIPPS_DIMES_results.pkl` instead of writing the default
  text, CSV, and XYZ outputs.
- `-q, --quiet`: Disable table output while retaining the progress bar.

#### COV quick start

COV is the noise-aware optimization method. Specify exactly one noise model:

- `--gaussian-noise-variance <v>` for homoskedastic absolute variance
  (`v_ij = v`).
- `--gaussian-noise-relative-std <w>` for heteroskedastic relative noise
  (`sigma_ij = w Dobs_ij`).

```bash
python -m hipps_dimes observed_ddmap.npy cov_fit \
  --input-type ddmap --input-format npy \
  --method COV --gaussian-noise-relative-std 0.1 \
  --use-gpu --iteration 10000 --no-xyzs --save-pickle
```

The default hybrid optimizer starts from a calibrated Rouse model, runs
variance-whitened PDHG, and hands the exact physical Gram matrix to monotone
FISTA for direct-Gram refinement. It stops
automatically when the built-in KKT convergence criterion is met or the
requested iteration limit is reached. Standalone PDHG remains available for
controlled comparisons or cases where its inverse-free iterations are
preferred:

```bash
python -m hipps_dimes observed_ddmap.npy cov_fit \
  --input-type ddmap --input-format npy \
  --method COV --gaussian-noise-relative-std 0.1 \
  --covariance-optimizer pdhg \
  --use-gpu --iteration 20000 --no-xyzs --save-pickle
```

An iteration limit is only a budget. A returned COV model is converged only
when its independently recomputed KKT certificate passes the requested
tolerance. The command-line program retains requested output files but exits
with status 1 when that certificate fails; the Python API emits a
`RuntimeWarning` and returns the partial result with
`results["covariance_optimization"]["converged"] == False`.

COV constructs an undirected observation graph from finite, positive
off-diagonal squared-distance constraints. Every retained locus must belong to
one connected component. Sparse observations are supported, but disconnected
clusters are rejected because their relative motion is unconstrained and the
noise-aware maximum-entropy objective has no finite optimum. When a locus has
no observations, choose exactly one explicit policy with
`--repair-fully-missing-loci` or `--remove-fully-missing-loci`, regardless of
the `--ignore-missing-data` setting. Repair or removal is applied first. The
remaining partial pairs are then either excluded with `--ignore-missing-data`
or interpolated without it. `DI` requires a complete target and therefore
cannot be combined with excluded pairs.

See [Noise-aware covariance optimization](doc/NOTE_ON_NOISE.md) for the noise
model, initialization, and scientific interpretation. See
[Variance-whitened PDHG and hybrid COV optimization](doc/COV_PDHG.md) for
solver details, convergence diagnostics, large-system guidance, and advanced
tuning.

### Examples

#### Example 1

First, download a cooler-format Hi-C contact map from
[here](https://drive.google.com/file/d/1eIxGv1JbIrEAVoUSQK_O_ebIjWo6toTJ/view?usp=sharing)
(**the file size is about 116 MB**). This map represents a mitotic chicken
chromosome and was originally retrieved from the
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740).
Rename it `hic_example.cool`, and then run:

```bash
HippsDimes hic_example.cool test --input-type cmap --input-format cooler -s chr7:10M-15M -i 10 -e 10
```

This command loads `hic_example.cool` and runs Iterative Scaling. The `test`
argument makes output filenames start with `test_`; `--input-type cmap`
identifies the input as a contact map; `--input-format cooler` identifies the
file format; and `-s chr7:10M-15M` selects the 10–15 Mb region of chromosome 7.
The input type, input format, and cooler region selection are required here.
Run `HippsDimes --help` for the complete option requirements.

With these defaults, the program writes `test.xyz`,
`test_connectivity_matrix.txt`, `test_dmap_final.txt`, `test_cmap_final.txt`,
`test_cmap_target.txt`, `test_run_parameters.csv`, and
`test_iteration_series.csv`. The `test.xyz` file contains 10 conformations and
can be viewed with VMD or other compatible visualization software.

#### Example 2

This example uses a Hi-C contact map for chromosome 14 in HeLa cells, measured
12 hours after release from prometaphase. Download the `.cool` file from
[here](https://drive.google.com/file/d/1j-zfDUP6LOZGCxz9uA3LaMI372ct1cU_/view?usp=sharing)
which was originally retrieved from the
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740)
under accession GSE102740. **The file size is about 655 MB.** Once downloaded,
run:

```bash
HippsDimes GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool::6 test --input-type cmap --input-format cooler -s chr14:20M-107M -i 10000 -e 10
```

This command loads group 6 of the multiresolution cooler file and runs
HIPPS-DIMES for 10,000 iterations. On an AMD Ryzen 5 3600 CPU, this example
takes approximately 3–4 minutes.

#### Example 3

A `.hic` example:

```bash
HippsDimes mydata.hic test \
  --input-type cmap --input-format hic \
  --selection chr1:31000000-41000000,chr1:31000000-41000000 \
  --binsize 25000 --norm KR --unit BP \
  -i 10000 -e 10
```

> For an example of applying HIPPS-DIMES directly to imaging data, see the
> corresponding notebook.

### Tips

- Computational cost and memory use grow rapidly with matrix size, but there
  is no universal matrix-size cutoff for convergence. The practical limit
  depends on the selected optimizer, hardware, number of observed pairs, and
  iteration budget. If a problem is too expensive, coarse-grain the input or
  analyze a smaller genomic region.
- Optimization tuning is method-specific. `--learning-rate` applies to IS and
  GD, while `--momentum` and `--nesterov` apply only to IS. Typical IS learning
  rates are between 1 and 30; GD generally requires a much smaller value, such
  as `1e-8`. These parameters are not used by DI or COV. Nesterov acceleration,
  for example `--momentum 0.95 --nesterov`, can accelerate IS, but the benefit
  is problem-dependent.
- For larger problems, consider `--use-gpu` when CuPy and an accessible CUDA
  GPU are available. The speedup depends on the matrix, optimizer, GPU, and CPU,
  so benchmark a representative case rather than assuming a fixed acceleration
  factor. COV uses float64 on the GPU and does not silently fall back to the CPU.
- Use `--save-steps 1000,5000,10000` to retain intermediate connectivity
  matrices. For IS, these checkpoints support the manual entropy-and-loss
  stopping decision described below. For COV, the built-in KKT test determines
  convergence; checkpoints are useful for diagnostics or restarting but do not
  replace the reported convergence certificate.
- Missing pairs are interpolated by default and the completed map becomes the
  target. Use `--ignore-missing-data` to exclude them instead. Non-finite input
  pairs are missing for every map type; nonpositive off-diagonal pairs are also
  missing for contact maps. If a locus has no observed off-diagonal pairs at
  all, select either `--repair-fully-missing-loci` to add its genomic-neighbor
  constraints or `--remove-fully-missing-loci` to remove it before the pair
  policy is applied. Repair and removal are mutually exclusive. COV
  additionally requires the retained observed-pair graph to be connected when
  pairs are excluded, including cases with no individually isolated locus.
- Contact-map inputs are normalized by their maximum entry by default. Use
  `--not-normalize` when the supplied contact map should retain its existing
  scale.
- A contact map alone does not define a physical length scale, so its inferred
  distances and structures are dimensionless. An external distance measurement,
  such as the mean distance between neighboring loci, can be used to rescale
  them.

### How to choose the optimization method, iteration count, and convergence criterion

The recommended method and stopping rule depend on the input and on whether
its uncertainty should be modeled:

| Input and interpretation | Recommended method | How to stop |
| --- | --- | --- |
| Experimental contact map, with uncertainty | COV with a user-specified noise model and noise level | Automatic KKT convergence test; `--iteration` is the maximum update budget |
| Experimental contact map, uncertainty unknown | IS is a practical alternative when only static 3D structures are needed | Manual judgment from the entropy and loss histories |
| Valid mean-squared-distance map, treated as exact | DI (preferred) or IS | DI is a one-step conversion; IS requires enough iterations to reproduce the target |
| Valid mean-squared-distance map, with uncertainty | COV with a user-specified noise model and noise level | Automatic KKT convergence test; `--iteration` is the maximum update budget |

#### Experimental contact maps

Converting an experimental contact map to pairwise distances does not guarantee
that the resulting mean-squared-distance map is a valid Euclidean distance
matrix. Consequently, the exact hard-constraint problem approached by
Iterative Scaling (IS) may have no feasible global solution: no Gaussian
ensemble may satisfy all inferred pair distances simultaneously.

If only static 3D structures are needed, IS remains a reasonable practical
choice. IS has no built-in convergence criterion and always runs for the number
of iterations requested by the user. The iteration-series file records both
`entropy` and `loss`; here, `loss` is the root-mean-square relative difference
between the model and target mean-squared distances over the constrained pairs.
As a default, select the iteration with the highest recorded entropy. In the
rarer case where the entropy curve has an early or local maximum while the loss
is still changing materially, do not treat that peak as convergence. Continue
until the loss begins to level off, and select a checkpoint using both curves.
This is necessarily an educated user judgment rather than a convergence
certificate supplied by HIPPS-DIMES.

> **Warning — IS with experimental contact maps:** Because the inferred target
> may have no feasible global solution, the high-frequency relaxation modes
> continue to drift as more IS iterations are performed. Interpret these modes
> with caution. This sensitivity is most relevant to the high-frequency regime
> of the predicted loss modulus, $G''(\omega)$. The low-frequency modes converge
> and, once converged, do not depend on the selected number of IS iterations.

For experimental contact maps, the preferred method is COV with an uncertainty
level supplied by the user. Either of the implemented noise models may be used:

- `--gaussian-noise-variance` specifies one homoskedastic absolute variance
  shared by all observed mean-squared-distance pairs.
- `--gaussian-noise-relative-std` specifies one shared relative standard
  deviation, `sigma_ij / Dobs_ij`; the resulting pair variances are
  heteroskedastic.

COV has a well-defined objective and built-in convergence tests. Its sole PDHG
implementation uses variance whitening and inverse-free runtime KKT checks.
The default hybrid PDHG-to-FISTA solver switches at relative KKT `1e-2`, uses
the PDHG Gram matrix directly without scalar recalibration, and stops
automatically only after the independently recomputed KKT residual meets the
requested tolerance. Direct-Gram FISTA uses dense centered matrices and a
matrix-free observed-pair operator, so its storage grows quadratically rather
than through a higher-order linear-system workspace. For COV, `--iteration` sets a
maximum update budget rather than a recommended stopping iteration. If the
budget is exhausted first, the returned convergence status is false and the
user should not describe the result as converged. Inspect
`results['covariance_optimization']['converged']`, `status`, and
`relative_eliminated_kkt_residual` for the final certificate. `iterations`
records the number of updates executed, while `returned_iteration` identifies
the model actually returned. The iteration-series column
`is_returned_iterate` marks that row, and the top-level COV `objective`, `loss`,
and `entropy` all describe the same returned model.

#### Valid mean-squared-distance maps

A mean-squared-distance map computed directly from an ensemble of 3D
coordinates is a valid Euclidean distance matrix, up to numerical precision.
If it is to be treated as exact, use Direct Inversion (`--method DI`). DI
performs the covariance-to-connectivity conversion in one step, so no iteration
count or iterative convergence decision is needed.

For a complete valid target, IS with suitable numerical settings converges to
the same solution as DI, but the required number of iterations depends on the
system size and optimization settings. DI is therefore preferred when the map
is valid and uncertainty is intentionally ignored. If uncertainty in the
mean-squared distances should be represented, use COV instead and specify
either the homoskedastic or heteroskedastic noise model and its level.

### Use HIPPS-DIMES as a Python library

HIPPS-DIMES can be used as either a command-line tool or a Python library,
making it straightforward to integrate into Python workflows, Jupyter
notebooks, and automated pipelines.

#### Main function: `run_optimization()`

The core functionality is available through the `run_optimization()` function:

```python
import numpy as np
import hipps_dimes as HD

# Load your contact map
cmap = np.loadtxt("contact_map.txt")

# Run optimization programmatically
results = HD.run_optimization(
    input_matrix=cmap,          # Provide matrix directly
    input_type="cmap",
    method="IS",
    iteration=10000,
    learning_rate=10.0,
    momentum=0.95,              # Use momentum for faster convergence
    nesterov=True,              # Enable Nesterov acceleration
    use_gpu=True,               # Use GPU if CuPy is available
    ensemble=1000,
    verbose=False,              # Suppress console output
)

# Access results
connectivity_matrix = results["connectivity_matrix"]
structures = results["xyzs"]          # (ensemble, n_beads, 3)
final_dmap = results["dmap_final"]
final_cmap = results["cmap_final"]
iteration_series = results["iteration_series"]
run_parameters = results["run_parameters"]
```

#### Return values

The `run_optimization()` function returns a dictionary with:

- `'connectivity_matrix'`: Final connectivity matrix (NumPy array)
- `'dmap_final'`: Final distance map (NumPy array)
- `'cmap_final'`: Final contact map (NumPy array, for `input_type="cmap"`)
- `'xyzs'`: Generated conformations (NumPy array, unless `no_xyzs=True`)
- `'iteration_series'`: Iteration-series scalar outputs (pandas DataFrame).
  Every method reports `iteration`, `loss`, and `entropy`; COV adds
  optimizer-specific diagnostics.
- `'run_parameters'`: Run parameters (pandas DataFrame with `parameter` and
  `value` columns)
- `'log'`: Alias for `'iteration_series'` (backward compatibility)
- `'rc_optimal'`: Optimal contact threshold (float, for `input_type="cmap"`)
- `'connectivity_matrix_at_steps'`: Saved intermediate connectivity matrices
  (dictionary, when `save_steps` is set)
- `'gram_matrix'` and `'covariance_optimization'`: Fitted Gram matrix and
  convergence diagnostics (for `method="COV"`)

#### Additional utility functions

The package also provides helper functions for direct use:

```python
import hipps_dimes as HD

# Generate structures from connectivity matrix
structures = HD.a2xyz_sample(connectivity_matrix, ensemble=1000)

# Convert connectivity matrix to distance map
dmap = HD.a2dmap_theory(connectivity_matrix)

# Convert connectivity matrix to contact map
cmap = HD.a2cmap_theory(connectivity_matrix, rc=5.0)

# Create a Rouse chain connectivity matrix
A = HD.construct_connectivity_matrix_rouse(n=100, k=1.0)
```

---

## Dynamics Prediction Functionality

In addition to reconstructing static 3D chromatin structures from contact or
distance maps, HIPPS-DIMES can simulate chromatin dynamics.[^3] The dynamics
are based on polymer physics and the Ornstein–Uhlenbeck process and provide
time-dependent observables such as autocorrelation functions (ACFs) and the
mean-square displacements (MSDs) of individual loci.

### Functions

- **`compute_acf_general_theory(i, j, t, a, zeta=1.0)`**: Numerically computes
  the time-dependent autocorrelation function between monomers *i* and *j* from
  a connectivity matrix `a`. It also returns the corresponding two-point MSD for
  every time in `t`.

- **`compute_m1_i(i, t, a, zeta=1.0)`**: Computes the single-locus MSD for
  monomer *i*. The returned two-dimensional array contains time in the first
  column and MSD in the second.

### `Dynamics` class

The `Dynamics` class provides trajectory simulation from a connectivity matrix
`a`.

#### Example code

```python
import hipps_dimes as HD

model = HD.Dynamics(a)  # a is the connectivity matrix
model.initialize(dt=1e-2, zeta=1.0, beta=1.0)

model.run(int(1e5), every=10)
model.resume(int(5e4), every=10)
```

Trajectory coordinates are available in `model.traj`, a `(T, N, 3)` NumPy
array, where `T` is the number of snapshots and `N` is the number of loci.
The reduced simulation time for each saved snapshot is stored in
`model.traj_time`, a length-`T` NumPy array.

Save both arrays with `model.save_traj("traj.npz")`. The `.npz` file contains
the `traj` and `traj_time` arrays.

`Dynamics.run(...)` starts a fresh trajectory and can be called only once per
simulation state. To continue an existing simulation, use
`Dynamics.resume(...)`. When arguments are omitted, `resume(...)` reuses the
previous passive simulation settings for `update`, `every`, `method`, and
`update_zero_modes`; any of them may still be overridden explicitly. To discard
the previous trajectory and start over on the same object, call
`Dynamics.reset()` before `run(...)`. By default, `run(...)` does not append the
post-integration final state to `model.traj`; set `include_final_state=True` to
include it.

### Dynamics under external force: `Dynamics.run_with_force(...)`

In addition to passive dynamics (`Dynamics.run`), you can simulate trajectories
with a constant external force applied to selected loci.

#### Key parameters

- `force_loci`: list of locus indices where the force is applied
- `force_amplitude`: force magnitude
- `force_direction`: `(3,)` direction vector (normalized internally)
- `force_duration`: optional number of time steps to apply the force; if `None`,
  the force is applied for the entire run

#### Example: forced dynamics

```python
import numpy as np
import hipps_dimes as HD

# Load a connectivity matrix or obtain one from HD.run_optimization().
a = np.loadtxt("my_connectivity_matrix.txt")
model = HD.Dynamics(a)
model.initialize(dt=1e-2, zeta=1.0, beta=1.0)

# Pull locus 10 along +x for the first 2e4 steps (then release)
model.run_with_force(
    T=int(1e5),
    force_loci=[10],
    force_amplitude=1.0,
    force_direction=[1.0, 0.0, 0.0],
    force_duration=int(2e4),
    every=10,
)

traj = model.traj  # shape: (n_snapshots, N, 3)
```

## Linear mechanical response

HIPPS-DIMES provides utilities to compute system-level **linear viscoelastic
moduli** and per-locus **mechanical susceptibilities** from a connectivity matrix
`a`. These routines decompose the polymer into normal modes and exclude the
zero, or center-of-mass, mode.

> **Note on units:** `freq` is interpreted as angular frequency $\omega$.
> The returned response functions are in the model's internal units and depend
> on the friction coefficient `zeta` used to define relaxation times.

### `compute_modulus(a, freq, zeta=1.0)`

Computes *system-level* moduli by summing contributions from all nonzero normal
modes.

- **Inputs**
  - `a`: `(N, N)` symmetric connectivity matrix
  - `freq`: `(n_freq,)` array of angular frequencies $\omega$
  - `zeta`: friction coefficient (default `1.0`)
- **Returns**
  - `G_storage`: `(n_freq, 2)` array with columns $[\omega, G'(\omega)]$
  - `G_loss`: `(n_freq, 2)` array with columns $[\omega, G''(\omega)]$

### `compute_monomer_mechanical_susceptibility(a, freq, zeta=1.0)`

Computes the real and imaginary parts of the *per-locus mechanical susceptibility*:

$$
\chi_i'(\omega) = \frac{1}{\zeta}\sum_{p>0}v_{pi}^2
\frac{\tau_p}{1+(\omega\tau_p)^2},
\qquad
\chi_i''(\omega) = \frac{1}{\zeta}\sum_{p>0}v_{pi}^2
\frac{\omega\tau_p^2}{1+(\omega\tau_p)^2}.
$$

Here, $\tau_p=-\zeta/\lambda_p$. There is no factor of two in these
monomer-level response functions.

- **Returns**
  - `freq`: `(n_freq,)`
  - `chi_prime_i`: `(n_freq, N)` array containing $\chi_i'(\omega)$
  - `chi_double_prime_i`: `(n_freq, N)` array containing $\chi_i''(\omega)$

#### Example: compute mechanical response

```python
import numpy as np
import hipps_dimes as HD

# Example: obtain a connectivity matrix 'a' from HIPPS-DIMES
# (you can also load a saved matrix from disk with np.loadtxt)
results = HD.run_optimization(
    input_matrix=np.loadtxt("contact_map.txt"),
    input_type="cmap",
    iteration=10000,
    verbose=False,
)
a = results["connectivity_matrix"]

# (1) Compute bulk moduli G'(ω), G''(ω)
freq = np.logspace(-3, 3, 200)  # angular frequencies ω
G_storage, G_loss = HD.compute_modulus(a, freq, zeta=1.0)

# (2) Compute per-locus mechanical susceptibilities
freq_out, chi_prime_i, chi_double_prime_i = (
    HD.compute_monomer_mechanical_susceptibility(a, freq, zeta=1.0)
)
```

## How to cite

If you use this program in a publication, please cite the following references:

* _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051._

* _Shi, G., Thirumalai, D. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nat Commun 14, 1150 (2023)._

* _Shi, G., Shin, S., and Thirumalai, D. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advances 11.31 (2025): eadx1763._

[^1]: _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051._
[^2]: _Shi, G., Thirumalai, D. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nat Commun 14, 1150 (2023)._
[^3]: _Shi, G., Shin, S., and Thirumalai, D. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advances 11.31 (2025): eadx1763._
