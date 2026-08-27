![movie](doc/source/movie.gif)

[**Installation**](#install) | [**Quick Start**](#quick-start) | [**API**](#api) | [**Chromatin Dynamics**](#dynamics-prediction-functionality) | [**Chromatin Mechanics**](#modulus-calculation)

This python program is the implementation of the **HIPPS-DIMES** method[^1][^2][^3]. HIPPS-DIMES is a computational method based on the maximum entropy principle, with experimental measured contact map or pair-wise distances as constraints, to generate a unique ensemble of <ins>3D chromatin structures</ins>. In a nutshell, this program accepts the input file of a mean spatial distance map (which can be measured in Multiplexed FISH experiment) or a Hi-C contact map (which is converted to distance map internally), and generates an ensemble of individual chromatin conformations that are consistent with the input.

In addition to reconstructing static 3D chromatin structures, the model is also able to predict **chromatin loci dynamics** as well as **chromatin loci mechanics** based on polymer physics and the Ornstein–Uhlenbeck process. This allows you to simulate time-dependent properties such as autocorrelation functions (ACF), mean-square displacement (MSD) of individual loci and modulus of individual loci, providing insights into the dynamic behavior of chromatin.

![schematic](doc/source/schematic.jpg)

The theory and applications of this method can be found in our work published:

- Shi, Guang, and D. Thirumalai. "Epigenetic state encodes locus-specific chromatin mechanics." bioRxiv (2025): 2025-12. [link](https://www.biorxiv.org/content/10.64898/2025.12.27.696709v1.abstract)
- Shi, Guang, and D. Thirumalai. "From Hi-C contact map to three-dimensional organization of interphase human chromosomes." Physical Review X 11.1 (2021): 011051. [link](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051)
- Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants." Nature Communications 14.1 (2023): 1150. [link](https://www.nature.com/articles/s41467-023-36412-4)
- Shi, Guang, Sucheol Shin, and D. Thirumalai. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advances 11.31 (2025): eadx1763. [link](https://www.science.org/doi/full/10.1126/sciadv.adx1763)

Other applications of this method can be found in:
- Dey, Atreya, et al. "Structural changes in chromosomes driven by multiple condensin motors during mitosis." Cell Reports 42.4 (2023).
- Jeong, Davin, et al. "Structural basis for the preservation of a subset of topologically associating domains in interphase chromosomes upon cohesin depletion." Elife 12 (2024): RP88564.

# Documentation

## Requirements

- Python 3.9+

## Install

First, download this repository using,

```bash
git clone https://github.com/anyuzx/HIPPS-DIMES
cd HIPPS-DIMES
```

#### Using pip
```bash
pip install -e .
```

#### Using UV
```bash
uv pip install -e .
```

This command will install the required packages from `pyproject.toml`, and install the script as a
python module. Once installed, you can call `HippsDimes` directly in the terminal to run the script.

### Dependencies

The package requires:
- Python 3.9+
- `Click` - Command line interface (fallback when rich-click is not installed)
- `Rich-Click` - Rich-styled CLI (used when installed; falls back to Click otherwise)
- `Numpy` - Numerical computing
- `Scipy` - Scientific computing
- `Pandas` - Data manipulation
- `Tqdm` - Progress bars
- `Cooler` - .mcool or .cool data format support
- `Rich` - Rich terminal output
- `hic-straw` - .hic format support

**Optional (for GPU acceleration):**
- `CuPy` - GPU-accelerated computing via CUDA

To install CuPy for GPU support:
```bash
conda install -c conda-forge cupy
# or
pip install cupy-cuda11x  # Replace with your CUDA version
```

## Quick start

The quickest way to get started is to run our example notebook in Google Colab. Here are some starter notebooks:

- [Quick start](https://colab.research.google.com/drive/1w7cK6S3z2_D5Mzgq2-SZgl4FoDjT9EC5?usp=sharing)
- [Chromatin Dynamics and Mechanics](https://colab.research.google.com/drive/1DsTIJTkiKc1vRq6lBncEg4nvPh4hQ4hp?usp=sharing)

In addition, it is helpful to view the help information for each arguments andoptions. To display help information, use

```bash
HippsDimes --help
```

## API

### Input files

This script accepts multiple input formats. If the input file is a Hi-C
contact map, it can be in `.cool` format (see https://github.com/open2c/cooler
for details of the `cooler` library), `.hic` format (via hicstraw), plain text
matrix format, or NumPy `.npy` format. If the input file is a mean spatial
distance map, the script accepts plain text matrix format or NumPy `.npy`
format. The text format for a matrix is the following: each row of the file
corresponds to the row of the matrix. Values are space-separated. The content
of the file should look like this,

```text
1  2  3
2  1  2
3  2  1
```

### Output files

This script will generate several files:

- A text file for the final simulated mean distance map
- A text file for the final simulated contact map (this is the best agreement contact map to the _normalized_ input contact map)
- For `cooler` and `.hic` contact-map inputs, a text file for the internal target contact map used to construct the target distance constraints: `{output_prefix}_cmap_target.txt`
- A text file for the connectivity matrix
- A `.xyz` formatted file for the ensemble of genome structures generated (can
  be turned off)
- A CSV file for run parameters: `{output_prefix}_run_parameters.csv` (can be turned off)
- A CSV file for iteration-series scalar data: `{output_prefix}_iteration_series.csv` (can be turned off)

### Explanation of the arguments and options

#### Arguments

- `INPUT`: File path for the input file. The input file can be a Hi-C contact
  map or a mean spatial distance map as measured in Multiplexed FISH experiment.
- `OUTPUT_PREFIX`: Prefix for output files. For instance, if one specifies it to be
  `TEST`, then all the output files will start with `TEST_`.

#### Options

- `-k, --connectivity-matrix`: Provide the path to the existing connectivity
  matrix one would like to use as initialization. Useful if restarting using the
  result from the previous run.
- `-e, --ensemble`: Number of individual conformations to be generated. This
  script will generate an ensemble of structures consistent with the input Hi-C
  contact map or the mean spatial distance map. Each individual conformations
  are different from each other. You can specify how many such individual
  conformations you want to generate. Default: 1000.
- `-a, --alpha`: Value of the contact map to distance map conversion
  exponent. If the input file is Hi-C contact map, the method first converts the
  contact map to a mean spatial distance map. The equation of the conversion is
  d_{ij} ~ c_{ij}^{1/\alpha}. The default value of α is 4.0, estimated in
  this work 10.1126/science.aaf8084. Default: 4.0.
- `-s, --selection`: Specify chromosome or region. This option is required
  when the input file has `cooler` or `.hic` format. For cooler files, the value is passed to the `cooler.Cooler.matrix().fetch()` method. For .hic files, use format "chr1:start1-end1,chr2:start2-end2". For details on cooler selectors, please refer to their [documentation](https://cooler.readthedocs.io/en/latest/concepts.html#matrix-selector).
- `-m, --method`: Select IS (Iterative Scaling, default), GD (Gradient Descent), DI (Direct Inversion), or COV (calibrated Gaussian covariance-cone optimization).
- `-l, --lamd`: Specify the weight for L1 or L2 regularization. Default value is 0.0, meaning no regularization. Regularization is typically used to avoid over-fitting.
- `-r, --reg`: Specify the type of regularization. Options: L1, L2 (default). This option should be used together with option `-l`.
- `--gaussian-noise-variance`: Positive scalar absolute variance on squared-distance constraints. COV only.
- `--gaussian-noise-relative-std`: Positive scalar relative standard deviation `sigma_ij / Dobs_ij`. COV converts it to pair variance `(value * Dobs_ij)^2` after preprocessing. Mutually exclusive with `--gaussian-noise-variance`.
- `--covariance-optimizer`: COV optimizer: `hybrid` (default), `pdhg`, or `newton`.
- `--covariance-relative-tolerance`: Relative COV KKT tolerance. Default: `1e-5`.
- `--covariance-absolute-tolerance`: Absolute internal COV KKT tolerance. Default: `1e-10`.
- `--covariance-handoff-relative-tolerance`: Relative KKT threshold for the
  default PDHG-to-Newton handoff. Default: `1e-3`.
- `-i, --iteration`: Maximum optimizer iterations. Default: 10000.
- `--learning-rate`: Learning rate. This hyperparameter controls the speed of convergence. If its value is too small, then convergence is very slow. If its value is too large, the program may never converge. Typically, learning rate can be set to be 1-30 if using Iterative scaling method. It should be a very small value (such as 1e-8) when using gradient descent optimization. Default: 10.0.
- `--momentum`: Momentum coefficient for IS method (0.0 to 1.0). Accelerates convergence by accumulating gradient history. **Recommended: Use 0.95 with `--nesterov` for fastest convergence (~50% faster).** Use 0.9 for more conservative settings. Only applies when method=IS. Default: 0.0.
- `--nesterov`: Use Nesterov Accelerated Gradient (NAG). Enables higher momentum values (0.95) without divergence. **Recommended: Use with `--momentum 0.95` for best performance.**
- `--use-gpu`: Enable GPU acceleration via CuPy. All COV optimizers use
  float64. COV fails rather than silently falling back when CUDA is
  unavailable.
- `--gpu-float32`: Use float32 for legacy GPU IS/GD. COV is float64-only.
- `--save-steps`: Comma-separated list of iteration steps at which to save the connectivity matrix. Example: `--save-steps 1000,5000,10000`. Files are saved as `{output_prefix}_connectivity_matrix_iter{step}.txt`. When used as a library (without `output_prefix`), connectivity matrices at these steps are still returned in `results['connectivity_matrix_at_steps']`.
- `--eigh-threads`: Number of threads for eigenvalue (eigh) and BLAS/LAPACK. If not set, the backend default is used. Set to 1 for single-threaded runs.
- `--input-type`: The type of the input file. To use the script, the type must be specified. Options: `cmap` (contact map) or `dmap` (distance map). This option is required.
- `--input-format`: The format of the input file. Options: `text`, `npy`, `cooler`, or `hic`. If the type of input file is Hi-C contact map, then the script supports `cooler` format Hi-C contact map file, `.hic` format, pure text-based file, or NumPy `.npy` file. In the text-based file, each line corresponds to the row of the contact map. If the type of input file is mean distance map, the script supports text-based matrix files and `.npy` files. This option is required.
- `--binsize`: Bin size (resolution) for .hic format in bp. Default: 25000.
- `--norm`: Normalization for .hic format. Options: KR, VC, NONE. Default: KR.
- `--unit`: Unit for .hic format. Options: BP, FRAG. Default: BP.
- `--no-log`: By default, the program writes two log files when `output-prefix` is provided: `{output_prefix}_run_parameters.csv` and `{output_prefix}_iteration_series.csv`. Use `--no-log` to disable writing both files.
- `--no-xyzs`: Turn off writing x,y,z coordinates of genome structures to files.
- `--ignore-missing-data`: Turn on this argument will let the program ignore the missing elements or infinite numbers in the contact map or distance map.
- `--balance`: Turn on the matrix balance for contact map. Only effective when `input_type == cmap` and `input_format == cooler`.
- `--neighbor-balance`: Turn on neighbor balancing for contact map. Only effective when `input_type == cmap`. Normalizes contact between i and j by dividing it by the geometric mean of neighbor contact for i and j. See Paggi, Zhang 2025 for method details.
- `--not-normalize`: Turn off the auto normalization of the contact map. Only effective when `input_type == cmap`.
- `--enforce-nonnegative-connectivity-matrix`: Constrain all the "spring constants" to be nonnegative.
- `-q, --quiet`: Quiet mode: disable fancy tables display, keep only the progress bar.

#### COV quick start

COV is the noise-aware optimization method. Specify exactly one noise model:

- `--gaussian-noise-variance <v>` for homoskedastic absolute variance
  (`v_ij = v`).
- `--gaussian-noise-relative-std <w>` for heteroskedastic relative noise
  (`sigma_ij = w D_ij`).

```bash
python -m hipps_dimes observed_ddmap.npy cov_fit \
  --input-type ddmap --input-format npy \
  --method COV --gaussian-noise-relative-std 0.1 \
  --use-gpu --iteration 10000 --no-xyzs --save-pickle
```

The default hybrid optimizer starts from a calibrated Rouse model, runs PDHG,
and hands off to Newton refinement. It stops automatically when the built-in
KKT convergence criterion is met or the requested iteration limit is reached.

See [Noise-aware covariance optimization](doc/NOTE_ON_NOISE.md) for the noise
model, initialization, and scientific interpretation. See
[PDHG and hybrid COV solver](doc/COV_PDHG.md) for solver details, convergence
diagnostics, and advanced tuning.

### Examples

#### Example 1

First, download a cooler format Hi-C contact map from
[here](https://drive.google.com/file/d/1eIxGv1JbIrEAVoUSQK_O_ebIjWo6toTJ/view?usp=sharing)
(**The file size is about 116 Mb**). This Hi-C contact map is for Chicken cell
mitotic chromosome, originally retrieved from
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740).
Rename it to `hic_example.cool`. Then execute the following command,

```bash
HippsDimes hic_example.cool test --input-type cmap --input-format cooler -s chr7:10M-15M -i 10 -e 10
```

This command tells the script to load the Hi-C contact map `hic_example.cool`
and perform the iterative scaling algorithm. The argument `test` instructs the
file names of output files start with `test_`. Option `--input-type cmap`
specifies that the input file is a contact map. Option `--input-format cooler`
specifies that the input file is a `cooler` file. Option `-s chr7:10M-15M`
specifies that the algorithm is performed on the region 10 Mbps - 15 Mbps on
Chromosome 7. Note that these three options are required and cannot be
neglected. **Some option arguments are optional, some are required. Please refer
to the section below and use `HippsDimes --help` for details**

When the program finishes, the script will generate several output files:
`test.xyz`, `test_connectivity_matrix.txt`, and `test_dmap_final.txt`.
`test.xyz` contains 10 sets of individual conformations of x, y, z coordinates
and can be viewed using `VMD` or other compatible visualization softwares.

#### Example 2

In this example, we use Hi-C contact map for HeLa cell line Chromosome 14 at
time point of 12 hours after the release from prometaphase. For the purpose of
demonstration, you can download the Hi-C `.cool` file from
[here](https://drive.google.com/file/d/1j-zfDUP6LOZGCxz9uA3LaMI372ct1cU_/view?usp=sharing)
which is originally retrieved from
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740)
under accession number GSE102740. **Before you download, note that the file has
size of about 655 Mb**. Once downloaded, execute the following command,

```bash
HippsDimes GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool::6 test --input-type cmap --input-format cooler -s chr14:20M-107M -i 10000 -e 10
```

Similar to the first example, this command tells the script to load the Hi-C
cooler file `GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool` and its group
6 (data for several different resolution is stored in different groups) and
perform the HIPPS/DIMES algorithm. In this example, we change the number of
iterations to be 10000 by using the option `-i 10000`. On a AMD Ryzen 5 3600 CPU
machine, it takes about 3-4 mins to finish the program. Once it is finished,
several output files are generated.

#### Example 3

A `.hic` format example:

```bash
HippsDimes mydata.hic test \
  --input-type cmap --input-format hic \
  --selection chr1:31000000-41000000,chr1:31000000-41000000 \
  --binsize 25000 --norm KR --unit BP \
  -i 10000 -e 10
```

> In particular, if you would like see an example of direct application of HIPPS-DIMES on imaging data, please go through the notebook.

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
- Use `--ignore-missing-data` to exclude missing or zero input pairs from the
  constraints. If a locus has no observed off-diagonal pairs at all, combine it
  with `--remove-fully-missing-loci` to remove that locus before optimization.
- Contact-map inputs are normalized by their maximum entry by default. Use
  `--not-normalize` when the supplied contact map should retain its existing
  scale.
- Note that when feeding the contact map, there is no physical length scale
  associated with it. Thus we cannot set a unit to the resulting distance matrix
  or the structures. In this sense, the structures generated are dimensionless.
  But one can use additional information to set the length scale of the problem.
  For instance, if you have a reasonable estimate of the average distance between
  the two nearest loci, then you can use this distance as the measure to rescale the
  structure to be consistent with it.

### How to choose optimization method, number of iteration and convergence

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

COV has a well-defined objective and built-in convergence tests. The default
hybrid PDHG-to-Newton solver stops automatically only after the independently
recomputed KKT residual meets the requested tolerance. For COV, `--iteration`
therefore sets a maximum update budget rather than a recommended stopping
iteration. If the budget is exhausted first, the returned convergence status
is false and the user should not describe the result as converged. Inspect
`results['covariance_optimization']['converged']`, `status`, and
`relative_eliminated_kkt_residual` for the final certificate.

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

### Use it as a Python Library

The HIPPS-DIMES code can be used both as a command-line tool and as a Python library. This makes it easy to integrate into your Python workflows, Jupyter notebooks, and automated pipelines.

#### Main Function: `run_optimization()`

The core functionality is available through the `run_optimization()` function:

```python
from HippsDimes import run_optimization
import numpy as np

# Load your contact map
cmap = np.loadtxt('contact_map.txt')

# Run optimization programmatically
results = run_optimization(
    input_matrix=cmap,          # Provide matrix directly
    input_type='cmap',
    method='IS',
    iteration=10000,
    learning_rate=10.0,
    momentum=0.95,              # Use momentum for faster convergence
    nesterov=True,              # Enable Nesterov acceleration
    use_gpu=True,               # Use GPU if CuPy is available
    ensemble=1000,
    verbose=False               # Suppress console output
)

# Access results
connectivity_matrix = results['connectivity_matrix']
structures = results['xyzs']          # (ensemble, n_beads, 3)
final_dmap = results['dmap_final']
final_cmap = results['cmap_final']
iteration_series = results['iteration_series']
run_parameters = results['run_parameters']
```

#### Return Values

The `run_optimization()` function returns a dictionary with:
- `'connectivity_matrix'`: Final connectivity matrix (numpy array)
- `'dmap_final'`: Final distance map (numpy array)
- `'cmap_final'`: Final contact map (numpy array, if input_type='cmap')
- `'xyzs'`: Generated conformations (numpy array)
- `'iteration_series'`: Iteration-series scalar outputs (pandas DataFrame; currently columns: iteration, loss, entropy)
- `'run_parameters'`: Run parameters (pandas DataFrame with columns: parameter, value)
- `'log'`: Alias for `'iteration_series'` (backward compatibility)
- `'rc_optimal'`: Optimal contact threshold (float, if input_type='cmap')

#### Additional Utility Functions

The package also provides helper functions for direct use:

```python
import HippsDimes as HD

# Generate structures from connectivity matrix
structures = HD.a2xyz_sample(connectivity_matrix, ensemble=1000)

# Convert connectivity matrix to distance map
dmap = HD.a2dmap_theory(connectivity_matrix)

# Convert connectivity matrix to contact map
cmap = HD.a2cmap_theory(connectivity_matrix, rc=5.0)

# Create a Rouse chain connectivity matrix
A = HD.construct_connectivity_matrix_rouse(n=100, k=1.0)
```

------

## Dynamics Prediction Functionality
In addition to reconstructing static 3D chromatin structures from contact or distance maps, the HIPPS-DIMES code now includes modules to simulate the dynamics of the chromatin[^3]. This new functionality is based on polymer physics and the Ornstein–Uhlenbeck process, which allows you to investigate time-dependent properties such as the autocorrelation function (ACF) and mean-square displacement (MSD) of individual locus.

### New Functions
- **`compute_acf_general_theory(i, j, t, a, zeta=1.0)`**: This function numerically computes the time-dependent autocorrelation function between monomers *i* and *j* using the connectivity matrix `a`. In addition to the ACF, it returns the corresponding two-point MSD for each time point in the provided time array `t`.

- **`compute_m1_i(i, t, a, zeta=1.0)`**: This function computes the single-loci mean-square displacement (MSD) for the i-th monomer as a function of time. The output is a 2D array where the first column is time and the second is the MSD of the monomer.

### New `Dynamics` Class
The new `Dynamics` class encapsulates the routines needed to run dynamic simulations provided the connectivity matrix `a`.

#### Example code
```python
from HippsDimes import Dynamics

model = Dynamics(a) # a is the connectivity matrix
model.initialize(dt=1e-2, zeta=1.0, beta=1.0)

model.run(int(1e5), every = 10)
model.resume(int(5e4), every = 10)
```

Trajectory coordinates can be accessed through `model.traj`. It is a `TxNx3` numpy array. `T` is the number of snapshots. `N` is the number of loci and 3 corresponds to coordinates at x, y, z dimensions. The corresponding reduced simulation time for each saved snapshot is stored in `model.traj_time`, a length-`T` numpy array.

You can save both arrays together with `model.save_traj("traj.npz")`. The `.npz` file stores two arrays: `traj` and `traj_time`.

`Dynamics.run(...)` starts a fresh trajectory and can be called only once per simulation state. To continue an existing simulation, use `Dynamics.resume(...)`. When omitted, `resume(...)` reuses the previous passive simulation settings for `update`, `every`, `method`, and `update_zero_modes`; you can still override any of them explicitly. To discard the previous trajectory and start over on the same object, call `Dynamics.reset()` before `run(...)`. By default, `run(...)` does not append the post-integration final state to `model.traj`; set `include_final_state=True` if you want that last state included.

### Dynamics under external force: `Dynamics.run_with_force(...)`

In addition to passive dynamics (`Dynamics.run`), you can simulate trajectories with a constant external force applied to selected loci.

- **Key parameters**
  - `force_loci`: list of locus indices where the force is applied
  - `force_amplitude`: force magnitude
  - `force_direction`: `(3,)` direction vector (it is normalized internally)
  - `force_duration`: optional number of timesteps to apply the force (if `None`, force is applied for the whole run)

#### Example: forced dynamics

```python
import numpy as np
from HippsDimes import Dynamics

# 'a' is the connectivity matrix (e.g., load from disk or obtain from run_optimization)
a = np.loadtxt("my_connectivity_matrix.txt")
model = Dynamics(a)  # a is the connectivity matrix
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

HIPPS-DIMES provides utilities to compute system-level **linear viscoelastic moduli** and per-locus **mechanical susceptibilities** from a connectivity matrix `a`. These routines decompose the polymer into normal modes, excluding the zero/center-of-mass mode.

> **Note on units**: `freq` is interpreted as **angular frequency** \(\omega\). The returned response functions are in the model’s internal units and depend on the friction coefficient `zeta` used to define relaxation times.

### `compute_modulus(a, freq, zeta=1.0)`

Computes *system-level* moduli by summing contributions from all non-zero normal modes.

- **Inputs**
  - `a`: `(N, N)` symmetric connectivity matrix
  - `freq`: `(n_freq,)` array of angular frequencies \(\omega\)
  - `zeta`: friction coefficient (default `1.0`)
- **Returns**
  - `(freq, G_storage)` as a `(n_freq, 2)` array (`[omega, G'(omega)]`)
  - `(freq, G_loss)` as a `(n_freq, 2)` array (`[omega, G''(omega)]`)

### `compute_monomer_mechanical_susceptibility(a, freq, zeta=1.0)`

Computes the real and imaginary parts of the *per-locus mechanical susceptibility*:

\[
\chi_i'(\omega) = \frac{1}{\zeta}\sum_{p>0}v_{pi}^2
\frac{\tau_p}{1+(\omega\tau_p)^2},
\qquad
\chi_i''(\omega) = \frac{1}{\zeta}\sum_{p>0}v_{pi}^2
\frac{\omega\tau_p^2}{1+(\omega\tau_p)^2}.
\]

Here, \(\tau_p=-\zeta/\lambda_p\). There is no factor of two in these monomer-level response functions.

- **Returns**
  - `freq`: `(n_freq,)`
  - `chi_prime_i`: `(n_freq, N)` array containing \(\chi_i'(\omega)\)
  - `chi_double_prime_i`: `(n_freq, N)` array containing \(\chi_i''(\omega)\)

#### Example: compute mechanical response

```python
import numpy as np
from HippsDimes import (
    run_optimization,
    compute_modulus,
    compute_monomer_mechanical_susceptibility,
)

# Example: obtain a connectivity matrix 'a' from HIPPS-DIMES
# (you can also load a saved matrix from disk with np.loadtxt)
results = run_optimization(input_matrix=np.loadtxt("contact_map.txt"), input_type="cmap", iteration=10000, verbose=False)
a = results["connectivity_matrix"]

# (1) Compute bulk moduli G'(ω), G''(ω)
freq = np.logspace(-3, 3, 200)  # angular frequencies ω
G_storage, G_loss = compute_modulus(a, freq, zeta=1.0)

# (2) Compute per-locus mechanical susceptibilities
freq_out, chi_prime_i, chi_double_prime_i = (
    compute_monomer_mechanical_susceptibility(a, freq, zeta=1.0)
)
```

## How to cite

If you used this program in your publication, please cite from the following
reference:

* _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051._

* _Shi, G., Thirumalai, D. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nat Commun 14, 1150 (2023)._

* _Shi, G., Shin, S., and Thirumalai, D. "Static three-dimensional structures determine fast dynamics between distal loci pairs in interphase chromosomes." Science Advances 11.31 (2025): eadx1763._

[^1]: _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional Organization of Interphase Human Chromosomes." Physical Review X 11.1 (2021): 011051._
[^2]: _Shi, G., Thirumalai, D. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nat Commun 14, 1150 (2023)._
[^3]: _Shi, G., Shin, S., and Thirumalai, D. Static Three-Dimensional Structures Determine Fast Dynamics Between Distal Loci Pairs in Interphase Chromosomes. bioRxiv (2025)._
