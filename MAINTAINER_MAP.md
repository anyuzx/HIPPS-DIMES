# Maintainer Map

## Purpose

This repository implements the HIPPS-DIMES method for:

- reconstructing 3D chromatin structure ensembles from contact or distance constraints
- converting contact maps into distance constraints
- predicting chromatin dynamics and mechanics from the inferred connectivity matrix

The supported user-facing entry surfaces today are:

- console script: `HippsDimes` via `hipps_dimes.cli:main`
- module execution: `python -m hipps_dimes`
- library imports: `import hipps_dimes` and legacy `import HippsDimes`

## Top-Level Layout

- `hipps_dimes/api.py`: workflow coordinator; input validation/loading, missing-data handling, solver dispatch, result packaging, artifact writing
- `hipps_dimes/models.py`: optimization engine plus the `Dynamics` class and progress-callback plumbing
- `hipps_dimes/numerics.py`: numerical kernels, matrix transforms, ensemble sampling, dynamics/mechanics math, optional GPU helpers
- `hipps_dimes/core.py`: compatibility aggregation module re-exporting numerics, models, and API surfaces
- `hipps_dimes/cli.py`: Click/rich-click CLI, option parsing, and CLI-only normalization such as `save_steps`
- `hipps_dimes/__init__.py`: package public API re-export
- `hipps_dimes/__main__.py`: `python -m hipps_dimes` entrypoint
- `HippsDimes.py`: backward-compatible module and script shim
- `tests/test_basic.py`: import compatibility and smoke-level numerical regression tests
- `tests/test_dynamics.py`: `Dynamics` lifecycle, progress callbacks, resume/reset behavior, and trajectory persistence
- `tests/test_entropy_and_api.py`: `run_optimization(...)`, CLI behavior, file outputs, save-pickle/save-steps paths, and missing-data handling
- `.github/workflows/ci.yml`: test, lint, and build automation
- `pyproject.toml`: packaging metadata, dependency declarations, tool config, and pytest defaults
- `Makefile`: local convenience commands for install, test, lint, format, build, and publish
- `.pre-commit-config.yaml`: local formatting and lint hook configuration
- `data/`: example/reference inputs
- `doc/`: notes and figures referenced by the README
- `README.md`: user-facing install, CLI, and API documentation

## Runtime Flow

### CLI path

1. `HippsDimes ...` resolves to `hipps_dimes.cli:main`.
2. `cli.py` parses flags and positional arguments, normalizes CLI-only values, and calls `run_optimization(...)`.
3. `api.py` validates inputs, loads matrices, handles missing-data repair/removal, and constructs `Optimize(...)`.
4. `models.py` runs the selected solver (`IS`, `GD`, or `DI`).
5. `api.py` builds the result dictionary, writes requested artifacts, and optionally generates XYZ ensembles.

### Library path

1. Users import `hipps_dimes` or the compatibility shim `HippsDimes`.
2. Public functions are re-exported through `hipps_dimes.__init__`, which pulls from `core.py`.
3. `run_optimization(...)` is the main orchestration entrypoint.
4. Lower-level numerical helpers and the `Dynamics` class are also available via the package exports.

## Common Change Guide

### Add or change a CLI flag

Touch:

- `hipps_dimes/cli.py`
- `hipps_dimes/api.py`
- `README.md`
- `tests/test_entropy_and_api.py`

Typical pattern:

1. Add the Click option in `cli.py`.
2. Thread the new argument into `run_optimization(...)`.
3. Validate or apply the behavior in `api.py` or `models.py`.
4. Add at least one CLI or API regression test.

### Change input loading or missing-data handling

Touch:

- `hipps_dimes/api.py`
- `README.md` if the user-facing behavior changes
- `tests/test_entropy_and_api.py`

Current high-level cases include:

- text and `.npy`
- `cooler`
- `.hic`
- direct NumPy arrays passed as `input_matrix`
- optional missing-data repair/removal when `ignore_missing_data` is enabled

### Change optimization behavior

Touch:

- `hipps_dimes/models.py`
- `hipps_dimes/api.py` if the result contract, progress reporting, or orchestration changes
- `hipps_dimes/numerics.py` if the change affects kernels or objective calculations
- `tests/test_entropy_and_api.py`

Main solver surface:

- `Optimize.run(...)`
- `Optimize.run_noisy(...)`
- internal update methods in `Optimize`

### Change output files or returned results

Touch:

- `hipps_dimes/api.py`
- `README.md`
- `tests/test_entropy_and_api.py`

Current result/output surface includes:

- `iteration_series` and backward-compatible alias `log`
- `run_parameters`
- final distance map
- final contact map for `cmap` inputs
- connectivity matrix
- optional `connectivity_matrix_at_steps`
- optional XYZ ensemble
- optional pickle payload via `save_pickle`
- optional metadata about removed/retained loci for missing-data workflows

### Change dynamics or mechanics functionality

Touch:

- `hipps_dimes/models.py`
- `hipps_dimes/numerics.py`
- `tests/test_dynamics.py`
- `tests/test_basic.py` for smoke-level numerics coverage

Look at surfaces such as:

- `Dynamics.run(...)`
- `Dynamics.resume(...)`
- `Dynamics.run_with_force(...)`
- `Dynamics.run_breakable_bond(...)`
- `compute_acf_general_theory`
- `compute_m1_all`
- `compute_stress_relaxation`
- `compute_modulus`
- `compute_monomer_modulus`
- `a2dmap_theory_with_force_applied`

### Change compatibility or public import behavior

Touch:

- `hipps_dimes/__init__.py`
- `hipps_dimes/__main__.py`
- `hipps_dimes/core.py`
- `HippsDimes.py`
- `pyproject.toml`
- `tests/test_basic.py`

Keep in mind:

- some tests still import `HippsDimes`
- the console script name is still `HippsDimes`
- `core.py` is an export shim, not the primary implementation site

## Test Map

### `tests/test_basic.py`

Use this file when touching:

- package import compatibility
- base matrix constructors and conversions
- stress relaxation and modulus calculations
- smoke-level `Dynamics` helpers
- force-applied distance-map behavior
- small utility regressions

### `tests/test_dynamics.py`

Use this file when touching:

- `Dynamics.run(...)` and `resume(...)`
- reset semantics
- trajectory time bookkeeping
- progress callback payloads
- `show_progress=False` behavior
- forced and breakable-bond dynamics
- `save_traj(...)`

### `tests/test_entropy_and_api.py`

Use this file when touching:

- `run_optimization(...)`
- CLI-to-API wiring
- output file writes
- log/run-parameter DataFrame contract
- `save_steps` and `save_pickle`
- progress callback behavior during optimization
- `cmap`/`dmap`/`ddmap` inputs
- missing-data repair/removal and related metadata
- entropy calculation

## Practical Commands

Run the full test suite without requiring `pytest-cov` locally:

```bash
pytest -q -o addopts=''
```

Run the most relevant focused test files:

```bash
pytest -q -o addopts='' tests/test_basic.py tests/test_dynamics.py
pytest -q -o addopts='' tests/test_entropy_and_api.py
```

Install the package in editable mode:

```bash
pip install -e .
```

Install development dependencies:

```bash
pip install -e ".[dev]"
```

## Known Repo Quirks

- `pyproject.toml` hardcodes coverage-related `pytest` addopts. If `pytest-cov` is not installed, plain `pytest` fails unless you override `addopts`.
- The current `Makefile` and CI lint job only target `HippsDimes.py`, not the full `hipps_dimes/` package. Package modules can therefore drift from lint coverage unless those commands are expanded.
- The historical `HippsDimes` import path is still part of the supported surface and is covered by tests.
- `README.md` documents the main CLI/API workflow, but newer programmatic features such as `save_pickle`, `save_steps`, missing-data locus removal, and progress callbacks are easier to verify from `api.py` and `tests/test_entropy_and_api.py`.

## Safe Mental Model

- `api.py` is the workflow and output-contract owner.
- `models.py` owns iterative optimization state and simulation lifecycle behavior.
- `numerics.py` owns scientific kernels and matrix math.
- `cli.py` should stay thin and mostly parse/forward.
- `core.py`, `__init__.py`, and `HippsDimes.py` are compatibility/export layers and should contain very little business logic.
