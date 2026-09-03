# Repository Guidelines

Use this file as the short operational guide for Codex and contributors. For the full module/test ownership map, detailed change routing, and repo quirks, read `MAINTAINER_MAP.md`.

## Project Structure

- `hipps_dimes/api.py`: workflow coordinator, input handling, result packaging, artifact writing
- `hipps_dimes/models.py`: optimization engine and `Dynamics` class
- `hipps_dimes/numerics.py`: scientific kernels and matrix math
- `hipps_dimes/covariance_pdhg.py`: noise-aware COV PDHG and hybrid solver
- `hipps_dimes/cli.py`: CLI parsing only; keep it thin
- `hipps_dimes/core.py`, `hipps_dimes/__init__.py`, `HippsDimes.py`: export and compatibility layers
- `tests/`: regression coverage split across `test_basic.py`, `test_dynamics.py`, `test_entropy_and_api.py`, and `test_covariance_pdhg.py`

## Working Rules

- Preserve both import surfaces: `import hipps_dimes` and legacy `import HippsDimes`.
- Prefer implementing behavior in `api.py`, `models.py`, or `numerics.py`, not in shims such as `core.py` or `HippsDimes.py`.
- Keep `cli.py` focused on argument parsing and forwarding into `run_optimization(...)`.
- When changing public behavior, update tests in the matching test file and update `README.md` if the user-facing contract changes.

## Test Commands

Local `pytest` is configured with coverage flags in `pyproject.toml`. If `pytest-cov` is unavailable, use:

```bash
pytest -q -o addopts=''
```

Focused runs:

```bash
pytest -q -o addopts='' tests/test_basic.py tests/test_dynamics.py
pytest -q -o addopts='' tests/test_entropy_and_api.py
```

## Change Routing

- CLI/options/output-contract changes: `api.py`, `cli.py`, `tests/test_entropy_and_api.py`
- COV optimization changes: `covariance_pdhg.py`, possibly `numerics.py`, plus `tests/test_covariance_pdhg.py`
- IS/GD/DI optimization changes: `models.py`, possibly `numerics.py`, plus API tests
- Dynamics/mechanics changes: `models.py`, `numerics.py`, `tests/test_dynamics.py`, and `tests/test_basic.py`
- Compatibility/export changes: `__init__.py`, `core.py`, `HippsDimes.py`, and import tests

## Repo Quirks

- The current `Makefile` only lints `HippsDimes.py`; CI additionally checks the COV solver and its dedicated tests.
