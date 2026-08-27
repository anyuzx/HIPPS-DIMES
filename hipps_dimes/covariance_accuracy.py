"""Accuracy guard for the public globalized COV entry point.

The covariance-gradient norm and the pairwise stationary residual are equivalent
at the exact optimum, but their finite numerical tolerances have different
scales.  The public solver therefore uses a modestly tighter internal first-order
tolerance than the user-facing tolerance.  This preserves the requested contract
for pairwise stationarity without changing the objective or returned solution.
"""

from __future__ import annotations

from functools import wraps

import numpy as np

from . import covariance_solver as _solver
from . import numerics as _num


_UNWRAPPED_SOLVER = _solver.fit_gaussian_noise_covariance
_INTERNAL_TOLERANCE_FACTOR = 0.05
_MINIMUM_INTERNAL_RELATIVE_TOLERANCE = 1e-13
_MINIMUM_INTERNAL_ABSOLUTE_TOLERANCE = 1e-14


@wraps(_UNWRAPPED_SOLVER)
def fit_gaussian_noise_covariance(*args, **kwargs):
    """Run COV with enough internal accuracy to satisfy pair stationarity."""
    requested_relative = float(kwargs.get("relative_tolerance", 1e-8))
    requested_absolute = float(kwargs.get("absolute_tolerance", 1e-10))

    if requested_relative > 0.0:
        internal_relative = max(
            _MINIMUM_INTERNAL_RELATIVE_TOLERANCE,
            _INTERNAL_TOLERANCE_FACTOR * requested_relative,
        )
    else:
        internal_relative = 0.0
    if requested_absolute > 0.0:
        internal_absolute = max(
            _MINIMUM_INTERNAL_ABSOLUTE_TOLERANCE,
            _INTERNAL_TOLERANCE_FACTOR * requested_absolute,
        )
    else:
        internal_absolute = 0.0

    kwargs["relative_tolerance"] = internal_relative
    kwargs["absolute_tolerance"] = internal_absolute
    fitted, gram, connectivity, info = _UNWRAPPED_SOLVER(*args, **kwargs)
    info["requested_relative_tolerance"] = requested_relative
    info["requested_absolute_tolerance"] = requested_absolute
    info["internal_relative_tolerance"] = internal_relative
    info["internal_absolute_tolerance"] = internal_absolute
    info["stationarity_accuracy_guard_factor"] = (
        _INTERNAL_TOLERANCE_FACTOR
    )
    return fitted, gram, connectivity, info


def install_covariance_accuracy_guard() -> None:
    """Install the guarded solver in both the module and public aggregator."""
    if not hasattr(_num, "_legacy_fit_gaussian_noise_covariance"):
        _num._legacy_fit_gaussian_noise_covariance = (
            _num.fit_gaussian_noise_covariance
        )
    _solver.fit_gaussian_noise_covariance = fit_gaussian_noise_covariance
    _num.fit_gaussian_noise_covariance = fit_gaussian_noise_covariance
