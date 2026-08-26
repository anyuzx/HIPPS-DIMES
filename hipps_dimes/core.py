"""Compatibility aggregation module for HIPPS-DIMES internals."""

from . import numerics as _numerics
from .covariance_globalization import install_covariance_globalization as _install_covariance_globalization

# Install the globalized covariance-cone solver before importing the public API.
# This preserves the existing API surface while keeping the legacy implementation
# available as ``numerics._legacy_fit_gaussian_noise_covariance`` for regression
# tests and controlled comparisons.
_install_covariance_globalization()

from .numerics import *  # noqa: F401,F403,E402
from .models import *  # noqa: F401,F403,E402
from .api import *  # noqa: F401,F403,E402
