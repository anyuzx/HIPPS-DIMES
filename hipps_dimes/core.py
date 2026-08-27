"""Compatibility aggregation module for HIPPS-DIMES internals."""

from . import numerics as _numerics
from .covariance_accuracy import (
    install_covariance_accuracy_guard as _install_covariance_accuracy_guard,
)

# Install the globalized covariance-cone solver before importing the public API.
# The original implementation remains available as
# ``numerics._legacy_fit_gaussian_noise_covariance`` for regression comparisons.
_install_covariance_accuracy_guard()

from .numerics import *  # noqa: F401,F403,E402
from .models import *  # noqa: F401,F403,E402
from .api import *  # noqa: F401,F403,E402
