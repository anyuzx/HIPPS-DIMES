"""Backward-compatible module shim for HIPPS-DIMES.

This module preserves the historical import path (`import HippsDimes`) while
internally delegating implementation to the `hipps_dimes` package.
"""

from hipps_dimes import *  # noqa: F401,F403
from hipps_dimes.cli import main

if __name__ == "__main__":
    main()
