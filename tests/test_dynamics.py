"""Regression tests for the Dynamics class."""

import numpy as np

import HippsDimes


def test_dynamics_copies_input_connectivity_matrix():
    """External mutations of the input matrix should not affect the model."""
    a = np.array(
        [
            [2.0, -1.0, -1.0],
            [-1.0, 2.0, -1.0],
            [-1.0, -1.0, 2.0],
        ]
    )
    a_initial = a.copy()

    model = HippsDimes.Dynamics(a)

    a[0, 1] = 123.0
    a[1, 0] = 123.0

    assert not np.shares_memory(model.A, a)
    assert np.allclose(model.A, a_initial)
    assert np.allclose(model.eigvalue, np.linalg.eigh(a_initial)[0])
