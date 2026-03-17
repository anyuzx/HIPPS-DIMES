"""Regression tests for the Dynamics class."""

import numpy as np
import pytest

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


def test_dynamics_run_twice_raises_with_resume_guidance():
    """A second fresh run should fail with an explicit resume/reset message."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(2, every=1, method="exact", update_zero_modes=False)

    with pytest.raises(RuntimeError, match="call resume\\(\\) to continue or reset\\(\\)"):
        model.run(1, every=1, method="exact", update_zero_modes=False)


def test_dynamics_resume_appends_from_current_state():
    """resume should continue from the current state and append to traj."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(3, every=1, method="exact", update_zero_modes=False)

    initial_traj = model.traj.copy()
    continuation_start = model.xyz.copy()

    model.resume(2, every=1, method="exact", update_zero_modes=False)

    assert model.traj.shape == (5, 4, 3)
    assert np.allclose(model.traj[:3], initial_traj)
    assert np.allclose(model.traj[3], continuation_start)
    assert np.allclose(model.traj_time, [0.0, 0.01, 0.02, 0.03, 0.04])


def test_dynamics_reset_allows_new_fresh_run():
    """reset should clear run history so run() can be called again."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(2, every=1, method="exact", update_zero_modes=False)

    model.reset()
    model.run(1, every=1, method="exact", update_zero_modes=False)

    assert model.traj.shape == (1, 4, 3)


def test_dynamics_run_excludes_final_state_by_default():
    """run should keep the historical default unless include_final_state is requested."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(5, every=2, method="exact", update_zero_modes=False)

    assert model.traj.shape == (3, 4, 3)
    assert np.allclose(model.traj_time, [0.0, 0.02, 0.04])


def test_dynamics_run_can_append_final_state_to_traj():
    """run should optionally append the post-integration state to traj."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(5, every=2, method="exact", update_zero_modes=False, include_final_state=True)

    assert model.traj.shape == (4, 4, 3)
    assert np.allclose(model.traj[-1], model.xyz)
    assert np.allclose(model.traj_time, [0.0, 0.02, 0.04, 0.05])
    assert model.time == pytest.approx(0.05)


def test_dynamics_reset_clears_traj_time():
    """reset should also clear saved trajectory times."""
    model = HippsDimes.Dynamics(4, k=1.0, model="rouse")
    model.initialize(dt=1e-2, zeta=1.0, beta=1.0)
    model.run(2, every=1, method="exact", update_zero_modes=False)

    model.reset()

    assert model.traj_time.shape == (0,)
    assert model.time == pytest.approx(0.0)
