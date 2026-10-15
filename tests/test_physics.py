"""Physics regression tests for the universal renormalization engine.

Run with:  python -m pytest tests/ -q
or standalone: python tests/test_physics.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from dvdp import renormalization as rg          # noqa: E402
from dvdp import integrators as it              # noqa: E402
from dvdp import metrics as mt                  # noqa: E402


def test_duffing_frequency_shift():
    """For Duffing at omega=1, the per-unit-eps frequency shift must be 3/8.

    This reproduces the Bender-Orszag / multiple-scales result quoted in the
    AFOSR proposal (alpha_omega = 3/8 at A=1, omega=1).
    """
    A_grid = np.array([1.0])
    _, beta_theta = rg.extract_generators("duffing", A_grid, omega=1.0)
    assert abs(beta_theta[0] - 0.375) < 5e-3, beta_theta[0]


def test_duffing_amplitude_preserved():
    """Duffing should renormalize frequency, not amplitude: beta_A ~ 0."""
    A_grid = np.array([1.0])
    beta_A, _ = rg.extract_generators("duffing", A_grid, omega=1.0)
    assert abs(beta_A[0]) < 1e-2, beta_A[0]


def test_vdp_amplitude_growth():
    """Van der Pol should grow the amplitude toward the limit cycle (A->2)."""
    A_grid = np.array([1.0])
    beta_A, _ = rg.extract_generators("vdp", A_grid, omega=1.0)
    # Weak-coupling theory: dA/dt = (alpha/2) A (1 - A^2/4); at A=1 this is +3/8.
    assert abs(beta_A[0] - 0.375) < 5e-3, beta_A[0]


def test_renorm_beats_naive_duffing():
    """Renormalized Duffing must track the full solution far better than naive."""
    t = np.linspace(0, 100, 4000)
    eps = 0.1
    phi_full, _ = it.solve_full(t, omega=1.0, eps=eps, A0=1.0)
    phi_naive, _ = it.solve_naive(t, omega=1.0, eps=eps, A0=1.0)
    res = rg.renormalize(t, omega=1.0, eps=eps, A0=1.0)
    e_rg = mt.rms_error(res.phi_R, phi_full)
    e_naive = mt.rms_error(phi_naive, phi_full)
    assert e_rg < e_naive, (e_rg, e_naive)
    assert e_rg < 0.1, e_rg


def test_combined_limits_reduce():
    """alpha=0 -> pure Duffing generators; eps=0 -> pure VdP generators."""
    t = np.linspace(0, 50, 2000)
    res_d = rg.renormalize(t, omega=1.0, eps=0.1, alpha=0.0, A0=1.0)
    res_v = rg.renormalize(t, omega=1.0, eps=0.0, alpha=0.1, A0=1.0)
    assert "vdp" not in res_d.components
    assert "duffing" not in res_v.components


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS  {name}")
            except AssertionError as e:
                print(f"FAIL  {name}: {e}")
