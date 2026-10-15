"""Time integrators: a correct fixed-step RK4 and scipy solve_ivp wrappers.

Used both to produce the reference "full" solution and for the numerical-method
benchmark (Phase 5).  The hand-written RK4 here fixes the step-size bug present
in the legacy notebooks, where ``h`` was derived from a ``t_span`` inconsistent
with ``t_eval``.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .systems import full_rhs
from . import background as bg
from . import perturbation as pert


def rk4(rhs, t_eval, y0, args=()):
    """Fixed-step classical RK4 on a uniform grid ``t_eval``.

    The step is taken directly between consecutive grid points, so the grid
    and the integration step are always consistent.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    y = np.asarray(y0, dtype=float)
    out = np.empty((len(y0), len(t_eval)))
    out[:, 0] = y
    for k in range(1, len(t_eval)):
        t = t_eval[k - 1]
        h = t_eval[k] - t_eval[k - 1]
        k1 = np.asarray(rhs(t, y, *args))
        k2 = np.asarray(rhs(t + h / 2, y + h / 2 * k1, *args))
        k3 = np.asarray(rhs(t + h / 2, y + h / 2 * k2, *args))
        k4 = np.asarray(rhs(t + h, y + h * k3, *args))
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[:, k] = y
    return out


def solve_full(t_eval, omega=1.0, eps=0.0, alpha=0.0, A0=1.0, theta0=0.0,
               method="DOP853", rtol=1e-10, atol=1e-12):
    """Reference ("full", exact) solution of the combined system."""
    y0 = [A0 * np.cos(theta0), -A0 * omega * np.sin(theta0)]
    sol = solve_ivp(full_rhs, (t_eval[0], t_eval[-1]), y0,
                    args=(omega, eps, alpha), t_eval=t_eval,
                    method=method, rtol=rtol, atol=atol)
    return sol.y[0], sol.y[1]


def solve_naive(t_eval, omega=1.0, eps=0.0, alpha=0.0, A0=1.0, theta0=0.0,
                method="DOP853", rtol=1e-10, atol=1e-12):
    """Naive first-order perturbative solution ``phi0 + eps phi1_D + alpha phi1_V``.

    Exhibits secular growth; included as a baseline for the comparison plots.
    """
    phi = bg.phi0(t_eval, A0, theta0, omega)
    pi = bg.Pi0(t_eval, A0, theta0, omega)
    if eps != 0.0:
        _, p1, q1 = pert.solve_first_order(
            pert.source_duffing, A0, omega=omega, theta=theta0,
            t_eval=t_eval, method=method, rtol=rtol, atol=atol)
        phi = phi + eps * p1
        pi = pi + eps * q1
    if alpha != 0.0:
        _, p1, q1 = pert.solve_first_order(
            pert.source_vdp, A0, omega=omega, theta=theta0,
            t_eval=t_eval, method=method, rtol=rtol, atol=atol)
        phi = phi + alpha * p1
        pi = pi + alpha * q1
    return phi, pi
