"""Full nonlinear right-hand sides for the Duffing-Van der Pol family.

The master equation (written as a first-order system for phi and its
conjugate momentum Pi = dphi/dt) is

    dphi/dt = Pi
    dPi/dt  = -omega**2 * phi - eps * phi**3 + alpha * (1 - phi**2) * Pi

with two *independent* deformation parameters:

* ``eps``   -- Duffing cubic stiffness  (-eps*phi**3)
* ``alpha`` -- Van der Pol nonlinear damping  (+alpha*(1 - phi**2)*Pi)

Setting ``eps = 0`` gives pure Van der Pol; ``alpha = 0`` gives pure Duffing.

Note (erratum vs. the AFOSR proposal): the proposal's Eq. (17b) prints the
Van der Pol term as ``(1 + phi**2)``; the standard Van der Pol oscillator and
all working notebooks use ``(1 - phi**2)``, which is what we use here.
"""

from __future__ import annotations

import numpy as np


def full_rhs(t, y, omega=1.0, eps=0.0, alpha=0.0):
    """Right-hand side of the combined Duffing-Van der Pol system.

    Parameters
    ----------
    t : float
        Time (unused; system is autonomous, kept for solve_ivp signature).
    y : array_like, shape (2,)
        State ``[phi, Pi]``.
    omega : float
        Natural angular frequency of the undeformed harmonic oscillator.
    eps : float
        Duffing deformation strength.
    alpha : float
        Van der Pol deformation strength.

    Returns
    -------
    list
        ``[dphi/dt, dPi/dt]``.
    """
    phi, Pi = y
    dphi = Pi
    dPi = -omega**2 * phi - eps * phi**3 + alpha * (1.0 - phi**2) * Pi
    return [dphi, dPi]


def duffing_rhs(t, y, omega=1.0, eps=1.0):
    """Pure Duffing oscillator (``alpha = 0``)."""
    return full_rhs(t, y, omega=omega, eps=eps, alpha=0.0)


def vdp_rhs(t, y, omega=1.0, alpha=1.0):
    """Pure Van der Pol oscillator (``eps = 0``)."""
    return full_rhs(t, y, omega=omega, eps=0.0, alpha=alpha)


def dvdp_rhs(t, y, omega=1.0, eps=1.0, alpha=1.0):
    """Combined Duffing-Van der Pol oscillator (both deformations on)."""
    return full_rhs(t, y, omega=omega, eps=eps, alpha=alpha)
