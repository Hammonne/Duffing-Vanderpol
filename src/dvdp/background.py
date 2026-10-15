"""Harmonic background solution and its parameter derivatives.

The background (order eps^0, alpha^0) is the harmonic oscillator written with
an *amplitude* ``A`` and a *phase* ``theta`` as the renormalizable parameters:

    phi0(t) =  A cos(omega t + theta)
    Pi0(t)  = -A omega sin(omega t + theta)

The renormalization group promotes ``A`` and ``theta`` to slowly varying
functions ``A(t)``, ``theta(t)``.  The generators are obtained by projecting
the secular growth of the first-order perturbation onto the *secular tangent
directions* of the background:

    P_A   = t * d phi0 / d A      (amplitude direction, weighted by t)
    P_th  = t * d phi0 / d theta  (phase direction, weighted by t)

and the analogous combinations for the momentum Pi0.  The explicit t-weight
encodes the secular (linear-in-time) character of the resonant growth that the
flow equations ``dA/dt = beta_A`` and ``dtheta/dt = beta_theta`` reabsorb.
"""

from __future__ import annotations

import numpy as np


def phi0(t, A, theta=0.0, omega=1.0):
    """Background displacement ``A cos(omega t + theta)``."""
    return A * np.cos(omega * t + theta)


def Pi0(t, A, theta=0.0, omega=1.0):
    """Background momentum ``-A omega sin(omega t + theta)``."""
    return -A * omega * np.sin(omega * t + theta)


# ---------------------------------------------------------------------------
# Parameter derivatives of the background (the tangent plane T_phi S).
# ---------------------------------------------------------------------------

def dphi0_dA(t, A, theta=0.0, omega=1.0):
    return np.cos(omega * t + theta)


def dPi0_dA(t, A, theta=0.0, omega=1.0):
    return -omega * np.sin(omega * t + theta)


def dphi0_dtheta(t, A, theta=0.0, omega=1.0):
    return -A * np.sin(omega * t + theta)


def dPi0_dtheta(t, A, theta=0.0, omega=1.0):
    return -A * omega * np.cos(omega * t + theta)


def secular_basis(t, A, theta=0.0, omega=1.0):
    """Return the secular tangent directions used for generator fitting.

    Returns
    -------
    (P_A, P_theta) : tuple of ndarray pairs
        Each element is the pair ``(phi-component, Pi-component)`` of the
        t-weighted background derivative, i.e.

            P_A     = ( t * dphi0/dA    , t * dPi0/dA    )
            P_theta = ( t * dphi0/dtheta, t * dPi0/dtheta )
    """
    P_A = (t * dphi0_dA(t, A, theta, omega), t * dPi0_dA(t, A, theta, omega))
    P_theta = (t * dphi0_dtheta(t, A, theta, omega),
               t * dPi0_dtheta(t, A, theta, omega))
    return P_A, P_theta
