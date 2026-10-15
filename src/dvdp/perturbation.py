"""First-order perturbation equations with background source terms.

Substituting the ansatz ``phi ~ phi0 + eps*phi1_D + alpha*phi1_V`` into the
master equation and collecting the terms linear in the couplings yields a
*linear* driven oscillator for each correction:

    dphi1/dt = Pi1
    dPi1/dt  = -omega**2 * phi1 + S(t)

where the source ``S`` is evaluated on the *background* (phi0, Pi0):

    Duffing     :  S_D(t) = -(phi0)**3
    Van der Pol :  S_V(t) =  (1 - phi0**2) * Pi0

Because the equations are linear, the correction for the combined system is
the linear superposition ``phi1 = eps*phi1_D + alpha*phi1_V``.  This is the
key fact that makes a single renormalization engine *universal*: the Duffing
and Van der Pol generators can be extracted independently and then combined.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import background as bg


def source_duffing(t, A, theta=0.0, omega=1.0):
    """Duffing source term ``-(phi0)**3`` evaluated on the background."""
    p0 = bg.phi0(t, A, theta, omega)
    return -p0**3


def source_vdp(t, A, theta=0.0, omega=1.0):
    """Van der Pol source term ``(1 - phi0**2) * Pi0`` on the background."""
    p0 = bg.phi0(t, A, theta, omega)
    pi0 = bg.Pi0(t, A, theta, omega)
    return (1.0 - p0**2) * pi0


#: registry mapping a deformation name to its background source function.
SOURCES = {"duffing": source_duffing, "vdp": source_vdp}


def solve_first_order(source, A, omega=1.0, theta=0.0, t_eval=None,
                      n_periods=6, n_points=2000, method="DOP853",
                      rtol=1e-9, atol=1e-12):
    """Integrate the linear first-order perturbation for a given source.

    Parameters
    ----------
    source : callable
        ``source(t, A, theta, omega) -> ndarray`` background forcing term.
    A : float
        Background amplitude the correction is computed about.
    omega, theta : float
        Background frequency and phase.
    t_eval : ndarray, optional
        Time grid.  If ``None``, a grid of ``n_periods`` periods with
        ``n_points`` samples is used (enough to expose the secular slope).
    n_periods, n_points : int
        Used to build the default grid.
    method, rtol, atol :
        Passed to :func:`scipy.integrate.solve_ivp`.

    Returns
    -------
    t : ndarray
    phi1 : ndarray
    pi1 : ndarray
    """
    if t_eval is None:
        T = 2.0 * np.pi / omega
        t_eval = np.linspace(0.0, n_periods * T, n_points)

    def rhs(t, y):
        phi1, pi1 = y
        return [pi1, -omega**2 * phi1 + source(t, A, theta, omega)]

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), [0.0, 0.0],
                    t_eval=t_eval, method=method, rtol=rtol, atol=atol)
    return sol.t, sol.y[0], sol.y[1]
