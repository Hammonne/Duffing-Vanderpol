"""The universal numerical-renormalization engine.

Given the two independent deformation strengths ``eps`` (Duffing) and
``alpha`` (Van der Pol), this module builds a renormalized approximation to
the full solution of

    dphi/dt = Pi,  dPi/dt = -omega^2 phi - eps phi^3 + alpha (1 - phi^2) Pi.

Algorithm (mirrors Sec. II.B of the AFOSR proposal and Galvez-Ghersi & Stein
2021, generalised to two couplings):

1.  For each deformation D in {duffing, vdp} and each amplitude A on a grid,
    solve the *linear* first-order perturbation driven by the background
    source ``S_D`` and project its secular growth onto the tangent directions
    ``t * d phi0/d A`` and ``t * d phi0/d theta`` by minimising the quadratic
    cost ``I = integral |phi1_perp|^2 dt``.  This yields the per-unit-coupling
    generators ``beta_A^D(A)`` and ``beta_theta^D(A)``.

2.  Because the first-order equation is linear, the generators of the combined
    system are the linear superposition

        beta_A(A)     = eps beta_A^D(A)     + alpha beta_A^V(A)
        beta_theta(A) = eps beta_theta^D(A) + alpha beta_theta^V(A).

3.  Integrate the (generally nonlinear) parameter-flow ODE

        dA/dt = beta_A(A),   dtheta/dt = beta_theta(A)

    and reconstruct the renormalized solution
    ``phi_R(t) = A(t) cos(omega t + theta(t))``.

Solving the flow as a full ODE (rather than truncating to a linear-in-Delta t
shift) is what lets the method remain accurate for *strong* deformations, not
only the infinitesimal regime.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from . import background as bg
from . import perturbation as pert
from .quadrature import RULES


# ---------------------------------------------------------------------------
# Step 1-2: generator extraction
# ---------------------------------------------------------------------------

def _project_generators(t, phi1, pi1, A, omega, theta, integrate=None):
    """Extract the secular generators (beta_A, beta_theta) robustly.

    Rather than integrating against t-weighted directions over a finite window
    (which lets the *bounded* oscillatory part of the perturbation contaminate
    the result), we first decompose the perturbation *pointwise* onto the
    instantaneous tangent directions

        v_A     = ( d phi0/d A    , d Pi0/d A    )
        v_theta = ( d phi0/d theta, d Pi0/d theta )

    giving coefficient time series ``c_A(t)``, ``c_theta(t)``.  A resonant
    (secular) generator makes these grow linearly in time, while any bounded
    response contributes only a bounded offset.  The generators are therefore
    the *slopes* of ``c_A(t)`` and ``c_theta(t)``, obtained by a straight-line
    fit that discards the constant offset.

    The ``integrate`` argument is accepted for API compatibility but unused.
    """
    vA = np.array([bg.dphi0_dA(t, A, theta, omega),
                   bg.dPi0_dA(t, A, theta, omega)])          # shape (2, N)
    vth = np.array([bg.dphi0_dtheta(t, A, theta, omega),
                    bg.dPi0_dtheta(t, A, theta, omega)])     # shape (2, N)

    # Pointwise solve [vA vth] @ [c_A, c_theta]^T = [phi1, pi1]^T for every t.
    # Build the 2x2 system per sample and invert analytically.
    a11, a21 = vA[0], vA[1]
    a12, a22 = vth[0], vth[1]
    det = a11 * a22 - a12 * a21
    # Guard the (measure-zero) instants where the tangent vectors are parallel.
    eps_det = 1e-12 * (np.abs(a11 * a22) + np.abs(a12 * a21) + 1.0)
    det = np.where(np.abs(det) < eps_det, np.sign(det) * eps_det + eps_det, det)
    c_A = (a22 * phi1 - a12 * pi1) / det
    c_theta = (-a21 * phi1 + a11 * pi1) / det

    # Slope of each coefficient series (linear fit, intercept absorbs offset).
    # Use the second half of the window to avoid the initial transient.
    half = len(t) // 2
    sl = slice(half, None)
    beta_A = np.polyfit(t[sl], c_A[sl], 1)[0]
    beta_theta = np.polyfit(t[sl], c_theta[sl], 1)[0]
    return beta_A, beta_theta


def extract_generators(deform, A_grid, omega=1.0, theta=0.0, rule="gauss",
                        n_periods=6, n_points=2000, pert_method="DOP853"):
    """Extract per-unit-coupling generators for one deformation.

    Parameters
    ----------
    deform : {"duffing", "vdp"}
        Which background source to use.
    A_grid : ndarray
        Amplitudes at which the generators are sampled.
    omega, theta : float
        Background frequency and phase.
    rule : {"trapz", "simpson", "gauss"}
        Quadrature rule for the cost-function inner products.
    n_periods, n_points, pert_method :
        Passed to the first-order solver.

    Returns
    -------
    beta_A : ndarray, beta_theta : ndarray
        Generator samples on ``A_grid`` (per unit of the deformation coupling).
    """
    source = pert.SOURCES[deform]
    integrate = RULES[rule]

    beta_A = np.empty_like(A_grid, dtype=float)
    beta_theta = np.empty_like(A_grid, dtype=float)

    for i, A in enumerate(A_grid):
        t, phi1, pi1 = pert.solve_first_order(
            source, A, omega=omega, theta=theta,
            n_periods=n_periods, n_points=n_points, method=pert_method)
        bA, bth = _project_generators(t, phi1, pi1, A, omega, theta, integrate)
        beta_A[i] = bA
        beta_theta[i] = bth

    return beta_A, beta_theta


# ---------------------------------------------------------------------------
# Step 3: parameter flow + reconstruction
# ---------------------------------------------------------------------------

@dataclass
class RenormResult:
    """Container for a renormalization run."""

    t: np.ndarray
    phi_R: np.ndarray
    pi_R: np.ndarray
    A_t: np.ndarray
    theta_t: np.ndarray
    A_grid: np.ndarray
    beta_A: np.ndarray           # combined amplitude generator on A_grid
    beta_theta: np.ndarray       # combined phase (frequency-shift) generator
    components: dict = field(default_factory=dict)  # per-deformation samples

    def frequency_shift(self, A=None):
        """Phase generator beta_theta (= instantaneous frequency shift)."""
        if A is None:
            return self.beta_theta
        return interp1d(self.A_grid, self.beta_theta, kind="cubic",
                        fill_value="extrapolate")(A)


def renormalize(t_eval, omega=1.0, eps=0.0, alpha=0.0, A0=1.0, theta0=0.0,
                A_grid=None, rule="gauss", n_periods=6, n_points=2000,
                flow_method="DOP853", pert_method="DOP853"):
    """Build the renormalized solution of the combined system.

    Parameters
    ----------
    t_eval : ndarray
        Output time grid for the renormalized solution.
    omega : float
        Natural frequency.
    eps, alpha : float
        Duffing and Van der Pol deformation strengths (independent).
    A0, theta0 : float
        Initial amplitude and phase (``phi(0)=A0 cos(theta0)``).
    A_grid : ndarray, optional
        Amplitude grid for generator extraction.  Defaults to a grid spanning
        a little below ``A0`` up to the Van der Pol limit-cycle radius (~2).
    rule, n_periods, n_points, flow_method, pert_method :
        Numerical knobs (see :func:`extract_generators`).

    Returns
    -------
    RenormResult
    """
    if A_grid is None:
        A_hi = max(4.0, 1.5 * A0)
        A_grid = np.linspace(0.05, A_hi, 30)

    components = {}
    beta_A = np.zeros_like(A_grid, dtype=float)
    beta_theta = np.zeros_like(A_grid, dtype=float)

    if eps != 0.0:
        bA_D, bth_D = extract_generators("duffing", A_grid, omega=omega,
                                         rule=rule, n_periods=n_periods,
                                         n_points=n_points,
                                         pert_method=pert_method)
        components["duffing"] = (bA_D, bth_D)
        beta_A += eps * bA_D
        beta_theta += eps * bth_D

    if alpha != 0.0:
        bA_V, bth_V = extract_generators("vdp", A_grid, omega=omega,
                                         rule=rule, n_periods=n_periods,
                                         n_points=n_points,
                                         pert_method=pert_method)
        components["vdp"] = (bA_V, bth_V)
        beta_A += alpha * bA_V
        beta_theta += alpha * bth_V

    beta_A_fn = interp1d(A_grid, beta_A, kind="cubic", fill_value="extrapolate")
    beta_theta_fn = interp1d(A_grid, beta_theta, kind="cubic",
                             fill_value="extrapolate")

    # Parameter flow: dA/dt = beta_A(A), dtheta/dt = beta_theta(A).
    def flow(t, y):
        A, _theta = y
        return [float(beta_A_fn(A)), float(beta_theta_fn(A))]

    sol = solve_ivp(flow, (t_eval[0], t_eval[-1]), [A0, theta0],
                    t_eval=t_eval, method=flow_method, rtol=1e-9, atol=1e-12)
    A_t = sol.y[0]
    theta_t = sol.y[1]

    # Reconstruct phi_R and its exact time derivative pi_R.
    arg = omega * t_eval + theta_t
    dA_dt = beta_A_fn(A_t)
    dtheta_dt = beta_theta_fn(A_t)
    phi_R = A_t * np.cos(arg)
    pi_R = dA_dt * np.cos(arg) - A_t * (omega + dtheta_dt) * np.sin(arg)

    return RenormResult(t=t_eval, phi_R=phi_R, pi_R=pi_R, A_t=A_t,
                        theta_t=theta_t, A_grid=A_grid, beta_A=beta_A,
                        beta_theta=beta_theta, components=components)
