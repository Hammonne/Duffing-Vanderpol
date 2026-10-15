"""dvdp: Numerical renormalization for the Duffing-Van der Pol family.

A single, *universal* renormalization engine that approximates the full
(nonlinear) solution of

    dphi/dt = Pi
    dPi/dt  = -omega^2 phi - eps*phi^3 + alpha*(1 - phi^2)*Pi

for arbitrary deformation strengths ``eps`` (Duffing) and ``alpha`` (Van der
Pol), independently.  Setting ``eps=0`` recovers pure Van der Pol, ``alpha=0``
recovers pure Duffing, and both nonzero gives the combined system.

Modules
-------
systems         : full nonlinear right-hand sides.
background       : harmonic background phi0 = A cos(omega t + theta) and its
                   parameter derivatives.
perturbation     : linear first-order perturbation equations with background
                   source terms (separable into Duffing and VdP pieces).
renormalization  : the universal engine -- extracts beta generators and
                   integrates the parameter flow to build phi_R, Pi_R.
integrators      : fixed-step RK4 plus scipy solve_ivp wrappers (benchmark).
quadrature       : trapezoid / composite Simpson / Gauss-Legendre rules.
metrics          : RMS error and related diagnostics.
"""

from .systems import full_rhs, duffing_rhs, vdp_rhs, dvdp_rhs
from .renormalization import (
    RenormResult,
    renormalize,
    extract_generators,
)

__all__ = [
    "full_rhs",
    "duffing_rhs",
    "vdp_rhs",
    "dvdp_rhs",
    "RenormResult",
    "renormalize",
    "extract_generators",
]
