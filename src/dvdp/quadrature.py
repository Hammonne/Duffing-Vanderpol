"""Numerical quadrature rules used by the cost-function inner products.

The renormalization cost function is ``I = sum_A integral |phi1_perp|^2 dt``.
Minimising it requires inner products of the form ``integral f(t) g(t) dt``.
We expose three interchangeable rules so their effect on the extracted
generators can be compared (they agree to many digits for smooth integrands):

* :func:`trapezoid`        -- fast, second order.
* :func:`composite_simpson`-- fourth order (needs an even number of intervals).
* :func:`gauss_legendre`   -- spectral accuracy via cubic interpolation.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d


def trapezoid(y, x):
    """Composite trapezoidal rule."""
    return np.trapz(y, x)


def composite_simpson(y, x):
    """Classical composite Simpson rule (requires even number of intervals)."""
    n = len(x) - 1
    if n % 2 != 0:
        raise ValueError("composite Simpson requires an even number of intervals")
    h = (x[-1] - x[0]) / n
    S = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2])
    return S * h / 3.0


def gauss_legendre(y, x, n_nodes=64):
    """Gauss-Legendre quadrature using cubic interpolation of the samples."""
    interp = interp1d(x, y, kind="cubic", bounds_error=False,
                      fill_value="extrapolate")
    nodes, weights = leggauss(n_nodes)
    a, b = x[0], x[-1]
    pts = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    return 0.5 * (b - a) * np.dot(weights, interp(pts))


#: registry so callers can select a rule by name
RULES = {
    "trapz": trapezoid,
    "simpson": composite_simpson,
    "gauss": lambda y, x: gauss_legendre(y, x, n_nodes=64),
}
