"""Error metrics for comparing approximate solutions against the full one."""

from __future__ import annotations

import numpy as np


def rms_error(approx, reference):
    """Root-mean-square error ``||approx - reference|| / sqrt(N)``."""
    approx = np.asarray(approx)
    reference = np.asarray(reference)
    return np.linalg.norm(approx - reference) / np.sqrt(len(reference))


def max_error(approx, reference):
    """Maximum absolute (sup-norm) error."""
    return float(np.max(np.abs(np.asarray(approx) - np.asarray(reference))))


def relative_rms(approx, reference):
    """RMS error normalised by the RMS amplitude of the reference."""
    reference = np.asarray(reference)
    denom = np.linalg.norm(reference) / np.sqrt(len(reference))
    if denom == 0:
        return np.nan
    return rms_error(approx, reference) / denom
