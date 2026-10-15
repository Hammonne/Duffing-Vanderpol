# Numerical-method benchmark (Phase 5)

Reference: DOP853 at rtol=atol=1e-13. Grid t in [0,50], 2001 points.

| System | Method | RMS error | RHS evals | wall time (s) |
|---|---|---|---|---|
| Duffing (eps=1) | RK45 | 1.83e-09 | 13976 | 0.718 |
| Duffing (eps=1) | DOP853 | 5.63e-10 | 6023 | 0.164 |
| Duffing (eps=1) | Radau | 2.41e-12 | 63990 | 4.256 |
| Duffing (eps=1) | RK4(fixed) | 7.08e-07 | 8004 | 0.134 |
| Van der Pol (alpha=1) | RK45 | 4.66e-10 | 15872 | 0.626 |
| Van der Pol (alpha=1) | DOP853 | 3.83e-11 | 9479 | 0.368 |
| Van der Pol (alpha=1) | Radau | 5.18e-12 | 86566 | 7.734 |
| Van der Pol (alpha=1) | RK4(fixed) | 1.91e-06 | 8004 | 0.099 |
| DVdP (eps=alpha=1) | RK45 | 1.20e-09 | 21620 | 0.841 |
| DVdP (eps=alpha=1) | DOP853 | 4.19e-10 | 13163 | 0.425 |
| DVdP (eps=alpha=1) | Radau | 5.02e-12 | 145952 | 10.641 |
| DVdP (eps=alpha=1) | RK4(fixed) | 5.38e-06 | 8004 | 0.101 |

## Decision

DOP853 gives the best accuracy-per-evaluation on the non-stiff Duffing and combined cases and matches the limit cycle of Van der Pol; Radau is the safe fallback if stiffness appears at very large alpha. The fixed-step RK4 is retained only for pedagogical comparison. **Default: DOP853** (rtol=1e-10, atol=1e-12), used throughout `dvdp.integrators.solve_full`.
