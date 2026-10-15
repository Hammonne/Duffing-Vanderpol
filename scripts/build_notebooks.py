"""Generate the clean, reproducible notebooks in notebooks/ from source.

Each notebook is thin: it imports the ``dvdp`` package and exercises one case,
so the physics lives in ``src/`` (single source of truth) and the notebooks are
narrative + figures.  Run:  python scripts/build_notebooks.py
"""

import os

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(__file__)
NBDIR = os.path.join(HERE, "..", "notebooks")
os.makedirs(NBDIR, exist_ok=True)

BOOT = (
    "import sys, os\n"
    "sys.path.insert(0, os.path.abspath('../src'))\n"
    "import numpy as np\n"
    "import matplotlib.pyplot as plt\n"
    "from dvdp import renormalization as rg, integrators as it, "
    "perturbation as pert, metrics as mt\n"
)


def write(name, cells):
    nb = new_notebook(cells=cells, metadata={
        "kernelspec": {"name": "python3", "display_name": "Python 3",
                       "language": "python"},
        "language_info": {"name": "python"}})
    path = os.path.join(NBDIR, name)
    with open(path, "w", encoding="utf-8") as fh:
        nbf.write(nb, fh)
    print("wrote", os.path.relpath(path))


# 01 - perturbation vs exact / secular growth
write("01_perturbation_vs_exact.ipynb", [
    new_markdown_cell(
        "# 01 - Perturbation theory vs. the exact solution\n\n"
        "We show that the naive first-order perturbative correction to the "
        "harmonic oscillator develops **secular growth** (unbounded amplitude) "
        "for both the Duffing and Van der Pol deformations, motivating the "
        "renormalization procedure. Mirrors Sec. II.A of the proposal."),
    new_code_cell(BOOT),
    new_code_cell(
        "t = np.linspace(0, 80, 3000)\n"
        "_, phi1_D, _ = pert.solve_first_order(pert.source_duffing, 1.0, t_eval=t)\n"
        "_, phi1_V, _ = pert.solve_first_order(pert.source_vdp, 1.0, t_eval=t)\n"
        "plt.figure(figsize=(8,4.5))\n"
        "plt.plot(t, phi1_D, label='Duffing $\\\\varphi^{(1)}$')\n"
        "plt.plot(t, phi1_V, label='Van der Pol $\\\\varphi^{(1)}$')\n"
        "plt.xlabel('t'); plt.ylabel('$\\\\varphi^{(1)}$'); plt.legend(); plt.grid(alpha=.3)\n"
        "plt.title('Secular growth of the first-order correction'); plt.show()"),
    new_markdown_cell(
        "The corrections grow without bound: adding `eps*phi1` to the background "
        "cannot approximate the full solution for long times. Renormalization "
        "reabsorbs this growth into slowly varying amplitude/phase."),
])

# 02 - Duffing renormalization
write("02_duffing_renormalization.ipynb", [
    new_markdown_cell(
        "# 02 - Duffing renormalization (frequency shift)\n\n"
        "Duffing deforms the oscillator with a cubic term. The renormalization "
        "reabsorbs the secular growth into a **frequency shift** "
        "$\\beta_\\theta(A)=\\tfrac{3}{8}A^2$ (Bender-Orszag), leaving the "
        "amplitude unchanged ($\\beta_A\\approx0$). Mirrors Sec. II.B / Figs. 1, 3."),
    new_code_cell(BOOT),
    new_code_cell(
        "A = np.linspace(0.05, 3.0, 40)\n"
        "bA, bth = rg.extract_generators('duffing', A)\n"
        "plt.figure(figsize=(7,4.5))\n"
        "plt.plot(A, bth, 'C3', label=r'$\\beta_\\theta(A)$ (numerical)')\n"
        "plt.plot(A, 0.375*A**2, 'k--', label=r'$\\frac38 A^2$ (Bender-Orszag)')\n"
        "plt.plot(A, bA, 'C0', label=r'$\\beta_A(A)\\approx0$')\n"
        "plt.xlabel('A'); plt.legend(); plt.grid(alpha=.3)\n"
        "plt.title('Duffing generators'); plt.show()\n"
        "print('beta_theta/A^2 at A=1:', np.interp(1.0, A, bth))"),
    new_code_cell(
        "t = np.linspace(0, 100, 4000); eps = 0.1\n"
        "pf,_ = it.solve_full(t, eps=eps, A0=1.0)\n"
        "pn,_ = it.solve_naive(t, eps=eps, A0=1.0)\n"
        "res = rg.renormalize(t, eps=eps, A0=1.0)\n"
        "plt.figure(figsize=(9,4.5))\n"
        "plt.plot(t, pn, 'k--', lw=1, label='naive')\n"
        "plt.plot(t, pf, 'C0', lw=1.3, label='full')\n"
        "plt.plot(t, res.phi_R, 'C3', lw=1.3, label='renormalized')\n"
        "plt.ylim(-4,4); plt.xlabel('t'); plt.ylabel('$\\\\varphi$'); plt.legend(); plt.grid(alpha=.3)\n"
        "plt.title(f'Duffing eps={eps}'); plt.show()\n"
        "print('RMS  RG=%.4f  naive=%.4f' % (mt.rms_error(res.phi_R,pf), mt.rms_error(pn,pf)))"),
])

# 03 - Van der Pol renormalization
write("03_vanderpol_renormalization.ipynb", [
    new_markdown_cell(
        "# 03 - Van der Pol renormalization (amplitude flow)\n\n"
        "Van der Pol requires an **amplitude** flow $\\beta_A(A)=\\tfrac12A(1-A^2/4)$ "
        "(zero at the limit cycle $A=2$), with negligible frequency shift. "
        "Mirrors Sec. II.C / Fig. 4."),
    new_code_cell(BOOT),
    new_code_cell(
        "A = np.linspace(0.05, 3.0, 40)\n"
        "bA, bth = rg.extract_generators('vdp', A)\n"
        "plt.figure(figsize=(7,4.5))\n"
        "plt.plot(A, bA, 'C0', label=r'$\\beta_A(A)$ (numerical)')\n"
        "plt.plot(A, 0.5*A*(1-A**2/4), 'k--', label=r'$\\frac12A(1-A^2/4)$')\n"
        "plt.axhline(0, color='gray', lw=.5); plt.axvline(2, color='gray', ls=':', lw=.8)\n"
        "plt.xlabel('A'); plt.legend(); plt.grid(alpha=.3); plt.title('Van der Pol generators'); plt.show()"),
    new_code_cell(
        "t = np.linspace(0, 100, 4000); alpha = 0.1\n"
        "pf,qf = it.solve_full(t, alpha=alpha, A0=1.0)\n"
        "res = rg.renormalize(t, alpha=alpha, A0=1.0)\n"
        "fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5))\n"
        "a1.plot(t, pf,'C0',label='full'); a1.plot(t,res.phi_R,'C3',label='renorm'); a1.legend(); a1.grid(alpha=.3); a1.set_xlabel('t'); a1.set_title(f'VdP alpha={alpha}')\n"
        "a2.plot(pf,qf,'C0',lw=1,label='full'); a2.plot(res.phi_R,res.pi_R,'C3',label='renorm'); a2.axis('equal'); a2.legend(); a2.grid(alpha=.3); a2.set_title('phase space (limit cycle)')\n"
        "plt.show(); print('A final =', res.A_t[-1])"),
])

# 04 - combined DVdP
write("04_combined_dvdp.ipynb", [
    new_markdown_cell(
        "# 04 - Combined Duffing-Van der Pol (the universal method)\n\n"
        "The master equation\n\n"
        "$$\\dot\\varphi=\\Pi,\\quad \\dot\\Pi=-\\omega^2\\varphi-\\epsilon\\varphi^3"
        "+\\alpha(1-\\varphi^2)\\Pi$$\n\n"
        "is renormalized with a **single universal engine**. Because the "
        "first-order equation is linear, the generators superpose:\n"
        "$\\beta=\\epsilon\\,\\beta^{\\rm Duf}+\\alpha\\,\\beta^{\\rm VdP}$. "
        "Duffing contributes a frequency shift, Van der Pol an amplitude flow; "
        "both act simultaneously."),
    new_code_cell(BOOT),
    new_code_cell(
        "def show(eps, alpha, tmax=100):\n"
        "    t = np.linspace(0, tmax, 4000)\n"
        "    pf,qf = it.solve_full(t, eps=eps, alpha=alpha, A0=1.0)\n"
        "    pn,_  = it.solve_naive(t, eps=eps, alpha=alpha, A0=1.0)\n"
        "    res   = rg.renormalize(t, eps=eps, alpha=alpha, A0=1.0)\n"
        "    fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5))\n"
        "    a1.plot(t,pn,'k--',lw=.8,alpha=.6,label='naive'); a1.plot(t,pf,'C0',lw=1.2,label='full'); a1.plot(t,res.phi_R,'C3',lw=1.2,label='renorm')\n"
        "    a1.set_ylim(-3.2,3.2); a1.legend(); a1.grid(alpha=.3); a1.set_xlabel('t'); a1.set_title(f'DVdP eps={eps}, alpha={alpha}')\n"
        "    a2.plot(pf,qf,'C0',lw=1,label='full'); a2.plot(res.phi_R,res.pi_R,'C3',label='renorm'); a2.axis('equal'); a2.legend(); a2.grid(alpha=.3); a2.set_title('phase space')\n"
        "    plt.show(); print('RMS  RG=%.4f  naive=%.4f' % (mt.rms_error(res.phi_R,pf), mt.rms_error(pn,pf)))\n"
        "\n"
        "show(0.1, 0.1)   # soft coupling\n"
        "show(1.0, 1.0)   # strong coupling"),
    new_markdown_cell(
        "`eps=0` recovers pure Van der Pol and `alpha=0` recovers pure Duffing, "
        "so this one engine solves all three problems. For soft coupling the "
        "approximation is excellent; for strong coupling it captures the "
        "amplitude envelope and limit-cycle size, with residual phase drift that "
        "is the expected first-order limitation (the true strong-coupling limit "
        "cycle is anharmonic)."),
])

# 05 - benchmark
write("05_numerical_method_benchmark.ipynb", [
    new_markdown_cell(
        "# 05 - Numerical-method benchmark\n\n"
        "Choice of the integrator for the full (reference) solution. We compare "
        "fixed-step RK4 against scipy's RK45, DOP853 and Radau on Duffing, Van "
        "der Pol and the combined system. See `docs/benchmark.md` for the full "
        "table; **DOP853** is adopted as the default."),
    new_code_cell(BOOT + "from dvdp.systems import full_rhs\nfrom scipy.integrate import solve_ivp\n"),
    new_code_cell(
        "t = np.linspace(0, 50, 2001)\n"
        "ref = solve_ivp(full_rhs,(0,50),[1,0],args=(1.0,1.0,1.0),t_eval=t,method='DOP853',rtol=1e-13,atol=1e-13).y[0]\n"
        "for m in ['RK45','DOP853','Radau']:\n"
        "    s = solve_ivp(full_rhs,(0,50),[1,0],args=(1.0,1.0,1.0),t_eval=t,method=m,rtol=1e-10,atol=1e-12)\n"
        "    print(f'{m:8s} rms={np.linalg.norm(s.y[0]-ref)/np.sqrt(len(ref)):.2e}  nfev={s.nfev}')"),
])

print("done")
