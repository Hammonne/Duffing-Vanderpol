"""Generate every figure used in the paper (Phases 4 and 6).

Reproduces Figs. 1-4 of the AFOSR proposal (Duffing and Van der Pol) and
produces the new combined Duffing-Van der Pol figures, plus the extracted
beta-generator functions.  All figures are written to figures/ as PDF.

Run:  python scripts/make_figures.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from dvdp import renormalization as rg       # noqa: E402
from dvdp import integrators as it           # noqa: E402
from dvdp import perturbation as pert        # noqa: E402
from dvdp import metrics as mt               # noqa: E402

FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 11, "figure.dpi": 120})
C_FULL, C_RG, C_NAIVE, C_LIN = "C0", "C3", "k", "C2"


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print("wrote", os.path.relpath(path))


# ---------------------------------------------------------------------------
def fig_secular_growth():
    """Naive first-order perturbation develops secular growth."""
    t = np.linspace(0, 80, 3000)
    _, phi1_D, _ = pert.solve_first_order(pert.source_duffing, 1.0, t_eval=t)
    _, phi1_V, _ = pert.solve_first_order(pert.source_vdp, 1.0, t_eval=t)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t, phi1_D, label=r"Duffing $\varphi^{(1)}$")
    ax.plot(t, phi1_V, label=r"Van der Pol $\varphi^{(1)}$")
    ax.set_xlabel("t"); ax.set_ylabel(r"$\varphi^{(1)}(t)$")
    ax.set_title("First-order corrections exhibit secular growth")
    ax.grid(alpha=0.3); ax.legend()
    save(fig, "fig_secular_growth.pdf")


def _compare_panel(ax, t, phi_full, phi_naive, phi_rg, title, show_lin=True):
    if show_lin:
        ax.plot(t, np.cos(t), ":", color=C_LIN, lw=1.2, label=r"$\varphi^{(0)}$ ($\epsilon=0$)")
    ax.plot(t, phi_naive, "--", color=C_NAIVE, lw=1, label=r"$\varphi^{(0)}+\varphi^{(1)}$ (naive)")
    ax.plot(t, phi_full, "-", color=C_FULL, lw=1.4, label=r"$\varphi_{\rm full}$")
    ax.plot(t, phi_rg, "-", color=C_RG, lw=1.4, label=r"$\varphi_R$ (renorm.)")
    ax.set_xlabel("t"); ax.set_ylabel(r"$\varphi(t)$"); ax.set_title(title)
    ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=8)


def fig_duffing():
    """Fig. 1 analogue: Duffing full vs naive vs renormalized vs linear."""
    t = np.linspace(0, 100, 4000)
    eps = 0.1
    pf, _ = it.solve_full(t, eps=eps, A0=1.0)
    pn, _ = it.solve_naive(t, eps=eps, A0=1.0)
    res = rg.renormalize(t, eps=eps, A0=1.0)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    _compare_panel(ax, t, pf, pn, res.phi_R,
                   rf"Duffing, $\epsilon={eps}$  "
                   rf"(RMS$_R$={mt.rms_error(res.phi_R,pf):.3f}, "
                   rf"RMS$_{{\rm naive}}$={mt.rms_error(pn,pf):.2f})")
    ax.set_ylim(-4, 4)
    save(fig, "fig_duffing_compare.pdf")


def fig_vdp():
    """Fig. 4 analogue: Van der Pol time series + phase space."""
    t = np.linspace(0, 100, 4000)
    alpha = 0.1
    pf, qf = it.solve_full(t, alpha=alpha, A0=1.0)
    pn, qn = it.solve_naive(t, alpha=alpha, A0=1.0)
    res = rg.renormalize(t, alpha=alpha, A0=1.0)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    _compare_panel(a1, t, pf, pn, res.phi_R,
                   rf"Van der Pol, $\alpha={alpha}$ "
                   rf"(RMS$_R$={mt.rms_error(res.phi_R,pf):.3f})",
                   show_lin=False)
    a1.set_ylim(-6, 6)
    a2.plot(pn, qn, "--", color=C_NAIVE, alpha=0.4, lw=0.8, label="naive")
    a2.plot(pf, qf, "-", color=C_FULL, lw=1.0, label="full")
    a2.plot(res.phi_R, res.pi_R, "-", color=C_RG, lw=1.4, label="renorm.")
    a2.set_xlabel(r"$\varphi$"); a2.set_ylabel(r"$\Pi$")
    a2.set_title("Phase space (limit cycle)"); a2.axis("equal")
    a2.grid(alpha=0.3); a2.legend(fontsize=8)
    save(fig, "fig_vdp_compare.pdf")


def fig_combined():
    """Combined DVdP: soft and strong coupling, time series + phase space."""
    for tag, (eps, alpha) in {"soft": (0.1, 0.1), "strong": (1.0, 1.0)}.items():
        t = np.linspace(0, 100, 4000)
        pf, qf = it.solve_full(t, eps=eps, alpha=alpha, A0=1.0)
        pn, qn = it.solve_naive(t, eps=eps, alpha=alpha, A0=1.0)
        res = rg.renormalize(t, eps=eps, alpha=alpha, A0=1.0)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
        _compare_panel(a1, t, pf, pn, res.phi_R,
                       rf"DVdP, $\epsilon=\alpha={eps}$ "
                       rf"(RMS$_R$={mt.rms_error(res.phi_R,pf):.3f}, "
                       rf"RMS$_{{\rm naive}}$={mt.rms_error(pn,pf):.2f})",
                       show_lin=False)
        ymax = max(3.0, 1.2 * np.max(np.abs(pf)))
        a1.set_ylim(-ymax, ymax)
        a2.plot(pf, qf, "-", color=C_FULL, lw=1.0, label="full")
        a2.plot(res.phi_R, res.pi_R, "-", color=C_RG, lw=1.4, label="renorm.")
        a2.set_xlabel(r"$\varphi$"); a2.set_ylabel(r"$\Pi$")
        a2.set_title("Phase space"); a2.axis("equal")
        a2.grid(alpha=0.3); a2.legend(fontsize=8)
        save(fig, f"fig_dvdp_{tag}.pdf")


def fig_generators():
    """Extracted beta-generators for Duffing and Van der Pol."""
    A = np.linspace(0.05, 3.0, 40)
    bA_D, bth_D = rg.extract_generators("duffing", A)
    bA_V, bth_V = rg.extract_generators("vdp", A)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.plot(A, bth_D, "C3", label=r"$\beta_\theta^{\rm Duf}(A)$")
    a1.plot(A, 0.375 * A**2, "k--", lw=1, label=r"$\frac{3}{8}A^2$ (Bender-Orszag)")
    a1.plot(A, bA_D, "C0", label=r"$\beta_A^{\rm Duf}(A)\approx 0$")
    a1.set_title("Duffing generators"); a1.set_xlabel("A"); a1.grid(alpha=0.3); a1.legend()
    a2.plot(A, bA_V, "C0", label=r"$\beta_A^{\rm VdP}(A)$")
    a2.plot(A, 0.5 * A * (1 - A**2 / 4), "k--", lw=1, label=r"$\frac{1}{2}A(1-A^2/4)$")
    a2.plot(A, bth_V, "C3", label=r"$\beta_\theta^{\rm VdP}(A)\approx 0$")
    a2.axhline(0, color="gray", lw=0.5)
    a2.set_title("Van der Pol generators"); a2.set_xlabel("A"); a2.grid(alpha=0.3); a2.legend()
    save(fig, "fig_generators.pdf")


def fig_error_scaling():
    """RMS error of RG vs naive as a function of coupling strength."""
    ts = np.linspace(0, 100, 4000)
    strengths = np.linspace(0.05, 1.0, 12)
    for label, kw in {"Duffing": dict(eps=None), "VdP": dict(alpha=None)}.items():
        pass
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, mode in [("Duffing", "eps"), ("VdP", "alpha"), ("DVdP", "both")]:
        e_rg, e_nv = [], []
        for s in strengths:
            kw = {"eps": s} if mode == "eps" else {"alpha": s} if mode == "alpha" else {"eps": s, "alpha": s}
            pf, _ = it.solve_full(ts, A0=1.0, **kw)
            pn, _ = it.solve_naive(ts, A0=1.0, **kw)
            res = rg.renormalize(ts, A0=1.0, **kw)
            e_rg.append(mt.rms_error(res.phi_R, pf))
            e_nv.append(mt.rms_error(pn, pf))
        ax.semilogy(strengths, e_rg, "o-", label=f"{label} RG")
        ax.semilogy(strengths, e_nv, "x--", alpha=0.5, label=f"{label} naive")
    ax.set_xlabel("deformation strength"); ax.set_ylabel("RMS error vs full")
    ax.set_title("Renormalized vs naive accuracy"); ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, ncol=3)
    save(fig, "fig_error_scaling.pdf")


if __name__ == "__main__":
    fig_secular_growth()
    fig_duffing()
    fig_vdp()
    fig_combined()
    fig_generators()
    fig_error_scaling()
    print("All figures done.")
