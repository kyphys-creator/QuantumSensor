"""Compare tie-break objectives on the degenerate chi^2 = 0 face.

The inverse problem is under-determined: n_vmin (~10^2-10^3) >> n_ebins (5).
The set of monotone non-negative fluxes that reproduce the observed counts
EXACTLY is a high-dimensional face

    F = { x : M x = y,  x_i >= x_{i+1},  x >= 0 },   y = M @ eta

(eta itself lies in F, so F is non-empty). The "best fit" is whichever single
point of F a tie-break objective selects. Production currently uses the
minimal column-norm flux (a vertex), which collapses weakly-constrained
high-v_min steps to zero (collapse) or piles flux onto one step (spike).

This script holds the face fixed and only swaps the tie-break objective, to
see which one tracks eta without collapse/spike:

  uniform   min  sum_j x_j                     (original minimal flux)
  colnorm   min  sum_j ||M[:,j]|| x_j          (production)
  smooth    min  sum_j (x_j - x_{j+1})^2       (L2 total variation)
  curv      min  sum_j (x_{j-1}-2x_j+x_{j+1})^2 (L2 second difference)
  maxent    max  sum_j entr(x_j)               (maximum entropy)

For each config we report the last-step ratio x[-1]/eta[-1] (1 = no
collapse/spike at the window edge) and the relative L2 error ||x-eta||/||eta||.
"""
import os
import sys

import numpy as np
import cvxpy as cp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, material, q, mass, nbins) -- 3 problem configs + 2 controls.
CONFIGS = [
    ("Al q2 M3  (collapse)", "Al", "2", "3", 5),
    ("TiN q2 M2 (collapse)", "TiN", "2", "2", 5),
    ("TiN q0 M3 (spike)",    "TiN", "0", "3", 5),
    ("Al q0 M1  (control)",  "Al", "0", "1", 5),
    ("TiN q0 M1 (control)",  "TiN", "0", "1", 5),
]

SOLVE = dict(solver=cp.CLARABEL, verbose=False)


def solve_face(M, y, eps, objective):
    """min objective(x) over F = {M x = y, x_i>=x_{i+1}, x>=eps}."""
    n = M.shape[1]
    x = cp.Variable(n)
    cons = [M @ x == y, x[:-1] - x[1:] >= 0, x >= eps]
    prob = cp.Problem(cp.Minimize(objective(x)), cons)
    prob.solve(**SOLVE)
    if x.value is None:
        return None
    return np.asarray(x.value)


def objectives(M):
    w = np.sqrt((M * M).sum(axis=0))  # column norms ||M[:,j]||

    def diff1(x):
        return cp.sum_squares(x[:-1] - x[1:])

    def diff2(x):
        return cp.sum_squares(x[:-2] - 2 * x[1:-1] + x[2:])

    return {
        "uniform": lambda x: cp.sum(x),
        "colnorm": lambda x: w @ x,
        "smooth":  diff1,
        "curv":    diff2,
        "maxent":  lambda x: -cp.sum(cp.entr(x)),
    }


def run():
    fig, axes = plt.subplots(len(CONFIGS), 1, figsize=(8, 3.0 * len(CONFIGS)))
    if len(CONFIGS) == 1:
        axes = [axes]

    hdr = f"{'config':>22} | " + " ".join(f"{k:>8}" for k in
                                          ["uniform", "colnorm", "smooth", "curv", "maxent"])
    print("\nLAST-STEP RATIO  x[-1]/eta[-1]  (1.00 = perfect, 0 = collapse, >>1 = spike)")
    print(hdr)
    print("-" * len(hdr))
    err_lines = []

    for ax, (label, mat, q, mass, nb) in zip(axes, CONFIGS):
        a = DarkMatterQuantumAnalysis(
            RunConfig(mat, q, mass, nb, "Halo", background="none"))
        M = a.m_phys
        eta = a.eta
        y = M @ eta
        v = a.vmin_mid
        eps = 1e-12 * float(eta.max())  # tiny positive floor (maxent needs x>0)

        objs = objectives(M)
        last = {}
        errs = {}
        ax.plot(v, eta, "k-", lw=2.0, label="eta (truth)", zorder=10)
        for name, obj in objs.items():
            x = solve_face(M, y, eps, obj)
            if x is None:
                last[name] = float("nan")
                errs[name] = float("nan")
                continue
            last[name] = x[-1] / eta[-1]
            errs[name] = np.linalg.norm(x - eta) / np.linalg.norm(eta)
            ax.plot(v, x, lw=1.0, alpha=0.85, label=name)

        print(f"{label:>22} | " +
              " ".join(f"{last[k]:8.2f}" for k in
                       ["uniform", "colnorm", "smooth", "curv", "maxent"]))
        err_lines.append(f"{label:>22} | " +
                         " ".join(f"{errs[k]:8.2e}" for k in
                                  ["uniform", "colnorm", "smooth", "curv", "maxent"]))

        ax.set_title(label)
        ax.set_xlabel("v_min [km/s]")
        ax.set_ylabel("flux (natural)")
        ax.set_yscale("log")
        ax.legend(fontsize=7, ncol=3)

    print("\nRELATIVE L2 ERROR  ||x-eta||/||eta||")
    print(hdr)
    print("-" * len(hdr))
    for ln in err_lines:
        print(ln)

    fig.tight_layout()
    out = os.path.join(HERE, "tiebreak_compare.pdf")
    fig.savefig(out)
    print(f"\nsaved {out}")

    edge_figure()


def edge_figure():
    """Edge-zoom: linear-y, last ~30% of v_min, only eta / current(colnorm) /
    curv / smooth so the collapse/spike vs recovery is unmistakable."""
    fig, axes = plt.subplots(len(CONFIGS), 1, figsize=(7, 2.6 * len(CONFIGS)))
    if len(CONFIGS) == 1:
        axes = [axes]
    show = {
        "colnorm (current)": ("colnorm", "tab:orange", "--"),
        "curv (proposed)":   ("curv",    "tab:red",    "-"),
        "smooth":            ("smooth",  "tab:green",  ":"),
    }
    for ax, (label, mat, q, mass, nb) in zip(axes, CONFIGS):
        a = DarkMatterQuantumAnalysis(
            RunConfig(mat, q, mass, nb, "Halo", background="none"))
        M, eta, v = a.m_phys, a.eta, a.vmin_mid
        y = M @ eta
        eps = 1e-12 * float(eta.max())
        objs = objectives(M)

        i0 = int(0.68 * len(v))  # last ~32%
        ax.plot(v[i0:], eta[i0:], "k-", lw=2.4, label="eta (truth)", zorder=10)
        for lbl, (name, col, ls) in show.items():
            x = solve_face(M, y, eps, objs[name])
            if x is not None:
                ax.plot(v[i0:], x[i0:], ls, color=col, lw=1.6, label=lbl)
        ax.set_title(label)
        ax.set_xlabel("v_min [km/s]")
        ax.set_ylabel("flux (natural)")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, "tiebreak_edge.pdf")
    fig.savefig(out)
    print(f"saved {out}")


if __name__ == "__main__":
    run()
