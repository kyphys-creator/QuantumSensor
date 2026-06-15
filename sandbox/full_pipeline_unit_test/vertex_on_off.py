"""Vertex selection ON vs OFF, central-area sweep, in physical units.

Isolated sandbox. The minimal-total-flux vertex (ON) gives the staircase with
the last-step shave / collapse; the OSQP interior solution (OFF) is the smooth
point on the same chi2=0 face. This overlays the two for M1/M2/M3 at a couple
of alphas, in physical cm^-1. (The on/off difference is an estimator/tie-break
question, independent of the unit scale, so it is computed with the production
solver at its native GeV where it is well-conditioned.)

    python vertex_on_off.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quantum_sensor.optimizer import run_optimize_qp
from quantum_sensor.model import condition
from quantum_sensor.config import CONDITION_C
from quantum_sensor.constants import CM as CM_1E9

from pipeline_invariance import eta_halo_at, GEV_NATIVE
from alpha_sweep_full import (central_area_cut, trap_weights, CSV_DIR,
                              MASS_LABEL, TARGET_COUNTS)

HERE = Path(__file__).resolve().parent
ALPHAS = [0.30, 0.50, 0.70]


def solve(M0, eta0, vertex):
    """Physical recovered flux [cm^-1]; vertex=True staircase, False interior."""
    exposure = TARGET_COUNTS / float((M0 @ eta0).max())
    m_phys = exposure * M0
    signal = m_phys @ eta0
    mc, dc, bc, un = condition(m_phys, signal, np.zeros_like(signal), c=CONDITION_C)
    x = un(run_optimize_qp(mc, dc, bc, M0.shape[1], vertex_select=vertex).x)
    return x * CM_1E9


def run_mass(mass_tag):
    d = np.genfromtxt(CSV_DIR / f"Al_q0M{mass_tag}_R5.csv", delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)
    vg = np.logspace(0, np.log10(v.max()), 400)
    eta_ref = eta_halo_at(GEV_NATIVE, mass_tag, vg) * CM_1E9

    fig, axes = plt.subplots(1, len(ALPHAS), figsize=(6 * len(ALPHAS), 5), sharey=True)
    print(f"\n=== M{mass_tag} ({MASS_LABEL[mass_tag]}) ===  last x/eta:")
    for ax, alpha in zip(axes, ALPHAS):
        Rc = np.array([central_area_cut(v, R[i], alpha) for i in range(len(R))])
        M0 = (Rc * w[None, :])
        keep = M0.any(axis=0)
        M0k, vk = M0[:, keep], v[keep]
        eta0 = eta_halo_at(GEV_NATIVE, mass_tag, vk)
        eta0_phys = eta0 * CM_1E9
        x_on = solve(M0k, eta0, vertex=True)
        x_off = solve(M0k, eta0, vertex=False)
        print(f"  alpha={alpha:.2f}:  ON = {x_on[-1]/eta0_phys[-1]:.2f}   "
              f"OFF = {x_off[-1]/eta0_phys[-1]:.2f}")

        ax.plot(vg, eta_ref, color="red", lw=2, label="input eta (SHM)")
        ax.step(vk, x_on, where="mid", color="C0", lw=1.8,
                label="vertex ON (min-flux staircase)")
        ax.step(vk, x_off, where="mid", color="C1", lw=1.4, ls="--",
                label="vertex OFF (interior smooth)")
        ax.set_xscale("log"); ax.set_xlim(v.min(), v.max()); ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        ax.set_title(rf"$\alpha$={alpha:.2f}")
        ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=8)
    axes[0].set_ylabel(r"flux $\tilde{\eta}$ [cm$^{-1}$]")
    fig.suptitle(f"Vertex ON vs OFF  -  M{mass_tag} ({MASS_LABEL[mass_tag]}), "
                 f"TES Al q0 R5 SHM", fontsize=13)
    fig.tight_layout()
    fig.savefig(HERE / f"vertex_on_off_M{mass_tag}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    for m in ("1", "2", "3"):
        run_mass(m)
    print("\nsaved vertex_on_off_M{1,2,3}.pdf")


if __name__ == "__main__":
    main()
