"""Sweep the central-area fraction alpha and watch the best-fit recovery.

Isolated sandbox. For each alpha, each bin's response keeps the highest-density
region around the peak covering `alpha` of the area (both tails cut); the matrix
is rebuilt and the production best-fit run (TES Al q0 M1 R5 SHM). Saves one
figure per alpha under by_alpha/, plus a single overlaid figure.

    python alpha_sweep.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from area_cut_test import (central_area_cut, trap_weights, best_fit,
                           CSV, MASS_TAG)
from quantum_sensor.eta_models import eta as eta_model
from quantum_sensor.constants import DM_MASS

ALPHAS = [0.30, 0.50, 0.70, 0.90, 0.95, 0.99]
HERE = Path(__file__).resolve().parent

# Vertex selection: True = production minimal-total-flux staircase; False =
# OSQP interior solution (no shaving of the lowest-weight high-v column).
VERTEX_SELECT = False
SUFFIX = "" if VERTEX_SELECT else "_novertex"
SUBDIR = "by_alpha" if VERTEX_SELECT else "by_alpha_novertex"


def cut_matrix(v, R, w, alpha):
    """Central-area-cut matrix and the kept-column mask."""
    R_cut = np.array([central_area_cut(v, R[i], alpha)[0] for i in range(len(R))])
    M = R_cut * w[None, :]
    keep = M.any(axis=0)
    return M[:, keep], keep


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)
    eta = eta_model("Halo", DM_MASS[MASS_TAG], v)

    vg = np.logspace(0, np.log10(v.max()), 400)
    eta_ref = eta_model("Halo", DM_MASS[MASS_TAG], vg)

    (HERE / SUBDIR).mkdir(exist_ok=True)
    print(f"vertex_select = {VERTEX_SELECT}\n")
    runs = []
    for alpha in ALPHAS:
        M, keep = cut_matrix(v, R, w, alpha)
        x, res = best_fit(M, eta[keep], vertex_select=VERTEX_SELECT)
        vk, ek = v[keep], eta[keep]
        # collapse metric: min of x/eta over the last 20% of the window
        tail = slice(int(0.8 * len(vk)), None)
        collapse = float(np.nanmin(np.where(ek[tail] > 0, x[tail] / ek[tail], np.nan)))
        runs.append((alpha, vk, x, ek, res))
        print(f"alpha={alpha:.2f}  {M.shape[1]:3d} cols  "
              f"window [{vk.min():5.0f},{vk.max():5.0f}] km/s "
              f"(width {vk.max()-vk.min():4.0f})  "
              f"tail x/eta min = {collapse:.2f}"
              f"{'  <- collapse' if collapse < 0.5 else ''}")

        # --- separate figure ---
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(vg, eta_ref, color="red", lw=2, label="input eta (SHM)")
        ax.step(vk, x, where="mid", color="C0", lw=1.5, label="best-fit")
        ax.axvspan(vk.min(), vk.max(), color="C0", alpha=0.08,
                   label=f"window {vk.min():.0f}-{vk.max():.0f} km/s")
        ax.set_xscale("log")
        ax.set_xlim(v.min(), v.max())
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        ax.set_ylabel(r"$\tilde{\eta}$ (natural units)")
        vs_lbl = "vertex" if VERTEX_SELECT else "interior (no vertex)"
        ax.set_title(rf"central-area cut $\alpha$ = {alpha:.2f}  "
                     f"(TES Al q0 M1 R5, SHM, {vs_lbl})")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=9)
        fig.savefig(HERE / SUBDIR / f"alpha_{alpha:.2f}.pdf", bbox_inches="tight")
        plt.close(fig)

    # --- overlaid figure ---
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(vg, eta_ref, color="red", lw=2.5, label="input eta (SHM)", zorder=5)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(runs)))
    for (alpha, vk, x, ek, res), col in zip(runs, colors):
        ax.step(vk, x, where="mid", color=col, lw=1.5,
                label=rf"$\alpha$={alpha:.2f}  ({vk.min():.0f}-{vk.max():.0f} km/s)")
    ax.set_xscale("log")
    ax.set_xlim(v.min(), v.max())
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$v_{min}$ [km/s]")
    ax.set_ylabel(r"$\tilde{\eta}$ (natural units)")
    vs_lbl = "vertex" if VERTEX_SELECT else "interior / no vertex"
    ax.set_title("Best-fit recovery vs central-area fraction "
                 rf"$\alpha$ (TES Al q0 M1 R5, SHM, {vs_lbl})")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=9, loc="upper right")
    fig.savefig(HERE / f"alpha_sweep_overlay{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"\nsaved {SUBDIR}/alpha_*.pdf  and  alpha_sweep_overlay{SUFFIX}.pdf")


if __name__ == "__main__":
    main()
