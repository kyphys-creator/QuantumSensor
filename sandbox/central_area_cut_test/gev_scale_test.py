"""Does the natural-unit GeV magnitude (not column scaling) set the conditioning?

Isolated sandbox. Hypothesis (user): the best-fit needs the column-scaling
band-aid only because GeV = 1e9 puts eta ~ 1e-31, far below the solver's
absolute tolerances. Choosing GeV so eta ~ O(1) should let a PLAIN solve (no
column scaling, no CONDITION_C) recover the clean staircase by itself.

A change of the base unit GeV -> G is a similarity transform: the dimensionless
physics is unchanged, only the magnitudes move. Here
    eta  proportional to GeV^1   (RHO_DM*SIGMA_E/mchi = GeV^4 * GeV^-2 * GeV^-1;
                                   the velocity factors are GeV^0)
    M    proportional to GeV^-1  (so the physical rate M @ eta is GeV-invariant).
The Mathematica CSV is left untouched; we rescale its values in-memory.

We scan G, rebuild (eta_G, M_G) at that unit, fix the exposure so the counts
stay ~1000 (a GeV-invariant choice), and solve the minimal-total-flux monotone
fit with a PLAIN LP (scipy HiGHS, no column scaling). The production solver
(column scaling ON) is the GeV-invariant reference.

    python gev_scale_test.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linprog

from area_cut_test import central_area_cut, trap_weights, best_fit, CSV, MASS_TAG
from quantum_sensor.eta_models import eta as eta_model
from quantum_sensor.constants import DM_MASS

P_ETA = 1.0                  # eta proportional to GeV^P_ETA (derived above)
ALPHA = 0.50                 # tight central cut -> ideal answer is a clean staircase
TARGET_COUNTS = 1000.0
GEV_NATIVE = 1e9
HERE = Path(__file__).resolve().parent


def plain_min_flux_lp(M_phys, data):
    """Minimal-total-flux monotone exact-fit, PLAIN (no column scaling).

    minimize sum(x)  s.t.  M_phys x = data,  x_j >= x_{j+1},  x >= 0.
    Same estimator as the production vertex, but solved on the raw matrix so
    its numerical quality is exposed to the absolute magnitude of x."""
    n = M_phys.shape[1]
    # monotone: x_{j+1} - x_j <= 0
    rows = np.zeros((n - 1, n))
    for j in range(n - 1):
        rows[j, j + 1] = 1.0
        rows[j, j] = -1.0
    res = linprog(c=np.ones(n), A_ub=rows, b_ub=np.zeros(n - 1),
                  A_eq=M_phys, b_eq=data, bounds=(0, None), method="highs")
    return (res.x if res.success else None), res


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)
    eta0 = eta_model("Halo", DM_MASS[MASS_TAG], v)

    R_cut = np.array([central_area_cut(v, R[i], ALPHA)[0] for i in range(len(R))])
    M0 = (R_cut * w[None, :])
    keep = M0.any(axis=0)
    M0, vk, ek0 = M0[:, keep], v[keep], eta0[keep]
    # GeV-invariant exposure that puts the native counts at TARGET
    exposure = TARGET_COUNTS / float((M0 @ ek0).max())
    print(f"native: |eta| ~ {np.abs(ek0).max():.1e}  (GeV={GEV_NATIVE:.0e})  "
          f"window {vk.min():.0f}-{vk.max():.0f} km/s, {M0.shape[1]} cols\n")

    GEVS = 10.0 ** np.arange(9, 46, 3)
    rel_err, last_ratio, ok = [], [], []
    curves = {}
    for G in GEVS:
        k = (G / GEV_NATIVE) ** P_ETA
        ek = ek0 * k                      # eta proportional to GeV^1
        M = M0 / k                        # M proportional to GeV^-1 (rate invariant)
        m_phys = exposure * M
        data = m_phys @ ek                # counts: GeV-invariant (~TARGET)
        x, res = plain_min_flux_lp(m_phys, data)
        if x is None:
            rel_err.append(np.nan); last_ratio.append(np.nan); ok.append(False)
            print(f"GeV={G:.0e}  |eta|~{np.abs(ek).max():.0e}  LP FAILED ({res.message[:40]})")
        else:
            err = np.linalg.norm(x - ek) / np.linalg.norm(ek)
            lr = x[-1] / ek[-1]
            rel_err.append(err); last_ratio.append(lr); ok.append(True)
            curves[G] = x / ek            # recovered / truth, shape over window
            print(f"GeV={G:.0e}  |eta|~{np.abs(ek).max():.0e}  "
                  f"||x-eta||/||eta|| = {err:.2e}  last x/eta = {lr:.3f}")

    # production reference (column scaling ON) at native GeV
    xprod, _ = best_fit(M0, ek0, vertex_select=True)
    err_prod = np.linalg.norm(xprod - ek0) / np.linalg.norm(ek0)
    print(f"\nproduction (column scaling, native GeV): "
          f"||x-eta||/||eta|| = {err_prod:.2e}")

    # --- figure 1: plain-LP solve status + quality vs GeV ---
    GEVS = np.asarray(GEVS, float)
    err = np.asarray(rel_err, float)
    okm = np.asarray(ok, bool)
    SENT = 5.0                       # sentinel height for FAILED solves
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    if okm.any():
        ax.scatter(GEVS[okm], err[okm], s=70, color="C0", zorder=5,
                   label="plain LP solved (no column scaling)")
    if (~okm).any():
        ax.scatter(GEVS[~okm], np.full((~okm).sum(), SENT), s=90, marker="x",
                   color="C3", zorder=5, label="plain LP FAILED (model error)")
    ax.axhline(err_prod, color="C2", ls="--",
               label=f"production (column scaling) err = {err_prod:.2e}")
    ax.axvline(1e9, color="0.5", ls=":", label="native GeV = 1e9")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(err_prod * 0.5, SENT * 2)
    ax.set_xlabel("GeV (natural-unit base);  larger GeV -> larger |eta|")
    ax.set_ylabel(r"recovery error $\|x-\eta\|/\|\eta\|$  (failed = top)")
    ax.set_title("Plain best-fit needs the right GeV magnitude; once big enough\n"
                 "it equals the column-scaled result (TES Al q0 M1 R5, alpha=0.5)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=9, loc="center right")
    fig.savefig(HERE / "gev_scale_quality.pdf", bbox_inches="tight")
    print("saved gev_scale_quality.pdf")

    # --- figure 2: recovered/truth shape at a few GeV ---
    fig, ax = plt.subplots(figsize=(8, 5.5))
    show = [g for g in (1e9, 1e18, 1e30, 1e39, 1e45) if g in curves]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    for g, col in zip(show, colors):
        ax.plot(vk, curves[g], "o-", ms=3, color=col, label=f"GeV={g:.0e}")
    ax.axhline(1.0, color="red", lw=2, label=r"perfect ($x=\eta$)")
    ax.set_xscale("log")
    ax.set_ylim(-0.1, 2.0)
    ax.set_xlabel(r"$v_{min}$ [km/s]")
    ax.set_ylabel(r"recovered / truth  $x/\eta$")
    ax.set_title("Plain-LP recovery shape vs GeV (perfect = flat 1.0)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=9)
    fig.savefig(HERE / "gev_scale_shape.pdf", bbox_inches="tight")
    print("saved gev_scale_shape.pdf")


if __name__ == "__main__":
    main()
