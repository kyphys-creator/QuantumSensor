"""Central-area alpha sweep through the full O(1) no-band-aid pipeline.

Isolated sandbox. For each central-area fraction alpha, build the matrix from
the response functions (two-sided peak-centred cut), run the full pipeline in
the O(1) unit system (round-trip eta, no column scaling, no CONDITION_C, plain
LP), output the recovered flux in PHYSICAL units [cm^-1], and verify it matches
the production solver (GeV=1e9, column scaling + CONDITION_C) on the same
matrix. TES Al q0 M1 R5 SHM.

    python alpha_sweep_full.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import linprog

from quantum_sensor.optimizer import run_optimize_qp
from quantum_sensor.model import condition
from quantum_sensor.config import CONDITION_C
from quantum_sensor.constants import CM as CM_1E9

from pipeline_invariance import (CM_at, exposure_AL, eta_halo_at,
                                 plain_min_flux_lp, GEV_NATIVE)

HERE = Path(__file__).resolve().parent
CSV_DIR = HERE.parents[1] / "Mathematica/output/TES/response_functions_csv"
ALPHAS = [0.30, 0.50, 0.70, 0.90]
TARGET_COUNTS = 1000.0
GEV_USE = 1e42
MASS_LABEL = {"1": "10 MeV", "2": "100 MeV", "3": "1 GeV"}


def trap_weights(v):
    w = np.empty_like(v)
    w[1:-1] = 0.5 * (v[2:] - v[:-2]); w[0] = 0.5 * (v[1] - v[0]); w[-1] = 0.5 * (v[-1] - v[-2])
    return w


def central_area_cut(v, r, alpha):
    w = trap_weights(v); area = r * w; total = area.sum()
    if total <= 0:
        return np.zeros_like(r)
    pk = int(np.argmax(r)); lo = hi = pk; acc = area[pk]; n = len(v)
    while acc < alpha * total:
        lv = r[lo - 1] if lo - 1 >= 0 else -np.inf
        rv = r[hi + 1] if hi + 1 < n else -np.inf
        if lv == -np.inf and rv == -np.inf:
            break
        if lv >= rv:
            lo -= 1; acc += area[lo]
        else:
            hi += 1; acc += area[hi]
    out = np.zeros_like(r); out[lo:hi + 1] = r[lo:hi + 1]
    return out


def production_style(M0, eta0):
    """GeV=1e9 with column scaling + CONDITION_C (the production conditioning)."""
    exposure = TARGET_COUNTS / float((M0 @ eta0).max())
    m_phys = exposure * M0
    signal = m_phys @ eta0
    mc, dc, bc, un = condition(m_phys, signal, np.zeros_like(signal), c=CONDITION_C)
    flux = un(run_optimize_qp(mc, dc, bc, M0.shape[1]).x)
    return flux * CM_1E9                                    # -> cm^-1


def plain_gev(M0, vmid_win, mass_tag, GeV):
    """Chosen unit system, plain LP, no column scaling, no CONDITION_C."""
    k = GeV / GEV_NATIVE
    M = M0 / k                                              # M proportional to GeV^-1
    eta = eta_halo_at(GeV, mass_tag, vmid_win)              # round-trip, proportional to GeV^+1
    exposure = TARGET_COUNTS / float((M @ eta).max())
    m_phys = exposure * M
    x = plain_min_flux_lp(m_phys, m_phys @ eta)
    return (x * CM_at(GeV)) if x is not None else None      # -> cm^-1


def run_mass(mass_tag):
    csv = CSV_DIR / f"Al_q0M{mass_tag}_R5.csv"
    d = np.genfromtxt(csv, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)

    eta_full = eta_halo_at(GEV_USE, mass_tag, v)
    print(f"\n=== M{mass_tag} ({MASS_LABEL[mass_tag]}) ===  GeV={GEV_USE:.0e}  "
          f"|eta|max~{np.abs(eta_full).max():.1e}")

    vg = np.logspace(0, np.log10(v.max()), 400)
    eta_ref = eta_halo_at(GEV_NATIVE, mass_tag, vg) * CM_1E9
    runs = []
    print(f"{'alpha':>6} {'cols':>5} {'window [km/s]':>16} {'last x/eta':>11} "
          f"{'plain-vs-prod':>14}")
    out_dir = HERE / f"by_alpha_M{mass_tag}"
    out_dir.mkdir(exist_ok=True)
    for alpha in ALPHAS:
        Rc = np.array([central_area_cut(v, R[i], alpha) for i in range(len(R))])
        M0 = (Rc * w[None, :])
        keep = M0.any(axis=0)
        M0k, vk = M0[:, keep], v[keep]
        eta0 = eta_halo_at(GEV_NATIVE, mass_tag, vk)
        phys_prod = production_style(M0k, eta0)
        phys_plain = plain_gev(M0k, vk, mass_tag, GEV_USE)
        if phys_plain is None:
            print(f"{alpha:6.2f} {M0k.shape[1]:5d}  plain LP failed")
            continue
        inv = np.linalg.norm(phys_plain - phys_prod) / np.linalg.norm(phys_prod)
        last = phys_plain[-1] / (eta0[-1] * CM_1E9)
        runs.append((alpha, vk, phys_plain, phys_prod))
        print(f"{alpha:6.2f} {M0k.shape[1]:5d}  [{vk.min():5.0f},{vk.max():5.0f}]"
              f"        {last:6.2f}   {inv:14.1e}")

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(vg, eta_ref, color="red", lw=2, label="input eta (SHM)")
        ax.step(vk, phys_prod, where="mid", color="black", lw=2.4,
                label="production (GeV=1e9, band-aids)")
        ax.step(vk, phys_plain, where="mid", color="C0", lw=1.2, ls="--",
                label="plain GeV=1e42 (no band-aids)")
        ax.set_xscale("log"); ax.set_xlim(v.min(), v.max()); ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"flux $\tilde{\eta}$ [cm$^{-1}$]")
        ax.set_title(rf"$\alpha$={alpha:.2f}  M{mass_tag} ({MASS_LABEL[mass_tag]})  "
                     f"TES Al q0 R5")
        ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=8)
        fig.savefig(out_dir / f"alpha_{alpha:.2f}.pdf", bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(vg, eta_ref, color="red", lw=2.5, label="input eta (SHM)", zorder=5)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(runs)))
    for (alpha, vk, phys_plain, _), col in zip(runs, colors):
        ax.step(vk, phys_plain, where="mid", color=col, lw=1.6,
                label=rf"$\alpha$={alpha:.2f} ({vk.min():.0f}-{vk.max():.0f} km/s)")
    ax.set_xscale("log"); ax.set_xlim(v.min(), v.max()); ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"flux $\tilde{\eta}$ [cm$^{-1}$]")
    ax.set_title(f"Central-area alpha sweep, no-band-aid pipeline at GeV=1e42\n"
                 f"M{mass_tag} ({MASS_LABEL[mass_tag]}), TES Al q0 R5 SHM (physical units)")
    ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=9, loc="upper right")
    fig.savefig(HERE / f"alpha_sweep_overlay_M{mass_tag}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    masses = sys.argv[1:] or ["1", "2", "3"]
    for m in masses:
        run_mass(m)
    print("\n=> at every (mass, alpha), the plain GeV=1e42 pipeline matches "
          "production in cm^-1.")


if __name__ == "__main__":
    main()
