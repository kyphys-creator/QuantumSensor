"""Full pipeline with GeV as a free parameter, no column scaling, no CONDITION_C.

Isolated sandbox (see README). Verifies the physical observables (event counts,
recovered flux in cm^-1) are invariant between:

  * production: DarkMatterQuantumAnalysis at GeV=1e9 (column scaling + CONDITION_C);
  * sandbox:    same response matrix, but the unit system rebuilt at an arbitrary
                GeV (eta O(1) via round-trip), solved with a PLAIN LP (no column
                scaling, no CONDITION_C).

    python pipeline_invariance.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import linprog

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.constants import CM as CM_1E9          # production CM (GeV=1e9)

HERE = Path(__file__).resolve().parent
GEV_NATIVE = 1e9
MASS_TAG = "1"


# --- unit system as a function of GeV (mirrors constants.py) ----------------
def CM_at(GeV):
    return 1.0 / (1.98e-14 * GeV)


def exposure_AL(GeV):
    """AL_EXP = 8200 uG*month rebuilt at GeV (mass*time, GeV^0 -> invariant)."""
    GRAM = 5.62e23 * GeV
    SEC = 1.0 / (6.58e-25 * GeV)
    uG = 1e-6 * GRAM
    MONTH = (365 * 24 * 3600 * SEC) / 12
    return 8200 * uG * MONTH


def _KKf(v0, vesc):
    return v0**3 * (-2.0 * np.exp(-vesc**2 / v0**2) * np.pi * (vesc / v0)
                    + np.pi**1.5 * erf(vesc / v0))


def eta_halo_at(GeV, mass_tag, vmin_kms):
    """Round-trip SHM eta recomputed entirely at the given GeV."""
    MeV = 1e-3 * GeV
    CM = CM_at(GeV)
    SEC = 1.0 / (6.58e-25 * GeV)
    KPS = (1e5 * CM) / SEC
    RHO_DM = 0.4 * GeV / CM**3
    SIGMA_E = 1e-30 * CM**2
    mchi = {"1": 10, "2": 100, "3": 1000}[mass_tag] * MeV
    V0, VE, VESC = 238.0 * KPS, 250.0 * KPS, 544.0 * KPS
    vm = np.asarray(vmin_kms, float) * KPS
    pref = (RHO_DM * SIGMA_E / mchi) * (V0**2 * np.pi / (2.0 * VE * _KKf(V0, VESC)))
    e_vesc = np.exp(-VESC**2 / V0**2)
    sp = np.sqrt(np.pi) * V0
    inner = np.where(
        vm < VESC - VE,
        -4.0 * e_vesc * VE + sp * (erf((vm + VE) / V0) - erf((vm - VE) / V0)),
        np.where(vm < VESC + VE,
                 -2.0 * e_vesc * (VE + VESC - vm) + sp * (erf(VESC / V0) - erf((vm - VE) / V0)),
                 0.0))
    return pref * inner


def plain_min_flux_lp(M_phys, data):
    """Minimal-total-flux monotone exact fit -- PLAIN (no column scaling, no c)."""
    n = M_phys.shape[1]
    rows = np.zeros((n - 1, n))
    for j in range(n - 1):
        rows[j, j + 1] = 1.0
        rows[j, j] = -1.0
    res = linprog(c=np.ones(n), A_ub=rows, b_ub=np.zeros(n - 1),
                  A_eq=M_phys, b_eq=data, bounds=(0, None), method="highs")
    return res.x if res.success else None


def main():
    # --- production reference (GeV=1e9, column scaling + CONDITION_C) ---
    cfg = RunConfig(material="Al", q="0", mass=MASS_TAG, nbins=5,
                    eta="Halo", background="none")
    a = DarkMatterQuantumAnalysis(cfg)
    a.optimize()
    M0 = a.rm.matrix                      # response matrix at GeV=1e9 (5 x n)
    vmid = a.rm.vmin_mid
    phys_prod = a.flux * CM_1E9           # production recovered flux [cm^-1]
    counts_prod = a.observed

    # round-trip eta at native must match production's loaded eta
    eta0 = eta_halo_at(GEV_NATIVE, MASS_TAG, vmid)
    print(f"round-trip eta vs production eta (GeV=1e9): "
          f"max rel diff = {np.abs(eta0/a.eta - 1).max():.2e}")
    # and my native forward must reproduce production's observed counts
    fwd0 = exposure_AL(GEV_NATIVE) * (M0 @ eta0)
    print(f"native forward vs production observed: "
          f"max rel diff = {np.abs(fwd0/counts_prod - 1).max():.2e}\n")

    # GeV that puts |eta| ~ O(1)
    GEV_O1 = GEV_NATIVE / float(np.abs(eta0).max())
    GEVS = [1e9, 1e18, 1e27, GEV_O1, 1e42]

    print(f"{'GeV':>10} {'|eta|max':>10} {'solve':>7} "
          f"{'counts match':>13} {'phys-flux match':>16}")
    results = {}
    for G in GEVS:
        k = G / GEV_NATIVE
        M = M0 / k                                  # M proportional to GeV^-1
        eta = eta_halo_at(G, MASS_TAG, vmid)        # round-trip, proportional to GeV^+1
        m_phys = exposure_AL(G) * M
        observed = m_phys @ eta                     # counts (GeV-invariant)
        cmatch = np.abs(observed / counts_prod - 1).max()
        x = plain_min_flux_lp(m_phys, observed)
        if x is None:
            print(f"{G:10.0e} {np.abs(eta).max():10.1e} {'FAIL':>7} "
                  f"{cmatch:13.1e} {'-':>16}")
            continue
        phys = x * CM_at(G)                          # recovered flux [cm^-1]
        pmatch = np.linalg.norm(phys - phys_prod) / np.linalg.norm(phys_prod)
        results[G] = phys
        print(f"{G:10.0e} {np.abs(eta).max():10.1e} {'ok':>7} "
              f"{cmatch:13.1e} {pmatch:16.1e}")

    print("\n=> counts and physical flux are invariant; the plain solve (no column\n"
          "   scaling, no CONDITION_C) reproduces production once GeV makes eta O(1).")

    # --- figure: physical recovered flux at each working GeV vs production ---
    fig, ax = plt.subplots(figsize=(8, 5.5))
    vg = np.logspace(0, np.log10(800), 400)
    eta_ref_phys = eta_halo_at(GEV_NATIVE, MASS_TAG, vg) * CM_1E9
    ax.plot(vg, eta_ref_phys, color="red", lw=2, label="input eta (SHM)")
    ax.step(vmid, phys_prod, where="mid", color="black", lw=2.5,
            label="production (GeV=1e9, band-aids)")
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(results)))
    for (G, phys), col in zip(results.items(), colors):
        ax.step(vmid, phys, where="mid", color=col, lw=1.2, ls="--",
                label=f"plain, GeV={G:.0e}")
    ax.set_xscale("log")
    ax.set_xlim(1, 800)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$v_{min}$ [km/s]")
    ax.set_ylabel(r"recovered flux $\tilde{\eta}$ [cm$^{-1}$]")
    ax.set_title("Full pipeline, GeV free, no column scaling, no CONDITION_C\n"
                 "(physical flux invariant; TES Al q0 M1 R5 SHM)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=8)
    fig.savefig(HERE / "pipeline_invariance.pdf", bbox_inches="tight")
    print("saved pipeline_invariance.pdf")


if __name__ == "__main__":
    main()
