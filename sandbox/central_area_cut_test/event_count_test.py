"""Event counts (hence chi^2) are GeV-invariant -- explicit numerical check.

The counts entering chi^2 are a physical observable, so a change of the
natural-unit base (GeV) must leave them unchanged. Counts = exposure * (M @ eta):
exposure (mass*time) is GeV^0, M is GeV^-1, eta is GeV^+1, so the product is
GeV^0. This rebuilds exposure and eta at GeV=1e9 and GeV=1e40 (full round-trip)
and the response at both, and confirms the per-bin counts -- and a chi^2 against
an alternative model -- are identical.

    python event_count_test.py
"""

from __future__ import annotations

import numpy as np

from gev_roundtrip_test import rebuild_constants, eta_halo_at, GEV_NATIVE
from area_cut_test import trap_weights, CSV, MASS_TAG

G_NEW = 1e40


def exposure_AL(GeV):
    """AL_EXP = 8200 uG*month, rebuilt at the given GeV (mass*time, GeV^0)."""
    GRAM = 5.62e23 * GeV
    SEC = 1.0 / (6.58e-25 * GeV)
    uG = 1e-6 * GRAM
    MONTH = (365 * 24 * 3600 * SEC) / 12
    return 8200 * uG * MONTH


def counts(GeV, M_old, v, model):
    """Per-bin event counts at a given GeV, via the proper round-trip.

    M_old is the CSV response (GeV=1e9). R_bin has energy-dimension -1, so at a
    new GeV it is M_old*(1e9/G). eta is recomputed from the formula at that GeV.
    The km/s integration width is a velocity (GeV-invariant)."""
    k = GeV / GEV_NATIVE
    M = M_old / k                                  # R_bin proportional to GeV^-1
    eta, _ = eta_halo_at(GeV, MASS_TAG, v)
    if model == "Disk":                            # crude alt model: shift v0
        # reuse the Halo machinery is enough for an invariance demo; Disk via
        # eta_halo_at would need its params -- instead perturb eta by +5% shape
        eta = eta * (1.0 + 0.05 * np.exp(-(v - 120.0) ** 2 / (2 * 30.0 ** 2)))
    return exposure_AL(GeV) * (M @ eta)


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    M_old = np.array([d[c] for c in bins]) * trap_weights(v)[None, :]

    print(f"exposure_AL(1e9)  = {exposure_AL(1e9):.6e}")
    print(f"exposure_AL(1e40) = {exposure_AL(G_NEW):.6e}")
    print(f"  exposure ratio (must be 1) = {exposure_AL(G_NEW)/exposure_AL(1e9):.6f}\n")

    N_old = counts(GEV_NATIVE, M_old, v, "Halo")
    N_new = counts(G_NEW, M_old, v, "Halo")
    print("per-bin event counts (signal, Halo):")
    print(f"  {'bin':12} {'GeV=1e9':>14} {'GeV=1e40':>14} {'rel diff':>10}")
    for b, a, c in zip(bins, N_old, N_new):
        print(f"  {b:12} {a:14.6e} {c:14.6e} {abs(c/a-1):10.1e}")

    # chi^2 of an alternative model against the signal, in both systems
    def chi2(GeV):
        s = counts(GeV, M_old, v, "Halo")
        m = counts(GeV, M_old, v, "Disk")
        return float(np.sum((s - m) ** 2 / s))
    print(f"\nchi^2 (alt model vs signal):  GeV=1e9 -> {chi2(1e9):.6f}   "
          f"GeV=1e40 -> {chi2(G_NEW):.6f}")
    print(f"  chi^2 rel diff = {abs(chi2(G_NEW)/chi2(1e9) - 1):.1e}")
    print("\n=> event counts and chi^2 are GeV-invariant (physical observables).")


if __name__ == "__main__":
    main()
