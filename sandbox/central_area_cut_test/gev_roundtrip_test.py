"""Is the GeV change a naive rescale, or must we go through physical units?

User's concern: the Mathematica numbers were all computed in the GeV=1e9 natural
system. To "switch to GeV=1e40" one should arguably convert each quantity back
to PHYSICAL units and then re-express it in the new natural system -- not just
multiply by a guessed power of (G/1e9).

This test does the proper round-trip explicitly: it rebuilds ALL constants at a
new GeV (mirroring constants.py / eta_models.py) and recomputes eta from the
formula at that GeV. It then checks this against

  (a) the shortcut  eta_new = eta_old * (G/1e9)^1     used in gev_scale_test, and
  (b) the physical invariant  eta * CM  [cm^-1], which must be identical in
      every unit system (that IS the round-trip: ->physical->new natural).

Isolated sandbox; nothing here is imported by production.

    python gev_roundtrip_test.py
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from quantum_sensor.eta_models import eta as eta_model_native
from quantum_sensor.constants import DM_MASS as DM_MASS_NATIVE, CM as CM_NATIVE

GEV_NATIVE = 1e9


def rebuild_constants(GeV):
    """Recompute the constants needed for the SHM eta from a GeV value,
    verbatim from constants.py (so this IS the physical<->natural conversion)."""
    eV = 1e-9 * GeV
    MeV = 1e-3 * GeV
    CM = 1.0 / (1.98e-14 * GeV)
    SEC = 1.0 / (6.58e-25 * GeV)
    KM = 1e5 * CM
    KPS = KM / SEC
    RHO_DM = 0.4 * GeV / CM**3
    SIGMA_E = 1e-30 * CM**2
    return dict(GeV=GeV, MeV=MeV, CM=CM, KPS=KPS, RHO_DM=RHO_DM, SIGMA_E=SIGMA_E)


def _KKf(v0, vesc):
    return v0**3 * (-2.0 * np.exp(-vesc**2 / v0**2) * np.pi * (vesc / v0)
                    + np.pi**1.5 * erf(vesc / v0))


def eta_halo_at(GeV, mass_tag, vmin_kms):
    """SHM eta(v_min) recomputed entirely at the given GeV (full round-trip)."""
    c = rebuild_constants(GeV)
    mchi = {"1": 10, "2": 100, "3": 1000}[mass_tag] * c["MeV"]
    V0, VE, VESC = 238.0 * c["KPS"], 250.0 * c["KPS"], 544.0 * c["KPS"]
    vm = np.asarray(vmin_kms, float) * c["KPS"]
    pref = (c["RHO_DM"] * c["SIGMA_E"] / mchi) * (V0**2 * np.pi / (2.0 * VE * _KKf(V0, VESC)))
    e_vesc = np.exp(-VESC**2 / V0**2)
    sp = np.sqrt(np.pi) * V0
    inner = np.where(
        vm < VESC - VE,
        -4.0 * e_vesc * VE + sp * (erf((vm + VE) / V0) - erf((vm - VE) / V0)),
        np.where(vm < VESC + VE,
                 -2.0 * e_vesc * (VE + VESC - vm) + sp * (erf(VESC / V0) - erf((vm - VE) / V0)),
                 0.0))
    return pref * inner, c["CM"]


def main():
    vk = np.logspace(np.log10(50), np.log10(300), 12)   # a few v_min points
    G_NEW = 1e40

    # native (GeV=1e9), three ways:
    eta_native_model = eta_model_native("Halo", DM_MASS_NATIVE["1"], vk)   # production
    eta_native_rb, CM_old = eta_halo_at(GEV_NATIVE, "1", vk)               # our rebuild
    print("rebuild matches production eta_model at GeV=1e9:  "
          f"max rel diff = {np.abs(eta_native_rb/eta_native_model - 1).max():.2e}\n")

    # new system (GeV=1e40): full round-trip recompute vs the shortcut
    eta_new_roundtrip, CM_new = eta_halo_at(G_NEW, "1", vk)
    k = (G_NEW / GEV_NATIVE) ** 1                          # eta proportional to GeV^1
    eta_new_shortcut = eta_native_rb * k

    print(f"switch GeV: {GEV_NATIVE:.0e} -> {G_NEW:.0e}   (k = G/1e9 = {k:.0e})")
    print(f"  |eta| native     ~ {np.abs(eta_native_rb).max():.3e}")
    print(f"  |eta| roundtrip   ~ {np.abs(eta_new_roundtrip).max():.3e}")
    print(f"  |eta| shortcut    ~ {np.abs(eta_new_shortcut).max():.3e}")
    print(f"  roundtrip vs shortcut: max rel diff = "
          f"{np.abs(eta_new_roundtrip/eta_new_shortcut - 1).max():.2e}")

    # the physical observable (eta in cm^-1) must be identical in both systems
    phys_old = eta_native_rb * CM_old
    phys_new = eta_new_roundtrip * CM_new
    print(f"\nphysical eta [cm^-1] invariance (the actual round-trip ->phys->):")
    print(f"  CM(1e9) = {CM_old:.3e}   CM(1e40) = {CM_new:.3e}   ratio = {CM_old/CM_new:.3e}")
    print(f"  max rel diff between eta*CM in the two systems = "
          f"{np.abs(phys_new/phys_old - 1).max():.2e}")

    print("\nconclusion:")
    print("  - full round-trip recompute == shortcut eta*(G/1e9)^1  (to machine eps)")
    print("  - because eta = [RHO_DM*SIGMA_E/mchi] (GeV^1) x [velocity ratios] (GeV^0),")
    print("    so the GeV dependence factorises cleanly; the velocities (KPS) are")
    print("    GeV-invariant. The shortcut used the correct energy-dimension (+1),")
    print("    hence it already WAS the ->physical->new-natural round-trip.")


if __name__ == "__main__":
    main()
