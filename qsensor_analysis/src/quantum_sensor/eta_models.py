"""Analytic velocity-distribution integrals eta(v_min), in natural units.

A Python port of the Mathematica eta definitions (03_functions_response.wl's
etath and 12_eta.wl's bound model), used ONLY for plotting: it lets the
reference eta curves be drawn over the full v_min axis (1-800 km/s), whereas the
``eta_<model>.csv`` files written by 12_eta.wl are sampled only on each matrix's
narrower windowed grid. The fit itself still uses those CSV values; these
functions reproduce them (validated to ~1e-12 on the window) and simply extend
them to any v_min.

All velocities are passed in km/s; the returned eta is in the same natural units
as the CSV (GeV-based), so the caller multiplies by ETA_TO_CM_INV for display.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from .constants import GeV, CM, KPS, RHO_DM, SIGMA_E

# ----- velocity-model parameters (mirror Mathematica 01_setup.wl / 12_eta.wl) --
V0, VE, VESC = 238.0 * KPS, 250.0 * KPS, 544.0 * KPS          # standard halo (SHM)
V0_DD, VE_DD, VESC_DD = 70.0 * KPS, 100.0 * KPS, 694.0 * KPS  # pure dark disk
KEL = 8.62e-14 * GeV                                          # Kelvin in GeV
TEB, VESC_EB = 300.0 * KEL, 11.2 * KPS                        # Earth-bound
RHO_B = 1e14 * GeV / CM**3                                    # bound DM density


def _KKf(v0: float, vesc: float) -> float:
    """Truncated-Maxwellian normalisation (03_functions_response.wl: KKf)."""
    return v0**3 * (-2.0 * np.exp(-vesc**2 / v0**2) * np.pi * (vesc / v0)
                    + np.pi**1.5 * erf(vesc / v0))


def _eta_th(mchi: float, vm: np.ndarray, v0: float, ve: float, vesc: float) -> np.ndarray:
    """SHM-type halo speed integral etath[mchi][vm][v0, ve, vesc] (natural units)."""
    vm = np.asarray(vm, dtype=float)
    pref = (RHO_DM * SIGMA_E / mchi) * (v0**2 * np.pi / (2.0 * ve * _KKf(v0, vesc)))
    e_vesc = np.exp(-vesc**2 / v0**2)
    sqrtpi_v0 = np.sqrt(np.pi) * v0
    inner = np.where(
        vm < vesc - ve,
        -4.0 * e_vesc * ve + sqrtpi_v0 * (erf((vm + ve) / v0) - erf((vm - ve) / v0)),
        np.where(
            vm < vesc + ve,
            -2.0 * e_vesc * (ve + vesc - vm) + sqrtpi_v0 * (erf(vesc / v0) - erf((vm - ve) / v0)),
            0.0,
        ),
    )
    return pref * inner


def _eta_bound(mchi: float, vm: np.ndarray) -> np.ndarray:
    """Earth-bound isotropic thermal-DM speed integral (12_eta.wl: etabound)."""
    vm = np.asarray(vm, dtype=float)
    v0 = np.sqrt(2.0 * TEB / mchi)
    val = ((RHO_B * SIGMA_E / mchi) * (2.0 * np.pi * v0**2 / _KKf(v0, VESC_EB))
           * (np.exp(-vm**2 / v0**2) - np.exp(-VESC_EB**2 / v0**2)))
    return np.where(vm < VESC_EB, val, 0.0)


def eta(model: str, mchi: float, vmin_kms: np.ndarray) -> np.ndarray:
    """eta(v_min) for ``model`` in {"Halo", "Disk", "Bound"} at v_min [km/s]."""
    vm = np.asarray(vmin_kms, dtype=float) * KPS
    if model == "Halo":
        return _eta_th(mchi, vm, V0, VE, VESC)
    if model == "Disk":
        return _eta_th(mchi, vm, V0_DD, VE_DD, VESC_DD)
    if model == "Bound":
        return _eta_bound(mchi, vm)
    raise ValueError(f"unknown eta model {model!r}")
