"""Typed configuration objects.

All run parameters and the background-model numbers live here as named
dataclass fields / a single scenario table -- never as bare literals scattered
through the code.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RunConfig:
    """One analysis configuration.

    material : detector material -- "Al" (TES) or "TiN" (MKID). Selects which
               Mathematica output tree the response matrix is read from
               (see data_loader.DETECTOR_OF) and which exposure is applied.
    q        : FDM/mediator tag -- "0" heavy (n=0), "2" light (n=2).
    mass     : DM-mass tag -- "1" (10 MeV), "2" (100 MeV), "3" (1 GeV).
    nbins    : number of observed-energy bins. Al (TES): 5 or 10. TiN (MKID):
               5 or 9 (threshold 0.2 eV, so the fine binning is R9, not R10).
    eta      : velocity-distribution model used when disk_fraction is None --
               "Halo" (SHM) or "Disk" (pure dark disk).
    disk_fraction : if set (0..1), the eta driving the fit is the MIXTURE
               (1 - p) * eta_Halo + p * eta_Disk with p = disk_fraction, instead
               of the single `eta` model. Both components are normalised to the
               same local DM density, so the mixture's density equals 100% SHM
               for any p. None (default) = use the single `eta` model.
    background : background scenario name (see BACKGROUND_SCENARIOS); "none"
                 fits the signal-only counts.
    run      : optional substring to pick a specific matrix folder/run.
    """

    material: str = "Al"
    q: str = "0"
    mass: str = "3"
    nbins: int = 5
    eta: str = "Halo"
    disk_fraction: float | None = None
    background: str = "none"
    run: str | None = None


@dataclass(frozen=True)
class BackgroundModel:
    """Background spectrum R(E) = A*exp(-E/B) + C, integrated over each E bin."""

    A: float
    B: float
    C: float


# The only place the background numbers live. Energy in eV, B in eV.
BACKGROUND_SCENARIOS: dict[str, BackgroundModel] = {
    "none": BackgroundModel(A=0.0, B=10.0, C=0.0),    # signal only
    "a":    BackgroundModel(A=0.0, B=10.0, C=0.0),    # null
    "c":    BackgroundModel(A=50.0, B=10.0, C=20.0),  # medium
    "b":    BackgroundModel(A=460.0, B=10.0, C=100.0),
    "b2":   BackgroundModel(A=920.0, B=10.0, C=200.0),
    "flat": BackgroundModel(A=0.0, B=10.0, C=100.0),
}


# --- Solver conditioning constant (single, common) --------------------------
# Overall reparametrisation of the unknown, ``x = CONDITION_C * u``. It exists
# only to bring the unknown to O(1) for the solver. With the unit base raised
# (``constants.GeV = 1e42``) the natural-unit eta is already O(1e2), so the
# reparametrisation is the identity: CONDITION_C = 1. (At the old GeV=1e9,
# eta ~ 1e-31 forced CONDITION_C ~ 1e-30; the constant is scale-specific, which
# is why it moves with GeV.) This will be removed entirely once column scaling
# is dropped (B4); kept as 1.0 here so the change is a pure no-op.
CONDITION_C: float = 1.0