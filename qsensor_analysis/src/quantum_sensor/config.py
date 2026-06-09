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
# Overall reparametrisation of the unknown, ``x = CONDITION_C * u`` (the DATA is
# never scaled, so the minimised objective stays the TRUE Neyman chi^2 and its
# value -- hence Delta-chi^2 intervals -- is preserved for statistics).
#
# The real conditioner is now the *per-column* scaling done inside the solver
# (``optimizer._column_scale``, the neutrinoAnalysis method): it normalises each
# design-matrix column to unit norm, which divides ``CONDITION_C`` out entirely.
# So the recovered flux is independent of this constant across any practical
# value (~1e-31 .. 1e-25); it survives only as a single, common, mass- and
# material-independent number (no per-quantity ``cons1``/``cons2``, no per-mass
# table). Keep it roughly O(1/||M column||) so OSQP's constraint matrix is well
# scaled.
CONDITION_C: float = 1e-27