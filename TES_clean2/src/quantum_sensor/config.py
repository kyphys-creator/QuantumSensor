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

    material : detector material ("Al").
    q        : FDM/mediator tag -- "0" heavy (n=0), "2" light (n=2).
    mass     : DM-mass tag -- "1" (10 MeV), "2" (100 MeV), "3" (1 GeV).
    nbins    : number of observed-energy bins (5 or 10).
    eta      : halo velocity-distribution model ("Halo", "Disk", "Bound").
    background : background scenario name (see BACKGROUND_SCENARIOS); "none"
                 fits the signal-only counts.
    run      : optional substring to pick a specific matrix folder/run.
    """

    material: str = "Al"
    q: str = "0"
    mass: str = "3"
    nbins: int = 5
    eta: str = "Halo"
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
