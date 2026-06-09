"""Background spectra, integrated onto the analysis' energy bins.

The bin edges come from the loaded response matrix (``bins.csv``), not a
hard-coded ``BINS`` array. The scenario amplitudes live in
``config.BACKGROUND_SCENARIOS``.
"""

from __future__ import annotations

import numpy as np

from .config import BackgroundModel, BACKGROUND_SCENARIOS
from .constants import eV, keV

# The legacy background amplitudes are quoted per keV per day; this is the same
# unit factor the original code applied (eV/keV * days-per-year).
DAYS_PER_YEAR = 365.0
AMP_SCALE = (eV / keV) * DAYS_PER_YEAR


def _integrate(model: BackgroundModel, e_lo: float, e_hi: float) -> float:
    """Integral of A*exp(-E/B) + C over [e_lo, e_hi] (eV), with AMP_SCALE."""
    expo = -model.A * model.B * (np.exp(-e_hi / model.B) - np.exp(-e_lo / model.B))
    const = model.C * (e_hi - e_lo)
    return AMP_SCALE * (expo + const)


def background_counts(scenario: str, ebin_low: np.ndarray,
                      ebin_high: np.ndarray) -> np.ndarray:
    """Per-energy-bin background counts for a named scenario.

    Returns zeros for the 'none'/'a' (null) scenarios. The absolute background
    normalisation relative to the new-pipeline signal may need recalibration;
    the shape and scenario ratios follow the legacy model.
    """
    if scenario not in BACKGROUND_SCENARIOS:
        raise ValueError(
            f"unknown background scenario {scenario!r}; "
            f"choose from {sorted(BACKGROUND_SCENARIOS)}")
    model = BACKGROUND_SCENARIOS[scenario]
    return np.array([_integrate(model, lo, hi)
                     for lo, hi in zip(ebin_low, ebin_high)])
