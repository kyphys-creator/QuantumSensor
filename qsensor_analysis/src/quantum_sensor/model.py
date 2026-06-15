"""Forward model: turn a response matrix + velocity distribution into expected
event counts, and assemble the linear operator the optimiser inverts.

Design note (no magic numbers):

* The exposure/unit scale ``K`` is derived from named constants (the Al
  exposure), not a tuned literal. It only sets the absolute count level
  (relevant for Poisson noise and the background comparison); the recovered
  flux *shape* is independent of ``K`` because it cancels in
  ``argmin_x |K M eta - K M x|`` = eta.
* The numerical conditioning scale (the old ``cons2``) is computed from the
  data itself, not hard-coded.
* The v_min grid and energy bins come from the loaded matrix, never from
  literals here.
"""

from __future__ import annotations

import numpy as np

from .constants import AL_EXP, TIN_EXP
from .data_loader import ResponseMatrix

# Exposure factor [natural units]. Multiplying the raw v_min-integrated response
# by the detector exposure gives expected event counts. Named/derived, not tuned.
EXPOSURE = {"Al": AL_EXP, "TiN": TIN_EXP}


def exposure_factor(material: str) -> float:
    try:
        return EXPOSURE[material]
    except KeyError as exc:
        raise ValueError(f"no exposure defined for material {material!r}") from exc


def response_operator(rm: ResponseMatrix, material: str) -> np.ndarray:
    """The physical forward operator M_phys = exposure * raw matrix.

    ``M_phys @ eta`` gives expected counts per energy bin.
    """
    return exposure_factor(material) * rm.matrix


def expected_counts(rm: ResponseMatrix, eta_aligned: np.ndarray,
                    material: str) -> np.ndarray:
    """Forward model: expected event counts per energy bin for flux ``eta``."""
    return response_operator(rm, material) @ eta_aligned


def condition(m_phys: np.ndarray, data: np.ndarray,
              background: np.ndarray | None = None):
    """Prepare the solver arrays: float data and a zero-filled background.

    Formerly this also applied an overall ``x = c * u`` reparametrisation (the
    ``config.CONDITION_C`` constant) to bring the ~1e-31 unknown up to O(1) for
    the solver. With ``constants.GeV`` raised so the natural-unit eta is already
    O(1e2), that reparametrisation is unnecessary and has been removed; column
    scaling inside the solver (:func:`optimizer._column_scale`) is the remaining,
    scale-adaptive conditioner. The matrix and data pass through untouched, so
    the minimised objective is the true Neyman chi^2.

    Returns ``(m_phys, data_cond, bkg_cond)``.
    """
    data_cond = np.asarray(data, dtype=float)
    if background is None:
        bkg_cond = np.zeros_like(data_cond)
    else:
        bkg_cond = np.asarray(background, dtype=float)
    return m_phys, data_cond, bkg_cond
