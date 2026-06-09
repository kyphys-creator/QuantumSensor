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

from .constants import AL_EXP
from .data_loader import ResponseMatrix

# Exposure factor [natural units]. Multiplying the raw v_min-integrated response
# by the detector exposure gives expected event counts. Named/derived, not tuned.
EXPOSURE = {"Al": AL_EXP}


def exposure_factor(material: str) -> float:
    try:
        return EXPOSURE[material]
    except KeyError as exc:
        raise ValueError(f"no exposure defined for material {material!r}") from exc


def align_eta(eta_raw: np.ndarray, vmin_mid: np.ndarray,
              v_lo: float = 1.0, v_hi: float = 800.0) -> np.ndarray:
    """Resample a tabulated eta(v_min) onto the matrix's v_min interval mid-points.

    The legacy eta files hold ``len(eta_raw)`` samples spanning ``[v_lo, v_hi]``
    km/s; values outside the table extrapolate to its endpoints (eta -> 0 at the
    high end). The exact source grid only affects how the test flux is built,
    not the self-consistency of the inversion.
    """
    v_src = np.linspace(v_lo, v_hi, len(eta_raw))
    return np.interp(vmin_mid, v_src, eta_raw)


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
    """Precondition the linear inverse so the unknown is O(1) for the solver.

    Replaces the hand-tuned ``cons1``/``cons2``: both scales are derived from
    the data. We change variables ``x = flux_scale * u`` and divide the
    equations by ``data_scale`` so that, in the scaled problem
    ``M_cond @ u = data_cond - bkg_cond``, both ``M_cond`` entries and ``u``
    are O(1).

    * ``data_scale``  = median positive count  (rows -> O(1))
    * ``flux_scale``  = data_scale / median(M_phys @ 1)
                        (the flux that a uniform unit vector would need to
                        produce the data; sets the natural magnitude of x)

    Returns ``(M_cond, data_cond, bkg_cond, unscale)`` where
    ``unscale(u) = flux_scale*u`` recovers the physical flux.
    """
    positive = data[data > 0]
    data_scale = float(np.median(positive)) if positive.size else 1.0
    ref_counts = m_phys @ np.ones(m_phys.shape[1])
    ref = float(np.median(ref_counts[ref_counts > 0])) if np.any(ref_counts > 0) else 1.0
    flux_scale = data_scale / ref if ref > 0 else 1.0

    m_cond = (flux_scale / data_scale) * m_phys
    data_cond = data / data_scale
    if background is None:
        bkg_cond = np.zeros_like(data_cond)
    else:
        bkg_cond = background / data_scale

    def unscale(u: np.ndarray) -> np.ndarray:
        return flux_scale * np.asarray(u)

    return m_cond, data_cond, bkg_cond, unscale
