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

    LEGACY / cross-check only. The main pipeline no longer uses this: it loads a
    natural-units eta already sampled on the matrix grid (``load_natural_eta`` /
    ``12_eta.wl``). This resamples the legacy *physical-units* eta files, which
    are NOT unit-consistent with the natural-units matrix.

    The legacy eta files hold ``len(eta_raw)`` samples spanning ``[v_lo, v_hi]``
    km/s; values outside the table extrapolate to its endpoints (eta -> 0 at the
    high end).
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
              background: np.ndarray | None = None,
              c: float = 1.0):
    """Recondition the linear inverse with a SINGLE common constant ``c``.

    The unknown is reparametrised ``x = c * u`` (column scaling only); the data
    and background are left **untouched**. Hence the objective the optimiser
    minimises is the *true* Neyman chi^2 -- its value, and so Delta-chi^2
    confidence intervals, are preserved. (A data-scaling conditioner, like the
    old ``cons2``, would multiply chi^2 by an arbitrary factor and break the
    statistical interpretation.)

    ``c`` is one common numerical knob (``config.CONDITION_C``), the same for
    all masses/materials. It only shifts which minimiser OSQP lands on for the
    underdetermined chi^2=0 face; it cancels out of the recovered flux
    (``unscale`` undoes it), so it never changes chi^2 or the physical answer.

    Returns ``(M_cond, data_cond, bkg_cond, unscale)`` with
    ``M_cond = c*m_phys``, ``data_cond = data``, ``bkg_cond = background`` and
    ``unscale(u) = c*u``.
    """
    data_cond = np.asarray(data, dtype=float)
    m_cond = c * m_phys
    if background is None:
        bkg_cond = np.zeros_like(data_cond)
    else:
        bkg_cond = np.asarray(background, dtype=float)

    def unscale(u: np.ndarray) -> np.ndarray:
        return c * np.asarray(u)

    return m_cond, data_cond, bkg_cond, unscale
