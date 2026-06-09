"""Top-level pipeline: dark-matter quantum-sensor response analysis.

Self-consistent forward + inverse model built on the new Mathematica response
matrices:

    1. load the response matrix M and its v_min / energy grids,
    2. align the chosen halo velocity distribution eta onto the matrix's v_min
       grid,
    3. forward model  observed = exposure * M @ eta  (+ background),
    4. invert: recover a monotone, non-negative flux x that reproduces the
       observed energy-bin counts (the under-determined inverse, regularised by
       the staircase/vertex selection).

No tuned scale factors: the exposure comes from constants, the numerical
conditioning is derived from the data, and the grids come from the data files.
"""

from __future__ import annotations

import numpy as np

from .config import RunConfig
from .data_loader import load_response_matrix, load_eta
from .model import align_eta, response_operator, condition
from .backgrounds import background_counts
from .optimizer import run_optimize, run_optimize_qp


class DarkMatterQuantumAnalysis:
    """Run the forward+inverse analysis for one :class:`RunConfig`."""

    def __init__(self, config: RunConfig):
        self.config = config

        # 1. response matrix + its grids (shapes/bins/v_min all come from here)
        self.rm = load_response_matrix(
            config.material, config.q, config.mass, config.nbins, run=config.run)

        # 2. velocity distribution aligned to the matrix v_min grid
        self.eta = align_eta(load_eta(config.eta, config.mass), self.rm.vmin_mid)

        # 3. forward model -> observed counts
        self.m_phys = response_operator(self.rm, config.material)
        self.signal = self.m_phys @ self.eta
        self.background = background_counts(
            config.background, self.rm.ebin_low, self.rm.ebin_high)
        self.observed = self.signal + self.background

        self.result = None
        self.flux = None  # physical recovered flux x(v_min)

    # -- convenience read-only grids ----------------------------------------
    @property
    def vmin_mid(self) -> np.ndarray:
        return self.rm.vmin_mid

    @property
    def n_vmin(self) -> int:
        return self.rm.n_vmin

    @property
    def n_ebins(self) -> int:
        return self.rm.n_ebins

    # -- inversion -----------------------------------------------------------
    def optimize(self, solver: str = "osqp", x0=None, display: bool = False,
                 fix=None, vertex_select: bool = True):
        """Recover the monotone non-negative flux from the observed counts."""
        m_cond, data_cond, bkg_cond, unscale = condition(
            self.m_phys, self.observed, self.background)

        if solver in ("osqp", "qp", "clarabel") or fix is not None:
            res = run_optimize_qp(
                m_cond, data_cond, bkg_cond, self.n_vmin,
                eps=0.0, fix=fix, vertex_select=vertex_select, verbose=display)
        else:
            res = run_optimize(
                m_cond, data_cond, bkg_cond, self.n_vmin,
                eps=0.0, x0=x0, display=display)

        self.result = res
        self.flux = unscale(res.x)        # back to physical units
        return self.flux

    # -- output (matplotlib imported lazily so it stays an optional dep) ------
    def plot(self, save: bool = True, ax=None):
        from .plotting import plot_flux_comparison
        return plot_flux_comparison(self, save=save, ax=ax)

    def save_flux(self, out_dir=None):
        from .plotting import save_flux
        return save_flux(self, out_dir=out_dir)
