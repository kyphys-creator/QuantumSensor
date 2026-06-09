"""Plot and save analysis results.

The v_min axis and the staircase step edges come straight from the response
matrix's ``vmin.csv`` grid (via the analysis object), so there are no
hard-coded index ranges.

Units: the whole pipeline runs in natural units (eta and the recovered flux
are natural units, dimension 1/length). Conversion to the physical unit
cm^-1 happens ONLY here, at plot time, via ``ETA_TO_CM_INV`` -- the stored
CSV (``save_flux``) stays in natural units.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from .constants import CM
from .data_loader import DETECTOR_OF

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"

# eta / flux have natural dimension 1/length, so 1 in natural units equals
# ``CM`` cm^-1 (cm = CM in natural units). Multiply to display in cm^-1.
ETA_TO_CM_INV = CM


def _config_tag(analysis) -> str:
    """The per-run config tag (no background -- that is a parent folder)."""
    c = analysis.config
    return f"{c.material}_q{c.q}_M{c.mass}_R{c.nbins}_{c.eta}"


def _label(analysis) -> str:
    """Full human-readable label including the detector and background."""
    c = analysis.config
    det = DETECTOR_OF.get(c.material, c.material)
    return f"{det} {_config_tag(analysis)} bkg-{c.background}"


def run_dir(analysis) -> Path:
    """Output folder for one run: results/<DET>/bkg-<background>/<config tag>/.

    Each run gets its own folder so the flux CSV and the figure live together,
    grouped by detector (TES/MKID) and background scenario.
    """
    c = analysis.config
    det = DETECTOR_OF.get(c.material, c.material)
    return RESULTS_DIR / det / f"bkg-{c.background}" / _config_tag(analysis)


def save_flux(analysis, out_dir: Path | None = None) -> Path:
    """Save the recovered flux next to its v_min grid as ``flux.csv``."""
    if analysis.flux is None:
        raise RuntimeError("run optimize() first")
    out_dir = out_dir or run_dir(analysis)
    out_dir.mkdir(parents=True, exist_ok=True)
    rm = analysis.rm
    table = np.column_stack([rm.vmin_low, rm.vmin_high, rm.vmin_mid, analysis.flux])
    path = out_dir / "flux.csv"
    np.savetxt(path, table, delimiter=",",
               header="vmin_low,vmin_high,vmin_mid,flux", comments="")
    return path


def plot_flux_comparison(analysis, save: bool = True, ax=None, out_dir: Path | None = None):
    """Overlay the input eta(v_min) and the recovered staircase flux."""
    if analysis.flux is None:
        raise RuntimeError("run optimize() first")
    rm = analysis.rm
    own_fig = ax is None
    if own_fig:
        _, ax = plt.subplots(figsize=(8, 6))

    # natural units -> physical cm^-1, only for display
    eta_phys = analysis.eta * ETA_TO_CM_INV
    flux_phys = analysis.flux * ETA_TO_CM_INV

    ax.plot(rm.vmin_mid, eta_phys, color="red", lw=2, label=r"input $\eta(v_{min})$")
    ax.hlines(flux_phys, rm.vmin_low, rm.vmin_high,
              color="C0", lw=1.5, label="Best-Fit")

    ax.set_xscale("log")
    ax.set_xlim(rm.vmin_low[0], rm.vmin_high[-1])
    ax.set_xlabel(r"$v_{min}$ [km/s]", fontsize=18)
    ax.set_ylabel(r"$\tilde{\eta}$  [cm$^{-1}$]", fontsize=18)
    ax.set_title(_label(analysis), fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(True, which="both", ls="--", alpha=0.4)

    if save:
        out_dir = out_dir or run_dir(analysis)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "flux.pdf"
        plt.savefig(path, bbox_inches="tight")
        print(f"saved {path}")
    return ax
