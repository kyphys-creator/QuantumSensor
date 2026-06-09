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
    eta_tag = c.eta if c.disk_fraction is None else f"mix{round(c.disk_fraction * 100)}"
    return f"{c.material}_q{c.q}_M{c.mass}_R{c.nbins}_{eta_tag}"


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

    p = analysis.config.disk_fraction
    is_bound = analysis.config.eta == "Bound" and p is None
    if p is not None:
        eta_label = r"fit $\eta$ (mix " + f"{round(p * 100)}% disk)"
    elif is_bound:
        eta_label = r"fit $\eta$ (Bound + 100% SHM)"
    else:
        eta_label = r"input $\eta(v_{min})$"
    ax.plot(rm.vmin_mid, eta_phys, color="red", lw=2, label=eta_label)
    ax.hlines(flux_phys, rm.vmin_low, rm.vmin_high,
              color="C0", lw=1.5, label="Best-Fit")

    # Show the components that make up the fit eta:
    #   mixture -> 100% SHM (gray) and p x pure disk (green);
    #   Bound   -> 100% SHM (gray) and the bound population (green).
    if p is not None:
        if analysis.eta_halo is not None:
            ax.plot(rm.vmin_mid, analysis.eta_halo * ETA_TO_CM_INV,
                    color="gray", lw=1.5, ls="--", label="100% SHM")
            # the SHM component actually in the mixture, (1-p) * SHM
            ax.plot(rm.vmin_mid, (1.0 - p) * analysis.eta_halo * ETA_TO_CM_INV,
                    color="orange", lw=1.5, ls="-.",
                    label=f"{round((1 - p) * 100)}% SHM")
        if analysis.eta_disk is not None:
            ax.plot(rm.vmin_mid, p * analysis.eta_disk * ETA_TO_CM_INV,
                    color="green", lw=1.5, ls=":",
                    label=f"{round(p * 100)}% pure disk")
    elif is_bound:
        if analysis.eta_halo is not None:
            ax.plot(rm.vmin_mid, analysis.eta_halo * ETA_TO_CM_INV,
                    color="gray", lw=1.5, ls="--", label="100% SHM")
        if analysis.eta_bound is not None:
            ax.plot(rm.vmin_mid, analysis.eta_bound * ETA_TO_CM_INV,
                    color="green", lw=1.5, ls=":", label="Bound")

    ax.set_xscale("log")
    # Bound eta falls many decades over a few km/s, so show it on a log y-axis
    # (non-positive entries -- the zero tail -- are simply dropped by matplotlib).
    if analysis.config.eta == "Bound" and analysis.config.disk_fraction is None:
        ax.set_yscale("log")
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
