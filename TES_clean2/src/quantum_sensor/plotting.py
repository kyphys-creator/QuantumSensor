"""Plot and save analysis results.

The v_min axis and the staircase step edges come straight from the response
matrix's ``vmin.csv`` grid (via the analysis object), so there are no
hard-coded index ranges. The recovered flux is already in the same units as
the input eta, so no extra unit factor is applied.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"


def _stem(analysis) -> str:
    c = analysis.config
    return f"{c.material}_q{c.q}_M{c.mass}_R{c.nbins}_{c.eta}_bkg-{c.background}"


def save_flux(analysis, out_dir: Path | None = None) -> Path:
    """Save the recovered flux next to its v_min grid as CSV."""
    if analysis.flux is None:
        raise RuntimeError("run optimize() first")
    out_dir = out_dir or RESULTS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    rm = analysis.rm
    table = np.column_stack([rm.vmin_low, rm.vmin_high, rm.vmin_mid, analysis.flux])
    path = out_dir / f"flux_{_stem(analysis)}.csv"
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

    ax.plot(rm.vmin_mid, analysis.eta, color="red", lw=2, label=r"input $\eta(v_{min})$")
    ax.hlines(analysis.flux, rm.vmin_low, rm.vmin_high,
              color="C0", lw=1.5, label="recovered flux")

    ax.set_xscale("log")
    ax.set_xlim(rm.vmin_low[0], rm.vmin_high[-1])
    ax.set_xlabel(r"$v_{min}$ [km/s]", fontsize=18)
    ax.set_ylabel(r"$\tilde{\eta}$  [cm$^{-1}$]", fontsize=18)
    ax.set_title(_stem(analysis).replace("_", " "), fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(True, which="both", ls="--", alpha=0.4)

    if save:
        out_dir = out_dir or (RESULTS_DIR / f"scenario_bkg_{analysis.config.background}")
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"flux_{_stem(analysis)}.pdf"
        plt.savefig(path, bbox_inches="tight")
        print(f"saved {path}")
    return ax
