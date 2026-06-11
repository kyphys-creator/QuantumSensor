"""Point-wise profile confidence bands (neutrinoAnalysis construction).

For each fine-binning Halo configuration, scan ~12 log-spaced v_min indices
over the populated window and locate the 68% / 95.4% point-wise profile
confidence interval of each flux step (Feldman-Cousins-style: observed
profile Delta-chi^2 against an MC-calibrated cutoff, edges by bracketing +
geometric bisection -- see ``quantum_sensor.statistics.find_confidence_band``).

Writes, next to each run's ``flux.csv``:

    flux_profile_band.csv    vmin_mid, best_fit, lo68, hi68, lo95, hi95
    band/band_idx*.json      one full record per scanned point (cutoffs, evals)
    flux_profile_band.pdf    eta + staircase + shaded point-wise bands

    python examples/run_pointwise_bands.py                    # 12-config grid
    python examples/run_pointwise_bands.py Al 0 2 10          # one config
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.statistics import pointwise_band, save_pointwise_band
from quantum_sensor.plotting import plot_flux_with_pointwise_bands, run_dir

LEVELS = (0.68, 0.954)
N_INDICES = 12


def run_one(material, q, mass, nbins):
    cfg = RunConfig(material=material, q=q, mass=mass, nbins=nbins,
                    eta="Halo", background="none")
    a = DarkMatterQuantumAnalysis(cfg)
    a.optimize()
    t0 = time.time()
    bands = pointwise_band(a, n_indices=N_INDICES, levels=LEVELS,
                           num_pseudo=30, n_pseudo_edge=200, rel_tol=0.02)
    print(f"{material} q{q} M{mass} R{nbins}: "
          f"{len(bands)} indices in {time.time() - t0:.0f}s")

    out = run_dir(a)
    out.mkdir(parents=True, exist_ok=True)

    table = np.array([[b["vmin_mid"], b["best_fit"],
                       b["band"][LEVELS[0]][0], b["band"][LEVELS[0]][1],
                       b["band"][LEVELS[1]][0], b["band"][LEVELS[1]][1]]
                      for b in bands])
    np.savetxt(out / "flux_profile_band.csv", table, delimiter=",", comments="",
               header="vmin_mid,best_fit,lo68,hi68,lo95,hi95")

    save_pointwise_band(a, bands)

    plot_flux_with_pointwise_bands(a, bands)
    plt.close("all")


if __name__ == "__main__":
    if len(sys.argv) == 5:
        run_one(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        fine = {"Al": 10, "TiN": 9}
        for material, nbins in fine.items():
            for q in ("0", "2"):
                for mass in ("1", "2", "3"):
                    try:
                        run_one(material, q, mass, nbins)
                    except Exception as exc:
                        print(f"skip {material} q{q} M{mass} R{nbins}: {exc}")
