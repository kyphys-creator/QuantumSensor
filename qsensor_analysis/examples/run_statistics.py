"""Full statistical analysis: discrimination grid + toy-ensemble flux bands.

Two products, written under ``results/stats/``:

1. ``discrimination.csv`` (+ ``discrimination_<bkg>.pdf`` heatmaps) -- the
   median expected significance (Asimov Poisson LLR, arXiv:1007.1727) for
   rejecting the pure-SHM Halo null when the alternative (dark-disk mixtures,
   pure Disk, Halo+Bound) is true, over every available configuration and
   background scenario, with a Monte-Carlo cross-check where p >~ 1/n_toys.

2. per-run ``flux_band.csv`` / ``flux_band.pdf`` (next to the existing
   ``flux.csv``) -- 68/95% percentile bands of the recovered staircase flux
   over Poisson pseudo-experiments, for the fine-binning Halo configurations.

    python examples/run_statistics.py            # both products
    python examples/run_statistics.py --no-toys  # significance grid only

Skips configurations whose response matrices / eta files are missing (same
behaviour as save_all.py).
"""

from __future__ import annotations

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.statistics import (ALTERNATIVES, discrimination_row,
                                       fit_toys, flux_band)
from quantum_sensor.plotting import RESULTS_DIR, plot_flux_band, run_dir
from quantum_sensor.data_loader import DETECTOR_OF

MATERIAL_GRID = {
    "Al":  [("0", 5), ("0", 10), ("2", 5), ("2", 10)],
    "TiN": [("0", 5), ("0", 9), ("2", 5), ("2", 9)],
}
MASSES = ("1", "2", "3")
BACKGROUNDS = ("none", "c", "b")
MC_TOYS = 200_000          # LLR Monte-Carlo cross-check sample
BAND_TOYS = 500            # pseudo-experiments per flux band
STATS_DIR = RESULTS_DIR / "stats"


def build(material, q, mass, nbins, background, **eta_kw):
    cfg = RunConfig(material=material, q=q, mass=mass, nbins=nbins,
                    background=background, **eta_kw)
    return DarkMatterQuantumAnalysis(cfg)


def significance_grid():
    rows = []
    for material, combos in MATERIAL_GRID.items():
        for q, nbins in combos:
            for mass in MASSES:
                for background in BACKGROUNDS:
                    try:
                        null = build(material, q, mass, nbins, background)
                    except Exception as exc:
                        print(f"skip {material} q{q} M{mass} R{nbins}: {exc}")
                        break
                    for name, delta in ALTERNATIVES.items():
                        try:
                            alt = build(material, q, mass, nbins, background,
                                        **delta)
                        except Exception:
                            continue          # e.g. Bound eta only for M3
                        row = discrimination_row(null, alt, name,
                                                 mc_toys=MC_TOYS, seed=42)
                        rows.append(row)
                        mc = ("inf" if row.sigma_mc is None
                              or not np.isfinite(row.sigma_mc)
                              else f"{row.sigma_mc:.2f}")
                        print(f"{material:3s} q{q} M{mass} R{nbins:<2d} "
                              f"bkg-{background:4s} vs {name:5s}: "
                              f"{row.sigma_asimov:6.2f} sigma (MC {mc})")
    return rows


def save_grid(rows):
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    path = STATS_DIR / "discrimination.csv"
    with open(path, "w") as f:
        f.write("material,q,mass,nbins,background,alternative,"
                "sigma_asimov,q_asimov,chi2_pearson,sigma_mc,p_value_mc\n")
        for r in rows:
            f.write(f"{r.material},{r.q},{r.mass},{r.nbins},{r.background},"
                    f"{r.alternative},{r.sigma_asimov:.6g},{r.q_asimov:.6g},"
                    f"{r.chi2_pearson:.6g},"
                    f"{'' if r.sigma_mc is None else f'{r.sigma_mc:.6g}'},"
                    f"{'' if r.p_value_mc is None else f'{r.p_value_mc:.6g}'}\n")
    print(f"saved {path}")
    return path


def heatmaps(rows):
    """One heatmap per background: rows = configs, cols = alternatives."""
    alts = list(ALTERNATIVES)
    for background in BACKGROUNDS:
        sub = [r for r in rows if r.background == background]
        if not sub:
            continue
        keys = sorted({(r.material, r.q, r.mass, r.nbins) for r in sub})
        grid = np.full((len(keys), len(alts)), np.nan)
        for r in sub:
            i = keys.index((r.material, r.q, r.mass, r.nbins))
            grid[i, alts.index(r.alternative)] = r.sigma_asimov

        fig, ax = plt.subplots(figsize=(6, 0.42 * len(keys) + 1.8))
        im = ax.imshow(grid, aspect="auto", cmap="viridis",
                       norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=100))
        ax.set_xticks(range(len(alts)), alts)
        ax.set_yticks(range(len(keys)),
                      [f"{DETECTOR_OF[m]} q{q} M{ms} R{nb}"
                       for m, q, ms, nb in keys], fontsize=8)
        for i in range(len(keys)):
            for j in range(len(alts)):
                if np.isfinite(grid[i, j]):
                    ax.text(j, i, f"{grid[i, j]:.1f}", ha="center",
                            va="center", fontsize=7,
                            color="white" if grid[i, j] < 10 else "black")
        ax.set_title(f"Median discrimination significance vs Halo "
                     f"[$\\sigma$], bkg-{background}", fontsize=11)
        fig.colorbar(im, ax=ax, label=r"$\sigma$ (Asimov)")
        path = STATS_DIR / f"discrimination_{background}.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        print(f"saved {path}")


def toy_bands():
    """68/95% flux bands for the fine-binning Halo configs, all backgrounds=none."""
    fine = {"Al": 10, "TiN": 9}
    for material, nbins in fine.items():
        for q in ("0", "2"):
            for mass in MASSES:
                try:
                    a = build(material, q, mass, nbins, "none")
                except Exception as exc:
                    print(f"skip band {material} q{q} M{mass}: {exc}")
                    continue
                a.optimize()
                ens = fit_toys(a, n_toys=BAND_TOYS, seed=7)
                band = flux_band(ens)
                plot_flux_band(a, band, ens.n_toys)
                plt.close("all")

                out = run_dir(a) / "flux_band.csv"
                lo68, hi68 = band[68.0]
                lo95, hi95 = band[95.0]
                table = np.column_stack([a.rm.vmin_mid, band["median"],
                                         lo68, hi68, lo95, hi95])
                np.savetxt(out, table, delimiter=",", comments="",
                           header="vmin_mid,median,lo68,hi68,lo95,hi95")
                print(f"saved {out}  (chi2 median {np.median(ens.chi2):.1f} "
                      f"for {a.n_ebins} bins, {ens.n_failed} failed)")


if __name__ == "__main__":
    rows = significance_grid()
    save_grid(rows)
    heatmaps(rows)
    if "--no-toys" not in sys.argv:
        toy_bands()
