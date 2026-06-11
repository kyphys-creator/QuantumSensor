"""Generate the production confidence-interval notebooks.

One notebook per (material, mediator, DM mass) — 12 in total — each with one
subsection per halo model (SHM / Disk / mixtures, plus Bound for the 1 GeV
notebooks). Rerun this script to regenerate all notebooks after a template
change (existing notebooks are overwritten, so keep results in results/, not
in the notebooks).

    python _generate_notebooks.py
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

DETECTOR = {"Al": "TES", "TiN": "MKID"}
MEDIATOR = {"0": "heavy mediator (F_DM = 1)", "2": "light mediator (F_DM ∝ 1/q²)"}
MASS_LABEL = {"1": "10 MeV", "2": "100 MeV", "3": "1 GeV"}
NBINS = {"Al": 5, "TiN": 5}          # production binning: 5 bins


def md(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": [l + "\n" for l in text.splitlines()]}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": [l + "\n" for l in text.splitlines()]}


def build_notebook(material, q, mass):
    det, med, mlab = DETECTOR[material], MEDIATOR[q], MASS_LABEL[mass]
    nbins = NBINS[material]
    cells = [

md(f"""# Confidence intervals — {det} ({material}), {med}, m_DM = {mlab}

Point-wise profile confidence bands of the recovered flux x(v_min), one
subsection per halo model. Method and conventions:
[`quantum_sensor.statistics.find_confidence_band`](../src/quantum_sensor/statistics.py)
(neutrinoAnalysis construction: observed profile Δχ² vs an MC-calibrated
cutoff; Poisson pseudo-experiments; CLARABEL solver, thread-parallel).

- Configuration: `material='{material}', q='{q}', mass='{mass}', nbins={nbins}`, background `none`
- Outputs per model: `results/{det}/bkg-none/<config>/` —
  `band/band_idx*.json` (primary, one per point), `flux_profile_band.csv`
  (summary), `flux_profile_band.pdf` (figure)
- **Each `run_model` cell takes ~30–60 min** (production MC settings).
  Already computed? Use `show_saved` in §Saved results instead."""),

code(f"""import numpy as np
import matplotlib.pyplot as plt

from quantum_sensor import (DarkMatterQuantumAnalysis, RunConfig,
                            pointwise_band, load_pointwise_band,
                            band_table, save_band_products)
from quantum_sensor.plotting import plot_flux_with_pointwise_bands, run_dir

MATERIAL, Q, MASS, NBINS = '{material}', '{q}', '{mass}', {nbins}
LEVELS = (0.68, 0.954)                                   # 1 sigma, 2 sigma
MC = dict(num_pseudo=50, n_pseudo_edge=500, rel_tol=0.02)  # production settings
N_INDICES = 30


def build(eta='Halo', disk_fraction=None):
    cfg = RunConfig(material=MATERIAL, q=Q, mass=MASS, nbins=NBINS,
                    eta=eta, disk_fraction=disk_fraction, background='none')
    a = DarkMatterQuantumAnalysis(cfg)
    a.optimize()
    return a


def run_model(eta='Halo', disk_fraction=None):
    \"\"\"Compute, save and display the band for one halo model (heavy).\"\"\"
    a = build(eta, disk_fraction)
    print(f'{{a.config}}')
    print(f'counts per bin: {{np.array2string(a.observed, precision=1)}}')
    bands = pointwise_band(a, n_indices=N_INDICES, levels=LEVELS, **MC)
    save_band_products(a, bands)
    display(band_table(bands))
    return a, bands


def show_saved(eta='Halo', disk_fraction=None):
    \"\"\"Reload a saved band (band/*.json) — no recomputation.\"\"\"
    a = build(eta, disk_fraction)
    bands = load_pointwise_band(run_dir(a))
    print(f'{{len(bands)}} points from {{run_dir(a) / "band"}}')
    display(band_table(bands))
    plot_flux_with_pointwise_bands(a, bands, save=False)
    return a, bands"""),

md("""## 1. Standard halo (SHM)"""),
code("""a_halo, b_halo = run_model(eta='Halo')"""),

md("""## 2. Pure dark disk"""),
code("""a_disk, b_disk = run_model(eta='Disk')"""),

md("""## 3. Halo + dark-disk mixtures

The fit eta is `(1-p)·Halo + p·Disk` with the same total local DM density."""),
code("""a_mix5, b_mix5 = run_model(disk_fraction=0.05)    # 5% disk"""),
code("""a_mix25, b_mix25 = run_model(disk_fraction=0.25)  # 25% disk"""),
    ]

    if mass == "3":
        cells += [
md("""## 4. Earth-bound population (Halo + Bound)

The bound population sits on top of the full SHM halo; it is nonzero only
below v_esc = 11.2 km/s, which only the 1 GeV window reaches."""),
code("""a_bound, b_bound = run_model(eta='Bound')"""),
        ]

    cells += [
md("""## Saved results (no recomputation)

Reload any model computed above (or in a previous session)."""),
code("""show_saved(eta='Halo');
# show_saved(eta='Disk');
# show_saved(disk_fraction=0.05);
# show_saved(disk_fraction=0.25);"""
     + ("\n# show_saved(eta='Bound');" if mass == "3" else "")),
    ]

    return {"cells": cells,
            "metadata": {"kernelspec": {"display_name": "Python 3",
                                        "language": "python",
                                        "name": "python3"},
                         "language_info": {"name": "python", "version": "3"}},
            "nbformat": 4, "nbformat_minor": 5}


if __name__ == "__main__":
    for material in ("Al", "TiN"):
        for q in ("0", "2"):
            for mass in ("1", "2", "3"):
                nb = build_notebook(material, q, mass)
                name = f"CI_{DETECTOR[material]}_{material}_q{q}_M{mass}.ipynb"
                with open(HERE / name, "w") as f:
                    json.dump(nb, f, indent=1, ensure_ascii=False)
                print("wrote", name)
