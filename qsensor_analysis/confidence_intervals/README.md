# Confidence-interval notebooks

Production notebooks for the point-wise profile confidence bands of the
recovered flux x(v_min) — one notebook per **(detector/material, mediator,
DM mass)**, with one subsection per **halo model** inside each notebook.

| | heavy mediator (q0) | light mediator (q2) |
|---|---|---|
| **TES (Al)** | `CI_TES_Al_q0_M{1,2,3}.ipynb` | `CI_TES_Al_q2_M{1,2,3}.ipynb` |
| **MKID (TiN)** | `CI_MKID_TiN_q0_M{1,2,3}.ipynb` | `CI_MKID_TiN_q2_M{1,2,3}.ipynb` |

Mass tags: M1 = 10 MeV, M2 = 100 MeV, M3 = 1 GeV. Binning: 5 bins (R5).

Subsections per notebook: **SHM (Halo)** / **pure Disk** / **mixtures
(5% / 25% disk)**, plus **Halo+Bound** in the M3 notebooks (the Earth-bound
population only reaches the 1 GeV window).

## Usage

- Each `run_model(...)` cell is a full computation (**~30–60 min**, production
  MC settings `num_pseudo=50, n_pseudo_edge=500, rel_tol=0.02`, 30 scan
  points) and saves everything under
  `../results/<DET>/bkg-none/<config>/` (per-point JSONs in `band/`, summary
  CSV, figure PDF).
- `show_saved(...)` reloads a saved band — table + figure, no recomputation.
- Method documentation: `quantum_sensor.statistics.find_confidence_band`.
  Interactive walkthrough / sandbox: `../confidence_bands_test.ipynb`.

## Regenerating the notebooks

The notebooks are generated from a single template:

```
python _generate_notebooks.py
```

**This overwrites all 12 notebooks** — keep results in `results/` (they are),
not inside the notebooks, and put template changes into
`_generate_notebooks.py` rather than editing notebooks by hand.
