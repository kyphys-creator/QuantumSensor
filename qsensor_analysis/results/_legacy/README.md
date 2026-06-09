# results/_legacy

Outputs from **before** the `results/` reorganization. None of these are
produced by the current code; they are kept only for reference.

The current layout writes one folder per run:

```
results/<DET>/bkg-<background>/<material>_q<q>_M<mass>_R<nbins>_<eta>/
    flux.csv   # recovered flux + v_min grid (natural units)
    flux.pdf   # flux-vs-eta figure (gitignored)
```

Regenerate everything with `python examples/save_all.py`.

Contents here:

- `flux_*.csv`, `scenario_bkg_none/flux_*.pdf` — the old flat / `scenario_bkg_*`
  naming this folder replaced.
- `TES_*.pdf`, `MKID_*.pdf`, `Res5*.dat`, `data/`, `events/` — artifacts from
  earlier (pre-package) pipelines.
