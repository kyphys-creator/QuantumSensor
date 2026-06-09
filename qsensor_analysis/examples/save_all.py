"""Regenerate the organized outputs for every available configuration.

Writes one folder per run under ``results/`` following the layout

    results/<DET>/bkg-<background>/<material>_q<q>_M<mass>_R<nbins>_<eta>/
        flux.csv   (recovered flux + v_min grid, natural units; tracked in git)
        flux.pdf   (flux-vs-eta comparison figure; gitignored)

``<DET>`` is TES (material Al) or MKID (material TiN). Only configurations whose
response matrices have actually been generated on the Mathematica side are run;
missing ones are skipped with a note (e.g. MKID light mediator q2, which needs
the MKID q2 matrices built first).

    python examples/save_all.py

Pass ``--no-pdf`` to write only the CSVs (skip the figures).
"""

from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")          # headless: save figures without a display
import matplotlib.pyplot as plt

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.data_loader import DETECTOR_OF


# (material, q, nbins) combinations that exist as response matrices.
# TiN (MKID) threshold is 0.2 eV, so its fine binning is R9 (not R10), and only
# the heavy mediator (q0) has been generated so far.
MATERIAL_GRID = {
    "Al":  [("0", 5), ("0", 10), ("2", 5), ("2", 10)],
    "TiN": [("0", 5), ("0", 9)],
}
MASSES = ("1", "2", "3")
# Pure velocity models rendered as references (need eta_<model>.csv from 12_eta.wl).
ETAS = ("Halo", "Disk")
# Dark-disk mixing fractions p for the best-fit eta = (1-p)*Halo + p*Disk.
DISK_FRACTIONS = (0.05, 0.25)
# Background scenarios to render. "none" is signal-only; add others as needed.
BACKGROUNDS = ("none",)

# eta specs as (eta_model, disk_fraction): pure references + mixture best-fits.
ETA_SPECS = [(m, None) for m in ETAS] + [("Halo", p) for p in DISK_FRACTIONS]


def main(write_pdf: bool = True) -> None:
    n_ok = n_skip = 0
    for material, combos in MATERIAL_GRID.items():
        det = DETECTOR_OF.get(material, material)
        for q, nbins in combos:
            for mass in MASSES:
                for eta, disk_fraction in ETA_SPECS:
                    for background in BACKGROUNDS:
                        cfg = RunConfig(material=material, q=q, mass=mass,
                                        nbins=nbins, eta=eta,
                                        disk_fraction=disk_fraction, background=background)
                        try:
                            a = DarkMatterQuantumAnalysis(cfg)
                        except FileNotFoundError as exc:
                            tag = eta if disk_fraction is None else f"mix{round(disk_fraction*100)}"
                            print(f"skip {det} {material} q{q} M{mass} R{nbins} {tag}: {exc}")
                            n_skip += 1
                            continue
                        a.optimize(solver="osqp")
                        csv = a.save_flux()
                        if write_pdf:
                            fig, ax = plt.subplots(figsize=(8, 6))
                            a.plot(save=True, ax=ax)
                            plt.close(fig)
                        print(f"ok   {csv.parent.relative_to(csv.parents[3])}")
                        n_ok += 1
    print(f"\n{n_ok} runs written, {n_skip} skipped.")


if __name__ == "__main__":
    main(write_pdf="--no-pdf" not in sys.argv[1:])
