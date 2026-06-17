"""Overlay plots: smooth response functions BEFORE vs AFTER the production
stage-10 cut (alpha high-side window + analytic low-v floor).

Reproduces the exact cut that Mathematica/src/{TES,MKID}/10_response_matrix.wl
applies, but on the already-exported smooth responses in
output/<DET>/response_functions_csv/*.csv (the .wdx are NOT touched). For each
config it writes <config>_cut.pdf into output/<DET>/response_function_plots/,
in the style of sandbox/central_area_cut_test/response_cut.pdf: one panel per
energy bin, original (grey) overlaid with the kept [floor, b] response (red).

Cut definition (peak-value cut, outermost crossings):
  * window [a, b]      : a is the first v_min (scanning up from the low end)
    where R(v_min) >= `frac` of the peak value, b is the last (scanning down
    from the high end); default frac = 0.20. Only the outer sub-frac tails are
    cut -- an interior dip below the threshold (the valley of a double-peaked
    response) is kept, so neither peak is removed. Found on a 4000-interval fine
    grid over [vminLo, vminHi].
  * kept response      : original on [a, b], zero elsewhere.

Usage:  python _make_cut_overlays.py [frac]   (default frac = 0.20)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
NFINE = 4000                             # matches stage-10 fine grid

DETECTORS = {
    "TES":  HERE / "TES",
    "MKID": HERE / "MKID",
}


def mass_tag(name: str) -> str:
    m = re.search(r"M(\d+)", name)
    return f"M{m.group(1)}" if m else "other"


def peak_window(v_fine, r_fine, frac):
    """Peak-value window [a, b] = the OUTERMOST crossings of frac * peak:
    a is the first v_min (scanning up from the low end) where R >= frac * peak,
    b is the last (scanning down from the high end). Only the outer sub-frac
    tails are cut; an interior dip below the threshold (e.g. the valley of a
    double-peaked response) is kept, so neither peak is removed."""
    if r_fine.max() <= 0:
        return v_fine[0], v_fine[0]      # empty response
    thr = frac * r_fine.max()
    idx = np.nonzero(r_fine >= thr)[0]
    return v_fine[idx[0]], v_fine[idx[-1]]


def process(csv_path: Path, out_dir: Path, frac: float):
    name = csv_path.stem
    raw = np.genfromtxt(csv_path, delimiter=",", names=True, deletechars="",
                        encoding="utf-8")
    cols = list(raw.dtype.names)
    vcol = cols[0]
    bins = cols[1:]
    v = np.asarray(raw[vcol], dtype=float)
    R = np.array([np.asarray(raw[b], dtype=float) for b in bins])

    mt = mass_tag(name)

    # fine grid for the window edges (matches stage-10)
    v_fine = np.linspace(v.min(), v.max(), NFINE + 1)

    fig, axes = plt.subplots(1, len(bins), figsize=(4.0 * len(bins), 4.0),
                             sharex=True)
    if len(bins) == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        r = R[i]
        r_fine = np.interp(v_fine, v, r)          # order-1, like the IF
        if frac > 0:
            a, b = peak_window(v_fine, r_fine, frac)
        else:
            a, b = v.min(), v.max()
        keep = (v >= a) & (v <= b)
        r_cut = np.where(keep, r, 0.0)

        ax.plot(v, r, color="0.7", lw=1.2, label="original")
        ax.plot(v, r_cut, color="C3", lw=2.0,
                label=f"cut: keep [{a:.0f}, {b:.0f}]")
        ax.axvspan(a, b, color="C3", alpha=0.10)
        ax.set_xscale("log")
        ax.set_xlim(v.min(), v.max())
        ax.set_ylim(bottom=0)
        ax.set_title(bins[i].replace("_", " "))
        ax.set_xlabel(r"$v_{\min}$ [km/s]")
        if i == 0:
            ax.set_ylabel(r"$R_{\rm bin}(v_{\min})$  [kg$^{-1}$]")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=8)

    fig.suptitle(
        f"{name} ({mt}): response functions, original vs cut "
        f"(both-side cut at {int(round(frac*100))}% of peak value)",
        fontsize=13)
    fig.tight_layout()
    out = out_dir / f"{name}_cut.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out.name


def main():
    frac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.20
    for det, base in DETECTORS.items():
        csv_dir = base / "response_functions_csv"
        out_dir = base / "response_function_plots"
        out_dir.mkdir(parents=True, exist_ok=True)
        if not csv_dir.is_dir():
            print(f"[{det}] no response_functions_csv -- skip")
            continue
        print(f"[{det}] both-side cut at {int(round(frac*100))}% of peak")
        for csv_path in sorted(csv_dir.glob("*.csv")):
            fname = process(csv_path, out_dir, frac)
            print(f"  saved {fname}")


if __name__ == "__main__":
    main()
