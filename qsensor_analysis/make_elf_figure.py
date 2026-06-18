"""Paper figure: energy-loss function Im[-1/eps(E', q)] vs momentum transfer q.

Self-contained reproduction. Exports the ELF from the project's own dielectric
definitions via wolframscript, then plots:

  * TES (Al)  -- Mermin tabulated eps (Al_mermin.dat): ELF = eps2/(eps1^2+eps2^2),
                 at E' = 0.1 eV (threshold) and 1.0 eV.
  * MKID (TiN) -- analytic Lindhard ImepsLTiN[w,q] (= Im[-1/eps] directly),
                 at E' = 0.2 eV (threshold) and 1.0 eV.

The loss function is positive (passivity) for both -- see the values printed.

Usage:  python make_elf_figure.py
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "Mathematica" / "src"
OUT = HERE / "results" / "final" / "elf"

WOLFRAM = shutil.which("wolframscript") or \
    "/Applications/Wolfram.app/Contents/MacOS/wolframscript"

# q grid [eV, natural units]: within the Al Mermin interpolation domain.
QGRID = "Table[10.^e, {e, Log10[37.3], Log10[37289.], 0.005}]"

TES_WLS = f'''
src = "{SRC / 'TES'}";
Get[FileNameJoin[{{src, "01_setup.wl"}}]];
Get[FileNameJoin[{{src, "04_material.wl"}}]];
elf[w_, q_] := eps2Alf[w, q]/(eps1Alf[w, q]^2 + eps2Alf[w, q]^2);
qs = {QGRID};
Export["{{csv}}", Prepend[Table[{{q, elf[0.1, q], elf[1.0, q]}}, {{q, qs}}],
  {{"q", "ELF_lo", "ELF_hi"}}]];
'''

MKID_WLS = f'''
src = "{SRC / 'MKID'}";
Get[FileNameJoin[{{src, "01_setup.wl"}}]];
Get[FileNameJoin[{{src, "04_material.wl"}}]];
qs = {QGRID};
Export["{{csv}}", Prepend[Table[{{q, Re@N@ImepsLTiN[0.2, q], Re@N@ImepsLTiN[1.0, q]}}, {{q, qs}}],
  {{"q", "ELF_lo", "ELF_hi"}}]];
'''


def export(wls_template: str) -> np.ndarray:
    with tempfile.TemporaryDirectory() as d:
        csv = Path(d) / "elf.csv"
        wls = Path(d) / "elf.wls"
        wls.write_text(wls_template.replace("{csv}", str(csv)))
        subprocess.run([WOLFRAM, "-file", str(wls)], check=True,
                       capture_output=True, text=True)
        return np.loadtxt(csv, delimiter=",", skiprows=1)


def main():
    tes = export(TES_WLS)
    mkid = export(MKID_WLS)
    for name, d in (("TES (Al)", tes), ("MKID (TiN)", mkid)):
        print(f"{name}: ELF min lo/hi = {d[:, 1].min():.2e} / {d[:, 2].min():.2e} "
              f"(>0 => Im[-1/eps] positive)")

    try:
        plt.style.use(str(HERE / "physrev.mplstyle"))
    except Exception:
        pass
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    specs = [(axes[0], tes, "TES (Al)", "tab:red", "0.1"),
             (axes[1], mkid, "MKID (TiN)", "tab:blue", "0.2")]
    for ax, d, title, color, elo in specs:
        ax.plot(d[:, 0], d[:, 1], "-", color=color, lw=2,
                label=rf"$E'={elo}$ eV")
        ax.plot(d[:, 0], d[:, 2], "--", color=color, lw=2, label=r"$E'=1.0$ eV")
        ax.set_title(title, fontsize=20)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$q$ [eV]", fontsize=22)
        ax.legend(fontsize=15, frameon=False)
        ax.grid(True, which="both", ls="--", alpha=0.3)
        ax.tick_params(labelsize=16)
    axes[0].set_ylabel(r"$\mathrm{Im}[-1/\varepsilon(E',q)]$", fontsize=22)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "elf_vs_q.pdf"
    fig.savefig(out, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
