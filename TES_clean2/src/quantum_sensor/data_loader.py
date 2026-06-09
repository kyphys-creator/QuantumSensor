"""Load the inputs of the analysis.

Two data sources:

* The new Mathematica response matrices, written by 10_response_matrix.wl as a
  per-run folder ``response_matrix/M<mass>/<name>/`` holding ``matrix.csv``,
  ``vmin.csv`` and ``bins.csv``. Shapes, energy bins and the v_min grid are all
  read from these files -- nothing about the grid is hard-coded here.
* The legacy ``data/`` folder (eta velocity distributions, observed Ratebin
  counts, and the old exposure-weighted ``curlyRplotting`` matrices kept for
  unit/exposure calibration).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

# quantum_sensor/ -> src/ -> TES_clean2/ -> QuantumSensor/
_PKG_DIR = Path(__file__).resolve()
LEGACY_DATA_DIR = _PKG_DIR.parents[2] / "data"
MATRIX_ROOT = _PKG_DIR.parents[3] / "Mathematica" / "output" / "TES" / "response_matrix"


@dataclass(frozen=True)
class ResponseMatrix:
    """A response matrix plus the grids that label its axes.

    ``matrix`` has shape (n_ebins, n_vmin). Rows are observed-energy bins,
    columns are v_min intervals. The matrix is the raw v_min integral of the
    response (natural units); exposure and unit factors are applied later in
    :mod:`model`.
    """

    matrix: np.ndarray      # (n_ebins, n_vmin)
    vmin_low: np.ndarray    # (n_vmin,)  [km/s]
    vmin_high: np.ndarray   # (n_vmin,)  [km/s]
    vmin_mid: np.ndarray    # (n_vmin,)  [km/s]
    ebin_low: np.ndarray    # (n_ebins,) [eV]
    ebin_high: np.ndarray   # (n_ebins,) [eV]
    name: str
    path: Path

    @property
    def n_ebins(self) -> int:
        return self.matrix.shape[0]

    @property
    def n_vmin(self) -> int:
        return self.matrix.shape[1]


def find_matrix_dir(material: str, q: str, mass: str, nbins: int,
                    run: str | None = None) -> Path:
    """Locate the response-matrix folder for one configuration.

    ``run`` optionally selects a specific run (a substring of the folder name,
    e.g. ``"v1-800_N1000"``) when several exist for the same configuration.
    """
    stem = f"{material}_q{q}M{mass}_R{nbins}"
    mass_dir = MATRIX_ROOT / f"M{mass}"
    matches = sorted(mass_dir.glob(f"{stem}_*"))
    if run is not None:
        matches = [m for m in matches if run in m.name]
    if not matches:
        raise FileNotFoundError(
            f"no response-matrix folder for {stem}"
            + (f" matching {run!r}" if run else "")
            + f" under {mass_dir}"
        )
    if len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        raise ValueError(
            f"{len(matches)} matrix folders match {stem}; pass run= to disambiguate: {names}"
        )
    return matches[0]


def load_response_matrix(material: str, q: str, mass: str, nbins: int,
                         run: str | None = None) -> ResponseMatrix:
    """Load a new Mathematica response matrix with its v_min and energy grids."""
    folder = find_matrix_dir(material, q, mass, nbins, run=run)
    matrix = np.atleast_2d(np.genfromtxt(folder / "matrix.csv", delimiter=","))
    vmin = np.atleast_2d(np.genfromtxt(folder / "vmin.csv", delimiter=","))
    ebins = np.atleast_2d(np.genfromtxt(folder / "bins.csv", delimiter=","))
    return ResponseMatrix(
        matrix=matrix,
        vmin_low=vmin[:, 0], vmin_high=vmin[:, 1], vmin_mid=vmin[:, 2],
        ebin_low=ebins[:, 0], ebin_high=ebins[:, 1],
        name=folder.name, path=folder,
    )


def load_natural_eta(rm: ResponseMatrix, model: str) -> np.ndarray:
    """Natural-units eta(v_min) sampled on the matrix's own v_min grid.

    Read from ``eta_<model>.csv`` inside the matrix folder, produced by the
    Mathematica stage ``12_eta.wl``. Because it is generated on the same v_min
    interval mid-points that label the matrix columns and in the same natural
    units as the matrix, ``M_phys @ eta`` is a real expected event count -- no
    interpolation (``align_eta``) and no unit conversion are needed or applied.

    Only the Halo (SHM) model is currently generated; see ``12_eta.wl``.
    """
    path = rm.path / f"eta_{model}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"natural-units eta not found: {path}\n"
            f"Generate it first with the Mathematica stage, e.g.:\n"
            f"  wolframscript -file Mathematica/src/TES/12_eta.wl {rm.name} {model}\n"
            f"(only the 'Halo' model is implemented)."
        )
    eta = np.atleast_1d(np.genfromtxt(path, delimiter=","))
    if eta.shape[0] != rm.n_vmin:
        raise ValueError(
            f"eta length {eta.shape[0]} != matrix n_vmin {rm.n_vmin} for {path}; "
            f"regenerate eta_{model}.csv with 12_eta.wl after rebuilding the matrix."
        )
    return eta


# ----------------------------------------------------------------- legacy data

def load_eta(eta: str, mass: str) -> np.ndarray:
    """Tabulated velocity distribution eta(v_min) (legacy, 800 samples).

    Physical units (~cm^-1), NOT consistent with the new natural-units
    matrices -- kept only for cross-checks. The analysis uses
    :func:`load_natural_eta` instead.
    """
    return np.genfromtxt(LEGACY_DATA_DIR / "Eta_data" / f"Eta{eta}M{mass}_Ko.csv",
                         delimiter=",")


def load_ratebin(eta: str, material: str, q: str, mass: str) -> np.ndarray:
    """Observed event counts per energy bin (legacy, R5 / q0 only)."""
    return np.genfromtxt(
        LEGACY_DATA_DIR / "Ratebin" / f"Event5Eta{eta}Mat{material}q{q}M{mass}.csv",
        delimiter=",")


def load_legacy_response_matrix(material: str, q: str, mass: str) -> np.ndarray:
    """Old exposure-weighted response matrix (curlyRplotting). Kept only to
    calibrate the exposure/unit factor against the new raw matrix."""
    path = (LEGACY_DATA_DIR / "curlyRplotting"
            / f"{material}_R5_q{q}M{mass}_Ko_thresholded_refined_exposure.csv")
    return np.genfromtxt(path, delimiter=",")
