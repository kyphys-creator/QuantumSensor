import numpy as np
import pandas as pd

from .constants import BINS, M_DETECTOR


def compute_bi(A: float, B: float, C: float) -> np.ndarray:
    """Integral of R(E) = A*exp(-E/B) + C over each bin in BINS."""
    bi = []
    for i in range(len(BINS) - 1):
        Ei, Ei1 = BINS[i], BINS[i + 1]
        integral = -A * B * (np.exp(-Ei1 / B) - np.exp(-Ei / B)) + C * (Ei1 - Ei)
        bi.append(M_DETECTOR * integral)
    return np.array(bi)


def prepare_backgrounds() -> pd.DataFrame:
    eV = 1
    keV = 1000
    days_in_year = 365
    scale = eV / keV * days_in_year

    ba = compute_bi(A=0 * scale, B=10, C=0 * scale)
    bb = compute_bi(A=460 * scale, B=10, C=100 * scale)
    b2 = compute_bi(A=920 * scale, B=10, C=200 * scale)
    bc = compute_bi(A=50 * scale, B=10, C=20 * scale)
    bflat = compute_bi(A=0 * scale, B=10, C=100 * scale)

    return pd.DataFrame({
        "Bin Start [eV]": BINS[:-1],
        "Bin End [eV]": BINS[1:],
        "b_i (a)": ba,
        "b_i (b)": bb,
        "b_i (c)": bc,
        "b_i (flat)": bflat,
        "b_i (b2)": b2,
    })


BKG_COLUMN = {
    'a': 'b_i (a)',
    'b': 'b_i (b)',
    'c': 'b_i (c)',
    'flat': 'b_i (flat)',
    'none': 'b_i (a)',
    'b2': 'b_i (b2)',
}
