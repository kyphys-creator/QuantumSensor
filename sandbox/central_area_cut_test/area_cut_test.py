"""Best-fit with a two-sided (peak-centred, 50%-area) response-function cut.

Isolated sandbox (see README). TES Al q0 M1 R5 SHM. For each bin's original
response function R_bin(v_min) we keep only the highest-density region around
the peak covering alpha=50% of the area (cutting BOTH the low-v rise and the
high-v tail), build the matrix from the cut responses, and run the production
best-fit. The full uncut matrix is run alongside as the reference.

    python area_cut_test.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# production pieces (read-only)
from quantum_sensor.optimizer import run_optimize_qp
from quantum_sensor.model import condition
from quantum_sensor.config import CONDITION_C
from quantum_sensor.eta_models import eta as eta_model
from quantum_sensor.constants import DM_MASS

ALPHA = 0.50                 # keep this fraction of the area around the peak
MASS_TAG = "1"               # M1 = 10 MeV
CSV = (Path(__file__).resolve().parents[2]
       / "Mathematica/output/TES/response_functions_csv/Al_q0M1_R5.csv")


def trap_weights(v):
    """Trapezoidal integration weights on a (possibly log-spaced) grid."""
    w = np.empty_like(v)
    w[1:-1] = 0.5 * (v[2:] - v[:-2])
    w[0] = 0.5 * (v[1] - v[0])
    w[-1] = 0.5 * (v[-1] - v[-2])
    return w


def central_area_cut(v, r, alpha=ALPHA):
    """Keep the highest-density region around the peak covering `alpha` of the
    area; zero the rest (both tails). Returns (r_cut, (v_lo, v_hi), kept_frac).

    Grows a contiguous window outward from the peak, each step annexing the
    neighbouring sample with the larger response value, until the accumulated
    area reaches `alpha` of the total."""
    w = trap_weights(v)
    area = r * w
    total = area.sum()
    if total <= 0:
        return np.zeros_like(r), (np.nan, np.nan), 0.0
    pk = int(np.argmax(r))
    lo = hi = pk
    acc = area[pk]
    n = len(v)
    while acc < alpha * total:
        left, right = lo - 1, hi + 1
        lv = r[left] if left >= 0 else -np.inf
        rv = r[right] if right < n else -np.inf
        if lv == -np.inf and rv == -np.inf:
            break
        if lv >= rv:
            lo = left
            acc += area[lo]
        else:
            hi = right
            acc += area[hi]
    out = np.zeros_like(r)
    out[lo:hi + 1] = r[lo:hi + 1]
    return out, (v[lo], v[hi]), acc / total


def best_fit(M, eta, target_counts=1000.0, vertex_select=True):
    """Production best-fit on noiseless self-consistent data; scale arbitrary.

    ``vertex_select=True`` is the production minimal-total-flux staircase
    vertex; ``False`` returns the OSQP interior solution (the smooth solution
    on the chi2=0 face), which does not shave the lowest-weight high-v column.
    """
    raw = M @ eta
    exposure = target_counts / float(raw.max())
    m_phys = exposure * M
    signal = m_phys @ eta
    m_cond, data, bkg, unscale = condition(m_phys, signal,
                                            np.zeros_like(signal), c=CONDITION_C)
    res = run_optimize_qp(m_cond, data, bkg, M.shape[1],
                          vertex_select=vertex_select)
    return unscale(res.x), res


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])               # (n_bins, n_v) original
    w = trap_weights(v)
    eta = eta_model("Halo", DM_MASS[MASS_TAG], v)    # SHM on the same grid

    # --- central area cut per bin ---
    R_cut = np.zeros_like(R)
    spans = []
    for i in range(len(bins)):
        R_cut[i], span, frac = central_area_cut(v, R[i], ALPHA)
        spans.append(span)
        print(f"{bins[i]:10s} peak@{v[np.argmax(R[i])]:6.1f}  "
              f"kept v in [{span[0]:6.1f}, {span[1]:6.1f}] km/s  "
              f"(area {frac*100:.0f}%)")

    # --- figure 1: original vs cut response functions ---
    fig, axes = plt.subplots(1, len(bins), figsize=(20, 4), sharex=True)
    for i, ax in enumerate(axes):
        ax.plot(v, R[i], color="0.7", lw=1.2, label="original")
        ax.plot(v, R_cut[i], color="C3", lw=2.0, label=f"central {int(ALPHA*100)}% area")
        ax.axvspan(spans[i][0], spans[i][1], color="C3", alpha=0.10)
        ax.set_xscale("log")
        ax.set_xlim(v.min(), v.max())
        ax.set_title(bins[i].replace("_", " "))
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        if i == 0:
            ax.set_ylabel(r"$R_{bin}(v_{min})$")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle("Response functions: original vs central-50%-area cut "
                 "(TES Al q0 M1 R5)", fontsize=13)
    fig.tight_layout()
    fig.savefig("response_cut.pdf", bbox_inches="tight")
    print("saved response_cut.pdf")

    # --- build matrices (grid columns, trapezoidal weights) and best-fit ---
    def matrix(Rmat):
        M = Rmat * w[None, :]                         # M[i,j] = R_i(v_j) * w_j
        keep = M.any(axis=0)                          # trim all-zero columns
        return M[:, keep], keep

    results = {}
    for name, Rmat in (("full (uncut)", R), (f"central {int(ALPHA*100)}% cut", R_cut)):
        M, keep = matrix(Rmat)
        x, res = best_fit(M, eta[keep])
        results[name] = (v[keep], x, eta[keep], res)
        rel_last = x[-3:] / eta[keep][-3:]
        print(f"[{name}] {M.shape[0]}x{M.shape[1]}  backend={res.backend}  "
              f"kept v [{v[keep].min():.0f},{v[keep].max():.0f}]  "
              f"last3 x/eta = {np.array2string(rel_last, precision=2, floatmode='fixed')}")

    # --- figure 2: best-fit vs eta ---
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), sharey=True)
    vg = np.logspace(0, np.log10(v.max()), 400)
    eta_ref = eta_model("Halo", DM_MASS[MASS_TAG], vg)
    for ax, (name, (vk, x, ek, res)) in zip(axes, results.items()):
        ax.plot(vg, eta_ref, color="red", lw=2, label="input eta (SHM)")
        ax.step(vk, x, where="mid", color="C0", lw=1.5, label="best-fit")
        ax.axvspan(vk.min(), vk.max(), color="C0", alpha=0.07,
                   label="matrix v_min window")
        ax.set_xscale("log")
        ax.set_xlim(v.min(), v.max())
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        ax.set_title(f"{name}\n(backend: {res.backend})")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=9)
    axes[0].set_ylabel(r"$\tilde{\eta}$ (natural units)")
    fig.suptitle("Best-fit recovery: full matrix vs central-area-cut matrix "
                 "(TES Al q0 M1 R5, SHM)", fontsize=13)
    fig.tight_layout()
    fig.savefig("best_fit.pdf", bbox_inches="tight")
    print("saved best_fit.pdf")


if __name__ == "__main__":
    main()
