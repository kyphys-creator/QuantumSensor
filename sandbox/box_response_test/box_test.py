"""Best-fit estimator test on tail-free 'box' response matrices.

Isolated sandbox (see README). Imports only the production *solver* path
(condition + run_optimize_qp) and the SHM halo model, read-only. Builds its own
synthetic response matrices, forward-models noiseless counts from a known
monotone eta, runs the best-fit inverse, and checks whether the recovered
staircase meets eta at the high-v_min edge of the window.

    python box_test.py
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- production pieces under test (read-only use) --------------------------
from quantum_sensor.optimizer import run_optimize_qp
from quantum_sensor.model import condition
from quantum_sensor.config import CONDITION_C
from quantum_sensor.eta_models import eta as eta_model
from quantum_sensor.constants import DM_MASS


# ---------------------------------------------------------------------------
# Geometry: rows = energy bins, cols = v_min intervals.
# The grid deliberately ENDS while eta is still positive (v_hi < v_esc+v_e),
# so the "last column" is a genuine test of the high-v_min edge -- exactly the
# situation created by the production per-row 1-alpha window cut.
# ---------------------------------------------------------------------------
N_EBINS = 5
N_VMIN = 60
V_LO, V_HI = 1.0, 400.0          # ends mid-eta (SHM reaches ~770 km/s)
MASS_TAG = "2"                    # 100 MeV


def vmin_grid(n=N_VMIN, lo=V_LO, hi=V_HI):
    edges = np.linspace(lo, hi, n + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    width = np.diff(edges)
    return mid, width


def thresholds(n_ebins=N_EBINS, lo=V_LO, hi=V_HI):
    """Per-energy-bin kinematic v_min threshold, spread across the grid.

    Higher energy deposition -> higher threshold, so the highest bin only sees
    the high-v_min part of eta (where the undershoot is observed)."""
    return np.linspace(lo, 0.75 * hi, n_ebins)


def build_matrix(kind, vmid, width, thr, tail_len=40.0, box_w=120.0):
    """Synthetic response matrix M[i,j] = R_i(vmid_j) * width_j.

    kind:
      'box_flat'   R_i = 1 for v >= thr_i, hard 0 past the grid edge (NO decay)
      'box_window' R_i = 1 on [thr_i, thr_i + box_w], else 0
      'tail'       R_i = exp(-(v-thr_i)/tail_len) for v >= thr_i  (control)
    """
    M = np.zeros((len(thr), len(vmid)))
    for i, t in enumerate(thr):
        above = vmid >= t
        if kind == "box_flat":
            r = np.where(above, 1.0, 0.0)
        elif kind == "box_window":
            r = np.where(above & (vmid <= t + box_w), 1.0, 0.0)
        elif kind == "tail":
            r = np.where(above, np.exp(-(vmid - t) / tail_len), 0.0)
        else:
            raise ValueError(kind)
        M[i] = r * width
    return M


def true_eta(vmid, mass_tag=MASS_TAG):
    """SHM eta on the grid (production model), monotone non-increasing."""
    return eta_model("Halo", DM_MASS[mass_tag], vmid)


def best_fit(M, eta, target_counts=1000.0):
    """The production best-fit: exact-fit minimal-total-flux monotone vertex.

    Noiseless self-consistent data (data = exposure * M @ eta) so this
    exercises the Stage-0 exact-fit vertex -- i.e. the best-fit estimator
    itself, with no Poisson noise in the way. An exposure scales the raw
    response to realistic counts (~target_counts max), as the production
    forward model does; the recovered flux is independent of it."""
    exposure = target_counts / float((M @ eta).max())
    m_phys = exposure * M
    signal = m_phys @ eta
    m_cond, data, bkg, unscale = condition(m_phys, signal,
                                           np.zeros_like(signal), c=CONDITION_C)
    res = run_optimize_qp(m_cond, data, bkg, M.shape[1])
    return unscale(res.x), res


def plot_matrices(vmid, width, thr, kinds):
    """Show the test response matrices actually fed to the best-fit.

    Top row: per-energy-bin response curves R_i(v_min) = M[i,j]/width_j (the
    box vs tail shape). Bottom row: the matrix M[i,j] as a heatmap."""
    fig, axes = plt.subplots(2, len(kinds), figsize=(18, 8),
                             gridspec_kw={"height_ratios": [1, 1.1]})
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(thr)))
    for c, kind in enumerate(kinds):
        M = build_matrix(kind, vmid, width, thr)
        R = M / width                       # back out R_i(v_min) (drop integ. weight)

        ax = axes[0, c]
        for i, t in enumerate(thr):
            ax.plot(vmid, R[i], color=colors[i], lw=1.8,
                    label=f"bin {i} (thr {t:.0f})")
        ax.set_title(f"{kind}: response functions $R_i(v_{{min}})$")
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        if c == 0:
            ax.set_ylabel(r"$R_i(v_{min})$")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(fontsize=8)

        ax = axes[1, c]
        im = ax.imshow(M, aspect="auto", origin="lower", cmap="magma",
                       extent=[vmid[0], vmid[-1], -0.5, len(thr) - 0.5])
        ax.set_title(f"{kind}: matrix $M[i,j]$")
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        if c == 0:
            ax.set_ylabel("energy bin $i$")
        ax.set_yticks(range(len(thr)))
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Test response matrices fed to the best-fit", fontsize=13)
    fig.tight_layout()
    out = "box_test_matrices.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"saved {out}")


# ---------------------------------------------------------------------------
def main():
    vmid, width = vmin_grid()
    thr = thresholds()
    eta = true_eta(vmid)

    kinds = ["box_flat", "box_window", "tail"]
    plot_matrices(vmid, width, thr, kinds)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), sharey=True)
    print(f"grid: {N_VMIN} cols, v_min {V_LO:.0f}-{V_HI:.0f} km/s  "
          f"(eta(v_hi) = {eta[-1]:.3e}, eta(v_lo) = {eta[0]:.3e})")
    print(f"thresholds (km/s): {np.array2string(thr, precision=0)}\n")

    for ax, kind in zip(axes, kinds):
        M = build_matrix(kind, vmid, width, thr)
        x, res = best_fit(M, eta)

        # last-step report + tie-break diagnosis: the true eta reproduces the
        # data exactly (it IS on the chi2=0 face); the best-fit is chosen over
        # it only because it has smaller total flux. So any undershoot is a
        # tie-break artifact, not a data/fit limitation.
        k = 5
        rel = np.where(eta[-k:] > 0, x[-k:] / eta[-k:], np.nan)
        residual = np.linalg.norm(M @ x - M @ eta) / np.linalg.norm(M @ eta)
        print(f"[{kind}] backend={res.backend}  "
              f"||M x - M eta|| / ||M eta|| = {residual:.1e}  "
              f"(eta is on the chi2=0 face)")
        print(f"    total flux  sum(x)/sum(eta) = {x.sum() / eta.sum():.3f}"
              f"   last {k} steps x/eta = "
              f"{np.array2string(rel, precision=2, floatmode='fixed')}")

        ax.plot(vmid, eta, color="red", lw=2, label="input eta (SHM)")
        ax.step(vmid, x, where="mid", color="C0", lw=1.5, label="best-fit")
        ax.axvline(thr[-1], color="0.6", ls=":", lw=1,
                   label="highest threshold")
        ax.set_xlabel(r"$v_{min}$ [km/s]")
        ax.set_title(f"{kind}\n(backend: {res.backend})")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(fontsize=9)
    axes[0].set_ylabel(r"$\tilde{\eta}$ (natural units)")
    fig.suptitle("Best-fit on box vs tail response matrices "
                 "(noiseless, grid ends mid-eta)", fontsize=13)
    fig.tight_layout()
    out = "box_test.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
