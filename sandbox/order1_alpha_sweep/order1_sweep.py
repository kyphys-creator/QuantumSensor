"""Central-area alpha sweep done in an O(1) unit system via the round-trip.

Isolated, self-contained sandbox (see README). Recomputes eta at a GeV chosen so
|eta| ~ O(1) (the round-trip function method), rescales the response matrix by
its energy-dimension, and solves the minimal-total-flux best-fit with a PLAIN LP
(no column scaling). TES Al q0 M1 R5 SHM.

    python order1_sweep.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import linprog

# production package (read-only) -- only for the column-scaled cross-check
from quantum_sensor.optimizer import run_optimize_qp
from quantum_sensor.model import condition
from quantum_sensor.config import CONDITION_C

HERE = Path(__file__).resolve().parent
CSV = (HERE.parents[1]
       / "Mathematica/output/TES/response_functions_csv/Al_q0M1_R5.csv")
GEV_NATIVE = 1e9
MASS_TAG = "1"
ALPHAS = [0.30, 0.50, 0.70, 0.90, 0.99]
TARGET_COUNTS = 1000.0


# --- round-trip unit machinery (mirrors constants.py / eta_models.py) -------
def rebuild_constants(GeV):
    CM = 1.0 / (1.98e-14 * GeV)
    SEC = 1.0 / (6.58e-25 * GeV)
    return dict(GeV=GeV, MeV=1e-3 * GeV, CM=CM, KPS=(1e5 * CM) / SEC,
                RHO_DM=0.4 * GeV / CM**3, SIGMA_E=1e-30 * CM**2)


def _KKf(v0, vesc):
    return v0**3 * (-2.0 * np.exp(-vesc**2 / v0**2) * np.pi * (vesc / v0)
                    + np.pi**1.5 * erf(vesc / v0))


def eta_halo_at(GeV, mass_tag, vmin_kms):
    """SHM eta(v_min) recomputed entirely at the given GeV (full round-trip)."""
    c = rebuild_constants(GeV)
    mchi = {"1": 10, "2": 100, "3": 1000}[mass_tag] * c["MeV"]
    V0, VE, VESC = 238.0 * c["KPS"], 250.0 * c["KPS"], 544.0 * c["KPS"]
    vm = np.asarray(vmin_kms, float) * c["KPS"]
    pref = (c["RHO_DM"] * c["SIGMA_E"] / mchi) * (V0**2 * np.pi / (2.0 * VE * _KKf(V0, VESC)))
    e_vesc = np.exp(-VESC**2 / V0**2)
    sp = np.sqrt(np.pi) * V0
    inner = np.where(
        vm < VESC - VE,
        -4.0 * e_vesc * VE + sp * (erf((vm + VE) / V0) - erf((vm - VE) / V0)),
        np.where(vm < VESC + VE,
                 -2.0 * e_vesc * (VE + VESC - vm) + sp * (erf(VESC / V0) - erf((vm - VE) / V0)),
                 0.0))
    return pref * inner


# --- response cut + matrix (self-contained) ---------------------------------
def trap_weights(v):
    w = np.empty_like(v)
    w[1:-1] = 0.5 * (v[2:] - v[:-2]); w[0] = 0.5 * (v[1] - v[0]); w[-1] = 0.5 * (v[-1] - v[-2])
    return w


def central_area_cut(v, r, alpha):
    w = trap_weights(v); area = r * w; total = area.sum()
    if total <= 0:
        return np.zeros_like(r)
    pk = int(np.argmax(r)); lo = hi = pk; acc = area[pk]; n = len(v)
    while acc < alpha * total:
        lv = r[lo - 1] if lo - 1 >= 0 else -np.inf
        rv = r[hi + 1] if hi + 1 < n else -np.inf
        if lv == -np.inf and rv == -np.inf:
            break
        if lv >= rv:
            lo -= 1; acc += area[lo]
        else:
            hi += 1; acc += area[hi]
    out = np.zeros_like(r); out[lo:hi + 1] = r[lo:hi + 1]
    return out


def plain_min_flux_lp(M_phys, data):
    """Minimal-total-flux monotone exact-fit, PLAIN (no column scaling)."""
    n = M_phys.shape[1]
    rows = np.zeros((n - 1, n))
    for j in range(n - 1):
        rows[j, j + 1] = 1.0; rows[j, j] = -1.0
    res = linprog(c=np.ones(n), A_ub=rows, b_ub=np.zeros(n - 1),
                  A_eq=M_phys, b_eq=data, bounds=(0, None), method="highs")
    return (res.x if res.success else None)


def production_fit(M_phys, eta, c):
    """Column-scaled production vertex with conditioning constant ``c``."""
    signal = M_phys @ eta
    mc, dc, bc, unscale = condition(M_phys, signal, np.zeros_like(signal), c=c)
    return unscale(run_optimize_qp(mc, dc, bc, M_phys.shape[1]).x)


def main():
    d = np.genfromtxt(CSV, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)

    # pick GeV so that |eta| ~ O(1) over the full grid (round-trip recompute)
    eta_native = eta_halo_at(GEV_NATIVE, MASS_TAG, v)
    scale = float(np.abs(eta_native).max())
    GEV_O1 = GEV_NATIVE / scale
    eta_o1_full = eta_halo_at(GEV_O1, MASS_TAG, v)         # recomputed at O(1) GeV
    print(f"native |eta| max = {scale:.3e}  ->  GeV_O1 = {GEV_O1:.3e}  "
          f"(|eta| max now {np.abs(eta_o1_full).max():.3f})\n")

    k = GEV_O1 / GEV_NATIVE                                # R_bin proportional to GeV^-1
    runs = []
    for alpha in ALPHAS:
        R_cut = np.array([central_area_cut(v, R[i], alpha) for i in range(len(R))])
        M = (R_cut * w[None, :]) / k                       # matrix at the O(1) GeV
        keep = M.any(axis=0)
        Mk, vk, ek = M[:, keep], v[keep], eta_o1_full[keep]
        exposure = TARGET_COUNTS / float((Mk @ ek).max())
        m_phys = exposure * Mk
        data = m_phys @ ek
        x = plain_min_flux_lp(m_phys, data)
        # production with c matched to the O(1) x-scale vs the native default
        c_match = float(np.abs(ek).max())
        xp = production_fit(m_phys, ek, c_match)
        xp_def = production_fit(m_phys, ek, CONDITION_C)   # native c=1e-30
        if x is None:
            print(f"alpha={alpha:.2f}  PLAIN LP FAILED")
            continue
        last = x[-1] / ek[-1]
        m_match = np.linalg.norm(x - xp) / np.linalg.norm(x)
        m_def = np.linalg.norm(x - xp_def) / np.linalg.norm(x)
        runs.append((alpha, vk, x, ek, xp))
        print(f"alpha={alpha:.2f}  {Mk.shape[1]:3d} cols  window "
              f"[{vk.min():5.0f},{vk.max():5.0f}]  last x/eta = {last:5.2f}   "
              f"prod(c~1) err={m_match:.0e}  prod(c=1e-30) err={m_def:.1f}")

        # per-alpha figure (y in O(1) natural units)
        fig, ax = plt.subplots(figsize=(7, 5))
        vg = np.logspace(0, np.log10(v.max()), 400)
        ax.plot(vg, eta_halo_at(GEV_O1, MASS_TAG, vg), color="red", lw=2,
                label="input eta (SHM, O(1))")
        ax.step(vk, x, where="mid", color="C0", lw=1.6, label="plain LP (no col-scaling)")
        ax.step(vk, xp, where="mid", color="C1", lw=1.0, ls="--",
                label="production (col-scaling, c~O(1))")
        ax.axvspan(vk.min(), vk.max(), color="C0", alpha=0.07)
        ax.set_xscale("log"); ax.set_xlim(v.min(), v.max()); ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"$\tilde{\eta}$ (O(1) units)")
        ax.set_title(rf"$\alpha$={alpha:.2f}  eta~O(1) via round-trip  (TES Al q0 M1 R5)")
        ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=8)
        (HERE / "by_alpha").mkdir(exist_ok=True)
        fig.savefig(HERE / "by_alpha" / f"alpha_{alpha:.2f}.pdf", bbox_inches="tight")
        plt.close(fig)

    # overlay
    fig, ax = plt.subplots(figsize=(9, 6))
    vg = np.logspace(0, np.log10(v.max()), 400)
    ax.plot(vg, eta_halo_at(GEV_O1, MASS_TAG, vg), color="red", lw=2.5,
            label="input eta (SHM, O(1))", zorder=5)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(runs)))
    for (alpha, vk, x, ek, xp), col in zip(runs, colors):
        ax.step(vk, x, where="mid", color=col, lw=1.6,
                label=rf"$\alpha$={alpha:.2f} ({vk.min():.0f}-{vk.max():.0f} km/s)")
    ax.set_xscale("log"); ax.set_xlim(v.min(), v.max()); ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"$\tilde{\eta}$ (O(1) units)")
    ax.set_title("Central-area alpha sweep, plain LP in O(1) units "
                 "(TES Al q0 M1 R5, SHM)")
    ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=9, loc="upper right")
    fig.savefig(HERE / "overlay.pdf", bbox_inches="tight")
    plt.close(fig)
    print("\nsaved by_alpha/alpha_*.pdf and overlay.pdf")


if __name__ == "__main__":
    main()
