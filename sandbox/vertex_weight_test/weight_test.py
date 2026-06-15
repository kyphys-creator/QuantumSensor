"""Vary the minimal-flux objective weighting (vertex ON) and check eta recovery.

Isolated sandbox (see README). Plain weighted-flux LP in an O(1) unit system:
    minimize sum_j w_j x_j  s.t.  M_phys x = data, x_j >= x_{j+1}, x >= 0.
Compares uniform (production), width, colnorm, invcolnorm weights for the
central-area-cut matrices of TES Al q0 R5 SHM (M1/M2/M3) at a fixed alpha.

    python weight_test.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
CSV_DIR = HERE.parents[1] / "Mathematica/output/TES/response_functions_csv"
GEV_NATIVE = 1e9
ALPHA = 0.50
MASS_LABEL = {"1": "10 MeV", "2": "100 MeV", "3": "1 GeV"}
TARGET = 1000.0


# --- round-trip unit machinery (eta proportional to GeV^1) ------------------
def CM_at(GeV):
    return 1.0 / (1.98e-14 * GeV)


def _KKf(v0, vesc):
    return v0**3 * (-2.0 * np.exp(-vesc**2 / v0**2) * np.pi * (vesc / v0)
                    + np.pi**1.5 * erf(vesc / v0))


def eta_halo_at(GeV, mass_tag, vmin_kms):
    MeV = 1e-3 * GeV
    CM = CM_at(GeV)
    SEC = 1.0 / (6.58e-25 * GeV)
    KPS = (1e5 * CM) / SEC
    RHO_DM = 0.4 * GeV / CM**3
    SIGMA_E = 1e-30 * CM**2
    mchi = {"1": 10, "2": 100, "3": 1000}[mass_tag] * MeV
    V0, VE, VESC = 238.0 * KPS, 250.0 * KPS, 544.0 * KPS
    vm = np.asarray(vmin_kms, float) * KPS
    pref = (RHO_DM * SIGMA_E / mchi) * (V0**2 * np.pi / (2.0 * VE * _KKf(V0, VESC)))
    e_vesc = np.exp(-VESC**2 / V0**2)
    sp = np.sqrt(np.pi) * V0
    inner = np.where(
        vm < VESC - VE,
        -4.0 * e_vesc * VE + sp * (erf((vm + VE) / V0) - erf((vm - VE) / V0)),
        np.where(vm < VESC + VE,
                 -2.0 * e_vesc * (VE + VESC - vm) + sp * (erf(VESC / V0) - erf((vm - VE) / V0)),
                 0.0))
    return pref * inner


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


def weighted_min_flux_lp(M_phys, data, w):
    """minimize sum_j w_j x_j s.t. M_phys x = data, monotone, x >= 0 (plain)."""
    n = M_phys.shape[1]
    rows = np.zeros((n - 1, n))
    for j in range(n - 1):
        rows[j, j + 1] = 1.0; rows[j, j] = -1.0
    res = linprog(c=w, A_ub=rows, b_ub=np.zeros(n - 1),
                  A_eq=M_phys, b_eq=data, bounds=(0, None), method="highs")
    return res.x if res.success else None


def weights_for(M_phys, vk):
    cn = np.linalg.norm(M_phys, axis=0)
    return {
        "uniform":    np.ones(M_phys.shape[1]),
        "width":      trap_weights(vk),
        "colnorm":    cn / cn.max(),
        "invcolnorm": (cn.max() / np.where(cn > 0, cn, cn.max())),
    }


def run_mass(mass_tag, GEV_O1):
    d = np.genfromtxt(CSV_DIR / f"Al_q0M{mass_tag}_R5.csv", delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)

    Rc = np.array([central_area_cut(v, R[i], ALPHA) for i in range(len(R))])
    M0 = (Rc * w[None, :])
    keep = M0.any(axis=0)
    M0k, vk = M0[:, keep], v[keep]
    k = GEV_O1 / GEV_NATIVE
    M = M0k / k
    eta = eta_halo_at(GEV_O1, mass_tag, vk)
    exposure = TARGET / float((M @ eta).max())
    m_phys = exposure * M
    data = m_phys @ eta

    print(f"\n=== M{mass_tag} ({MASS_LABEL[mass_tag]}) ===  alpha={ALPHA}  "
          f"window [{vk.min():.0f},{vk.max():.0f}] km/s")
    print(f"  {'weight':>11} {'last x/eta':>11} {'||x-eta||/||eta||':>18}")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(vk, eta, color="red", lw=2.5, label="input eta", zorder=5)
    colors = plt.cm.viridis(np.linspace(0, 0.85, 4))
    for (name, wj), col in zip(weights_for(m_phys, vk).items(), colors):
        x = weighted_min_flux_lp(m_phys, data, wj)
        if x is None:
            print(f"  {name:>11}   LP failed"); continue
        last = x[-1] / eta[-1]
        err = np.linalg.norm(x - eta) / np.linalg.norm(eta)
        print(f"  {name:>11} {last:11.2f} {err:18.2e}")
        ax.step(vk, x, where="mid", color=col, lw=1.6, label=f"{name}")
    ax.set_xscale("log"); ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"$\tilde{\eta}$ (O(1) units)")
    ax.set_title(f"Minimal-flux objective weighting, M{mass_tag} "
                 f"({MASS_LABEL[mass_tag]}), alpha={ALPHA}")
    ax.grid(True, which="both", ls="--", alpha=0.35); ax.legend(fontsize=9)
    fig.savefig(HERE / f"weights_M{mass_tag}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    # one O(1) GeV for all masses (eta ~ O(1) for M1; M2/M3 a bit larger -- fine)
    eta1 = eta_halo_at(GEV_NATIVE, "1", np.array([200.0]))
    GEV_O1 = GEV_NATIVE / float(np.abs(eta_halo_at(GEV_NATIVE, "1",
              np.genfromtxt(CSV_DIR / "Al_q0M1_R5.csv", delimiter=",", names=True)["vmin"])).max())
    print(f"GeV_O1 = {GEV_O1:.2e}")
    for m in ("1", "2", "3"):
        run_mass(m, GEV_O1)
    print("\nsaved weights_M{1,2,3}.pdf")


if __name__ == "__main__":
    main()
