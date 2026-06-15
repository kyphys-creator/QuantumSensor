"""Robustness of the colnorm weighting vs uniform across detector/mass/alpha.

Isolated sandbox. For TES(Al) and MKID(TiN), q0, M1/M2/M3, and a sweep of the
central-area alpha, compare the minimal-flux vertex under the production
`uniform` weight vs the `colnorm` weight (w_j = ||M[:,j]||). Reports last-step
x/eta and overall ||x-eta||/||eta|| for both, to confirm colnorm never does
worse and removes the last-step shave where uniform has it.

    python colnorm_sweep.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from weight_test import (central_area_cut, trap_weights, eta_halo_at,
                         weighted_min_flux_lp, GEV_NATIVE, MASS_LABEL, TARGET)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DETECTORS = {
    "TES":  ("Al",  ROOT / "Mathematica/output/TES/response_functions_csv"),
    "MKID": ("TiN", ROOT / "Mathematica/output/MKID/response_functions_csv"),
}
MASSES = ("1", "2", "3")
ALPHAS = [0.30, 0.50, 0.70, 0.90]


def solve_one(csv, mass_tag, alpha, GEV_O1):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    v = d["vmin"]
    bins = [c for c in d.dtype.names if c != "vmin"]
    R = np.array([d[c] for c in bins])
    w = trap_weights(v)
    Rc = np.array([central_area_cut(v, R[i], alpha) for i in range(len(R))])
    M0 = (Rc * w[None, :])
    keep = M0.any(axis=0)
    M0k, vk = M0[:, keep], v[keep]
    M = M0k / (GEV_O1 / GEV_NATIVE)
    eta = eta_halo_at(GEV_O1, mass_tag, vk)
    m_phys = (TARGET / float((M @ eta).max())) * M
    data = m_phys @ eta
    cn = np.linalg.norm(m_phys, axis=0)
    out = {}
    for name, wj in (("uniform", np.ones(m_phys.shape[1])),
                     ("colnorm", cn / cn.max())):
        x = weighted_min_flux_lp(m_phys, data, wj)
        if x is None:
            out[name] = (np.nan, np.nan)
        else:
            out[name] = (x[-1] / eta[-1],
                         np.linalg.norm(x - eta) / np.linalg.norm(eta))
    return out


def main():
    GEV_O1 = GEV_NATIVE / float(np.abs(eta_halo_at(
        GEV_NATIVE, "1",
        np.genfromtxt(DETECTORS["TES"][1] / "Al_q0M1_R5.csv",
                      delimiter=",", names=True)["vmin"])).max())

    for det, (mat, csvdir) in DETECTORS.items():
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharey=True)
        print(f"\n=== {det} ({mat}) q0 R5 ===")
        print(f"  {'mass':>5} {'alpha':>6} "
              f"{'uniform last':>13} {'colnorm last':>13} "
              f"{'unif err':>10} {'coln err':>10}")
        for ax, mass in zip(axes, MASSES):
            csv = csvdir / f"{mat}_q0M{mass}_R5.csv"
            lu, lc = [], []
            for alpha in ALPHAS:
                r = solve_one(csv, mass, alpha, GEV_O1)
                lu.append(r["uniform"][0]); lc.append(r["colnorm"][0])
                print(f"  M{mass:>4} {alpha:6.2f} "
                      f"{r['uniform'][0]:13.2f} {r['colnorm'][0]:13.2f} "
                      f"{r['uniform'][1]:10.2e} {r['colnorm'][1]:10.2e}")
            ax.plot(ALPHAS, lu, "o-", color="C1", label="uniform")
            ax.plot(ALPHAS, lc, "s-", color="C0", label="colnorm")
            ax.axhline(1.0, color="red", ls="--", lw=1, label="perfect")
            ax.set_xlabel(r"$\alpha$"); ax.set_title(f"M{mass} ({MASS_LABEL[mass]})")
            ax.set_ylim(-0.1, 1.3); ax.grid(True, ls="--", alpha=0.35); ax.legend(fontsize=8)
        axes[0].set_ylabel(r"last-step $x/\eta$")
        fig.suptitle(f"colnorm vs uniform last-step recovery -- {det} ({mat}) q0 R5 SHM",
                     fontsize=13)
        fig.tight_layout()
        fig.savefig(HERE / f"colnorm_sweep_{det}.pdf", bbox_inches="tight")
        plt.close(fig)
    print("\nsaved colnorm_sweep_{TES,MKID}.pdf")


if __name__ == "__main__":
    main()
