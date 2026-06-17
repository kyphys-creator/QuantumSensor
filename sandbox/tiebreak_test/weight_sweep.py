"""What does changing the vertex-LP WEIGHTING do? (on the current 20%-peak matrices)

The staircase best-fit is the vertex of the chi^2=0 face
   { x : M x = y,  x_i >= x_{i+1},  x >= 0 }
that minimises a LINEAR objective  min sum_j w_j x_j. Production uses
w_j = ||M[:,j]|| (colnorm). This sweeps several weight profiles over all 12
configs and reports the last-step ratio x[-1]/eta[-1] (1 = no collapse/spike)
and the relative L2 error ||x-eta||/||eta||.

    uniform   w_j = 1                  (original minimal total flux)
    colnorm   w_j = ||M[:,j]||         (production)
    colnorm2  w_j = ||M[:,j]||^2
    colsqrt   w_j = sqrt(||M[:,j]||)
    colsum    w_j = sum_i M[i,j]
    invcol    w_j = 1 / ||M[:,j]||     (penalise low-response columns)
"""
import numpy as np
import warnings
from scipy.optimize import linprog

warnings.filterwarnings("ignore")
from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

CONFIGS = [(m, q, ms, 5) for m in ("Al", "TiN") for q in ("0", "2")
           for ms in ("1", "2", "3")]
WEIGHTS = ["uniform", "colnorm", "colnorm2", "colsqrt", "colsum", "invcol"]


def weights_for(M):
    cn = np.sqrt((M * M).sum(axis=0))
    cn = np.where(cn > 0, cn, cn[cn > 0].min() if (cn > 0).any() else 1.0)
    return {
        "uniform":  np.ones(M.shape[1]),
        "colnorm":  cn,
        "colnorm2": cn ** 2,
        "colsqrt":  np.sqrt(cn),
        "colsum":   M.sum(axis=0),
        "invcol":   1.0 / cn,
    }


def vertex(M, y, w):
    n = M.shape[1]
    A = np.zeros((n - 1, n))
    for i in range(n - 1):
        A[i, i], A[i, i + 1] = 1.0, -1.0
    r = linprog(w, A_ub=-A, b_ub=np.zeros(n - 1), A_eq=M, b_eq=y,
                bounds=[(0.0, None)] * n, method="highs-ds",
                options={"presolve": True})
    return np.asarray(r.x) if r.success else None


def main():
    print("LAST-STEP RATIO  x[-1]/eta[-1]   (1.00 = perfect, 0 = collapse, >>1 = spike)")
    hdr = f"{'config':>12} " + " ".join(f"{w:>9}" for w in WEIGHTS)
    print(hdr); print("-" * len(hdr))
    errs = {}
    for mat, q, ms, nb in CONFIGS:
        a = DarkMatterQuantumAnalysis(RunConfig(mat, q, ms, nb, "Halo", background="none"))
        M, eta = a.m_phys, a.eta
        y = M @ eta
        ws = weights_for(M)
        row, erow = [], []
        for name in WEIGHTS:
            x = vertex(M, y, ws[name])
            if x is None:
                row.append(float("nan")); erow.append(float("nan")); continue
            row.append(x[-1] / eta[-1])
            erow.append(np.linalg.norm(x - eta) / np.linalg.norm(eta))
        errs[f"{mat} q{q}M{ms}"] = erow
        print(f"{mat+' q'+q+'M'+ms:>12} " + " ".join(f"{v:9.2f}" for v in row))

    print("\nRELATIVE L2 ERROR  ||x-eta||/||eta||")
    print(hdr); print("-" * len(hdr))
    for k, erow in errs.items():
        print(f"{k:>12} " + " ".join(f"{v:9.2e}" for v in erow))


if __name__ == "__main__":
    main()
