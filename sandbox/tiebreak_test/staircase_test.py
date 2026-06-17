"""Find a STAIRCASE (piecewise-constant vertex) tie-break that does not
collapse/spike.

Hard requirement: the recovered flux MUST be a staircase (a vertex of the
monotone polytope), as in neutrinoAnalysis. The smooth min-curvature solution
is therefore disqualified (it is the interior of the chi^2=0 face, not a
vertex). Note the smooth solution is itself feasible (monotone, non-negative,
M x = y), so projecting onto the polytope just returns it -- only a LINEAR
objective min sum_j w_j x_j is guaranteed to land on a vertex.

Two routes are compared, both on the exact-fit face (mu = y = M @ eta):

  (A) LINEAR-WEIGHT VERTICES  min sum_j w_j x_j   (LP -> guaranteed vertex)
        colnorm    w_j = ||M[:,j]||            (current production)
        colnorm2   w_j = ||M[:,j]||^2
        colsqrt    w_j = ||M[:,j]||^0.5
        colsum     w_j = sum_i M[i,j]
      Each is a true staircase; we look for one that tracks eta to the edge.

  (B) PC-APPROX OF SMOOTH  : fit the smooth min-curvature curve x_smooth with a
      K = n_ebins level monotone piecewise-constant staircase (exact DP), then
      report the chi^2 residual ||M x - y|| / ||y||. This GUARANTEES a clean
      K-step staircase that tracks eta; the residual says how far off the
      exact fit it is.

Reported per config: number of distinct plateaus (staircase check), last-step
ratio x[-1]/eta[-1], relative L2 error to eta, and (B only) the fit residual.
"""
import os

import numpy as np
import cvxpy as cp
from scipy.optimize import linprog
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

HERE = os.path.dirname(os.path.abspath(__file__))

CONFIGS = [
    ("Al q2 M3  (collapse)", "Al", "2", "3", 5),
    ("TiN q2 M2 (collapse)", "TiN", "2", "2", 5),
    ("TiN q0 M3 (spike)",    "TiN", "0", "3", 5),
    ("Al q0 M1  (control)",  "Al", "0", "1", 5),
    ("TiN q0 M1 (control)",  "TiN", "0", "1", 5),
]


def vertex_lp(M, y, w, eps):
    """min w . x  s.t.  M x = y,  x_i >= x_{i+1},  x >= eps   (LP -> vertex)."""
    n = M.shape[1]
    A_ord = np.zeros((n - 1, n))
    for i in range(n - 1):
        A_ord[i, i] = 1.0
        A_ord[i, i + 1] = -1.0
    res = linprog(w, A_ub=-A_ord, b_ub=np.zeros(n - 1),
                  A_eq=M, b_eq=y, bounds=[(eps, None)] * n,
                  method="highs-ds", options={"presolve": True})
    return np.asarray(res.x) if res.success else None


def smooth_curv(M, y, eps):
    """Min second-difference^2 over the face -> the smooth reference curve."""
    n = M.shape[1]
    x = cp.Variable(n)
    cons = [M @ x == y, x[:-1] - x[1:] >= 0, x >= eps]
    obj = cp.Minimize(cp.sum_squares(x[:-2] - 2 * x[1:-1] + x[2:]))
    cp.Problem(obj, cons).solve(solver=cp.CLARABEL)
    return np.asarray(x.value) if x.value is not None else None


def pc_fit_monotone(target, K):
    """Exact K-segment piecewise-constant L2 fit of `target` (DP).

    target is (near) monotone decreasing, so segment means are automatically
    non-increasing -> the result is a valid monotone staircase. Returns the
    fitted piecewise-constant vector.
    """
    n = len(target)
    pre = np.concatenate([[0.0], np.cumsum(target)])
    pre2 = np.concatenate([[0.0], np.cumsum(target * target)])

    def seg_cost(i, j):  # SSE of target[i:j] about its mean
        cnt = j - i
        s = pre[j] - pre[i]
        s2 = pre2[j] - pre2[i]
        return s2 - s * s / cnt

    INF = float("inf")
    # dp[k][i] = min cost to fit target[:i] with k segments
    dp = [[INF] * (n + 1) for _ in range(K + 1)]
    bk = [[0] * (n + 1) for _ in range(K + 1)]
    dp[0][0] = 0.0
    for k in range(1, K + 1):
        for i in range(1, n + 1):
            for j in range(k - 1, i):
                if dp[k - 1][j] == INF:
                    continue
                c = dp[k - 1][j] + seg_cost(j, i)
                if c < dp[k][i]:
                    dp[k][i] = c
                    bk[k][i] = j
    # backtrack segment boundaries
    bounds = []
    i = n
    for k in range(K, 0, -1):
        j = bk[k][i]
        bounds.append((j, i))
        i = j
    bounds.reverse()
    out = np.empty(n)
    for (a, b) in bounds:
        out[a:b] = target[a:b].mean()
    return out


def segment_bounds(target, K):
    """DP breakpoints: list of (a, b) segments for K-piecewise-constant fit."""
    n = len(target)
    pre = np.concatenate([[0.0], np.cumsum(target)])
    pre2 = np.concatenate([[0.0], np.cumsum(target * target)])

    def seg_cost(i, j):
        cnt = j - i
        s = pre[j] - pre[i]
        s2 = pre2[j] - pre2[i]
        return s2 - s * s / cnt

    INF = float("inf")
    dp = [[INF] * (n + 1) for _ in range(K + 1)]
    bk = [[0] * (n + 1) for _ in range(K + 1)]
    dp[0][0] = 0.0
    for k in range(1, K + 1):
        for i in range(1, n + 1):
            for j in range(k - 1, i):
                if dp[k - 1][j] == INF:
                    continue
                c = dp[k - 1][j] + seg_cost(j, i)
                if c < dp[k][i]:
                    dp[k][i] = c
                    bk[k][i] = j
    bounds = []
    i = n
    for k in range(K, 0, -1):
        j = bk[k][i]
        bounds.append((j, i))
        i = j
    bounds.reverse()
    return bounds


def refit_levels(M, y, bounds, eps=0.0):
    """Given fixed plateau segments, re-fit the level VALUES to the data by
    chi^2 (monotone, non-negative) -> an exact/near-exact-fit staircase that is
    a genuine vertex with the smooth curve's plateau structure."""
    K = len(bounds)
    n = M.shape[1]
    # group matrix: column k sums the M columns inside segment k
    G = np.zeros((M.shape[0], K))
    for k, (a, b) in enumerate(bounds):
        G[:, k] = M[:, a:b].sum(axis=1)
    lv = cp.Variable(K)
    cons = [lv[:-1] - lv[1:] >= 0, lv >= eps]
    cp.Problem(cp.Minimize(cp.sum_squares(G @ lv - y)), cons).solve(solver=cp.CLARABEL)
    levels = np.asarray(lv.value)
    x = np.empty(n)
    for k, (a, b) in enumerate(bounds):
        x[a:b] = levels[k]
    return x


def n_plateaus(x, rtol=1e-4):
    """Count distinct constant levels (staircase steps)."""
    scale = max(abs(x).max(), 1e-300)
    jumps = np.abs(np.diff(x)) > rtol * scale
    return int(jumps.sum()) + 1


WEIGHTS = ["colnorm", "colnorm2", "colsqrt", "colsum"]


def weights_for(M):
    cn = np.sqrt((M * M).sum(axis=0))
    return {
        "colnorm":  cn,
        "colnorm2": cn ** 2,
        "colsqrt":  np.sqrt(cn),
        "colsum":   M.sum(axis=0),
    }


def run():
    fig, axes = plt.subplots(len(CONFIGS), 1, figsize=(8, 3.0 * len(CONFIGS)))
    if len(CONFIGS) == 1:
        axes = [axes]

    cols = WEIGHTS + ["PCfit(B)"]
    print("\nSTAIRCASE steps (#plateaus) | last-ratio x[-1]/eta[-1] | relErr | PCfit residual")
    for ax, (label, mat, q, mass, nb) in zip(axes, CONFIGS):
        a = DarkMatterQuantumAnalysis(
            RunConfig(mat, q, mass, nb, "Halo", background="none"))
        M, eta, v = a.m_phys, a.eta, a.vmin_mid
        y = M @ eta
        eps = 0.0
        ny = np.linalg.norm(y)

        ax.plot(v, eta, "k-", lw=2.2, label="eta", zorder=10)
        print(f"\n{label}   (n_vmin={len(v)}, n_ebins={M.shape[0]})")

        ws = weights_for(M)
        for name in WEIGHTS:
            x = vertex_lp(M, y, ws[name], eps)
            if x is None:
                print(f"  {name:9} LP failed")
                continue
            npl = n_plateaus(x)
            lr = x[-1] / eta[-1]
            re = np.linalg.norm(x - eta) / np.linalg.norm(eta)
            print(f"  {name:9} steps={npl:3d}  last={lr:6.2f}  relErr={re:.2e}")
            ax.step(v, x, where="post", lw=1.0, alpha=0.8, label=name)

        # (B) PC approx of smooth curve
        xs = smooth_curv(M, y, eps)
        if xs is not None:
            bounds = segment_bounds(xs, M.shape[0])
            xpc = pc_fit_monotone(xs, M.shape[0])
            npl = n_plateaus(xpc)
            lr = xpc[-1] / eta[-1]
            re = np.linalg.norm(xpc - eta) / np.linalg.norm(eta)
            resid = np.linalg.norm(M @ xpc - y) / ny
            print(f"  PCfit(B)  steps={npl:3d}  last={lr:6.2f}  relErr={re:.2e}  "
                  f"chi-residual={resid:.2e}")
            ax.step(v, xpc, where="post", lw=1.4, color="tab:orange",
                    label="PCfit smooth (B)")

            # (C) refit the level VALUES to data on the smooth-curve segments
            xrf = refit_levels(M, y, bounds, eps)
            npl = n_plateaus(xrf)
            lr = xrf[-1] / eta[-1]
            re = np.linalg.norm(xrf - eta) / np.linalg.norm(eta)
            resid = np.linalg.norm(M @ xrf - y) / ny
            print(f"  refit(C)  steps={npl:3d}  last={lr:6.2f}  relErr={re:.2e}  "
                  f"chi-residual={resid:.2e}")
            ax.step(v, xrf, where="post", lw=2.2, color="tab:red",
                    label="refit levels (C)")

        ax.set_title(label)
        ax.set_xlabel("v_min [km/s]")
        ax.set_ylabel("flux")
        ax.set_yscale("log")
        ax.legend(fontsize=7, ncol=3)

    fig.tight_layout()
    out = os.path.join(HERE, "staircase_compare.pdf")
    fig.savefig(out)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    run()
