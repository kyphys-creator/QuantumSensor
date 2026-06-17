"""Method C (smooth-curvature -> DP segments -> level refit) staircase tie-break.

Reusable, solver-agnostic implementation plus two test regimes:

  PART 1  noiseless exact fit (mu = M @ eta), all 12 configs
          -> #steps, last-step ratio x[-1]/eta[-1], relErr, chi-residual
  PART 2  Poisson noise, representative configs, many seeds
          -> distribution of last-step ratio + relErr for method C vs the
             production colnorm vertex (does method C stay collapse/spike-free
             and unbiased under noise?)

Method C, given fitted counts mu (= y for exact fit, = M @ x_chi2min under
noise) and K = #energy bins:
  1. x_smooth = argmin sum (2nd diff)^2   s.t. M x = mu, x_i>=x_{i+1}, x>=0
  2. segment x_smooth into K monotone piecewise-constant pieces (exact DP)
  3. refit the K level values: argmin ||M x - mu||^2 on those segments
     (monotone, >=0)  -> a genuine K-step staircase reproducing mu.
"""
import numpy as np
import warnings
import cvxpy as cp
from scipy.optimize import linprog

warnings.filterwarnings("ignore")
from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

CONFIGS = [(m, q, ms, 5) for m in ("Al", "TiN") for q in ("0", "2")
           for ms in ("1", "2", "3")]


# ----------------------------- method C -----------------------------
def smooth_curv(M, mu, eps=0.0):
    n = M.shape[1]
    x = cp.Variable(n)
    cons = [M @ x == mu, x[:-1] - x[1:] >= 0, x >= eps]
    cp.Problem(cp.Minimize(cp.sum_squares(x[:-2] - 2 * x[1:-1] + x[2:])),
               cons).solve(solver=cp.CLARABEL)
    return None if x.value is None else np.asarray(x.value)


def dp_segments(target, K):
    """Exact K-segment piecewise-constant L2 partition of `target` (DP)."""
    n = len(target)
    pre = np.concatenate([[0.0], np.cumsum(target)])
    pre2 = np.concatenate([[0.0], np.cumsum(target * target)])

    def cost(i, j):
        c = j - i
        s = pre[j] - pre[i]
        return (pre2[j] - pre2[i]) - s * s / c

    INF = float("inf")
    dp = [[INF] * (n + 1) for _ in range(K + 1)]
    bk = [[0] * (n + 1) for _ in range(K + 1)]
    dp[0][0] = 0.0
    for k in range(1, K + 1):
        for i in range(1, n + 1):
            for j in range(k - 1, i):
                if dp[k - 1][j] == INF:
                    continue
                c = dp[k - 1][j] + cost(j, i)
                if c < dp[k][i]:
                    dp[k][i], bk[k][i] = c, j
    bnds, i = [], n
    for k in range(K, 0, -1):
        j = bk[k][i]
        bnds.append((j, i)); i = j
    return bnds[::-1]


def refit_levels(M, mu, bounds, eps=0.0):
    K, n = len(bounds), M.shape[1]
    G = np.zeros((M.shape[0], K))
    for k, (a, b) in enumerate(bounds):
        G[:, k] = M[:, a:b].sum(axis=1)
    lv = cp.Variable(K)
    cp.Problem(cp.Minimize(cp.sum_squares(G @ lv - mu)),
               [lv[:-1] - lv[1:] >= 0, lv >= eps]).solve(solver=cp.CLARABEL)
    levels = np.asarray(lv.value)
    x = np.empty(n)
    for k, (a, b) in enumerate(bounds):
        x[a:b] = levels[k]
    return x


def method_c(M, mu, K, eps=0.0):
    xs = smooth_curv(M, mu, eps)
    if xs is None:
        return None
    return refit_levels(M, mu, dp_segments(xs, K), eps)


# --------------------- baselines / helpers ---------------------
def colnorm_vertex(M, mu, eps=0.0):
    n = M.shape[1]
    w = np.sqrt((M * M).sum(axis=0))
    A = np.zeros((n - 1, n))
    for i in range(n - 1):
        A[i, i], A[i, i + 1] = 1.0, -1.0
    r = linprog(w, A_ub=-A, b_ub=np.zeros(n - 1), A_eq=M, b_eq=mu,
                bounds=[(eps, None)] * n, method="highs-ds")
    return np.asarray(r.x) if r.success else None


def chi2min_interior(M, data):
    """Neyman chi^2 interior fit (monotone, >=0) -> fitted counts mu = M x*."""
    n = M.shape[1]
    x = cp.Variable(n)
    mask = data > 0
    res = cp.multiply((data[mask] - M[mask] @ x), 1.0 / np.sqrt(data[mask]))
    cp.Problem(cp.Minimize(cp.sum_squares(res)),
               [x[:-1] - x[1:] >= 0, x >= 0]).solve(solver=cp.CLARABEL)
    return None if x.value is None else M @ np.asarray(x.value)


def n_steps(x, rtol=1e-4):
    s = max(abs(x).max(), 1e-300)
    return int((np.abs(np.diff(x)) > rtol * s).sum()) + 1


# ----------------------------- PART 1 -----------------------------
def part1():
    print("=" * 78)
    print("PART 1  noiseless exact fit (mu = M @ eta), all 12 configs")
    print(f"{'config':>12} {'steps':>6} {'last':>7} {'relErr':>9} {'residual':>10}")
    print("-" * 78)
    for mat, q, ms, nb in CONFIGS:
        a = DarkMatterQuantumAnalysis(RunConfig(mat, q, ms, nb, "Halo", background="none"))
        M, eta = a.m_phys, a.eta
        y = M @ eta
        x = method_c(M, y, M.shape[0])
        if x is None:
            print(f"{mat+' q'+q+'M'+ms:>12}   solve failed"); continue
        last = x[-1] / eta[-1]
        re = np.linalg.norm(x - eta) / np.linalg.norm(eta)
        resid = np.linalg.norm(M @ x - y) / np.linalg.norm(y)
        print(f"{mat+' q'+q+'M'+ms:>12} {n_steps(x):6d} {last:7.2f} {re:9.2e} {resid:10.2e}")


# ----------------------------- PART 2 -----------------------------
def part2(seeds=200, target=200.0):
    cfgs = [("Al", "2", "3", 5), ("TiN", "2", "2", 5), ("TiN", "0", "3", 5)]
    print("\n" + "=" * 78)
    print(f"PART 2  Poisson noise ({seeds} seeds, ~{target:.0f} peak counts)")
    print(f"        median [16-84%] of last-step ratio and relErr")
    print(f"{'config':>12} {'method':>9} {'last (med [16,84])':>26} {'relErr (med [16,84])':>26}")
    print("-" * 78)
    rng = np.random.default_rng(0)
    for mat, q, ms, nb in cfgs:
        a = DarkMatterQuantumAnalysis(RunConfig(mat, q, ms, nb, "Halo", background="none"))
        M, eta = a.m_phys, a.eta
        mu0 = M @ eta
        scale = target / mu0.max()
        Ms = M * scale
        mu0s = Ms @ eta
        rec = {"methodC": ([], []), "colnorm": ([], [])}
        for s in range(seeds):
            data = rng.poisson(mu0s).astype(float)
            mu = chi2min_interior(Ms, data)
            if mu is None:
                continue
            for name, fn in (("methodC", method_c), ("colnorm", colnorm_vertex)):
                x = fn(Ms, mu, Ms.shape[0]) if name == "methodC" else fn(Ms, mu)
                if x is None:
                    continue
                rec[name][0].append(x[-1] / eta[-1])
                rec[name][1].append(np.linalg.norm(x - eta) / np.linalg.norm(eta))
        for name in ("methodC", "colnorm"):
            la = np.array(rec[name][0]); er = np.array(rec[name][1])
            lq = np.percentile(la, [16, 50, 84]); eq = np.percentile(er, [16, 50, 84])
            print(f"{mat+' q'+q+'M'+ms:>12} {name:>9} "
                  f"{lq[1]:7.2f} [{lq[0]:5.2f},{lq[2]:5.2f}]      "
                  f"{eq[1]:7.2e} [{eq[0]:.2e},{eq[2]:.2e}]")
        print()


if __name__ == "__main__":
    part1()
    part2()
