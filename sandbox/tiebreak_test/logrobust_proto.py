"""Prototype: log-scale-robust method C smooth step (RELATIVE curvature).

The Bound eta (M3) is a ~1e14-dynamic-range cliff; plain ||D2 x||^2 curvature
blows up at the cliff and the QP fails. Relative curvature ||(D2 x)/s||^2
(s = local magnitude) divides the cliff's huge second difference by its huge
magnitude -> O(1), so it stays conditioned. For a gentle eta s ~ flat and it
reduces to the plain curvature, so the working configs must be preserved.

We need a magnitude profile s that (a) is robust to the dynamic range and (b)
does NOT itself collapse at the edge (or it would drag x down). Test choices:
  s_lp   : colnorm minimal-flux LP vertex (robust, but collapses at the edge)
  s_env  : running-max envelope of s_lp from the high-v end (non-increasing,
           never below a fraction of the local level) -> kills the edge collapse
Compared against plain linear curvature (baseline) on Bound + Halo.
"""
import numpy as np
import warnings
import cvxpy as cp
from scipy.optimize import linprog

warnings.filterwarnings("ignore")
from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.optimizer import _build_qp, _dp_segments, _refit_levels


def colnorm_lp(M, mu):
    n = M.shape[1]
    w = np.sqrt((M * M).sum(0))
    A = np.zeros((n - 1, n))
    for i in range(n - 1):
        A[i, i], A[i, i + 1] = 1.0, -1.0
    r = linprog(w, A_ub=-A, b_ub=np.zeros(n - 1), A_eq=M, b_eq=mu,
                bounds=[(0, None)] * n, method="highs-ds")
    return np.asarray(r.x) if r.success else None


def smooth(M, mu, s=None):
    """min ||(D2 x)/s_mid||^2 s.t. M x = mu, monotone, x>=0. s=None -> plain."""
    n = M.shape[1]
    x = cp.Variable(n)
    d2 = x[:-2] - 2 * x[1:-1] + x[2:]
    if s is None:
        obj = cp.sum_squares(d2)
    else:
        obj = cp.sum_squares(cp.multiply(1.0 / s[1:-1], d2))
    cons = [M @ x == mu, x[:-1] - x[1:] >= 0, x >= 0]
    cp.Problem(cp.Minimize(obj), cons).solve(solver=cp.CLARABEL)
    return None if x.value is None else np.asarray(x.value)


def smooth_relcurv_wsub(M, mu, D):
    """Conditioned RELATIVE-curvature smooth step.

    Objective is relative curvature ||(D2 x)/D_mid||^2 (s = D as divisor): where
    D is small the curvature is penalised hard -> x forced flat (NOT collapsed),
    so it is anti-collapse and reduces to plain curvature for a flat D (gentle
    configs preserved). The variable is rescaled x = D (.) w purely for
    numerical conditioning (exact change of variables, objective unchanged), so
    the 1e14 cliff no longer breaks the solver; the equality is row-scaled by
    |mu|. D is a non-collapsing magnitude envelope of a robust LP vertex."""
    n = M.shape[1]
    w = cp.Variable(n)
    MD = M * D[None, :]
    r = np.where(np.abs(mu) > 0, np.abs(mu), 1.0)        # equality row scale
    x = cp.multiply(D, w)
    d2x = x[:-2] - 2 * x[1:-1] + x[2:]
    obj = cp.sum_squares(cp.multiply(1.0 / D[1:-1], d2x))    # relative curvature
    cons = [cp.multiply(1.0 / r, MD @ w) == mu / r,
            D[:-1] * w[:-1] - D[1:] * w[1:] >= 0,
            w >= 0]
    cp.Problem(cp.Minimize(obj), cons).solve(solver=cp.CLARABEL)
    return None if w.value is None else D * np.asarray(w.value)


def refit_cond(M, mu, bounds, xs):
    """Level refit conditioned by the smooth solution's per-segment magnitude:
    levels l = L (.) lt with L_k = mean(xs in segment k), lt ~ O(1)."""
    K = len(bounds)
    G = np.zeros((M.shape[0], K))
    L = np.zeros(K)
    for k, (a, b) in enumerate(bounds):
        G[:, k] = M[:, a:b].sum(1)
        L[k] = max(xs[a:b].mean(), 1e-300)
    lt = cp.Variable(K)
    r = np.where(np.abs(mu) > 0, np.abs(mu), 1.0)
    lev = cp.multiply(L, lt)
    cons = [L[:-1] * lt[:-1] - L[1:] * lt[1:] >= 0, lt >= 0]
    cp.Problem(cp.Minimize(cp.sum_squares(cp.multiply(1.0 / r, G @ lev - mu))),
               cons).solve(solver=cp.CLARABEL)
    if lt.value is None:
        return None
    levels = L * np.asarray(lt.value)
    x = np.empty(M.shape[1])
    for k, (a, b) in enumerate(bounds):
        x[a:b] = levels[k]
    return x


def hold_floor(s, floor_frac=1e-3):
    """Non-increasing magnitude from the LP vertex, holding the last positive
    plateau through any edge collapse. Sweep low->high: a genuine decrease
    (incl. a cliff step, ~7x) updates the running level; an edge collapse to ~0
    (drop below floor_frac * level) holds the previous level. So a flat LP stays
    flat (gentle eta preserved), and a cliff is followed but the high-v edge is
    floored at the last resolved level instead of collapsing to ~0."""
    s = np.maximum(np.asarray(s, float), 0.0)
    D = np.empty_like(s)
    level = max(s[0], 1e-300)
    for j in range(len(s)):
        if s[j] >= floor_frac * level:
            level = s[j]
        D[j] = level
    return np.maximum(D, 1e-300)


def full_method_c(M, mu, K, xs, cond_refit=False):
    if xs is None:
        return None, "smooth failed"
    bounds = _dp_segments(xs, K)
    x = refit_cond(M, mu, bounds, xs) if cond_refit else \
        _refit_levels(M, mu, bounds, M.shape[1], 0.0)
    return (x, "ok") if x is not None else (None, "refit failed")


CONFIGS = [("Al", "0", "3", 5), ("TiN", "0", "3", 5),
           ("Al", "0", "1", 5), ("TiN", "0", "3", 5)]


def run():
    for mat, q, ms, nb in [("Al", "0", "3", 5), ("TiN", "0", "3", 5)]:
        for et in ("Halo", "Bound"):
            a = DarkMatterQuantumAnalysis(RunConfig(mat, q, ms, nb, et, background="none"))
            M, eta = a.m_phys, a.eta
            mu = M @ eta
            qp = _build_qp(a.m_phys, a.observed, np.zeros_like(a.observed), a.n_vmin, 0.0)
            K = qp["M_active"].shape[0]
            slp = colnorm_lp(M, mu)
            senv = hold_floor(slp) if slp is not None else None
            print(f"\n=== {mat} q{q}M{ms} {et}  (eta range {eta.min():.2g}..{eta.max():.2g}) ===")
            variants = [("plain", smooth(M, mu, None), False),
                        ("relcurv_wsub", smooth_relcurv_wsub(M, mu, senv) if senv is not None else None, True)]
            for tag, xs, cond in variants:
                x, msg = full_method_c(M, mu, K, xs, cond_refit=cond)
                if x is None:
                    print(f"  {tag:14} {msg}"); continue
                last = x[-1] / eta[-1]
                re = np.linalg.norm(x - eta) / np.linalg.norm(eta)
                print(f"  {tag:14} last={last:10.2f}  relErr={re:.2e}")


if __name__ == "__main__":
    run()
