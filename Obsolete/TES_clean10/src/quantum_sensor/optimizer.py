import numpy as np
from numba import njit
from scipy.optimize import minimize, linprog, LinearConstraint
from scipy.sparse import csc_matrix, eye as sp_eye, vstack as sp_vstack


@njit
def chi_neyman(xvec, M, data, Bkg, eps):
    """Neyman chi-squared: sum((data - Bkg - M@x)^2 / data) over bins with data > 0."""
    mod = M.dot(xvec)
    mask = data > 0
    return np.sum(((data[mask] - Bkg[mask] - mod[mask]) ** 2) / data[mask])


def ordering_constraint(n: int) -> LinearConstraint:
    """x[i] >= x[i+1] for all i — non-increasing flux."""
    A = np.zeros((n - 1, n))
    for i in range(n - 1):
        A[i, i] = 1
        A[i, i + 1] = -1
    return LinearConstraint(csc_matrix(A), lb=0, ub=np.inf)


def _ordering_csc(n: int):
    A = np.zeros((n - 1, n))
    for i in range(n - 1):
        A[i, i] = 1
        A[i, i + 1] = -1
    return csc_matrix(A)


def run_optimize(M_matrix, data, Bkg_vector, n, eps=0.0, x0=None, display=True):
    if x0 is None:
        x0 = np.ones(n)

    bounds = [(eps, None)] * n
    constraints = [ordering_constraint(n)]
    options = {'maxiter': 100000, 'xtol': 1e-9, 'gtol': 1e-9}
    if display:
        options['verbose'] = 3

    return minimize(
        lambda x: chi_neyman(x, M_matrix, data, Bkg_vector, eps),
        x0,
        method='trust-constr',
        bounds=bounds,
        constraints=constraints,
        options=options,
        jac='3-point',
    )


def _trim_bkg(Bkg_vector, m):
    if len(Bkg_vector) >= m:
        return Bkg_vector[:m]
    return np.pad(Bkg_vector, (0, m - len(Bkg_vector)))


def _build_qp(M_matrix, data, Bkg_vector, n, eps):
    """Build χ² as a QP: (1/2) xᵀ P x + qᵀ x + const, with mask on data>0."""
    m = len(data)
    Bkg = _trim_bkg(Bkg_vector, m)
    y = data - Bkg
    mask = data > 0
    M = M_matrix[mask]
    d = data[mask]
    yv = y[mask]
    inv_d = 1.0 / d

    P = 2.0 * (M.T * inv_d) @ M
    q = -2.0 * M.T @ (yv * inv_d)
    const = float((yv * yv * inv_d).sum())

    A_ord = _ordering_csc(n)
    I = sp_eye(n, format='csc')

    return {
        'P': P, 'q': q, 'const': const,
        'M_active': M, 'y_active': yv, 'inv_d': inv_d,
        'A_ord': A_ord, 'I': I, 'n': n, 'm_active': int(mask.sum()),
        'eps': eps, 'mask': mask,
    }


class _Result:
    __slots__ = ('x', 'fun', 'success', 'message', 'nit', 'solve_time', 'backend', 'staircase')


def _staircase_vertex(qp, x_ref=None):
    """HiGHS-simplex LP: pick a vertex on the optimal face.

    For underdetermined problems, the QP optimum is the affine set
        {x : M_active x = y_active, A_ord x >= 0, x >= eps}
    (χ²_min = 0). We solve a feasibility LP on it; the simplex returns a
    basic vertex with at most m_active distinct non-degenerate coordinates
    → staircase with ≤ m_active + 1 plateaus.

    Objective: minimize Σ x (least-flux staircase that matches the data).

    `x_ref` is accepted for API compatibility; the equality target is the
    true y_active so the residual is machine-precision regardless.
    """
    n = qp['n']
    M = qp['M_active']
    y_eq = qp['y_active']

    # HiGHS rejects ill-scaled constraints; normalize each equality row.
    row_scale = np.abs(M).max(axis=1)
    row_scale[row_scale == 0] = 1.0
    M_n = M / row_scale[:, None]
    y_n = y_eq / row_scale

    A_ord_dense = -qp['A_ord'].toarray()  # ordering: -(x_i - x_{i+1}) <= 0
    b_ub = np.zeros(n - 1)
    bounds = [(qp['eps'], None)] * n
    c = np.ones(n)

    res = linprog(
        c, A_ub=A_ord_dense, b_ub=b_ub,
        A_eq=M_n, b_eq=y_n,
        bounds=bounds, method='highs-ds',
        options={'presolve': True},
    )
    if not res.success:
        return None
    return np.asarray(res.x)


class _OSQPBackend:
    """OSQP-based QP solver with optional simplex vertex selection.

    For underdetermined problems (n > m_active), the QP minimu m is a face of
    the polytope rather than a point. `vertex_select=True` (default) reroutes
    to a HiGHS dual-simplex LP that picks a basic feasible vertex on the
    optimal face, which is a piecewise-constant (staircase) flux with at
    most m_active distinct levels.
    """

    def __init__(self, vertex_select=True, eps_abs=1e-10, eps_rel=1e-10,
                 max_iter=200000, polish=True, verbose=False):
        self.vertex_select = vertex_select
        self.eps_abs = eps_abs
        self.eps_rel = eps_rel
        self.max_iter = max_iter
        self.polish = polish
        self.verbose = verbose

    def solve(self, qp):
        n = qp['n']

        # Underdetermined: skip QP, go straight to LP vertex (χ²_min = 0 anyway).
        if self.vertex_select and n > qp['m_active']:
            x_v = _staircase_vertex(qp)
            if x_v is not None:
                out = _Result()
                out.x = x_v
                out.fun = float(0.5 * x_v @ qp['P'] @ x_v + qp['q'] @ x_v + qp['const'])
                out.success = True
                out.message = 'staircase vertex (HiGHS dual simplex)'
                out.nit = 0
                out.solve_time = 0.0
                out.backend = 'highs-vertex'
                out.staircase = True
                return out

        import osqp
        P, q = qp['P'], qp['q']

        # OSQP needs absolute-scale conditioning — rescale (NOT Tikhonov).
        p_scale = max(np.abs(P).max(), 1e-300)
        P_s = csc_matrix(P / p_scale)
        q_s = q / p_scale

        A = sp_vstack([qp['A_ord'], qp['I']], format='csc')
        l = np.concatenate([np.zeros(n - 1), np.full(n, qp['eps'])])
        u = np.full(2 * n - 1, np.inf)

        solver = osqp.OSQP()
        solver.setup(P_s, q_s, A, l, u,
                     eps_abs=self.eps_abs, eps_rel=self.eps_rel,
                     max_iter=self.max_iter, polish=self.polish,
                     verbose=self.verbose)
        res = solver.solve()

        if res.x is None or not np.all(np.isfinite(res.x)):
            raise RuntimeError(f"OSQP did not return usable x: {res.info.status}")

        x = np.asarray(res.x)
        out = _Result()
        out.backend = 'osqp'
        out.staircase = False
        out.nit = res.info.iter
        out.solve_time = res.info.run_time
        out.message = res.info.status
        out.success = True

        out.x = x
        chi2 = float(0.5 * x @ qp['P'] @ x + qp['q'] @ x + qp['const'])
        out.fun = chi2
        return out


class _CLARABELBackend:
    """CLARABEL-based QP solver — used for fixed-parameter / well-conditioned QPs."""

    def __init__(self, eps_abs=1e-10, eps_rel=1e-10, max_iter=200, verbose=False):
        self.eps_abs = eps_abs
        self.eps_rel = eps_rel
        self.max_iter = max_iter
        self.verbose = verbose

    def solve(self, qp, fix=None):
        import clarabel
        from scipy.sparse import csc_matrix as csc

        n = qp['n']
        p_scale = max(np.abs(qp['P']).max(), 1e-300)
        P = csc(qp['P'] / p_scale)
        q = qp['q'].copy() / p_scale

        # Constraints: ordering (NonnegativeCone via -A_ord x in NonnegativeCone),
        # bounds x >= eps, optional equality for fixed indices.
        rows = []
        cones = []

        # x[i] - x[i+1] >= 0  →  -(x[i] - x[i+1]) <= 0  →  Ax in Nonneg means Ax >= 0
        # clarabel uses Ax + s = b with s in cone; for NonnegativeCone: Ax <= b, residual b-Ax >= 0
        # Equivalent: -A_ord x <= 0
        rows.append((-qp['A_ord'], np.zeros(n - 1), clarabel.NonnegativeConeT(n - 1)))
        # x >= eps  →  -x <= -eps
        rows.append((-qp['I'], np.full(n, -qp['eps']), clarabel.NonnegativeConeT(n)))

        if fix is not None:
            idx = np.asarray(list(fix.keys()), dtype=int)
            vals = np.asarray([fix[i] for i in idx], dtype=float)
            E = np.zeros((len(idx), n))
            for r, i in enumerate(idx):
                E[r, i] = 1.0
            rows.append((csc(E), vals, clarabel.ZeroConeT(len(idx))))

        A = sp_vstack([r[0] for r in rows], format='csc')
        b = np.concatenate([r[1] for r in rows])
        cones = [r[2] for r in rows]

        settings = clarabel.DefaultSettings()
        settings.tol_gap_abs = self.eps_abs
        settings.tol_gap_rel = self.eps_rel
        settings.tol_feas = self.eps_abs
        settings.max_iter = self.max_iter
        settings.verbose = self.verbose

        solver = clarabel.DefaultSolver(P, q, A, b, cones, settings)
        sol = solver.solve()

        x = np.asarray(sol.x)
        out = _Result()
        out.backend = 'clarabel'
        out.staircase = False
        out.x = x
        # P, q passed to clarabel were rescaled by 1/p_scale; recover physical chi^2 from original qp.
        out.fun = float(0.5 * x @ qp['P'] @ x + qp['q'] @ x + qp['const'])
        out.success = str(sol.status) == 'Solved'
        out.message = str(sol.status)
        out.nit = sol.iterations
        out.solve_time = sol.solve_time
        return out


def run_optimize_qp(M_matrix, data, Bkg_vector, n, eps=0.0,
                    fix=None, vertex_select=True, verbose=False):
    """Solve χ²-minimization as a QP.

    Routing:
      - fix is None (free solve)        → OSQP. Underdetermined ⇒ HiGHS staircase vertex.
      - fix is dict {idx: value}        → CLARABEL (auto-route for fixed-parameter solves).
    """
    qp = _build_qp(M_matrix, data, Bkg_vector, n, eps)
    if fix:
        backend = _CLARABELBackend(verbose=verbose)
        return backend.solve(qp, fix=fix)
    backend = _OSQPBackend(vertex_select=vertex_select, verbose=verbose)
    return backend.solve(qp)


# Back-compat alias
def run_optimize_osqp(M_matrix, data, Bkg_vector, n, eps=0.0, **kw):
    return run_optimize_qp(M_matrix, data, Bkg_vector, n, eps=eps,
                            vertex_select=kw.pop('vertex_select', True),
                            verbose=kw.pop('verbose', False))
