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


def _column_scale(M: np.ndarray) -> np.ndarray:
    """Per-column reciprocal norm ``D_j = 1/‖M[:,j]‖`` (unit-norm columns).

    Deterministic conditioning (the key trick in kyphys-creator/neutrinoAnalysis):
    the optimal flux and the design-matrix columns otherwise span many orders of
    magnitude, and OSQP/CLARABEL terminate early at a wrong point while reporting
    success. Reparametrising ``x = D ⊙ y`` makes every column O(1). Because ``D``
    divides out the overall conditioning constant ``c`` (``config.CONDITION_C``)
    too, the recovered flux is independent of ``c`` -- no per-mass tuning. ``D``
    is data-independent, so it can be reused across pseudo-data (Monte Carlo).
    """
    cn = np.sqrt((M * M).sum(axis=0))
    cn[cn == 0] = 1.0
    return 1.0 / cn


def _vertex_select(qp, D: np.ndarray, mu: np.ndarray):
    """HiGHS-simplex vertex on the QP's optimal face (neutrinoAnalysis method).

    The χ² minimiser is a whole face of the monotone polytope whenever
    ``n > rank(M)`` (here ~948 params vs 5 bins): OSQP returns a smooth interior
    *ramp*, but the physically meaningful estimate is a *vertex*
    (piecewise-constant flux). Once the QP fixes the fitted values
    ``mu = M_active @ x``, we pick a vertex of
    ``{x : M_active x = mu, x_i ≥ x_{i+1}, x ≥ eps}`` with a simplex method.

    Objective: minimise the total flux ``Σ x_j`` -- the minimal-norm monotone
    flux that still reproduces ``mu``. Minimising the total pushes every step
    down to the smallest value the monotone + count constraints allow, so the
    under-constrained ends both collapse: the high-v_min *tail* (where the small
    response barely affects the counts) goes to zero, and the low-v_min *head*
    is not inflated above what the counts require -- it just meets eta instead of
    overshooting it.

    (Earlier this minimised a geometric tail weight ``Σ R^(i/(n-1)) x_i``. The
    ``R`` factor under-weighted the head (``R^0 = 1``), letting the first step
    overshoot eta by a few percent, while adding nothing the plain total ``Σ x``
    doesn't already give for the tail.)

    Solved in the column-scaled variable ``x = D ⊙ z`` (``M_s = M_active·D`` has
    unit-norm, ``c``-independent columns), so HiGHS sees a well-conditioned
    system and the picked vertex does not depend on the conditioning constant.
    """
    n = qp['n']
    M_s = qp['M_active'] * D[None, :]          # unit-norm columns, c-free
    mu = np.asarray(mu)

    # Minimise total physical flux Σ x_j. In the scaled variable z (x = D⊙z)
    # that is Σ D_j z_j, i.e. weight D on z. Normalise by the max to strip the
    # overall 1/c factor in D -> a c-independent, well-scaled LP.
    weight = D / D.max()
    # ordering on physical x = D⊙z:  D_i z_i - D_{i+1} z_{i+1} >= 0
    A_ord_z = -(qp['A_ord'].toarray() * D[None, :])
    bounds = [(qp['eps'] / d if d else 0.0, None) for d in D]

    res = linprog(
        weight, A_ub=A_ord_z, b_ub=np.zeros(n - 1),
        A_eq=M_s, b_eq=mu,
        bounds=bounds, method='highs-ds',
        options={'presolve': True},
    )
    if not res.success:
        return None
    return D * np.asarray(res.x)


class _OSQPBackend:
    """Column-scaled OSQP QP + simplex vertex selection (neutrinoAnalysis method).

    Two stages, mirroring kyphys-creator/neutrinoAnalysis:

      1. solve the Neyman-χ² QP with OSQP in the *column-scaled* variable
         ``x = D ⊙ y`` (``D`` from :func:`_column_scale`) -> a smooth interior
         solution and the fitted values ``mu = M_active @ x``;
      2. if ``vertex_select`` (default), replace that ramp with a
         piecewise-constant *vertex* reproducing the same ``mu`` via a
         tail-weighted HiGHS simplex LP (:func:`_vertex_select`). χ² is
         unchanged because ``mu`` (hence the residual) is identical.

    Column scaling makes the result independent of the conditioning constant
    ``c`` and removes the need for any per-mass tuning.
    """

    def __init__(self, vertex_select=True, eps_abs=1e-10, eps_rel=1e-10,
                 max_iter=200000, polish=True, verbose=False):
        self.vertex_select = vertex_select
        self.eps_abs = eps_abs
        self.eps_rel = eps_rel
        self.max_iter = max_iter
        self.polish = polish
        self.verbose = verbose

    def _osqp_interior(self, qp, D):
        """Solve the χ² QP in the column-scaled variable z (x = D ⊙ z).

        Build the QP from the unit-norm matrix ``M_s = M_active·D`` directly so
        the conditioning constant ``c`` cancels exactly (no ``c²``-scaled
        intermediate), giving a well-scaled, ``c``-independent problem.
        """
        import osqp
        from scipy.sparse import diags as sp_diags

        n = qp['n']
        M_s = qp['M_active'] * D[None, :]         # unit-norm columns, c-free
        inv_d = qp['inv_d']
        yv = qp['y_active']

        P_z = 2.0 * (M_s.T * inv_d) @ M_s
        q_z = -2.0 * M_s.T @ (yv * inv_d)
        p_scale = max(np.abs(P_z).max(), 1e-300)

        # Constraints in z: ordering A_ord·(D⊙z) ≥ 0 and x = D⊙z ≥ eps.
        Dm = sp_diags(D)
        A = sp_vstack([qp['A_ord'] @ Dm, Dm], format='csc')
        l = np.concatenate([np.zeros(n - 1), np.full(n, qp['eps'])])
        u = np.full(2 * n - 1, np.inf)

        solver = osqp.OSQP()
        solver.setup(csc_matrix(P_z / p_scale), q_z / p_scale, A, l, u,
                     eps_abs=self.eps_abs, eps_rel=self.eps_rel,
                     max_iter=self.max_iter, polish=self.polish,
                     verbose=self.verbose)
        res = solver.solve()
        if res.x is None or not np.all(np.isfinite(res.x)):
            raise RuntimeError(f"OSQP did not return usable x: {res.info.status}")
        return D * np.asarray(res.x), res        # x in m_cond space, info

    def solve(self, qp):
        D = _column_scale(qp['M_active'])
        out = _Result()

        # Stage 0: exact-fit vertex. The chi^2 = 0 minimiser exists whenever the
        # signal y = data - background is reproducible by a monotone non-negative
        # flux -- always the case for the self-consistent forward model, where
        # y = M @ eta and eta is itself monotone and >= 0. Solving the vertex LP
        # directly at mu = y bypasses the QP interior, which is badly conditioned
        # when the counts span many orders of magnitude (the dense, steeply
        # falling Bound eta makes OSQP report "dual infeasible"). It reduces to
        # the same vertex the QP path would pick when that path converges, so it
        # is the default fast/robust route, not a special case.
        if self.vertex_select and qp['n'] > qp['m_active']:
            y = qp['y_active']
            x_exact = _vertex_select(qp, D, y)
            if x_exact is not None and np.allclose(qp['M_active'] @ x_exact, y,
                                                   rtol=1e-6, atol=0):
                out.backend = 'highs-vertex-exact'
                out.staircase = True
                out.success = True
                out.message = 'exact-fit vertex'
                out.nit = 0
                out.solve_time = 0.0
                out.x = x_exact
                out.fun = float(0.5 * x_exact @ qp['P'] @ x_exact
                                + qp['q'] @ x_exact + qp['const'])
                return out

        # Stage 1 (fallback): QP interior solution (smooth ramp) + fitted mu.
        # Needed when no exact monotone fit exists (e.g. Poisson pseudo-data),
        # so the chi^2 minimum is > 0 and mu must come from the QP.
        x_int, res = self._osqp_interior(qp, D)

        out.backend = 'osqp'
        out.staircase = False
        out.nit = res.info.iter
        out.solve_time = res.info.run_time
        out.message = res.info.status
        out.success = True
        out.x = x_int

        # Stage 2: vertex selection on mu = M_active @ x_interior.
        if self.vertex_select and qp['n'] > qp['m_active']:
            mu = qp['M_active'] @ x_int
            x_v = _vertex_select(qp, D, mu)
            if x_v is not None:
                out.x = x_v
                out.staircase = True
                out.backend = 'osqp+highs-vertex'

        x = out.x
        out.fun = float(0.5 * x @ qp['P'] @ x + qp['q'] @ x + qp['const'])
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
