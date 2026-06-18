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


def run_optimize(M_matrix, data, Bkg_vector, n, eps=0.0, x0=None, display=True,
                 vertex_select=True):
    """SciPy (trust-constr) Neyman-χ² solve, with optional staircase vertex.

    ``trust-constr`` is an interior-point method: it returns a *smooth ramp* on
    the interior of the optimal face, not a vertex, so the recovered flux is not
    piecewise-constant. When ``vertex_select`` (default) and the problem is
    under-determined (``n > #active bins``), we keep the converged fitted values
    ``mu = M_active @ x`` and pick a piecewise-constant *staircase* reproducing
    the same ``mu`` -- the same method-C step the OSQP backend uses -- so both
    solvers yield the same staircase estimate (χ² unchanged).
    """
    if x0 is None:
        x0 = np.ones(n)

    bounds = [(eps, None)] * n
    constraints = [ordering_constraint(n)]
    options = {'maxiter': 100000, 'xtol': 1e-9, 'gtol': 1e-9}
    if display:
        options['verbose'] = 3

    res = minimize(
        lambda x: chi_neyman(x, M_matrix, data, Bkg_vector, eps),
        x0,
        method='trust-constr',
        bounds=bounds,
        constraints=constraints,
        options=options,
        jac='3-point',
    )

    out = _Result()
    out.backend = 'scipy'
    out.staircase = False
    out.x = np.asarray(res.x)
    out.fun = float(res.fun)
    out.success = bool(res.success)
    out.message = str(res.message)
    out.nit = int(getattr(res, 'nit', 0))
    out.solve_time = 0.0

    # Vertex selection: collapse the smooth ramp to a piecewise-constant vertex
    # reproducing the same fitted values mu (χ² unchanged), matching OSQP.
    qp = _build_qp(M_matrix, data, Bkg_vector, n, eps)
    if vertex_select and qp['n'] > qp['m_active']:
        mu = qp['M_active'] @ out.x
        x_v = _vertex_select(qp, mu)
        if x_v is not None:
            out.x = x_v
            out.staircase = True
            out.backend = 'scipy+method-c-vertex'
            out.fun = float(0.5 * x_v @ qp['P'] @ x_v + qp['q'] @ x_v + qp['const'])
    return out


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


def _second_diff_operator(n: int) -> np.ndarray:
    """(n-2) x n second-difference matrix D2; ``||D2 x||^2`` is the discrete
    curvature (sum of squared second differences) of the staircase x."""
    D2 = np.zeros((max(n - 2, 0), n))
    for i in range(n - 2):
        D2[i, i], D2[i, i + 1], D2[i, i + 2] = 1.0, -2.0, 1.0
    return D2


def _clarabel_qp(P, q, blocks):
    """Solve  min 1/2 xᵀP x + qᵀx  s.t. the given (A, b, cone) blocks.

    Each block is ``(A, b, cone)`` in CLARABEL's ``A x + s = b, s in cone``
    convention (ZeroConeT -> equality, NonnegativeConeT -> A x <= b). Returns
    the solution x, or None if the problem is not solved (e.g. infeasible)."""
    import clarabel
    from scipy.sparse import csc_matrix as csc

    A = sp_vstack([blk[0] for blk in blocks], format='csc')
    b = np.concatenate([blk[1] for blk in blocks])
    cones = [blk[2] for blk in blocks]
    s = clarabel.DefaultSettings()
    s.verbose = False
    sol = clarabel.DefaultSolver(csc(P), np.asarray(q, float), A, b, cones, s).solve()
    return np.asarray(sol.x) if str(sol.status) == 'Solved' else None


def _curvature_min(M, mu, A_ord, eps):
    """Smoothest monotone non-negative flux reproducing the fitted counts mu:

        argmin ||D2 x||^2   s.t.  M x = mu,  x_i >= x_{i+1},  x >= eps

    The minimum-curvature element of the chi^2-optimal face is unbiased at the
    window edge (it does not drive the tail to zero the way a minimal-flux /
    tail-weighted objective does), so it tracks eta where eta != 0. Returns None
    if ``M x = mu`` is infeasible under monotone/non-negativity (e.g. noisy data
    with chi^2 > 0) -- the caller then supplies a reproducible ``mu`` instead."""
    import clarabel
    from scipy.sparse import csc_matrix as csc

    n = M.shape[1]
    D2 = _second_diff_operator(n)
    P = 2.0 * (D2.T @ D2) + 1e-10 * np.eye(n)        # tiny ridge -> strictly PD
    P = P / max(np.abs(P).max(), 1e-300)
    blocks = [
        (csc(np.asarray(M)), np.asarray(mu, float),
         clarabel.ZeroConeT(M.shape[0])),                       # M x = mu
        (-A_ord, np.zeros(n - 1), clarabel.NonnegativeConeT(n - 1)),  # x_i>=x_{i+1}
        (-sp_eye(n, format='csc'), np.full(n, -eps),
         clarabel.NonnegativeConeT(n)),                         # x >= eps
    ]
    return _clarabel_qp(P, np.zeros(n), blocks)


def _colnorm_vertex(M, mu, A_ord, eps):
    """Column-norm-weighted minimal-flux staircase: argmin sum_j ||M[:,j]|| x_j
    over {M x = mu, x_i >= x_{i+1}, x >= eps}, via a HiGHS simplex LP.

    Used as the fallback when the curvature QP (:func:`_curvature_min`) is too
    ill-conditioned to solve -- i.e. a steeply-falling (Bound) eta with a ~1e14
    dynamic-range cliff. The simplex is scale-invariant, so it solves the cliff
    trivially and reproduces the counts to machine precision. For such a cliff
    the counts are dominated by the low-v_min spike, so the weighted minimal-flux
    is forced to put the large flux there (it captures the cliff) and only
    collapses the count-irrelevant high-v_min halo floor -- a clean monotone
    staircase that tracks the cliff (validated in sandbox/tiebreak_test; for the
    cliff this is identical to the plain uniform Sum x_j vertex). It is NOT used
    for a gentle eta, where minimal-flux would collapse the informative window
    edge; there the curvature QP solves and method C is used instead. Returns
    None if the LP fails (e.g. no exact monotone fit for noisy data)."""
    n = M.shape[1]
    weight = np.sqrt((M * M).sum(axis=0))               # ||M[:,j]||
    res = linprog(weight, A_ub=-A_ord.toarray(), b_ub=np.zeros(n - 1),
                  A_eq=M, b_eq=np.asarray(mu, float),
                  bounds=[(eps, None)] * n, method='highs-ds')
    return np.asarray(res.x) if res.success else None


def _dp_segments(target, K):
    """Exact K-segment piecewise-constant L2 partition of ``target`` (dynamic
    programming). Returns the K segment bounds [(a0,b0), ...]. The (near-)
    monotone smooth ``target`` yields non-increasing segment means, so the
    staircase built on these bounds stays monotone. Vectorised over the inner
    split index so it is cheap even at n ~ 10^3 (used per Poisson toy)."""
    target = np.asarray(target, float)
    n = len(target)
    K = min(K, n)
    pre = np.concatenate([[0.0], np.cumsum(target)])
    pre2 = np.concatenate([[0.0], np.cumsum(target * target)])
    dp = np.full((K + 1, n + 1), np.inf)
    bk = np.zeros((K + 1, n + 1), dtype=int)
    dp[0, 0] = 0.0
    for k in range(1, K + 1):
        for i in range(k, n + 1):
            js = np.arange(k - 1, i)
            seg = (pre2[i] - pre2[js]) - (pre[i] - pre[js]) ** 2 / (i - js)
            vals = dp[k - 1, js] + seg
            t = int(np.argmin(vals))
            dp[k, i], bk[k, i] = vals[t], js[t]
    bnds, i = [], n
    for k in range(K, 0, -1):
        j = int(bk[k, i])
        bnds.append((j, i))
        i = j
    return bnds[::-1]


def _refit_levels(M, mu, bounds, n, eps):
    """Refit the K plateau VALUES to the fitted counts mu on fixed segment
    ``bounds``:  argmin ||G l - mu||^2  s.t.  l_k >= l_{k+1}, l >= eps, where
    column k of G sums the M columns inside segment k. Returns the assembled
    n-vector staircase (a genuine vertex with the smooth curve's plateaus)."""
    import clarabel
    from scipy.sparse import csc_matrix as csc

    K = len(bounds)
    G = np.zeros((M.shape[0], K))
    for k, (a, b) in enumerate(bounds):
        G[:, k] = np.asarray(M[:, a:b]).sum(axis=1)
    P = 2.0 * (G.T @ G) + 1e-10 * np.eye(K)
    sc = max(np.abs(P).max(), 1e-300)
    P, q = P / sc, (-2.0 * G.T @ np.asarray(mu, float)) / sc
    blocks = [(-sp_eye(K, format='csc'), np.full(K, -eps),
               clarabel.NonnegativeConeT(K))]                   # l >= eps
    if K >= 2:
        A_ord = np.zeros((K - 1, K))
        for i in range(K - 1):
            A_ord[i, i], A_ord[i, i + 1] = 1.0, -1.0
        blocks.insert(0, (csc(-A_ord), np.zeros(K - 1),
                          clarabel.NonnegativeConeT(K - 1)))      # l_k >= l_{k+1}
    lv = _clarabel_qp(P, q, blocks)
    if lv is None:
        return None
    x = np.empty(n)
    for k, (a, b) in enumerate(bounds):
        x[a:b] = lv[k]
    return x


def _vertex_select(qp, mu: np.ndarray):
    """Method C staircase tie-break on the QP's optimal face.

    The chi^2 minimiser is a whole face of the monotone polytope whenever
    ``n > rank(M)`` (here ~10^2-10^3 params vs ~5 bins). Among the many staircases
    that reproduce the fitted counts ``mu = M_active @ x`` we want the one that
    tracks eta WITHOUT collapsing/spiking at the window edge -- but eta does not
    vanish there, so the minimal-flux / tail-weighted vertex (neutrinoAnalysis,
    appropriate when the true tail vanishes) drives the edge to zero. Instead:

      1. ``_curvature_min`` -- the smoothest (min ||D2 x||^2) flux reproducing
         mu: the unbiased, weight-free representative of the optimal face;
      2. ``_dp_segments``   -- segment it into K = #active-bins monotone plateaus;
      3. ``_refit_levels``  -- refit those K plateau values to mu.

    The result is a genuine K-step staircase reproducing mu (chi^2 unchanged to
    a negligible DP-segmentation residual) that follows eta to the edge. No
    contrived weight is added -- the weighted-LP vertex is replaced by snapping
    the natural smooth solution to its nearest K-step staircase.

    Two tiers. For a gentle eta the plain curvature step solves and this is the
    whole method. For a steeply-falling eta (the Bound population, ~1e14 dynamic
    range) the curvature QP is too ill-conditioned and returns None; we then fall
    back to the column-norm-weighted minimal-flux simplex vertex
    (:func:`_colnorm_vertex`), which is scale-invariant so it solves the cliff
    trivially. There minimal-flux is the right choice -- the counts are dominated
    by the low-v_min cliff, so it captures the cliff and only collapses the
    count-irrelevant halo floor (it would, however, collapse a gentle window edge,
    which is why method C is used whenever the curvature QP solves). Returns None
    only if ``M x = mu`` is genuinely infeasible (noisy data) -- the OSQP path
    then supplies a reproducible mu.
    """
    M = qp['M_active']
    mu = np.asarray(mu)
    # Tier 1: method C (gentle eta) -- smooth curvature, DP segments, level refit.
    x_smooth = _curvature_min(M, mu, qp['A_ord'], qp['eps'])
    if x_smooth is not None:
        bounds = _dp_segments(x_smooth, M.shape[0])
        return _refit_levels(M, mu, bounds, qp['n'], qp['eps'])
    # Tier 2: column-norm-weighted minimal-flux vertex (steeply-falling Bound cliff).
    return _colnorm_vertex(M, mu, qp['A_ord'], qp['eps'])


class _OSQPBackend:
    """OSQP QP + method-C staircase selection.

    Two stages:

      1. solve the Neyman-χ² QP with OSQP directly in physical ``x`` -> a smooth
         interior solution and the fitted values ``mu = M_active @ x``;
      2. if ``vertex_select`` (default), replace that ramp with a
         piecewise-constant *staircase* reproducing the same ``mu`` via method C
         (:func:`_vertex_select`: min-curvature solution -> DP segments -> level
         refit). χ² is unchanged because ``mu`` (hence the residual) is identical.

    No column scaling: with the raised-GeV unit base the design-matrix columns
    are O(1), so the QP is well-conditioned in physical units directly.
    """

    def __init__(self, vertex_select=True, eps_abs=1e-10, eps_rel=1e-10,
                 max_iter=200000, polish=True, verbose=False):
        self.vertex_select = vertex_select
        self.eps_abs = eps_abs
        self.eps_rel = eps_rel
        self.max_iter = max_iter
        self.polish = polish
        self.verbose = verbose

    def _osqp_interior(self, qp):
        """Solve the χ² QP directly in physical ``x`` (no column scaling).

        With the raised-GeV unit base the design-matrix columns are O(1), so
        OSQP is well-conditioned without the unit-norm reparametrisation.
        """
        import osqp

        n = qp['n']
        M = qp['M_active']
        inv_d = qp['inv_d']
        yv = qp['y_active']

        P_x = 2.0 * (M.T * inv_d) @ M
        q_x = -2.0 * M.T @ (yv * inv_d)
        p_scale = max(np.abs(P_x).max(), 1e-300)

        # Constraints in x: ordering A_ord·x ≥ 0 and x ≥ eps.
        A = sp_vstack([qp['A_ord'], qp['I']], format='csc')
        l = np.concatenate([np.zeros(n - 1), np.full(n, qp['eps'])])
        u = np.full(2 * n - 1, np.inf)

        solver = osqp.OSQP()
        solver.setup(csc_matrix(P_x / p_scale), q_x / p_scale, A, l, u,
                     eps_abs=self.eps_abs, eps_rel=self.eps_rel,
                     max_iter=self.max_iter, polish=self.polish,
                     verbose=self.verbose)
        res = solver.solve()
        if res.x is None or not np.all(np.isfinite(res.x)):
            raise RuntimeError(f"OSQP did not return usable x: {res.info.status}")
        return np.asarray(res.x), res

    def solve(self, qp):
        out = _Result()

        # Stage 0: exact-fit staircase. The chi^2 = 0 minimiser exists whenever
        # the signal y = data - background is reproducible by a monotone
        # non-negative flux -- always the case for the self-consistent forward
        # model, where y = M @ eta and eta is itself monotone and >= 0. Building
        # the method-C staircase directly at mu = y bypasses the QP interior,
        # which is badly conditioned when the counts span many orders of
        # magnitude (the dense, steeply falling Bound eta makes OSQP report
        # "dual infeasible"). The K-step staircase reproduces y up to a small
        # DP-segmentation residual (the staircase's own fit), so the acceptance
        # tolerance is relative, not machine-precision. For genuinely noisy data
        # M x = y is infeasible, _vertex_select returns None, and we fall through
        # to the OSQP interior, which supplies a reproducible mu for Stage 2.
        if self.vertex_select and qp['n'] > qp['m_active']:
            y = qp['y_active']
            x_exact = _vertex_select(qp, y)
            if x_exact is not None and (
                    np.linalg.norm(qp['M_active'] @ x_exact - y)
                    <= 1e-2 * max(np.linalg.norm(y), 1e-300)):
                out.backend = 'method-c-exact'
                out.staircase = True
                out.success = True
                out.message = 'exact-fit staircase (method C)'
                out.nit = 0
                out.solve_time = 0.0
                out.x = x_exact
                out.fun = float(0.5 * x_exact @ qp['P'] @ x_exact
                                + qp['q'] @ x_exact + qp['const'])
                return out

        # Stage 1 (fallback): QP interior solution (smooth ramp) + fitted mu.
        # Needed when no exact monotone fit exists (e.g. Poisson pseudo-data),
        # so the chi^2 minimum is > 0 and mu must come from the QP.
        x_int, res = self._osqp_interior(qp)

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
            x_v = _vertex_select(qp, mu)
            if x_v is not None:
                out.x = x_v
                out.staircase = True
                out.backend = 'osqp+method-c-vertex'

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
      - fix is None (free solve)        → OSQP. Underdetermined ⇒ method-C staircase.
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
