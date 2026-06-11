"""Full statistical layer: Poisson pseudo-experiments and model discrimination.

Two independent questions, two tools:

1. **How well is the recovered flux determined?** -- ``fit_toys`` draws Poisson
   pseudo-data from the forward model's expected counts and re-runs the full
   monotone non-negative inversion on every toy. The ensemble of recovered
   staircases gives per-``v_min`` percentile *bands* (``flux_band``) around the
   self-consistent (Asimov) recovery.

2. **Can two velocity models be told apart?** -- the inversion itself cannot
   discriminate (any monotone eta fits some monotone flux), so discrimination
   is a *hypothesis test between the two fully-specified count vectors*.
   ``asimov_significance`` gives the median expected significance from the
   Poisson log-likelihood ratio evaluated on the Asimov dataset (Cowan et al.,
   arXiv:1007.1727); ``mc_significance`` validates it by Monte Carlo where the
   significance is small enough to resolve.

Everything reuses the existing forward model and solver: toys go through the
same ``condition`` + ``run_optimize_qp`` path as the real fit, so the bands
inherit the staircase/vertex estimator exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm

from .config import RunConfig, CONDITION_C
from .model import condition
from .optimizer import run_optimize_qp


# ---------------------------------------------------------------------------
# 1. Pseudo-experiment ensemble for the recovered flux
# ---------------------------------------------------------------------------

@dataclass
class ToyEnsemble:
    """Result of ``fit_toys``: one inversion per Poisson pseudo-experiment."""

    fluxes: np.ndarray          # (n_toys, n_vmin) recovered staircases
    chi2: np.ndarray            # (n_toys,) true Neyman chi^2 at the optimum
    observed: np.ndarray        # (n_toys, n_ebins) the Poisson draws
    expected: np.ndarray        # (n_ebins,) the Asimov expectation used
    n_failed: int = 0           # toys whose solve did not return a usable flux
    seed: int | None = None
    config: RunConfig | None = None

    @property
    def n_toys(self) -> int:
        return self.fluxes.shape[0]


def fit_toys(analysis, n_toys: int = 300, seed: int = 0,
             progress: bool = False) -> ToyEnsemble:
    """Poisson-fluctuate the expected counts and re-run the inversion per toy.

    Each toy draws ``observed ~ Poisson(signal + background)`` and solves the
    same conditioned QP as :meth:`DarkMatterQuantumAnalysis.optimize` (OSQP
    interior + HiGHS vertex). The exact-fit fast path does not apply to noisy
    data, so every toy exercises the genuine chi^2 > 0 route.

    Bins drawn at 0 are dropped by the Neyman mask inside the solver (they
    carry no Neyman weight); with the count levels of this analysis that is
    rare outside the heaviest-mass configurations.
    """
    rng = np.random.default_rng(seed)
    expected = analysis.observed                    # Asimov: signal + background
    n_vmin = analysis.n_vmin

    draws = rng.poisson(expected, size=(n_toys, expected.size)).astype(float)

    fluxes = np.full((n_toys, n_vmin), np.nan)
    chi2 = np.full(n_toys, np.nan)
    n_failed = 0
    for i, data in enumerate(draws):
        m_cond, data_cond, bkg_cond, unscale = condition(
            analysis.m_phys, data, analysis.background, c=CONDITION_C)
        try:
            res = run_optimize_qp(m_cond, data_cond, bkg_cond, n_vmin)
            fluxes[i] = unscale(res.x)
            chi2[i] = res.fun
        except Exception:
            n_failed += 1
        if progress and (i + 1) % 50 == 0:
            print(f"  toys {i + 1}/{n_toys}")

    ok = ~np.isnan(chi2)
    return ToyEnsemble(fluxes=fluxes[ok], chi2=chi2[ok], observed=draws[ok],
                       expected=expected, n_failed=n_failed, seed=seed,
                       config=analysis.config)


def flux_band(ensemble: ToyEnsemble,
              levels: tuple[float, ...] = (68.0, 95.0)) -> dict:
    """Per-``v_min`` percentile bands of the toy fluxes.

    Returns ``{"median": (n_vmin,), <level>: (lo, hi)}`` where ``lo``/``hi``
    are the symmetric percentile envelopes (e.g. 16/84 for the 68 band).
    """
    out = {"median": np.nanmedian(ensemble.fluxes, axis=0)}
    for lv in levels:
        half = (100.0 - lv) / 2.0
        lo = np.nanpercentile(ensemble.fluxes, half, axis=0)
        hi = np.nanpercentile(ensemble.fluxes, 100.0 - half, axis=0)
        out[lv] = (lo, hi)
    return out


# ---------------------------------------------------------------------------
# 2. Model discrimination (simple vs simple Poisson hypothesis test)
# ---------------------------------------------------------------------------

def _llr(n: np.ndarray, mu_null: np.ndarray, mu_alt: np.ndarray) -> np.ndarray:
    """q = -2 ln(L_null / L_alt) for Poisson counts ``n`` (per toy row).

    Only bins where both expectations are positive contribute; both hypotheses
    here always include the same background, so a zero-expectation bin in one
    is zero in the other too.
    """
    mask = (mu_null > 0) & (mu_alt > 0)
    a, b = mu_null[mask], mu_alt[mask]
    n = np.atleast_2d(n)[:, mask]
    return 2.0 * ((a - b).sum() + n @ np.log(b / a))


def asimov_significance(mu_null: np.ndarray, mu_alt: np.ndarray) -> dict:
    """Median expected significance for rejecting ``mu_null`` when ``mu_alt``
    is true, from the Asimov dataset ``n = mu_alt`` (arXiv:1007.1727).

    Returns ``q_asimov`` (the LLR on the Asimov data), ``sigma`` (its square
    root, the median significance), and ``chi2_pearson``
    (``sum (mu_alt - mu_null)^2 / mu_null``, the Gaussian-limit check).
    """
    mu_null = np.asarray(mu_null, float)
    mu_alt = np.asarray(mu_alt, float)
    q_asimov = float(_llr(mu_alt, mu_null, mu_alt)[0])
    mask = mu_null > 0
    chi2 = float((((mu_alt - mu_null) ** 2)[mask] / mu_null[mask]).sum())
    return {"q_asimov": q_asimov,
            "sigma": float(np.sqrt(max(q_asimov, 0.0))),
            "chi2_pearson": chi2}


def mc_significance(mu_null: np.ndarray, mu_alt: np.ndarray,
                    n_toys: int = 100_000, seed: int = 0) -> dict:
    """Monte-Carlo check of the Asimov significance.

    Draws the LLR distribution under the null, takes the median LLR under the
    alternative, and converts the null-tail p-value to one-sided sigma. With
    ``n_toys`` toys the smallest resolvable p is ~1/n_toys (~4.3 sigma at 1e5);
    beyond that ``sigma`` is reported as ``inf`` and ``p_value`` as 0 -- quote
    the Asimov number there.
    """
    rng = np.random.default_rng(seed)
    mu_null = np.asarray(mu_null, float)
    mu_alt = np.asarray(mu_alt, float)

    q_null = _llr(rng.poisson(mu_null, size=(n_toys, mu_null.size)),
                  mu_null, mu_alt)
    q_alt = _llr(rng.poisson(mu_alt, size=(n_toys, mu_alt.size)),
                 mu_null, mu_alt)
    q_med = float(np.median(q_alt))
    p = float((q_null >= q_med).mean())
    return {"p_value": p,
            "sigma": float(norm.isf(p)) if p > 0 else float("inf"),
            "q_median_alt": q_med}


# ---------------------------------------------------------------------------
# 3. Discrimination grid over configurations
# ---------------------------------------------------------------------------

#: Alternatives tested against the pure-SHM (Halo) null, as RunConfig deltas.
ALTERNATIVES: dict[str, dict] = {
    "mix5":  {"disk_fraction": 0.05},
    "mix25": {"disk_fraction": 0.25},
    "Disk":  {"eta": "Disk"},
    "Bound": {"eta": "Bound"},       # Halo + Earth-bound population (M3 only)
}


@dataclass
class DiscriminationRow:
    material: str
    q: str
    mass: str
    nbins: int
    background: str
    alternative: str
    sigma_asimov: float
    q_asimov: float
    chi2_pearson: float
    sigma_mc: float | None = None
    p_value_mc: float | None = None


def discrimination_row(analysis_null, analysis_alt, alternative: str,
                       mc_toys: int | None = 100_000,
                       seed: int = 0) -> DiscriminationRow:
    """Significance for one (null, alternative) pair of built analyses."""
    res = asimov_significance(analysis_null.observed, analysis_alt.observed)
    row = DiscriminationRow(
        material=analysis_null.config.material, q=analysis_null.config.q,
        mass=analysis_null.config.mass, nbins=analysis_null.config.nbins,
        background=analysis_null.config.background, alternative=alternative,
        sigma_asimov=res["sigma"], q_asimov=res["q_asimov"],
        chi2_pearson=res["chi2_pearson"])
    if mc_toys:
        mc = mc_significance(analysis_null.observed, analysis_alt.observed,
                             n_toys=mc_toys, seed=seed)
        row.sigma_mc = mc["sigma"]
        row.p_value_mc = mc["p_value"]
    return row


# ---------------------------------------------------------------------------
# 4. Point-wise profile confidence band (kyphys-creator/neutrinoAnalysis)
# ---------------------------------------------------------------------------
#
# Feldman-Cousins-style construction, per flux parameter x_j:
#
#   1. observed Delta-chi^2: fit the real data free (chi2_free) and with x_j
#      fixed at a trial value v (chi2_fixed); dchi2_obs = chi2_fixed - chi2_free.
#   2. calibration: generate pseudo-data from the FIXED-fit model
#      (mu = M @ x_fixed + bkg) and refit each pseudo both free and fixed at
#      the same v; the CL-quantile of the |Delta-chi^2| distribution is the
#      cutoff at this v (no chi^2(1) assumption -- the monotone non-negative
#      constraints make the true distribution non-trivial).
#   3. v is inside the band iff dchi2_obs < cutoff.
#   4. the band edge is located by outward geometric bracketing from the best
#      fit, then geometric bisection, with more pseudo-data during refinement.
#
# Differences from the neutrino code: pseudo-data is drawn Poisson (their
# Gaussian sqrt(N) approximation is unsafe at the O(10) counts of the heaviest
# masses here), and the fixed-parameter solve routes to CLARABEL through
# ``run_optimize_qp(fix=...)``. As there, the SAME seed is reused for every
# trial value (common random numbers), so the inside/outside decision varies
# smoothly along the scan and the bisection does not jitter.


def _profile_solve(analysis, data, fix: dict | None = None):
    """One conditioned solve; returns (true chi^2, physical x).

    Both the free and the fixed fits go through CLARABEL: its interior-point
    method solves this QP in ~10 iterations where OSQP's ADMM needs ~10^4
    (measured 27 ms vs 1.1 s at n_vmin=286, chi^2 agreement ~3e-5). Only the
    chi^2 value matters here (Delta-chi^2), which is invariant under where on
    the optimal face the solver lands, so no vertex selection is run.
    """
    from .optimizer import _build_qp, _CLARABELBackend

    m_cond, data_cond, bkg_cond, unscale = condition(
        analysis.m_phys, data, analysis.background, c=CONDITION_C)
    fix_cond = ({i: v / CONDITION_C for i, v in fix.items()}
                if fix is not None else None)
    qp = _build_qp(m_cond, data_cond, bkg_cond, analysis.n_vmin, 0.0)
    res = _CLARABELBackend(eps_abs=1e-9, eps_rel=1e-9).solve(qp, fix=fix_cond)
    return float(res.fun), unscale(res.x)


def _band_eval(analysis, index: int, v: float, levels, n_pseudo: int,
               seed: int, cache: dict, chi2_free_obs: float,
               n_jobs: int = 8) -> dict:
    """Inside/outside decision at one trial value ``v`` (cached by value).

    The pseudo-experiment fits run on a thread pool: CLARABEL releases the
    GIL during its solve, giving ~3.5x on Apple Silicon (processes do no
    better -- the solves saturate memory bandwidth, not the GIL). All Poisson
    draws are made serially up front, so the result is independent of
    ``n_jobs``.
    """
    key = float(f"{v:.6e}")           # relative rounding (values are ~1e-32)
    if key in cache:
        return cache[key]

    chi2_fixed_obs, x_fixed = _profile_solve(analysis, analysis.observed,
                                             fix={index: v})
    dchi2_obs = abs(chi2_fixed_obs - chi2_free_obs)

    mu = analysis.m_phys @ x_fixed + analysis.background
    rng = np.random.default_rng(seed)       # common random numbers across v
    draws = [rng.poisson(mu).astype(float) for _ in range(n_pseudo)]

    def one(pseudo):
        try:
            c_free, _ = _profile_solve(analysis, pseudo)
            c_fix, _ = _profile_solve(analysis, pseudo, fix={index: v})
            return abs(c_fix - c_free)
        except Exception:
            return None

    if n_jobs > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(n_jobs) as ex:
            dchi2_mc = [d for d in ex.map(one, draws) if d is not None]
    else:
        dchi2_mc = [d for d in map(one, draws) if d is not None]

    ds = np.sort(dchi2_mc)
    included, cutoff = {}, {}
    for lv in levels:
        cut = ds[max(int(lv * len(ds)) - 1, 0)] if len(ds) else np.inf
        cutoff[lv] = float(cut)
        included[lv] = bool(dchi2_obs < cut)

    res = {"v": float(v), "dchi2_obs": float(dchi2_obs),
           "included": included, "cutoff": cutoff, "n_mc": len(ds)}
    cache[key] = res
    return res


def _bisect_edge(analysis, index, level, v_in, v_out, levels, n_pseudo_edge,
                 rel_tol, seed, cache, chi2_free_obs, n_jobs=8):
    """Geometric bisection between an included and an excluded value."""
    a, b = float(v_in), float(v_out)
    for _ in range(40):
        if abs(b - a) <= rel_tol * max(abs(b), abs(a)):
            break
        m = np.sqrt(a * b) if (a > 0 and b > 0) else 0.5 * (a + b)
        r = _band_eval(analysis, index, m, levels, n_pseudo_edge,
                       seed, cache, chi2_free_obs, n_jobs=n_jobs)
        if r["included"][level]:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def find_confidence_band(analysis, index: int,
                         levels: tuple[float, ...] = (0.68, 0.954),
                         num_pseudo: int = 30, n_pseudo_edge: int = 100,
                         step: float = 1.5, rel_tol: float = 0.05,
                         max_bracket: int = 40, seed: int = 42,
                         n_jobs: int = 8, verbose: bool = False) -> dict:
    """Point-wise profile confidence interval for the flux step at ``index``.

    Mirrors ``neutrinoAnalysis.find_confidence_band``: bracket outward from the
    best fit with factor ``step`` until the widest level excludes, then pin each
    level's edge by geometric bisection (``n_pseudo_edge`` pseudo-experiments
    during refinement, ``num_pseudo`` while bracketing; the per-value cache is
    shared). Requires ``analysis.optimize()`` to have been run.
    """
    if analysis.flux is None:
        raise RuntimeError("run optimize() first")
    levels = tuple(sorted(levels))
    widest = levels[-1]
    cache: dict = {}
    chi2_free_obs = float(analysis.result.fun)

    v0 = float(analysis.flux[index])
    # A zero best-fit step (the collapsed tail) still has a nonzero upper
    # limit; start the upward bracket from a small positive scale instead.
    v_start = v0 if v0 > 0 else float(analysis.flux.max()) * 1e-4

    r0 = _band_eval(analysis, index, max(v0, v_start * 1e-12), levels,
                    num_pseudo, seed, cache, chi2_free_obs, n_jobs=n_jobs)
    if verbose and not r0["included"][widest]:
        print(f"[warn] best fit at index {index} already outside the "
              f"{widest} band")

    def bracket(direction):
        v = v_start
        for _ in range(max_bracket):
            v = v * step if direction > 0 else v / step
            if v <= 0:
                return None
            r = _band_eval(analysis, index, v, levels, num_pseudo,
                           seed, cache, chi2_free_obs, n_jobs=n_jobs)
            if not r["included"][widest]:
                return v
        return None

    up_out = bracket(+1)
    lo_out = bracket(-1) if v0 > 0 else None    # x >= 0: zero step -> lower = 0

    band = {}
    for lv in levels:
        upper = (_bisect_edge(analysis, index, lv, max(v0, v_start), up_out,
                              levels, n_pseudo_edge, rel_tol, seed, cache,
                              chi2_free_obs, n_jobs=n_jobs)
                 if up_out is not None else np.inf)
        lower = (_bisect_edge(analysis, index, lv, v0, lo_out, levels,
                              n_pseudo_edge, rel_tol, seed, cache,
                              chi2_free_obs, n_jobs=n_jobs)
                 if lo_out is not None else 0.0)
        band[lv] = (float(lower), float(upper))
        if verbose:
            print(f"  idx {index} level {lv}: [{lower:.4g}, {upper:.4g}]")

    return {"index": int(index),
            "vmin_mid": float(analysis.vmin_mid[index]),
            "best_fit": v0,
            "levels": levels,
            "band": band,
            "n_evaluations": len(cache)}


def save_pointwise_band(analysis, bands: list[dict], out_dir=None) -> "Path":
    """Write one JSON per scanned point under ``<run_dir>/band/``.

    Mirrors neutrinoAnalysis's ``bands/band_*idx*.json`` layout: each scanned
    v_min index gets its own ``band_idx<index>.json``, so points can be added,
    recomputed or diffed individually. Returns the ``band/`` directory.
    """
    import json
    from pathlib import Path
    from .plotting import run_dir

    out = Path(out_dir) if out_dir is not None else run_dir(analysis)
    out = out / "band"
    out.mkdir(parents=True, exist_ok=True)
    for b in bands:
        rec = {**b, "levels": list(b["levels"]),
               "band": {str(k): list(v) for k, v in b["band"].items()}}
        with open(out / f"band_idx{b['index']:04d}.json", "w") as f:
            json.dump(rec, f, indent=1)
    return out


def load_pointwise_band(run_directory) -> list[dict]:
    """Read ``band/band_idx*.json`` back into a ``pointwise_band``-style list
    (sorted by index; levels/band keys restored to floats)."""
    import json
    from pathlib import Path

    band_dir = Path(run_directory) / "band"
    bands = []
    for path in sorted(band_dir.glob("band_idx*.json")):
        with open(path) as f:
            b = json.load(f)
        b["levels"] = tuple(float(l) for l in b["levels"])
        b["band"] = {float(k): tuple(v) for k, v in b["band"].items()}
        bands.append(b)
    if not bands:
        raise FileNotFoundError(f"no band_idx*.json under {band_dir}")
    return sorted(bands, key=lambda b: b["index"])


def band_table(bands: list[dict]):
    """Readable pandas table of ``pointwise_band`` results: physical units
    (cm^-1, via the same conversion the plots use) plus relative 68% errors."""
    import pandas as pd
    from .constants import CM

    levels = bands[0]["levels"]
    rows = []
    for b in bands:
        f = b["best_fit"]
        row = {"index": b["index"], "vmin [km/s]": round(b["vmin_mid"], 1),
               "best fit [cm^-1]": f * CM}
        for lv in levels:
            lo, hi = b["band"][lv]
            label = _SIGMA_LABEL.get(round(lv, 3), f"{lv:g}")
            row[f"lo {label} [cm^-1]"] = lo * CM
            row[f"hi {label} [cm^-1]"] = hi * CM
        lo, hi = b["band"][levels[0]]
        row["-68%"] = f"-{(1 - lo / f) * 100:.0f}%" if f > 0 else "-"
        row["+68%"] = (f"+{(hi / f - 1) * 100:.0f}%"
                       if f > 0 and np.isfinite(hi) else "-")
        row["evals"] = b["n_evaluations"]
        rows.append(row)
    return pd.DataFrame(rows)


def save_band_products(analysis, bands: list[dict]):
    """Save all band products into the run folder: the summary CSV
    (``flux_profile_band.csv``, natural units), the per-point JSONs under
    ``band/`` and the figure (``flux_profile_band.pdf``)."""
    from .plotting import plot_flux_with_pointwise_bands, run_dir

    out = run_dir(analysis)
    out.mkdir(parents=True, exist_ok=True)
    levels = bands[0]["levels"]
    table = np.array([[b["vmin_mid"], b["best_fit"]]
                      + [x for lv in levels for x in b["band"][lv]]
                      for b in bands])
    header = "vmin_mid,best_fit," + ",".join(
        f"lo{round(lv * 100)},hi{round(lv * 100)}" for lv in levels)
    np.savetxt(out / "flux_profile_band.csv", table, delimiter=",",
               comments="", header=header)
    band_dir = save_pointwise_band(analysis, bands)
    plot_flux_with_pointwise_bands(analysis, bands)
    print(f"saved: flux_profile_band.csv / {len(bands)} json under "
          f"{band_dir.name}/ / pdf  ({out})")
    return out


def pointwise_band(analysis, indices=None, n_indices: int = 12,
                   verbose: bool = True, **kwargs) -> list[dict]:
    """Profile band at a set of v_min indices (default: log-spaced over the
    populated window). Returns one ``find_confidence_band`` dict per index."""
    if analysis.flux is None:
        raise RuntimeError("run optimize() first")
    if indices is None:
        pos = np.flatnonzero(analysis.flux > 0)
        lo, hi = (0, analysis.n_vmin - 1) if pos.size == 0 else (pos[0], pos[-1])
        grid = np.unique(np.geomspace(lo + 1, hi + 1, n_indices).astype(int) - 1)
        indices = grid.tolist()
    out = []
    for k, idx in enumerate(indices, 1):
        b = find_confidence_band(analysis, idx, **kwargs)
        out.append(b)
        if verbose:
            print(f"[{k}/{len(indices)}] idx {idx} "
                  f"(v={b['vmin_mid']:.0f} km/s)  {_band_summary(b)}  "
                  f"({b['n_evaluations']} evals)")
    return out


#: sigma label for the conventional confidence levels (else the raw number).
_SIGMA_LABEL = {0.68: "1sigma", 0.683: "1sigma", 0.90: "90%",
                0.95: "2sigma", 0.954: "2sigma"}


def _band_summary(b: dict) -> str:
    """One-line physical-units (cm^-1) summary of a band record."""
    from .constants import CM           # eta/flux are 1/length: x*CM = cm^-1
    parts = [f"best {b['best_fit'] * CM:.3g} cm^-1"]
    for lv in b["levels"]:
        lo, hi = b["band"][lv]
        label = _SIGMA_LABEL.get(round(lv, 3), f"{lv:g}")
        hi_s = f"{hi * CM:.3g}" if np.isfinite(hi) else "inf"
        parts.append(f"{label} [{lo * CM:.3g}, {hi_s}]")
    return "  ".join(parts)
