# box_response_test — isolated best-fit sandbox

**Not production.** Nothing here is imported by `quantum_sensor`, writes into
`qsensor_analysis/results/`, or is touched by the CI notebooks. It only
*reads* the production solver (`run_optimize_qp`), conditioner (`condition`)
and halo model (`eta_model`) to test the **best-fit estimator itself** on
synthetic response matrices.

## Question

After the response matrix is windowed (the per-row `1-alpha` tail cut), the
**last best-fit step undershoots the input η** — it does not meet the halo
curve at the high-`v_min` edge of the window. Is that undershoot

1. **intrinsic to the best-fit** (the minimal-total-flux monotone vertex), or
2. **caused by the decaying tail** of the response (the windowed columns carry
   tiny, decaying weight, so the minimal-total-flux objective collapses them)?

## Test

`box_test.py` builds tail-free **box** response matrices (response is flat
above each energy bin's threshold, then a hard cut — no decay) and compares
the recovered staircase to a known monotone η on a grid that **ends while η is
still positive** (mirroring the windowed matrix). A **decaying-tail** matrix of
the same geometry is run alongside as the control.

- `box_flat`  — R_i(v) = 1 for v ≥ thr_i, hard 0 past the grid edge
- `box_window`— R_i(v) = 1 on a finite box [thr_i, thr_i+W]
- `tail`      — R_i(v) = exp(-(v-thr_i)/L) for v ≥ thr_i  (control: the real shape)

Run:

```
python box_test.py
```

Writes `box_test.pdf` and a printed last-steps report here in the sandbox only.

## Findings

In all three matrices the true η reproduces the data to machine precision
(`||M x − M η|| / ||M η|| ~ 1e-16`): **η is always on the χ²=0 face**, so the
last-step undershoot is never a fit-quality problem — it is which point on the
degenerate optimal face the **minimal-total-flux** vertex tie-break picks.

| matrix | last steps x/η | Σx / Ση | behaviour |
|---|---|---|---|
| `box_window` (localized) | ≈ 1 | 0.98 | **last step recovered** — tracks η to the edge |
| `box_flat` (flat, extended) | → 0 | 1.00 | collapses (degenerate optimum, equal total flux) |
| `tail` (decaying) | → 0 | 0.95 | collapses — **strictly lower total flux** |

**Conclusion.** The undershoot is **not intrinsic to the best-fit** — with a
response that localizes the high-`v_min` information (`box_window`) the
minimal-total-flux vertex tracks η to the window edge. It appears when the
highest-`v_min` columns are not independently pinned:

- `tail` (the production case): the decaying columns carry almost no count
  information, so **zeroing them strictly lowers the total flux** (Σx/Ση =
  0.95). Minimal-total-flux actively prefers the collapsed vertex over η.
- `box_flat`: equal total flux (1.00) — a degenerate optimum; the simplex just
  picks a sparse vertex arbitrarily.

So the production undershoot is a **minimal-total-flux regularization artifact
on the windowed (decaying) tail**, not a failure to fit. Cure directions (for
production, not done here): keep the highest-`v_min` response localized rather
than a decaying tail, or replace the minimal-total-flux tie-break with one that
does not reward zeroing the weakly-constrained tail (anchoring / smoothness /
max-entropy).
