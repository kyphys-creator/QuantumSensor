# central_area_cut_test — two-sided (peak-centred) response cut

**Not production.** Isolated sandbox: reads the original full-domain response
functions and the SHM model read-only, builds its own matrix, runs the
production best-fit. Writes only inside this folder.

## Idea

The production matrix windows each row **one-sided** (from its kinematic
threshold up to where the cumulative integral reaches `1-alpha`), which keeps a
**decaying high-v_min tail** — and that tail is what the minimal-total-flux
best-fit collapses (see `../box_response_test`).

Here we cut **both ends instead**: for each bin's response function
`R_bin(v_min)` keep only the **highest-density region around the peak that
covers `alpha = 50%` of the area**, zero everything else (both the low-v rise
and the high-v tail). Then build the matrix from the cut response functions and
run the best-fit.

## Config

TES (Al), heavy mediator (q0), M1 (10 MeV), SHM (Halo), 5 bins —
`Mathematica/output/TES/response_functions_csv/Al_q0M1_R5.csv`.

## Run

```
python area_cut_test.py
```

Writes `response_cut.pdf` (original vs central-50%-area-cut response functions)
and `best_fit.pdf` (eta vs recovered staircase, central-cut matrix vs the full
uncut matrix) here in the sandbox only.

## Follow-ups in this folder

- `alpha_sweep.py` — sweeps the area fraction alpha; saves per-alpha figures
  (`by_alpha/`, `by_alpha_novertex/`) and an overlay. `VERTEX_SELECT` toggles
  the minimal-total-flux vertex vs the OSQP interior solution.
- `gev_scale_test.py` — tests whether the natural-unit **GeV magnitude**, not
  the column-scaling band-aid, sets the conditioning.

### gev_scale findings

`eta` is proportional to GeV^1 (RHO_DM·SIGMA_E/mchi), `M` to GeV^-1, so a change
of the base unit GeV→G is a similarity transform (physics unchanged, magnitudes
move). Solving the minimal-total-flux fit with a **plain LP (no column scaling,
no CONDITION_C)** across GeV:

| GeV | \|eta\| | plain LP |
|---|---|---|
| 1e9 (native) … 1e24 | 8e-31 … 8e-16 | **fails** (HiGHS model error) |
| 1e27 and up | ≥ 8e-13 | **solves**, error 2.09e-2 |

So **the user's hypothesis is correct about conditioning**: at the native
GeV=1e9 the numbers (\|eta\|~1e-31) are too small for a plain solve; the fit is
only possible because column scaling rescues it. Choosing GeV so \|eta\| ≳ 1e-13
(GeV ≳ ~1e26) makes a plain solve work with **no scaling tricks at all**.

But two things are *not* a scale problem:
1. Once big enough, the plain-LP answer is **identical** to the column-scaled
   production answer (error 2.09e-2) — column scaling was doing its job, not
   distorting anything.
2. The **last-step undershoot (x/eta = 0.855) is scale-invariant** — it persists
   at every working GeV (`gev_scale_shape.pdf`: x/eta = 1.0 across the window,
   0.855 only on the last column). It is the minimal-total-flux tie-break, not
   conditioning, and cannot be fixed by changing GeV.

Bottom line: GeV magnitude sets *whether you can solve cleanly without
band-aids*; it does **not** change *which* monotone vertex is chosen, so it does
not remove the last-step undershoot (that needs a different tie-break).
