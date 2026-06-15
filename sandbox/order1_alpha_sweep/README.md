# order1_alpha_sweep — alpha sweep in an O(1) unit system (no column scaling)

**Not production.** Isolated sandbox. Combines two earlier sandbox results:

1. the **round-trip** unit change (`../central_area_cut_test/gev_roundtrip_test.py`):
   recompute every constant at a new GeV so that `eta ~ O(1)` instead of ~1e-31,
   which makes the problem well-conditioned **without** the column-scaling
   band-aid;
2. the **central-area cut** alpha sweep (`../central_area_cut_test/alpha_sweep.py`):
   keep the highest-density region around each bin's peak covering area `alpha`.

Here the matrix is built at the O(1) GeV (eta recomputed via the round-trip
function, M rescaled by its energy-dimension -1), and the minimal-total-flux
best-fit is solved with a **PLAIN LP** (scipy HiGHS, no column scaling, no
CONDITION_C). The production column-scaled solver is overlaid as a cross-check.

Config: TES (Al), heavy mediator (q0), M1 (10 MeV), SHM.

## Run

```
python order1_sweep.py
```

Writes `by_alpha/alpha_*.pdf` (per-alpha, y-axis in O(1) units) and
`overlay.pdf` here in the sandbox only.

## What it tests

Whether moving to the principled O(1) unit system (right GeV) and dropping the
column-scaling band-aid reproduces the same alpha behaviour seen before:
clean staircase tracking eta at small alpha, collapse re-appearing at large
alpha, and the minimal-total-flux last-step shave. (It does — the shape is
unit-system-invariant; O(1) just lets a plain solver reach it.)

## Findings

GeV chosen so |eta|_max = 1.0 (GeV_O1 = 1.21e39, via the round-trip recompute).

| alpha | window [km/s] | plain LP last x/eta | behaviour |
|---|---|---|---|
| 0.30 | 92-143 | 0.87 | tracks eta, last-step shave |
| 0.50 | 79-145 | 0.86 | tracks eta, last-step shave |
| 0.70 | 59-159 | 0.00 | collapse |
| 0.90 | 52-350 | 0.00 | collapse |
| 0.99 | 49-518 | 0.00 | collapse |

Identical to the native-scale sweep (`../central_area_cut_test/alpha_sweep.py`):
the alpha behaviour is **unit-system-invariant**. The plain LP (no column
scaling, no CONDITION_C) solves cleanly at every alpha because eta ~ O(1).

**Bonus — CONDITION_C is scale-specific.** Running the production column-scaled
solver on the same O(1) problem:

- conditioning constant `c ~ O(1)` (matched to the O(1) x-scale): matches the
  plain LP to machine precision (err ~ 1e-16 ... 1e-14);
- the native default `c = CONDITION_C = 1e-30`: **100% wrong** (err = 1.0) at
  every alpha.

So `CONDITION_C = 1e-30` is tuned to the GeV=1e9 magnitude (x ~ 1e-31), not a
scale-free choice. In a different unit system the conditioning constant must
move with it. This is the concrete sense in which "the GeV choice, not the
column-scaling band-aid, sets the scale".
