# full_pipeline_unit_test — GeV as a free parameter, no conditioning band-aids

**Not production.** Isolated sandbox. Reruns the full forward+inverse pipeline
in a unit system where **GeV is an external parameter**, with eta set to O(1)
via the round-trip, and with **both** numerical band-aids removed:

- no column scaling,
- no `CONDITION_C`.

It then checks that the **physical observables are invariant** vs the
production pipeline (`DarkMatterQuantumAnalysis`, which runs at GeV=1e9 with
column scaling + `CONDITION_C=1e-30`):

1. event counts entering chi^2 (dimensionless) — identical at every GeV;
2. the recovered flux in physical units, `x * CM` [cm^-1] — identical at every
   GeV where the plain solve is numerically viable.

The production response matrix (`a.rm.matrix`, GeV=1e9) is reused unchanged; the
sandbox only (a) rescales it by its energy-dimension (M ∝ GeV^-1), (b) recomputes
eta at the chosen GeV via the round-trip function (eta ∝ GeV^+1), (c) rebuilds
the exposure at that GeV, and (d) solves with a plain minimal-total-flux LP.

Config: TES (Al), heavy mediator (q0), M1 (10 MeV), SHM, background none.

## Run

```
python pipeline_invariance.py
```

Prints a per-GeV table (counts match, physical-flux match, solve status) and
writes `pipeline_invariance.pdf` (physical recovered flux over GeV) here only.

## Expected result

- GeV=1e9 (eta~1e-31): the plain solve **fails** (numbers below solver
  tolerance) — this is exactly why production needs the band-aids.
- GeV chosen so eta~O(1) (and above): the plain solve **succeeds** and returns
  the **same physical flux** as production, to machine precision — confirming
  the pipeline is physics-invariant and the band-aids are only a substitute for
  a well-chosen unit scale.
