# vertex_weight_test — the minimal-flux objective weighting

**Not production.** Isolated sandbox. Keeps vertex selection ON but varies the
**weight** in its LP objective

    minimize  sum_j w_j x_j     s.t.  M x = mu,  x_j >= x_{j+1},  x >= 0

The production weight is **uniform** (`w_j = 1`, plain total flux Sum x_j),
which shaves the highest-v_min / lowest-response step (see
`../full_pipeline_unit_test`). This tries other weights and measures how well
the recovered staircase tracks eta, especially the last step.

Weights tested (in physical x-space):

| name | w_j | idea |
|---|---|---|
| `uniform` | 1 | production default |
| `width` | Δv_j | minimise the integrated flux ∫x dv |
| `colnorm` | ‖M[:,j]‖ | weight by response strength |
| `invcolnorm` | 1/‖M[:,j]‖ | make low-response columns expensive to shave |

Solved with a plain LP (no column scaling, no CONDITION_C) in an O(1) unit
system via the round-trip, so the comparison is purely about the objective
weight. TES Al q0 R5 SHM, M1/M2/M3.

## Run

```
python weight_test.py
```

Writes `weights_M{1,2,3}.pdf` (overlay of the weightings vs eta) and prints a
last-step / overall-error table here in the sandbox only.
