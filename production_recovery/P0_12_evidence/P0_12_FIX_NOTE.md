# P0_12 fix note

This patch fixes the P0_12 zero-force proxy validation failure.

Observed failure:

```text
P0_12_SYNTHETIC_CLOSED_LOOP_CHECK FAIL: zero-force proxy run failed
```

Root cause:

`fibre_prod_synthetic_closed_loop_run` treated every `lambda_fsi > 0` run as requiring a nonzero RHS increment. That is correct for the normal small-lambda path with a nonzero force buffer, but it is wrong for the zero-force proxy case used by P0_12. In that proxy case `beta_hydro=0`, so the structure input, reaction force, Eulerian force buffer, and RHS increment must all remain zero while the run still returns `status=0`.

Fix:

The RHS increment check now first verifies the exact force-buffer scaling relation

```text
rhs_increment = lambda_fsi * penalty_beta * force_buffer
```

and only requires a nonzero RHS increment when `lambda_fsi > 0` and the corresponding force buffer is nonzero. A zero force buffer with zero RHS response is now accepted as the intended fail-closed/no-response proof.

No production DNS path, pressure/projection/RK3/channel-forcing logic, RHS adapter formula, force-buffer spreading formula, or previously closed P0_0-P0_11 evidence is modified.
