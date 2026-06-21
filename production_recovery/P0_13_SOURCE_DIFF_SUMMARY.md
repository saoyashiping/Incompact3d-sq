# P0_13 Source Diff Summary

## Modified files

- `src/fibre_prod_synthetic_closed_loop.f90`
- `production_recovery/P0_13_evidence/P0_13_FIX_NOTE.md`
- `production_recovery/P0_13_SOURCE_DIFF_SUMMARY.md`

## Summary

P0_13 failed only during the synthetic np=1/2/4 closure rerun because `fibre_prod_synthetic_closed_loop.f90` still required a nonzero RHS increment whenever `lambda_fsi > 0`.  That incorrectly rejects the intentional zero-force proxy case where `beta_hydro=0` makes the force buffer zero and the correct RHS response is exactly zero.

The check now requires nonzero RHS increment only when `lambda_fsi > 0` and the checked Eulerian force-buffer component is nonzero.  The lambda-gated equality `rhs_increment_x = lambda_fsi * penalty_beta * force_buffer%fx` remains enforced.

## Production status

This fix does not make paper-scale DNS production ready.  It only repairs the P0_13 closure validation predicate so that the already-designed zero-force no-response case can pass correctly.
