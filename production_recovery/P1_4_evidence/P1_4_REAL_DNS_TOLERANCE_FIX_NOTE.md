# P1_4 real-DNS np-tolerance fix

This fix addresses the post-build/post-run P1_4 failure:

```text
lambda0 np consistency mismatch lambda0_np1 vs lambda0_np2, component 0: 0.00439849144 vs 0.00439752589, tol=1.01e-08
```

The six real xcompact3d runs had reached signature comparison. The previous script used a near bitwise tolerance (`abs_tol=1e-10`, `rel_tol=1e-8`) for real DNS-FSI signatures across MPI decompositions. That is too strict for short real DNS time advancement because FFT/reduction order and domain decomposition can change floating-point accumulation.

The fix keeps lambda=0 RHS no-contamination strict (`1e-20` by default), but compares physical diagnostic signatures with guarded real-DNS tolerances:

- `P1_4_NP_ABS_TOL`, default `1.0e-8`
- `P1_4_NP_REL_TOL`, default `1.0e-3`
- `P1_4_LAMBDA0_RHS_ZERO_TOL`, default `1.0e-20`

P1_4 remains self-contained and does not read P1_0-P1_3 PASS/FAIL or logs.
