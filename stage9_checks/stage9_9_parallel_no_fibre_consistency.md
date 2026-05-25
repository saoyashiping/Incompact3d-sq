# Stage 9.9 parallel no-fibre DNS consistency

Stage 9.9 checks short no-fibre/no-coupling production channel DNS across MPI decompositions (`np=1,2,4`) using a deterministic Stage 9.9-only initial condition.

- Stage 9.9 comes after Stage 9.8 restart I/O checks.
- Stage 9.9 remains no-fibre/no-coupling/no-IBM-injection.
- It runs real production time advance for `X3D_STAGE9_9_MAX_STEPS` complete outer steps (default `3`) and stops only after completed outer steps.

## Recorded global signatures
For each run, Stage 9.9 records and writes:

- `stage9_9_signature_sum_ux`, `stage9_9_signature_sum_uy`, `stage9_9_signature_sum_uz`
- `stage9_9_signature_max_ux`, `stage9_9_signature_max_uy`, `stage9_9_signature_max_uz`
- `stage9_9_signature_l2_ux`, `stage9_9_signature_l2_uy`, `stage9_9_signature_l2_uz`

Plus finite statuses for velocity/pressure/divergence/CFL/massflux and final local status.

## Tolerances
- `X3D_STAGE9_9_SIGNATURE_TOL` (default `1.0e-8`) for **initial** signatures (strict absolute check)
- `X3D_STAGE9_9_FINAL_SIGNATURE_ABS_TOL` (default `1.0e-6`) for **final** signatures
- `X3D_STAGE9_9_FINAL_SIGNATURE_REL_TOL` (default `1.0e-12`) for **final** signatures
- `X3D_STAGE9_9_DIVERGENCE_TOL` (default `1.0e-8`)
- `X3D_STAGE9_9_MASSFLUX_TOL` (default `1.0e-6`)

Raw fresh channel initialization can be decomposition-dependent, so Stage 9.9 does **not** rely on checkpoint-based cross-decomposition restart.

- Stage 9.9 applies a deterministic analytic channel profile from global y-index:
  - `ux = 4*eta*(1-eta)`, `uy = 0`, `uz = 0`, `p = 0`
  - `eta = (global_j - 1)/max(1,ny-1)`
- Dat key `stage9_9_decomposition_invariant_initial_state_status` is set only when this deterministic field is actually applied.
- Initial signatures are compared across `np=1/2/4` first.
- If initial signatures differ, Stage 9.9 fails before judging time-advance consistency.
- If initial signatures match, final signatures are compared with mixed tolerance:
  - pass iff `abs(delta) <= max(abs_tol, rel_tol*max(1.0,abs(reference)))`.
- This mixed final criterion accounts for floating-point reduction/order differences after real time advancement and does **not** relax initial-state consistency.

## Not tested yet
- Long-time statistical equivalence beyond short smoke steps.
- Detailed turbulence-spectrum equivalence.

## Manual run
`bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh`

Expected output includes:
- `STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS`
- `STAGE 9.9 FINAL VERDICT: PASS`
