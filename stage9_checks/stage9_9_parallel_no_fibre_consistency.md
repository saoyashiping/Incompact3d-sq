# Stage 9.9 parallel no-fibre DNS consistency

Stage 9.9 checks short no-fibre/no-coupling production channel DNS across MPI decompositions (`np=1,2,4`) with methodology guards.

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
- `X3D_STAGE9_9_SIGNATURE_TOL` (default `1.0e-8`)
- `X3D_STAGE9_9_DIVERGENCE_TOL` (default `1.0e-8`)
- `X3D_STAGE9_9_MASSFLUX_TOL` (default `1.0e-6`)

Raw cross-np signature comparison is only valid when initialization is decomposition-invariant.

- Dat key `stage9_9_decomposition_invariant_initial_state_status` controls whether cross-np signature deltas are compared.
- If this status is `0`, the gate intentionally fails with a clear diagnostic instead of reporting misleading signature-delta failures.
- This avoids claiming parallel inconsistency when the initial perturbation/state itself is decomposition-dependent.

## Not tested yet
- Long-time statistical equivalence beyond short smoke steps.
- Detailed turbulence-spectrum equivalence.

## Manual run
`bash stage9_checks/run_stage9_9_parallel_no_fibre_consistency.sh`

Expected output includes:
- `STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS`
- `STAGE 9.9 FINAL VERDICT: PASS`
