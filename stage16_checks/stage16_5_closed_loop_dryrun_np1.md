# Stage 16.5: np=1 first closed-loop one-fibre dry run

Stage 16.5 adds a tightly guarded standalone diagnostic dry run for the first
np=1 one-fibre closed-loop path. It does **not** insert a production hook into
`xcompact3d.f90`, does **not** modify fluid RHS arrays directly, and does
**not** touch pressure/projection/Poisson/RK3/channel-forcing code.

## Manual command

```bash
bash stage16_checks/run_stage16_5_closed_loop_dryrun_np1.sh
```

The wrapper configures `BUILD_DIR` if needed, builds only
`fibre_stage16_closed_loop_dryrun_check`, runs the standalone check, and then
uses the Stage 16.5 helper to audit runtime diagnostics and static regression
evidence.

## Safe defaults

The wrapper supports these variables with fail-closed defaults:

- `BUILD_DIR=build_stage9`
- `MPIEXEC=mpirun`
- `MPIEXEC_FLAGS` may be empty
- `STAGE16_5_RUN_STAGE16_4=0`
- `STAGE16_5_REQUIRE_STAGE14_CLOSED=1`
- `STAGE16_5_REQUIRE_STAGE15_CLOSED=1`
- `STAGE16_5_REQUIRE_STAGE16_4=1`
- `STAGE16_5_ACCEPT_STAGE16_4_CLOSED_EVIDENCE=1`
- `STAGE16_5_ENABLE=1`
- `STAGE16_5_ONE_FIBRE_FSI_ENABLE=1`
- `STAGE16_5_STRUCTURE_UPDATE_ENABLE=1` for this controlled dry run only
- `STAGE16_5_TWO_WAY_RHS_ENABLE=1` for diagnostics only
- `STAGE16_5_DIAGNOSTIC_ONLY=1`
- `STAGE16_5_NP=1`
- `STAGE16_5_NPTS=8`
- `STAGE16_5_FEEDBACK_ALPHA=1.0`
- `STAGE16_5_LAMBDA=1.0e-8`
- `STAGE16_5_DT=1.0e-4`
- `STAGE16_5_RHO_TILDE=1.0`
- `STAGE16_5_MAX_FORCE_INPUT=1.0e-6`
- `STAGE16_5_MAX_STRUCTURE_UPDATE=1.0e-12`
- `STAGE16_5_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_5_MAX_FLUID_DELTA=1.0e-8`

## Diagnostic path

The standalone check exercises a controlled analytic np=1 path:

1. deterministic Stage 11-style sampled `U_f` and structure `V_f` values,
2. Stage 12-style slip and feedback candidate,
3. Stage 16.4 structure-side fluid-on-fibre force input,
4. action-reaction-compatible fluid-side force signature,
5. a bounded diagnostic structure-update increment,
6. Stage 13-style force-density signature,
7. Stage 14-style RHS diagnostic increment, and
8. a no-contamination fluid signature.

The Stage 14 RHS value is a scalar diagnostic increment only. It is not a direct
fluid RHS write and is not connected to the production time loop.

## Output file

The check/helper writes:

```text
stage16_outputs/fibre_stage16_5_closed_loop_dryrun_np1.dat
```

Required PASS evidence includes `final_status 1` plus these status keys:

- `stage16_5_requested_status`
- `np`
- `one_fibre_count_status`
- `closed_loop_path_status`
- `stage11_sampling_status`
- `stage12_feedback_status`
- `stage16_4_force_input_status`
- `structure_force_input_finite_status`
- `structure_force_input_bounded_status`
- `controlled_structure_update_status`
- `controlled_structure_update_bounded_status`
- `stage13_force_density_status`
- `stage14_rhs_status`
- `rhs_increment_bounded_status`
- `fluid_signature_delta`
- `fluid_signature_status`
- `approved_stage12_13_14_chain_status`
- `no_direct_rhs_injection_status`
- `no_production_hook_status`
- `no_pressure_projection_modification_status`
- `no_poisson_modification_status`
- `no_rk3_channel_forcing_modification_status`
- `no_channel_forcing_modification_status`
- `no_wall_contact_status`
- `no_multifibre_status`
- `no_legacy_ibm_forcing_status`
- `no_nan_inf_status`
- `stage14_regression_status`
- `stage15_regression_status`
- `stage16_1_regression_status`
- `stage16_2_regression_status`
- `stage16_3_regression_status`
- `stage16_4_regression_status`
- `final_status`

The wrapper/helper prints:

```text
STAGE 16.5 CLOSED-LOOP DRY RUN NP1 VERDICT: PASS
STAGE 16.5 FINAL VERDICT: PASS
```

or `FAIL` with explicit failure reasons.

## Static-audit false-positive policy

The Stage 16.5 helper intentionally reuses the corrected Stage 16.4 audit
pattern. It scans real `.f90` and shell-script logic, but does not scan markdown
documentation as executable regression evidence. It does not treat negative-check
strings or regex literals as production behavior, and it allows the legitimate
Stage 13.5 conservation/sign audit while still requiring current Stage 13.6
production force-density diagnostic evidence.
