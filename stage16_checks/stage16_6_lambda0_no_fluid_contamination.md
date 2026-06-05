# Stage 16.6: lambda=0 closed-loop no-fluid-contamination validation

Stage 16.6 runs the already validated one-fibre closed-loop diagnostic path with
`lambda = 0`. The feedback, structure-force-input, controlled structure-update,
Stage 13 force-density, and Stage 14 RHS diagnostic statuses must be active, but
the final RHS increment and fluid signature delta must remain zero or within the
strict Stage 16.6 no-contamination tolerances.

Stage 16.6 is diagnostic-only. It does not insert a production hook into
`xcompact3d.f90`, does not directly modify fluid RHS arrays, and does not touch
pressure/projection/Poisson/RK3/channel-forcing code.

## Manual command

```bash
bash stage16_checks/run_stage16_6_lambda0_no_fluid_contamination.sh
```

The wrapper configures `BUILD_DIR` if needed, builds only
`fibre_stage16_lambda0_no_contamination_check`, runs the standalone check, and
then invokes the Stage 16.6 helper to audit runtime diagnostics and static
regression evidence.

## Safe defaults

The wrapper supports these variables with fail-closed defaults:

- `BUILD_DIR=build_stage9`
- `MPIEXEC=mpirun`
- `MPIEXEC_FLAGS` may be empty
- `STAGE16_6_RUN_STAGE16_5=0`
- `STAGE16_6_REQUIRE_STAGE14_CLOSED=1`
- `STAGE16_6_REQUIRE_STAGE15_CLOSED=1`
- `STAGE16_6_REQUIRE_STAGE16_5=1`
- `STAGE16_6_ACCEPT_STAGE16_5_CLOSED_EVIDENCE=1`
- `STAGE16_6_ENABLE=1`
- `STAGE16_6_ONE_FIBRE_FSI_ENABLE=1`
- `STAGE16_6_STRUCTURE_UPDATE_ENABLE=1` for controlled diagnostics only
- `STAGE16_6_TWO_WAY_RHS_ENABLE=1` for RHS hook diagnostics only
- `STAGE16_6_DIAGNOSTIC_ONLY=1`
- `STAGE16_6_NP=1`
- `STAGE16_6_NPTS=8`
- `STAGE16_6_FEEDBACK_ALPHA=1.0`
- `STAGE16_6_LAMBDA=0.0`
- `STAGE16_6_DT=1.0e-4`
- `STAGE16_6_RHO_TILDE=1.0`
- `STAGE16_6_MAX_FORCE_INPUT=1.0e-6`
- `STAGE16_6_MAX_STRUCTURE_UPDATE=1.0e-12`
- `STAGE16_6_MAX_RHS_INCREMENT=1.0e-14`
- `STAGE16_6_MAX_FLUID_DELTA=1.0e-14`

## Diagnostic path

The standalone check exercises this controlled path:

1. deterministic Stage 11-style sampled `U_f` and structure `V_f` values,
2. Stage 12-style slip and feedback candidate,
3. Stage 16.4 structure-side fluid-on-fibre force input,
4. action-reaction-compatible fluid-side force signature,
5. bounded controlled structure-state diagnostic response,
6. Stage 13-style force-density candidate signature,
7. Stage 14-style RHS diagnostic path with `lambda = 0`, and
8. pre/post fluid-signature comparison proving no contamination.

The Stage 14 RHS value is a scalar diagnostic increment only. Since `lambda=0`,
`rhs_increment_value` must be zero or within `STAGE16_6_MAX_RHS_INCREMENT`; the
fluid signature delta must be zero or within `STAGE16_6_MAX_FLUID_DELTA`.

## Output file

The check/helper writes:

```text
stage16_outputs/fibre_stage16_6_lambda0_no_fluid_contamination.dat
```

Required PASS evidence includes `final_status 1` plus these status keys:

- `stage16_6_requested_status`
- `np`
- `lambda_value`
- `lambda_zero_status`
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
- `rhs_increment_value`
- `rhs_increment_zero_status`
- `rhs_increment_bounded_status`
- `fluid_signature_delta`
- `fluid_signature_unchanged_status`
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
- `stage16_5_regression_status`
- `final_status`

The wrapper/helper prints:

```text
STAGE 16.6 LAMBDA0 NO-FLUID-CONTAMINATION VERDICT: PASS
STAGE 16.6 FINAL VERDICT: PASS
```

or `FAIL` with explicit failure reasons.

## Static-audit false-positive policy

The Stage 16.6 helper intentionally reuses the corrected Stage 16.5 / Stage 16.4
audit pattern. It scans real `.f90` and shell-script logic, but does not scan
markdown documentation as executable regression evidence. It does not treat
negative-check strings or regex literals as production behavior, and it allows
the legitimate Stage 13.5 conservation/sign audit while still requiring current
Stage 13.6 production force-density diagnostic evidence.
