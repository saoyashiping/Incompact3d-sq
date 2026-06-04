# Stage 15.5: lambda=0 / structure-no-op production no-contamination validation

## Scope

Stage 15.5 validates a production no-contamination case for the guarded Stage 15 production structure hook. The hook is enabled in diagnostic-only/no-op mode while Stage 14 RHS injection uses `lambda=0`, and the wrapper compares the resulting production signatures against the corresponding no-op baseline.

This stage is a validation wrapper and documentation stage only. It does not add or activate structure advance, bending, tension, wall/contact handling, multi-fibre logic, or any solver-side coupling.

## Manual command

```bash
bash stage15_checks/run_stage15_5_structure_noop_invariance.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds `xcompact3d`, runs or verifies the baseline and no-op production cases, and writes Stage 15.5 diagnostics under `stage15_outputs/`.

## Environment variables

The wrapper supports these variables:

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher forwarded to production smoke runs. |
| `MPIEXEC_FLAGS` | empty | Optional launcher flags. |
| `STAGE15_5_RUN_STAGE15_4` | `0` | Optionally run the Stage 15.4 wrapper as prerequisite evidence. |
| `STAGE15_5_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_5_REQUIRE_STAGE15_4` | `1` | Require Stage 15.4 source/check/doc evidence. |
| `STAGE15_5_ENABLE` | `1` | Enable the Stage 15 hook in the no-op production case. |
| `STAGE15_5_STRUCTURE_ADVANCE_ENABLE` | `0` | Keep production structure advance disabled. |
| `STAGE15_5_DIAGNOSTIC_ONLY` | `1` | Require diagnostic-only Stage 15 behavior. |
| `STAGE15_5_LAMBDA` | `0.0` | Stage 14 RHS injection gain for the invariance case. |
| `STAGE15_5_NP` | `1` | Recorded Stage 15.5 process-count intent. The inherited Stage 9.9 smoke still checks its configured production decompositions. |
| `STAGE15_5_MAX_FLUID_DELTA` | `0.0` | Maximum baseline/no-op fluid signature delta. |
| `STAGE15_5_MAX_RHS_DELTA` | `0.0` | Maximum lambda-zero RHS increment/signature delta. |
| `STAGE15_5_MAX_STRUCTURE_DELTA` | `0.0` | Maximum structure-state delta inferred from no update counters. |
| `STAGE15_5_RUN_BASELINE` | `1` | Run the baseline case; when set to `0`, existing Stage 15.5 baseline evidence must already exist. |
| `STAGE15_5_RUN_PRODUCTION_SMOKE` | `1` | Run the Stage 15.5 no-op production smoke; when set to `0`, existing no-op evidence must already exist. |

## Runtime cases

### Baseline case

The baseline case enables the already closed Stage 11, Stage 12, Stage 13, and Stage 14 diagnostic chain with Stage 14 injection gain set to `STAGE15_5_LAMBDA=0.0`, while leaving the Stage 15 structure hook disabled or unset. The wrapper records the final production fluid signatures from `stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat` into:

```text
stage15_outputs/stage15_5_baseline_np1.dat
```

It also captures the Stage 14 lambda-zero RHS diagnostic as:

```text
stage15_outputs/stage15_5_baseline_stage14_rhs_hook.dat
```

### Stage 15.5 structure no-op case

The no-op case uses the same Stage 11-14 diagnostic chain with `lambda=0`, and additionally enables the Stage 15 production structure hook in diagnostic-only mode with structure advance disabled. The wrapper requires the production hook to be active and to report a positive apply count while all structure update/solve counters remain zero.

The wrapper records no-op evidence in:

```text
stage15_outputs/stage15_5_noop_np1.dat
stage15_outputs/stage15_5_noop_stage14_rhs_hook.dat
stage15_outputs/stage15_5_noop_production_structure_hook.dat
```

## Static audits

The Stage 15.5 wrapper fails closed if any required static evidence is missing or regressed. It checks that:

- the Stage 15.5 wrapper and documentation exist;
- the Stage 15.4 hook target and source are still registered;
- `xcompact3d.f90` retains the guarded Stage 15 diagnostic/no-op hook insertion and has no Stage 15 production structure-advance call;
- the old Stage 14 lambda-zero hook-registration gate is absent;
- Stage 11.5, Stage 13.6, Stage 14.5, and Stage 15.1-15.4 diagnostics still exist;
- rank0-safe diagnostic writing markers are preserved;
- Stage 13 force-density sampling repair evidence is preserved;
- `stage14_checks/STAGE14_CLOSED.md` exists when required.

The wrapper uses portable `grep`/`awk` checks and does not depend on any optional search utility.

## Required diagnostic output

The required Stage 15.5 diagnostic file is:

```text
stage15_outputs/fibre_stage15_5_structure_noop_invariance.dat
```

It contains at least:

```text
stage15_5_requested_status
baseline_run_status
noop_run_status
hook_active_status
hook_apply_count
diagnostic_only_status
noop_status
lambda_value
rhs_increment_zero_status
fluid_signature_delta
rhs_signature_delta
structure_state_delta
structure_advance_count
x_position_update_count
v_velocity_update_count
a_acceleration_update_count
bending_solve_count
tension_solve_count
wall_contact_count
multifibre_count
no_fluid_rhs_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
final_status
```

The wrapper also writes a gate summary:

```text
stage15_outputs/stage15_5_structure_noop_invariance.dat
```

## PASS evidence

A successful run prints:

```text
STAGE 15.5 STRUCTURE NO-OP INVARIANCE VERDICT: PASS
STAGE 15.5 FINAL VERDICT: PASS
```

and the required diagnostic file reports:

- `baseline_run_status 1`
- `noop_run_status 1`
- `hook_active_status 1`
- `hook_apply_count` greater than zero when the production smoke is run
- `diagnostic_only_status 1`
- `noop_status 1`
- `rhs_increment_zero_status 1`
- `fluid_signature_delta` less than or equal to `STAGE15_5_MAX_FLUID_DELTA`
- `rhs_signature_delta` less than or equal to `STAGE15_5_MAX_RHS_DELTA`
- `structure_state_delta` less than or equal to `STAGE15_5_MAX_STRUCTURE_DELTA`
- all structure advance, position, velocity, acceleration, bending, tension, wall/contact, and multi-fibre counters equal to zero
- `final_status 1`

## Assumptions and risks

- The wrapper inherits the existing Stage 9.9 production smoke driver to produce deterministic fluid signatures and existing Stage 14 lambda-zero RHS diagnostics.
- With the strict default tolerances set to zero, any non-bitwise production signature drift is treated as contamination and fails closed.
- If `STAGE15_5_RUN_BASELINE=0` or `STAGE15_5_RUN_PRODUCTION_SMOKE=0`, the wrapper does not silently skip validation; it requires existing Stage 15.5 evidence files and fails if they are missing.

## Stage 15.5 static hook-guard audit note

The xcompact3d hook-guard audit must distinguish executable `call stage15_production_structure_hook_apply(...)` statements from the `use fibre_stage15_production_structure_hook, only: ...` import list. Only real runtime `call` statements require same-line `stage15_structure_hook_reg` guard evidence. The import-only symbol list is not an unguarded hook invocation and must not trigger `xcompact3d_stage15_hook_apply_not_guarded`.
