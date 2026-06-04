# Stage 15.6: controlled single-step structure advance at np=1

## Scope

Stage 15.6 adds a standalone, tightly guarded controlled single-step structure-state update check at `np=1`. It uses the already checked Stage 15.3 explicit candidate formula:

```text
A_f_cand = F_test / rho_tilde
V_f_new  = V_f_old + dt * A_f_cand
X_f_new  = X_f_old + dt * V_f_new
```

The stage verifies that a very small test force produces finite, bounded, sign-consistent, diagnostically auditable `X_f`/`V_f`/`A_f` candidate updates and exactly one controlled update count. It does not activate full production flexible-fibre physics, bending, tension, wall/contact handling, multi-fibre handling, Stage 14 RHS injection coupling, or solver-path modifications.

## Files

Stage 15.6 adds:

```text
src/fibre_stage15_controlled_structure_step.f90
src/fibre_stage15_controlled_structure_step_check.f90
stage15_checks/run_stage15_6_controlled_structure_step_np1.sh
stage15_checks/stage15_6_controlled_structure_step_np1.md
```

The build registration is the standalone CMake target:

```text
fibre_stage15_controlled_structure_step_check
```

No Stage 10-15.5 source, script, or documentation behavior is intentionally changed except for minimal build registration of the new standalone check target.

## Manual command

```bash
bash stage15_checks/run_stage15_6_controlled_structure_step_np1.sh
```

The wrapper configures `BUILD_DIR` if missing, builds only the Stage 15.6 check target, runs the standalone check through `MPIEXEC -n 1`, audits required regression markers, and validates the Stage 15.6 diagnostic file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher for the standalone np=1 check. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_6_RUN_STAGE15_5` | `0` | Optionally run the Stage 15.5 wrapper as prerequisite evidence. |
| `STAGE15_6_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_6_REQUIRE_STAGE15_5` | `1` | Require Stage 15.5 wrapper/documentation evidence. |
| `STAGE15_6_ENABLE` | `1` | Enable Stage 15 for the standalone controlled-step check. |
| `STAGE15_6_CONTROLLED_STEP_ENABLE` | `1` | Enable the controlled single-step path. |
| `STAGE15_6_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only this controlled Stage 15.6 single-step update. |
| `STAGE15_6_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_6_NP` | `1` | Required process count for Stage 15.6. |
| `STAGE15_6_NPTS` | `8` | Number of controlled points in the standalone check. |
| `STAGE15_6_DT` | `1.0e-4` | Controlled time-step size. |
| `STAGE15_6_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_6_TEST_FORCE` | `1.0e-6` | Small test-force magnitude. |
| `STAGE15_6_MAX_POSITION_UPDATE` | `1.0e-12` | Maximum allowed position update. |
| `STAGE15_6_MAX_VELOCITY_UPDATE` | `1.0e-9` | Maximum allowed velocity update. |
| `STAGE15_6_MAX_ACCELERATION` | `1.0e-5` | Maximum allowed candidate acceleration. |
| `STAGE15_6_MAX_ZERO_FORCE_DRIFT` | `1.0e-14` | Maximum allowed zero-force drift. |
| `STAGE15_6_MAX_FLUID_DELTA` | `0.0` | Recorded no-fluid-modification tolerance. |
| `STAGE15_6_RUN_PRODUCTION_SMOKE` | `0` | Production smoke is deferred by default; standalone controlled-step validation is reported explicitly. |

## Static audits

The wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- Stage 15.6 source files call production bending, tension, wall/contact, multi-fibre, Stage 14 RHS injection, fluid RHS, pressure/projection/Poisson/RK3, or channel-forcing routines;
- `xcompact3d.f90` loses the approved Stage 15 guarded diagnostic hook insertion;
- a Stage 15.6 production-time-loop connection appears without an approved guard;
- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, Stage 14.5, or Stage 15.1-15.5 diagnostics are missing;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 force-density sampling repair evidence is missing;
- `stage14_checks/STAGE14_CLOSED.md` is missing when required.

## Runtime diagnostics

The required Stage 15.6 diagnostic file is:

```text
stage15_outputs/fibre_stage15_6_controlled_structure_step_np1.dat
```

It contains at least:

```text
stage15_6_requested_status
controlled_step_enabled_status
structure_advance_enable_status
diagnostic_only_status
np
npts
dt
rho_tilde
test_force_magnitude
zero_force_drift
zero_force_status
small_force_status
max_acceleration
max_velocity_update
max_position_update
forced_component_nonzero_status
sign_consistency_status
bounded_update_status
controlled_update_count
production_full_structure_advance_count
bending_solve_count
tension_solve_count
wall_contact_count
multifibre_count
rhs_injection_connection_count
no_fluid_rhs_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
final_status
```

The wrapper also writes a gate summary:

```text
stage15_outputs/stage15_6_controlled_structure_step_np1.dat
```

## PASS evidence

A successful run prints:

```text
STAGE 15.6 CONTROLLED STRUCTURE STEP NP1 VERDICT: PASS
STAGE 15.6 FINAL VERDICT: PASS
```

and the required diagnostic file reports:

- `stage15_6_requested_status 1`
- `controlled_step_enabled_status 1`
- `structure_advance_enable_status 1`
- `diagnostic_only_status 1`
- `np 1`
- `zero_force_status 1`
- `small_force_status 1`
- `forced_component_nonzero_status 1`
- `sign_consistency_status 1`
- `bounded_update_status 1`
- `controlled_update_count 1`
- all full-production advance, bending, tension, wall/contact, multi-fibre, and RHS-injection connection counts equal to zero
- no fluid/RHS, pressure/projection, Poisson, or RK3/channel-forcing modification statuses equal to one
- `final_status 1`

## Assumptions and risks

- The controlled single-step update is standalone and diagnostic-only; it is not wired into the production time loop in this stage.
- `STAGE15_6_RUN_PRODUCTION_SMOKE=0` is the safe default because the Stage 15.5 production no-contamination path is still the authoritative production no-op validation; the wrapper records this as deferred rather than silently pretending a production smoke was run.
- The update magnitudes are intentionally tiny under the default force and time step, so tightening any bound below the expected formula result will correctly fail the stage.
