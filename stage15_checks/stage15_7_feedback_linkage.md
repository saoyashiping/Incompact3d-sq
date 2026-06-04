# Stage 15.7: structure advance and Stage 12 feedback-force linkage validation

## Scope

Stage 15.7 validates, in a standalone diagnostic check, that a controlled Stage 15 structure-owned velocity update changes the Stage 12 feedback-force candidate through the expected slip-velocity linkage:

```text
F_fs_cand = alpha * (U_f - V_f)
Delta slip = -Delta V_f
Delta F_fs_cand = -alpha * Delta V_f
```

This is a linkage-validation stage only. It does not activate full production FSI physics, bending, tension, wall/contact logic, multi-fibre logic, direct Stage 14 RHS injection, pressure/projection, Poisson, RK3, or channel-forcing changes.

## Files

Stage 15.7 adds:

```text
src/fibre_stage15_feedback_linkage_check.f90
stage15_checks/run_stage15_7_feedback_linkage.sh
stage15_checks/stage15_7_feedback_linkage.md
```

The only build-system change is the standalone check target:

```text
fibre_stage15_feedback_linkage_check
```

The check links the existing Stage 12 feedback formula module and calls `stage12_feedback_formula_compute_controlled` directly, so the feedback-force response is not a divergent reimplementation.

## Manual command

```bash
bash stage15_checks/run_stage15_7_feedback_linkage.sh
```

The wrapper configures `BUILD_DIR` if needed, builds only the Stage 15.7 feedback-linkage check target, runs it with `MPIEXEC -n 1`, performs static regression audits, and validates the generated diagnostic file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher for the standalone check. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_7_RUN_STAGE15_6` | `0` | Optionally run the Stage 15.6 wrapper as prerequisite evidence. |
| `STAGE15_7_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_7_REQUIRE_STAGE15_6` | `1` | Require Stage 15.6 wrapper/check/doc evidence. |
| `STAGE15_7_ENABLE` | `1` | Enable the standalone Stage 15.7 linkage check. |
| `STAGE15_7_CONTROLLED_STEP_ENABLE` | `1` | Enable the controlled velocity update. |
| `STAGE15_7_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only this controlled linkage-check update. |
| `STAGE15_7_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_7_NP` | `1` | Required process count for this stage. |
| `STAGE15_7_NPTS` | `8` | Number of Lagrangian/control points. |
| `STAGE15_7_DT` | `1.0e-4` | Controlled step size. |
| `STAGE15_7_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_7_TEST_FORCE` | `1.0e-6` | Small controlled force magnitude. |
| `STAGE15_7_FEEDBACK_ALPHA` | `1.0` | Stage 12 feedback gain used in the standalone formula call. |
| `STAGE15_7_LAMBDA` | `1.0e-8` | Small nonzero Stage 14 diagnostic response scale. |
| `STAGE15_7_MAX_VELOCITY_UPDATE` | `1.0e-9` | Maximum allowed controlled velocity update. |
| `STAGE15_7_MAX_SLIP_ERROR` | `1.0e-14` | Maximum allowed error in `Delta slip = -Delta V_f`. |
| `STAGE15_7_MAX_FORCE_ERROR` | `1.0e-14` | Maximum allowed error in `Delta F = -alpha*Delta V_f`. |
| `STAGE15_7_MAX_FORCE_RESPONSE` | `1.0e-8` | Maximum allowed feedback-force response. |
| `STAGE15_7_MAX_RHS_RESPONSE` | `1.0e-12` | Maximum allowed bounded RHS diagnostic response. |
| `STAGE15_7_RUN_PRODUCTION_SMOKE` | `0` | Production smoke is deferred by default; standalone linkage validation is reported explicitly. |

## Static audits

The Stage 15.7 wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- Stage 15.7 source calls production bending, tension, wall/contact, multi-fibre, pressure/projection/Poisson/RK3, channel-forcing, or production time-loop routines;
- a direct Stage 14 module/source connection is introduced outside the approved Stage 12-13-14 diagnostic chain;
- `xcompact3d.f90` loses the approved guarded Stage 15 diagnostic hook insertion;
- a Stage 15.7 production-time-loop linkage appears without an approved guard;
- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, Stage 14.5, or Stage 15.1-15.6 diagnostics are missing;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 sampling-repair evidence is missing;
- `stage14_checks/STAGE14_CLOSED.md` is missing when required.

## Runtime diagnostics

The required Stage 15.7 diagnostic file is:

```text
stage15_outputs/fibre_stage15_7_feedback_linkage.dat
```

It contains at least:

```text
stage15_7_requested_status
controlled_step_enabled_status
structure_advance_enable_status
diagnostic_only_status
np
npts
dt
rho_tilde
test_force_magnitude
feedback_alpha
lambda_value
max_velocity_update
velocity_update_nonzero_status
max_slip_change
slip_change_nonzero_status
slip_error
slip_consistency_status
max_feedback_force_change
feedback_force_change_nonzero_status
feedback_force_error
feedback_force_consistency_status
force_response_bounded_status
rhs_response_bounded_status
controlled_update_count
production_full_structure_advance_count
bending_solve_count
tension_solve_count
wall_contact_count
multifibre_count
rhs_injection_connection_count
approved_stage12_13_14_chain_status
no_fluid_rhs_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
final_status
```

The wrapper also writes:

```text
stage15_outputs/stage15_7_feedback_linkage.dat
```

## PASS evidence

A successful run prints:

```text
STAGE 15.7 FEEDBACK LINKAGE VERDICT: PASS
STAGE 15.7 FINAL VERDICT: PASS
```

and the required diagnostic file reports:

- `stage15_7_requested_status 1`
- `controlled_step_enabled_status 1`
- `structure_advance_enable_status 1`
- `diagnostic_only_status 1`
- `np 1`
- `velocity_update_nonzero_status 1`
- `slip_change_nonzero_status 1`
- `slip_consistency_status 1`
- `feedback_force_change_nonzero_status 1`
- `feedback_force_consistency_status 1`
- `force_response_bounded_status 1`
- `rhs_response_bounded_status 1`
- `controlled_update_count 1`
- all full-production advance, bending, tension, wall/contact, multi-fibre, and direct RHS-injection connection counts equal to zero
- no fluid/RHS, pressure/projection, Poisson, or RK3/channel-forcing modification statuses equal to one
- `final_status 1`

## Assumptions and risks

- Stage 15.7 intentionally uses a standalone diagnostic linkage check and does not add production time-loop coupling.
- Stage 13 and Stage 14 responses are represented as bounded diagnostic responses from the approved Stage 12 feedback-force change and the small nonzero `lambda` scale; direct Stage 14 injection is not called in this standalone stage.
- `STAGE15_7_RUN_PRODUCTION_SMOKE=0` is the safe default. The wrapper records production smoke as deferred instead of silently claiming that a production smoke was run.
