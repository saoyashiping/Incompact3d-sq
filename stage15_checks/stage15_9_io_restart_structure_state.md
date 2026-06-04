# Stage 15.9: restart / stats / visu / coarse I/O compatibility with controlled structure state
## Stage 15.9 RHS tolerance separation

Stage 15.9 intentionally keeps two RHS-related tolerances separate:

- `STAGE15_9_MAX_RHS_RESPONSE` is the strict Stage 15.7 diagnostic linkage bound for `rhs_response = lambda * Delta F`.
- `STAGE15_9_MAX_STAGE14_RHS_INCREMENT` is the inherited Stage 14.9 production small-lambda RHS-increment bound used when Stage 15.9 reuses the Stage 14.9 restart/stats/visu/coarse-I/O smoke. Its default is `1.0e-4`, matching the existing Stage 14.9 controlled small-lambda production bound.

Do not pass the much tighter Stage 15.7 diagnostic RHS-response tolerance into Stage 14.9 as `STAGE14_9_MAX_RHS_INCREMENT`; doing so falsely fails valid Stage 14.9 small-lambda production I/O evidence with `stage15_9_rhs_response_exceeds_tolerance`.

## Scope

Stage 15.9 validates I/O and restart compatibility for the controlled Stage 15 structure-state diagnostics. It combines two already passed diagnostic paths without enabling full production FSI physics:

1. the Stage 15.7 feedback-linkage standalone check, which exercises the Stage 15.6 controlled `X_f/V_f/A_f` update and verifies the Stage 12-style feedback relation; and
2. the Stage 14.9 production I/O/restart/stats/visu/coarse compatibility path, which verifies the closed Stage 11→13→14 diagnostic/RHS chain through the existing production I/O and restart machinery.

The controlled structure update remains the Stage 15.6 candidate formula:

```text
A_f_cand = F_test / rho_tilde
V_f_new  = V_f_old + dt * A_f_cand
X_f_new  = X_f_old + dt * V_f_new
```

The feedback-force linkage remains:

```text
F_fs_cand = alpha * (U_f - V_f)
Delta slip = -Delta V_f
Delta F_fs_cand = -alpha * Delta V_f
```

Stage 15.9 is an I/O/restart compatibility validation only. It does not activate full production bending, tension, wall/contact handling, multi-fibre logic, direct Stage 14 RHS injection outside the approved Stage 12→13→14 chain, pressure/projection, Poisson, RK3, channel forcing, or unrelated DNS numerics.

## Files

Stage 15.9 adds:

```text
stage15_checks/run_stage15_9_io_restart_structure_state.sh
stage15_checks/stage15_9_io_restart_structure_state.md
```

No new Fortran source is required. The wrapper reuses the already registered Stage 15.7 check target:

```text
fibre_stage15_feedback_linkage_check
```

and the existing Stage 14.9 production I/O/restart wrapper:

```text
stage14_checks/run_stage14_9_io_restart_stats_visu_rhs_injection.sh
```

## Manual command

```bash
bash stage15_checks/run_stage15_9_io_restart_structure_state.sh
```

The wrapper configures `BUILD_DIR` if needed, builds the Stage 15.7 feedback-linkage target, runs the controlled linkage check under MPI, runs the production I/O/restart/statistics/visualization/coarse compatibility path when requested, copies supporting diagnostics into `stage15_outputs/`, and writes the Stage 15.9 summary file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_9_RUN_STAGE15_8` | `0` | Optionally run the Stage 15.8 prerequisite wrapper. |
| `STAGE15_9_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_9_REQUIRE_STAGE15_8` | `1` | Require Stage 15.8 wrapper/documentation evidence. |
| `STAGE15_9_ENABLE` | `1` | Enable Stage 15.9 controlled compatibility validation. |
| `STAGE15_9_CONTROLLED_STEP_ENABLE` | `1` | Enable the Stage 15 controlled step in the standalone linkage check. |
| `STAGE15_9_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only the controlled Stage 15.9 compatibility check. |
| `STAGE15_9_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_9_NP` | `2` | MPI process count used by the Stage 15 linkage launch and Stage 14.9 I/O/restart path. |
| `STAGE15_9_NPTS` | `8` | Number of controlled structure points for the Stage 15 linkage check. |
| `STAGE15_9_DT` | `1.0e-4` | Controlled step size. |
| `STAGE15_9_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_9_TEST_FORCE` | `1.0e-6` | Small controlled force magnitude. |
| `STAGE15_9_FEEDBACK_ALPHA` | `1.0` | Feedback-force gain. |
| `STAGE15_9_LAMBDA` | `1.0e-8` | Small nonzero Stage 14 diagnostic response scale. |
| `STAGE15_9_MAX_STRUCTURE_RESTART_DELTA` | `1.0e-14` | Maximum allowed deterministic structure-state restart delta. |
| `STAGE15_9_MAX_FLUID_RESTART_DELTA` | `1.0e-8` | Maximum allowed inherited fluid restart delta. |
| `STAGE15_9_MAX_IO_SIGNATURE_DELTA` | `1.0e-8` | Maximum allowed inherited I/O signature delta. |
| `STAGE15_9_MAX_FORCE_RESPONSE` | `1.0e-8` | Maximum allowed controlled feedback-force response. |
| `STAGE15_9_MAX_RHS_RESPONSE` | `1.0e-12` | Maximum allowed controlled RHS response. |
| `STAGE15_9_RUN_RESTART` | `1` | Require restart write/read/continuation evidence. |
| `STAGE15_9_RUN_STATS_VISU` | `1` | Require statistics and visualization evidence. |
| `STAGE15_9_RUN_COARSE_IO` | `1` | Require coarse I/O evidence when supported by the inherited Stage 14.9 path. |
| `STAGE15_9_RUN_PRODUCTION_SMOKE` | `1` | Run the inherited production I/O/restart smoke path. If set to `0`, the wrapper reports diagnostic-only deferral instead of pretending production evidence was collected. |

## Static audits

The wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- the Stage 15.9 wrapper or documentation is missing;
- the Stage 15.7 feedback-linkage check target is missing;
- the Stage 15.8 wrapper/documentation evidence is missing when required;
- `xcompact3d.f90` loses the approved guarded Stage 15 diagnostic hook insertion or gains an unapproved Stage 15.9/full-advance production call;
- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, Stage 14.5, or Stage 15.1-15.8 diagnostic markers are missing;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 force-density sampling repair evidence is missing;
- `stage14_checks/STAGE14_CLOSED.md` is missing when required;
- the wrapper appears to require `rg` or recursive `grep`.

## Runtime and I/O evidence

The wrapper writes controlled Stage 15 linkage diagnostics to:

```text
stage15_outputs/stage15_9_feedback_linkage.dat
```

It copies the inherited Stage 14.9 I/O/restart evidence to:

```text
stage15_outputs/stage15_9_stage14_9_io_restart_stats_visu_rhs_injection.dat
```

It then writes the required Stage 15.9 summary file:

```text
stage15_outputs/fibre_stage15_9_io_restart_structure_state.dat
```

The summary contains at least:

```text
stage15_9_requested_status
np
npts
dt
rho_tilde
test_force_magnitude
feedback_alpha
lambda_value
controlled_update_status
controlled_update_count
structure_state_finite_status
feedback_linkage_status
stage13_force_density_status
stage14_rhs_status
restart_write_status
restart_file_status
restart_read_status
restart_continuation_status
structure_restart_delta
structure_restart_status
fluid_restart_delta
fluid_restart_status
stats_output_status
stats_nonempty_status
visu_output_status
visu_nonempty_status
coarse_io_output_status
coarse_io_nonempty_status
io_signature_delta
io_signature_status
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

## PASS evidence

A successful run prints:

```text
STAGE 15.9 IO/RESTART STRUCTURE STATE VERDICT: PASS
STAGE 15.9 FINAL VERDICT: PASS
```

and writes `final_status 1` in the Stage 15.9 summary file.

## Failure behavior

The final verdict fails if required restart write/read/continuation evidence, restart files, stats/visu/coarse I/O evidence, controlled structure update diagnostics, Stage 12/13/14 chain diagnostics, static no-regression evidence, or no-contamination evidence are missing or outside tolerance. The wrapper reports explicit reasons and never silently treats missing evidence as PASS.
