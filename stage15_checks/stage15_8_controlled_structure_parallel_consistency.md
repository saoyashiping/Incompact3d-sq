# Stage 15.8: np=1/2/4 controlled structure-advance parallel consistency

## Scope

Stage 15.8 validates that the controlled Stage 15 structure-state update and the Stage 15.7 feedback-linkage diagnostics are decomposition-consistent when the standalone linkage executable is launched with MPI process counts:

```text
np = 1, 2, 4
```

The controlled update remains the Stage 15.6 candidate formula:

```text
A_f_cand = F_test / rho_tilde
V_f_new  = V_f_old + dt * A_f_cand
X_f_new  = X_f_old + dt * V_f_new
```

The feedback linkage remains the Stage 15.7 relation:

```text
F_fs_cand = alpha * (U_f - V_f)
Delta slip = -Delta V_f
Delta F_fs_cand = -alpha * Delta V_f
```

Stage 15.8 is a parallel-consistency validation stage only. It does not activate full production FSI physics, bending, tension, wall/contact logic, multi-fibre logic, direct Stage 14 RHS injection, pressure/projection, Poisson, RK3, or channel-forcing changes.

## Files

Stage 15.8 adds:

```text
stage15_checks/run_stage15_8_controlled_structure_parallel_consistency.sh
stage15_checks/stage15_8_controlled_structure_parallel_consistency.md
```

No new Fortran source is required. The wrapper reuses the Stage 15.7 standalone check target:

```text
fibre_stage15_feedback_linkage_check
```

## Manual command

```bash
bash stage15_checks/run_stage15_8_controlled_structure_parallel_consistency.sh
```

The wrapper configures `BUILD_DIR` if needed, builds the Stage 15.7 feedback-linkage target, runs it under `MPIEXEC -n 1`, `MPIEXEC -n 2`, and `MPIEXEC -n 4`, stores per-np diagnostics, compares the np=2 and np=4 results against np=1, and writes a Stage 15.8 summary file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher used for np=1/2/4 runs. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_8_RUN_STAGE15_7` | `0` | Optionally run the Stage 15.7 wrapper as prerequisite evidence. |
| `STAGE15_8_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_8_REQUIRE_STAGE15_7` | `1` | Require Stage 15.7 wrapper/check/doc evidence. |
| `STAGE15_8_ENABLE` | `1` | Enable the standalone controlled linkage validation. |
| `STAGE15_8_CONTROLLED_STEP_ENABLE` | `1` | Enable the controlled velocity-update path. |
| `STAGE15_8_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only this controlled parallel-consistency check. |
| `STAGE15_8_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_8_NP_LIST` | `1 2 4` | Required MPI process-count list. |
| `STAGE15_8_NPTS` | `8` | Number of Lagrangian/control points. |
| `STAGE15_8_DT` | `1.0e-4` | Controlled step size. |
| `STAGE15_8_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_8_TEST_FORCE` | `1.0e-6` | Small controlled force magnitude. |
| `STAGE15_8_FEEDBACK_ALPHA` | `1.0` | Stage 12 feedback gain. |
| `STAGE15_8_LAMBDA` | `1.0e-8` | Small nonzero Stage 14 diagnostic response scale. |
| `STAGE15_8_MAX_VELOCITY_UPDATE` | `1.0e-9` | Maximum allowed controlled velocity update per run. |
| `STAGE15_8_MAX_SLIP_ERROR` | `1.0e-14` | Maximum allowed slip-consistency error per run. |
| `STAGE15_8_MAX_FORCE_ERROR` | `1.0e-14` | Maximum allowed feedback-force consistency error per run. |
| `STAGE15_8_MAX_FORCE_RESPONSE` | `1.0e-8` | Maximum allowed force response per run. |
| `STAGE15_8_MAX_RHS_RESPONSE` | `1.0e-12` | Maximum allowed RHS response per run. |
| `STAGE15_8_MAX_PARALLEL_VELOCITY_DIFF` | `1.0e-14` | Maximum allowed velocity-response difference versus np=1. |
| `STAGE15_8_MAX_PARALLEL_SLIP_DIFF` | `1.0e-14` | Maximum allowed slip-response difference versus np=1. |
| `STAGE15_8_MAX_PARALLEL_FORCE_DIFF` | `1.0e-14` | Maximum allowed force-response difference versus np=1. |
| `STAGE15_8_MAX_PARALLEL_RHS_DIFF` | `1.0e-14` | Maximum allowed RHS-response difference versus np=1. |
| `STAGE15_8_RUN_PRODUCTION_SMOKE` | `0` | Production smoke is deferred by default; MPI-launched standalone linkage validation is reported explicitly. |

## Static audits

The Stage 15.8 wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- the Stage 15.8 wrapper or documentation is missing;
- the required Stage 15.7 feedback-linkage check target or source is missing;
- `xcompact3d.f90` loses the approved guarded Stage 15 diagnostic hook insertion;
- a Stage 15.8 production linkage appears without an approved guard;
- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, Stage 14.5, or Stage 15.1-15.7 diagnostics are missing;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 force-density sampling repair evidence is missing;
- `stage14_checks/STAGE14_CLOSED.md` is missing when required;
- `STAGE15_8_NP_LIST` differs from the required `1 2 4` list.

## Runtime diagnostics

The wrapper stores per-np diagnostics as:

```text
stage15_outputs/stage15_8_np1_feedback_linkage.dat
stage15_outputs/stage15_8_np2_feedback_linkage.dat
stage15_outputs/stage15_8_np4_feedback_linkage.dat
```

It writes the required Stage 15.8 summary file:

```text
stage15_outputs/fibre_stage15_8_controlled_structure_parallel_consistency.dat
```

The summary contains at least:

```text
stage15_8_requested_status
np_list
np1_run_status
np2_run_status
np4_run_status
np1_final_status
np2_final_status
np4_final_status
np1_max_velocity_update
np2_max_velocity_update
np4_max_velocity_update
np1_max_slip_change
np2_max_slip_change
np4_max_slip_change
np1_max_feedback_force_change
np2_max_feedback_force_change
np4_max_feedback_force_change
parallel_velocity_diff_np2
parallel_velocity_diff_np4
parallel_slip_diff_np2
parallel_slip_diff_np4
parallel_force_diff_np2
parallel_force_diff_np4
parallel_rhs_diff_np2
parallel_rhs_diff_np4
parallel_velocity_status
parallel_slip_status
parallel_force_status
parallel_rhs_status
controlled_update_count_consistency_status
no_rank_corruption_status
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

The wrapper also writes a gate-status file:

```text
stage15_outputs/stage15_8_controlled_structure_parallel_consistency.dat
```

## PASS evidence

A successful run prints:

```text
STAGE 15.8 CONTROLLED STRUCTURE PARALLEL CONSISTENCY VERDICT: PASS
STAGE 15.8 FINAL VERDICT: PASS
```

and the summary file reports:

- `np1_run_status 1`, `np2_run_status 1`, and `np4_run_status 1`
- `np1_final_status 1`, `np2_final_status 1`, and `np4_final_status 1`
- `parallel_velocity_status 1`
- `parallel_slip_status 1`
- `parallel_force_status 1`
- `parallel_rhs_status 1`
- `controlled_update_count_consistency_status 1`
- `no_rank_corruption_status 1`
- all full-production advance, bending, tension, wall/contact, multi-fibre, and direct RHS-injection connection counts equal to zero
- no fluid/RHS, pressure/projection, Poisson, or RK3/channel-forcing modification statuses equal to one
- `final_status 1`

## Assumptions and risks

- Stage 15.8 intentionally reuses the Stage 15.7 standalone check and launches it under MPI process counts 1, 2, and 4. The logical controlled-linkage check remains deterministic and should produce identical diagnostics across launches.
- The wrapper detects obvious rank-corrupted diagnostics by requiring a single `final_status` entry and all required keys in each per-np file.
- `STAGE15_8_RUN_PRODUCTION_SMOKE=0` is the safe default. The wrapper records production smoke as deferred rather than silently claiming that a production smoke was run.
