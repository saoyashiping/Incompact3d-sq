# Stage 15.10: RHS / IBM / structure contamination audit under controlled Stage 15 update

## Scope

Stage 15.10 is a static plus runtime contamination audit for the controlled Stage 15 structure update and Stage 12 feedback-force linkage. It verifies that the only approved coupling route remains:

```text
Stage 15 controlled structure-state update
-> Stage 12 feedback-force candidate
-> Stage 13 Eulerian force-density candidate
-> Stage 14 controlled RHS diagnostic/injection path
```

It does not add new physics and does not activate full production bending, tension, wall/contact handling, multi-fibre logic, legacy production IBM forcing, pressure/projection, Poisson, RK3, channel forcing, or unrelated DNS numerics.

The controlled structure update remains:

```text
A_f_cand = F_test / rho_tilde
V_f_new  = V_f_old + dt * A_f_cand
X_f_new  = X_f_old + dt * V_f_new
```

The feedback-force linkage remains:

```text
F_fs_cand = alpha * (U_f - V_f)
```

## Files

Stage 15.10 adds:

```text
stage15_checks/run_stage15_10_rhs_ibm_structure_contamination_audit.sh
stage15_checks/stage15_10_rhs_ibm_structure_contamination_audit.md
```

No new Fortran source is required. The wrapper reuses the existing Stage 15.7 feedback-linkage check target:

```text
fibre_stage15_feedback_linkage_check
```

and the existing Stage 14.10 contamination audit wrapper:

```text
stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh
```

## Manual command

```bash
bash stage15_checks/run_stage15_10_rhs_ibm_structure_contamination_audit.sh
```

The wrapper configures `BUILD_DIR` if needed, builds the Stage 15.7 linkage target, performs static anti-contamination audits, runs the controlled Stage 15 linkage audit under MPI, runs the inherited Stage 14.10 production contamination audit when production smoke is enabled, and writes a Stage 15.10 summary file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_10_RUN_STAGE15_9` | `0` | Optionally run the Stage 15.9 prerequisite wrapper. |
| `STAGE15_10_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_10_REQUIRE_STAGE15_9` | `1` | Require Stage 15.9 wrapper/documentation evidence. |
| `STAGE15_10_ENABLE` | `1` | Enable the Stage 15.10 controlled contamination audit. |
| `STAGE15_10_CONTROLLED_STEP_ENABLE` | `1` | Enable the controlled Stage 15 step in the standalone linkage audit. |
| `STAGE15_10_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only this controlled contamination audit. |
| `STAGE15_10_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_10_NP` | `2` | MPI process count for the Stage 15 linkage launch and Stage 14.10 audit path. |
| `STAGE15_10_NPTS` | `8` | Number of controlled structure points. |
| `STAGE15_10_DT` | `1.0e-4` | Controlled step size. |
| `STAGE15_10_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_10_TEST_FORCE` | `1.0e-6` | Small controlled test force. |
| `STAGE15_10_FEEDBACK_ALPHA` | `1.0` | Feedback-force gain. |
| `STAGE15_10_LAMBDA` | `1.0e-8` | Small nonzero Stage 14 diagnostic response scale. |
| `STAGE15_10_MAX_FORCE_RESPONSE` | `1.0e-8` | Maximum allowed Stage 15 feedback-force response. |
| `STAGE15_10_MAX_RHS_RESPONSE` | `1.0e-12` | Maximum allowed Stage 15 linkage RHS response. |
| `STAGE15_10_MAX_STAGE14_RHS_INCREMENT` | `1.0e-4` | Maximum allowed inherited Stage 14 RHS increment. |
| `STAGE15_10_MAX_FLUID_DELTA` | `1.0e-8` | Maximum allowed fluid signature delta from the production smoke audit. |
| `STAGE15_10_RUN_PRODUCTION_SMOKE` | `1` | Run the inherited Stage 14.10 production contamination audit. If disabled, the wrapper reports that production evidence was deferred instead of silently passing. |

## Static audits

The wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- the Stage 15.10 wrapper or documentation is missing;
- the Stage 15.9 wrapper/documentation evidence is missing when required;
- the Stage 15.7 feedback-linkage check target is missing;
- the Stage 14.10 contamination wrapper is missing;
- `xcompact3d.f90` loses the approved guarded Stage 15 hook or gains unapproved Stage 15.10/full-advance/IBM calls;
- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, Stage 14.5, or Stage 15.1-15.9 diagnostic markers are missing;
- old Stage 13.5 force-density diagnostic names reappear;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 force-density sampling repair evidence is missing;
- active Stage 15 source code calls or uses production pressure/projection/Poisson/RK3/channel forcing, production/legacy IBM forcing, bending/tension, wall/contact, multi-fibre, implicit/full structure solve, or direct Stage 14 RHS injection outside the approved chain.

## Runtime evidence

The wrapper writes the controlled Stage 15 linkage audit copy to:

```text
stage15_outputs/stage15_10_feedback_linkage.dat
```

It copies inherited Stage 14.10 contamination evidence to:

```text
stage15_outputs/stage15_10_stage14_rhs_ibm_structure_contamination_audit.dat
```

It writes the Stage 15.10 summary file:

```text
stage15_outputs/fibre_stage15_10_rhs_ibm_structure_contamination_audit.dat
```

The summary contains at least:

```text
stage15_10_requested_status
np
controlled_update_status
feedback_linkage_status
stage13_force_density_status
stage14_rhs_status
force_response_bounded_status
rhs_response_bounded_status
stage14_rhs_increment_bounded_status
fluid_signature_delta
fluid_signature_status
approved_stage12_13_14_chain_status
direct_rhs_injection_connection_count
legacy_ibm_forcing_count
production_ibm_forcing_outside_approved_chain_count
production_full_structure_advance_count
controlled_structure_update_count
bending_solve_count
tension_solve_count
implicit_structure_solve_count
wall_contact_count
multifibre_count
no_fluid_rhs_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
no_channel_forcing_modification_status
no_nan_inf_status
final_status
```

## PASS evidence

A successful run prints:

```text
STAGE 15.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: PASS
STAGE 15.10 FINAL VERDICT: PASS
```

and writes `final_status 1` in the Stage 15.10 summary file.

## Failure behavior

The wrapper fails if required static evidence, build evidence, controlled update evidence, Stage 12/13/14 chain evidence, Stage 14 RHS increment bounds, fluid signature bounds, no-NaN/Inf evidence, no-production-IBM evidence, no-full-structure-advance evidence, or no-solver-contamination evidence is missing or outside tolerance. Missing evidence is reported explicitly in the reasons block.
