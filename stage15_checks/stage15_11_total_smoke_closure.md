# Stage 15.11: Stage 15 total smoke and closure

## Closed early-stage evidence handling

Stage 15.11 treats Stage 15.0-15.7 as closed historical stages. In a fresh unzip tree, their runtime `stage15_outputs/*.dat` files may be absent even though the scripts, documentation, source modules, and check targets are present and the user has already manually passed those stages. Therefore `STAGE15_11_ACCEPT_CLOSED_EARLY_EVIDENCE=1` accepts static closed-stage evidence for Stage 15.0-15.7 and does not re-run or re-fail them just to reconstruct historical output.

Stage 15.8, Stage 15.9, and Stage 15.10 remain the runtime closure prerequisites. If their `.dat` evidence is missing or invalid and `STAGE15_11_AUTO_RUN_MISSING_PREREQS=1`, Stage 15.11 regenerates those prerequisite checks. Stage 15.9 must receive the separated `STAGE15_9_MAX_STAGE14_RHS_INCREMENT` tolerance so that the Stage 15.7 diagnostic RHS-response bound is not incorrectly applied to the Stage 14.9 production RHS increment.

## Scope

Stage 15.11 is the Stage 15 closure stage. It verifies the Stage 15.0-15.10 controlled structure-state evidence, preserves the Stage 14 and Stage 15 anti-regression protections, and writes the final closure file only after full PASS evidence:

```text
stage15_checks/STAGE15_CLOSED.md
```

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

The only approved coupling chain remains:

```text
Stage 15 controlled structure-state update
-> Stage 12 feedback-force candidate
-> Stage 13 Eulerian force-density candidate
-> Stage 14 controlled RHS diagnostic/injection path
```

Stage 15.11 does not introduce new physics and does not activate full production bending, tension, wall/contact handling, multi-fibre logic, legacy production IBM forcing, pressure/projection, Poisson, RK3, channel forcing, or unrelated DNS numerics.

## Files

Stage 15.11 adds:

```text
stage15_checks/run_stage15_11_total_smoke_closure.sh
stage15_checks/stage15_11_total_smoke_closure.md
```

The wrapper generates this file only on full pass:

```text
stage15_checks/STAGE15_CLOSED.md
```

The closure file is deleted/not generated on partial pass, missing evidence, skipped required evidence, or failed regression checks.

## Manual command

```bash
bash stage15_checks/run_stage15_11_total_smoke_closure.sh
```

The wrapper configures `BUILD_DIR` if needed, builds the required Stage 15 feedback-linkage target, performs static anti-regression audits, verifies or regenerates Stage 15.0-15.10 evidence according to the configured prerequisite flags, writes a Stage 15.11 summary file, and generates `STAGE15_CLOSED.md` only if every required gate passes.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE15_11_RUN_STAGE15_8` | `0` | Force-run the Stage 15.8 prerequisite wrapper. |
| `STAGE15_11_RUN_STAGE15_9` | `0` | Force-run the Stage 15.9 prerequisite wrapper. |
| `STAGE15_11_RUN_STAGE15_10` | `0` | Force-run the Stage 15.10 prerequisite wrapper. |
| `STAGE15_11_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE15_11_REQUIRE_STAGE15_10` | `1` | Require Stage 15.10 evidence. |
| `STAGE15_11_ENABLE` | `1` | Enable Stage 15.11 closure checks. |
| `STAGE15_11_CONTROLLED_STEP_ENABLE` | `1` | Enable controlled-step prerequisite paths. |
| `STAGE15_11_STRUCTURE_ADVANCE_ENABLE` | `1` | Enable only the controlled closure validation. |
| `STAGE15_11_DIAGNOSTIC_ONLY` | `1` | Keep diagnostic/safety mode active. |
| `STAGE15_11_NP` | `2` | MPI process count for np=2 prerequisite paths. |
| `STAGE15_11_NP_LIST` | `1 2 4` | Required parallel consistency list. |
| `STAGE15_11_NPTS` | `8` | Number of controlled structure points. |
| `STAGE15_11_DT` | `1.0e-4` | Controlled step size. |
| `STAGE15_11_RHO_TILDE` | `1.0` | Controlled density ratio. |
| `STAGE15_11_TEST_FORCE` | `1.0e-6` | Small controlled test force. |
| `STAGE15_11_FEEDBACK_ALPHA` | `1.0` | Feedback-force gain. |
| `STAGE15_11_LAMBDA` | `1.0e-8` | Small nonzero Stage 14 diagnostic response scale. |
| `STAGE15_11_MAX_FORCE_RESPONSE` | `1.0e-8` | Maximum allowed feedback-force response. |
| `STAGE15_11_MAX_RHS_RESPONSE` | `1.0e-12` | Maximum allowed linkage RHS response. |
| `STAGE15_11_MAX_STAGE14_RHS_INCREMENT` | `1.0e-4` | Maximum allowed inherited Stage 14 RHS increment. |
| `STAGE15_11_MAX_FLUID_DELTA` | `1.0e-8` | Maximum allowed fluid signature delta. |
| `STAGE15_11_MAX_STRUCTURE_RESTART_DELTA` | `1.0e-14` | Maximum allowed structure restart delta. |
| `STAGE15_11_MAX_IO_SIGNATURE_DELTA` | `1.0e-8` | Maximum allowed I/O signature delta. |
| `STAGE15_11_AUTO_RUN_MISSING_PREREQS` | `1` | Auto-run missing Stage 15.1-15.10 prerequisite wrappers. |

## Static audits

The wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It fails closed if:

- the forbidden Stage 14 lambda-zero hook-registration gate reappears;
- Stage 11.5, Stage 13.6, or Stage 14.5 diagnostic markers are missing;
- old Stage 13.5 force-density diagnostic names reappear;
- rank0-safe diagnostic-writing markers regress;
- Stage 13 force-density sampling repair evidence is missing;
- `stage14_checks/STAGE14_CLOSED.md` is missing when required;
- Stage 15.1-15.11 wrapper/doc markers are missing;
- Stage 15.0 configuration evidence is missing;
- the Stage 15 feedback-linkage check target is missing;
- the wrapper appears to require `rg` or recursive `grep`.

## Runtime and closure evidence

The wrapper verifies or regenerates Stage 15.0-15.10 diagnostics, then verifies:

- controlled structure update status;
- feedback linkage status;
- np=1/2/4 parallel consistency status;
- restart/I/O compatibility status;
- contamination-audit status;
- force/RHS/Stage14 RHS-increment bounded statuses;
- approved Stage 12→13→14 chain status;
- no full structure advance, bending, tension, implicit structure solve, wall/contact, or multi-fibre activation;
- no direct RHS injection outside the approved chain;
- no legacy IBM forcing;
- no pressure/projection/Poisson/RK3/channel-forcing contamination;
- no NaN/Inf in diagnostics.

The Stage 15.11 summary file is:

```text
stage15_outputs/fibre_stage15_11_total_smoke_closure.dat
```

It contains at least:

```text
stage15_11_requested_status
stage15_0_evidence_status
stage15_1_evidence_status
stage15_2_evidence_status
stage15_3_evidence_status
stage15_4_evidence_status
stage15_5_evidence_status
stage15_6_evidence_status
stage15_7_evidence_status
stage15_8_evidence_status
stage15_9_evidence_status
stage15_10_evidence_status
controlled_update_status
feedback_linkage_status
parallel_consistency_status
restart_io_status
contamination_audit_status
stage14_regression_status
stage15_regression_status
approved_stage12_13_14_chain_status
force_response_bounded_status
rhs_response_bounded_status
stage14_rhs_increment_bounded_status
rank0_safe_diagnostic_status
stage13_sampling_repair_status
no_full_structure_advance_status
no_bending_solve_status
no_tension_solve_status
no_wall_contact_status
no_multifibre_status
no_direct_rhs_injection_status
no_legacy_ibm_forcing_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
no_nan_inf_status
stage15_closed_file_status
final_status
```

## Closure file contents

On full PASS only, `stage15_checks/STAGE15_CLOSED.md` records:

- `STAGE15_CLOSED=YES`;
- Stage 15 objective;
- controlled structure update formula;
- feedback-force linkage formula;
- completed substeps Stage 15.0-15.11;
- np=1/2/4 evidence;
- restart/I/O evidence;
- contamination-audit evidence;
- explicit statement that full bending/tension/wall/contact/multi-fibre production physics remains inactive;
- explicit statement that the approved coupling route remains Stage 15 → Stage 12 → Stage 13 → Stage 14;
- next recommended stage: Stage 16 first controlled one-flexible-fibre channel DNS FSI.

## PASS evidence

A successful run prints:

```text
STAGE 15.11 TOTAL SMOKE VERDICT: PASS
STAGE 15.11 FINAL VERDICT: PASS
STAGE15_CLOSED=YES
```

and writes `final_status 1` plus `stage15_closed_file_status 1` in the Stage 15.11 summary file.

## Failure behavior

The wrapper fails if any Stage 15.0-15.10 evidence is missing and cannot be regenerated, if any Stage 14 or Stage 15 regression is detected, if parallel/restart/I/O/contamination evidence fails, if force/RHS/Stage14 increment bounds are exceeded, if solver contamination is detected, if NaN/Inf appears, or if `STAGE15_CLOSED.md` would be generated without full PASS evidence.


## Stage 15.11 closure-prerequisite evidence handling repair

Stage 15.11 accepts Stage 15.0--15.7 as closed early-stage evidence in a fresh unzip tree. When those historical runtime `stage15_outputs/*.dat` files are absent, the closure script must not leave `controlled_update_status` or `feedback_linkage_status` unset; those closure properties are inherited from the stronger Stage 15.8 parallel-consistency evidence and Stage 15.10 contamination-audit evidence.

The wrapper now also emits explicit unmet-status reasons before failure, so `unknown_stage15_11_failure` cannot mask a missing closure gate.


Stage 15.11 review note: `approved_stage12_13_14_chain_status` must be inherited from the stronger Stage 15.8/15.10 closure evidence when those diagnostics pass, so a fresh unzip tree does not fail after accepting Stage 15.0-15.7 closed static evidence. The final wrapper must also emit explicit failure reasons for any unmet closure gate.
