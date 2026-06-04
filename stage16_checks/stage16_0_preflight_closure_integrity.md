# Stage 16.0 preflight / closure integrity

Stage 16.0 is a diagnostic-only preflight audit before any Stage 16 one-fibre
channel-DNS FSI production physics is introduced. It verifies that the closed
Stage 10-15 evidence remains intact and that the approved Stage 15 -> Stage 12
-> Stage 13 -> Stage 14 coupling route is still the only validated route into
controlled RHS response.

## Manual command

Run exactly:

```bash
bash stage16_checks/run_stage16_0_preflight_closure_integrity.sh
```

The wrapper writes:

```text
stage16_outputs/fibre_stage16_0_preflight_closure_integrity.dat
```

A successful run prints both required PASS lines:

```text
STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS
STAGE 16.0 FINAL VERDICT: PASS
```

A failed run prints `FAIL` and explicit reasons. Unknown-only failure is not an
acceptable Stage 16.0 outcome.

## Environment variables

The wrapper supports these controls with safe defaults:

| Variable | Default | Purpose |
| --- | ---: | --- |
| `BUILD_DIR` | `build_stage9` | CMake build directory used for the preflight build-dir check. |
| `DECOMP2D_ROOT` | empty | Passed through `-DCMAKE_PREFIX_PATH` when configuring a missing build directory. |
| `MPIEXEC` | `mpirun` | Preserved for consistency with closed-stage wrappers. Stage 16.0 itself does not run MPI validation. |
| `MPIEXEC_FLAGS` | empty | Preserved for consistency with closed-stage wrappers. |
| `STAGE16_0_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md` to exist and be nonempty. |
| `STAGE16_0_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md` to exist and be nonempty. |
| `STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE` | `1` | Accept Stage 15.1-15.11 evidence through `STAGE15_CLOSED.md` instead of requiring fresh output files. |
| `STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING` | `0` | If explicitly enabled, run Stage 15.11 to regenerate missing Stage 15 closure evidence. |
| `STAGE16_0_ENABLE` | `1` | Stage 16.0 request gate. |
| `STAGE16_0_DIAGNOSTIC_ONLY` | `1` | Must remain diagnostic-only; Stage 16.0 must not activate new physics. |

## Audit coverage

Stage 16.0 checks these preflight categories:

1. **Closure-file integrity**
   - `stage14_checks/STAGE14_CLOSED.md` exists and is nonempty.
   - `stage15_checks/STAGE15_CLOSED.md` exists and is nonempty.
   - Stage 15 closure evidence explicitly states Stage 15 is closed.
   - Stage 15 closure evidence states production bending, tension, wall/contact,
     and multi-fibre physics remain inactive.
   - Stage 15 closure evidence preserves the Stage 15 -> Stage 12 -> Stage 13
     -> Stage 14 approved coupling route.

2. **Stage 14 anti-regression protections**
   - The forbidden `stage14_get_injection_gain() == 0.0` hook-registration gate
     is absent.
   - Small nonzero lambda Stage 14 RHS diagnostics remain present and are not
     blocked.
   - Stage 11.5, Stage 13.6, and Stage 14.5 production diagnostics still exist.
   - Stage 13 production force-density diagnostics use the current Stage 13.6
     names, not old Stage 13.5 production names.
   - Rank0-safe diagnostic writing is preserved for Stage 11, Stage 13,
     Stage 14, and Stage 15 diagnostics.
   - The Stage 13 force-density sampling repair remains in place and local
     subdomain-centre sampling is rejected.

3. **Stage 15 anti-regression protections**
   - Stage 15.1 through Stage 15.11 scripts and documentation exist.
   - Stage 15.1 through Stage 15.11 evidence is present or is accepted through
     the closed Stage 15 closure file.
   - Stage 15 scripts must not depend on ripgrep without a grep fallback.
   - Stage 15.11 must provide specific failure reasons and must not rely on an
     unknown-only failure.
   - Stage 15.11 closure generation remains guarded by full PASS evidence.

4. **Stage 16.0 boundary protections**
   - No Stage 16 source file is allowed yet.
   - No Stage 16 production hook may appear in `src/xcompact3d.f90`.
   - Stage 16.0 must not introduce direct RHS injection, legacy IBM forcing,
     bending solve, tension solve, wall/contact handling, multi-fibre logic,
     pressure/projection changes, Poisson changes, RK3 changes, or channel
     forcing changes.
   - Stage 16.0 scripts must not modify closed Stage 10-15 files.

5. **Approved coupling-chain audit**
   - Stage 15 controlled structure update evidence remains present.
   - Stage 12 feedback linkage remains the approved feedback route.
   - Stage 13 force-density candidate evidence remains present.
   - Stage 14 controlled RHS path remains present.
   - Direct Stage 14 RHS injection outside the approved chain and legacy IBM
     forcing outside the approved chain remain inactive.

## Expected summary keys

The Stage 16.0 summary file contains these status keys, where `1` means PASS
and `0` means FAIL:

```text
stage16_0_requested_status
stage14_closed_file_status
stage15_closed_file_status
stage15_closed_content_status
stage14_regression_status
stage15_regression_status
stage13_6_diagnostic_name_status
stage13_sampling_repair_status
rank0_safe_diagnostic_status
no_rg_only_dependency_status
no_unknown_failure_status
stage16_boundary_status
approved_stage12_13_14_chain_status
no_direct_rhs_injection_status
no_legacy_ibm_forcing_status
no_bending_solve_status
no_tension_solve_status
no_wall_contact_status
no_multifibre_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
final_status
```

## Stage 16.0 boundary

This stage intentionally adds only Stage 16.0 audit scripts and documentation.
It does not add one-fibre FSI physics, production structure advance, bending or
tension solves, wall/contact handling, multi-fibre logic, fluid RHS changes,
pressure/projection/Poisson changes, RK3 changes, or channel-forcing changes.

## Stage 16.0 repair notes

This Stage 16.0 audit must treat Stage 15 as a closed stage when
`stage15_checks/STAGE15_CLOSED.md` is present and complete. In a fresh unpacked
source tree, historical `stage15_outputs/*.dat` runtime files may be absent; the
closure file is the versioned evidence that Stage 15.1--15.11 were closed.

The preflight audit intentionally scans production/source evidence for the
forbidden Stage 14 lambda-zero registration gate and obsolete Stage 13.5
production diagnostic names. Documentation may mention those old strings as
negative anti-regression examples, and those negative mentions must not be
misclassified as active regressions.

The `rg`/ripgrep audit checks actual command invocations only. Regex strings or
comments that merely name `rg` as a forbidden dependency are not treated as a
runtime dependency.
