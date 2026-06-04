# Stage 16.0: preflight / closure integrity

## Scope

Stage 16.0 is a preflight audit for beginning Stage 16. It verifies that prerequisite closure evidence from Stage 14 and Stage 15 is intact before any one-flexible-fibre channel-DNS FSI physics is introduced.

Stage 16.0 does **not** add production FSI physics. It does not activate structure advance beyond the already closed Stage 15 controlled evidence, full bending, full tension, wall/contact handling, multi-fibre logic, fluid RHS changes, pressure/projection, Poisson, RK3, channel forcing, or legacy IBM forcing outside the approved Stage 11-14 chain.

The approved Stage 15 closure chain remains:

```text
Stage 15 controlled structure-state update
-> Stage 12 feedback-force candidate
-> Stage 13 Eulerian force-density candidate
-> Stage 14 controlled RHS diagnostic/injection path
```

## Files

Stage 16.0 adds:

```text
stage16_checks/run_stage16_0_preflight_closure_integrity.sh
stage16_checks/stage16_0_preflight_closure_integrity.md
```

No Stage 16 source file or production hook is added in this preflight stage.

## Manual command

```bash
bash stage16_checks/run_stage16_0_preflight_closure_integrity.sh
```

The wrapper configures `BUILD_DIR` if needed, audits Stage 14 and Stage 15 closure evidence, verifies Stage 16.0 has not introduced production physics, and writes a Stage 16.0 summary file.

## Environment variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `DECOMP2D_ROOT` | `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4` | Prefix path used when configuring a missing build directory. |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `MPIEXEC` | `mpirun` | MPI launcher forwarded only if Stage 15.11 auto-run is explicitly enabled. |
| `MPIEXEC_FLAGS` | empty | Optional MPI launcher flags. |
| `STAGE16_0_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE16_0_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md`. |
| `STAGE16_0_ACCEPT_CLOSED_STAGE15_EVIDENCE` | `1` | Accept valid `STAGE15_CLOSED.md` as closed Stage 15 evidence. |
| `STAGE16_0_AUTO_RUN_STAGE15_11_IF_MISSING` | `0` | Optionally run Stage 15.11 if `STAGE15_CLOSED.md` is missing. Disabled by default to preserve closed files. |
| `STAGE16_0_ENABLE` | `1` | Enable Stage 16.0 preflight checks. |
| `STAGE16_0_DIAGNOSTIC_ONLY` | `1` | Keep Stage 16.0 diagnostic/preflight-only. |

## Required closure checks

The wrapper fails closed unless:

- `stage14_checks/STAGE14_CLOSED.md` exists and is nonempty when required;
- `stage15_checks/STAGE15_CLOSED.md` exists and is nonempty when required;
- `STAGE15_CLOSED.md` includes `STAGE15_CLOSED=YES`;
- `STAGE15_CLOSED.md` states that Stage 15 is closed;
- `STAGE15_CLOSED.md` preserves the inactive full bending/tension/wall/contact/multi-fibre production-physics statement;
- `STAGE15_CLOSED.md` preserves the approved `Stage 15 -> Stage 12 -> Stage 13 -> Stage 14` coupling route;
- `STAGE15_CLOSED.md` preserves the Stage 16 recommendation.

## Static anti-regression audits

The wrapper uses portable `grep`/`awk` checks and does not require ripgrep. It verifies:

- the forbidden `stage14_get_injection_gain() == 0.0` registration gate is absent;
- Stage 11.5, Stage 13.6, and Stage 14.5 diagnostic markers remain present;
- old Stage 13.5 production force-density diagnostic names do not reappear;
- rank0-safe diagnostic-writing markers remain present;
- Stage 13 force-density sampling repair evidence remains present;
- Stage 15.1-15.11 wrapper/doc markers remain present;
- Stage 15.11 has a reasons-file failure path rather than an unknown-only failure mode;
- no Stage 16 source files are present yet;
- no Stage 16 hook/call is inserted into `xcompact3d.f90`;
- Stage 16 scripts do not call direct RHS injection, IBM forcing, pressure/projection/Poisson/RK3/channel forcing, bending/tension, wall/contact, multi-fibre, or structure-advance routines;
- Stage 16 scripts do not write into closed Stage 10-15 directories or `src/`.

## Diagnostic output

The wrapper writes:

```text
stage16_outputs/fibre_stage16_0_preflight_closure_integrity.dat
```

The summary contains at least:

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

## PASS evidence

A successful run prints:

```text
STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS
STAGE 16.0 FINAL VERDICT: PASS
```

and writes `final_status 1` in the Stage 16.0 summary file.

## Failure behavior

The wrapper prints explicit reasons on failure. It fails if closure files are missing or incomplete, Stage 14 or Stage 15 regression markers are detected, Stage 13.6 diagnostics regress, rank0-safe diagnostics regress, Stage 13 local-subdomain-center sampling evidence reappears or the repair is missing, any Stage 16 production physics appears, direct RHS injection or legacy IBM forcing is detected, pressure/projection/Poisson/RK3/channel-forcing contamination is detected, or required docs/scripts/build evidence is missing.
