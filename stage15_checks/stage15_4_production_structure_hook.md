# Stage 15.4 Production Structure-Advance Diagnostic Hook Skeleton

Stage 15.4 inserts a guarded production Stage 15 structure hook into `xcompact3d` in diagnostic-only/no-op mode.  The hook records request, registration, apply, finalize, time-step context, structure-state availability, and no-op diagnostics.  It does not advance the fibre, update `X_f`, update `V_f`, update `A_f`, solve bending/tension, modify Stage 14 RHS injection, or alter fluid/solver behavior.

The eventual structure equation remains:

```text
rho_tilde * X_tt = d_s(T d_s X) - gamma * d_ssss X - F_fs
```

with inextensibility:

```text
d_s X · d_s X = 1
```

Stage 15.4 does not solve that equation.  It only adds the production diagnostic hook skeleton required for future stages.

## Added source/check scope

Stage 15.4 owns these new files:

- `src/fibre_stage15_production_structure_hook.f90`
- `src/fibre_stage15_production_structure_hook_check.f90`
- `stage15_checks/run_stage15_4_production_structure_hook.sh`
- `stage15_checks/stage15_4_production_structure_hook.md`

Stage 15.4 also makes minimal production hook insertion changes to:

- `src/xcompact3d.f90`
- `src/CMakeLists.txt`

The build-system registration adds only the standalone check target:

```text
fibre_stage15_production_structure_hook_check
```

The main `xcompact3d` target includes the new diagnostic hook module only so the guarded no-op production hook can be called when Stage 15 is explicitly requested.

## Production hook behavior

The hook may record:

- `stage15_4_requested_status`
- hook registration status
- hook apply count
- hook finalize status
- diagnostic-only status
- no-op status
- structure-state availability
- production time-loop hook status
- last time-step/RK-substep context

The hook records zero counts for:

- production structure advance
- `X_f` position update
- `V_f` velocity update
- `A_f` production acceleration update
- bending solve
- tension solve
- wall/contact logic
- multi-fibre logic
- Stage 14 RHS injection connection

## Wrapper environment variables

The Stage 15.4 wrapper supports:

```bash
DECOMP2D_ROOT
BUILD_DIR
MPIEXEC
MPIEXEC_FLAGS
STAGE15_4_RUN_STAGE15_3
STAGE15_4_REQUIRE_STAGE14_CLOSED
STAGE15_4_REQUIRE_STAGE15_3
STAGE15_4_ENABLE
STAGE15_4_STRUCTURE_ADVANCE_ENABLE
STAGE15_4_DIAGNOSTIC_ONLY
STAGE15_4_NP
STAGE15_4_MAX_FLUID_DELTA
STAGE15_4_RUN_PRODUCTION_SMOKE
```

Safe defaults are:

```bash
BUILD_DIR=build_stage9
MPIEXEC=mpirun
MPIEXEC_FLAGS=
STAGE15_4_RUN_STAGE15_3=0
STAGE15_4_REQUIRE_STAGE14_CLOSED=1
STAGE15_4_REQUIRE_STAGE15_3=1
STAGE15_4_ENABLE=1
STAGE15_4_STRUCTURE_ADVANCE_ENABLE=0
STAGE15_4_DIAGNOSTIC_ONLY=1
STAGE15_4_NP=1
STAGE15_4_MAX_FLUID_DELTA=0.0
STAGE15_4_RUN_PRODUCTION_SMOKE=1
```

If `BUILD_DIR` has no `CMakeCache.txt`, the wrapper configures it with CMake and `DECOMP2D_ROOT`, then builds only `fibre_stage15_production_structure_hook_check`.

The wrapper intentionally uses POSIX `grep`/`awk` static checks and does not require `ripgrep`.

## Production smoke status

The wrapper always runs the standalone hook check.  A lightweight production smoke is recorded as deferred in this stage because the existing safe Stage 15.4 evidence is the standalone hook check plus static verification of the guarded `xcompact3d` insertion.  The wrapper writes `stage15_4_production_smoke_deferred_status 1` in its gate summary rather than silently skipping production-smoke evidence.

## Static audit gates

The wrapper fails closed if it detects Stage 15.4 source calls or links to:

- production structure advance
- bending solve
- tension solve
- wall/contact logic
- multi-fibre logic
- Stage 14 RHS injection
- fluid RHS modification
- pressure/projection, Poisson, RK3, or channel forcing

It also checks that:

- `xcompact3d.f90` contains only the approved guarded Stage 15.4 diagnostic/no-op hook insertion.
- `xcompact3d.f90` does not call Stage 15 production structure advance.
- The forbidden Stage 14 hook-registration gate `stage14_get_injection_gain() == 0.0` remains absent.
- Stage 11.5, Stage 13.6, Stage 14, Stage 15.1, Stage 15.2, and Stage 15.3 diagnostics remain present.
- Rank0-safe Stage 11/13/14/15 diagnostic writing guards remain present.
- The Stage 13 force-density sampling repair remains present.
- `stage14_checks/STAGE14_CLOSED.md` exists when required.

## Runtime diagnostics

The standalone check writes:

```text
stage15_outputs/fibre_stage15_4_production_structure_hook.dat
```

Required fields include:

- `stage15_4_requested_status`
- `hook_registered_status`
- `hook_apply_count`
- `hook_finalize_status`
- `diagnostic_only_status`
- `noop_status`
- `structure_state_available_status`
- `production_time_loop_hook_status`
- `production_time_loop_connection_count`
- `production_structure_advance_count`
- `x_position_update_count`
- `v_velocity_update_count`
- `a_acceleration_update_count`
- `bending_solve_count`
- `tension_solve_count`
- `wall_contact_count`
- `multifibre_count`
- `rhs_injection_connection_count`
- `no_fluid_rhs_modification_status`
- `no_pressure_projection_modification_status`
- `no_poisson_modification_status`
- `no_rk3_channel_forcing_modification_status`
- `final_status`

The wrapper also writes a gate summary:

```text
stage15_outputs/stage15_4_production_structure_hook.dat
```

## Manual command

Run Stage 15.4 manually with:

```bash
bash stage15_checks/run_stage15_4_production_structure_hook.sh
```

Expected successful terminal evidence:

```text
STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: PASS
STAGE 15.4 FINAL VERDICT: PASS
```
