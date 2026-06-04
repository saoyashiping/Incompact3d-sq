# Stage 15.2 Prescribed Fibre Velocity to Structure-Owned Velocity Adapter

Stage 15.2 replaces the standalone prescribed/placeholder fibre-velocity source used in controlled Stage 12 feedback-force diagnostics with a Stage 15 structure-owned velocity source, while preserving identical numerical values.

The Stage 12 feedback candidate remains:

```text
F_fs_cand = alpha * (U_f - V_f)
```

Stage 15.2 changes only the source interface for `V_f` in the standalone adapter check:

```text
V_f_stage15_structure_owned == V_f_prescribed
```

within strict diagnostic tolerance.  It does not advance the fibre, does not update `X_f` or `V_f` in time, does not solve bending or tension, and does not connect directly to Stage 14 RHS injection except through the already validated Stage 12 feedback-formula diagnostic chain.

## Added source/check scope

Stage 15.2 owns these new files:

- `src/fibre_stage15_velocity_source_adapter.f90`
- `src/fibre_stage15_velocity_source_adapter_check.f90`
- `stage15_checks/run_stage15_2_velocity_source_adapter.sh`
- `stage15_checks/stage15_2_velocity_source_adapter.md`

It also uses a minimal read/write accessor extension in `src/fibre_stage15_structure_state.f90` so the new adapter can initialize and read the Stage 15.1 structure-owned `V_f` buffer without changing Stage 15.1 default behavior.

The only build-system registration is the standalone check target:

```text
fibre_stage15_velocity_source_adapter_check
```

## Adapter behavior

The adapter:

1. Receives a prescribed velocity reference from the current Stage 12 convention.
2. Copies that finite reference into the Stage 15.1 structure-owned `V_f` buffer.
3. Reads `V_f` back from the Stage 15 structure state.
4. Compares `V_f_stage15` against `V_f_prescribed` with `STAGE15_2_MAX_VELOCITY_DIFF`.
5. Uses the existing Stage 12 feedback-formula routine to compare `F_fs_cand(U_f,V_f_stage15)` with `F_fs_cand(U_f,V_f_prescribed)` using `STAGE15_2_MAX_FORCE_DIFF`.
6. Confirms zero-slip behavior remains zero.

No structure advance, bending solve, tension solve, fibre position time update, or fibre velocity time update routine is introduced.

## Wrapper environment variables

The Stage 15.2 wrapper supports:

```bash
DECOMP2D_ROOT
BUILD_DIR
MPIEXEC
MPIEXEC_FLAGS
STAGE15_2_RUN_STAGE15_1
STAGE15_2_REQUIRE_STAGE14_CLOSED
STAGE15_2_REQUIRE_STAGE15_1
STAGE15_2_NPTS
STAGE15_2_ENABLE
STAGE15_2_STRUCTURE_ADVANCE_ENABLE
STAGE15_2_DIAGNOSTIC_ONLY
STAGE15_2_MAX_VELOCITY_DIFF
STAGE15_2_MAX_FORCE_DIFF
```

Safe defaults are:

```bash
BUILD_DIR=build_stage9
MPIEXEC=mpirun
MPIEXEC_FLAGS=
STAGE15_2_RUN_STAGE15_1=0
STAGE15_2_REQUIRE_STAGE14_CLOSED=1
STAGE15_2_REQUIRE_STAGE15_1=1
STAGE15_2_NPTS=8
STAGE15_2_ENABLE=1
STAGE15_2_STRUCTURE_ADVANCE_ENABLE=0
STAGE15_2_DIAGNOSTIC_ONLY=1
STAGE15_2_MAX_VELOCITY_DIFF=1.0e-14
STAGE15_2_MAX_FORCE_DIFF=1.0e-14
```

If `BUILD_DIR` has no `CMakeCache.txt`, the wrapper configures it with CMake and `DECOMP2D_ROOT`, then builds only `fibre_stage15_velocity_source_adapter_check`.

## Static audit gates

The wrapper fails closed if it detects Stage 15.2 source calls to:

- Structure advance.
- Bending solve.
- Tension solve.
- Fibre position or velocity time update.
- Wall/contact logic.
- Multi-fibre logic.
- Pressure/projection, Poisson, RK3, or channel-forcing code paths.

It also checks that:

- The forbidden Stage 14 hook-registration gate `stage14_get_injection_gain() == 0.0` remains absent.
- Stage 11.5, Stage 13, Stage 14, and Stage 15.1 production diagnostics remain present.
- Rank0-safe Stage 11/13/14 diagnostic writing guards remain present.
- The Stage 13 force-density sampling repair remains present.
- `stage14_checks/STAGE14_CLOSED.md` exists when required.
- Stage 15.1 evidence exists, or the Stage 15.1 wrapper is explicitly requested with `STAGE15_2_RUN_STAGE15_1=1`.

## Runtime diagnostics

The standalone check writes:

```text
stage15_outputs/fibre_stage15_2_velocity_source_adapter.dat
```

Required fields include:

- `stage15_2_requested_status`
- `stage15_structure_owned_velocity_status`
- `prescribed_velocity_reference_status`
- `velocity_source_adapter_status`
- `npts`
- `max_velocity_source_diff`
- `max_feedback_force_diff`
- `velocity_equivalence_status`
- `feedback_force_equivalence_status`
- `zero_slip_status`
- `finite_value_status`
- `structure_advance_count`
- `bending_solve_count`
- `tension_solve_count`
- `position_time_update_count`
- `velocity_time_update_count`
- `no_fluid_rhs_modification_status`
- `no_pressure_projection_modification_status`
- `no_poisson_modification_status`
- `no_rk3_channel_forcing_modification_status`
- `final_status`

The wrapper also writes a gate summary:

```text
stage15_outputs/stage15_2_velocity_source_adapter.dat
```

## Manual command

Run Stage 15.2 manually with:

```bash
bash stage15_checks/run_stage15_2_velocity_source_adapter.sh
```

Expected successful terminal evidence:

```text
STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: PASS
STAGE 15.2 FINAL VERDICT: PASS
```


Stage 15.11 review note: the wrapper must not require ripgrep; static searches use grep-compatible helpers so the check runs on systems without `rg`.
