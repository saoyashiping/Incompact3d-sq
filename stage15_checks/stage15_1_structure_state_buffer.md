# Stage 15.1 Production Flexible-Fibre Structure-State Buffer Skeleton

Stage 15.1 introduces a guarded production-owned flexible-fibre structure-state buffer skeleton for future structure advance work.  It allocates, initializes, clears/resets, validates, and diagnoses finite storage only.

Stage 15.1 does **not** solve the eventual flexible-fibre equation

```text
rho_tilde * X_tt = d_s(T d_s X) - gamma * d_ssss X - F_fs
```

and does **not** enforce the future inextensibility constraint

```text
d_s X · d_s X = 1
```

Those equations remain deferred to later stages.

## Added source/check scope

Stage 15.1 owns these new files:

- `src/fibre_stage15_structure_state.f90`
- `src/fibre_stage15_structure_state_check.f90`
- `stage15_checks/run_stage15_1_structure_state_buffer.sh`
- `stage15_checks/stage15_1_structure_state_buffer.md`

The only build-system registration is the standalone check target:

```text
fibre_stage15_structure_state_check
```

The structure-state module is not connected to the production DNS time loop, Stage 12 feedback force, Stage 14 RHS injection, wall/contact handling, or multi-fibre logic.

## Buffer contents

The Stage 15.1 structure-state module provides guarded storage for one future flexible fibre:

- `X_f(npts,3)` position buffer.
- `V_f(npts,3)` velocity buffer.
- `A_f(npts,3)` optional acceleration/RHS placeholder buffer, initialized to zero.
- `T_f(npts)` optional tension-candidate placeholder buffer, initialized to zero.
- Allocation, initialization, clear/reset, validation, finite-value, and diagnostic status fields.
- Zero-valued counters proving that no structure advance, bending solve, tension solve, position time update, or velocity time update has occurred.

`X_f` is initialized as a finite straight reference fibre along the first coordinate.  `V_f`, `A_f`, and `T_f` initialize to zero.  The clear/reset path zeros all buffers and revalidates finite storage.

## Wrapper environment variables

The Stage 15.1 wrapper supports:

```bash
DECOMP2D_ROOT
BUILD_DIR
MPIEXEC
MPIEXEC_FLAGS
STAGE15_1_RUN_STAGE15_0
STAGE15_1_REQUIRE_STAGE14_CLOSED
STAGE15_1_REQUIRE_STAGE15_0
STAGE15_1_NPTS
STAGE15_1_ENABLE
STAGE15_1_STRUCTURE_ADVANCE_ENABLE
STAGE15_1_DIAGNOSTIC_ONLY
```

Safe defaults are:

```bash
BUILD_DIR=build_stage9
MPIEXEC=mpirun
MPIEXEC_FLAGS=
STAGE15_1_RUN_STAGE15_0=0
STAGE15_1_REQUIRE_STAGE14_CLOSED=1
STAGE15_1_REQUIRE_STAGE15_0=1
STAGE15_1_NPTS=8
STAGE15_1_ENABLE=0
STAGE15_1_STRUCTURE_ADVANCE_ENABLE=0
STAGE15_1_DIAGNOSTIC_ONLY=1
```

If `BUILD_DIR` has no `CMakeCache.txt`, the wrapper configures it with CMake and `DECOMP2D_ROOT`, then builds only `fibre_stage15_structure_state_check`.

## Static audit gates

The wrapper fails closed if it detects any Stage 15.1 source connection to:

- Structure advance.
- Bending solve.
- Tension solve.
- Fibre position or velocity time update.
- Stage 12 feedback force.
- Stage 14 RHS injection.
- Wall/contact logic.
- Multi-fibre logic.
- Fluid RHS, pressure/projection, Poisson, RK3, or channel-forcing code paths.

It also checks that the Stage 14 forbidden hook-registration gate

```text
stage14_get_injection_gain() == 0.0
```

is absent; that Stage 11, Stage 13, and Stage 14 production diagnostics remain present; that rank0-safe diagnostic writing remains present; that the Stage 13 force-density sampling repair remains present; and that `stage14_checks/STAGE14_CLOSED.md` exists when required.

## Runtime diagnostics

The standalone check writes:

```text
stage15_outputs/fibre_stage15_1_structure_state_buffer.dat
```

Required fields include:

- `stage15_1_requested_status`
- `allocation_status`
- `initialization_status`
- `clear_status`
- `validation_status`
- `npts`
- `x_finite_status`
- `v_finite_status`
- `optional_a_or_rhs_finite_status`
- `optional_tension_finite_status`
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
stage15_outputs/stage15_1_structure_state_buffer.dat
```

## Manual command

Run Stage 15.1 manually with:

```bash
bash stage15_checks/run_stage15_1_structure_state_buffer.sh
```

Expected successful terminal evidence:

```text
STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: PASS
STAGE 15.1 FINAL VERDICT: PASS
```


## 2026-06-03 robustness note

The Stage 15.1 wrapper must not depend exclusively on `rg`/ripgrep being installed. Static audits use a grep-compatible fallback so that absence of `rg` cannot mask or fabricate Stage 11/13/14 regression evidence. The Stage 13 production force-density diagnostic evidence is keyed to the Stage 13.6 production diagnostic file and status fields (`fibre_stage13_6_production_force_density_candidate.dat`, `stage13_6_production_force_density_candidate_status`).


## 2026-06-03 static-audit correction

The Stage 15.1 static audit must distinguish real production-fluid/RHS access from negative diagnostic status fields.  Identifiers such as `no_fluid_rhs_modification_status` are required diagnostics and must not be treated as actual RHS or fluid-array access.  The audit therefore checks real `call` sites, production velocity/RHS array references, and forbidden production-fluid module imports, while allowing diagnostic names that document the absence of fluid/RHS modification.
