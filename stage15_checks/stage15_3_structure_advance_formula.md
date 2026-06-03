# Stage 15.3 Standalone Controlled Structure-Advance Formula Checks

Stage 15.3 adds a standalone, diagnostic-only controlled formula check for future flexible-fibre structure updates.  It does not connect any structure update to the production `xcompact3d` time loop and does not solve the full constrained structure equation.

The eventual equation remains:

```text
rho_tilde * X_tt = d_s(T d_s X) - gamma * d_ssss X - F_fs
```

with inextensibility:

```text
d_s X · d_s X = 1
```

Stage 15.3 does not solve that system.  It checks only a controlled explicit candidate formula in standalone diagnostics:

```text
A_f_cand = F_test / rho_tilde
V_f_new  = V_f_old + dt * A_f_cand
X_f_new  = X_f_old + dt * V_f_new
```

## Added source/check scope

Stage 15.3 owns these new files:

- `src/fibre_stage15_structure_advance_formula.f90`
- `src/fibre_stage15_structure_advance_formula_check.f90`
- `stage15_checks/run_stage15_3_structure_advance_formula.sh`
- `stage15_checks/stage15_3_structure_advance_formula.md`

The only build-system registration is the standalone check target:

```text
fibre_stage15_structure_advance_formula_check
```

The formula module is not added to `xcompact3d`, Stage 14 RHS injection, fluid RHS, pressure/projection, Poisson, RK3, channel forcing, wall/contact, or multi-fibre code paths.

## Controlled checks

The standalone check verifies:

- Zero-force response has no drift beyond `STAGE15_3_MAX_ZERO_FORCE_DRIFT`.
- Small-force acceleration, velocity update, and position update are finite and bounded by `STAGE15_3_MAX_SMALL_FORCE_UPDATE`.
- Acceleration/update signs are consistent with the applied test force.
- Force scaling is linear for a doubled small test force.
- `dt` scaling is consistent: velocity update scales by two and position update scales by four when `dt` is doubled for the controlled formula.
- Production connection counters remain zero.
- Solver/RHS modification status fields remain safe.

## Wrapper environment variables

The Stage 15.3 wrapper supports:

```bash
DECOMP2D_ROOT
BUILD_DIR
MPIEXEC
MPIEXEC_FLAGS
STAGE15_3_RUN_STAGE15_2
STAGE15_3_REQUIRE_STAGE14_CLOSED
STAGE15_3_REQUIRE_STAGE15_2
STAGE15_3_NPTS
STAGE15_3_DT
STAGE15_3_RHO_TILDE
STAGE15_3_TEST_FORCE
STAGE15_3_MAX_ZERO_FORCE_DRIFT
STAGE15_3_MAX_SMALL_FORCE_UPDATE
STAGE15_3_MAX_SIGN_ERROR
STAGE15_3_MAX_SCALING_ERROR
STAGE15_3_ENABLE
STAGE15_3_STRUCTURE_ADVANCE_ENABLE
STAGE15_3_DIAGNOSTIC_ONLY
```

Safe defaults are:

```bash
BUILD_DIR=build_stage9
MPIEXEC=mpirun
MPIEXEC_FLAGS=
STAGE15_3_RUN_STAGE15_2=0
STAGE15_3_REQUIRE_STAGE14_CLOSED=1
STAGE15_3_REQUIRE_STAGE15_2=1
STAGE15_3_NPTS=8
STAGE15_3_DT=1.0e-4
STAGE15_3_RHO_TILDE=1.0
STAGE15_3_TEST_FORCE=1.0e-6
STAGE15_3_MAX_ZERO_FORCE_DRIFT=1.0e-14
STAGE15_3_MAX_SMALL_FORCE_UPDATE=1.0e-6
STAGE15_3_MAX_SIGN_ERROR=1.0e-14
STAGE15_3_MAX_SCALING_ERROR=1.0e-12
STAGE15_3_ENABLE=1
STAGE15_3_STRUCTURE_ADVANCE_ENABLE=0
STAGE15_3_DIAGNOSTIC_ONLY=1
```

If `BUILD_DIR` has no `CMakeCache.txt`, the wrapper configures it with CMake and `DECOMP2D_ROOT`, then builds only `fibre_stage15_structure_advance_formula_check`.

The wrapper intentionally uses POSIX `grep`/`awk` static checks and does not require `ripgrep`.

## Static audit gates

The wrapper fails closed if it detects Stage 15.3 source calls or links to:

- Production structure advance.
- Bending solve.
- Tension solve.
- Wall/contact logic.
- Multi-fibre logic.
- `xcompact3d` production time loop.
- Stage 14 RHS injection.
- Fluid RHS, pressure/projection, Poisson, RK3, or channel forcing.

It also checks that:

- The forbidden Stage 14 hook-registration gate `stage14_get_injection_gain() == 0.0` remains absent.
- Stage 11.5, Stage 13, Stage 14, Stage 15.1, and Stage 15.2 diagnostics remain present.
- Rank0-safe Stage 11/13/14 diagnostic writing guards remain present.
- The Stage 13 force-density sampling repair remains present.
- `stage14_checks/STAGE14_CLOSED.md` exists when required.
- Stage 15.2 evidence exists, or the Stage 15.2 wrapper is explicitly requested with `STAGE15_3_RUN_STAGE15_2=1`.

## Runtime diagnostics

The standalone check writes:

```text
stage15_outputs/fibre_stage15_3_structure_advance_formula.dat
```

Required fields include:

- `stage15_3_requested_status`
- `npts`
- `dt`
- `rho_tilde`
- `test_force_magnitude`
- `zero_force_drift`
- `zero_force_status`
- `small_force_max_acceleration`
- `small_force_max_velocity_update`
- `small_force_max_position_update`
- `small_force_status`
- `sign_consistency_status`
- `force_scaling_error`
- `force_scaling_status`
- `dt_scaling_error`
- `dt_scaling_status`
- `finite_value_status`
- `production_time_loop_connection_count`
- `production_structure_advance_count`
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
stage15_outputs/stage15_3_structure_advance_formula.dat
```

## Manual command

Run Stage 15.3 manually with:

```bash
bash stage15_checks/run_stage15_3_structure_advance_formula.sh
```

Expected successful terminal evidence:

```text
STAGE 15.3 STRUCTURE ADVANCE FORMULA VERDICT: PASS
STAGE 15.3 FINAL VERDICT: PASS
```


## Stage 15.3 static-audit compatibility note

Stage 13 production force-density diagnostics use the closed Stage 13.6 production diagnostic names:

- `stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat`
- `stage13_6_production_force_density_candidate_status`

Stage 15.3 must audit those current Stage 13.6 names. It must not regress to obsolete Stage 13.5 diagnostic names, because that causes a false `stage13_production_force_density_diagnostics_missing` failure even when the Stage 13 production diagnostic path is present and preserved.
