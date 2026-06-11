# Stage 18.11: restart / I/O compatibility for physical structure state

Stage 18.11 is a diagnostic-only restart/I/O compatibility gate for the
single-fibre physical structure state.  It writes helper-local JSON snapshots
only under `stage18_outputs` and never modifies production restart, statistics,
visualisation, RHS, IBM, DNS-core, or runtime structure-state paths.

## Snapshot model

The state vector is

```text
S = {X, V, A, F_b, F_T, F_h, F_total, X_next, V_next}
```

The diagnostic vector is

```text
D = {E_k, E_b, E_s, P_h, max_response, max_velocity, max_acceleration}
```

A helper-local snapshot is

```text
R = serialize(S, D, metadata)
S_loaded, D_loaded, metadata_loaded = deserialize(R)
```

Roundtrip consistency requires

```text
||S_loaded - S|| <= tol
|D_loaded - D| <= tol
metadata_loaded == metadata
```

## Required schema fields

The snapshot contains `schema_name`, `schema_version`, `stage_id`, `npts`,
`component_dim`, `fibre_length`, `ds`, `dt_structure`, `rho_l`, `rho_tilde`,
`bending_stiffness`, `gamma`, `np_list`, all state arrays, all scalar
diagnostics, exact array shape metadata, and a deterministic digest.

The required values are:

```text
schema_name = stage18_11_physical_structure_state_snapshot
stage_id = 18.11
component_dim = 3
```

## Restart recompute equivalence

Given the loaded state, Stage 18.11 recomputes

```text
A_loaded = F_total_loaded / rho_l
X_next_loaded = X_loaded + dt*V_loaded + 0.5*dt^2*A_loaded
V_next_loaded = V_loaded + dt*A_loaded
E_loaded = energy(S_loaded)
P_loaded = power(F_h_loaded, V_loaded)
```

and compares the recomputed values with the serialized values.

## Partition snapshot compatibility

Stage 18.11 reuses Stage 18.10 deterministic partition/reduction emulation for
`np=1`, `np=2`, and `np=4`.  Partition snapshot chunks are helper-local JSON
objects, not MPI I/O and not production parallel restart.  Reloading partition
chunks must reconstruct the global arrays and reduce scalar diagnostics back to
the `np=1` helper truth.

## Diagnostic-only boundary

**RESTART/I/O COMPATIBILITY DIAGNOSTIC** means helper-local JSON snapshot
roundtrip under `stage18_outputs`.  **PRODUCTION RESTART I/O** means runtime
restart files or Fortran restart module changes; Stage 18.11 does not do this.
**LOCAL SNAPSHOT** means temporary helper output for verification only.
**PARTITION SNAPSHOT** means helper-local emulated chunks, not MPI I/O and not
production parallel restart.

Protective references to `src/fibre_stage14_production_rhs_injection.f90`,
`src/fibre_stage13_production_force_density_candidate.f90`, `src/xcompact3d.f90`,
`statistics.f90`, `visu.f90`, restart I/O files, `MPIEXEC`, or `MPIEXEC_FLAGS`
are compatibility/audit text only and are not runtime activation.

## Defaults

```text
STAGE18_11_ENABLE=1
STAGE18_11_RESTART_IO_COMPATIBILITY_ENABLE=1
STAGE18_11_SINGLE_FIBRE_ONLY=1
STAGE18_11_DIAGNOSTIC_ONLY=1
STAGE18_11_NP_LIST=1,2,4
STAGE18_11_NPTS=64
STAGE18_11_NSTEPS=5
STAGE18_11_FIBRE_LENGTH=1.0
STAGE18_11_COMPONENT_DIM=3
STAGE18_11_RHO_L=1.0
STAGE18_11_RHO_TILDE=1.0
STAGE18_11_BENDING_STIFFNESS=1.0e-3
STAGE18_11_GAMMA=1.0e-3
STAGE18_11_DT_STRUCTURE=1.0e-4
STAGE18_11_SINE_EPS=1.0e-3
STAGE18_11_SINE_MODE=1
STAGE18_11_FLUID_FORCE_MAG=1.0e-3
STAGE18_11_INITIAL_VELOCITY_MAG=0.0
STAGE18_11_RESPONSE_BOUND=1.0e-4
STAGE18_11_VELOCITY_BOUND=1.0e-3
STAGE18_11_ACCELERATION_BOUND=1.0e-2
STAGE18_11_ZERO_TOL=1.0e-14
STAGE18_11_FORMULA_TOL=1.0e-12
STAGE18_11_RESTART_TOL=1.0e-12
STAGE18_11_REDUCTION_TOL=1.0e-12
STAGE18_11_ENERGY_TOL=1.0e-12
STAGE18_11_TEST_CASE=snapshot_roundtrip_restart_equivalence_partition_io
```

## Manual command and expected PASS evidence

Run:

```bash
stage18_checks/run_stage18_11_restart_io_physical_structure_state.sh
```

Expected terminal evidence:

```text
STAGE 18.11 RESTART IO PHYSICAL STRUCTURE STATE VERDICT: PASS
STAGE 18.11 FINAL VERDICT: PASS
```

The wrapper writes only helper-local diagnostics under `stage18_outputs`:

```text
stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state.dat
stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_snapshot.json
stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_partition_snapshot.json
```
