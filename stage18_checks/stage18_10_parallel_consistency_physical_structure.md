# Stage 18.10: parallel consistency for physical structure dynamics

Stage 18.10 is a diagnostic-only, MPI-free partition/reduction emulation gate
for the single-fibre physical-structure response.  It uses Stage 18.9-style
helper-local arrays as the `np=1` reference truth and checks deterministic
contiguous partitions for `np=1`, `np=2`, and `np=4`.

## Partition consistency model

**Global truth** is computed on all Lagrangian fibre points:

```text
X_ref, V_ref, A_ref
F_b_ref, F_T_ref, F_h_ref, F_total_ref
X_next_ref, V_next_ref
E_k_ref, E_b_ref, E_s_ref, P_h_ref
max_response_ref, max_velocity_ref, max_acceleration_ref
```

**Partition emulation** splits the point index range into deterministic
contiguous partitions.  Each partition owns a unique slice; all points are
covered exactly once with no gap and no overlap.  Stage 18.10 uses the preferred
Option A from the stage plan: compute global diagnostic arrays first, then slice
and reduce them.  This validates partition/reduction consistency without
emulating production halo exchange.

## Physical response equations

The dimensional local diagnostic equation is

```text
rho_l A = F_T_candidate + F_b_candidate + F_h_candidate
F_b_candidate = -B X_ssss
```

The nondimensional local diagnostic equation is

```text
rho_tilde A = F_T_candidate + F_b_candidate + F_h_candidate
```

The helper-local candidate step is

```text
A_candidate^n = F_total^n / rho_l
V_candidate^{n+1} = V_candidate^n + dt * A_candidate^n
X_candidate^{n+1} = X_candidate^n + dt * V_candidate^n + 0.5 * dt^2 * A_candidate^n
```

## Energy and power reductions

The helper uses uniform trapezoidal quadrature:

```text
w_0 = ds/2
w_{n-1} = ds/2
w_i = ds for interior points
```

and reduces

```text
E_k = 1/2 int rho_l |V|^2 ds
E_b = 1/2 int B |X_ss|^2 ds
E_s = E_k + E_b
P_h = int F_h dot V ds
```

from local Python partitions back to the `np=1` reference values.

## Diagnostic checks

1. `partition_coverage_case`: confirms coverage once, no gap, no overlap, and
   balanced partition sizes for `np=1`, `np=2`, and `np=4`.
2. `np1_reference_case`: builds the Stage 18.9-style `np=1` truth and checks
   finite arrays, finite nonnegative energies, and bounded response.
3. `partition_reconstruction_case`: reconstructs `X/V/A/F_b/F_T/F_h/F_total`
   and `X_next/V_next` from partitions and compares against `np=1`.
4. `partition_scalar_reduction_case`: deterministically reduces `E_k`, `E_b`,
   `E_s`, `P_h`, and response maxima from partitions.
5. `np2_np4_consistency_case`: requires `np=2` and `np=4` arrays/scalars to
   match the `np=1` reference.
6. `controlled_forced_response_parallel_case`: verifies the controlled `F_h`
   response remains identical under partitioning and is not spread to fluid RHS.
7. `dry_response_parallel_case`: verifies `F_h=0` dry response consistency and
   dry `P_h=0` across partition layouts.
8. `no_production_parallel_contamination_case`: confirms no actual MPI command,
   no production DNS execution, no Stage 14 RHS call, no RHS/IBM/DNS-core change,
   no stats/visu/restart I/O change, no production `X/V/A` update, no production
   structure hook, no contact/collision force, and no production multi-fibre
   logic.

## Diagnostic-only boundary

**PARALLEL CONSISTENCY DIAGNOSTIC** means local helper-only partition/reduction
emulation.  **PRODUCTION PARALLEL RUN** means actual MPI, `mpirun`, `mpiexec`, or
production DNS execution; Stage 18.10 does not do this.  **LOCAL PARTITION**
means temporary Python-side slices of one fibre; they are not written into
production state, RHS, IBM, or runtime output.  **NP=1 REFERENCE** means global
helper truth, not production serial DNS.

Protective references to `src/fibre_stage14_production_rhs_injection.f90`,
`src/fibre_stage13_production_force_density_candidate.f90`, `src/xcompact3d.f90`,
`statistics.f90`, `visu.f90`, restart I/O files, `MPIEXEC`, or `MPIEXEC_FLAGS`
are compatibility/audit text only and are not runtime activation.

## Defaults

```text
STAGE18_10_ENABLE=1
STAGE18_10_PARALLEL_CONSISTENCY_ENABLE=1
STAGE18_10_SINGLE_FIBRE_ONLY=1
STAGE18_10_DIAGNOSTIC_ONLY=1
STAGE18_10_NP_LIST=1,2,4
STAGE18_10_NPTS=64
STAGE18_10_NSTEPS=5
STAGE18_10_FIBRE_LENGTH=1.0
STAGE18_10_COMPONENT_DIM=3
STAGE18_10_RHO_L=1.0
STAGE18_10_RHO_TILDE=1.0
STAGE18_10_BENDING_STIFFNESS=1.0e-3
STAGE18_10_GAMMA=1.0e-3
STAGE18_10_USE_DIMENSIONAL_RESPONSE=1
STAGE18_10_USE_NONDIMENSIONAL_RESPONSE=1
STAGE18_10_DT_STRUCTURE=1.0e-4
STAGE18_10_SINE_EPS=1.0e-3
STAGE18_10_SINE_MODE=1
STAGE18_10_FLUID_FORCE_MAG=1.0e-3
STAGE18_10_INITIAL_VELOCITY_MAG=0.0
STAGE18_10_RESPONSE_BOUND=1.0e-4
STAGE18_10_VELOCITY_BOUND=1.0e-3
STAGE18_10_ACCELERATION_BOUND=1.0e-2
STAGE18_10_ZERO_TOL=1.0e-14
STAGE18_10_FORMULA_TOL=1.0e-12
STAGE18_10_REDUCTION_TOL=1.0e-12
STAGE18_10_ENERGY_TOL=1.0e-12
STAGE18_10_BOUNDED_TOL=1.0e-8
STAGE18_10_TEST_CASE=np1_np2_np4_partition_reduction_consistency
```

## Manual command and expected PASS evidence

Run:

```bash
stage18_checks/run_stage18_10_parallel_consistency_physical_structure.sh
```

Expected terminal evidence:

```text
STAGE 18.10 PARALLEL CONSISTENCY PHYSICAL STRUCTURE VERDICT: PASS
STAGE 18.10 FINAL VERDICT: PASS
```

The wrapper writes
`stage18_outputs/fibre_stage18_10_parallel_consistency_physical_structure.dat`
with all required Stage 18.10 `*_status` fields and `final_status PASS`.
