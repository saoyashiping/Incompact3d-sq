# Stage 18.9: controlled one-fibre physical response np=1

Stage 18.9 is a diagnostic-only, helper-local controlled response benchmark for
one physical fibre in serial (`np=1`).  It combines the already closed Stage
18.3--18.8 diagnostic concepts without activating production physics.

## Local diagnostic model

The dimensional helper equation is

```text
rho_l A = F_T_candidate + F_b_candidate + F_h_candidate
F_b_candidate = -B X_ssss
```

The nondimensional helper equation is

```text
rho_tilde A = F_T_candidate + F_b_candidate + F_h_candidate
```

The Stage 18.9 candidate step is local to the Python helper:

```text
A_candidate^n = F_total^n / rho_l
V_candidate^{n+1} = V_candidate^n + dt * A_candidate^n
X_candidate^{n+1} = X_candidate^n + dt * V_candidate^n + 0.5 * dt^2 * A_candidate^n
```

The same formula may be checked with `rho_tilde` for nondimensional validation.

## Energy and power bookkeeping

The helper evaluates local diagnostic quantities only:

```text
E_k = 1/2 int rho_l |V|^2 ds
E_b = 1/2 int B |X_ss|^2 ds
E_s = E_k + E_b
P_h = int F_h dot V ds
```

Uniform trapezoidal weights are used:

```text
w_0 = ds/2
w_{n-1} = ds/2
w_i = ds for interior points
```

## Controlled checks

1. `np1_serial_mode_case` verifies `STAGE18_9_NP=1`; no MPI command is needed.
2. `dry_response_control_case` uses `F_h=0`, local candidate arrays, finite
   energies, nonnegative energies, bounded response, and dry `P_h=0`.
3. `forced_response_direction_case` applies a controlled fluid-on-fibre force in
   the positive y direction and checks positive acceleration, velocity, and
   displacement signs without spreading force to the Eulerian RHS.
4. `bending_restoring_response_case` uses
   `X(s)=(s, eps*sin(k*s), 0)` and confirms the bending force and acceleration
   oppose the sine displacement.
5. `controlled_energy_power_case` checks finite, nonnegative energy values and
   local power sign bookkeeping.
6. `bounded_response_window_case` verifies configured displacement, velocity,
   and acceleration bounds across a small local diagnostic window.
7. `no_production_contamination_case` confirms no Stage 14 RHS call, no force
   spreading, no production DNS/MPI execution, no production X/V/A update, and
   no IBM, DNS-core, statistics, visualisation, or restart I/O modification.

## Diagnostic-only boundary

**CONTROLLED RESPONSE DIAGNOSTIC** means helper-local arrays for one fibre and
`np=1`.  **PRODUCTION STRUCTURE RESPONSE** means persistent Fortran/runtime
`X/V/A` update; Stage 18.9 does not do this.  **LOCAL CANDIDATE FORCE/STEP**
means temporary arrays inside the helper; they are not written into production
state, RHS, IBM, or runtime output.  **NP=1 DIAGNOSTIC** means serial helper
execution, not MPI production execution.

Protective references to `src/fibre_stage14_production_rhs_injection.f90`,
`src/fibre_stage13_production_force_density_candidate.f90`, `src/xcompact3d.f90`,
`statistics.f90`, `visu.f90`, restart I/O files, `MPIEXEC`, or `MPIEXEC_FLAGS`
are compatibility/audit text only and are not runtime activation.

## Defaults

The wrapper supplies safe defaults:

```text
STAGE18_9_ENABLE=1
STAGE18_9_CONTROLLED_RESPONSE_ENABLE=1
STAGE18_9_SINGLE_FIBRE_ONLY=1
STAGE18_9_DIAGNOSTIC_ONLY=1
STAGE18_9_NP=1
STAGE18_9_NPTS=64
STAGE18_9_NSTEPS=5
STAGE18_9_FIBRE_LENGTH=1.0
STAGE18_9_COMPONENT_DIM=3
STAGE18_9_RHO_L=1.0
STAGE18_9_RHO_TILDE=1.0
STAGE18_9_BENDING_STIFFNESS=1.0e-3
STAGE18_9_GAMMA=1.0e-3
STAGE18_9_USE_DIMENSIONAL_RESPONSE=1
STAGE18_9_USE_NONDIMENSIONAL_RESPONSE=1
STAGE18_9_DT_STRUCTURE=1.0e-4
STAGE18_9_SINE_EPS=1.0e-3
STAGE18_9_SINE_MODE=1
STAGE18_9_FLUID_FORCE_MAG=1.0e-3
STAGE18_9_INITIAL_VELOCITY_MAG=0.0
STAGE18_9_RESPONSE_BOUND=1.0e-4
STAGE18_9_VELOCITY_BOUND=1.0e-3
STAGE18_9_ACCELERATION_BOUND=1.0e-2
STAGE18_9_ZERO_TOL=1.0e-14
STAGE18_9_FORMULA_TOL=1.0e-12
STAGE18_9_ENERGY_TOL=1.0e-12
STAGE18_9_BOUNDED_TOL=1.0e-8
STAGE18_9_TEST_CASE=np1_dry_forced_bending_energy_bounded
```

## Manual command and expected PASS evidence

Run:

```bash
stage18_checks/run_stage18_9_controlled_one_fibre_physical_response_np1.sh
```

Expected terminal evidence:

```text
STAGE 18.9 CONTROLLED ONE FIBRE PHYSICAL RESPONSE NP1 VERDICT: PASS
STAGE 18.9 FINAL VERDICT: PASS
```

The wrapper writes
`stage18_outputs/fibre_stage18_9_controlled_one_fibre_physical_response_np1.dat`
with all required Stage 18.9 `*_status` fields and `final_status PASS`.
