# Stage 18.7: structure energy and power diagnostics

Stage 18.7 is a **diagnostic-only** standalone gate for local single-fibre
energy and power bookkeeping.  It validates kinetic energy, bending energy, total
structure energy, fluid-on-fibre power, and an energy-residual formula using only
helper-local arrays.

## Energy and power equations

Dimensional kinetic energy:

```text
E_k = 1/2 * int_0^L rho_l * |V|^2 ds
```

Nondimensional kinetic energy:

```text
E_k = 1/2 * int_0^L rho_tilde * |V|^2 ds
```

Dimensional and nondimensional bending energy:

```text
E_b = 1/2 * int_0^L B * |X_ss|^2 ds
E_b = 1/2 * int_0^L gamma * |X_ss|^2 ds
```

Total structure energy and fluid-on-fibre power:

```text
E_s = E_k + E_b
P_h = int_0^L F_h dot V ds
```

Trapezoidal quadrature weights are:

```text
w_0 = ds/2
w_{n-1} = ds/2
w_i = ds for interior points
```

The discrete diagnostic formulas are:

```text
E_k ~= 1/2 * sum_i rho_l * |V_i|^2 * w_i
E_b ~= 1/2 * sum_i B * |X_ss_i|^2 * w_i
E_s = E_k + E_b
P_h ~= sum_i F_h_i dot V_i * w_i
R_E = (E_s^{n+1} - E_s^n)/dt - P_h
```

`R_E` is a local diagnostic residual formula only.  Stage 18.7 does not add
damping, runtime energy output, production energy streams, or full production
energy-budget logic.

## Diagnostic cases

The helper validates:

1. `straight_rest_zero_energy_case`: straight fibre, zero velocity, zero
   curvature, and zero fluid force give `E_k = 0`, `E_b = 0`, `E_s = 0`, and
   `P_h = 0`.
2. `straight_uniform_velocity_case`: constant velocity gives
   `E_k = 1/2*rho_l*U0^2*L`, zero bending energy, and `E_s = E_k`.
3. `sine_bending_energy_case`: sinusoidal centreline curvature gives zero kinetic
   energy and positive bending energy.
4. `fluid_power_sign_case`: `F_h` parallel to `V` gives positive power,
   anti-parallel gives negative power, and perpendicular gives zero power.
5. `energy_residual_formula_case`: local scalars validate
   `R_E = (E_s^{n+1} - E_s^n)/dt - P_h`.
6. `candidate_update_energy_case`: Stage 18.5-style candidate arrays are used only
   locally to confirm energy values are finite and nonnegative.

## Configuration defaults

```text
STAGE18_7_ENABLE=1
STAGE18_7_ENERGY_POWER_DIAGNOSTIC_ENABLE=1
STAGE18_7_SINGLE_FIBRE_ONLY=1
STAGE18_7_DIAGNOSTIC_ONLY=1
STAGE18_7_NPTS=64
STAGE18_7_FIBRE_LENGTH=1.0
STAGE18_7_COMPONENT_DIM=3
STAGE18_7_RHO_L=1.0
STAGE18_7_RHO_TILDE=1.0
STAGE18_7_BENDING_STIFFNESS=1.0e-3
STAGE18_7_GAMMA=1.0e-3
STAGE18_7_USE_DIMENSIONAL_ENERGY=1
STAGE18_7_USE_NONDIMENSIONAL_ENERGY=1
STAGE18_7_VELOCITY_MAG=1.0e-3
STAGE18_7_FLUID_FORCE_MAG=1.0e-3
STAGE18_7_SINE_EPS=1.0e-3
STAGE18_7_SINE_MODE=1
STAGE18_7_DT_STRUCTURE=1.0e-4
STAGE18_7_ZERO_TOL=1.0e-14
STAGE18_7_FORMULA_TOL=1.0e-12
STAGE18_7_POWER_TOL=1.0e-12
STAGE18_7_ENERGY_TOL=1.0e-12
STAGE18_7_TEST_CASE=straight_rest_uniform_sine_power_residual
```

`DECOMP2D_ROOT` remains interface compatibility only.  The wrapper infers the
repository root from its own script path and does not `cd` into `DECOMP2D_ROOT`.

## Non-activation guarantees

Stage 18.7 must not and does not:

- write production energy diagnostic output;
- modify production `X`, `V`, or `A` state;
- insert production structure hooks;
- modify Stage 16 structure-force code;
- modify `src/fibre_stage13_production_force_density_candidate.f90`;
- modify `src/fibre_stage14_production_rhs_injection.f90` or call Stage 14 RHS
  injection;
- spread force to fluid RHS or modify fluid RHS;
- modify IBM, pressure projection, Poisson, RK3/channel forcing, or DNS-core;
- modify statistics, visualisation, restart, or runtime output I/O;
- run structure time integration runtime;
- apply bending or tension forces to production state;
- project or repair `X` or `V` for inextensibility;
- introduce contact/collision/wall/penalty/repulsive/lubrication/friction/
  adhesion/contact-damping forces;
- introduce production multi-fibre logic;
- build targets, run MPI, or run production validation.

Protective documentation references to `src/fibre_stage14_production_rhs_injection.f90`,
`src/fibre_stage13_production_force_density_candidate.f90`, and `src/xcompact3d.f90`
are not runtime activation.

## False-positive-safe audit policy

The Stage 18.7 helper continues the corrected Stage 18.6 / Stage 18.5 / Stage
18.4 / Stage 18.3 / Stage 18.2 / Stage 18.1 / Stage 18.0 / Stage 17.11 / Stage
17.10 / Stage 17.6 / Stage 16 policy:

- targeted structural checks only;
- no broad repository-wide code scans;
- Markdown examples and negative-check strings are not executable regressions;
- regex literals such as `rg[[:space:]]` are not treated as `rg` usage;
- no mandatory `rg` dependency;
- source-only archives without `.git` metadata are accepted as non-contamination;
- existing Stage 2 or legacy structure files are not Stage 18.7 activation;
- local diagnostic `E_k`, `E_b`, `E_s`, `P_h`, and `R_E` values are not production
  energy output unless written into runtime output, stats/visu/restart I/O,
  production modules, hooks, RHS, IBM, or Fortran production paths;
- local diagnostic `F_h dot V` power bookkeeping is not RHS contamination;
- only `*_status` fields, not `*_value`, `*_formula_value`, `*_shape_value`, or
  `*_case_value` fields, control `final_status`.

## Manual command and expected evidence

Run:

```bash
stage18_checks/run_stage18_7_structure_energy_power_diagnostics.sh
```

Expected wrapper evidence:

```text
STAGE 18.7 STRUCTURE ENERGY POWER DIAGNOSTICS VERDICT: PASS
STAGE 18.7 FINAL VERDICT: PASS
```

The helper writes:

```text
stage18_outputs/fibre_stage18_7_structure_energy_power_diagnostics.dat
```

with all required Stage 18.7 `*_status` fields and `final_status PASS`.
