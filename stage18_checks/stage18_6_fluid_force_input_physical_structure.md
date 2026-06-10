# Stage 18.6: fluid-on-fibre force input into physical structure

Stage 18.6 is a **diagnostic-only** standalone gate for the structure-side use of
an already available fluid-on-fibre force candidate.  It validates sign
convention, acceleration contribution, split-force bookkeeping, and local power
bookkeeping for a single fibre without modifying production state or fluid RHS.

## Force-input equation boundary

`F_h` is the force of fluid on fibre.  In the physical structure equation it
enters with positive structure-side sign:

```text
rho_l A = d_s(T X_s) - B X_ssss + F_h
rho_tilde A = d_s(T X_s) - gamma X_ssss + F_h
```

The local diagnostic acceleration candidates are:

```text
A_h_candidate = F_h / rho_l
A_h_candidate = F_h / rho_tilde
A_total_candidate = (F_T_candidate + F_b_candidate + F_h_candidate) / rho_l
A_total_candidate = (F_T_candidate + F_b_candidate + F_h_candidate) / rho_tilde
```

These arrays are local diagnostic candidates only.  They are not written to
production `X`, `V`, `A`, structure integrators, RHS, IBM, or DNS-core.

## Sign convention

The structure-side force is `+F_h`.  The equal-and-opposite bookkeeping force
that would correspond to fibre-on-fluid is:

```text
F_fibre_on_fluid = -F_h
F_h + F_fibre_on_fluid = 0
```

Stage 18.6 validates this sign convention only.  It does not spread
`F_fibre_on_fluid` to Eulerian RHS and does not call Stage 14 RHS injection.

## Power bookkeeping candidate

Stage 18.6 performs a local sign check for:

```text
P_h = int F_h dot V ds
P_h ~= sum_l F_h,l dot V_l * w_l
```

with uniform trapezoidal weights:

```text
w_0 = ds/2
w_{n-1} = ds/2
w_l = ds for interior points
```

This is only power bookkeeping.  Stage 18.6 does not activate runtime
energy/power production diagnostics; that boundary is reserved for later stages.

## Diagnostic cases

The helper validates:

1. `zero_fluid_force_case`: `F_h = 0`, so `A_h_candidate = 0` and `P_h = 0`.
2. `uniform_fluid_force_case`: constant `F_h` gives `A_h_candidate = F_h/rho_l`
   and, when enabled, `F_h/rho_tilde`.
3. `sign_convention_case`: structure-side sign is `+F_h`, fluid-side bookkeeping
   sign is `-F_h`, and action-reaction sum is zero.
4. `split_structure_rhs_case`: `F_total = F_T_candidate + F_b_candidate +
   F_h_candidate` and `A_total_candidate = F_total/rho`.
5. `power_sign_case`: parallel velocity gives positive power, anti-parallel gives
   negative power, and perpendicular gives zero power within tolerance.

## Wrapper defaults

```text
STAGE18_6_ENABLE=1
STAGE18_6_FLUID_FORCE_INPUT_ENABLE=1
STAGE18_6_SINGLE_FIBRE_ONLY=1
STAGE18_6_DIAGNOSTIC_ONLY=1
STAGE18_6_NPTS=16
STAGE18_6_FIBRE_LENGTH=1.0
STAGE18_6_COMPONENT_DIM=3
STAGE18_6_RHO_L=1.0
STAGE18_6_RHO_TILDE=1.0
STAGE18_6_USE_DIMENSIONAL_MASS=1
STAGE18_6_USE_NONDIMENSIONAL_MASS=1
STAGE18_6_FLUID_FORCE_MAG=1.0e-3
STAGE18_6_VELOCITY_MAG=1.0e-3
STAGE18_6_ZERO_TOL=1.0e-14
STAGE18_6_FORMULA_TOL=1.0e-12
STAGE18_6_POWER_TOL=1.0e-12
STAGE18_6_TEST_CASE=zero_uniform_sign_split_power
```

`DECOMP2D_ROOT` remains an interface-compatibility variable only.  The wrapper
infers the repository root from its own script path and does not `cd` into
`DECOMP2D_ROOT`.

## Non-activation guarantees

Stage 18.6 must not and does not:

- modify production structure force input or production `X`, `V`, `A` state;
- insert a production structure hook;
- modify Stage 16 structure-force code;
- modify Stage 13 force-density code;
- modify Stage 14 RHS code or call Stage 14 RHS injection;
- spread force to fluid RHS;
- modify fluid RHS, IBM, pressure projection, Poisson, RK3/channel forcing, or
  DNS-core;
- run structure time integration;
- apply bending or tension forces to production state;
- project or repair `X` or `V` for inextensibility;
- activate runtime energy/power production diagnostics;
- introduce contact/collision/wall/penalty/repulsive/lubrication/friction/
  adhesion/contact-damping forces;
- introduce production multi-fibre logic;
- build targets, run MPI, or run production validation.

Protective documentation references to `src/fibre_stage14_production_rhs_injection.f90`
and `src/xcompact3d.f90` are not runtime activation.

## False-positive-safe audit policy

The Stage 18.6 helper continues the corrected Stage 18.5 / Stage 18.4 / Stage
18.3 / Stage 18.2 / Stage 18.1 / Stage 18.0 / Stage 17.11 / Stage 17.10 / Stage
17.6 / Stage 16 policy:

- targeted structural checks only;
- no broad repository-wide code scans;
- Markdown examples and negative-check strings are not executable regressions;
- regex literals such as `rg[[:space:]]` are not treated as `rg` usage;
- no mandatory `rg` dependency;
- source-only archives without `.git` metadata are accepted as non-contamination;
- existing Stage 2 or legacy structure files are not Stage 18.6 activation;
- local diagnostic `F_h` and `-F_h` arrays are not production force input or RHS
  contamination unless injected into production paths;
- only `*_status` fields, not `*_value`, `*_formula_value`, `*_shape_value`, or
  `*_case_value` fields, control `final_status`.

## Manual command and expected evidence

Run:

```bash
stage18_checks/run_stage18_6_fluid_force_input_physical_structure.sh
```

Expected wrapper evidence:

```text
STAGE 18.6 FLUID FORCE INPUT PHYSICAL STRUCTURE VERDICT: PASS
STAGE 18.6 FINAL VERDICT: PASS
```

The helper writes:

```text
stage18_outputs/fibre_stage18_6_fluid_force_input_physical_structure.dat
```

with all required Stage 18.6 `*_status` fields and `final_status PASS`.
