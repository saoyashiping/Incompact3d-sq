# Stage 18.8: dry physical structure benchmark

Stage 18.8 is a **diagnostic-only** standalone dry benchmark for one fibre with
`F_h = 0`.  It validates dry no-fluid bookkeeping, straight-rest/no-drift,
uniform-translation/no-force, bending-only sine candidate acceleration, bounded
local candidate steps, and dry energy checks using helper-local arrays only.

## Dry structure equations

Dry condition:

```text
F_h = 0
```

Dimensional and nondimensional dry structure equations:

```text
rho_l X_tt = d_s(T X_s) - B X_ssss
rho_tilde X_tt = d_s(T X_s) - gamma X_ssss
```

Local dry acceleration candidates:

```text
A_dry_candidate = (F_T_candidate + F_b_candidate) / rho_l
A_dry_candidate = (F_T_candidate + F_b_candidate) / rho_tilde
```

Bending-only diagnostic candidates use:

```text
F_T_candidate = 0
F_b_candidate = -B X_ssss
A_b_candidate = -B X_ssss / rho_l
A_b_candidate = -gamma X_ssss / rho_tilde
```

Dry energy and no-fluid power are:

```text
E_k = 1/2 int rho_l |V|^2 ds
E_b = 1/2 int B |X_ss|^2 ds
E_s = E_k + E_b
P_h = int F_h dot V ds = 0
R_E_dry = (E_s^{n+1} - E_s^n)/dt
```

`R_E_dry` is a local benchmark value only; it is not production runtime output.

## Diagnostic cases

The helper validates:

1. `dry_straight_rest_no_drift_case`: straight fibre, zero velocity, zero
   candidate forces, zero acceleration, no position/velocity drift, zero energy,
   and zero no-fluid power.
2. `dry_uniform_translation_case`: straight fibre with constant velocity, zero dry
   acceleration, velocity preservation, `X_next = X + dt * V`, correct kinetic
   energy, zero bending energy, and zero no-fluid power.
3. `dry_sine_bending_candidate_case`: sinusoidal fibre with `X_ssss = (0,
   eps*k^4*sin(k*s), 0)`, `F_b_candidate = -B X_ssss`, `A_b_candidate =
   F_b_candidate/rho_l`, acceleration opposing displacement, positive bending
   energy, and zero no-fluid power.
4. `dry_candidate_time_step_bounded_case`: Stage 18.5-style local candidate
   updates have finite arrays, bounded displacement increments, and bounded
   velocity increments.
5. `dry_energy_bounded_case`: local dry energies are finite and nonnegative.
6. `dry_no_fluid_contamination_case`: `F_h = 0`, `P_h = 0`, no Stage 14 RHS call,
   no force spreading, and no RHS/IBM/DNS-core modification.

## Configuration defaults

```text
STAGE18_8_ENABLE=1
STAGE18_8_DRY_STRUCTURE_BENCHMARK_ENABLE=1
STAGE18_8_SINGLE_FIBRE_ONLY=1
STAGE18_8_DIAGNOSTIC_ONLY=1
STAGE18_8_NPTS=64
STAGE18_8_FIBRE_LENGTH=1.0
STAGE18_8_COMPONENT_DIM=3
STAGE18_8_RHO_L=1.0
STAGE18_8_RHO_TILDE=1.0
STAGE18_8_BENDING_STIFFNESS=1.0e-3
STAGE18_8_GAMMA=1.0e-3
STAGE18_8_USE_DIMENSIONAL_DRY=1
STAGE18_8_USE_NONDIMENSIONAL_DRY=1
STAGE18_8_VELOCITY_MAG=1.0e-3
STAGE18_8_SINE_EPS=1.0e-3
STAGE18_8_SINE_MODE=1
STAGE18_8_DT_STRUCTURE=1.0e-4
STAGE18_8_MAX_DISPLACEMENT_INCREMENT=1.0e-4
STAGE18_8_MAX_VELOCITY_INCREMENT=1.0e-3
STAGE18_8_ZERO_TOL=1.0e-14
STAGE18_8_FORMULA_TOL=1.0e-12
STAGE18_8_ENERGY_TOL=1.0e-12
STAGE18_8_BOUNDED_TOL=1.0e-8
STAGE18_8_TEST_CASE=dry_straight_translation_sine_bounded_energy
```

`DECOMP2D_ROOT` remains interface compatibility only.  The wrapper infers the
repository root from its own script path and does not `cd` into `DECOMP2D_ROOT`.

## Non-activation guarantees

Stage 18.8 must not and does not:

- write production dry benchmark output;
- modify production `X`, `V`, or `A` state;
- insert production structure hooks;
- modify Stage 16 structure-force code;
- modify `src/fibre_stage13_production_force_density_candidate.f90`;
- modify `src/fibre_stage14_production_rhs_injection.f90` or call Stage 14 RHS
  injection;
- spread force to fluid RHS or modify fluid RHS;
- modify IBM, pressure projection, Poisson, RK3/channel forcing, or DNS-core;
- modify `statistics.f90`, `visu.f90`, restart I/O, or runtime output I/O;
- run production structure time integration;
- apply bending or tension forces to production state;
- project or repair `X` or `V` for inextensibility;
- introduce contact/collision/wall/penalty/repulsive/lubrication/friction/
  adhesion/contact-damping forces;
- introduce production multi-fibre logic;
- build targets, run MPI, or run production validation.

Protective documentation references to `src/fibre_stage14_production_rhs_injection.f90`,
`src/fibre_stage13_production_force_density_candidate.f90`, `src/xcompact3d.f90`,
`statistics.f90`, `visu.f90`, or restart I/O files are not runtime activation.

## False-positive-safe audit policy

The Stage 18.8 helper continues the corrected Stage 18.7 / Stage 18.6 / Stage
18.5 / Stage 18.4 / Stage 18.3 / Stage 18.2 / Stage 18.1 / Stage 18.0 / Stage
17.11 / Stage 17.10 / Stage 17.6 / Stage 16 policy:

- targeted structural checks only;
- no broad repository-wide code scans;
- Markdown examples and negative-check strings are not executable regressions;
- regex literals such as `rg[[:space:]]` are not treated as `rg` usage;
- no mandatory `rg` dependency;
- source-only archives without `.git` metadata are accepted as non-contamination;
- existing Stage 2 or legacy structure files are not Stage 18.8 activation;
- local diagnostic dry `X_candidate`, `V_candidate`, and `A_candidate` arrays are
  not production updates unless written into runtime output, production modules,
  hooks, RHS, IBM, or Fortran production paths;
- local dry energy computations are not production diagnostics or stats/visu/
  restart output;
- only `*_status` fields, not `*_value`, `*_formula_value`, `*_shape_value`, or
  `*_case_value` fields, control `final_status`.

## Manual command and expected evidence

Run:

```bash
stage18_checks/run_stage18_8_dry_physical_structure_benchmark.sh
```

Expected wrapper evidence:

```text
STAGE 18.8 DRY PHYSICAL STRUCTURE BENCHMARK VERDICT: PASS
STAGE 18.8 FINAL VERDICT: PASS
```

The helper writes:

```text
stage18_outputs/fibre_stage18_8_dry_physical_structure_benchmark.dat
```

with all required Stage 18.8 `*_status` fields and `final_status PASS`.
