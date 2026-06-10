# Stage 18.5: structure time integration core

Stage 18.5 is a **diagnostic-only** standalone gate for the single-fibre
structure time-integration core.  It validates local candidate arrays and update
formulas that later Stage 18 production work may use, while preserving all
closed Stage 10--17 and Stage 18.0--18.4 files.

## Scope boundary

Stage 18.5 defines the following local diagnostic state variables only:

- `X^n`: centreline position at the diagnostic time level.
- `V^n`: centreline velocity at the diagnostic time level.
- `A_candidate^n`: diagnostic acceleration candidate.
- `X_candidate^{n+1}`: diagnostic next-position candidate.
- `V_candidate^{n+1}`: diagnostic next-velocity candidate.
- `dt_structure`, `rho_l`, `rho_tilde`, `npts`, and `component_dim = 3`.

The diagnostic structure equation boundary is:

```text
X_t = V
rho_l V_t = F_total
rho_tilde V_t = F_total
F_total = F_T_candidate + F_b_candidate + F_h_candidate
```

where `F_T_candidate`, `F_b_candidate`, and `F_h_candidate` are local diagnostic
candidate-force arrays only.  The candidate force sum is never spread to fluid
RHS, IBM, DNS-core, production structure modules, or persistent production state.

## Candidate update formulas

Stage 18.5 validates a velocity-Verlet / constant-acceleration diagnostic update:

```text
A_candidate^n = F_total^n / rho_l
X_candidate^{n+1} = X^n + dt * V^n + 0.5 * dt^2 * A_candidate^n
V_candidate^{n+1} = V^n + dt * A_candidate^n
```

and the nondimensional companion:

```text
A_candidate^n = F_total^n / rho_tilde
X_candidate^{n+1} = X^n + dt * V^n + 0.5 * dt^2 * A_candidate^n
V_candidate^{n+1} = V^n + dt * A_candidate^n
```

This is a local **candidate update** inside the Stage 18.5 helper only.

## Diagnostic cases

The helper validates these cases without running production physics:

1. `zero_force_rest_case`: straight fibre, zero velocity, zero force, zero drift.
2. `uniform_velocity_no_force_case`: zero acceleration, velocity preserved, and
   `X_candidate^{n+1} = X^n + dt * V^n`.
3. `constant_force_case`: `A_candidate = F_total / rho_l`, with velocity and
   position formulas checked directly.
4. `split_force_sum_case`: `F_total = F_T_candidate + F_b_candidate +
   F_h_candidate`, checked before computing diagnostic acceleration.
5. `dt_refinement_case`: for constant acceleration, two half steps match one full
   step within tolerance.

## Configuration interface

The wrapper supports these environment variables with safe defaults:

```text
STAGE18_5_ENABLE=1
STAGE18_5_TIME_INTEGRATION_CORE_ENABLE=1
STAGE18_5_SINGLE_FIBRE_ONLY=1
STAGE18_5_DIAGNOSTIC_ONLY=1
STAGE18_5_NPTS=16
STAGE18_5_FIBRE_LENGTH=1.0
STAGE18_5_COMPONENT_DIM=3
STAGE18_5_DT_STRUCTURE=1.0e-4
STAGE18_5_RHO_L=1.0
STAGE18_5_RHO_TILDE=1.0
STAGE18_5_USE_DIMENSIONAL_MASS=1
STAGE18_5_USE_NONDIMENSIONAL_MASS=1
STAGE18_5_UNIFORM_VELOCITY=1.0e-3
STAGE18_5_CONSTANT_FORCE=1.0e-3
STAGE18_5_ZERO_TOL=1.0e-14
STAGE18_5_FORMULA_TOL=1.0e-12
STAGE18_5_DT_REFINEMENT_TOL=1.0e-12
STAGE18_5_TEST_CASE=zero_uniform_constant_split_dt_refinement
```

`DECOMP2D_ROOT` remains an interface-compatibility variable only.  The wrapper
infers the repository root from its own script path and does not `cd` into
`DECOMP2D_ROOT`.

## Non-activation guarantees

Stage 18.5 must not and does not:

- modify production `X`, `V`, or `A` state;
- insert a production structure hook;
- apply bending or tension force to production state;
- project or repair `X` or `V` for inextensibility;
- activate fluid-force physical-structure integration;
- call Stage 14 RHS injection;
- modify fluid RHS, IBM, pressure projection, Poisson, RK3/channel forcing, or
  DNS-core;
- introduce contact, collision, wall, penalty, repulsive, lubrication, friction,
  adhesion, or damping forces;
- introduce production multi-fibre logic;
- build targets, run MPI, or run production validation.

## False-positive-safe audit policy

The Stage 18.5 helper continues the corrected Stage 18.4 / Stage 18.3 / Stage
18.2 / Stage 18.1 / Stage 18.0 / Stage 17.11 / Stage 17.10 / Stage 17.6 / Stage
16 policy:

- targeted structural checks only;
- no broad repository-wide code scans;
- Markdown examples are documentation, not executable regression evidence;
- negative-check strings are not treated as activation;
- regex literals such as `rg[[:space:]]` are not treated as `rg` usage;
- no mandatory `rg` dependency;
- source-only archives without `.git` metadata are accepted as non-contamination;
- existing Stage 2 or legacy structure files are not treated as Stage 18.5
  activation;
- only `*_status` fields, not `*_value`, `*_formula_value`, `*_shape_value`, or
  `*_case_value` fields, control `final_status`.

## Manual command and expected evidence

Run:

```bash
stage18_checks/run_stage18_5_structure_time_integration_core.sh
```

Expected wrapper evidence:

```text
STAGE 18.5 STRUCTURE TIME INTEGRATION CORE VERDICT: PASS
STAGE 18.5 FINAL VERDICT: PASS
```

The helper writes:

```text
stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat
```

with all required Stage 18.5 `*_status` fields and `final_status PASS`.

## Stage 18.5 false-positive-safe runtime-audit note

Stage 18.5 is allowed to mention protected production filenames only inside
read-only evidence-preservation checks, for example when confirming that the
previous Stage 14 small-lambda hook evidence still points to the correct files.
Such literals are not runtime activation.  Runtime activation is defined only as
executing a build/MPI/production command, inserting a production hook, writing to
protected source/CMake/closed-stage paths, or feeding candidate updates into
production X/V/A, RHS, IBM, or DNS-core state.

Therefore the Stage 18.5 helper audits executable behaviour and changed paths,
not every protected filename literal in its own source.  This preserves the
Stage 17.6 / 17.10 / 17.11 / 18.0 false-positive-safe source-only audit policy.

