# Stage 18.4: tension solve and inextensibility constraint

Stage 18.4 is a diagnostic-only standalone tension-solve and inextensibility-constraint gate for one fibre.  It validates constraint diagnostics and a scalar tension equation in isolation, but it does not apply tension force, repair inextensibility, advance structure state, modify fluid RHS, touch IBM, or change DNS-core logic.

## Geometry and state inputs

Stage 18.4 uses diagnostic inputs only:

- `X(s)`
- `V(s)`
- `A(s)`
- `X_s`
- `V_s`
- `A_s`
- `X_ss`
- `X_sssss` as a diagnostic derivative for manufactured tension-equation checks only
- `ds`
- `npts`
- `fibre_length`
- `component_dim = 3`

## Constraint diagnostics

Inextensibility constraint:

```text
X_s dot X_s = 1
```

Arclength error:

```text
arc_error = X_s dot X_s - 1
```

Velocity-level constraint:

```text
X_s dot V_s = 0
```

Acceleration-level constraint:

```text
X_s dot A_s + V_s dot V_s = 0
```

Arclength-error detection is diagnostic-only.  Stage 18.4 does not project, repair, or otherwise change `X` or `V`.

## Structure-equation boundary

Dimensional structure-equation boundary for later stages:

```text
rho_l A = d_s(T X_s) - B X_ssss + F_h
```

Nondimensional boundary:

```text
rho_tilde A = d_s(T X_s) - gamma X_ssss + F_h
```

These equations are documented only; Stage 18.4 does not solve a coupled physical structure step.

## Standalone tension equation

Dimensional diagnostic tension equation:

```text
T_ss - |X_ss|^2 T = B X_s dot X_sssss - X_s dot d_s(F_h) - rho_l |V_s|^2
```

Nondimensional diagnostic tension equation:

```text
T_ss - |X_ss|^2 T = gamma X_s dot X_sssss - X_s dot d_s(F_h) - rho_tilde |V_s|^2
```

Stage 18.4 may solve this scalar equation as a standalone diagnostic linear system and check residuals.  It must not apply the resulting tension force:

```text
F_T = d_s(T X_s)
```

to acceleration, velocity, position, a structure integrator, fluid RHS, IBM, or production state.

## Diagnostic cases

### `straight_zero_tension`

```text
X(s) = (s, 0, 0)
V(s) = 0
F_h(s) = 0
```

Expected diagnostic evidence:

```text
X_s dot X_s = 1
X_s dot V_s = 0
V_s dot V_s = 0
|X_ss| = 0
RHS = 0
T = 0 with Dirichlet T(0)=T(L)=0
residual small
```

### `manufactured_tension_straight`

```text
X(s) = (s, 0, 0)
V(s) = 0
F_h(s) = 0
T_exact(s) = sin(pi*s/L)
T_ss = -(pi/L)^2 sin(pi*s/L)
```

The helper solves `T_ss = RHS` with Dirichlet endpoints and verifies bounded residual and bounded solution error.

### `velocity_constraint_case`

```text
X(s) = (s, 0, 0)
V(s) = (0, eps*sin(pi*s/L), 0)
```

Expected: `X_s dot V_s = 0`.  The acceleration-level source term `|V_s|^2` is finite and nonnegative, but Stage 18.4 does not update `A`.

### `arclength_error_detection`

```text
X(s) = ((1+eps)*s, 0, 0)
```

Expected: `X_s dot X_s - 1` is nonzero and detected.  No repair is attempted.

## Diagnostic solve is not activation

Stage 18.4 explicitly distinguishes:

- `TENSION SOLVE DIAGNOSTIC` means solving the scalar tension equation in isolation and checking residual.
- `TENSION FORCE ACTIVATED` means `d_s(T X_s)` has been applied to acceleration, velocity, position, a structure integrator, fluid RHS, IBM, or production state.
- `INEXTENSIBILITY DIAGNOSTIC` means measuring `X_s dot X_s - 1`.
- `INEXTENSIBILITY ENFORCEMENT` means changing `X` or `V` to satisfy the constraint.

Stage 18.4 must not:

- apply tension force to acceleration, velocity, or position;
- project or repair `X` or `V` for inextensibility;
- perform structure time integration or structure state update;
- activate fluid-force physical-structure integration;
- activate runtime structure energy/power implementation;
- apply bending force at runtime;
- modify any fluid RHS;
- introduce wall-contact, fibre-fibre collision, penalty, repulsive, lubrication, friction, adhesion, or contact-damping forces;
- introduce production multi-fibre logic; or
- introduce direct RHS injection, unapproved Stage 14 RHS calls, legacy IBM forcing, unapproved production IBM forcing, pressure/projection/Poisson/RK3/channel-forcing changes, or DNS-core changes.

## Wrapper defaults

`stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh` infers the repository root from its own script path.  `DECOMP2D_ROOT` is retained only as an interface compatibility variable and is not used as the repository root.  The wrapper creates `stage18_outputs/`, invokes only the Stage 18.4 helper, builds nothing, runs no MPI, and runs no production validation.

Safe defaults are:

```text
STAGE18_4_ENABLE=1
STAGE18_4_TENSION_CONSTRAINT_ENABLE=1
STAGE18_4_SINGLE_FIBRE_ONLY=1
STAGE18_4_DIAGNOSTIC_ONLY=1
STAGE18_4_NPTS=64
STAGE18_4_FIBRE_LENGTH=1.0
STAGE18_4_COMPONENT_DIM=3
STAGE18_4_RHO_L=1.0
STAGE18_4_RHO_TILDE=1.0
STAGE18_4_BENDING_STIFFNESS=1.0e-3
STAGE18_4_GAMMA=1.0e-3
STAGE18_4_VELOCITY_EPS=1.0e-3
STAGE18_4_ARCLENGTH_STRETCH_EPS=1.0e-3
STAGE18_4_ZERO_TOL=1.0e-14
STAGE18_4_FORMULA_TOL=1.0e-10
STAGE18_4_SOLVE_TOL=5.0e-3
STAGE18_4_ARC_ERROR_TOL=1.0e-6
STAGE18_4_TEST_CASE=straight_manufactured_velocity_arclength
```

## Closed-stage preservation

Stage 18.4 adds only new Stage 18.4 files.  It does not modify closed Stage 10--17 files, Stage 18.0 files, Stage 18.1 files, Stage 18.2 files, or Stage 18.3 files.  It preserves the Stage 18.0 wrapper-root fix, Stage 18.1 physical-structure configuration evidence, Stage 18.2 geometry-operator evidence, Stage 18.3 bending-operator evidence, Stage 17.6 source-only/static-audit fix, Stage 17.10 evidence/static-audit fix, and Stage 17.11 total audit/closure fix.

## False-positive-safe audit policy

Stage 18.4 reuses the corrected Stage 18.3, Stage 18.2, Stage 18.1, Stage 18.0, Stage 17.11, Stage 17.10, Stage 17.6, and Stage 16 false-positive-safe audit pattern without editing closed files:

- use targeted structural checks rather than broad repository-wide scans;
- do not scan Markdown documentation as real code-regression evidence;
- do not treat documentation text, negative-check strings, or old failure labels as real regressions;
- do not treat regex literals such as `rg[[:space:]]` as actual `rg` command usage;
- do not require `rg`; any future shell use of `rg` must include a `grep` fallback;
- do not classify a source-only archive without `.git` metadata as DNS-core contamination or closed-stage modification;
- do not treat Stage 2 or legacy structure files as Stage 18.4 activation evidence;
- do not treat a standalone diagnostic tension solve as tension-force application;
- do not treat arclength-error detection as inextensibility enforcement; and
- do not treat `final_status` as missing before it is written.

Numeric and string evidence fields such as `*_value`, `*_formula_value`, `*_shape_value`, and `*_case_value` are not boolean PASS/FAIL fields.  Only `*_status` fields and explicit boolean fields control `final_status`.

## Expected PASS verdict

With Stage 18.3, Stage 18.2, Stage 18.1, and Stage 18.0 evidence intact; Stage 17 closure evidence intact or structurally accepted; corrected preservation fixes intact; and default valid tension/inextensibility parameters, the wrapper should print:

```text
STAGE 18.4 TENSION INEXTENSIBILITY CONSTRAINT VERDICT: PASS
STAGE 18.4 FINAL VERDICT: PASS
```

The helper writes `stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat` with all required Stage 18.4 summary fields and `final_status PASS` on the passing path.
