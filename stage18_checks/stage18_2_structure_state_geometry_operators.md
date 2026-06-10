# Stage 18.2: structure state and geometry operators

Stage 18.2 is a diagnostic-only structure-state and geometry-operator gate for later single-fibre physical structure dynamics.  It defines and validates in-memory state concepts and analytic geometry-operator diagnostics, but it does not activate physical force, tension, inextensibility, time integration, structure-state update, fluid-force physical-structure coupling, RHS/IBM logic, contact/collision logic, or DNS-core changes.

## State concepts

Stage 18.2 defines the single-fibre state concepts used by later Stage 18 work:

- `X(s,t)`: fibre centreline position
- `V(s,t)`: fibre velocity
- `A(s,t)`: fibre acceleration
- `npts`: number of Lagrangian points
- `ds`: nominal arclength spacing
- `component_dim = 3`
- single-fibre state only

Velocity and acceleration are state concepts only in Stage 18.2:

```text
V = X_t
A = V_t
```

No time integration is performed and no `X`, `V`, or `A` state update is applied.

## Geometry operators

Stage 18.2 validates the following diagnostic geometry operators:

```text
X_s    = dX/ds
X_ss   = d^2X/ds^2
X_sss  = d^3X/ds^3
X_ssss = d^4X/ds^4
```

Tangent:

```text
tau = X_s
```

Curvature vector:

```text
kappa_vec = X_ss
```

Curvature magnitude:

```text
kappa = |X_ss|
```

Arclength constraint diagnostic:

```text
arc_error = X_s dot X_s - 1
```

`max_abs_arc_error` is a diagnostic reduction only.  Endpoint/free-free geometry metadata is diagnostic-only and applies no endpoint or boundary force.

## Straight-fibre diagnostic

For a straight fibre aligned with `x`:

```text
X(s) = (s, 0, 0)
```

Expected diagnostic geometry:

```text
X_s = (1,0,0)
X_ss = 0
X_sss = 0
X_ssss = 0
kappa = 0
arc_error = 0
```

## Sinusoidal diagnostic

For a small-amplitude sinusoidal diagnostic fibre:

```text
X(s) = (s, eps*sin(k*s), 0)
```

Expected analytic derivatives:

```text
X_s = (1, eps*k*cos(k*s), 0)
X_ss = (0, -eps*k^2*sin(k*s), 0)
X_sss = (0, -eps*k^3*cos(k*s), 0)
X_ssss = (0, eps*k^4*sin(k*s), 0)
```

These derivatives are geometry-operator diagnostics only.

## Geometry is not physical-force activation

Stage 18.2 explicitly distinguishes:

- `GEOMETRY OPERATOR AVAILABLE` means diagnostic derivative, curvature, and arclength calculations exist.
- `PHYSICAL FORCE ACTIVATED` means a force or state update has been formed or applied.

Stage 18.2 may validate `X_ssss` as a geometric derivative, but it must not form or apply:

```text
F_b = -B X_ssss
```

Stage 18.2 must not:

- compute a real bending force;
- solve tension;
- enforce inextensibility;
- advance `X`, `V`, or `A`;
- insert any production hook;
- modify any fluid RHS;
- implement structure energy or power diagnostics;
- introduce wall-contact, fibre-fibre collision, penalty, repulsive, lubrication, friction, adhesion, or contact-damping forces;
- introduce production multi-fibre logic;
- introduce collision-induced RHS or collision-induced structure update;
- introduce direct RHS injection, unapproved Stage 14 RHS calls, legacy IBM forcing, or unapproved production IBM forcing; or
- modify pressure, projection, Poisson, RK3, or channel-forcing logic.

## Wrapper defaults

`stage18_checks/run_stage18_2_structure_state_geometry_operators.sh` infers the repository root from its own script path.  `DECOMP2D_ROOT` is retained only as an interface compatibility variable and is not used as the repository root.  The wrapper creates `stage18_outputs/`, invokes only the Stage 18.2 helper, builds nothing, runs no MPI, and runs no production validation.

Safe defaults are:

```text
STAGE18_2_ENABLE=1
STAGE18_2_STRUCTURE_STATE_GEOMETRY_ENABLE=1
STAGE18_2_SINGLE_FIBRE_ONLY=1
STAGE18_2_DIAGNOSTIC_ONLY=1
STAGE18_2_NPTS=16
STAGE18_2_FIBRE_LENGTH=1.0
STAGE18_2_COMPONENT_DIM=3
STAGE18_2_SINE_EPS=1.0e-3
STAGE18_2_SINE_MODE=1
STAGE18_2_ZERO_TOL=1.0e-14
STAGE18_2_FORMULA_TOL=1.0e-12
STAGE18_2_DERIVATIVE_TOL=5.0e-2
STAGE18_2_ARC_ERROR_TOL=1.0e-2
STAGE18_2_TEST_CASE=straight_and_sine_geometry
```

## Closed-stage preservation

Stage 18.2 adds only new Stage 18.2 files.  It does not modify closed Stage 10--17 files, Stage 18.0 files, or Stage 18.1 files.  Stage 18.2 preserves the Stage 18.0 wrapper-root fix, Stage 18.1 physical-structure configuration evidence, Stage 17.6 source-only/static-audit fix, Stage 17.10 evidence/static-audit fix, and Stage 17.11 total audit/closure fix.

## False-positive-safe audit policy

Stage 18.2 reuses the corrected Stage 18.1, Stage 18.0, Stage 17.11, Stage 17.10, Stage 17.6, and Stage 16 false-positive-safe audit pattern without editing closed files:

- use targeted structural checks rather than broad repository-wide scans;
- do not scan Markdown documentation as real code-regression evidence;
- do not treat documentation text, negative-check strings, or old failure labels as real regressions;
- do not treat regex literals such as `rg[[:space:]]` as actual `rg` command usage;
- do not require `rg`; any future shell use of `rg` must include a `grep` fallback;
- do not classify a source-only archive without `.git` metadata as DNS-core contamination or closed-stage modification;
- do not treat Stage 2 or legacy structure files as Stage 18.2 activation evidence;
- do not treat the diagnostic computation of `X_ssss` as bending-force activation; and
- do not treat `final_status` as missing before it is written.

Numeric and string evidence fields such as `*_value`, `*_formula_value`, `*_shape_value`, and `*_case_value` are not boolean PASS/FAIL fields.  Only `*_status` fields and explicit boolean fields control `final_status`.

## Expected PASS verdict

With Stage 18.1 evidence intact, Stage 18.0 evidence intact, Stage 17 closure evidence intact or structurally accepted, corrected preservation fixes intact, and default valid geometry parameters, the wrapper should print:

```text
STAGE 18.2 STRUCTURE STATE GEOMETRY OPERATORS VERDICT: PASS
STAGE 18.2 FINAL VERDICT: PASS
```

The helper writes `stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat` with all required Stage 18.2 summary fields and `final_status PASS` on the passing path.
