# Stage 18.3: physical bending force operator

Stage 18.3 is a diagnostic-only physical bending-force operator gate for a single fibre.  It validates the bending-force candidate formulas `F_b = -B * X_ssss` and `F_b = -gamma * X_ssss` using Stage 18.2 geometry-operator concepts, and it validates the bending-energy formula.  The candidate force is never applied to acceleration, velocity, position, a structure integrator, fluid RHS, IBM, or production state.

## Geometry inputs

Stage 18.3 uses diagnostic geometry inputs only:

- `X(s)`
- `X_ss`
- `X_ssss`
- `ds`
- `npts`
- `fibre_length`
- `component_dim = 3`

## Bending parameters

Stage 18.3 validates these bending parameters and switches:

- `bending_stiffness B`
- `gamma`: nondimensional bending stiffness
- `use_dimensional_bending`
- `use_nondimensional_bending`

At least one of dimensional or nondimensional bending validation must be enabled.

## Bending-force candidate formulas

Dimensional bending-force candidate:

```text
F_b = -B * X_ssss
```

Nondimensional bending-force candidate:

```text
F_b = -gamma * X_ssss
```

Stage 18.3 may compute these candidate arrays for diagnostics.  It must not apply them.

## Bending-force norm diagnostics

The helper records only diagnostic candidate-array checks:

- `max_abs_bending_force_candidate`
- `l2_bending_force_candidate`
- `bending_force_candidate_finite`

These are diagnostic reductions, not force application.

## Bending energy

Bending energy boundary:

```text
E_b = 1/2 * int_0^L B * |X_ss|^2 ds
```

Discrete Stage 18.3 diagnostic approximation:

```text
E_b ~= 1/2 * sum_l B * |X_ss,l|^2 * w_l
```

with uniform trapezoidal weights:

```text
w_0 = ds/2
w_{n-1} = ds/2
w_l = ds for interior points
```

## Analytic validation cases

### `straight_fibre`

```text
X(s) = (s, 0, 0)
```

Expected:

```text
X_ssss = 0
F_b = 0
X_ss = 0
E_b = 0
```

### `sine_y_fibre`

```text
X(s) = (s, eps * sin(k*s), 0)
```

Expected:

```text
X_ss = (0, -eps*k^2*sin(k*s), 0)
X_ssss = (0, eps*k^4*sin(k*s), 0)
F_b = (0, -B*eps*k^4*sin(k*s), 0)
```

or nondimensional:

```text
F_b = (0, -gamma*eps*k^4*sin(k*s), 0)
```

### `quadratic_y_fibre`

```text
X(s) = (s, eps*s^2, 0)
```

Expected:

```text
X_ss = (0, 2*eps, 0)
X_ssss = 0
F_b = 0
E_b > 0
```

This case proves that bending energy can be nonzero even when the fourth-derivative bending-force candidate is zero.

## Candidate is not activation

Stage 18.3 explicitly distinguishes:

- `BENDING FORCE CANDIDATE / OPERATOR` means diagnostic computation of `-B X_ssss` or `-gamma X_ssss`.
- `BENDING FORCE ACTIVATED` means the force has been applied to acceleration, velocity, position, a structure integrator, fluid RHS, IBM, or production state.

Stage 18.3 must not:

- apply bending force to acceleration;
- update `V`;
- update `X`;
- insert the candidate force into any structure integrator;
- insert the candidate force into fluid RHS;
- solve tension;
- enforce inextensibility;
- perform structure time integration or structure state update;
- activate fluid-force physical-structure integration;
- activate runtime structure energy/power coupling;
- introduce wall-contact, fibre-fibre collision, penalty, repulsive, lubrication, friction, adhesion, or contact-damping forces;
- introduce production multi-fibre logic;
- introduce direct RHS injection, unapproved Stage 14 RHS calls, legacy IBM forcing, or unapproved production IBM forcing; or
- modify pressure, projection, Poisson, RK3, or channel-forcing logic.

## Wrapper defaults

`stage18_checks/run_stage18_3_physical_bending_force_operator.sh` infers the repository root from its own script path.  `DECOMP2D_ROOT` is retained only as an interface compatibility variable and is not used as the repository root.  The wrapper creates `stage18_outputs/`, invokes only the Stage 18.3 helper, builds nothing, runs no MPI, and runs no production validation.

Safe defaults are:

```text
STAGE18_3_ENABLE=1
STAGE18_3_BENDING_OPERATOR_ENABLE=1
STAGE18_3_SINGLE_FIBRE_ONLY=1
STAGE18_3_DIAGNOSTIC_ONLY=1
STAGE18_3_NPTS=32
STAGE18_3_FIBRE_LENGTH=1.0
STAGE18_3_COMPONENT_DIM=3
STAGE18_3_BENDING_STIFFNESS=1.0e-3
STAGE18_3_GAMMA=1.0e-3
STAGE18_3_USE_DIMENSIONAL_BENDING=1
STAGE18_3_USE_NONDIMENSIONAL_BENDING=1
STAGE18_3_SINE_EPS=1.0e-3
STAGE18_3_SINE_MODE=1
STAGE18_3_QUADRATIC_EPS=1.0e-3
STAGE18_3_ZERO_TOL=1.0e-14
STAGE18_3_FORMULA_TOL=1.0e-12
STAGE18_3_ENERGY_TOL=1.0e-12
STAGE18_3_TEST_CASE=straight_sine_quadratic_bending
```

## Closed-stage preservation

Stage 18.3 adds only new Stage 18.3 files.  It does not modify closed Stage 10--17 files, Stage 18.0 files, Stage 18.1 files, or Stage 18.2 files.  Stage 18.3 preserves the Stage 18.0 wrapper-root fix, Stage 18.1 physical-structure configuration evidence, Stage 18.2 structure-state/geometry-operator evidence, Stage 17.6 source-only/static-audit fix, Stage 17.10 evidence/static-audit fix, and Stage 17.11 total audit/closure fix.

## False-positive-safe audit policy

Stage 18.3 reuses the corrected Stage 18.2, Stage 18.1, Stage 18.0, Stage 17.11, Stage 17.10, Stage 17.6, and Stage 16 false-positive-safe audit pattern without editing closed files:

- use targeted structural checks rather than broad repository-wide scans;
- do not scan Markdown documentation as real code-regression evidence;
- do not treat documentation text, negative-check strings, or old failure labels as real regressions;
- do not treat regex literals such as `rg[[:space:]]` as actual `rg` command usage;
- do not require `rg`; any future shell use of `rg` must include a `grep` fallback;
- do not classify a source-only archive without `.git` metadata as DNS-core contamination or closed-stage modification;
- do not treat Stage 2 or legacy structure files as Stage 18.3 activation evidence;
- do not treat diagnostic computation of `F_b` as application of bending force; and
- do not treat `final_status` as missing before it is written.

Numeric and string evidence fields such as `*_value`, `*_formula_value`, `*_shape_value`, and `*_case_value` are not boolean PASS/FAIL fields.  Only `*_status` fields and explicit boolean fields control `final_status`.

## Expected PASS verdict

With Stage 18.2, Stage 18.1, and Stage 18.0 evidence intact; Stage 17 closure evidence intact or structurally accepted; corrected preservation fixes intact; and default valid bending parameters, the wrapper should print:

```text
STAGE 18.3 PHYSICAL BENDING FORCE OPERATOR VERDICT: PASS
STAGE 18.3 FINAL VERDICT: PASS
```

The helper writes `stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat` with all required Stage 18.3 summary fields and `final_status PASS` on the passing path.
