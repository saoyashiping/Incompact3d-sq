# Stage 18.1: physical structure dynamics configuration

Stage 18.1 is a diagnostic-only configuration gate for later single-fibre physical structure dynamics.  It defines and validates physical and nondimensional parameter relations, but it does not activate structure dynamics, advance structure state, compute forces, solve constraints, modify the fluid RHS, or touch DNS-core logic.

## Configuration enable flags

Stage 18.1 distinguishes configuration availability from physics activation:

- `CONFIG ENABLED` means the named configuration group is present and validated.
- `PHYSICS ACTIVATED` means an actual solver, force, constraint enforcement, state update, RHS update, or diagnostic implementation has been applied.

Stage 18.1 may set these configuration gates to enabled:

- `stage18_enable`
- `stage18_structure_dynamics_config_enable`
- `stage18_single_fibre_only`
- `stage18_physical_structure_boundary`
- `stage18_bending_config_enable`
- `stage18_tension_config_enable`
- `stage18_inextensibility_config_enable`
- `stage18_time_integration_config_enable`
- `stage18_energy_diagnostic_config_enable`
- `stage18_diagnostic_only`

All physics activation remains absent in Stage 18.1.

## Physical and nondimensional parameters

Stage 18.1 validates these scalar configuration concepts for later stages:

- `rho_s`: fibre material density
- `fibre_radius`: fibre radius `a`
- `fibre_area`: cross-sectional area `A`
- `young_modulus`: Young's modulus `E`
- `second_moment_area`: second moment of area `I`
- `bending_stiffness`: `B = E * I`
- `mass_per_length`: `rho_l = rho_s * A`
- `fibre_length`: `L_f`
- `npts`: number of Lagrangian points
- `ds`: nominal arclength spacing
- `dt_structure`: structure time step for later stages
- `rho_tilde`: nondimensional mass parameter
- `gamma`: nondimensional bending stiffness parameter

## Formula boundary validated in Stage 18.1

Cross-sectional area:

```text
A = pi * a^2
```

Second moment of area for a circular cross-section:

```text
I = pi * a^4 / 4
```

Dimensional bending stiffness:

```text
B = E * I
```

Mass per unit length:

```text
rho_l = rho_s * A
```

Nominal arclength spacing:

```text
ds = L_f / (npts - 1)
```

Stage 18.1 validates only these parameter relations and finite/positive constraints.  It does not discretize or apply any force operator.

## Later Stage 18 physical equation boundary

Later Stage 18 substages target the single-fibre structure equation:

```text
rho_l * X_tt = d_s(T X_s) - B X_ssss + F_h
```

and the nondimensional form:

```text
rho_tilde * X_tt = d_s(T X_s) - gamma X_ssss + F_h
```

Later stages also target the inextensibility constraint:

```text
X_s dot X_s = 1
```

Later structure energy diagnostics target:

```text
E_k = 1/2 int rho_l |V|^2 ds
E_b = 1/2 int B |X_ss|^2 ds
P_h = int F_h dot V ds
```

These equations and diagnostics are boundary declarations in Stage 18.1, not implementations.

## Stage 18.1 non-activation guarantees

Stage 18.1 must not:

- implement the structure equation;
- compute a real bending force;
- solve tension;
- enforce inextensibility;
- advance `X`, `V`, or `A`;
- modify any fluid RHS;
- implement structure energy or power diagnostics;
- introduce wall-contact, fibre-fibre collision, penalty, repulsive, lubrication, friction, adhesion, or contact-damping forces;
- introduce production multi-fibre logic;
- introduce collision-induced RHS or collision-induced structure updates;
- introduce direct RHS injection, unapproved Stage 14 RHS calls, legacy IBM forcing, or unapproved production IBM forcing; or
- modify pressure, projection, Poisson, RK3, or channel-forcing logic.

## Defaults and wrapper interface

`stage18_checks/run_stage18_1_physical_structure_config.sh` infers the repository root from its own script path.  `DECOMP2D_ROOT` is retained only as an interface compatibility variable and is not used as the repository root.  The wrapper creates `stage18_outputs/`, invokes only the Stage 18.1 helper, builds nothing, runs no MPI, and runs no production validation.

Safe defaults are:

```text
STAGE18_1_ENABLE=1
STAGE18_1_STRUCTURE_DYNAMICS_CONFIG_ENABLE=1
STAGE18_1_SINGLE_FIBRE_ONLY=1
STAGE18_1_PHYSICAL_STRUCTURE_BOUNDARY=1
STAGE18_1_BENDING_CONFIG_ENABLE=1
STAGE18_1_TENSION_CONFIG_ENABLE=1
STAGE18_1_INEXTENSIBILITY_CONFIG_ENABLE=1
STAGE18_1_TIME_INTEGRATION_CONFIG_ENABLE=1
STAGE18_1_ENERGY_DIAGNOSTIC_CONFIG_ENABLE=1
STAGE18_1_DIAGNOSTIC_ONLY=1
STAGE18_1_RHO_S=1.0
STAGE18_1_FIBRE_RADIUS=1.0e-3
STAGE18_1_YOUNG_MODULUS=1.0
STAGE18_1_FIBRE_LENGTH=1.0
STAGE18_1_NPTS=8
STAGE18_1_DT_STRUCTURE=1.0e-4
STAGE18_1_RHO_TILDE=1.0
STAGE18_1_GAMMA=1.0e-3
STAGE18_1_ZERO_TOL=1.0e-14
STAGE18_1_FORMULA_TOL=1.0e-12
```

## False-positive-safe audit policy

Stage 18.1 reuses the corrected Stage 18.0, Stage 17.11, Stage 17.10, Stage 17.6, and Stage 16 false-positive-safe audit pattern without editing closed files:

- use targeted structural checks rather than broad repository-wide scans;
- do not scan Markdown documentation as real code-regression evidence;
- do not treat documentation text, negative-check strings, or old failure labels as real regressions;
- do not treat regex literals such as `rg[[:space:]]` as actual `rg` command usage;
- do not require `rg`; any future shell use of `rg` must include a `grep` fallback;
- do not classify a source-only archive without `.git` metadata as DNS-core contamination or closed-stage modification;
- do not treat Stage 2 or legacy structure files as Stage 18.1 activation evidence;
- accept `VALUE_SUFFIXES`, `VALUE_KEYS`, or `pass_fail_keys` style exclusion logic for non-boolean numeric/string fields; and
- do not treat `final_status` as missing before it is written.

## Expected PASS verdict

With Stage 18.0 evidence intact, Stage 17 closure evidence intact or structurally accepted, the corrected Stage 18.0 wrapper root fix preserved, and default valid physical parameters, the wrapper should print:

```text
STAGE 18.1 PHYSICAL STRUCTURE CONFIG VERDICT: PASS
STAGE 18.1 FINAL VERDICT: PASS
```

The helper writes `stage18_outputs/fibre_stage18_1_physical_structure_config.dat` with all required Stage 18.1 summary fields and `final_status PASS` on the passing path.
