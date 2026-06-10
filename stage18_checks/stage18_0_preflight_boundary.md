# Stage 18.0: Stage 17 closure and Stage 18 preflight boundary

Stage 18.0 is a diagnostic-only preflight layer.  It verifies that Stage 17 is closed, or that fresh-archive Stage 17.11 structural closure evidence can be safely accepted when `STAGE18_0_ACCEPT_STAGE17_CLOSED_EVIDENCE=1`, and it declares the strict Stage 18 boundary.

## Stage 18 boundary

Stage 18 is limited to **single-fibre physical structure dynamics enhancement**.  Later Stage 18 substages may introduce and validate physical single-fibre structure dynamics, but Stage 18.0 does not implement any of that physics.

Stage 18 later targets the single-fibre structure equation:

```text
rho_l X_tt = d_s(T X_s) - B X_ssss + F_h
```

or nondimensional:

```text
rho_tilde X_tt = d_s(T X_s) - gamma X_ssss + F_h
```

Stage 18 later targets the inextensibility constraint:

```text
X_s dot X_s = 1
```

Stage 18 later targets structure energy diagnostics:

```text
E_k = 1/2 int rho_l |V|^2 ds
E_b = 1/2 int B |X_ss|^2 ds
P_h = int F_h dot V ds
```

Stage 18.0 itself must not implement these equations yet; it only declares the boundary.

## Stage 18.0 non-implementation guarantees

The Stage 18.0 wrapper and helper:

- add no Fortran source;
- modify no CMake file;
- build no target;
- run no MPI command;
- run no production physics or production validation;
- introduce no bending operator activation;
- introduce no tension solve activation;
- introduce no inextensibility activation;
- introduce no structure time-integration core;
- introduce no structure state update;
- introduce no fluid-on-fibre force integration into the physical structure equation;
- introduce no structure energy/power diagnostic implementation;
- introduce no real wall-contact, fibre-fibre collision, penalty, repulsive, lubrication, friction, adhesion, or contact-damping force;
- introduce no production multi-fibre logic;
- introduce no collision-induced RHS or collision-induced structure update;
- introduce no direct RHS injection, unapproved Stage 14 RHS call, legacy IBM forcing, or unapproved production IBM forcing; and
- modify no pressure, projection, Poisson, RK3, or channel-forcing logic.

## Closed-stage preservation policy

Stage 18.0 does not modify closed Stage 10--17 files.  Closed files are inspected only as read-only evidence.  The helper preserves the corrected false-positive-safe audit policy from Stage 17.11, Stage 17.10, Stage 17.6, and Stage 16:

- use targeted structural checks rather than broad repository-wide scans;
- do not scan Markdown documentation as real code-regression evidence;
- do not treat documentation text, negative-check strings, or old failure labels as real regressions;
- do not treat regex literals such as `rg[[:space:]]` as actual `rg` command usage;
- do not require `rg`; any future shell use of `rg` must include a `grep` fallback;
- do not classify a source-only archive without `.git` metadata as DNS-core contamination or closed-stage modification;
- accept `VALUE_SUFFIXES`, `VALUE_KEYS`, or `pass_fail_keys` style exclusion logic for non-boolean numeric/string fields; and
- do not treat `final_status` as missing before it is written.

## Wrapper interface

`stage18_checks/run_stage18_0_preflight_boundary.sh` supports these environment variables with safe defaults:

- `DECOMP2D_ROOT=$(pwd)`
- `BUILD_DIR=build_stage9`
- `MPIEXEC=mpirun`
- `MPIEXEC_FLAGS=`
- `STAGE18_0_REQUIRE_STAGE17_CLOSED=1`
- `STAGE18_0_ACCEPT_STAGE17_CLOSED_EVIDENCE=1`
- `STAGE18_0_ENABLE=1`
- `STAGE18_0_SINGLE_FIBRE_STRUCTURE_DYNAMICS_BOUNDARY=1`
- `STAGE18_0_DIAGNOSTIC_ONLY=1`

The wrapper creates `stage18_outputs/`, invokes only the Stage 18.0 helper, and writes `stage18_outputs/fibre_stage18_0_preflight_boundary.dat`.

## Expected PASS verdict

With Stage 17.11 closure evidence intact and no unapproved closed-stage modifications, the wrapper should print:

```text
STAGE 18.0 PREFLIGHT BOUNDARY VERDICT: PASS
STAGE 18.0 FINAL VERDICT: PASS
```

The summary output contains all Stage 18.0 required status fields, including `stage17_closed_file_status`, `stage17_closed_evidence_status`, all no-contamination statuses, and `final_status`.
