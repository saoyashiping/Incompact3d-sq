# P0_13 Fix Note: synthetic zero-force proxy no-response

## Problem

P0_13 static audits, CMake configuration, target builds, and P0_2-P0_11/P0_13 closure micro-checks passed, but the validation failed when running `fibre_prod_synthetic_closed_loop_check` under np=1/2/4.

The failing diagnostic was:

```text
P0_12_SYNTHETIC_CLOSED_LOOP_CHECK FAIL: zero-force proxy run failed
```

## Root cause

`src/fibre_prod_synthetic_closed_loop.f90` still contained the older P0_12 small-lambda assertion:

```fortran
if (lambda_fsi > 0.0_dp .and. maxval(abs(rhs_increment_x)) <= 0.0_dp) status = 41
```

This is correct for a nonzero force buffer, but it is wrong for the deliberate zero-force proxy case used by the closure check:

```text
beta_hydro = 0
  -> hydro_force_candidate = 0
  -> structure_input_force = 0
  -> reaction_force_candidate = 0
  -> force_buffer = 0
  -> lambda_fsi > 0
  -> RHS increment must remain 0
```

The old assertion incorrectly treated this valid no-response proof as a failure.

## Fix

The nonzero RHS increment requirement is now gated by a nonzero force-buffer component:

```fortran
if (status == 0 .and. lambda_fsi > 0.0_dp .and. &
    maxval(abs(force_buffer%fx)) > 0.0_dp .and. &
    maxval(abs(rhs_increment_x)) <= 0.0_dp) status = 41
```

The exact force-buffer-to-RHS relation remains enforced first:

```fortran
rhs_increment_x == lambda_fsi * penalty_beta * force_buffer%fx
```

## Scope

Only the synthetic closed-loop validation predicate was corrected.  No pressure/projection/RK3/channel-forcing logic was modified.  No production DNS path was enabled.
