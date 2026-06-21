# P0.1 Source Diff Summary

## Modified files

- `src/CMakeLists.txt`
- `src/fibre_prod_rhs_adapter.f90`
- `src/fibre_prod_main_hook.f90`
- `src/fibre_prod_main_hook_check.f90`

## Added evidence files

- `production_recovery/P0_1_PLAN.md`
- `production_recovery/P0_1_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/P0_1_PASS_FAIL.md`
- `production_recovery/P0_1_evidence/P0_1_VALIDATION_COMMAND.sh`

## Main changes

### 1. Main target module integration

`xcompact3d` now includes the full production fibre helper chain in its source list, including production state, grid adapter, IBM delta/interpolation/spreading, force buffer, structure solver, FSI coupling, wall contact, and fibre collision helper modules.

### 2. RHS adapter safety gate

The previous uniform scalar perturbation was removed:

```fortran
contribution = config%lambda_fsi * config%penalty_beta * config%dt
rhs_x(i, j, k) = rhs_x(i, j, k) + contribution
```

The new nonzero-lambda path requires explicit finite force-density buffers:

```fortran
call fibre_prod_rhs_adapter_apply(..., force_x, force_y, force_z)
```

If lambda is nonzero but no force-density buffer is supplied, the adapter returns status `13` and leaves the RHS unchanged.

### 3. Lambda-zero safety

Disabled mode and `lambda_fsi = 0` remain strict no-ops.

### 4. Main hook API

`fibre_prod_main_hook_apply` keeps the old call signature compatible with `xcompact3d`, and adds optional force-density buffers for the future physical force path.

## Remaining blocker after P0.1

P0.1 does not yet make `xcompact3d` produce the runtime Eulerian force-density buffer from Lagrangian fibre dynamics. That is the required P0.2/P0.3 work before production DNS-FSI runs.
