# R1 Baseline Run Plan

## Candidate baseline case/input

Preferred candidate if build prerequisites become available: `tests/TGV-Taylor-Green-vortex/reference_input.i3d`, because it is an existing test-oriented baseline input and does not represent a fibre/FSI production validation case.

## R1 input copy

No R1 input copy was created in this attempt because configure failed before any executable was built. Without a built executable and MPI launcher, no runnable baseline case could be exercised.

## Fibre/FSI disabled assurance

R1 would require the selected input copy to keep fibre/FSI disabled or absent. In this blocked attempt, no runtime confirmation was possible.

## Short, low-cost sanity-run intent

If build and MPI tooling become available, the selected R1 input copy should be reduced only in the R1 evidence directory to a short-step, small-cost sanity run. Original input files must not be modified.

## Intended np=1/2/4 commands

```sh
mpirun -np 1 <path-to-built-xcompact3d> <R1-input>
mpirun -np 2 <path-to-built-xcompact3d> <R1-input>
mpirun -np 4 <path-to-built-xcompact3d> <R1-input>
```

## Actual run-plan status

BLOCKED: `cmake` configure failed because no Fortran compiler was found, and `mpirun` is not available in the environment.
