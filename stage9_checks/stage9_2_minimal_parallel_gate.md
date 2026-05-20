# Stage 9.2 minimal parallel gate after real 2decomp interface cleanup

## Goal

Stage 9.2 establishes a minimal MPI/decomposition gate after Stage 9.1 real interface cleanup. It validates that:
- `xcompact3d` builds;
- `fibre_stage9_dependency_gate_check` builds;
- Stage 9.1 static interface gate still passes;
- a new minimal MPI decomposition program passes for `np=1/2/4`.

## Method and checks

The Stage 9.2 program `fibre_stage9_2_minimal_parallel_gate` uses real `decomp_2d`, `decomp_2d_mpi`, and `decomp_2d_constants` and performs:
- MPI init/finalize and 2D pencil decomposition init/finalize;
- per-rank reporting of local decomposition metadata (`xsize/ysize/zsize`, start/end indices);
- local consistency checks for non-negative sizes and valid start/end/index relations;
- coarse/stat/visu bounds validity checks (`xstS/xenS/xstV/xenV` non-negative and non-inverted);
- MPI reductions for total local cell count and fail-flag aggregation.

Rank 0 prints final single-run verdict:
- `STAGE 9.2 PROGRAM VERDICT: PASS` or
- `STAGE 9.2 PROGRAM VERDICT: FAIL`.

## Full gate pass criteria

`stage9_checks/run_stage9_2_minimal_parallel_gate.sh` passes only if all succeed:
1. `cmake --build build_stage9 --target xcompact3d -j`
2. `cmake --build build_stage9 --target fibre_stage9_dependency_gate_check -j`
3. `cmake --build build_stage9 --target fibre_stage9_2_minimal_parallel_gate -j`
4. `bash stage9_checks/run_stage9_1_interface_consistency.sh`
5. `${MPIEXEC:-mpirun} -np 1 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate`
6. `${MPIEXEC:-mpirun} -np 2 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate`
7. `${MPIEXEC:-mpirun} -np 4 ./build_stage9/bin/fibre_stage9_2_minimal_parallel_gate`

Script final verdict:
- `STAGE 9.2 FINAL VERDICT: PASS` or
- `STAGE 9.2 FINAL VERDICT: FAIL` (with failed steps listed).

## Manual run

```bash
cmake --build build_stage9 --target xcompact3d -j
cmake --build build_stage9 --target fibre_stage9_dependency_gate_check -j
cmake --build build_stage9 --target fibre_stage9_2_minimal_parallel_gate -j
MPIEXEC=mpirun bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
```

## Scope guard

Stage 9.2 does **not** change physics models, IBM numerical algorithms, RHS/projection/pressure solver behavior, or production two-way coupling logic. It only validates parallel infrastructure consistency after Stage 9.1 interface cleanup.
