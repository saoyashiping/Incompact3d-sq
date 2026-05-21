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

## Driver behavior (build/run control)

`stage9_checks/run_stage9_2_minimal_parallel_gate.sh` provides staged control:
1. configure stage;
2. build stage;
3. Stage 9.1 interface gate;
4. MPI run stage (`np=1/2/4`).

Key behavior:
- If build directory is missing, the script auto-configures with CMake.
- `BUILD_DIR` is configurable (default: `build_stage9`).
- `DECOMP2D_ROOT` (or `CMAKE_PREFIX_PATH` in environment) can be used to locate external 2decomp install.
- `MPIEXEC` is configurable (default: `mpirun`; can be `mpiexec`/`srun`).
- `MPIEXEC_FLAGS` is configurable for launcher flags (default empty).
- If configure fails: no build, no MPI runs.
- If any required build fails: no MPI runs.
- If Stage 9.1 interface gate fails: no MPI runs.
- Before MPI execution, script verifies stage9.2 executable exists and is executable; attempts `chmod +x` if needed.

Every run ends with final script verdict:
- `STAGE 9.2 FINAL VERDICT: PASS` or
- `STAGE 9.2 FINAL VERDICT: FAIL` (with failed checks listed).

## Full gate pass criteria

`stage9_checks/run_stage9_2_minimal_parallel_gate.sh` passes only if all succeed:
1. `cmake -S . -B ${BUILD_DIR}` (or with `-DCMAKE_PREFIX_PATH=${DECOMP2D_ROOT}` when set)
2. `cmake --build ${BUILD_DIR} --target xcompact3d -j`
3. `cmake --build ${BUILD_DIR} --target fibre_stage9_dependency_gate_check -j`
4. `cmake --build ${BUILD_DIR} --target fibre_stage9_2_minimal_parallel_gate -j`
5. `bash stage9_checks/run_stage9_1_interface_consistency.sh`
6. `${MPIEXEC:-mpirun} -np 1 ${BUILD_DIR}/bin/fibre_stage9_2_minimal_parallel_gate`
7. `${MPIEXEC:-mpirun} -np 2 ${BUILD_DIR}/bin/fibre_stage9_2_minimal_parallel_gate`
8. `${MPIEXEC:-mpirun} -np 4 ${BUILD_DIR}/bin/fibre_stage9_2_minimal_parallel_gate`

## Manual run

```bash
BUILD_DIR=build_stage9 MPIEXEC=mpirun bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
# optional external 2decomp path:
DECOMP2D_ROOT=/path/to/2decomp/install BUILD_DIR=build_stage9 MPIEXEC=mpiexec bash stage9_checks/run_stage9_2_minimal_parallel_gate.sh
```

## Scope guard

Stage 9.2 does **not** change physics models, IBM numerical algorithms, RHS/projection/pressure solver behavior, or production two-way coupling logic. It only validates parallel infrastructure consistency after Stage 9.1 interface cleanup.


## Stage 9.2 diagnosis and verdict details (latest)

- The Fortran gate now prints explicit per-rank `[PASS]/[FAIL]` sub-check lines for minimal MPI/decomposition consistency.
- Rank 0 prints aggregate failed sub-check count and program verdict.
- Minimal gate strong criteria are limited to MPI/decomposition consistency:
  - MPI init;
  - decomp init;
  - non-negative local sizes;
  - start/end consistency for non-empty local regions;
  - `size = end-start+1` for non-empty local regions;
  - global local-cell sum positive.
- Coarse/stat/visu bounds checks are treated as optional diagnostics unless real coarse initializers are explicitly called.
- 2decomp small-grid warnings and OpenMPI TCP reachable pairing warnings are not automatically treated as gate failures.
- Script supports `MPIEXEC_FLAGS` for launcher tuning, for example:
  - `MPIEXEC_FLAGS="--mca btl self,vader,tcp"`
  - `MPIEXEC_FLAGS="--mca btl self,vader"`
