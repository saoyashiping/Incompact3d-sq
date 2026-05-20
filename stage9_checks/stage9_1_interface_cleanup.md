# Stage 9.1 real 2decomp-fft interface cleanup

## Scope

This stage is limited to real 2decomp-fft interface takeover and compatibility cleanup:
- remove local compat/stub layers;
- route I/O interfaces to `decomp_2d_io`;
- route coarse conversion/bounds and `update_halo` to `decomp_2d`;
- restore compile-time dependency consistency.

No changes are introduced to physical models, RHS formulations, projection behavior, IBM numerical coupling logic, or Stage 7/8 verification logic.

## Stage 9.1 outcomes

- Compat layer removed:
  - `src/xcompact3d_decomp_io_compat.f90` deleted.
  - all active `use xcompact3d_decomp_io_compat` removed from source.
  - `src/CMakeLists.txt` no longer lists compat source entries.

- Stub removed:
  - `src/fibre_stage_decomp2d_constants_stub.f90` deleted.
  - no source/build references remain to this local constants stub.

- I/O imports unified to `decomp_2d_io`:
  - `decomp_2d_init_io`, `decomp_2d_register_variable`, `decomp_2d_open_io`,
    `decomp_2d_close_io`, `decomp_2d_start_io`, `decomp_2d_end_io`,
    `decomp_2d_write_one`, `decomp_2d_read_one`, `decomp_2d_write_plane`,
    `decomp_2d_write_mode`, `decomp_2d_read_mode`, `decomp_2d_append_mode`,
    `gen_iodir_name` imported from real `decomp_2d_io`.

- Coarse conversion/bounds unified to `decomp_2d`:
  - `fine_to_coarseS`/`fine_to_coarseV` imported from `decomp_2d`.
  - local coarse-bounds declarations and initializer removed from `variables.f90`.
  - coarse-bounds call path removed from `xcompact3d.f90` and `variables.f90`.

- `update_halo` source restored:
  - `use m_halo` removed.
  - `update_halo` imported from `decomp_2d` in affected files.

## Static checks

- Added: `stage9_checks/run_stage9_1_interface_consistency.sh`
- This script checks for:
  - compat/stub file removal;
  - no active compat/m_halo import residues;
  - no local coarse-bounds initializer residues;
  - active usage of real `decomp_2d_io` and `decomp_2d` symbols;
  - decomp2d target links retained in CMake.
