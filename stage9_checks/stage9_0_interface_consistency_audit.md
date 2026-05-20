# Stage 9.0 interface consistency audit (pre-Stage 9.1)

## 1) Scope

This is an **audit-only** report comparing the current branch against the original Xcompact3D/master interface expectations to detect interface drift before Stage 9.1.

The audit focuses on:
- original Xcompact3D public procedure signatures;
- module dependency/source-of-symbol changes;
- I/O interface consistency (`decomp_2d_io` / `decomp_2d`);
- decomposition-bound interface consistency (`xstS/xenS/xstV/xenV`, `ph*`).

No source behavior changes are made in this PR.

---

## 2) Original Xcompact3D files modified

Based on current branch delta against the immediate upstream baseline available in local history (`HEAD` change set), the following original files differ.

### Case files
- `src/Case-PTBL.f90`

### Core solver files
- `src/acl_source.f90`
- `src/mhd.f90`
- `src/navier.f90`
- `src/particle.f90`
- `src/statistics.f90`

### I/O and visualization files
- None in this change set.

### Variable/decomposition files
- `src/variables.f90`

### Build-system files
- None in this change set.

---

## 3) Procedure signature audit

Checked modified original files for changes to existing subroutine/function argument lists.

### Result
- **No existing original Xcompact3D subroutine/function signatures were changed** in the modified original files listed above.
- Changes are primarily in `use` dependencies and internal call-site/module-source adjustments.

### New fibre-only procedures/programs
- None added in the listed modified original files for this audit scope.

---

## 4) Module dependency audit

This section records dependency-source drift with interface risk potential.

### 4.1 `update_halo` import source changed

1. `src/acl_source.f90`
   - Original dependency (baseline): `use decomp_2d, only: ..., update_halo`
   - Current dependency: `use decomp_2d, only: ...` + `use m_halo, only: update_halo`
   - Risk: module-source drift can break compatibility if `m_halo::update_halo` interface/version diverges from original expected provider.

2. `src/particle.f90`
   - Original dependency (baseline): `use decomp_2d, only: ..., update_halo`
   - Current dependency: `use decomp_2d, only: ...` + `use m_halo, only: update_halo`
   - Risk: same as above; interface-provider substitution in original files.

### 4.2 Decomposition metadata imported from `var`

Multiple original files now import `ph1/ph2/ph3/ph4` via `use var` additions:
- `src/Case-PTBL.f90`
- `src/mhd.f90`
- `src/navier.f90`

Original expectation in master for decomposition objects is typically native `decomp_2d` ownership path; routing through `var` can create interface coupling to project-local declarations.

Risk:
- if `var` locally redefines or shadows native decomposition objects, interface consistency with upstream `decomp_2d` contracts can drift.

### 4.3 `decomp_2d_io` vs `decomp_2d` source shift for coarse-bound/coarsening symbols

`src/statistics.f90` changed imports:
- baseline used `decomp_2d_io` for `xstS/xenS` and `fine_to_coarseS` in affected scopes;
- current uses `decomp_2d` for these symbols.

Risk:
- if the source-matched library exposes these symbols in `decomp_2d`, this may be acceptable;
- however it is a module-source drift from baseline and should be validated as an intentional upstream-compatible change.

### 4.4 `decomp_2d_io -> xcompact3d_decomp_io_compat`

Audit result in current branch scope:
- no active replacement detected in currently modified original files.

---

## 5) Decomposition-bound audit

Target symbols:
- `xstS`, `xenS`, `xstV`, `xenV`
- `xcompact3d_coarse_bounds_initialized`
- `init_xcompact3d_coarse_bounds`

### Findings in current branch scope
- `xcompact3d_coarse_bounds_initialized`: not introduced in current modified set.
- `init_xcompact3d_coarse_bounds`: not introduced in current modified set.
- `xstS/xenS/xstV/xenV`: used/imported; this audit did **not** find current-branch additions of local coarse-bound state in the changed files list for this PR scope.
- `src/variables.f90` includes local `type(DECOMP_INFO), save :: ph1, ph2, ph3, ph4, phG` in the current file content, which may duplicate native ownership expectations depending on upstream master state.

Conclusion:
- no evidence in this PR scope of new local coarse-bound init hook;
- possible decomposition object ownership duplication risk remains via `var`-level `ph*` declarations and downstream imports.

---

## 6) I/O compatibility-layer audit

Searched for:
- `xcompact3d_decomp_io_compat`
- `xcompact3d_decomp_io_compat.f90`

### Findings in this branch scope
- No active use/remapping to `xcompact3d_decomp_io_compat` found in the currently changed files for this audit.
- `src/CMakeLists.txt` in this audit scope does not list compat module additions in this change set.

### Interface interpretation
- No direct evidence in this specific delta of replacing original `decomp_2d_io` with compat-layer imports.
- If compat layer exists in other branches/history, it is out-of-scope for this change set and should be re-audited at merge-time against `origin/main`.

---

## 7) Fibre-module dependency audit

Added/modified fibre-prefixed file in current branch scope:
- `src/fibre_stage9_build_regression_check.f90`

Dependency characteristics:
- reads text evidence files only;
- does not call solver/DNS/projection/RHS/fluid/fibre advance APIs;
- does not depend on `xcompact3d_decomp_io_compat`;
- does not depend on duplicated coarse-bound variables in `var`.

Risk:
- low interface risk to original solver interfaces (standalone audit executable behavior).

---

## 8) Risk classification

### A. Must fix before Stage 9.1
1. **Module-source drift for `update_halo` in original files** (`decomp_2d` -> `m_halo`) without explicit upstream interface-equivalence proof.
2. **Decomposition metadata source drift via `use var, only: ph*`** in original files (`Case-PTBL`, `mhd`, `navier`) that may tie original interfaces to project-local ownership.
3. **Potential `decomp_2d_io` -> `decomp_2d` source drift** for `xstS/xenS/fine_to_coarseS` in `statistics.f90` must be confirmed against intended upstream interface contract.

### B. Can defer
1. Standalone fibre audit executable (`fibre_stage9_build_regression_check`) that does not alter solver/public behavior.
2. Documentation/diagnostic/comment-level changes with no interface implications.

### C. Needs manual decision
1. Whether `decomp_2d` should now be canonical source for `xstS/xenS/fine_to_coarseS` under the chosen source-matched 2decomp package.
2. Whether `ph*` should remain `var`-sourced in modified original files or be restored to native decomp ownership/import pattern used by master.

---

## 9) Recommended next PR

**Recommended PR name:**

> Stage 9.0A restore native Xcompact3D interfaces

**Recommended scope (do not implement in this audit PR):**
- restore original `decomp_2d_io` imports in original Xcompact3D files where master expects them;
- remove dependency on compatibility-layer indirection from original files;
- avoid duplicated local decomposition bounds/state in `variables.f90`;
- keep fibre additions isolated from original solver interfaces;
- preserve original Xcompact3D public subroutine/function signatures.

