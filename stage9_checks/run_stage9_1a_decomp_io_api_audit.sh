#!/usr/bin/env bash
set -euo pipefail
: "${DECOMP2D_ROOT:=/home/sq/opt/2decomp-fft}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p stage9_outputs

OUT="stage9_outputs/fibre_stage9_1a_decomp_io_api_audit.dat"
REPORT="stage9_outputs/STAGE9_1A_DECOMP_IO_API_AUDIT.md"

: > "$OUT"
: > "$REPORT"

emit(){ echo "$1=$2" >> "$OUT"; }

get_key(){
  local f="$1" k="$2"
  awk -F= -v key="$k" '$1==key{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$f" 2>/dev/null || true
}

dep_file="stage9_outputs/fibre_stage9_dependency_gate_check.dat"
if [[ -f "$dep_file" ]]; then emit stage9_1a_stage9_0_output_exists 1; else emit stage9_1a_stage9_0_output_exists 0; fi
s90=$(get_key "$dep_file" stage9_dependency_gate_check_status); s90=${s90:-0}
a91=$(get_key "$dep_file" stage9_gate_stage9_1_allowed_flag); a91=${a91:-0}
pca=$(get_key "$dep_file" stage9_gate_production_coupling_allowed_flag); pca=${pca:-0}
emit stage9_1a_stage9_0_status "$s90"
emit stage9_1a_stage9_1_allowed_flag "$a91"
emit stage9_1a_production_coupling_allowed_flag "$pca"
if [[ "$s90" == "1" && "$a91" == "1" && "$pca" == "0" ]]; then emit stage9_1a_stage9_0_dependency_status 1; else emit stage9_1a_stage9_0_dependency_status 0; fi

[[ -d "$DECOMP2D_ROOT" ]] && emit stage9_1a_decomp2d_root_exists 1 || emit stage9_1a_decomp2d_root_exists 0
if [[ -d "$DECOMP2D_ROOT" ]] && [[ -n "$(find "$DECOMP2D_ROOT" -mindepth 1 -maxdepth 2 2>/dev/null | head -n1)" ]]; then emit stage9_1a_decomp2d_root_nonempty_flag 1; else emit stage9_1a_decomp2d_root_nonempty_flag 0; fi
emit stage9_1a_real_decomp2d_required_flag 1
emit stage9_1a_no_decomp2d_stub_policy_flag 1
if [[ "$(get_key "$OUT" stage9_1a_decomp2d_root_exists)" == "1" && "$(get_key "$OUT" stage9_1a_decomp2d_root_nonempty_flag)" == "1" ]]; then emit stage9_1a_decomp2d_presence_status 1; else emit stage9_1a_decomp2d_presence_status 0; fi

[[ -f src/xcompact3d_decomp_io_compat.f90 ]] && emit stage9_1a_compat_module_exists 1 || emit stage9_1a_compat_module_exists 0
compat_line=$(grep -n "xcompact3d_decomp_io_compat.f90" src/CMakeLists.txt | head -n1 | cut -d: -f1 || true)
[[ -n "$compat_line" ]] && emit stage9_1a_compat_in_cmake_flag 1 || emit stage9_1a_compat_in_cmake_flag 0
check_before(){ local tline; tline=$(grep -n "$1" src/CMakeLists.txt | head -n1 | cut -d: -f1 || true); [[ -n "$compat_line" && -n "$tline" && "$compat_line" -lt "$tline" ]] && echo 1 || echo 0; }
emit stage9_1a_compat_before_forces_flag "$(check_before 'forces.f90')"
emit stage9_1a_compat_before_statistics_flag "$(check_before 'statistics.f90')"
emit stage9_1a_compat_before_visu_flag "$(check_before 'visu.f90')"
emit stage9_1a_compat_before_tools_flag "$(check_before 'tools.f90')"
emit stage9_1a_compat_before_les_models_flag "$(check_before 'les_models.f90')"
case_line=$(grep -n "Case-" src/CMakeLists.txt | head -n1 | cut -d: -f1 || true)
[[ -n "$compat_line" && -n "$case_line" && "$compat_line" -lt "$case_line" ]] && emit stage9_1a_compat_before_case_files_flag 1 || emit stage9_1a_compat_before_case_files_flag 0
if [[ "$(get_key "$OUT" stage9_1a_compat_module_exists)" == "1" && "$(get_key "$OUT" stage9_1a_compat_in_cmake_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_forces_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_statistics_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_visu_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_tools_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_les_models_flag)" == "1" && "$(get_key "$OUT" stage9_1a_compat_before_case_files_flag)" == "1" ]]; then emit stage9_1a_compat_cmake_order_status 1; else emit stage9_1a_compat_cmake_order_status 0; fi

forbidden_pat='use[[:space:]]+decomp_2d_io[[:space:]]*,[[:space:]]*only[[:space:]]*:[[:space:]]*(decomp_2d_init_io|decomp_2d_register_variable|decomp_2d_write_one|decomp_2d_read_one|decomp_2d_write_plane|fine_to_coarseS|fine_to_coarseV|gen_iodir_name|decomp_2d_open_io|decomp_2d_close_io|decomp_2d_start_io|decomp_2d_end_io|decomp_2d_write_mode|decomp_2d_read_mode|decomp_2d_append_mode)'
forbidden_count=$( (grep -RInE "$forbidden_pat" src || true) | wc -l | tr -d ' ')
emit stage9_1a_forbidden_direct_old_io_import_count "$forbidden_count"
[[ "$forbidden_count" == "0" ]] && emit stage9_1a_forbidden_direct_old_io_imports_absent_flag 1 || emit stage9_1a_forbidden_direct_old_io_imports_absent_flag 0
xio_ok=0
if grep -q "use decomp_2d_io, only : decomp_2d_io_init" src/xcompact3d.f90 && grep -q "use decomp_2d_io, only : decomp_2d_io_finalise" src/xcompact3d.f90; then xio_ok=1; fi
emit stage9_1a_xcompact3d_real_io_init_finalise_allowed_flag "$xio_ok"
if [[ "$(get_key "$OUT" stage9_1a_forbidden_direct_old_io_imports_absent_flag)" == "1" && "$xio_ok" == "1" ]]; then emit stage9_1a_direct_import_policy_status 1; else emit stage9_1a_direct_import_policy_status 0; fi

compat=src/xcompact3d_decomp_io_compat.f90
has_all(){ local ok=1; for p in "$@"; do grep -q "$p" "$compat" || ok=0; done; echo "$ok"; }
emit stage9_1a_compat_metadata_wrappers_status "$(has_all 'public :: decomp_2d_init_io' 'public :: decomp_2d_register_variable' 'public :: decomp_2d_open_io' 'public :: decomp_2d_close_io' 'public :: decomp_2d_start_io' 'public :: decomp_2d_end_io' 'subroutine decomp_2d_init_io' 'subroutine decomp_2d_register_variable' 'subroutine decomp_2d_open_io' 'subroutine decomp_2d_close_io' 'subroutine decomp_2d_start_io' 'subroutine decomp_2d_end_io')"
emit stage9_1a_compat_write_one_wrapper_status "$(has_all 'public :: decomp_2d_write_one' 'interface decomp_2d_write_one' 'x3d_write_one_r3_simple' 'x3d_write_one_r3_legacy')"
emit stage9_1a_compat_read_one_wrapper_status "$(has_all 'public :: decomp_2d_read_one' 'interface decomp_2d_read_one' 'x3d_read_one_r3_simple' 'x3d_read_one_r3_legacy')"
emit stage9_1a_compat_write_plane_wrapper_status "$(has_all 'public :: decomp_2d_write_plane' 'interface decomp_2d_write_plane' 'x3d_write_plane_r3_simple' 'x3d_write_plane_r3_legacy')"
emit stage9_1a_compat_fine_to_coarse_wrapper_status "$(has_all 'public :: fine_to_coarseS' 'public :: fine_to_coarseV' 'interface fine_to_coarseS' 'interface fine_to_coarseV' 'fine_to_coarseS_r3' 'fine_to_coarseV_r3')"
emit stage9_1a_compat_gen_iodir_name_status "$(has_all 'public :: gen_iodir_name' 'function gen_iodir_name')"
if [[ "$(get_key "$OUT" stage9_1a_compat_metadata_wrappers_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_write_one_wrapper_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_read_one_wrapper_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_write_plane_wrapper_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_fine_to_coarse_wrapper_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_gen_iodir_name_status)" == "1" ]]; then emit stage9_1a_compat_coverage_status 1; else emit stage9_1a_compat_coverage_status 0; fi

legacy_calls=$(grep -RInE "call[[:space:]]+(decomp_2d_write_one|decomp_2d_read_one|decomp_2d_write_plane|fine_to_coarseS|fine_to_coarseV|decomp_2d_register_variable|decomp_2d_init_io|decomp_2d_open_io|decomp_2d_close_io|decomp_2d_start_io|decomp_2d_end_io)" src | cut -d: -f1 | sort -u || true)
legacy_count=$(echo "$legacy_calls" | sed '/^$/d' | wc -l | tr -d ' ')
uncovered=0
coverage_table=""
while IFS= read -r f; do
  [[ -z "$f" ]] && continue
  if grep -qi "use[[:space:]]\+xcompact3d_decomp_io_compat" "$f"; then cov=YES; else cov=NO; uncovered=$((uncovered+1)); fi
  coverage_table+="| $f | $cov |\n"
done <<< "$legacy_calls"
emit stage9_1a_legacy_io_call_file_count "$legacy_count"
emit stage9_1a_legacy_io_uncovered_file_count "$uncovered"
[[ "$legacy_count" -ge 1 && "$uncovered" -eq 0 ]] && emit stage9_1a_legacy_io_call_coverage_flag 1 || emit stage9_1a_legacy_io_call_coverage_flag 0
[[ "$(get_key "$OUT" stage9_1a_legacy_io_call_coverage_flag)" == "1" ]] && emit stage9_1a_legacy_call_coverage_status 1 || emit stage9_1a_legacy_call_coverage_status 0

grep -q "Stage 9.1 build-only compatibility path" "$compat" && emit stage9_1a_build_only_wrapper_marker_flag 1 || emit stage9_1a_build_only_wrapper_marker_flag 0
grep -q "Stage 9.10" "$compat" && grep -q "restart/output semantics" "$compat" && emit stage9_1a_stage9_10_output_semantics_marker_flag 1 || emit stage9_1a_stage9_10_output_semantics_marker_flag 0
[[ "$(get_key "$OUT" stage9_1a_build_only_wrapper_marker_flag)" == "1" && "$(get_key "$OUT" stage9_1a_stage9_10_output_semantics_marker_flag)" == "1" ]] && emit stage9_1a_wrapper_semantics_warning_status 1 || emit stage9_1a_wrapper_semantics_warning_status 0

if [[ "$(get_key "$OUT" stage9_1a_stage9_0_dependency_status)" == "1" && "$(get_key "$OUT" stage9_1a_decomp2d_presence_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_cmake_order_status)" == "1" && "$(get_key "$OUT" stage9_1a_direct_import_policy_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_coverage_status)" == "1" && "$(get_key "$OUT" stage9_1a_legacy_call_coverage_status)" == "1" && "$(get_key "$OUT" stage9_1a_wrapper_semantics_warning_status)" == "1" ]]; then
  emit stage9_1a_stage9_1b_allowed_flag 1
  emit stage9_1a_stage9_2_allowed_flag 0
  emit stage9_1a_production_coupling_allowed_flag 0
  emit stage9_1a_next_step_policy_status 1
else
  emit stage9_1a_stage9_1b_allowed_flag 0
  emit stage9_1a_stage9_2_allowed_flag 0
  emit stage9_1a_production_coupling_allowed_flag 0
  emit stage9_1a_next_step_policy_status 0
fi

if [[ "$(get_key "$OUT" stage9_1a_stage9_0_dependency_status)" == "1" && "$(get_key "$OUT" stage9_1a_decomp2d_presence_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_cmake_order_status)" == "1" && "$(get_key "$OUT" stage9_1a_direct_import_policy_status)" == "1" && "$(get_key "$OUT" stage9_1a_compat_coverage_status)" == "1" && "$(get_key "$OUT" stage9_1a_legacy_call_coverage_status)" == "1" && "$(get_key "$OUT" stage9_1a_wrapper_semantics_warning_status)" == "1" && "$(get_key "$OUT" stage9_1a_next_step_policy_status)" == "1" ]]; then emit stage9_1a_decomp_io_api_audit_status 1; else emit stage9_1a_decomp_io_api_audit_status 0; fi

allowed_direct=$(grep -RIn "use decomp_2d_io" src || true)
forbidden_direct=$(grep -RInE "$forbidden_pat" src || true)

cat > "$REPORT" <<MD
# Stage 9.1a Decomp IO API Compatibility Audit

## 1. Purpose
Static API compatibility audit for Xcompact3D production I/O usage vs installed decomp2d/decomp_2d_io prior to Stage 9.1b build regression.

## 2. DECOMP2D_ROOT

- audit value: $DECOMP2D_ROOT

## 3. Allowed remaining direct decomp_2d_io uses

\`\`\`
$allowed_direct
\`\`\`

## 4. Forbidden old direct imports count

- count: $(get_key "$OUT" stage9_1a_forbidden_direct_old_io_import_count)

\`\`\`
$forbidden_direct
\`\`\`

## 5. Compat wrapper coverage table

| Check | Flag |
|---|---|
| metadata wrappers | $(get_key "$OUT" stage9_1a_compat_metadata_wrappers_status) |
| write_one wrappers | $(get_key "$OUT" stage9_1a_compat_write_one_wrapper_status) |
| read_one wrappers | $(get_key "$OUT" stage9_1a_compat_read_one_wrapper_status) |
| write_plane wrappers | $(get_key "$OUT" stage9_1a_compat_write_plane_wrapper_status) |
| fine_to_coarse wrappers | $(get_key "$OUT" stage9_1a_compat_fine_to_coarse_wrapper_status) |
| gen_iodir_name | $(get_key "$OUT" stage9_1a_compat_gen_iodir_name_status) |

## 6. Legacy call coverage table

| File | compat imported |
|---|---|
$coverage_table

## 7. Build-only wrapper warning

Compat module markers indicate this is a build-only compatibility path for Stage 9.1.

## 8. Stage 9.10 requirement

Stage 9.10 restart/output semantics must audit real I/O behavior.

## 9. Next allowed step

Stage 9.1b only (Stage 9.2 disallowed).
MD
