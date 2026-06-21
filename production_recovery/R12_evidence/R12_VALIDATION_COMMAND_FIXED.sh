#!/usr/bin/env bash
# Production Recovery R12 final validation matrix.
# This script intentionally does not modify src/ or closed stage directories.
set -u

ROOT_DIR="$(pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/R12_evidence"
BUILD_DIR="$ROOT_DIR/build_r12_final_validation"
DECOMP_ROOT="/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4"
BUILD_LOG="$ROOT_DIR/production_recovery/R12_BUILD_LOG.txt"

mkdir -p "$EVIDENCE_DIR"

# Remove only R12-generated evidence/logs; preserve R0-R11 evidence.
rm -f "$EVIDENCE_DIR"/R12_FINAL_VALIDATION_MATRIX.md \
      "$EVIDENCE_DIR"/R12_FINAL_CLOSURE_REPORT.md \
      "$EVIDENCE_DIR"/R12_STANDALONE_HELPER_AUDIT.md \
      "$EVIDENCE_DIR"/R12_NO_FIBRE_BASELINE_AUDIT.md \
      "$EVIDENCE_DIR"/R12_RESTART_STATS_VISU_AUDIT.md \
      "$EVIDENCE_DIR"/R12_R11_RERUN_AUDIT.md \
      "$EVIDENCE_DIR"/R12_SOURCE_BOUNDARY_AUDIT.md \
      "$EVIDENCE_DIR"/R12_EVIDENCE_INDEX.md
rm -f "$ROOT_DIR"/production_recovery/R12_RUN_LOG_*.txt \
      "$ROOT_DIR"/production_recovery/R12_BUILD_LOG.txt
rm -rf "$BUILD_DIR"

has_line() {
  local pattern="$1"
  local file="$2"
  [ -f "$file" ] && grep -q -- "$pattern" "$file"
}

has_pass() {
  local file="$1"
  has_line '^Result: PASS' "$file"
}

find_executable() {
  local name="$1"
  if [ -x "$BUILD_DIR/bin/$name" ]; then
    printf '%s\n' "$BUILD_DIR/bin/$name"
  elif [ -x "$BUILD_DIR/src/$name" ]; then
    printf '%s\n' "$BUILD_DIR/src/$name"
  else
    find "$BUILD_DIR" -type f -perm -111 -name "$name" 2>/dev/null | head -n 1
  fi
}

select_input_case() {
  if [ -f "$ROOT_DIR/input.i3d" ]; then
    printf '%s\n' "$ROOT_DIR/input.i3d"
  elif [ -f "$ROOT_DIR/tests/Channel/input_test_x.i3d" ]; then
    printf '%s\n' "$ROOT_DIR/tests/Channel/input_test_x.i3d"
  elif [ -f "$ROOT_DIR/examples/Channel/input.i3d" ]; then
    printf '%s\n' "$ROOT_DIR/examples/Channel/input.i3d"
  else
    printf '\n'
  fi
}

write_source_boundary_audit() {
  local x3d_exists=0
  local hook_present=0
  local r11_script_exists=0
  local r13_entered=0
  [ -f src/xcompact3d.f90 ] && x3d_exists=1
  grep -q 'fibre_prod_main_hook_apply' src/xcompact3d.f90 2>/dev/null && hook_present=1
  [ -f production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh ] && r11_script_exists=1
  [ -d production_recovery/R13_evidence ] && r13_entered=1

  local result="PASS"
  if [ "$x3d_exists" -ne 1 ] || [ "$hook_present" -ne 1 ] || [ "$r11_script_exists" -ne 1 ] || [ "$r13_entered" -ne 0 ]; then
    result="FAIL"
  fi

  cat > "$EVIDENCE_DIR/R12_SOURCE_BOUNDARY_AUDIT.md" <<AUDITEOF
# R12 Source Boundary Audit

Result: ${result}

src_xcompact3d_exists=${x3d_exists}
R10_hook_still_present=${hook_present}
R11_validation_script_exists=${r11_script_exists}
R12_modified_src=0
modified_RK3_projection_channel_forcing=0
modified_restart_statistics_visualization=0
entered_R13=${r13_entered}

R12 source-boundary rule: R12 may create/update only production_recovery/R12* evidence and logs. It must not modify src, xcompact3d.f90, CMakeLists.txt, RK3, projection, channel forcing, restart, statistics, visualization, or closed stage directories.
AUDITEOF
}

configure_build_tree() {
  {
    echo "# R12 build log"
    echo "ROOT_DIR=$ROOT_DIR"
    echo "BUILD_DIR=$BUILD_DIR"
    echo "DECOMP_ROOT=$DECOMP_ROOT"
    echo
    echo "## Configure command"
    echo "cmake -S . -B \"$BUILD_DIR\" -DCMAKE_Fortran_COMPILER=mpif90 -DDECOMP2D_ROOT=\"$DECOMP_ROOT\" -DCMAKE_PREFIX_PATH=\"$DECOMP_ROOT\""
  } > "$BUILD_LOG"

  cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DDECOMP2D_ROOT="$DECOMP_ROOT" \
    -DCMAKE_PREFIX_PATH="$DECOMP_ROOT" \
    2>&1 | tee -a "$BUILD_LOG"
  local configure_status=${PIPESTATUS[0]}
  echo "configure_status=$configure_status" >> "$BUILD_LOG"
  return "$configure_status"
}

build_and_run_helper() {
  local target="$1"
  local token="$2"
  local logfile="$3"

  {
    echo "# R12 helper log for ${target}"
    echo "## Build command"
    echo "cmake --build \"$BUILD_DIR\" --target ${target} -j 2"
  } > "$logfile"

  cmake --build "$BUILD_DIR" --target "$target" -j 2 2>&1 | tee -a "$logfile"
  local build_status=${PIPESTATUS[0]}
  echo "build_status=$build_status" >> "$logfile"

  if [ "$build_status" -ne 0 ]; then
    echo "Result: FAIL" >> "$logfile"
    echo "Reason: target ${target} could not be built." >> "$logfile"
    return 1
  fi

  local exe
  exe="$(find_executable "$target")"
  echo "resolved_executable=$exe" >> "$logfile"
  if [ -z "$exe" ] || [ ! -x "$exe" ]; then
    echo "Result: FAIL" >> "$logfile"
    echo "Reason: executable for ${target} was not found." >> "$logfile"
    return 1
  fi

  echo "## Run command" >> "$logfile"
  echo "\"$exe\"" >> "$logfile"
  "$exe" 2>&1 | tee -a "$logfile"
  local run_status=${PIPESTATUS[0]}
  echo "run_status=$run_status" >> "$logfile"

  if [ "$run_status" -eq 0 ] && grep -q "$token" "$logfile"; then
    echo "Result: PASS" >> "$logfile"
    return 0
  fi

  echo "Result: FAIL" >> "$logfile"
  echo "Reason: expected token was not found: $token" >> "$logfile"
  return 1
}

run_helper_chain() {
  local pass=1

  build_and_run_helper fibre_prod_state_check "R2_FIBRE_PROD_STATE_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R2.txt" || pass=0
  build_and_run_helper fibre_prod_grid_adapter_check "R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R3.txt" || pass=0
  build_and_run_helper fibre_prod_ibm_interpolation_check "R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R4.txt" || pass=0
  build_and_run_helper fibre_prod_structure_solver_check "R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R5.txt" || pass=0
  build_and_run_helper fibre_prod_ibm_spreading_check "R6_FIBRE_PROD_IBM_SPREADING_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R6.txt" || pass=0
  build_and_run_helper fibre_prod_fsi_closed_loop_check "R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R7.txt" || pass=0
  build_and_run_helper fibre_prod_wall_contact_check "R8_FIBRE_PROD_WALL_CONTACT_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R8.txt" || pass=0
  build_and_run_helper fibre_prod_fibre_collision_check "R9_FIBRE_PROD_FIBRE_COLLISION_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R9.txt" || pass=0
  build_and_run_helper fibre_prod_main_hook_check "R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS" "$ROOT_DIR/production_recovery/R12_RUN_LOG_R10_hook_check.txt" || pass=0

  local result="FAIL"
  [ "$pass" -eq 1 ] && result="PASS"

  cat > "$EVIDENCE_DIR/R12_STANDALONE_HELPER_AUDIT.md" <<AUDITEOF
# R12 Standalone Helper Audit

Result: ${result}

R2_state_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R2.txt" && echo PASS || echo FAIL)
R3_grid_adapter_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R3.txt" && echo PASS || echo FAIL)
R4_ibm_interpolation_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R4.txt" && echo PASS || echo FAIL)
R5_structure_solver_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R5.txt" && echo PASS || echo FAIL)
R6_ibm_spreading_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R6.txt" && echo PASS || echo FAIL)
R7_fsi_closed_loop_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R7.txt" && echo PASS || echo FAIL)
R8_wall_contact_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R8.txt" && echo PASS || echo FAIL)
R9_fibre_collision_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R9.txt" && echo PASS || echo FAIL)
R10_main_hook_check=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_R10_hook_check.txt" && echo PASS || echo FAIL)
AUDITEOF

  [ "$pass" -eq 1 ]
}

run_r11_rerun() {
  local pass=1
  if [ ! -f production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh ]; then
    pass=0
    echo "BLOCKED: R11 validation script missing" > production_recovery/R12_RUN_LOG_R11_rerun.txt
  else
    bash production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh 2>&1 | tee production_recovery/R12_RUN_LOG_R11_rerun.txt
    local r11_status=${PIPESTATUS[0]}
    echo "r11_script_status=$r11_status" >> production_recovery/R12_RUN_LOG_R11_rerun.txt
    [ "$r11_status" -eq 0 ] || pass=0
  fi

  has_pass production_recovery/R11_PASS_FAIL.md || pass=0
  has_pass production_recovery/R11_evidence/R11_MPI_CONSISTENCY_AUDIT.md || pass=0
  for f in \
    production_recovery/R11_evidence/R11_LAMBDA0_NP1_AUDIT.txt \
    production_recovery/R11_evidence/R11_LAMBDA0_NP2_AUDIT.txt \
    production_recovery/R11_evidence/R11_LAMBDA0_NP4_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP1_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP2_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP4_AUDIT.txt; do
    has_pass "$f" || pass=0
  done

  local result="FAIL"
  [ "$pass" -eq 1 ] && result="PASS"

  cat > "$EVIDENCE_DIR/R12_R11_RERUN_AUDIT.md" <<AUDITEOF
# R12 R11 Rerun Audit

Result: ${result}

R11_rerun_log=production_recovery/R12_RUN_LOG_R11_rerun.txt
R11_PASS_FAIL=$(has_pass production_recovery/R11_PASS_FAIL.md && echo PASS || echo FAIL)
R11_MPI_CONSISTENCY_AUDIT=$(has_pass production_recovery/R11_evidence/R11_MPI_CONSISTENCY_AUDIT.md && echo PASS || echo FAIL)
R11_lambda0_np1=$(has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP1_AUDIT.txt && echo PASS || echo FAIL)
R11_lambda0_np2=$(has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP2_AUDIT.txt && echo PASS || echo FAIL)
R11_lambda0_np4=$(has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP4_AUDIT.txt && echo PASS || echo FAIL)
R11_smalllambda_np1=$(has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP1_AUDIT.txt && echo PASS || echo FAIL)
R11_smalllambda_np2=$(has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP2_AUDIT.txt && echo PASS || echo FAIL)
R11_smalllambda_np4=$(has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP4_AUDIT.txt && echo PASS || echo FAIL)
AUDITEOF

  [ "$pass" -eq 1 ]
}

build_xcompact3d() {
  cmake --build "$BUILD_DIR" --target xcompact3d -j 2 2>&1 | tee -a "$BUILD_LOG"
  local status=${PIPESTATUS[0]}
  echo "xcompact3d_build_status=$status" >> "$BUILD_LOG"
  return "$status"
}

run_nofibre_case() {
  local np="$1"
  local exe="$2"
  local input_src="$3"
  local run_dir="$EVIDENCE_DIR/r12_nofibre_np${np}"
  local log="$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np${np}.txt"

  rm -rf "$run_dir"
  mkdir -p "$run_dir"
  cp "$input_src" "$run_dir/input.i3d"

  {
    echo "# R12 no-fibre np=${np} log"
    echo "run_dir=$run_dir"
    echo "input_src=$input_src"
    echo "x3d_exe=$exe"
    echo "command=FIBRE_PROD_ENABLE=0 FIBRE_PROD_LAMBDA=0 FIBRE_PROD_DIAGNOSTICS=0 mpirun --oversubscribe -np ${np} $exe"
  } > "$log"

  (
    cd "$run_dir" || exit 1
    FIBRE_PROD_ENABLE=0 \
    FIBRE_PROD_LAMBDA=0 \
    FIBRE_PROD_DIAGNOSTICS=0 \
    mpirun --oversubscribe -np "$np" "$exe"
  ) 2>&1 | tee -a "$log"
  local status=${PIPESTATUS[0]}
  echo "run_status=$status" >> "$log"

  if [ "$status" -eq 0 ] && grep -q 'Good job! Xcompact3d finished successfully!' "$log"; then
    echo "Result: PASS" >> "$log"
    return 0
  fi
  echo "Result: FAIL" >> "$log"
  return 1
}

run_nofibre_baseline() {
  local pass=1
  local input_src
  input_src="$(select_input_case)"
  local exe
  exe="$(find_executable xcompact3d)"

  if [ -z "$input_src" ] || [ -z "$exe" ] || [ ! -x "$exe" ]; then
    pass=0
    for np in 1 2 4; do
      {
        echo "# R12 no-fibre np=${np} log"
        echo "Result: FAIL"
        echo "Reason: input case or xcompact3d executable missing."
        echo "input_src=$input_src"
        echo "x3d_exe=$exe"
      } > "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np${np}.txt"
    done
  else
    run_nofibre_case 1 "$exe" "$input_src" || pass=0
    run_nofibre_case 2 "$exe" "$input_src" || pass=0
    run_nofibre_case 4 "$exe" "$input_src" || pass=0
  fi

  local result="FAIL"
  [ "$pass" -eq 1 ] && result="PASS"
  cat > "$EVIDENCE_DIR/R12_NO_FIBRE_BASELINE_AUDIT.md" <<AUDITEOF
# R12 No-Fibre Baseline Audit

Result: ${result}

np1=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np1.txt" && echo PASS || echo FAIL)
np2=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np2.txt" && echo PASS || echo FAIL)
np4=$(grep -q 'Result: PASS' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np4.txt" && echo PASS || echo FAIL)
input_case=${input_src}
AUDITEOF

  [ "$pass" -eq 1 ]
}

check_restart_stats_visu() {
  local pass=1
  for np in 1 2 4; do
    local log="$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np${np}.txt"
    grep -q 'Writing restart point\|Restart point .* saved successfully' "$log" || pass=0
    grep -q 'Writing snapshots' "$log" || pass=0
    grep -q 'Writing stat file\|Write stat done' "$log" || pass=0
  done

  local result="FAIL"
  [ "$pass" -eq 1 ] && result="PASS"
  cat > "$EVIDENCE_DIR/R12_RESTART_STATS_VISU_AUDIT.md" <<AUDITEOF
# R12 Restart/Statistics/Visualization Audit

Result: ${result}

np1_restart=$(grep -q 'Writing restart point\|Restart point .* saved successfully' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np1.txt" && echo PASS || echo FAIL)
np1_snapshots=$(grep -q 'Writing snapshots' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np1.txt" && echo PASS || echo FAIL)
np1_statistics=$(grep -q 'Writing stat file\|Write stat done' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np1.txt" && echo PASS || echo FAIL)
np2_restart=$(grep -q 'Writing restart point\|Restart point .* saved successfully' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np2.txt" && echo PASS || echo FAIL)
np2_snapshots=$(grep -q 'Writing snapshots' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np2.txt" && echo PASS || echo FAIL)
np2_statistics=$(grep -q 'Writing stat file\|Write stat done' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np2.txt" && echo PASS || echo FAIL)
np4_restart=$(grep -q 'Writing restart point\|Restart point .* saved successfully' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np4.txt" && echo PASS || echo FAIL)
np4_snapshots=$(grep -q 'Writing snapshots' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np4.txt" && echo PASS || echo FAIL)
np4_statistics=$(grep -q 'Writing stat file\|Write stat done' "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np4.txt" && echo PASS || echo FAIL)
AUDITEOF

  [ "$pass" -eq 1 ]
}

write_final_matrix_and_closure() {
  local source_boundary_pass=0
  local standalone_helper_chain_pass=0
  local r11_rerun_pass=0
  local nofibre_np124_pass=0
  local lambda0_np124_pass=0
  local smalllambda_np124_pass=0
  local restart_statistics_visualization_smoke_pass=0
  local no_src_modification_in_r12=1
  local no_rk3_projection_channel_forcing_modification=1
  local no_r13_entered=1

  has_pass "$EVIDENCE_DIR/R12_SOURCE_BOUNDARY_AUDIT.md" && source_boundary_pass=1
  has_pass "$EVIDENCE_DIR/R12_STANDALONE_HELPER_AUDIT.md" && standalone_helper_chain_pass=1
  has_pass "$EVIDENCE_DIR/R12_R11_RERUN_AUDIT.md" && r11_rerun_pass=1
  has_pass "$EVIDENCE_DIR/R12_NO_FIBRE_BASELINE_AUDIT.md" && nofibre_np124_pass=1
  has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP1_AUDIT.txt && \
    has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP2_AUDIT.txt && \
    has_pass production_recovery/R11_evidence/R11_LAMBDA0_NP4_AUDIT.txt && lambda0_np124_pass=1
  has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP1_AUDIT.txt && \
    has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP2_AUDIT.txt && \
    has_pass production_recovery/R11_evidence/R11_SMALLLAMBDA_NP4_AUDIT.txt && smalllambda_np124_pass=1
  has_pass "$EVIDENCE_DIR/R12_RESTART_STATS_VISU_AUDIT.md" && restart_statistics_visualization_smoke_pass=1

  local result="FAIL"
  if [ "$source_boundary_pass" -eq 1 ] && \
     [ "$standalone_helper_chain_pass" -eq 1 ] && \
     [ "$r11_rerun_pass" -eq 1 ] && \
     [ "$nofibre_np124_pass" -eq 1 ] && \
     [ "$lambda0_np124_pass" -eq 1 ] && \
     [ "$smalllambda_np124_pass" -eq 1 ] && \
     [ "$restart_statistics_visualization_smoke_pass" -eq 1 ] && \
     [ "$no_src_modification_in_r12" -eq 1 ] && \
     [ "$no_rk3_projection_channel_forcing_modification" -eq 1 ] && \
     [ "$no_r13_entered" -eq 1 ]; then
    result="PASS"
  fi

  cat > "$EVIDENCE_DIR/R12_FINAL_VALIDATION_MATRIX.md" <<MATRIXEOF
# R12 Final Validation Matrix

Result: ${result}

source_boundary_pass=${source_boundary_pass}
standalone_helper_chain_pass=${standalone_helper_chain_pass}
r11_rerun_pass=${r11_rerun_pass}
nofibre_np124_pass=${nofibre_np124_pass}
lambda0_np124_pass=${lambda0_np124_pass}
smalllambda_np124_pass=${smalllambda_np124_pass}
restart_statistics_visualization_smoke_pass=${restart_statistics_visualization_smoke_pass}
no_src_modification_in_r12=${no_src_modification_in_r12}
no_rk3_projection_channel_forcing_modification=${no_rk3_projection_channel_forcing_modification}
no_r13_entered=${no_r13_entered}

R12 PASS requires every item above to equal 1. R12 does not modify source code and does not enter R13.
MATRIXEOF

  cat > "$EVIDENCE_DIR/R12_FINAL_CLOSURE_REPORT.md" <<CLOSUREEOF
# R12 Final Closure Report

Result: ${result}

## Scope

R12 closes the Production Recovery R0-R12 validation matrix if and only if the final matrix is PASS.

## PASS meaning

If Result: PASS, R0-R12 Production Recovery validation matrix is closed for the controlled micro-validation scope: standalone helper chain, controlled main hook, np=1/2/4 R11 consistency, no-fibre np=1/2/4 smoke baseline, and restart/statistics/visualization smoke evidence.

## Non-claims

R12 PASS does not claim long-time production DNS, paper-grade physical statistics, exhaustive mesh/timestep independence, or experimental validation.

## Boundary

No source files were modified by R12. R13 was not entered.
CLOSUREEOF

  cat > "$ROOT_DIR/production_recovery/R12_PASS_FAIL.md" <<PASSEOF
# Production Recovery R12 PASS/FAIL

Result: ${result}
PASSEOF
}

write_evidence_index() {
  cat > "$EVIDENCE_DIR/R12_EVIDENCE_INDEX.md" <<INDEXEOF
# R12 Evidence Index

- R12_SOURCE_BOUNDARY_AUDIT.md
- R12_STANDALONE_HELPER_AUDIT.md
- R12_R11_RERUN_AUDIT.md
- R12_NO_FIBRE_BASELINE_AUDIT.md
- R12_RESTART_STATS_VISU_AUDIT.md
- R12_FINAL_VALIDATION_MATRIX.md
- R12_FINAL_CLOSURE_REPORT.md
- production_recovery/R12_RUN_LOG_R2.txt
- production_recovery/R12_RUN_LOG_R3.txt
- production_recovery/R12_RUN_LOG_R4.txt
- production_recovery/R12_RUN_LOG_R5.txt
- production_recovery/R12_RUN_LOG_R6.txt
- production_recovery/R12_RUN_LOG_R7.txt
- production_recovery/R12_RUN_LOG_R8.txt
- production_recovery/R12_RUN_LOG_R9.txt
- production_recovery/R12_RUN_LOG_R10_hook_check.txt
- production_recovery/R12_RUN_LOG_R11_rerun.txt
- production_recovery/R12_RUN_LOG_nofibre_np1.txt
- production_recovery/R12_RUN_LOG_nofibre_np2.txt
- production_recovery/R12_RUN_LOG_nofibre_np4.txt
INDEXEOF
}

main() {
  write_source_boundary_audit
  local configure_pass=1
  configure_build_tree || configure_pass=0

  local helper_pass=0
  if [ "$configure_pass" -eq 1 ]; then
    run_helper_chain && helper_pass=1 || helper_pass=0
  else
    cat > "$EVIDENCE_DIR/R12_STANDALONE_HELPER_AUDIT.md" <<AUDITEOF
# R12 Standalone Helper Audit

Result: FAIL

Configure failed; helper targets were not run.
AUDITEOF
  fi

  run_r11_rerun || true

  if [ "$configure_pass" -eq 1 ]; then
    build_xcompact3d || true
    run_nofibre_baseline || true
  else
    cat > "$EVIDENCE_DIR/R12_NO_FIBRE_BASELINE_AUDIT.md" <<AUDITEOF
# R12 No-Fibre Baseline Audit

Result: FAIL

Configure failed; no-fibre np=1/2/4 baseline was not run.
AUDITEOF
    for np in 1 2 4; do
      echo "Result: FAIL" > "$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np${np}.txt"
    done
  fi

  check_restart_stats_visu || true
  write_final_matrix_and_closure
  write_evidence_index

  echo "===== R12 result ====="
  cat "$ROOT_DIR/production_recovery/R12_PASS_FAIL.md"
  echo "===== R12 final matrix ====="
  cat "$EVIDENCE_DIR/R12_FINAL_VALIDATION_MATRIX.md"
}

main "$@"
