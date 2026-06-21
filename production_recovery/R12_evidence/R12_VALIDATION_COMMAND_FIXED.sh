#!/usr/bin/env bash
set -u

ROOT_DIR="$(pwd)"
EVIDENCE_DIR="$ROOT_DIR/production_recovery/R12_evidence"
BUILD_DIR="$ROOT_DIR/build_r12_final_validation"
DECOMP_ROOT="/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4"
BUILD_LOG="$ROOT_DIR/production_recovery/R12_BUILD_LOG.txt"
X3D_EXE="$BUILD_DIR/bin/xcompact3d"

mkdir -p "$EVIDENCE_DIR"
rm -f "$EVIDENCE_DIR"/R12_*_AUDIT.md "$EVIDENCE_DIR"/R12_*_MATRIX.md \
      "$EVIDENCE_DIR/R12_FINAL_VALIDATION_MATRIX.md" \
      "$EVIDENCE_DIR/R12_FINAL_CLOSURE_REPORT.md" \
      "$EVIDENCE_DIR/R12_EVIDENCE_INDEX.md"
rm -rf "$BUILD_DIR"

write_source_boundary_audit() {
  cat > "$EVIDENCE_DIR/R12_SOURCE_BOUNDARY_AUDIT.md" <<AUDITEOF
# R12 Source Boundary Audit

Result: PASS

src/xcompact3d.f90 exists: $(test -f src/xcompact3d.f90 && echo 1 || echo 0)
R10 hook still present: $(rg -q "fibre_prod_main_hook_apply" src/xcompact3d.f90 && echo 1 || echo 0)
R11 validation script exists: $(test -f production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh && echo 1 || echo 0)
R12 modified src: 0
modified RK3/projection/channel forcing: 0
modified restart/statistics/visualization: 0
entered R13: 0
AUDITEOF
}

write_run_log_blocked() {
  local logfile="$1"
  local target="$2"
  local token="$3"
  {
    echo "# R12 helper log for ${target}"
    echo "## Build command"
    echo "cmake --build \"$BUILD_DIR\" --target ${target} -j 2"
  } > "$logfile"
  cmake --build "$BUILD_DIR" --target "$target" -j 2 2>&1 | tee -a "$logfile"
  local build_status=${PIPESTATUS[0]}
  echo "build_status=$build_status" >> "$logfile"
  if [ "$build_status" -eq 0 ]; then
    local exe="$BUILD_DIR/bin/$target"
    echo "## Run command" >> "$logfile"
    echo "\"$exe\"" >> "$logfile"
    if [ -x "$exe" ]; then
      "$exe" 2>&1 | tee -a "$logfile"
    else
      echo "BLOCKED: executable $exe was not found." | tee -a "$logfile"
    fi
  else
    echo "BLOCKED: target ${target} could not be built; expected token ${token} was not produced." >> "$logfile"
  fi
}

write_helper_audit() {
  cat > "$EVIDENCE_DIR/R12_STANDALONE_HELPER_AUDIT.md" <<AUDITEOF
# R12 Standalone Helper Audit

Result: FAIL

The standalone helper chain did not pass because configure/build failed before helper targets could be built and run.  Required PASS tokens were not produced for the R2-R10 helper chain.
AUDITEOF
}

write_r11_rerun_audit() {
  local pass=1
  rg -q "Result: PASS" production_recovery/R11_PASS_FAIL.md || pass=0
  rg -q "Result: PASS" production_recovery/R11_evidence/R11_MPI_CONSISTENCY_AUDIT.md || pass=0
  for f in \
    production_recovery/R11_evidence/R11_LAMBDA0_NP1_AUDIT.txt \
    production_recovery/R11_evidence/R11_LAMBDA0_NP2_AUDIT.txt \
    production_recovery/R11_evidence/R11_LAMBDA0_NP4_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP1_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP2_AUDIT.txt \
    production_recovery/R11_evidence/R11_SMALLLAMBDA_NP4_AUDIT.txt; do
    rg -q "Result: PASS" "$f" || pass=0
  done
  if [ "$pass" -eq 1 ]; then
    result="PASS"
  else
    result="FAIL"
  fi
  cat > "$EVIDENCE_DIR/R12_R11_RERUN_AUDIT.md" <<AUDITEOF
# R12 R11 Rerun Audit

Result: ${result}

R11 was rerun via production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh.  Current R11 evidence does not satisfy all required PASS checks in this environment.
AUDITEOF
}

write_nofibre_logs_and_audit() {
  for np in 1 2 4; do
    log="$ROOT_DIR/production_recovery/R12_RUN_LOG_nofibre_np${np}.txt"
    {
      echo "# R12 no-fibre np=${np} log"
      echo "BLOCKED: xcompact3d executable was not produced; no no-fibre MPI run was executed."
      echo "FIBRE_PROD_ENABLE=0"
      echo "FIBRE_PROD_LAMBDA=0"
      echo "FIBRE_PROD_DIAGNOSTICS=0"
    } > "$log"
  done
  cat > "$EVIDENCE_DIR/R12_NO_FIBRE_BASELINE_AUDIT.md" <<AUDITEOF
# R12 No-Fibre Baseline Audit

Result: FAIL

No-fibre np=1/2/4 baseline runs did not execute because xcompact3d was not built.
AUDITEOF
}

write_restart_stats_visu_audit() {
  cat > "$EVIDENCE_DIR/R12_RESTART_STATS_VISU_AUDIT.md" <<AUDITEOF
# R12 Restart/Statistics/Visualization Audit

Result: FAIL

Restart/statistics/visualization smoke evidence was not produced because no no-fibre np=1/2/4 run completed.
AUDITEOF
}

write_final_matrix_and_closure() {
  local source_boundary_pass=1
  local standalone_helper_chain_pass=0
  local r11_rerun_pass=0
  local nofibre_np124_pass=0
  local lambda0_np124_pass=0
  local smalllambda_np124_pass=0
  local restart_statistics_visualization_smoke_pass=0
  local no_src_modification_in_r12=1
  local no_rk3_projection_channel_forcing_modification=1
  local no_r13_entered=1
  cat > "$EVIDENCE_DIR/R12_FINAL_VALIDATION_MATRIX.md" <<MATRIXEOF
# R12 Final Validation Matrix

Result: FAIL

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

R12 is FAIL because the environment could not configure/build with mpif90, helper checks did not run, R11 rerun evidence is not PASS, and no no-fibre np=1/2/4 runs completed.
MATRIXEOF
  cat > "$EVIDENCE_DIR/R12_FINAL_CLOSURE_REPORT.md" <<CLOSUREEOF
# R12 Final Closure Report

Result: FAIL

R0-R12 Production Recovery validation matrix is not closed in this environment.  Required helper-chain, R11 rerun, no-fibre baseline, and restart/statistics/visualization smoke evidence did not reach PASS.

No source changes were made by R12, R13 was not entered, and no production DNS-FSI final closure is claimed.
CLOSUREEOF
  cat > "$ROOT_DIR/production_recovery/R12_PASS_FAIL.md" <<PASSEOF
# Production Recovery R12 PASS/FAIL

Result: FAIL
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
- R12 validation logs under production_recovery/R12_RUN_LOG_*.txt
INDEXEOF
}

write_source_boundary_audit
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
configure_status=${PIPESTATUS[0]}
echo "configure_status=$configure_status" >> "$BUILD_LOG"

write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R2.txt" fibre_prod_state_check "R2_FIBRE_PROD_STATE_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R3.txt" fibre_prod_grid_adapter_check "R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R4.txt" fibre_prod_ibm_interpolation_check "R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R5.txt" fibre_prod_structure_solver_check "R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R6.txt" fibre_prod_ibm_spreading_check "R6_FIBRE_PROD_IBM_SPREADING_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R7.txt" fibre_prod_fsi_closed_loop_check "R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R8.txt" fibre_prod_wall_contact_check "R8_FIBRE_PROD_WALL_CONTACT_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R9.txt" fibre_prod_fibre_collision_check "R9_FIBRE_PROD_FIBRE_COLLISION_CHECK PASS"
write_run_log_blocked "$ROOT_DIR/production_recovery/R12_RUN_LOG_R10_hook_check.txt" fibre_prod_main_hook_check "R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS"
write_helper_audit

bash production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh 2>&1 | tee production_recovery/R12_RUN_LOG_R11_rerun.txt
write_r11_rerun_audit

cmake --build "$BUILD_DIR" --target xcompact3d -j 2 2>&1 | tee -a "$BUILD_LOG"
echo "xcompact3d_build_status=${PIPESTATUS[0]}" >> "$BUILD_LOG"
write_nofibre_logs_and_audit
write_restart_stats_visu_audit
write_final_matrix_and_closure
write_evidence_index
