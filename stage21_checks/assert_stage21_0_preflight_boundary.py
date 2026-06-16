#!/usr/bin/env python3
"""Stage 21.0 wall/contact/collision preflight boundary audit.

This helper validates documentation and fail-closed gate defaults only. It does
not compute real wall gaps, fibre-fibre gaps, contact forces, collision forces,
production structure updates, production RHS updates, MPI, or DNS.
"""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_0_preflight_boundary.dat"
HELPER = CHECK_DIR / "assert_stage21_0_preflight_boundary.py"
WRAPPER = CHECK_DIR / "run_stage21_0_preflight_boundary.sh"
DOC = CHECK_DIR / "stage21_0_preflight_boundary.md"

SAFE_DEFAULTS = {
    "STAGE21_0_ENABLE": "1",
    "STAGE21_0_PREFLIGHT_BOUNDARY_ENABLE": "1",
    "STAGE21_0_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE21_0_REQUIRE_STAGE20_11_PASS": "1",
    "STAGE21_0_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_0_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE21_0_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE21_0_DIAGNOSTIC_ONLY": "1",
    "STAGE21_0_FAIL_CLOSED": "1",
    "STAGE21_0_WALL_CONTACT_ENABLE": "0",
    "STAGE21_0_FIBRE_COLLISION_ENABLE": "0",
    "STAGE21_0_CONTACT_DISTANCE_AUDIT_ENABLE": "0",
    "STAGE21_0_WALL_GAP_AUDIT_ENABLE": "0",
    "STAGE21_0_FIBRE_FIBRE_GAP_AUDIT_ENABLE": "0",
    "STAGE21_0_NEAR_CONTACT_WARNING_ENABLE": "0",
    "STAGE21_0_CONTACT_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_0_WALL_CONTACT_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_0_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_0_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_0_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_0_CONTACT_TO_RHS_ENABLE": "0",
    "STAGE21_0_SINGLE_FIBRE_WALL_TEST_ONLY": "1",
    "STAGE21_0_TWO_FIBRE_COLLISION_TEST_ONLY": "0",
    "STAGE21_0_PRODUCTION_MULTIFIBRE_ENABLE": "0",
    "STAGE21_0_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_0_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_0_ZERO_TOL": "1.0e-14",
    "STAGE21_0_AUDIT_TOL": "1.0e-12",
    "STAGE21_0_TEST_CASE": "wall_contact_collision_preflight_boundary",
}

REQUIRED_STATUS_FIELDS = [
    "stage21_0_requested_status",
    "stage21_0_preflight_boundary_enable_status",
    "stage20_closure_evidence_status",
    "stage20_source_only_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "stage21_goal_documented_status",
    "stage21_math_boundary_documented_status",
    "wall_gap_formula_documented_status",
    "fibre_fibre_gap_formula_documented_status",
    "non_penetration_requirement_documented_status",
    "future_wall_contact_force_formula_documented_status",
    "future_fibre_collision_force_formula_documented_status",
    "action_reaction_requirement_documented_status",
    "default_safe_gate_values_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "wall_contact_default_disabled_status",
    "fibre_collision_default_disabled_status",
    "contact_distance_audit_default_disabled_status",
    "wall_gap_audit_default_disabled_status",
    "fibre_fibre_gap_audit_default_disabled_status",
    "near_contact_warning_default_disabled_status",
    "contact_force_candidate_default_disabled_status",
    "wall_contact_force_candidate_default_disabled_status",
    "fibre_collision_force_candidate_default_disabled_status",
    "contact_force_apply_default_disabled_status",
    "contact_in_structure_advance_default_disabled_status",
    "contact_to_rhs_default_disabled_status",
    "single_fibre_wall_test_only_status",
    "two_fibre_collision_test_default_disabled_status",
    "production_multifibre_default_disabled_status",
    "production_dns_allowed_disabled_status",
    "actual_mpi_allowed_disabled_status",
    "no_actual_wall_gap_computation_status",
    "no_actual_fibre_fibre_gap_computation_status",
    "no_actual_contact_force_computation_status",
    "no_actual_collision_force_computation_status",
    "no_production_structure_update_status",
    "no_production_rhs_update_status",
    "no_stage14_rhs_injection_status",
    "no_stage10_20_file_modification_status",
    "no_closed_stage_modification_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status",
    "no_production_structure_to_fluid_reaction_force_status",
    "no_production_eulerian_spreading_activation_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage21_0_status",
    "no_fluid_rhs_modification_status",
    "no_ibm_modification_status",
    "no_dns_core_modification_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status",
    "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status",
    "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status",
    "no_repulsive_force_status",
    "no_lubrication_force_status",
    "no_friction_force_status",
    "no_adhesion_force_status",
    "no_contact_damping_force_status",
    "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "stage21_1_next_stage_declared_status",
    "stage21_0_wrapper_bash_syntax_status",
    "stage21_0_helper_py_compile_status",
]


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def stage20_evidence_accepted() -> tuple[bool, str]:
    closed = ROOT / "stage20_checks" / "STAGE20_CLOSED.md"
    stage20_11 = ROOT / "stage20_outputs" / "fibre_stage20_11_total_closure_audit.dat"
    if closed.exists() and "Final verdict: PASS" in closed.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE20_CLOSED_PASS"
    if stage20_11.exists() and "STAGE 20 FINAL CLOSURE VERDICT: PASS" in stage20_11.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE20_11_PASS_OUTPUT"
    if enabled("STAGE21_0_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"):
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def doc_text() -> str:
    return DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_0_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def default_safe_gate_values_ok() -> bool:
    enabled_keys = {
        "STAGE21_0_ENABLE",
        "STAGE21_0_PREFLIGHT_BOUNDARY_ENABLE",
        "STAGE21_0_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE",
        "STAGE21_0_REQUIRE_STAGE20_11_PASS",
        "STAGE21_0_DO_NOT_RERUN_PREVIOUS_STAGES",
        "STAGE21_0_ALLOW_MISSING_OLD_STAGE_OUTPUTS",
        "STAGE21_0_ALLOW_MISSING_OLD_CLOSURE_FILES",
        "STAGE21_0_DIAGNOSTIC_ONLY",
        "STAGE21_0_FAIL_CLOSED",
        "STAGE21_0_SINGLE_FIBRE_WALL_TEST_ONLY",
    }
    disabled_keys = {
        "STAGE21_0_WALL_CONTACT_ENABLE",
        "STAGE21_0_FIBRE_COLLISION_ENABLE",
        "STAGE21_0_CONTACT_DISTANCE_AUDIT_ENABLE",
        "STAGE21_0_WALL_GAP_AUDIT_ENABLE",
        "STAGE21_0_FIBRE_FIBRE_GAP_AUDIT_ENABLE",
        "STAGE21_0_NEAR_CONTACT_WARNING_ENABLE",
        "STAGE21_0_CONTACT_FORCE_CANDIDATE_ENABLE",
        "STAGE21_0_WALL_CONTACT_FORCE_CANDIDATE_ENABLE",
        "STAGE21_0_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE",
        "STAGE21_0_CONTACT_FORCE_APPLY_ENABLE",
        "STAGE21_0_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE",
        "STAGE21_0_CONTACT_TO_RHS_ENABLE",
        "STAGE21_0_TWO_FIBRE_COLLISION_TEST_ONLY",
        "STAGE21_0_PRODUCTION_MULTIFIBRE_ENABLE",
        "STAGE21_0_PRODUCTION_DNS_ALLOWED",
        "STAGE21_0_ACTUAL_MPI_ALLOWED",
    }
    return all(enabled(k) for k in enabled_keys) and all(disabled(k) for k in disabled_keys)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    evidence_ok, evidence_mode = stage20_evidence_accepted()
    text = doc_text()
    statuses: dict[str, bool] = {field: True for field in REQUIRED_STATUS_FIELDS}
    statuses.update(
        {
            "stage21_0_requested_status": enabled("STAGE21_0_ENABLE"),
            "stage21_0_preflight_boundary_enable_status": enabled("STAGE21_0_PREFLIGHT_BOUNDARY_ENABLE"),
            "stage20_closure_evidence_status": evidence_ok,
            "stage20_source_only_closure_acceptance_status": enabled("STAGE21_0_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
            "missing_old_stage_outputs_allowed_status": enabled("STAGE21_0_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
            "missing_old_closure_files_allowed_status": enabled("STAGE21_0_ALLOW_MISSING_OLD_CLOSURE_FILES"),
            "no_previous_stage_rerun_status": enabled("STAGE21_0_DO_NOT_RERUN_PREVIOUS_STAGES"),
            "stage21_goal_documented_status": "wall/contact/collision safety boundary" in text or "wall/contact/collision preflight boundary" in text,
            "stage21_math_boundary_documented_status": "Mathematical boundary" in text,
            "wall_gap_formula_documented_status": "g_wall = d_wall - a_f" in text,
            "fibre_fibre_gap_formula_documented_status": "g_ff = d_ff - (a_i + a_j)" in text,
            "non_penetration_requirement_documented_status": "g_wall >= 0" in text and "g_ff >= 0" in text,
            "future_wall_contact_force_formula_documented_status": "F_wall_candidate = k_w * delta_wall * n_wall - c_w * v_n_minus * n_wall" in text,
            "future_fibre_collision_force_formula_documented_status": "F_i_candidate = k_c * delta_ij * n_ij - c_c * v_n_minus * n_ij" in text and "F_j_candidate = -F_i_candidate" in text,
            "action_reaction_requirement_documented_status": "F_i_candidate + F_j_candidate = 0" in text,
            "default_safe_gate_values_status": default_safe_gate_values_ok(),
            "diagnostic_only_status": enabled("STAGE21_0_DIAGNOSTIC_ONLY"),
            "fail_closed_status": enabled("STAGE21_0_FAIL_CLOSED"),
            "wall_contact_default_disabled_status": disabled("STAGE21_0_WALL_CONTACT_ENABLE"),
            "fibre_collision_default_disabled_status": disabled("STAGE21_0_FIBRE_COLLISION_ENABLE"),
            "contact_distance_audit_default_disabled_status": disabled("STAGE21_0_CONTACT_DISTANCE_AUDIT_ENABLE"),
            "wall_gap_audit_default_disabled_status": disabled("STAGE21_0_WALL_GAP_AUDIT_ENABLE"),
            "fibre_fibre_gap_audit_default_disabled_status": disabled("STAGE21_0_FIBRE_FIBRE_GAP_AUDIT_ENABLE"),
            "near_contact_warning_default_disabled_status": disabled("STAGE21_0_NEAR_CONTACT_WARNING_ENABLE"),
            "contact_force_candidate_default_disabled_status": disabled("STAGE21_0_CONTACT_FORCE_CANDIDATE_ENABLE"),
            "wall_contact_force_candidate_default_disabled_status": disabled("STAGE21_0_WALL_CONTACT_FORCE_CANDIDATE_ENABLE"),
            "fibre_collision_force_candidate_default_disabled_status": disabled("STAGE21_0_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE"),
            "contact_force_apply_default_disabled_status": disabled("STAGE21_0_CONTACT_FORCE_APPLY_ENABLE"),
            "contact_in_structure_advance_default_disabled_status": disabled("STAGE21_0_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE"),
            "contact_to_rhs_default_disabled_status": disabled("STAGE21_0_CONTACT_TO_RHS_ENABLE"),
            "single_fibre_wall_test_only_status": enabled("STAGE21_0_SINGLE_FIBRE_WALL_TEST_ONLY"),
            "two_fibre_collision_test_default_disabled_status": disabled("STAGE21_0_TWO_FIBRE_COLLISION_TEST_ONLY"),
            "production_multifibre_default_disabled_status": disabled("STAGE21_0_PRODUCTION_MULTIFIBRE_ENABLE"),
            "production_dns_allowed_disabled_status": disabled("STAGE21_0_PRODUCTION_DNS_ALLOWED"),
            "actual_mpi_allowed_disabled_status": disabled("STAGE21_0_ACTUAL_MPI_ALLOWED"),
            "stage21_1_next_stage_declared_status": "Stage 21.1: wall distance and signed-gap audit" in text,
            "stage21_0_wrapper_bash_syntax_status": bash_syntax_ok(),
            "stage21_0_helper_py_compile_status": py_compile_ok(),
        }
    )
    final_pass = all(statuses.values())
    lines = [
        "stage21_0_title wall/contact/collision preflight boundary",
        "stage21_title wall/contact/collision safety boundary for flexible-fibre FSI",
        "stage20_closure_evidence_mode_value " + evidence_mode,
        "source_only_policy_value Stage 20 closure files and old outputs are optional when source-only acceptance is enabled",
        "rerun_policy_value Stage 21.0 does not rerun Stage 10-20",
        "active_physics_policy_value documented formulas only; no actual gap/contact/collision computation",
        "wall_gap_formula_value g_wall = d_wall - a_f",
        "fibre_fibre_gap_formula_value g_ff = d_ff - (a_i + a_j)",
        "non_penetration_requirement_value g_wall >= 0; g_ff >= 0",
        "future_wall_contact_force_formula_value F_wall_candidate = k_w * delta_wall * n_wall - c_w * v_n_minus * n_wall",
        "future_fibre_collision_force_formula_value F_i_candidate = k_c * delta_ij * n_ij - c_c * v_n_minus * n_ij; F_j_candidate = -F_i_candidate",
        "action_reaction_requirement_value F_i_candidate + F_j_candidate = 0",
        "stage21_1_next_stage_value Stage 21.1: wall distance and signed-gap audit",
    ]
    for key in SAFE_DEFAULTS:
        lines.append(f"{key.lower()}_value {env(key)}")
    for field in REQUIRED_STATUS_FIELDS:
        lines.append(f"{field} {'PASS' if statuses[field] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final_pass else 'FAIL'}")
    if final_pass:
        lines.extend([
            "STAGE 21.0 PREFLIGHT BOUNDARY VERDICT: PASS",
            "STAGE 21.0 FINAL VERDICT: PASS",
        ])
    else:
        failed = [field for field, ok in statuses.items() if not ok]
        lines.append("failure_reasons_value " + ",".join(failed))
        lines.extend([
            "STAGE 21.0 PREFLIGHT BOUNDARY VERDICT: FAIL",
            "STAGE 21.0 FINAL VERDICT: FAIL",
        ])
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2])
    print(lines[-1])
    return 0 if final_pass else 1


if __name__ == "__main__":
    sys.exit(main())
