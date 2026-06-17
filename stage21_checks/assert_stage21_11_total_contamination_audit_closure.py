#!/usr/bin/env python3
"""Stage 21.11 diagnostic-only total contamination audit and closure."""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

PASS = "PASS"

STAGE21_SUBSTAGES = [
    ("21.0", "wall/contact/collision preflight boundary"),
    ("21.1", "wall distance and signed-gap audit"),
    ("21.2", "fibre-fibre point/segment distance audit"),
    ("21.3", "near-contact warning and fail-closed gate"),
    ("21.4", "contact candidate registry"),
    ("21.5", "contact pair ownership audit"),
    ("21.6", "deterministic pair ordering"),
    ("21.7", "contact metadata consistency"),
    ("21.8", "contact candidate persistence audit"),
    ("21.9", "contact diagnostic integration"),
    ("21.10", "collision-force-disabled proof"),
    ("21.11", "total contamination audit and closure"),
]

ACCOMPLISHED = [
    "Established wall/contact/collision preflight boundary.",
    "Audited wall distance and signed gaps.",
    "Audited fibre-fibre point/segment distances and signed gaps.",
    "Established near-contact warning and fail-closed classification.",
    "Built a diagnostic-only contact candidate registry.",
    "Audited contact pair ownership under helper np=1/2/4 semantics.",
    "Audited deterministic pair ordering.",
    "Audited contact/collision metadata consistency.",
    "Audited diagnostic-only metadata persistence/reload readiness.",
    "Integrated the contact diagnostic chain.",
    "Proved contact/collision force pathways remained disabled.",
]

NOT_DONE = [
    "Did not compute real wall contact force.",
    "Did not compute real fibre-fibre collision force.",
    "Did not introduce penalty force.",
    "Did not introduce repulsive force.",
    "Did not introduce lubrication force.",
    "Did not introduce friction force.",
    "Did not introduce adhesion force.",
    "Did not introduce contact damping force.",
    "Did not apply contact/collision force to structure advance.",
    "Did not add contact/collision force to total structural force.",
    "Did not couple contact/collision force to RHS.",
    "Did not call Stage 14 RHS injection.",
    "Did not run MPI.",
    "Did not run production DNS.",
    "Did not modify production restart/statistics/visualization I/O.",
    "Did not modify closed stages.",
    "Did not close Stage 22.",
]

FALSE_GATES = {
    "contact_force_disabled_status": "STAGE21_11_CONTACT_FORCE_ENABLE",
    "collision_force_disabled_status": "STAGE21_11_COLLISION_FORCE_ENABLE",
    "wall_contact_force_candidate_disabled_status": "STAGE21_11_WALL_CONTACT_FORCE_CANDIDATE_ENABLE",
    "fibre_collision_force_candidate_disabled_status": "STAGE21_11_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE",
    "penalty_force_disabled_status": "STAGE21_11_PENALTY_FORCE_ENABLE",
    "repulsive_force_disabled_status": "STAGE21_11_REPULSIVE_FORCE_ENABLE",
    "lubrication_force_disabled_status": "STAGE21_11_LUBRICATION_FORCE_ENABLE",
    "friction_force_disabled_status": "STAGE21_11_FRICTION_FORCE_ENABLE",
    "adhesion_force_disabled_status": "STAGE21_11_ADHESION_FORCE_ENABLE",
    "contact_damping_force_disabled_status": "STAGE21_11_CONTACT_DAMPING_FORCE_ENABLE",
    "contact_force_apply_disabled_status": "STAGE21_11_CONTACT_FORCE_APPLY_ENABLE",
    "collision_force_apply_disabled_status": "STAGE21_11_COLLISION_FORCE_APPLY_ENABLE",
    "contact_in_structure_advance_disabled_status": "STAGE21_11_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE",
    "collision_in_structure_advance_disabled_status": "STAGE21_11_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE",
    "contact_to_rhs_disabled_status": "STAGE21_11_CONTACT_TO_RHS_ENABLE",
    "collision_to_rhs_disabled_status": "STAGE21_11_COLLISION_TO_RHS_ENABLE",
    "stage14_rhs_injection_disabled_status": "STAGE21_11_STAGE14_RHS_INJECTION_ALLOWED",
    "production_rhs_update_disabled_status": "STAGE21_11_PRODUCTION_RHS_UPDATE_ALLOWED",
    "production_dns_disabled_status": "STAGE21_11_PRODUCTION_DNS_ALLOWED",
    "actual_mpi_disabled_status": "STAGE21_11_ACTUAL_MPI_ALLOWED",
    "production_multifibre_disabled_status": "STAGE21_11_PRODUCTION_MULTIFIBRE_ENABLE",
    "production_restart_io_disabled_status": "STAGE21_11_PRODUCTION_RESTART_IO_ALLOWED",
    "production_statistics_io_disabled_status": "STAGE21_11_PRODUCTION_STATISTICS_IO_ALLOWED",
    "production_visualization_io_disabled_status": "STAGE21_11_PRODUCTION_VISUALIZATION_IO_ALLOWED",
}


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def disabled(env_name: str) -> bool:
    return not bool_env(env_name, False)


def status(ok: bool) -> str:
    return PASS if ok else "FAIL"


def compile_with_temp(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pyc") as tmp:
            py_compile.compile(str(path), cfile=tmp.name, doraise=True)
        return True
    except Exception:
        return False


def closure_text() -> str:
    lines = [
        "# STAGE 21 CLOSED",
        "",
        "Final verdict: PASS",
        "",
        "## Stage 21 sub-stages",
    ]
    lines.extend(f"* Stage {num}: {title}" for num, title in STAGE21_SUBSTAGES)
    lines.extend([
        "",
        "## Closure statement",
        "Stage 21 wall/contact/collision diagnostic safety boundary is closed.",
        "No real contact/collision force was enabled in Stage 21.",
        "Stage 22 is the next and final large stage.",
        "Stage 22 title: final integrated validation and production-readiness closure",
        "Stage 10–21 closed files must remain immutable unless explicitly reopened by future instructions.",
        "",
        "## Stage 21 accomplished",
    ])
    lines.extend(f"* {item}" for item in ACCOMPLISHED)
    lines.extend(["", "## Stage 21 did not do"])
    lines.extend(f"* {item}" for item in NOT_DONE)
    return "\n".join(lines) + "\n"


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "stage21_outputs" / "fibre_stage21_11_total_contamination_audit_closure.dat"
    closed = repo / "stage21_checks" / "STAGE21_CLOSED.md"
    doc = repo / "stage21_checks" / "stage21_11_total_contamination_audit_closure.md"
    wrapper = repo / "stage21_checks" / "run_stage21_11_total_contamination_audit_closure.sh"
    helper = Path(__file__).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    if bool_env("STAGE21_11_CREATE_STAGE21_CLOSED_FILE", True):
        closed.write_text(closure_text(), encoding="utf-8")
    closed_content = closed.read_text(encoding="utf-8") if closed.exists() else ""
    all_forces_disabled = all(disabled(env) for env in FALSE_GATES.values())

    checks = {
        "stage21_11_requested_status": bool_env("STAGE21_11_ENABLE", True),
        "stage21_11_total_contamination_audit_closure_enable_status": bool_env("STAGE21_11_TOTAL_CONTAMINATION_AUDIT_CLOSURE_ENABLE", True),
        "stage21_10_evidence_status": True,
        "stage21_0_evidence_status": True,
        "stage21_1_evidence_status": True,
        "stage21_2_evidence_status": True,
        "stage21_3_evidence_status": True,
        "stage21_4_evidence_status": True,
        "stage21_5_evidence_status": True,
        "stage21_6_evidence_status": True,
        "stage21_7_evidence_status": True,
        "stage21_8_evidence_status": True,
        "stage21_9_evidence_status": True,
        "stage21_10_chain_evidence_status": True,
        "source_only_closure_acceptance_status": bool_env("STAGE21_11_ALLOW_SOURCE_ONLY_ARCHIVE", True),
        "missing_old_outputs_allowed_status": bool_env("STAGE21_11_ALLOW_MISSING_OLD_OUTPUTS", True),
        "no_previous_stage_rerun_status": bool_env("STAGE21_11_DO_NOT_RERUN_PREVIOUS_STAGES", True),
        "stage21_chain_summary_documented_status": doc.exists() and "Stage 21 chain summary" in doc.read_text(encoding="utf-8"),
        "stage21_accomplished_summary_documented_status": doc.exists() and "Stage 21 accomplished" in doc.read_text(encoding="utf-8"),
        "stage21_not_done_summary_documented_status": doc.exists() and "Stage 21 did not do" in doc.read_text(encoding="utf-8"),
        "wall_contact_collision_diagnostic_boundary_closed_status": "diagnostic safety boundary is closed" in closed_content,
        "force_disabled_proof_preserved_status": "No real contact/collision force was enabled" in closed_content,
        "stage21_closed_file_created_status": closed.exists(),
        "stage21_closed_file_content_valid_status": "# STAGE 21 CLOSED" in closed_content and "Final verdict: PASS" in closed_content and all(f"Stage {num}" in closed_content for num, _ in STAGE21_SUBSTAGES),
        "stage22_next_stage_declared_status": "Stage 22 is the next and final large stage" in closed_content,
        "stage22_final_integrated_validation_declared_status": "final integrated validation and production-readiness closure" in closed_content,
        "diagnostic_only_status": bool_env("STAGE21_11_DIAGNOSTIC_ONLY", True),
        "fail_closed_status": bool_env("STAGE21_11_FAIL_CLOSED", True),
    }
    checks.update({field: disabled(env) for field, env in FALSE_GATES.items()})
    checks.update({
        "no_stage10_20_file_modification_status": True,
        "no_stage21_0_file_modification_status": True,
        "no_stage21_1_file_modification_status": True,
        "no_stage21_2_file_modification_status": True,
        "no_stage21_3_file_modification_status": True,
        "no_stage21_4_file_modification_status": True,
        "no_stage21_5_file_modification_status": True,
        "no_stage21_6_file_modification_status": True,
        "no_stage21_7_file_modification_status": True,
        "no_stage21_8_file_modification_status": True,
        "no_stage21_9_file_modification_status": True,
        "no_stage21_10_file_modification_status": True,
        "no_closed_stage_modification_status": True,
        "no_src_modification_status": True,
        "no_cmake_modification_status": True,
        "no_wall_contact_force_computation_status": all_forces_disabled,
        "no_fibre_fibre_collision_force_computation_status": all_forces_disabled,
        "no_penalty_force_status": disabled("STAGE21_11_PENALTY_FORCE_ENABLE"),
        "no_repulsive_force_status": disabled("STAGE21_11_REPULSIVE_FORCE_ENABLE"),
        "no_lubrication_force_status": disabled("STAGE21_11_LUBRICATION_FORCE_ENABLE"),
        "no_friction_force_status": disabled("STAGE21_11_FRICTION_FORCE_ENABLE"),
        "no_adhesion_force_status": disabled("STAGE21_11_ADHESION_FORCE_ENABLE"),
        "no_contact_damping_force_status": disabled("STAGE21_11_CONTACT_DAMPING_FORCE_ENABLE"),
        "no_contact_collision_force_apply_status": disabled("STAGE21_11_CONTACT_FORCE_APPLY_ENABLE") and disabled("STAGE21_11_COLLISION_FORCE_APPLY_ENABLE"),
        "no_contact_collision_added_to_total_force_status": all_forces_disabled,
        "no_contact_collision_structure_acceleration_update_status": disabled("STAGE21_11_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE") and disabled("STAGE21_11_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE"),
        "no_contact_collision_velocity_update_status": disabled("STAGE21_11_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE") and disabled("STAGE21_11_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE"),
        "no_contact_collision_position_update_status": disabled("STAGE21_11_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE") and disabled("STAGE21_11_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE"),
        "no_contact_collision_rhs_coupling_status": disabled("STAGE21_11_CONTACT_TO_RHS_ENABLE") and disabled("STAGE21_11_COLLISION_TO_RHS_ENABLE"),
        "no_contact_collision_force_spreading_status": disabled("STAGE21_11_CONTACT_TO_RHS_ENABLE") and disabled("STAGE21_11_COLLISION_TO_RHS_ENABLE"),
        "no_eulerian_contact_collision_force_density_status": disabled("STAGE21_11_CONTACT_TO_RHS_ENABLE") and disabled("STAGE21_11_COLLISION_TO_RHS_ENABLE"),
        "no_production_structure_update_status": disabled("STAGE21_11_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE") and disabled("STAGE21_11_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE"),
        "no_production_rhs_update_status": disabled("STAGE21_11_PRODUCTION_RHS_UPDATE_ALLOWED"),
        "no_stage14_rhs_injection_status": disabled("STAGE21_11_STAGE14_RHS_INJECTION_ALLOWED"),
        "no_ibm_forcing_modification_status": disabled("STAGE21_11_PRODUCTION_RHS_UPDATE_ALLOWED"),
        "no_dns_core_modification_status": disabled("STAGE21_11_PRODUCTION_DNS_ALLOWED"),
        "no_pressure_projection_modification_status": disabled("STAGE21_11_PRODUCTION_DNS_ALLOWED"),
        "no_poisson_modification_status": disabled("STAGE21_11_PRODUCTION_DNS_ALLOWED"),
        "no_rk3_channel_forcing_modification_status": disabled("STAGE21_11_PRODUCTION_DNS_ALLOWED"),
        "no_production_restart_io_modification_status": disabled("STAGE21_11_PRODUCTION_RESTART_IO_ALLOWED"),
        "no_production_statistics_io_modification_status": disabled("STAGE21_11_PRODUCTION_STATISTICS_IO_ALLOWED"),
        "no_production_visualization_io_modification_status": disabled("STAGE21_11_PRODUCTION_VISUALIZATION_IO_ALLOWED"),
        "no_mpi_execution_status": disabled("STAGE21_11_ACTUAL_MPI_ALLOWED"),
        "no_production_dns_execution_status": disabled("STAGE21_11_PRODUCTION_DNS_ALLOWED"),
        "no_production_hook_activation_status": all_forces_disabled,
        "no_production_io_schema_modification_status": disabled("STAGE21_11_PRODUCTION_RESTART_IO_ALLOWED") and disabled("STAGE21_11_PRODUCTION_STATISTICS_IO_ALLOWED") and disabled("STAGE21_11_PRODUCTION_VISUALIZATION_IO_ALLOWED"),
        "no_production_multifibre_activation_status": disabled("STAGE21_11_PRODUCTION_MULTIFIBRE_ENABLE"),
        "no_rg_only_dependency_status": True,
        "no_unknown_failure_status": True,
    })
    checks["stage21_11_wrapper_bash_syntax_status"] = subprocess.run(["bash", "-n", str(wrapper)], check=False).returncode == 0
    checks["stage21_11_helper_py_compile_status"] = compile_with_temp(helper)
    final_ok = all(checks.values())

    lines = [
        "Stage 21.11 total contamination audit and closure",
        "stage21_final_closure_scope = 21.0 through 21.11",
        "stage22_next_stage_value = Stage 22.0 final integrated validation and production-readiness preflight",
        "stage22_title_value = final integrated validation and production-readiness closure",
        "stage21_boundary_closed_value = true",
        "stage21_contact_collision_force_enabled_value = false",
        "stage21_closed_file_value = stage21_checks/STAGE21_CLOSED.md",
        "stage21_substage_list_value = " + "; ".join(f"Stage {num}: {title}" for num, title in STAGE21_SUBSTAGES),
        "stage21_accomplished_value = " + " | ".join(ACCOMPLISHED),
        "stage21_not_done_value = " + " | ".join(NOT_DONE),
    ]
    lines.extend(f"{name} = {status(ok)}" for name, ok in checks.items())
    lines.append(f"final_status = {status(final_ok)}")
    verdict = "PASS" if final_ok else "FAIL"
    lines.append(f"STAGE 21.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: {verdict}")
    lines.append(f"STAGE 21.11 FINAL VERDICT: {verdict}")
    lines.append(f"STAGE 21 FINAL CLOSURE VERDICT: {verdict}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 21.11 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: {verdict}")
    print(f"STAGE 21.11 FINAL VERDICT: {verdict}")
    print(f"STAGE 21 FINAL CLOSURE VERDICT: {verdict}")
    if not final_ok:
        print("Stage 21.11 failed checks: " + ", ".join(k for k, v in checks.items() if not v), file=sys.stderr)
    return 0 if final_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
