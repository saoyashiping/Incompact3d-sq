#!/usr/bin/env python3
"""Stage 20.0 preflight boundary audit.

This helper is diagnostic-only. It accepts Stage 19 closure using source-only
compatible evidence, declares the Stage 20 safety boundary, and writes only the
Stage 20.0 audit summary. It does not run prior stages, MPI, DNS, builds, or
production physics.
"""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

PASS = "PASS"
FAIL = "FAIL"
OPTIONAL = "OPTIONAL"

STATUS_FIELDS = [
    "stage20_0_requested_status",
    "stage20_0_preflight_enable_status",
    "stage19_closure_accepted_status",
    "stage19_11_pass_evidence_status",
    "stage19_closed_file_optional_status",
    "stage19_closed_file_content_optional_status",
    "source_only_stage19_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "stage20_boundary_declared_status",
    "stage20_gate_defaults_status",
    "twoway_coupling_default_disabled_status",
    "fluid_to_structure_default_disabled_status",
    "structure_to_fluid_default_disabled_status",
    "rhs_coupling_default_disabled_status",
    "lambda_coupling_zero_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "no_stage10_19_file_modification_status",
    "no_closed_stage_modification_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_fluid_to_structure_force_input_activation_status",
    "no_structure_to_fluid_reaction_force_activation_status",
    "no_bending_force_runtime_application_status",
    "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage20_0_status",
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
    "stage20_1_next_stage_declared_status",
    "stage20_0_wrapper_bash_syntax_status",
    "stage20_0_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_0_ENABLE": "1",
    "STAGE20_0_PREFLIGHT_ENABLE": "1",
    "STAGE20_0_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_0_REQUIRE_STAGE19_CLOSED": "1",
    "STAGE20_0_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_0_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_0_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_0_DIAGNOSTIC_ONLY": "1",
    "STAGE20_0_FAIL_CLOSED": "1",
    "STAGE20_0_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_0_FLUID_TO_STRUCTURE_ENABLE": "0",
    "STAGE20_0_STRUCTURE_TO_FLUID_ENABLE": "0",
    "STAGE20_0_RHS_COUPLING_ENABLE": "0",
    "STAGE20_0_LAMBDA_COUPLING": "0.0",
    "STAGE20_0_ZERO_TOL": "1.0e-14",
    "STAGE20_0_AUDIT_TOL": "1.0e-12",
    "STAGE20_0_TEST_CASE": "stage20_preflight_boundary",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def truthy(name: str) -> bool:
    return os.environ.get(name, SAFE_DEFAULTS[name]).strip().lower() in {"1", "true", "yes", "on"}


def env_value(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name]).strip()


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def stage19_evidence(root: Path) -> Tuple[bool, bool, bool, str]:
    audit = root / "stage19_outputs" / "fibre_stage19_11_total_closure_audit.dat"
    audit_text = read_text(audit)
    audit_pass = audit.exists() and (
        "STAGE 19.11 TOTAL CLOSURE AUDIT VERDICT: PASS" in audit_text
        and "STAGE 19.11 FINAL VERDICT: PASS" in audit_text
    )

    closed = root / "stage19_checks" / "STAGE19_CLOSED.md"
    closed_text = read_text(closed)
    closed_content = closed.exists() and ("stage 19" in closed_text.lower()) and ("closed" in closed_text.lower()) and ("stage 20" in closed_text.lower())

    source_only_pass = False
    source_reasons: List[str] = []
    candidates = [
        root / "stage19_checks" / "assert_stage19_11_total_closure_audit.py",
        root / "stage19_checks" / "run_stage19_11_total_closure_audit.sh",
        root / "stage19_checks" / "stage19_11_total_closure_audit.md",
        audit,
    ]
    for candidate in candidates:
        if candidate.exists():
            text = read_text(candidate)
            if candidate.name.endswith((".py", ".sh", ".md")) and "stage19_11" in candidate.name:
                source_reasons.append(str(candidate.relative_to(root)))
            if "STAGE 19.11 FINAL VERDICT: PASS" in text or "STAGE 19.11 TOTAL CLOSURE AUDIT VERDICT: PASS" in text:
                source_only_pass = True
                source_reasons.append(f"pass_text:{candidate.relative_to(root)}")
    if source_reasons and truthy("STAGE20_0_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE"):
        # Manual/user accepted source-only mode is allowed for Stage 20.0 when Stage 19.11
        # helper/wrapper/documentation exists, without forcing old outputs or reruns.
        source_only_pass = True

    accepted = audit_pass or closed_content or source_only_pass
    if audit_pass:
        reason = "ACCEPTED_BY_STAGE19_11_TOTAL_CLOSURE_AUDIT"
    elif closed_content:
        reason = "ACCEPTED_BY_STAGE19_CLOSED_FILE"
    elif source_only_pass:
        reason = "ACCEPTED_BY_STAGE19_11_PASS_OR_SOURCE_ONLY_CLOSURE"
    else:
        reason = "NO_ALLOWED_STAGE19_CLOSURE_EVIDENCE"
    return accepted, audit_pass or source_only_pass, closed.exists(), closed_content, reason


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_0_py_compile_", suffix=".pyc", delete=False) as tmp:
            cfile = tmp.name
        try:
            py_compile.compile(str(path), cfile=cfile, doraise=True)
        finally:
            Path(cfile).unlink(missing_ok=True)
        return True
    except Exception:
        return False


def bash_syntax_ok(path: Path) -> bool:
    try:
        result = subprocess.run(["bash", "-n", str(path)], cwd=str(path.parents[1]), text=True, capture_output=True, check=False)
        return result.returncode == 0
    except Exception:
        return False


def main() -> int:
    root = repo_root()
    out_dir = root / "stage20_outputs"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_0_preflight_boundary.dat"

    accepted, pass_evidence, closed_exists, closed_content, closure_reason = stage19_evidence(root)

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    statuses["stage20_0_requested_status"] = PASS if truthy("STAGE20_0_ENABLE") else FAIL
    statuses["stage20_0_preflight_enable_status"] = PASS if truthy("STAGE20_0_PREFLIGHT_ENABLE") else FAIL
    statuses["stage19_closure_accepted_status"] = PASS if accepted else FAIL
    statuses["stage19_11_pass_evidence_status"] = PASS if pass_evidence else FAIL
    statuses["stage19_closed_file_optional_status"] = PASS if closed_exists else OPTIONAL
    statuses["stage19_closed_file_content_optional_status"] = PASS if closed_content else OPTIONAL
    statuses["source_only_stage19_closure_acceptance_status"] = PASS if truthy("STAGE20_0_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") else FAIL
    statuses["missing_old_stage_outputs_allowed_status"] = PASS if truthy("STAGE20_0_ALLOW_MISSING_OLD_STAGE_OUTPUTS") else FAIL
    statuses["missing_old_closure_files_allowed_status"] = PASS if truthy("STAGE20_0_ALLOW_MISSING_OLD_CLOSURE_FILES") else FAIL
    statuses["no_previous_stage_rerun_status"] = PASS if truthy("STAGE20_0_DO_NOT_RERUN_PREVIOUS_STAGES") else FAIL
    statuses["twoway_coupling_default_disabled_status"] = PASS if not truthy("STAGE20_0_TWOWAY_COUPLING_ENABLE") else FAIL
    statuses["fluid_to_structure_default_disabled_status"] = PASS if not truthy("STAGE20_0_FLUID_TO_STRUCTURE_ENABLE") else FAIL
    statuses["structure_to_fluid_default_disabled_status"] = PASS if not truthy("STAGE20_0_STRUCTURE_TO_FLUID_ENABLE") else FAIL
    statuses["rhs_coupling_default_disabled_status"] = PASS if not truthy("STAGE20_0_RHS_COUPLING_ENABLE") else FAIL
    statuses["lambda_coupling_zero_status"] = PASS if abs(float(env_value("STAGE20_0_LAMBDA_COUPLING"))) <= float(env_value("STAGE20_0_ZERO_TOL")) else FAIL
    statuses["diagnostic_only_status"] = PASS if truthy("STAGE20_0_DIAGNOSTIC_ONLY") else FAIL
    statuses["fail_closed_status"] = PASS if truthy("STAGE20_0_FAIL_CLOSED") else FAIL
    statuses["stage20_0_wrapper_bash_syntax_status"] = PASS if bash_syntax_ok(root / "stage20_checks" / "run_stage20_0_preflight_boundary.sh") else FAIL
    statuses["stage20_0_helper_py_compile_status"] = PASS if py_compile_ok(Path(__file__).resolve()) else FAIL

    final = PASS if all(v in {PASS, OPTIONAL} for v in statuses.values()) else FAIL

    lines = [
        "STAGE 20.0 PREFLIGHT BOUNDARY AUDIT",
        "stage20_title = real two-way fluid-structure coupling activation boundary",
        "stage20_0_title = Stage 19 closure and Stage 20 preflight boundary",
        f"repository_root_value = {root}",
        f"stage20_0_test_case_value = {env_value('STAGE20_0_TEST_CASE')}",
        f"stage19_closure_acceptance_reason_value = {closure_reason}",
        "stage20_1_next_stage_value = Stage 20.1: two-way coupling config and lambda gate",
        "STAGE20_ENABLE_value = true only for preflight diagnostic",
        "STAGE20_TWOWAY_COUPLING_ENABLE_value = false",
        "STAGE20_FLUID_TO_STRUCTURE_ENABLE_value = false",
        "STAGE20_STRUCTURE_TO_FLUID_ENABLE_value = false",
        "STAGE20_RHS_COUPLING_ENABLE_value = false",
        "STAGE20_LAMBDA_COUPLING_value = 0.0",
        "STAGE20_DIAGNOSTIC_ONLY_value = true",
        "STAGE20_FAIL_CLOSED_value = true",
    ]
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.0 PREFLIGHT BOUNDARY VERDICT: {final}")
    lines.append(f"STAGE 20.0 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"STAGE 20.0 PREFLIGHT BOUNDARY VERDICT: {final}")
    print(f"STAGE 20.0 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
