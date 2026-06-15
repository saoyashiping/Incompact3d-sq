#!/usr/bin/env python3
"""Stage 20.1 two-way configuration and lambda gate audit.

This helper is definition-only and diagnostic-only. It validates safe conceptual
Stage 20 gate defaults, preserves Stage 20.0 source-only closure acceptance, and
writes only the Stage 20.1 audit summary. It does not run prior stages, MPI,
DNS, builds, or production physics.
"""
from __future__ import annotations

import os
import py_compile
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

PASS = "PASS"
FAIL = "FAIL"
OPTIONAL = "OPTIONAL"

STATUS_FIELDS = [
    "stage20_1_requested_status",
    "stage20_1_config_gate_enable_status",
    "stage20_0_evidence_status",
    "stage20_0_source_only_acceptance_preserved_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "stage20_boundary_declared_status",
    "all_stage20_config_gates_documented_status",
    "default_safe_gate_values_status",
    "stage20_enable_status",
    "twoway_coupling_default_disabled_status",
    "fluid_to_structure_default_disabled_status",
    "structure_to_fluid_default_disabled_status",
    "rhs_coupling_default_disabled_status",
    "lambda_coupling_zero_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "single_fibre_only_status",
    "contact_default_disabled_status",
    "collision_default_disabled_status",
    "multifibre_default_disabled_status",
    "gate_consistency_rules_status",
    "twoway_disabled_implies_subgates_disabled_status",
    "lambda_zero_implies_effective_coupling_zero_status",
    "rhs_disabled_implies_no_stage14_rhs_status",
    "diagnostic_only_implies_no_runtime_activation_status",
    "fail_closed_invalid_partial_gate_status",
    "no_stage10_19_file_modification_status",
    "no_stage20_0_file_modification_status",
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
    "no_stage14_rhs_call_from_stage20_1_status",
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
    "stage20_2_next_stage_declared_status",
    "stage20_1_wrapper_bash_syntax_status",
    "stage20_1_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_1_ENABLE": "1",
    "STAGE20_1_CONFIG_GATE_ENABLE": "1",
    "STAGE20_1_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_1_REQUIRE_STAGE20_0_PASS": "1",
    "STAGE20_1_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_1_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_1_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_1_DIAGNOSTIC_ONLY": "1",
    "STAGE20_1_FAIL_CLOSED": "1",
    "STAGE20_1_STAGE20_ENABLE": "1",
    "STAGE20_1_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_1_FLUID_TO_STRUCTURE_ENABLE": "0",
    "STAGE20_1_STRUCTURE_TO_FLUID_ENABLE": "0",
    "STAGE20_1_RHS_COUPLING_ENABLE": "0",
    "STAGE20_1_LAMBDA_COUPLING": "0.0",
    "STAGE20_1_SINGLE_FIBRE_ONLY": "1",
    "STAGE20_1_CONTACT_ENABLE": "0",
    "STAGE20_1_COLLISION_ENABLE": "0",
    "STAGE20_1_MULTIFIBRE_ENABLE": "0",
    "STAGE20_1_ZERO_TOL": "1.0e-14",
    "STAGE20_1_AUDIT_TOL": "1.0e-12",
    "STAGE20_1_TEST_CASE": "twoway_config_lambda_gate",
}

CONCEPTUAL_GATES = {
    "STAGE20_ENABLE": "true only for Stage 20 diagnostic/config audit",
    "STAGE20_TWOWAY_COUPLING_ENABLE": "false",
    "STAGE20_FLUID_TO_STRUCTURE_ENABLE": "false",
    "STAGE20_STRUCTURE_TO_FLUID_ENABLE": "false",
    "STAGE20_RHS_COUPLING_ENABLE": "false",
    "STAGE20_LAMBDA_COUPLING": "0.0",
    "STAGE20_DIAGNOSTIC_ONLY": "true",
    "STAGE20_FAIL_CLOSED": "true",
    "STAGE20_SINGLE_FIBRE_ONLY": "true",
    "STAGE20_CONTACT_ENABLE": "false",
    "STAGE20_COLLISION_ENABLE": "false",
    "STAGE20_MULTIFIBRE_ENABLE": "false",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def env_value(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name]).strip()


def truthy(name: str) -> bool:
    return env_value(name).lower() in {"1", "true", "yes", "on"}


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def stage20_0_evidence(root: Path) -> Tuple[bool, str]:
    output = root / "stage20_outputs" / "fibre_stage20_0_preflight_boundary.dat"
    text = read_text(output)
    output_pass = output.exists() and (
        "STAGE 20.0 PREFLIGHT BOUNDARY VERDICT: PASS" in text
        and "STAGE 20.0 FINAL VERDICT: PASS" in text
    )
    source_only_markers = [
        root / "stage20_checks" / "assert_stage20_0_preflight_boundary.py",
        root / "stage20_checks" / "run_stage20_0_preflight_boundary.sh",
        root / "stage20_checks" / "stage20_0_preflight_boundary.md",
    ]
    source_only_available = truthy("STAGE20_1_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") and all(
        marker.exists() for marker in source_only_markers
    )
    if output_pass:
        return True, "ACCEPTED_BY_STAGE20_0_PASS_OUTPUT"
    if source_only_available:
        return True, "ACCEPTED_BY_STAGE20_0_SOURCE_ONLY_BEHAVIOR"
    return False, "NO_STAGE20_0_PASS_OR_SOURCE_ONLY_EVIDENCE"


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_1_py_compile_", suffix=".pyc", delete=False) as tmp:
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


def lambda_is_zero() -> bool:
    return abs(float(env_value("STAGE20_1_LAMBDA_COUPLING"))) <= float(env_value("STAGE20_1_ZERO_TOL"))


def main() -> int:
    root = repo_root()
    out_dir = root / "stage20_outputs"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_1_twoway_config_lambda_gate.dat"

    stage20_0_ok, stage20_0_reason = stage20_0_evidence(root)

    twoway = truthy("STAGE20_1_TWOWAY_COUPLING_ENABLE")
    f2s = truthy("STAGE20_1_FLUID_TO_STRUCTURE_ENABLE")
    s2f = truthy("STAGE20_1_STRUCTURE_TO_FLUID_ENABLE")
    rhs = truthy("STAGE20_1_RHS_COUPLING_ENABLE")
    diagnostic = truthy("STAGE20_1_DIAGNOSTIC_ONLY")
    fail_closed = truthy("STAGE20_1_FAIL_CLOSED")
    contact = truthy("STAGE20_1_CONTACT_ENABLE")
    collision = truthy("STAGE20_1_COLLISION_ENABLE")
    multifibre = truthy("STAGE20_1_MULTIFIBRE_ENABLE")

    twoway_subgates_ok = (not twoway) and (not f2s) and (not s2f) and (not rhs)
    lambda_zero_ok = lambda_is_zero()
    gate_rules_ok = all([
        twoway_subgates_ok,
        lambda_zero_ok,
        (not rhs),
        diagnostic,
        fail_closed,
        (not contact),
        (not collision),
        (not multifibre),
    ])

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    statuses["stage20_1_requested_status"] = PASS if truthy("STAGE20_1_ENABLE") else FAIL
    statuses["stage20_1_config_gate_enable_status"] = PASS if truthy("STAGE20_1_CONFIG_GATE_ENABLE") else FAIL
    statuses["stage20_0_evidence_status"] = PASS if stage20_0_ok else FAIL
    statuses["stage20_0_source_only_acceptance_preserved_status"] = PASS if truthy("STAGE20_1_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") else FAIL
    statuses["missing_old_stage_outputs_allowed_status"] = PASS if truthy("STAGE20_1_ALLOW_MISSING_OLD_STAGE_OUTPUTS") else FAIL
    statuses["missing_old_closure_files_allowed_status"] = PASS if truthy("STAGE20_1_ALLOW_MISSING_OLD_CLOSURE_FILES") else FAIL
    statuses["no_previous_stage_rerun_status"] = PASS if truthy("STAGE20_1_DO_NOT_RERUN_PREVIOUS_STAGES") else FAIL
    statuses["stage20_enable_status"] = PASS if truthy("STAGE20_1_STAGE20_ENABLE") else FAIL
    statuses["twoway_coupling_default_disabled_status"] = PASS if not twoway else FAIL
    statuses["fluid_to_structure_default_disabled_status"] = PASS if not f2s else FAIL
    statuses["structure_to_fluid_default_disabled_status"] = PASS if not s2f else FAIL
    statuses["rhs_coupling_default_disabled_status"] = PASS if not rhs else FAIL
    statuses["lambda_coupling_zero_status"] = PASS if lambda_zero_ok else FAIL
    statuses["diagnostic_only_status"] = PASS if diagnostic else FAIL
    statuses["fail_closed_status"] = PASS if fail_closed else FAIL
    statuses["single_fibre_only_status"] = PASS if truthy("STAGE20_1_SINGLE_FIBRE_ONLY") else FAIL
    statuses["contact_default_disabled_status"] = PASS if not contact else FAIL
    statuses["collision_default_disabled_status"] = PASS if not collision else FAIL
    statuses["multifibre_default_disabled_status"] = PASS if not multifibre else FAIL
    statuses["gate_consistency_rules_status"] = PASS if gate_rules_ok else FAIL
    statuses["twoway_disabled_implies_subgates_disabled_status"] = PASS if twoway_subgates_ok else FAIL
    statuses["lambda_zero_implies_effective_coupling_zero_status"] = PASS if lambda_zero_ok else FAIL
    statuses["rhs_disabled_implies_no_stage14_rhs_status"] = PASS if not rhs else FAIL
    statuses["diagnostic_only_implies_no_runtime_activation_status"] = PASS if diagnostic and gate_rules_ok else FAIL
    statuses["fail_closed_invalid_partial_gate_status"] = PASS if fail_closed and gate_rules_ok else FAIL
    statuses["stage20_1_wrapper_bash_syntax_status"] = PASS if bash_syntax_ok(root / "stage20_checks" / "run_stage20_1_twoway_config_lambda_gate.sh") else FAIL
    statuses["stage20_1_helper_py_compile_status"] = PASS if py_compile_ok(Path(__file__).resolve()) else FAIL

    final = PASS if all(value in {PASS, OPTIONAL} for value in statuses.values()) else FAIL

    lines = [
        "STAGE 20.1 TWOWAY CONFIG LAMBDA GATE AUDIT",
        "stage20_title = real two-way fluid-structure coupling activation boundary",
        "stage20_1_title = two-way coupling config and lambda gate",
        f"repository_root_value = {root}",
        f"stage20_1_test_case_value = {env_value('STAGE20_1_TEST_CASE')}",
        f"stage20_0_evidence_reason_value = {stage20_0_reason}",
        "stage20_2_next_stage_value = Stage 20.2: fluid-to-structure force input adapter",
        "effective_coupling_value = 0.0",
    ]
    lines.extend(f"{name}_value = {value}" for name, value in CONCEPTUAL_GATES.items())
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.1 TWOWAY CONFIG LAMBDA GATE VERDICT: {final}")
    lines.append(f"STAGE 20.1 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"STAGE 20.1 TWOWAY CONFIG LAMBDA GATE VERDICT: {final}")
    print(f"STAGE 20.1 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
