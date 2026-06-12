#!/usr/bin/env python3
"""Stage 19.1 production physical structure configuration gate.

Pure-Python, configuration-only diagnostic for Stage 19.1.  It defines the
fail-closed Stage 19 physical-structure gate schema, validates safe defaults and
environment override consistency, accepts Stage 19.0 / Stage 18 closure evidence
without rerunning prior stages, and writes only a helper-local stage19_outputs
artifact.  It does not build, run MPI, run DNS, add Fortran, touch CMake, create
production X/V/A state, insert hooks, activate advance/commit paths, spread RHS
forces, modify IBM/DNS-core/projection/Poisson/RK3, or add contact/collision /
production multifibre logic.
"""
from __future__ import annotations

import argparse
import os
import py_compile
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS: List[str] = [
    "stage19_1_requested_status",
    "stage19_1_config_gate_enable_status",
    "stage19_0_evidence_status",
    "stage18_closure_evidence_status",
    "stage18_source_only_closure_acceptance_preserved_status",
    "no_stage10_18_file_modification_status",
    "no_stage19_0_file_modification_status",
    "no_closed_stage_modification_status",
    "config_schema_documented_status",
    "all_required_gates_present_status",
    "default_safe_values_status",
    "diagnostic_only_status",
    "single_fibre_only_status",
    "fail_closed_status",
    "physical_structure_default_disabled_status",
    "structure_state_default_disabled_status",
    "structure_init_default_disabled_status",
    "force_candidate_default_disabled_status",
    "advance_candidate_default_disabled_status",
    "commit_default_disabled_status",
    "hook_default_disabled_status",
    "fluid_force_input_default_disabled_status",
    "rhs_spreading_default_disabled_status",
    "stage14_rhs_injection_default_disabled_status",
    "restart_io_default_disabled_status",
    "statistics_io_default_disabled_status",
    "visualization_io_default_disabled_status",
    "contact_model_default_disabled_status",
    "fibre_fibre_collision_default_disabled_status",
    "multifibre_production_default_disabled_status",
    "diagnostic_only_consistency_status",
    "single_fibre_only_consistency_status",
    "fail_closed_consistency_status",
    "rhs_spreading_disabled_consistency_status",
    "stage14_rhs_injection_disabled_consistency_status",
    "contact_collision_disabled_consistency_status",
    "io_disabled_consistency_status",
    "stage19_1_wrapper_bash_syntax_status",
    "stage19_1_helper_py_compile_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_bending_force_runtime_application_status",
    "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage19_1_status",
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
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]

@dataclass(frozen=True)
class Gate:
    logical_name: str
    env_name: str
    default: bool

GATES: List[Gate] = [
    Gate("stage19_physical_structure_enable", "STAGE19_1_PHYSICAL_STRUCTURE_ENABLE", False),
    Gate("stage19_physical_structure_diagnostic_only", "STAGE19_1_DIAGNOSTIC_ONLY", True),
    Gate("stage19_physical_structure_single_fibre_only", "STAGE19_1_SINGLE_FIBRE_ONLY", True),
    Gate("stage19_physical_structure_fail_closed", "STAGE19_1_FAIL_CLOSED", True),
    Gate("stage19_physical_structure_state_enable", "STAGE19_1_STRUCTURE_STATE_ENABLE", False),
    Gate("stage19_physical_structure_init_enable", "STAGE19_1_STRUCTURE_INIT_ENABLE", False),
    Gate("stage19_physical_structure_force_candidate_enable", "STAGE19_1_FORCE_CANDIDATE_ENABLE", False),
    Gate("stage19_physical_structure_advance_candidate_enable", "STAGE19_1_ADVANCE_CANDIDATE_ENABLE", False),
    Gate("stage19_physical_structure_commit_enable", "STAGE19_1_COMMIT_ENABLE", False),
    Gate("stage19_physical_structure_hook_enable", "STAGE19_1_HOOK_ENABLE", False),
    Gate("stage19_fluid_force_input_enable", "STAGE19_1_FLUID_FORCE_INPUT_ENABLE", False),
    Gate("stage19_rhs_spreading_enable", "STAGE19_1_RHS_SPREADING_ENABLE", False),
    Gate("stage19_stage14_rhs_injection_enable", "STAGE19_1_STAGE14_RHS_INJECTION_ENABLE", False),
    Gate("stage19_restart_io_enable", "STAGE19_1_RESTART_IO_ENABLE", False),
    Gate("stage19_statistics_io_enable", "STAGE19_1_STATISTICS_IO_ENABLE", False),
    Gate("stage19_visualization_io_enable", "STAGE19_1_VISUALIZATION_IO_ENABLE", False),
    Gate("stage19_contact_model_enable", "STAGE19_1_CONTACT_MODEL_ENABLE", False),
    Gate("stage19_fibre_fibre_collision_enable", "STAGE19_1_FIBRE_FIBRE_COLLISION_ENABLE", False),
    Gate("stage19_multifibre_production_enable", "STAGE19_1_MULTIFIBRE_PRODUCTION_ENABLE", False),
]

ALLOWED_NEW_OR_MODIFIED = {
    "stage19_checks/run_stage19_1_physical_structure_config_gate.sh",
    "stage19_checks/assert_stage19_1_physical_structure_config_gate.py",
    "stage19_checks/stage19_1_physical_structure_config_gate.md",
    "stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat",
}

STAGE19_0_FILES = {
    "stage19_checks/run_stage19_0_preflight_boundary.sh",
    "stage19_checks/assert_stage19_0_preflight_boundary.py",
    "stage19_checks/stage19_0_preflight_boundary.md",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
}

ACCEPTED_UNTRACKED_EVIDENCE = {
    "stage17_checks/STAGE17_CLOSED.md",
    "stage18_checks/STAGE18_CLOSED.md",
    "stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
}

TRUE_VALUES = {"1", "true", "TRUE", "yes", "YES", "on", "ON", "t", "T"}
FALSE_VALUES = {"0", "false", "FALSE", "no", "NO", "off", "OFF", "f", "F"}


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def parse_bool_env(name: str, default: bool, invalid: List[str]) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    value = raw.strip()
    if value in TRUE_VALUES:
        return True
    if value in FALSE_VALUES:
        return False
    invalid.append(f"{name}={raw}")
    return default


def run_quiet(cmd: Sequence[str], cwd: Path) -> Tuple[int, str]:
    try:
        proc = subprocess.run(cmd, cwd=str(cwd), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        return proc.returncode, proc.stdout
    except OSError as exc:
        return 127, str(exc)


def git_status_entries(root: Path) -> List[Tuple[str, str]]:
    code, out = run_quiet(["git", "status", "--porcelain", "--untracked-files=all"], root)
    if code != 0:
        return []
    entries: List[Tuple[str, str]] = []
    for raw in out.splitlines():
        if not raw:
            continue
        xy = raw[:2]
        path = raw[3:] if len(raw) > 3 else raw
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        entries.append((xy, path.strip()))
    return entries


def outside_allowed_changes(root: Path) -> List[str]:
    outside: List[str] = []
    for xy, path in git_status_entries(root):
        if path in ALLOWED_NEW_OR_MODIFIED:
            continue
        if xy == "??" and path in ACCEPTED_UNTRACKED_EVIDENCE:
            continue
        outside.append(path)
    return outside


def evidence_has_final_pass(path: Path, verdict: str) -> bool:
    text = read_text(path)
    return "final_status PASS" in text and verdict in text


def stage18_closed_marker_valid(root: Path) -> bool:
    marker = root / "stage18_checks" / "STAGE18_CLOSED.md"
    text = read_text(marker)
    lowered = text.lower()
    stages = all(token in text for token in ["18.0", "18.1", "18.2", "18.3", "18.4", "18.5", "18.6", "18.7", "18.8", "18.9", "18.10", "18.11", "18.12"])
    closure = marker.is_file() and "stage 18" in lowered and any(word in lowered for word in ("closed", "closure"))
    no_contamination = all(term in text for term in ("RHS", "IBM", "DNS-core"))
    no_collision = all(term in lowered for term in ("contact", "collision", "multifibre"))
    return closure and stages and no_contamination and no_collision


def stage18_12_dat_valid(root: Path) -> bool:
    path = root / "stage18_outputs" / "fibre_stage18_12_total_contamination_audit_closure.dat"
    text = read_text(path)
    required = [
        "final_status PASS",
        "STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS",
        "STAGE 18.12 FINAL VERDICT: PASS",
        "stage18_closed_file_created_status PASS",
    ]
    return path.is_file() and all(item in text for item in required)


def stage19_0_source_acceptance_preserved(root: Path) -> bool:
    helper = read_text(root / "stage19_checks" / "assert_stage19_0_preflight_boundary.py")
    doc = read_text(root / "stage19_checks" / "stage19_0_preflight_boundary.md")
    return all(token in helper for token in (
        "stage18_closure_accepted_status",
        "prior_stage18_outputs_required_status",
        "stage18_closure_supersedes_individual_outputs_status",
        "ACCEPTED_BY_STAGE18_CLOSURE",
    )) and "must not force users to rerun Stage 18.0 through Stage 18.11" in doc


def stage19_0_evidence_valid(root: Path) -> bool:
    output = root / "stage19_outputs" / "fibre_stage19_0_preflight_boundary.dat"
    if evidence_has_final_pass(output, "STAGE 19.0 FINAL VERDICT: PASS"):
        return True
    return all((root / path).is_file() for path in STAGE19_0_FILES if path.startswith("stage19_checks/")) and stage19_0_source_acceptance_preserved(root)


def stage18_evidence_safely_accepted(root: Path) -> bool:
    return (stage18_closed_marker_valid(root) and stage18_12_dat_valid(root)) or stage19_0_source_acceptance_preserved(root)


def syntax_status(root: Path) -> Tuple[str, str]:
    wrapper = root / "stage19_checks" / "run_stage19_1_physical_structure_config_gate.sh"
    helper = root / "stage19_checks" / "assert_stage19_1_physical_structure_config_gate.py"
    bash_code, _ = run_quiet(["bash", "-n", str(wrapper)], root)
    try:
        with tempfile.TemporaryDirectory() as td:
            cfile = Path(td) / "stage19_1_helper.pyc"
            py_compile.compile(str(helper), cfile=str(cfile), doraise=True)
        py_status = "PASS"
    except py_compile.PyCompileError:
        py_status = "FAIL"
    return pass_fail(bash_code == 0), py_status


def all_required_gates_present(values: Dict[str, bool]) -> bool:
    return {gate.logical_name for gate in GATES} == set(values)


def defaults_are_safe(values: Dict[str, bool], env_overridden: bool) -> bool:
    if env_overridden:
        return True
    return all(values[gate.logical_name] is gate.default for gate in GATES)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 19.1 physical structure config gate")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    root = args.repo_root.resolve()
    output = args.output or (root / "stage19_outputs" / "fibre_stage19_1_physical_structure_config_gate.dat")
    statuses: Dict[str, str] = {key: "PASS" for key in SUMMARY_KEYS if key != "final_status"}
    failure_reasons: List[str] = []
    invalid_env: List[str] = []

    stage19_1_enable = parse_bool_env("STAGE19_1_ENABLE", True, invalid_env)
    config_gate_enable = parse_bool_env("STAGE19_1_CONFIG_GATE_ENABLE", True, invalid_env)
    require_stage19_0 = parse_bool_env("STAGE19_1_REQUIRE_STAGE19_0_PASS", True, invalid_env)
    require_stage18 = parse_bool_env("STAGE19_1_REQUIRE_STAGE18_CLOSED", True, invalid_env)

    values = {gate.logical_name: parse_bool_env(gate.env_name, gate.default, invalid_env) for gate in GATES}
    env_overridden = any(gate.env_name in os.environ for gate in GATES)

    statuses["stage19_1_requested_status"] = pass_fail(stage19_1_enable)
    statuses["stage19_1_config_gate_enable_status"] = pass_fail(config_gate_enable)
    statuses["stage19_0_evidence_status"] = pass_fail((not require_stage19_0) or stage19_0_evidence_valid(root))
    statuses["stage18_closure_evidence_status"] = pass_fail((not require_stage18) or stage18_evidence_safely_accepted(root))
    statuses["stage18_source_only_closure_acceptance_preserved_status"] = pass_fail(stage19_0_source_acceptance_preserved(root))

    changed = outside_allowed_changes(root)
    statuses["no_stage10_18_file_modification_status"] = pass_fail(not any(path.startswith(("stage10_", "stage11_", "stage12_", "stage13_", "stage14_", "stage15_", "stage16_", "stage17_", "stage18_")) for path in changed))
    statuses["no_stage19_0_file_modification_status"] = pass_fail(not any(path in STAGE19_0_FILES for path in changed))
    statuses["no_closed_stage_modification_status"] = pass_fail(statuses["no_stage10_18_file_modification_status"] == "PASS" and statuses["no_stage19_0_file_modification_status"] == "PASS")

    doc = read_text(root / "stage19_checks" / "stage19_1_physical_structure_config_gate.md")
    statuses["config_schema_documented_status"] = pass_fail(all(gate.logical_name in doc for gate in GATES))
    statuses["all_required_gates_present_status"] = pass_fail(all_required_gates_present(values))
    statuses["default_safe_values_status"] = pass_fail(defaults_are_safe(values, env_overridden))

    statuses["diagnostic_only_status"] = pass_fail(values["stage19_physical_structure_diagnostic_only"])
    statuses["single_fibre_only_status"] = pass_fail(values["stage19_physical_structure_single_fibre_only"])
    statuses["fail_closed_status"] = pass_fail(values["stage19_physical_structure_fail_closed"])
    statuses["physical_structure_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_enable"])
    statuses["structure_state_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_state_enable"])
    statuses["structure_init_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_init_enable"])
    statuses["force_candidate_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_force_candidate_enable"])
    statuses["advance_candidate_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_advance_candidate_enable"])
    statuses["commit_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_commit_enable"])
    statuses["hook_default_disabled_status"] = pass_fail(not values["stage19_physical_structure_hook_enable"])
    statuses["fluid_force_input_default_disabled_status"] = pass_fail(not values["stage19_fluid_force_input_enable"])
    statuses["rhs_spreading_default_disabled_status"] = pass_fail(not values["stage19_rhs_spreading_enable"])
    statuses["stage14_rhs_injection_default_disabled_status"] = pass_fail(not values["stage19_stage14_rhs_injection_enable"])
    statuses["restart_io_default_disabled_status"] = pass_fail(not values["stage19_restart_io_enable"])
    statuses["statistics_io_default_disabled_status"] = pass_fail(not values["stage19_statistics_io_enable"])
    statuses["visualization_io_default_disabled_status"] = pass_fail(not values["stage19_visualization_io_enable"])
    statuses["contact_model_default_disabled_status"] = pass_fail(not values["stage19_contact_model_enable"])
    statuses["fibre_fibre_collision_default_disabled_status"] = pass_fail(not values["stage19_fibre_fibre_collision_enable"])
    statuses["multifibre_production_default_disabled_status"] = pass_fail(not values["stage19_multifibre_production_enable"])

    runtime_activation = [
        "stage19_physical_structure_state_enable",
        "stage19_physical_structure_init_enable",
        "stage19_physical_structure_force_candidate_enable",
        "stage19_physical_structure_advance_candidate_enable",
        "stage19_physical_structure_commit_enable",
        "stage19_physical_structure_hook_enable",
        "stage19_fluid_force_input_enable",
        "stage19_rhs_spreading_enable",
        "stage19_stage14_rhs_injection_enable",
    ]
    physical_disabled_ok = values["stage19_physical_structure_enable"] or not any(values[name] for name in runtime_activation)
    diagnostic_ok = (not values["stage19_physical_structure_diagnostic_only"]) or not any(values[name] for name in (
        "stage19_physical_structure_commit_enable",
        "stage19_rhs_spreading_enable",
        "stage19_stage14_rhs_injection_enable",
        "stage19_restart_io_enable",
        "stage19_statistics_io_enable",
        "stage19_visualization_io_enable",
    ))
    single_fibre_ok = (not values["stage19_physical_structure_single_fibre_only"]) or not any(values[name] for name in (
        "stage19_multifibre_production_enable",
        "stage19_fibre_fibre_collision_enable",
    ))
    rhs_disabled_ok = values["stage19_rhs_spreading_enable"] or not any(values[name] for name in (
        "stage19_stage14_rhs_injection_enable",
        "stage19_fluid_force_input_enable",
    ))
    stage14_rhs_ok = not values["stage19_stage14_rhs_injection_enable"]
    contact_collision_ok = not any(values[name] for name in (
        "stage19_contact_model_enable",
        "stage19_fibre_fibre_collision_enable",
        "stage19_multifibre_production_enable",
    ))
    io_disabled_ok = not any(values[name] for name in (
        "stage19_restart_io_enable",
        "stage19_statistics_io_enable",
        "stage19_visualization_io_enable",
    ))
    fail_closed_ok = values["stage19_physical_structure_fail_closed"] and not invalid_env and physical_disabled_ok and diagnostic_ok and single_fibre_ok and rhs_disabled_ok and stage14_rhs_ok and contact_collision_ok and io_disabled_ok

    statuses["diagnostic_only_consistency_status"] = pass_fail(diagnostic_ok)
    statuses["single_fibre_only_consistency_status"] = pass_fail(single_fibre_ok)
    statuses["fail_closed_consistency_status"] = pass_fail(fail_closed_ok)
    statuses["rhs_spreading_disabled_consistency_status"] = pass_fail(rhs_disabled_ok)
    statuses["stage14_rhs_injection_disabled_consistency_status"] = pass_fail(stage14_rhs_ok)
    statuses["contact_collision_disabled_consistency_status"] = pass_fail(contact_collision_ok)
    statuses["io_disabled_consistency_status"] = pass_fail(io_disabled_ok)

    wrapper_syntax, helper_compile = syntax_status(root)
    statuses["stage19_1_wrapper_bash_syntax_status"] = wrapper_syntax
    statuses["stage19_1_helper_py_compile_status"] = helper_compile

    statuses["no_production_fortran_modification_status"] = pass_fail(not any(path.startswith("src/") and path.endswith((".f90", ".F90", ".f", ".F")) for path in changed))
    statuses["no_cmake_modification_status"] = pass_fail(not any(path == "CMakeLists.txt" or path.endswith("CMakeLists.txt") or path.startswith("cmake/") for path in changed))

    no_production_keys = [
        "no_production_structure_state_creation_status",
        "no_production_structure_buffer_creation_status",
        "no_production_structure_update_status",
        "no_production_structure_hook_status",
        "no_production_structure_advance_api_activation_status",
        "no_production_structure_commit_activation_status",
        "no_bending_force_runtime_application_status",
        "no_tension_force_runtime_application_status",
        "no_fluid_force_input_activation_status",
        "no_force_spreading_to_fluid_rhs_status",
        "no_stage14_rhs_call_from_stage19_1_status",
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
        "stage13_6_diagnostic_preserved_status",
        "stage13_no_local_subdomain_center_regression_status",
        "stage14_small_lambda_hook_status",
        "no_rg_only_dependency_status",
    ]
    for key in no_production_keys:
        statuses[key] = pass_fail(not changed)
    statuses["no_unknown_failure_status"] = pass_fail(not invalid_env)

    failing = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and statuses.get(key) != "PASS"]
    if invalid_env:
        failure_reasons.extend(f"invalid_boolean:{item}" for item in invalid_env)
    failure_reasons.extend(f"{key}={statuses.get(key, 'MISSING')}" for key in failing)
    final_status = "PASS" if not failing and not invalid_env else "FAIL"
    statuses["final_status"] = final_status

    output.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    lines.append("# Stage 19.1 physical structure config gate")
    lines.append("stage19_title production-side physical structure integration boundary")
    lines.append("stage19_1_title production physical structure config gate")
    lines.append(f"stage19_1_test_case {os.environ.get('STAGE19_1_TEST_CASE', 'production_physical_structure_config_gate')}")
    lines.append(f"stage19_1_zero_tol_value {os.environ.get('STAGE19_1_ZERO_TOL', '1.0e-14')}")
    lines.append(f"stage19_1_audit_tol_value {os.environ.get('STAGE19_1_AUDIT_TOL', '1.0e-12')}")
    lines.append("stage19_1_scope configuration-only; no production state, hook, advance, commit, RHS, IBM, DNS-core, production I/O, MPI, DNS, contact, collision, or multifibre activation")
    for gate in GATES:
        value = values[gate.logical_name]
        lines.append(f"{gate.logical_name}_value {str(value).lower()}")
        lines.append(f"{gate.logical_name}_default_value {str(gate.default).lower()}")
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses[key]}")
    if failure_reasons:
        lines.append("failure_reasons_begin")
        lines.extend(failure_reasons)
        lines.append("failure_reasons_end")
    lines.append(f"STAGE 19.1 PHYSICAL STRUCTURE CONFIG GATE VERDICT: {final_status}")
    lines.append(f"STAGE 19.1 FINAL VERDICT: {final_status}")
    output.write_text("\n".join(lines) + "\n")

    print(f"STAGE 19.1 PHYSICAL STRUCTURE CONFIG GATE VERDICT: {final_status}")
    print(f"STAGE 19.1 FINAL VERDICT: {final_status}")
    if failure_reasons:
        print("STAGE 19.1 FAILURE REASONS:")
        for reason in failure_reasons:
            print(f"  - {reason}")
    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
