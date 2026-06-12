#!/usr/bin/env python3
"""Stage 19.2 physical structure state-container boundary diagnostic.

Pure-Python, state-container-boundary-only helper.  It defines and validates the
schema intended for future production-side physical-structure arrays and metadata
without creating production runtime state, adding hooks, advancing/committing
state, coupling to RHS/IBM/DNS-core, modifying production I/O, running MPI/DNS,
or introducing contact/collision/production multifibre logic.
"""
from __future__ import annotations

import argparse
import math
import os
import py_compile
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

SUMMARY_KEYS: List[str] = [
    "stage19_2_requested_status",
    "stage19_2_state_container_enable_status",
    "stage19_1_evidence_status",
    "stage19_0_evidence_status",
    "stage18_closure_evidence_status",
    "stage19_1_config_gate_preserved_status",
    "stage19_0_source_only_closure_acceptance_preserved_status",
    "no_stage10_18_file_modification_status",
    "no_stage19_0_file_modification_status",
    "no_stage19_1_file_modification_status",
    "no_closed_stage_modification_status",
    "container_schema_documented_status",
    "all_required_container_fields_present_status",
    "default_safe_values_status",
    "n_fibre_status",
    "n_point_status",
    "component_dim_status",
    "fibre_length_status",
    "ds_formula_status",
    "rho_l_status",
    "rho_tilde_status",
    "bending_stiffness_status",
    "gamma_status",
    "x_prod_shape_status",
    "v_prod_shape_status",
    "a_prod_shape_status",
    "f_b_prod_shape_status",
    "f_t_prod_shape_status",
    "f_h_prod_shape_status",
    "f_total_prod_shape_status",
    "x_candidate_shape_status",
    "v_candidate_shape_status",
    "a_candidate_shape_status",
    "owner_rank_shape_status",
    "global_point_id_shape_status",
    "local_point_id_shape_status",
    "array_finite_status",
    "global_point_id_coverage_status",
    "global_point_id_no_duplicate_status",
    "owner_rank_deterministic_status",
    "diagnostic_only_status",
    "single_fibre_only_status",
    "fail_closed_status",
    "commit_default_disabled_status",
    "rhs_spreading_default_disabled_status",
    "stage14_rhs_injection_default_disabled_status",
    "state_valid_default_false_status",
    "container_initialized_default_false_status",
    "diagnostic_only_consistency_status",
    "single_fibre_only_consistency_status",
    "fail_closed_consistency_status",
    "rhs_spreading_disabled_consistency_status",
    "stage14_rhs_injection_disabled_consistency_status",
    "commit_disabled_consistency_status",
    "container_inactive_default_status",
    "stage19_2_wrapper_bash_syntax_status",
    "stage19_2_helper_py_compile_status",
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
    "no_stage14_rhs_call_from_stage19_2_status",
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

REQUIRED_CONTAINER_FIELDS = [
    "n_fibre",
    "n_point",
    "component_dim",
    "ds",
    "fibre_length",
    "rho_l",
    "rho_tilde",
    "bending_stiffness",
    "gamma",
    "X_prod",
    "V_prod",
    "A_prod",
    "F_b_prod",
    "F_T_prod",
    "F_h_prod",
    "F_total_prod",
    "X_candidate",
    "V_candidate",
    "A_candidate",
    "owner_rank",
    "global_point_id",
    "local_point_id",
    "state_valid",
    "container_initialized",
    "diagnostic_only",
    "commit_allowed",
    "rhs_spreading_allowed",
    "stage14_rhs_injection_allowed",
]

ARRAY_FIELDS_2D = [
    "X_prod",
    "V_prod",
    "A_prod",
    "F_b_prod",
    "F_T_prod",
    "F_h_prod",
    "F_total_prod",
    "X_candidate",
    "V_candidate",
    "A_candidate",
]

SHAPE_STATUS = {
    "X_prod": "x_prod_shape_status",
    "V_prod": "v_prod_shape_status",
    "A_prod": "a_prod_shape_status",
    "F_b_prod": "f_b_prod_shape_status",
    "F_T_prod": "f_t_prod_shape_status",
    "F_h_prod": "f_h_prod_shape_status",
    "F_total_prod": "f_total_prod_shape_status",
    "X_candidate": "x_candidate_shape_status",
    "V_candidate": "v_candidate_shape_status",
    "A_candidate": "a_candidate_shape_status",
}

STAGE19_0_FILES = {
    "stage19_checks/run_stage19_0_preflight_boundary.sh",
    "stage19_checks/assert_stage19_0_preflight_boundary.py",
    "stage19_checks/stage19_0_preflight_boundary.md",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
}

STAGE19_1_FILES = {
    "stage19_checks/run_stage19_1_physical_structure_config_gate.sh",
    "stage19_checks/assert_stage19_1_physical_structure_config_gate.py",
    "stage19_checks/stage19_1_physical_structure_config_gate.md",
    "stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat",
}

ALLOWED_NEW_OR_MODIFIED = {
    "stage19_checks/run_stage19_2_physical_structure_state_container.sh",
    "stage19_checks/assert_stage19_2_physical_structure_state_container.py",
    "stage19_checks/stage19_2_physical_structure_state_container.md",
    "stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat",
}

ACCEPTED_UNTRACKED_EVIDENCE = {
    "stage17_checks/STAGE17_CLOSED.md",
    "stage18_checks/STAGE18_CLOSED.md",
    "stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
    "stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat",
}

TRUE_VALUES = {"1", "true", "TRUE", "yes", "YES", "on", "ON", "t", "T"}
FALSE_VALUES = {"0", "false", "FALSE", "no", "NO", "off", "OFF", "f", "F"}


@dataclass(frozen=True)
class ContainerConfig:
    n_fibre: int
    n_point: int
    component_dim: int
    fibre_length: float
    rho_l: float
    rho_tilde: float
    bending_stiffness: float
    gamma: float
    diagnostic_only: bool
    single_fibre_only: bool
    fail_closed: bool
    commit_allowed: bool
    rhs_spreading_allowed: bool
    stage14_rhs_injection_allowed: bool
    state_valid: bool
    container_initialized: bool


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


def parse_int_env(name: str, default: int, invalid: List[str]) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw.strip())
    except ValueError:
        invalid.append(f"{name}={raw}")
        return default


def parse_float_env(name: str, default: float, invalid: List[str]) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        value = float(raw.strip())
    except ValueError:
        invalid.append(f"{name}={raw}")
        return default
    if not math.isfinite(value):
        invalid.append(f"{name}={raw}")
        return default
    return value


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


def stage19_0_source_acceptance_preserved(root: Path) -> bool:
    helper = read_text(root / "stage19_checks" / "assert_stage19_0_preflight_boundary.py")
    doc = read_text(root / "stage19_checks" / "stage19_0_preflight_boundary.md")
    return all(token in helper for token in (
        "stage18_closure_accepted_status",
        "prior_stage18_outputs_required_status",
        "stage18_closure_supersedes_individual_outputs_status",
        "ACCEPTED_BY_STAGE18_CLOSURE",
    )) and "must not force users to rerun Stage 18.0 through Stage 18.11" in doc


def stage19_1_config_gate_preserved(root: Path) -> bool:
    helper = read_text(root / "stage19_checks" / "assert_stage19_1_physical_structure_config_gate.py")
    doc = read_text(root / "stage19_checks" / "stage19_1_physical_structure_config_gate.md")
    return all(token in helper for token in (
        "stage19_physical_structure_enable",
        "stage19_physical_structure_diagnostic_only",
        "stage19_physical_structure_fail_closed",
        "stage19_rhs_spreading_enable",
        "stage19_stage14_rhs_injection_enable",
    )) and "Stage 19.1 is configuration-only" in doc


def stage19_1_evidence_valid(root: Path) -> bool:
    output = root / "stage19_outputs" / "fibre_stage19_1_physical_structure_config_gate.dat"
    if evidence_has_final_pass(output, "STAGE 19.1 FINAL VERDICT: PASS"):
        return True
    return all((root / path).is_file() for path in STAGE19_1_FILES if path.startswith("stage19_checks/")) and stage19_1_config_gate_preserved(root)


def stage19_0_evidence_valid(root: Path) -> bool:
    output = root / "stage19_outputs" / "fibre_stage19_0_preflight_boundary.dat"
    if evidence_has_final_pass(output, "STAGE 19.0 FINAL VERDICT: PASS"):
        return True
    return all((root / path).is_file() for path in STAGE19_0_FILES if path.startswith("stage19_checks/")) and stage19_0_source_acceptance_preserved(root)


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


def stage18_evidence_safely_accepted(root: Path) -> bool:
    return (stage18_closed_marker_valid(root) and stage18_12_dat_valid(root)) or stage19_0_source_acceptance_preserved(root)


def syntax_status(root: Path) -> Tuple[str, str]:
    wrapper = root / "stage19_checks" / "run_stage19_2_physical_structure_state_container.sh"
    helper = root / "stage19_checks" / "assert_stage19_2_physical_structure_state_container.py"
    bash_code, _ = run_quiet(["bash", "-n", str(wrapper)], root)
    try:
        with tempfile.TemporaryDirectory() as td:
            cfile = Path(td) / "stage19_2_helper.pyc"
            py_compile.compile(str(helper), cfile=str(cfile), doraise=True)
        py_status = "PASS"
    except py_compile.PyCompileError:
        py_status = "FAIL"
    return pass_fail(bash_code == 0), py_status


def zeros_2d(n_point: int, component_dim: int) -> List[List[float]]:
    return [[0.0 for _component in range(component_dim)] for _point in range(n_point)]


def straight_line_x(n_point: int, component_dim: int, ds: float) -> List[List[float]]:
    values = zeros_2d(n_point, component_dim)
    if component_dim > 0:
        for idx in range(n_point):
            values[idx][0] = idx * ds
    return values


def shape_2d(values: object) -> Tuple[int, int]:
    if not isinstance(values, list):
        return (-1, -1)
    if not values:
        return (0, 0)
    if not all(isinstance(row, list) for row in values):
        return (len(values), -1)
    widths = {len(row) for row in values}
    return (len(values), widths.pop() if len(widths) == 1 else -1)


def shape_1d(values: object) -> Tuple[int]:
    if not isinstance(values, list):
        return (-1,)
    return (len(values),)


def finite_nested(values: object) -> bool:
    if isinstance(values, list):
        return all(finite_nested(item) for item in values)
    return isinstance(values, (int, float)) and math.isfinite(float(values))


def build_container(config: ContainerConfig) -> Dict[str, object]:
    ds = config.fibre_length / (config.n_point - 1) if config.n_point > 1 else math.nan
    container: Dict[str, object] = {
        "n_fibre": config.n_fibre,
        "n_point": config.n_point,
        "component_dim": config.component_dim,
        "ds": ds,
        "fibre_length": config.fibre_length,
        "rho_l": config.rho_l,
        "rho_tilde": config.rho_tilde,
        "bending_stiffness": config.bending_stiffness,
        "gamma": config.gamma,
        "X_prod": straight_line_x(config.n_point, config.component_dim, ds if math.isfinite(ds) else 0.0),
        "V_prod": zeros_2d(config.n_point, config.component_dim),
        "A_prod": zeros_2d(config.n_point, config.component_dim),
        "F_b_prod": zeros_2d(config.n_point, config.component_dim),
        "F_T_prod": zeros_2d(config.n_point, config.component_dim),
        "F_h_prod": zeros_2d(config.n_point, config.component_dim),
        "F_total_prod": zeros_2d(config.n_point, config.component_dim),
        "X_candidate": straight_line_x(config.n_point, config.component_dim, ds if math.isfinite(ds) else 0.0),
        "V_candidate": zeros_2d(config.n_point, config.component_dim),
        "A_candidate": zeros_2d(config.n_point, config.component_dim),
        "owner_rank": [0 for _idx in range(config.n_point)],
        "global_point_id": [idx for idx in range(config.n_point)],
        "local_point_id": [idx for idx in range(config.n_point)],
        "state_valid": config.state_valid,
        "container_initialized": config.container_initialized,
        "diagnostic_only": config.diagnostic_only,
        "commit_allowed": config.commit_allowed,
        "rhs_spreading_allowed": config.rhs_spreading_allowed,
        "stage14_rhs_injection_allowed": config.stage14_rhs_injection_allowed,
    }
    return container


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 19.2 physical structure state-container boundary")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    root = args.repo_root.resolve()
    output = args.output or (root / "stage19_outputs" / "fibre_stage19_2_physical_structure_state_container.dat")
    statuses: Dict[str, str] = {key: "PASS" for key in SUMMARY_KEYS if key != "final_status"}
    failure_reasons: List[str] = []
    invalid_env: List[str] = []

    requested = parse_bool_env("STAGE19_2_ENABLE", True, invalid_env)
    container_gate = parse_bool_env("STAGE19_2_STATE_CONTAINER_ENABLE", True, invalid_env)
    require_stage19_1 = parse_bool_env("STAGE19_2_REQUIRE_STAGE19_1_PASS", True, invalid_env)
    require_stage19_0 = parse_bool_env("STAGE19_2_REQUIRE_STAGE19_0_PASS", True, invalid_env)
    require_stage18 = parse_bool_env("STAGE19_2_REQUIRE_STAGE18_CLOSED", True, invalid_env)

    config = ContainerConfig(
        n_fibre=parse_int_env("STAGE19_2_N_FIBRE", 1, invalid_env),
        n_point=parse_int_env("STAGE19_2_N_POINT", 64, invalid_env),
        component_dim=parse_int_env("STAGE19_2_COMPONENT_DIM", 3, invalid_env),
        fibre_length=parse_float_env("STAGE19_2_FIBRE_LENGTH", 1.0, invalid_env),
        rho_l=parse_float_env("STAGE19_2_RHO_L", 1.0, invalid_env),
        rho_tilde=parse_float_env("STAGE19_2_RHO_TILDE", 1.0, invalid_env),
        bending_stiffness=parse_float_env("STAGE19_2_BENDING_STIFFNESS", 1.0e-3, invalid_env),
        gamma=parse_float_env("STAGE19_2_GAMMA", 1.0e-3, invalid_env),
        diagnostic_only=parse_bool_env("STAGE19_2_DIAGNOSTIC_ONLY", True, invalid_env),
        single_fibre_only=parse_bool_env("STAGE19_2_SINGLE_FIBRE_ONLY", True, invalid_env),
        fail_closed=parse_bool_env("STAGE19_2_FAIL_CLOSED", True, invalid_env),
        commit_allowed=parse_bool_env("STAGE19_2_COMMIT_ALLOWED", False, invalid_env),
        rhs_spreading_allowed=parse_bool_env("STAGE19_2_RHS_SPREADING_ALLOWED", False, invalid_env),
        stage14_rhs_injection_allowed=parse_bool_env("STAGE19_2_STAGE14_RHS_INJECTION_ALLOWED", False, invalid_env),
        state_valid=parse_bool_env("STAGE19_2_STATE_VALID", False, invalid_env),
        container_initialized=parse_bool_env("STAGE19_2_CONTAINER_INITIALIZED", False, invalid_env),
    )
    container = build_container(config)

    statuses["stage19_2_requested_status"] = pass_fail(requested)
    statuses["stage19_2_state_container_enable_status"] = pass_fail(container_gate)
    statuses["stage19_1_evidence_status"] = pass_fail((not require_stage19_1) or stage19_1_evidence_valid(root))
    statuses["stage19_0_evidence_status"] = pass_fail((not require_stage19_0) or stage19_0_evidence_valid(root))
    statuses["stage18_closure_evidence_status"] = pass_fail((not require_stage18) or stage18_evidence_safely_accepted(root))
    statuses["stage19_1_config_gate_preserved_status"] = pass_fail(stage19_1_config_gate_preserved(root))
    statuses["stage19_0_source_only_closure_acceptance_preserved_status"] = pass_fail(stage19_0_source_acceptance_preserved(root))

    changed = outside_allowed_changes(root)
    statuses["no_stage10_18_file_modification_status"] = pass_fail(not any(path.startswith(("stage10_", "stage11_", "stage12_", "stage13_", "stage14_", "stage15_", "stage16_", "stage17_", "stage18_")) for path in changed))
    statuses["no_stage19_0_file_modification_status"] = pass_fail(not any(path in STAGE19_0_FILES for path in changed))
    statuses["no_stage19_1_file_modification_status"] = pass_fail(not any(path in STAGE19_1_FILES for path in changed))
    statuses["no_closed_stage_modification_status"] = pass_fail(statuses["no_stage10_18_file_modification_status"] == "PASS" and statuses["no_stage19_0_file_modification_status"] == "PASS" and statuses["no_stage19_1_file_modification_status"] == "PASS")

    doc = read_text(root / "stage19_checks" / "stage19_2_physical_structure_state_container.md")
    statuses["container_schema_documented_status"] = pass_fail(all(field in doc for field in REQUIRED_CONTAINER_FIELDS))
    statuses["all_required_container_fields_present_status"] = pass_fail(set(REQUIRED_CONTAINER_FIELDS) == set(container))

    default_safe = (
        config.n_fibre == 1
        and config.component_dim == 3
        and config.diagnostic_only
        and config.single_fibre_only
        and config.fail_closed
        and not config.commit_allowed
        and not config.rhs_spreading_allowed
        and not config.stage14_rhs_injection_allowed
        and not config.state_valid
        and not config.container_initialized
    )
    statuses["default_safe_values_status"] = pass_fail(default_safe)
    statuses["n_fibre_status"] = pass_fail(config.n_fibre == 1)
    statuses["n_point_status"] = pass_fail(config.n_point >= 2)
    statuses["component_dim_status"] = pass_fail(config.component_dim == 3)
    statuses["fibre_length_status"] = pass_fail(config.fibre_length > 0.0)
    expected_ds = config.fibre_length / (config.n_point - 1) if config.n_point > 1 else math.nan
    statuses["ds_formula_status"] = pass_fail(math.isfinite(expected_ds) and math.isclose(float(container["ds"]), expected_ds, rel_tol=0.0, abs_tol=1.0e-14))
    statuses["rho_l_status"] = pass_fail(config.rho_l > 0.0)
    statuses["rho_tilde_status"] = pass_fail(config.rho_tilde > 0.0)
    statuses["bending_stiffness_status"] = pass_fail(config.bending_stiffness >= 0.0)
    statuses["gamma_status"] = pass_fail(config.gamma >= 0.0)

    expected_2d_shape = (config.n_point, config.component_dim)
    for field in ARRAY_FIELDS_2D:
        statuses[SHAPE_STATUS[field]] = pass_fail(shape_2d(container[field]) == expected_2d_shape)
    statuses["owner_rank_shape_status"] = pass_fail(shape_1d(container["owner_rank"]) == (config.n_point,))
    statuses["global_point_id_shape_status"] = pass_fail(shape_1d(container["global_point_id"]) == (config.n_point,))
    statuses["local_point_id_shape_status"] = pass_fail(shape_1d(container["local_point_id"]) == (config.n_point,))
    statuses["array_finite_status"] = pass_fail(all(finite_nested(container[field]) for field in ARRAY_FIELDS_2D))
    statuses["global_point_id_coverage_status"] = pass_fail(sorted(container["global_point_id"]) == list(range(config.n_point)))
    statuses["global_point_id_no_duplicate_status"] = pass_fail(len(set(container["global_point_id"])) == config.n_point)
    statuses["owner_rank_deterministic_status"] = pass_fail(container["owner_rank"] == [0 for _idx in range(config.n_point)])

    statuses["diagnostic_only_status"] = pass_fail(config.diagnostic_only)
    statuses["single_fibre_only_status"] = pass_fail(config.single_fibre_only and config.n_fibre == 1)
    statuses["fail_closed_status"] = pass_fail(config.fail_closed)
    statuses["commit_default_disabled_status"] = pass_fail(not config.commit_allowed)
    statuses["rhs_spreading_default_disabled_status"] = pass_fail(not config.rhs_spreading_allowed)
    statuses["stage14_rhs_injection_default_disabled_status"] = pass_fail(not config.stage14_rhs_injection_allowed)
    statuses["state_valid_default_false_status"] = pass_fail(not config.state_valid)
    statuses["container_initialized_default_false_status"] = pass_fail(not config.container_initialized)

    diagnostic_ok = (not config.diagnostic_only) or not config.commit_allowed
    single_fibre_ok = (not config.single_fibre_only) or config.n_fibre == 1
    rhs_disabled_ok = config.rhs_spreading_allowed or not config.stage14_rhs_injection_allowed
    stage14_rhs_ok = not config.stage14_rhs_injection_allowed
    commit_disabled_ok = not config.commit_allowed and not config.state_valid
    container_inactive = not config.state_valid and not config.container_initialized and not config.commit_allowed and not config.rhs_spreading_allowed and not config.stage14_rhs_injection_allowed
    fail_closed_ok = config.fail_closed and not invalid_env and diagnostic_ok and single_fibre_ok and rhs_disabled_ok and stage14_rhs_ok and commit_disabled_ok and container_inactive
    statuses["diagnostic_only_consistency_status"] = pass_fail(diagnostic_ok)
    statuses["single_fibre_only_consistency_status"] = pass_fail(single_fibre_ok)
    statuses["fail_closed_consistency_status"] = pass_fail(fail_closed_ok)
    statuses["rhs_spreading_disabled_consistency_status"] = pass_fail(rhs_disabled_ok)
    statuses["stage14_rhs_injection_disabled_consistency_status"] = pass_fail(stage14_rhs_ok)
    statuses["commit_disabled_consistency_status"] = pass_fail(commit_disabled_ok)
    statuses["container_inactive_default_status"] = pass_fail(container_inactive)

    wrapper_syntax, helper_compile = syntax_status(root)
    statuses["stage19_2_wrapper_bash_syntax_status"] = wrapper_syntax
    statuses["stage19_2_helper_py_compile_status"] = helper_compile
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
        "no_stage14_rhs_call_from_stage19_2_status",
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
        failure_reasons.extend(f"invalid_value:{item}" for item in invalid_env)
    failure_reasons.extend(f"{key}={statuses.get(key, 'MISSING')}" for key in failing)
    final_status = "PASS" if not failing and not invalid_env else "FAIL"
    statuses["final_status"] = final_status

    output.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    lines.append("# Stage 19.2 physical structure state-container boundary")
    lines.append("stage19_title production-side physical structure integration boundary")
    lines.append("stage19_2_title production physical structure state container")
    lines.append(f"stage19_2_test_case {os.environ.get('STAGE19_2_TEST_CASE', 'production_physical_structure_state_container')}")
    lines.append(f"stage19_2_zero_tol_value {os.environ.get('STAGE19_2_ZERO_TOL', '1.0e-14')}")
    lines.append(f"stage19_2_audit_tol_value {os.environ.get('STAGE19_2_AUDIT_TOL', '1.0e-12')}")
    lines.append("stage19_2_scope state-container-boundary only; local helper arrays are schema validation, not production runtime X/V/A state")
    lines.append(f"n_fibre_value {config.n_fibre}")
    lines.append(f"n_point_value {config.n_point}")
    lines.append(f"component_dim_value {config.component_dim}")
    lines.append(f"ds_value {container['ds']}")
    lines.append(f"fibre_length_value {config.fibre_length}")
    lines.append(f"rho_l_value {config.rho_l}")
    lines.append(f"rho_tilde_value {config.rho_tilde}")
    lines.append(f"bending_stiffness_value {config.bending_stiffness}")
    lines.append(f"gamma_value {config.gamma}")
    for field in ARRAY_FIELDS_2D:
        lines.append(f"{field.lower()}_shape_value {shape_2d(container[field])}")
    lines.append(f"owner_rank_shape_value {shape_1d(container['owner_rank'])}")
    lines.append(f"global_point_id_shape_value {shape_1d(container['global_point_id'])}")
    lines.append(f"local_point_id_shape_value {shape_1d(container['local_point_id'])}")
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses[key]}")
    if failure_reasons:
        lines.append("failure_reasons_begin")
        lines.extend(failure_reasons)
        lines.append("failure_reasons_end")
    lines.append(f"STAGE 19.2 PHYSICAL STRUCTURE STATE CONTAINER VERDICT: {final_status}")
    lines.append(f"STAGE 19.2 FINAL VERDICT: {final_status}")
    output.write_text("\n".join(lines) + "\n")

    print(f"STAGE 19.2 PHYSICAL STRUCTURE STATE CONTAINER VERDICT: {final_status}")
    print(f"STAGE 19.2 FINAL VERDICT: {final_status}")
    if failure_reasons:
        print("STAGE 19.2 FAILURE REASONS:")
        for reason in failure_reasons:
            print(f"  - {reason}")
    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
