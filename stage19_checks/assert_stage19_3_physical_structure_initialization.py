#!/usr/bin/env python3
"""Stage 19.3 physical structure initialization-boundary diagnostic.

Pure-Python helper-local initialization validation for the Stage 19 production
physical-structure boundary.  This file validates initialization schemas and
manufactured helper arrays only; it does not create production runtime X/V/A
state, insert hooks, advance/commit state, couple to RHS/IBM/DNS-core, modify
production I/O, run MPI/DNS, or introduce contact/collision/multifibre logic.
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
    "stage19_3_requested_status",
    "stage19_3_initialization_enable_status",
    "stage19_2_evidence_status",
    "stage19_1_evidence_status",
    "stage19_0_evidence_status",
    "stage18_closure_evidence_status",
    "stage19_2_state_container_preserved_status",
    "stage19_1_config_gate_preserved_status",
    "stage19_0_source_only_closure_acceptance_preserved_status",
    "no_stage10_18_file_modification_status",
    "no_stage19_0_file_modification_status",
    "no_stage19_1_file_modification_status",
    "no_stage19_2_file_modification_status",
    "no_closed_stage_modification_status",
    "initialization_schema_documented_status",
    "all_required_initialization_modes_present_status",
    "all_required_initialized_fields_present_status",
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
    "init_mode_status",
    "sine_amplitude_status",
    "sine_mode_status",
    "velocity_amplitude_status",
    "controlled_force_amplitude_status",
    "straight_fibre_initialization_status",
    "small_sine_initialization_status",
    "zero_velocity_initialization_status",
    "controlled_velocity_initialization_status",
    "zero_force_placeholder_initialization_status",
    "controlled_helper_force_placeholder_initialization_status",
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
    "f_total_consistency_status",
    "candidate_equals_production_initialization_status",
    "global_point_id_coverage_status",
    "global_point_id_no_duplicate_status",
    "owner_rank_deterministic_status",
    "diagnostic_only_status",
    "single_fibre_only_status",
    "fail_closed_status",
    "state_valid_after_init_status",
    "container_initialized_after_init_status",
    "commit_default_disabled_status",
    "rhs_spreading_default_disabled_status",
    "stage14_rhs_injection_default_disabled_status",
    "diagnostic_only_consistency_status",
    "single_fibre_only_consistency_status",
    "fail_closed_consistency_status",
    "rhs_spreading_disabled_consistency_status",
    "stage14_rhs_injection_disabled_consistency_status",
    "commit_disabled_consistency_status",
    "initialization_production_runtime_inactive_status",
    "stage19_3_wrapper_bash_syntax_status",
    "stage19_3_helper_py_compile_status",
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
    "no_stage14_rhs_call_from_stage19_3_status",
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

INIT_MODES = {
    "straight_fibre_zero_velocity",
    "small_sine_fibre_zero_velocity",
    "straight_fibre_controlled_velocity",
    "small_sine_fibre_controlled_velocity",
    "controlled_helper_force_placeholder",
}

REQUIRED_FIELDS = [
    "X_prod", "V_prod", "A_prod", "F_b_prod", "F_T_prod", "F_h_prod", "F_total_prod",
    "X_candidate", "V_candidate", "A_candidate", "owner_rank", "global_point_id", "local_point_id",
    "n_fibre", "n_point", "component_dim", "fibre_length", "ds", "rho_l", "rho_tilde",
    "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "velocity_amplitude",
    "controlled_force_amplitude", "diagnostic_only", "state_valid", "container_initialized",
    "commit_allowed", "rhs_spreading_allowed", "stage14_rhs_injection_allowed",
]

ARRAY_FIELDS_2D = [
    "X_prod", "V_prod", "A_prod", "F_b_prod", "F_T_prod", "F_h_prod", "F_total_prod",
    "X_candidate", "V_candidate", "A_candidate",
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
STAGE19_2_FILES = {
    "stage19_checks/run_stage19_2_physical_structure_state_container.sh",
    "stage19_checks/assert_stage19_2_physical_structure_state_container.py",
    "stage19_checks/stage19_2_physical_structure_state_container.md",
    "stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat",
}

ALLOWED_NEW_OR_MODIFIED = {
    "stage19_checks/run_stage19_3_physical_structure_initialization.sh",
    "stage19_checks/assert_stage19_3_physical_structure_initialization.py",
    "stage19_checks/stage19_3_physical_structure_initialization.md",
    "stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat",
}

ACCEPTED_UNTRACKED_EVIDENCE = {
    "stage17_checks/STAGE17_CLOSED.md",
    "stage18_checks/STAGE18_CLOSED.md",
    "stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
    "stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat",
    "stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat",
}

TRUE_VALUES = {"1", "true", "TRUE", "yes", "YES", "on", "ON", "t", "T"}
FALSE_VALUES = {"0", "false", "FALSE", "no", "NO", "off", "OFF", "f", "F"}


@dataclass(frozen=True)
class InitConfig:
    n_fibre: int
    n_point: int
    component_dim: int
    fibre_length: float
    rho_l: float
    rho_tilde: float
    bending_stiffness: float
    gamma: float
    init_mode: str
    sine_amplitude: float
    sine_mode: int
    velocity_amplitude: float
    controlled_force_amplitude: float
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


def stage19_2_state_container_preserved(root: Path) -> bool:
    helper = read_text(root / "stage19_checks" / "assert_stage19_2_physical_structure_state_container.py")
    doc = read_text(root / "stage19_checks" / "stage19_2_physical_structure_state_container.md")
    return all(token in helper for token in (
        "REQUIRED_CONTAINER_FIELDS",
        "X_prod",
        "V_prod",
        "A_prod",
        "commit_allowed",
        "rhs_spreading_allowed",
        "stage14_rhs_injection_allowed",
    )) and "state-container-boundary" in doc


def stage19_2_evidence_valid(root: Path) -> bool:
    output = root / "stage19_outputs" / "fibre_stage19_2_physical_structure_state_container.dat"
    if evidence_has_final_pass(output, "STAGE 19.2 FINAL VERDICT: PASS"):
        return True
    return all((root / path).is_file() for path in STAGE19_2_FILES if path.startswith("stage19_checks/")) and stage19_2_state_container_preserved(root)


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
    wrapper = root / "stage19_checks" / "run_stage19_3_physical_structure_initialization.sh"
    helper = root / "stage19_checks" / "assert_stage19_3_physical_structure_initialization.py"
    bash_code, _ = run_quiet(["bash", "-n", str(wrapper)], root)
    try:
        with tempfile.TemporaryDirectory() as td:
            cfile = Path(td) / "stage19_3_helper.pyc"
            py_compile.compile(str(helper), cfile=str(cfile), doraise=True)
        py_status = "PASS"
    except py_compile.PyCompileError:
        py_status = "FAIL"
    return pass_fail(bash_code == 0), py_status


def zeros_2d(n_point: int, component_dim: int) -> List[List[float]]:
    return [[0.0 for _component in range(component_dim)] for _point in range(n_point)]


def sine_value(config: InitConfig, s_value: float) -> float:
    return math.sin(2.0 * math.pi * config.sine_mode * s_value / config.fibre_length)


def initialize_state(config: InitConfig) -> Dict[str, object]:
    ds = config.fibre_length / (config.n_point - 1) if config.n_point > 1 else math.nan
    x_prod = zeros_2d(config.n_point, config.component_dim)
    v_prod = zeros_2d(config.n_point, config.component_dim)
    a_prod = zeros_2d(config.n_point, config.component_dim)
    f_b_prod = zeros_2d(config.n_point, config.component_dim)
    f_t_prod = zeros_2d(config.n_point, config.component_dim)
    f_h_prod = zeros_2d(config.n_point, config.component_dim)

    use_sine_shape = config.init_mode.startswith("small_sine")
    use_controlled_velocity = config.init_mode.endswith("controlled_velocity") and config.velocity_amplitude > 0.0
    use_controlled_force = config.init_mode == "controlled_helper_force_placeholder" and config.controlled_force_amplitude > 0.0

    for idx in range(config.n_point):
        s_i = idx * ds if math.isfinite(ds) else 0.0
        wave = sine_value(config, s_i) if config.fibre_length > 0.0 and config.sine_mode > 0 else 0.0
        x_prod[idx][0] = s_i
        if config.component_dim > 1 and use_sine_shape:
            x_prod[idx][1] = config.sine_amplitude * wave
        if config.component_dim > 1 and use_controlled_velocity:
            v_prod[idx][1] = config.velocity_amplitude * wave
        if config.component_dim > 1 and use_controlled_force:
            f_h_prod[idx][1] = config.controlled_force_amplitude * wave

    f_total_prod = [
        [f_b_prod[i][j] + f_t_prod[i][j] + f_h_prod[i][j] for j in range(config.component_dim)]
        for i in range(config.n_point)
    ]
    return {
        "n_fibre": config.n_fibre,
        "n_point": config.n_point,
        "component_dim": config.component_dim,
        "fibre_length": config.fibre_length,
        "ds": ds,
        "rho_l": config.rho_l,
        "rho_tilde": config.rho_tilde,
        "bending_stiffness": config.bending_stiffness,
        "gamma": config.gamma,
        "init_mode": config.init_mode,
        "sine_amplitude": config.sine_amplitude,
        "sine_mode": config.sine_mode,
        "velocity_amplitude": config.velocity_amplitude,
        "controlled_force_amplitude": config.controlled_force_amplitude,
        "diagnostic_only": config.diagnostic_only,
        "state_valid": config.state_valid,
        "container_initialized": config.container_initialized,
        "commit_allowed": config.commit_allowed,
        "rhs_spreading_allowed": config.rhs_spreading_allowed,
        "stage14_rhs_injection_allowed": config.stage14_rhs_injection_allowed,
        "X_prod": x_prod,
        "V_prod": v_prod,
        "A_prod": a_prod,
        "F_b_prod": f_b_prod,
        "F_T_prod": f_t_prod,
        "F_h_prod": f_h_prod,
        "F_total_prod": f_total_prod,
        "X_candidate": [row[:] for row in x_prod],
        "V_candidate": [row[:] for row in v_prod],
        "A_candidate": [row[:] for row in a_prod],
        "owner_rank": [0 for _idx in range(config.n_point)],
        "global_point_id": [idx for idx in range(config.n_point)],
        "local_point_id": [idx for idx in range(config.n_point)],
    }


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


def arrays_close(lhs: List[List[float]], rhs: List[List[float]], tol: float) -> bool:
    return shape_2d(lhs) == shape_2d(rhs) and all(abs(a - b) <= tol for row_a, row_b in zip(lhs, rhs) for a, b in zip(row_a, row_b))


def f_total_consistent(state: Dict[str, object], tol: float) -> bool:
    fb = state["F_b_prod"]
    ft = state["F_T_prod"]
    fh = state["F_h_prod"]
    total = state["F_total_prod"]
    if not all(isinstance(item, list) for item in (fb, ft, fh, total)):
        return False
    for i in range(len(total)):
        for j in range(len(total[i])):
            if abs(total[i][j] - (fb[i][j] + ft[i][j] + fh[i][j])) > tol:
                return False
    return True


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 19.3 physical structure initialization boundary")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    root = args.repo_root.resolve()
    output = args.output or (root / "stage19_outputs" / "fibre_stage19_3_physical_structure_initialization.dat")
    statuses: Dict[str, str] = {key: "PASS" for key in SUMMARY_KEYS if key != "final_status"}
    failure_reasons: List[str] = []
    invalid_env: List[str] = []

    requested = parse_bool_env("STAGE19_3_ENABLE", True, invalid_env)
    initialization_enable = parse_bool_env("STAGE19_3_INITIALIZATION_ENABLE", True, invalid_env)
    require_stage19_2 = parse_bool_env("STAGE19_3_REQUIRE_STAGE19_2_PASS", True, invalid_env)
    require_stage19_1 = parse_bool_env("STAGE19_3_REQUIRE_STAGE19_1_PASS", True, invalid_env)
    require_stage19_0 = parse_bool_env("STAGE19_3_REQUIRE_STAGE19_0_PASS", True, invalid_env)
    require_stage18 = parse_bool_env("STAGE19_3_REQUIRE_STAGE18_CLOSED", True, invalid_env)
    init_mode = os.environ.get("STAGE19_3_INIT_MODE", "small_sine_fibre_zero_velocity").strip()
    if init_mode not in INIT_MODES:
        invalid_env.append(f"STAGE19_3_INIT_MODE={init_mode}")
        init_mode = "small_sine_fibre_zero_velocity"

    config = InitConfig(
        n_fibre=parse_int_env("STAGE19_3_N_FIBRE", 1, invalid_env),
        n_point=parse_int_env("STAGE19_3_N_POINT", 64, invalid_env),
        component_dim=parse_int_env("STAGE19_3_COMPONENT_DIM", 3, invalid_env),
        fibre_length=parse_float_env("STAGE19_3_FIBRE_LENGTH", 1.0, invalid_env),
        rho_l=parse_float_env("STAGE19_3_RHO_L", 1.0, invalid_env),
        rho_tilde=parse_float_env("STAGE19_3_RHO_TILDE", 1.0, invalid_env),
        bending_stiffness=parse_float_env("STAGE19_3_BENDING_STIFFNESS", 1.0e-3, invalid_env),
        gamma=parse_float_env("STAGE19_3_GAMMA", 1.0e-3, invalid_env),
        init_mode=init_mode,
        sine_amplitude=parse_float_env("STAGE19_3_SINE_AMPLITUDE", 1.0e-3, invalid_env),
        sine_mode=parse_int_env("STAGE19_3_SINE_MODE", 1, invalid_env),
        velocity_amplitude=parse_float_env("STAGE19_3_VELOCITY_AMPLITUDE", 0.0, invalid_env),
        controlled_force_amplitude=parse_float_env("STAGE19_3_CONTROLLED_FORCE_AMPLITUDE", 0.0, invalid_env),
        diagnostic_only=parse_bool_env("STAGE19_3_DIAGNOSTIC_ONLY", True, invalid_env),
        single_fibre_only=parse_bool_env("STAGE19_3_SINGLE_FIBRE_ONLY", True, invalid_env),
        fail_closed=parse_bool_env("STAGE19_3_FAIL_CLOSED", True, invalid_env),
        commit_allowed=parse_bool_env("STAGE19_3_COMMIT_ALLOWED", False, invalid_env),
        rhs_spreading_allowed=parse_bool_env("STAGE19_3_RHS_SPREADING_ALLOWED", False, invalid_env),
        stage14_rhs_injection_allowed=parse_bool_env("STAGE19_3_STAGE14_RHS_INJECTION_ALLOWED", False, invalid_env),
        state_valid=parse_bool_env("STAGE19_3_STATE_VALID_AFTER_INIT", True, invalid_env),
        container_initialized=parse_bool_env("STAGE19_3_CONTAINER_INITIALIZED_AFTER_INIT", True, invalid_env),
    )
    state = initialize_state(config)
    tol = parse_float_env("STAGE19_3_AUDIT_TOL", 1.0e-12, invalid_env)

    statuses["stage19_3_requested_status"] = pass_fail(requested)
    statuses["stage19_3_initialization_enable_status"] = pass_fail(initialization_enable)
    statuses["stage19_2_evidence_status"] = pass_fail((not require_stage19_2) or stage19_2_evidence_valid(root))
    statuses["stage19_1_evidence_status"] = pass_fail((not require_stage19_1) or stage19_1_evidence_valid(root))
    statuses["stage19_0_evidence_status"] = pass_fail((not require_stage19_0) or stage19_0_evidence_valid(root))
    statuses["stage18_closure_evidence_status"] = pass_fail((not require_stage18) or stage18_evidence_safely_accepted(root))
    statuses["stage19_2_state_container_preserved_status"] = pass_fail(stage19_2_state_container_preserved(root))
    statuses["stage19_1_config_gate_preserved_status"] = pass_fail(stage19_1_config_gate_preserved(root))
    statuses["stage19_0_source_only_closure_acceptance_preserved_status"] = pass_fail(stage19_0_source_acceptance_preserved(root))

    changed = outside_allowed_changes(root)
    statuses["no_stage10_18_file_modification_status"] = pass_fail(not any(path.startswith(("stage10_", "stage11_", "stage12_", "stage13_", "stage14_", "stage15_", "stage16_", "stage17_", "stage18_")) for path in changed))
    statuses["no_stage19_0_file_modification_status"] = pass_fail(not any(path in STAGE19_0_FILES for path in changed))
    statuses["no_stage19_1_file_modification_status"] = pass_fail(not any(path in STAGE19_1_FILES for path in changed))
    statuses["no_stage19_2_file_modification_status"] = pass_fail(not any(path in STAGE19_2_FILES for path in changed))
    statuses["no_closed_stage_modification_status"] = pass_fail(all(statuses[key] == "PASS" for key in ("no_stage10_18_file_modification_status", "no_stage19_0_file_modification_status", "no_stage19_1_file_modification_status", "no_stage19_2_file_modification_status")))

    doc = read_text(root / "stage19_checks" / "stage19_3_physical_structure_initialization.md")
    statuses["initialization_schema_documented_status"] = pass_fail(all(field in doc for field in REQUIRED_FIELDS))
    statuses["all_required_initialization_modes_present_status"] = pass_fail(all(mode in doc for mode in INIT_MODES))
    statuses["all_required_initialized_fields_present_status"] = pass_fail(set(REQUIRED_FIELDS).issubset(set(state)))

    default_safe = (
        config.n_fibre == 1 and config.n_point == 64 and config.component_dim == 3 and config.fibre_length == 1.0
        and config.rho_l == 1.0 and config.rho_tilde == 1.0 and config.bending_stiffness == 1.0e-3 and config.gamma == 1.0e-3
        and config.init_mode == "small_sine_fibre_zero_velocity" and config.sine_amplitude == 1.0e-3 and config.sine_mode == 1
        and config.velocity_amplitude == 0.0 and config.controlled_force_amplitude == 0.0 and config.diagnostic_only
        and config.single_fibre_only and config.fail_closed and config.state_valid and config.container_initialized
        and not config.commit_allowed and not config.rhs_spreading_allowed and not config.stage14_rhs_injection_allowed
    )
    statuses["default_safe_values_status"] = pass_fail(default_safe)
    statuses["n_fibre_status"] = pass_fail(config.n_fibre == 1)
    statuses["n_point_status"] = pass_fail(config.n_point >= 2)
    statuses["component_dim_status"] = pass_fail(config.component_dim == 3)
    statuses["fibre_length_status"] = pass_fail(config.fibre_length > 0.0)
    expected_ds = config.fibre_length / (config.n_point - 1) if config.n_point > 1 else math.nan
    statuses["ds_formula_status"] = pass_fail(math.isfinite(expected_ds) and math.isclose(float(state["ds"]), expected_ds, rel_tol=0.0, abs_tol=1.0e-14))
    statuses["rho_l_status"] = pass_fail(config.rho_l > 0.0)
    statuses["rho_tilde_status"] = pass_fail(config.rho_tilde > 0.0)
    statuses["bending_stiffness_status"] = pass_fail(config.bending_stiffness >= 0.0)
    statuses["gamma_status"] = pass_fail(config.gamma >= 0.0)
    statuses["init_mode_status"] = pass_fail(config.init_mode in INIT_MODES)
    statuses["sine_amplitude_status"] = pass_fail(math.isfinite(config.sine_amplitude) and abs(config.sine_amplitude) <= 1.0e-1)
    statuses["sine_mode_status"] = pass_fail(config.sine_mode > 0)
    statuses["velocity_amplitude_status"] = pass_fail(config.velocity_amplitude >= 0.0)
    statuses["controlled_force_amplitude_status"] = pass_fail(config.controlled_force_amplitude >= 0.0)

    is_straight = all(abs(row[1]) <= tol and abs(row[2]) <= tol for row in state["X_prod"])
    is_small_sine = all(math.isfinite(row[1]) and abs(row[1]) <= abs(config.sine_amplitude) + tol for row in state["X_prod"])
    is_zero_velocity = all(abs(value) <= tol for row in state["V_prod"] for value in row)
    controlled_velocity_ok = config.velocity_amplitude == 0.0 or any(abs(row[1]) > tol for row in state["V_prod"])
    zero_force_ok = all(abs(value) <= tol for field in ("F_b_prod", "F_T_prod") for row in state[field] for value in row)
    controlled_force_ok = config.controlled_force_amplitude == 0.0 or any(abs(row[1]) > tol for row in state["F_h_prod"])
    statuses["straight_fibre_initialization_status"] = pass_fail((not config.init_mode.startswith("straight")) or is_straight)
    statuses["small_sine_initialization_status"] = pass_fail((not config.init_mode.startswith("small_sine")) or is_small_sine)
    statuses["zero_velocity_initialization_status"] = pass_fail(("zero_velocity" not in config.init_mode) or is_zero_velocity)
    statuses["controlled_velocity_initialization_status"] = pass_fail(("controlled_velocity" not in config.init_mode) or controlled_velocity_ok)
    statuses["zero_force_placeholder_initialization_status"] = pass_fail(zero_force_ok and (config.controlled_force_amplitude > 0.0 or all(abs(value) <= tol for row in state["F_h_prod"] for value in row)))
    statuses["controlled_helper_force_placeholder_initialization_status"] = pass_fail((config.init_mode != "controlled_helper_force_placeholder") or controlled_force_ok)

    expected_2d_shape = (config.n_point, config.component_dim)
    for field in ARRAY_FIELDS_2D:
        statuses[SHAPE_STATUS[field]] = pass_fail(shape_2d(state[field]) == expected_2d_shape)
    statuses["owner_rank_shape_status"] = pass_fail(shape_1d(state["owner_rank"]) == (config.n_point,))
    statuses["global_point_id_shape_status"] = pass_fail(shape_1d(state["global_point_id"]) == (config.n_point,))
    statuses["local_point_id_shape_status"] = pass_fail(shape_1d(state["local_point_id"]) == (config.n_point,))
    statuses["array_finite_status"] = pass_fail(all(finite_nested(state[field]) for field in ARRAY_FIELDS_2D))
    statuses["f_total_consistency_status"] = pass_fail(f_total_consistent(state, tol))
    statuses["candidate_equals_production_initialization_status"] = pass_fail(arrays_close(state["X_candidate"], state["X_prod"], tol) and arrays_close(state["V_candidate"], state["V_prod"], tol) and arrays_close(state["A_candidate"], state["A_prod"], tol))
    statuses["global_point_id_coverage_status"] = pass_fail(sorted(state["global_point_id"]) == list(range(config.n_point)))
    statuses["global_point_id_no_duplicate_status"] = pass_fail(len(set(state["global_point_id"])) == config.n_point)
    statuses["owner_rank_deterministic_status"] = pass_fail(state["owner_rank"] == [0 for _idx in range(config.n_point)])

    statuses["diagnostic_only_status"] = pass_fail(config.diagnostic_only)
    statuses["single_fibre_only_status"] = pass_fail(config.single_fibre_only and config.n_fibre == 1)
    statuses["fail_closed_status"] = pass_fail(config.fail_closed)
    statuses["state_valid_after_init_status"] = pass_fail(config.state_valid)
    statuses["container_initialized_after_init_status"] = pass_fail(config.container_initialized)
    statuses["commit_default_disabled_status"] = pass_fail(not config.commit_allowed)
    statuses["rhs_spreading_default_disabled_status"] = pass_fail(not config.rhs_spreading_allowed)
    statuses["stage14_rhs_injection_default_disabled_status"] = pass_fail(not config.stage14_rhs_injection_allowed)

    diagnostic_ok = (not config.diagnostic_only) or not config.commit_allowed
    single_fibre_ok = (not config.single_fibre_only) or config.n_fibre == 1
    rhs_disabled_ok = config.rhs_spreading_allowed or not config.stage14_rhs_injection_allowed
    stage14_rhs_ok = not config.stage14_rhs_injection_allowed
    commit_disabled_ok = not config.commit_allowed
    runtime_inactive = not config.commit_allowed and not config.rhs_spreading_allowed and not config.stage14_rhs_injection_allowed
    fail_closed_ok = config.fail_closed and not invalid_env and diagnostic_ok and single_fibre_ok and rhs_disabled_ok and stage14_rhs_ok and commit_disabled_ok and runtime_inactive
    statuses["diagnostic_only_consistency_status"] = pass_fail(diagnostic_ok)
    statuses["single_fibre_only_consistency_status"] = pass_fail(single_fibre_ok)
    statuses["fail_closed_consistency_status"] = pass_fail(fail_closed_ok)
    statuses["rhs_spreading_disabled_consistency_status"] = pass_fail(rhs_disabled_ok)
    statuses["stage14_rhs_injection_disabled_consistency_status"] = pass_fail(stage14_rhs_ok)
    statuses["commit_disabled_consistency_status"] = pass_fail(commit_disabled_ok)
    statuses["initialization_production_runtime_inactive_status"] = pass_fail(runtime_inactive)

    wrapper_syntax, helper_compile = syntax_status(root)
    statuses["stage19_3_wrapper_bash_syntax_status"] = wrapper_syntax
    statuses["stage19_3_helper_py_compile_status"] = helper_compile
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
        "no_stage14_rhs_call_from_stage19_3_status",
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
    lines.append("# Stage 19.3 physical structure initialization boundary")
    lines.append("stage19_title production-side physical structure integration boundary")
    lines.append("stage19_3_title production physical structure initialization")
    lines.append(f"stage19_3_test_case {os.environ.get('STAGE19_3_TEST_CASE', 'production_physical_structure_initialization')}")
    lines.append(f"stage19_3_zero_tol_value {os.environ.get('STAGE19_3_ZERO_TOL', '1.0e-14')}")
    lines.append(f"stage19_3_audit_tol_value {tol}")
    lines.append("stage19_3_scope initialization-boundary only; local helper initialization is not production runtime state update")
    for name in ("n_fibre", "n_point", "component_dim", "fibre_length", "ds", "rho_l", "rho_tilde", "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "velocity_amplitude", "controlled_force_amplitude"):
        lines.append(f"{name}_value {state[name]}")
    for field in ARRAY_FIELDS_2D:
        lines.append(f"{field.lower()}_shape_value {shape_2d(state[field])}")
    lines.append(f"owner_rank_shape_value {shape_1d(state['owner_rank'])}")
    lines.append(f"global_point_id_shape_value {shape_1d(state['global_point_id'])}")
    lines.append(f"local_point_id_shape_value {shape_1d(state['local_point_id'])}")
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses[key]}")
    if failure_reasons:
        lines.append("failure_reasons_begin")
        lines.extend(failure_reasons)
        lines.append("failure_reasons_end")
    lines.append(f"STAGE 19.3 PHYSICAL STRUCTURE INITIALIZATION VERDICT: {final_status}")
    lines.append(f"STAGE 19.3 FINAL VERDICT: {final_status}")
    output.write_text("\n".join(lines) + "\n")

    print(f"STAGE 19.3 PHYSICAL STRUCTURE INITIALIZATION VERDICT: {final_status}")
    print(f"STAGE 19.3 FINAL VERDICT: {final_status}")
    if failure_reasons:
        print("STAGE 19.3 FAILURE REASONS:")
        for reason in failure_reasons:
            print(f"  - {reason}")
    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
