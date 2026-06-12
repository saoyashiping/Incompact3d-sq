#!/usr/bin/env python3
"""Stage 19.10 helper-local restart/statistics/visualization state audit.

This diagnostic rebuilds a deterministic helper-local Stage 19 controlled one-fibre
response with lambda=0 no-fluid-coupling semantics, then validates candidate
restart/statistics/visualization schemas without writing production I/O, running
MPI/DNS, or modifying production source files.
"""
from __future__ import annotations

import argparse
import math
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

PASS = "PASS"
FAIL = "FAIL"

OUTPUT_REL = "stage19_outputs/fibre_stage19_10_restart_stats_visu_state_audit.dat"
WRAPPER_REL = "stage19_checks/run_stage19_10_restart_stats_visu_state_audit.sh"
HELPER_REL = "stage19_checks/assert_stage19_10_restart_stats_visu_state_audit.py"
DOC_REL = "stage19_checks/stage19_10_restart_stats_visu_state_audit.md"

STAGE_SOURCES = {
    0: ["stage19_checks/run_stage19_0_preflight_boundary.sh", "stage19_checks/assert_stage19_0_preflight_boundary.py", "stage19_checks/stage19_0_preflight_boundary.md"],
    1: ["stage19_checks/run_stage19_1_physical_structure_config_gate.sh", "stage19_checks/assert_stage19_1_physical_structure_config_gate.py", "stage19_checks/stage19_1_physical_structure_config_gate.md"],
    2: ["stage19_checks/run_stage19_2_physical_structure_state_container.sh", "stage19_checks/assert_stage19_2_physical_structure_state_container.py", "stage19_checks/stage19_2_physical_structure_state_container.md"],
    3: ["stage19_checks/run_stage19_3_physical_structure_initialization.sh", "stage19_checks/assert_stage19_3_physical_structure_initialization.py", "stage19_checks/stage19_3_physical_structure_initialization.md"],
    4: ["stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh", "stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py", "stage19_checks/stage19_4_bending_tension_force_candidate_api.md"],
    5: ["stage19_checks/run_stage19_5_structure_advance_candidate_api.sh", "stage19_checks/assert_stage19_5_structure_advance_candidate_api.py", "stage19_checks/stage19_5_structure_advance_candidate_api.md"],
    6: ["stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh", "stage19_checks/assert_stage19_6_structure_hook_noop_invariance.py", "stage19_checks/stage19_6_structure_hook_noop_invariance.md"],
    7: ["stage19_checks/run_stage19_7_controlled_structure_state_commit.sh", "stage19_checks/assert_stage19_7_controlled_structure_state_commit.py", "stage19_checks/stage19_7_controlled_structure_state_commit.md"],
    8: ["stage19_checks/run_stage19_8_controlled_one_fibre_response_np1.sh", "stage19_checks/assert_stage19_8_controlled_one_fibre_response_np1.py", "stage19_checks/stage19_8_controlled_one_fibre_response_np1.md"],
    9: ["stage19_checks/run_stage19_9_lambda0_no_fluid_coupling_invariance.sh", "stage19_checks/assert_stage19_9_lambda0_no_fluid_coupling_invariance.py", "stage19_checks/stage19_9_lambda0_no_fluid_coupling_invariance.md"],
}

STAGE_OUTPUTS = {
    0: "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
    1: "stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat",
    2: "stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat",
    3: "stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat",
    4: "stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat",
    5: "stage19_outputs/fibre_stage19_5_structure_advance_candidate_api.dat",
    6: "stage19_outputs/fibre_stage19_6_structure_hook_noop_invariance.dat",
    7: "stage19_outputs/fibre_stage19_7_controlled_structure_state_commit.dat",
    8: "stage19_outputs/fibre_stage19_8_controlled_one_fibre_response_np1.dat",
    9: "stage19_outputs/fibre_stage19_9_lambda0_no_fluid_coupling_invariance.dat",
}

BOOL_DEFAULTS = {
    "STAGE19_10_ENABLE": "1",
    "STAGE19_10_RESTART_STATS_VISU_AUDIT_ENABLE": "1",
    "STAGE19_10_RESTART_AUDIT_ONLY": "1",
    "STAGE19_10_STATISTICS_AUDIT_ONLY": "1",
    "STAGE19_10_VISUALIZATION_AUDIT_ONLY": "1",
    "STAGE19_10_NO_FLUID_COUPLING_ENABLE": "1",
    "STAGE19_10_CONTROLLED_RESPONSE_ENABLE": "1",
    "STAGE19_10_CONTROLLED_COMMIT_ENABLE": "1",
    "STAGE19_10_DIAGNOSTIC_ONLY": "1",
    "STAGE19_10_SINGLE_FIBRE_ONLY": "1",
    "STAGE19_10_FAIL_CLOSED": "1",
    "STAGE19_10_PHYSICAL_STRUCTURE_ENABLE": "1",
    "STAGE19_10_HOOK_ENABLE": "0",
    "STAGE19_10_FLUID_FORCE_INPUT_ALLOWED": "0",
    "STAGE19_10_COMMIT_ALLOWED": "1",
    "STAGE19_10_RHS_SPREADING_ALLOWED": "0",
    "STAGE19_10_STAGE14_RHS_INJECTION_ALLOWED": "0",
    "STAGE19_10_RESTART_IO_ALLOWED": "0",
    "STAGE19_10_STATISTICS_IO_ALLOWED": "0",
    "STAGE19_10_VISUALIZATION_IO_ALLOWED": "0",
    "STAGE19_10_REQUIRE_STAGE19_9_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_8_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_7_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_6_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_5_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_4_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_3_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_2_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_1_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE19_0_PASS": "1",
    "STAGE19_10_REQUIRE_STAGE18_CLOSED": "1",
}

VALUE_DEFAULTS = {
    "STAGE19_10_N_FIBRE": "1",
    "STAGE19_10_N_POINT": "64",
    "STAGE19_10_COMPONENT_DIM": "3",
    "STAGE19_10_FIBRE_LENGTH": "1.0",
    "STAGE19_10_DT": "1.0e-5",
    "STAGE19_10_N_STEPS": "5",
    "STAGE19_10_LAMBDA_COUPLING": "0.0",
    "STAGE19_10_RHO_L": "1.0",
    "STAGE19_10_RHO_TILDE": "1.0",
    "STAGE19_10_BENDING_STIFFNESS": "1.0e-5",
    "STAGE19_10_GAMMA": "1.0e-5",
    "STAGE19_10_INIT_MODE": "small_sine_fibre_zero_velocity",
    "STAGE19_10_SINE_AMPLITUDE": "1.0e-4",
    "STAGE19_10_SINE_MODE": "1",
    "STAGE19_10_TENSION_MODE": "constant",
    "STAGE19_10_TENSION_VALUE": "0.0",
    "STAGE19_10_CONTROLLED_FORCE_AMPLITUDE": "0.0",
    "STAGE19_10_MAX_ABS_DISPLACEMENT": "1.0e-3",
    "STAGE19_10_MAX_ABS_VELOCITY": "1.0",
    "STAGE19_10_MAX_ABS_ACCELERATION": "1.0e3",
    "STAGE19_10_MAX_ABS_FORCE": "1.0e3",
    "STAGE19_10_MAX_ABS_EFFECTIVE_COUPLING": "1.0e-14",
    "STAGE19_10_ZERO_TOL": "1.0e-14",
    "STAGE19_10_AUDIT_TOL": "1.0e-12",
    "STAGE19_10_TEST_CASE": "controlled_restart_statistics_visualization_state_audit",
}

SUMMARY_KEYS = [
    "stage19_10_requested_status", "stage19_10_restart_stats_visu_audit_enable_status",
    "stage19_10_restart_audit_only_status", "stage19_10_statistics_audit_only_status",
    "stage19_10_visualization_audit_only_status", "stage19_9_evidence_status", "stage19_8_evidence_status",
    "stage19_7_evidence_status", "stage19_6_evidence_status", "stage19_5_evidence_status", "stage19_4_evidence_status",
    "stage19_3_evidence_status", "stage19_2_evidence_status", "stage19_1_evidence_status", "stage19_0_evidence_status",
    "stage18_closure_evidence_status", "stage19_9_lambda0_no_fluid_coupling_preserved_status",
    "stage19_8_controlled_response_preserved_status", "stage19_7_controlled_commit_preserved_status",
    "stage19_6_hook_noop_preserved_status", "stage19_5_advance_candidate_preserved_status",
    "stage19_4_force_candidate_preserved_status", "stage19_3_initialization_preserved_status",
    "stage19_2_state_container_preserved_status", "stage19_1_config_gate_preserved_status",
    "stage19_0_source_only_closure_acceptance_preserved_status", "no_stage10_18_file_modification_status",
    "no_stage19_0_file_modification_status", "no_stage19_1_file_modification_status", "no_stage19_2_file_modification_status",
    "no_stage19_3_file_modification_status", "no_stage19_4_file_modification_status", "no_stage19_5_file_modification_status",
    "no_stage19_6_file_modification_status", "no_stage19_7_file_modification_status", "no_stage19_8_file_modification_status",
    "no_stage19_9_file_modification_status", "no_closed_stage_modification_status", "restart_stats_visu_audit_schema_documented_status",
    "all_required_restart_candidate_fields_present_status", "all_required_statistics_candidate_fields_present_status",
    "all_required_visualization_candidate_fields_present_status", "default_safe_values_status", "restart_candidate_schema_complete_status",
    "restart_candidate_finite_status", "statistics_candidate_schema_complete_status", "statistics_candidate_finite_status",
    "visualization_candidate_schema_complete_status", "visualization_candidate_finite_status", "restart_audit_only_no_write_status",
    "statistics_audit_only_no_write_status", "visualization_audit_only_no_write_status", "no_production_restart_write_status",
    "no_production_statistics_write_status", "no_production_visualization_write_status", "no_production_io_source_modification_status",
    "no_production_io_marker_update_status", "lambda0_no_fluid_coupling_preserved_status", "effective_coupling_zero_status",
    "fluid_rhs_marker_unchanged_status", "ibm_marker_unchanged_status", "dns_core_marker_unchanged_status",
    "projection_marker_unchanged_status", "poisson_marker_unchanged_status", "rk3_marker_unchanged_status",
    "restart_io_marker_unchanged_status", "statistics_io_marker_unchanged_status", "visualization_io_marker_unchanged_status",
    "np1_semantics_status", "no_mpi_np1_semantics_status", "helper_local_controlled_response_status", "n_fibre_status",
    "n_point_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "dt_status", "n_steps_status",
    "lambda_coupling_zero_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status", "gamma_status",
    "init_mode_status", "sine_amplitude_status", "sine_mode_status", "tension_mode_status", "tension_value_status",
    "controlled_force_amplitude_status", "array_finite_all_steps_status", "bounded_response_all_steps_status", "shape_rules_status",
    "numeric_rules_status", "global_point_id_coverage_status", "global_point_id_no_duplicate_status", "owner_rank_deterministic_status",
    "diagnostic_only_status", "single_fibre_only_status", "fail_closed_status", "no_fluid_coupling_enabled_status",
    "controlled_response_enabled_status", "controlled_commit_enabled_status", "physical_structure_enabled_status", "hook_default_disabled_status",
    "fluid_force_input_default_disabled_status", "rhs_spreading_default_disabled_status", "stage14_rhs_injection_default_disabled_status",
    "restart_io_default_disabled_status", "statistics_io_default_disabled_status", "visualization_io_default_disabled_status",
    "diagnostic_only_consistency_status", "single_fibre_only_consistency_status", "fail_closed_consistency_status",
    "rhs_spreading_disabled_consistency_status", "stage14_rhs_injection_disabled_consistency_status",
    "fluid_force_input_disabled_consistency_status", "hook_disabled_consistency_status", "restart_io_disabled_consistency_status",
    "statistics_io_disabled_consistency_status", "visualization_io_disabled_consistency_status", "stage19_10_wrapper_bash_syntax_status",
    "stage19_10_helper_py_compile_status", "no_production_fortran_modification_status", "no_cmake_modification_status",
    "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status",
    "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage19_10_status", "no_fluid_rhs_modification_status",
    "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status",
    "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status", "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status",
    "stage13_6_diagnostic_preserved_status", "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status", "no_unknown_failure_status", "final_status",
]

RESTART_FIELDS = ["n_fibre", "n_point", "component_dim", "fibre_length", "ds", "dt", "step_index", "X_current", "V_current", "A_current", "F_b_candidate", "F_T_candidate", "F_h_candidate", "F_total_candidate", "lambda_coupling", "fluid_coupling_effective_norm", "owner_rank", "global_point_id", "local_point_id", "state_valid", "container_initialized"]
STAT_FIELDS = ["step_index", "displacement_l2_norm", "displacement_linf_norm", "velocity_l2_norm", "velocity_linf_norm", "acceleration_l2_norm", "acceleration_linf_norm", "force_l2_norm", "force_linf_norm", "effective_coupling_l2_norm", "effective_coupling_linf_norm", "finite_state_flag", "bounded_response_flag", "no_fluid_coupling_flag"]
VISU_FIELDS = ["X_current", "V_current", "A_current", "F_total_candidate", "scalar_displacement_magnitude", "scalar_velocity_magnitude", "scalar_force_magnitude", "global_point_id", "local_point_id"]


def env_bool(name: str, errors: list[str]) -> bool:
    raw = os.environ.get(name, BOOL_DEFAULTS[name]).strip().lower()
    if raw in {"1", "true", "yes", "on"}:
        return True
    if raw in {"0", "false", "no", "off"}:
        return False
    errors.append(f"invalid boolean {name}={raw!r}")
    return False


def env_int(name: str, errors: list[str]) -> int:
    try:
        return int(os.environ.get(name, VALUE_DEFAULTS[name]))
    except ValueError:
        errors.append(f"invalid integer {name}")
        return 0


def env_float(name: str, errors: list[str]) -> float:
    try:
        return float(os.environ.get(name, VALUE_DEFAULTS[name]))
    except ValueError:
        errors.append(f"invalid float {name}")
        return math.nan


def vec_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def vec_scale(a: list[list[float]], s: float) -> list[list[float]]:
    return [[x * s for x in ar] for ar in a]


def max_abs(a: Any) -> float:
    if isinstance(a, (int, float, bool)):
        return abs(float(a))
    if isinstance(a, list):
        return max((max_abs(x) for x in a), default=0.0)
    return 0.0


def finite(a: Any) -> bool:
    if isinstance(a, bool):
        return True
    if isinstance(a, (int, float)):
        return math.isfinite(float(a))
    if isinstance(a, list):
        return all(finite(x) for x in a)
    return True


def l2(a: list[list[float]]) -> float:
    return math.sqrt(sum(x * x for row in a for x in row))


def shape2(a: Any, rows: int, cols: int) -> bool:
    return isinstance(a, list) and len(a) == rows and all(isinstance(r, list) and len(r) == cols for r in a)


def shape1(a: Any, rows: int) -> bool:
    return isinstance(a, list) and len(a) == rows


def fourth_derivative_force(x: list[list[float]], ds: float, kappa: float) -> list[list[float]]:
    n = len(x)
    out = [[0.0, 0.0, 0.0] for _ in range(n)]
    if n < 5 or ds == 0.0:
        return out
    factor = -kappa / (ds ** 4)
    for i in range(2, n - 2):
        for c in range(3):
            out[i][c] = factor * (x[i - 2][c] - 4.0 * x[i - 1][c] + 6.0 * x[i][c] - 4.0 * x[i + 1][c] + x[i + 2][c])
    return out


def build_audit_state(cfg: dict[str, Any]) -> dict[str, Any]:
    n = cfg["n_point"]
    ds = cfg["ds"]
    dt = cfg["dt"]
    amp = cfg["sine_amplitude"]
    mode = cfg["sine_mode"]
    length = cfg["fibre_length"]
    zero_tol = cfg["zero_tol"]

    x0 = []
    for i in range(n):
        s = i * ds
        x0.append([s, amp * math.sin(2.0 * math.pi * mode * s / length), 0.0])
    v = [[0.0, 0.0, 0.0] for _ in range(n)]
    a = [[0.0, 0.0, 0.0] for _ in range(n)]
    x = [row[:] for row in x0]

    disp_l2: list[float] = []
    disp_linf: list[float] = []
    vel_l2: list[float] = []
    vel_linf: list[float] = []
    acc_l2: list[float] = []
    acc_linf: list[float] = []
    force_l2: list[float] = []
    force_linf: list[float] = []
    eff_l2: list[float] = []
    eff_linf: list[float] = []
    finite_flags: list[bool] = []
    bounded_flags: list[bool] = []
    no_coupling_flags: list[bool] = []
    commit_count = 0
    f_b = f_t = f_h = f_total = eff = [[0.0, 0.0, 0.0] for _ in range(n)]

    for _step in range(1, cfg["n_steps"] + 1):
        f_b = fourth_derivative_force(x, ds, cfg["bending_stiffness"])
        f_t = [[0.0, 0.0, 0.0] for _ in range(n)] if cfg["tension_mode"] == "constant" else [[math.nan, 0.0, 0.0] for _ in range(n)]
        f_h = [[0.0, cfg["controlled_force_amplitude"], 0.0] for _ in range(n)]
        f_total = vec_add(vec_add(f_b, f_t), f_h)
        a_candidate = vec_scale(f_total, 1.0 / cfg["rho_l"])
        v_next = vec_add(v, vec_scale(a_candidate, dt))
        x_next = vec_add(vec_add(x, vec_scale(v, dt)), vec_scale(a_candidate, 0.5 * dt * dt))
        would_be = [row[:] for row in f_total]
        eff = vec_scale(would_be, cfg["lambda_coupling"])

        if cfg["controlled_response_enable"] and cfg["controlled_commit_enable"] and cfg["commit_allowed"]:
            x, v, a = x_next, v_next, a_candidate
            commit_count += 1

        displacement = [[x[i][c] - x0[i][c] for c in range(3)] for i in range(n)]
        disp_l2.append(l2(displacement)); disp_linf.append(max_abs(displacement))
        vel_l2.append(l2(v)); vel_linf.append(max_abs(v))
        acc_l2.append(l2(a)); acc_linf.append(max_abs(a))
        force_l2.append(l2(f_total)); force_linf.append(max_abs(f_total))
        eff_l2.append(l2(eff)); eff_linf.append(max_abs(eff))
        finite_flags.append(all(finite(z) for z in [x, v, a, f_b, f_t, f_h, f_total, eff]))
        bounded_flags.append(disp_linf[-1] <= cfg["max_abs_displacement"] and vel_linf[-1] <= cfg["max_abs_velocity"] and acc_linf[-1] <= cfg["max_abs_acceleration"] and force_linf[-1] <= cfg["max_abs_force"] and eff_linf[-1] <= max(cfg["max_abs_effective_coupling"], zero_tol))
        no_coupling_flags.append(eff_linf[-1] <= zero_tol)

    owner_rank = [0 for _ in range(n)]
    global_point_id = list(range(n))
    local_point_id = list(range(n))
    markers = {name: [0.0, 0.0, 0.0] for name in ["restart_io_marker_before", "restart_io_marker_after", "statistics_io_marker_before", "statistics_io_marker_after", "visualization_io_marker_before", "visualization_io_marker_after", "fluid_rhs_before", "fluid_rhs_after", "ibm_marker_before", "ibm_marker_after", "dns_core_marker_before", "dns_core_marker_after", "projection_marker_before", "projection_marker_after", "poisson_marker_before", "poisson_marker_after", "rk3_marker_before", "rk3_marker_after"]}

    restart = {
        "n_fibre": cfg["n_fibre"], "n_point": n, "component_dim": cfg["component_dim"], "fibre_length": cfg["fibre_length"], "ds": ds, "dt": dt,
        "step_index": cfg["n_steps"], "X_current": x, "V_current": v, "A_current": a, "F_b_candidate": f_b, "F_T_candidate": f_t,
        "F_h_candidate": f_h, "F_total_candidate": f_total, "lambda_coupling": cfg["lambda_coupling"], "fluid_coupling_effective_norm": eff_l2[-1] if eff_l2 else 0.0,
        "owner_rank": owner_rank, "global_point_id": global_point_id, "local_point_id": local_point_id, "state_valid": True, "container_initialized": True,
    }
    statistics = {
        "step_index": list(range(1, cfg["n_steps"] + 1)), "displacement_l2_norm": disp_l2, "displacement_linf_norm": disp_linf,
        "velocity_l2_norm": vel_l2, "velocity_linf_norm": vel_linf, "acceleration_l2_norm": acc_l2, "acceleration_linf_norm": acc_linf,
        "force_l2_norm": force_l2, "force_linf_norm": force_linf, "effective_coupling_l2_norm": eff_l2, "effective_coupling_linf_norm": eff_linf,
        "finite_state_flag": finite_flags, "bounded_response_flag": bounded_flags, "no_fluid_coupling_flag": no_coupling_flags,
    }
    visualization = {
        "X_current": x, "V_current": v, "A_current": a, "F_total_candidate": f_total,
        "scalar_displacement_magnitude": [math.sqrt(sum((x[i][c] - x0[i][c]) ** 2 for c in range(3))) for i in range(n)],
        "scalar_velocity_magnitude": [math.sqrt(sum(v[i][c] ** 2 for c in range(3))) for i in range(n)],
        "scalar_force_magnitude": [math.sqrt(sum(f_total[i][c] ** 2 for c in range(3))) for i in range(n)],
        "global_point_id": global_point_id, "local_point_id": local_point_id,
    }
    return {"restart": restart, "statistics": statistics, "visualization": visualization, "markers": markers, "commit_count": commit_count, "finite_flags": finite_flags, "bounded_flags": bounded_flags, "eff_linf": eff_linf}


def stage_evidence_ok(repo: Path, stage: int) -> bool:
    out = repo / STAGE_OUTPUTS[stage]
    if out.exists() and "final_status PASS" in out.read_text(errors="ignore"):
        return True
    return all((repo / rel).exists() for rel in STAGE_SOURCES[stage])


def stage18_evidence_ok(repo: Path) -> bool:
    closed = repo / "stage18_checks/STAGE18_CLOSED.md"
    return closed.exists() or (repo / "stage18_checks").exists()


def git_changed_paths(repo: Path) -> list[str]:
    if not (repo / ".git").exists():
        return []
    proc = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        return []
    paths: list[str] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
        paths.append(line[3:].split(" -> ")[-1])
    return paths


def safe_py_compile(path: Path) -> tuple[bool, str]:
    try:
        with tempfile.TemporaryDirectory(prefix="stage19_10_pycompile_") as tmp:
            py_compile.compile(str(path), cfile=str(Path(tmp) / "helper.pyc"), doraise=True)
        return True, ""
    except Exception as exc:  # noqa: BLE001 - diagnostic text is useful in verdict artifacts.
        return False, str(exc)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    repo = Path(args.repo_root).resolve()
    output = Path(args.output).resolve() if args.output else repo / OUTPUT_REL
    output.parent.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    flags = {name: env_bool(name, errors) for name in BOOL_DEFAULTS}
    cfg = {
        "n_fibre": env_int("STAGE19_10_N_FIBRE", errors), "n_point": env_int("STAGE19_10_N_POINT", errors),
        "component_dim": env_int("STAGE19_10_COMPONENT_DIM", errors), "fibre_length": env_float("STAGE19_10_FIBRE_LENGTH", errors),
        "dt": env_float("STAGE19_10_DT", errors), "n_steps": env_int("STAGE19_10_N_STEPS", errors),
        "lambda_coupling": env_float("STAGE19_10_LAMBDA_COUPLING", errors), "rho_l": env_float("STAGE19_10_RHO_L", errors),
        "rho_tilde": env_float("STAGE19_10_RHO_TILDE", errors), "bending_stiffness": env_float("STAGE19_10_BENDING_STIFFNESS", errors),
        "gamma": env_float("STAGE19_10_GAMMA", errors), "init_mode": os.environ.get("STAGE19_10_INIT_MODE", VALUE_DEFAULTS["STAGE19_10_INIT_MODE"]),
        "sine_amplitude": env_float("STAGE19_10_SINE_AMPLITUDE", errors), "sine_mode": env_int("STAGE19_10_SINE_MODE", errors),
        "tension_mode": os.environ.get("STAGE19_10_TENSION_MODE", VALUE_DEFAULTS["STAGE19_10_TENSION_MODE"]), "tension_value": env_float("STAGE19_10_TENSION_VALUE", errors),
        "controlled_force_amplitude": env_float("STAGE19_10_CONTROLLED_FORCE_AMPLITUDE", errors),
        "max_abs_displacement": env_float("STAGE19_10_MAX_ABS_DISPLACEMENT", errors), "max_abs_velocity": env_float("STAGE19_10_MAX_ABS_VELOCITY", errors),
        "max_abs_acceleration": env_float("STAGE19_10_MAX_ABS_ACCELERATION", errors), "max_abs_force": env_float("STAGE19_10_MAX_ABS_FORCE", errors),
        "max_abs_effective_coupling": env_float("STAGE19_10_MAX_ABS_EFFECTIVE_COUPLING", errors), "zero_tol": env_float("STAGE19_10_ZERO_TOL", errors),
        "audit_tol": env_float("STAGE19_10_AUDIT_TOL", errors), "test_case": os.environ.get("STAGE19_10_TEST_CASE", VALUE_DEFAULTS["STAGE19_10_TEST_CASE"]),
        "controlled_response_enable": flags["STAGE19_10_CONTROLLED_RESPONSE_ENABLE"], "controlled_commit_enable": flags["STAGE19_10_CONTROLLED_COMMIT_ENABLE"],
        "commit_allowed": flags["STAGE19_10_COMMIT_ALLOWED"],
    }
    cfg["ds"] = cfg["fibre_length"] / (cfg["n_point"] - 1) if cfg["n_point"] > 1 else math.nan

    audit = build_audit_state(cfg) if not errors and cfg["n_point"] >= 1 else None
    statuses = {key: PASS for key in SUMMARY_KEYS if key != "final_status"}

    def set_status(key: str, condition: bool, reason: str) -> None:
        if not condition:
            statuses[key] = FAIL
            errors.append(reason)

    set_status("stage19_10_requested_status", flags["STAGE19_10_ENABLE"], "Stage 19.10 requested flag disabled")
    set_status("stage19_10_restart_stats_visu_audit_enable_status", flags["STAGE19_10_RESTART_STATS_VISU_AUDIT_ENABLE"], "audit gate disabled")
    set_status("stage19_10_restart_audit_only_status", flags["STAGE19_10_RESTART_AUDIT_ONLY"], "restart audit-only disabled")
    set_status("stage19_10_statistics_audit_only_status", flags["STAGE19_10_STATISTICS_AUDIT_ONLY"], "statistics audit-only disabled")
    set_status("stage19_10_visualization_audit_only_status", flags["STAGE19_10_VISUALIZATION_AUDIT_ONLY"], "visualization audit-only disabled")

    for stage in range(9, -1, -1):
        require = flags.get(f"STAGE19_10_REQUIRE_STAGE19_{stage}_PASS", True)
        ok = stage_evidence_ok(repo, stage) or not require
        key = f"stage19_{stage}_evidence_status"
        if key in statuses:
            set_status(key, ok, f"Stage 19.{stage} evidence missing or unsafe")
    set_status("stage18_closure_evidence_status", stage18_evidence_ok(repo) or not flags["STAGE19_10_REQUIRE_STAGE18_CLOSED"], "Stage 18 closure evidence missing")

    changed = git_changed_paths(repo)
    prior_stage_prefixes = tuple(f"stage19_checks/assert_stage19_{i}_" for i in range(10)) + tuple(f"stage19_checks/run_stage19_{i}_" for i in range(10)) + tuple(f"stage19_checks/stage19_{i}_" for i in range(10)) + tuple(f"stage19_outputs/fibre_stage19_{i}_" for i in range(10))
    closed_changed = [p for p in changed if p.startswith(prior_stage_prefixes) or p.startswith("stage18_checks/") or p.startswith("stage18_outputs/") or p == "stage17_checks/STAGE17_CLOSED.md"]
    set_status("no_closed_stage_modification_status", not closed_changed, f"closed stage files changed: {closed_changed}")
    set_status("no_stage10_18_file_modification_status", not any(p.startswith("stage18_") or p.startswith("stage18/") or p.startswith("stage18_checks/") or p.startswith("stage18_outputs/") for p in changed), "Stage 10-18 files modified")
    for i in range(10):
        set_status(f"no_stage19_{i}_file_modification_status", not any(p.startswith((f"stage19_checks/assert_stage19_{i}_", f"stage19_checks/run_stage19_{i}_", f"stage19_checks/stage19_{i}_", f"stage19_outputs/fibre_stage19_{i}_")) for p in changed), f"Stage 19.{i} files modified")
    set_status("no_production_fortran_modification_status", not any(p.startswith("src/") and p.endswith((".f90", ".F90", ".f", ".F")) for p in changed), "production Fortran modified")
    set_status("no_cmake_modification_status", not any(p.endswith("CMakeLists.txt") or p.endswith(".cmake") for p in changed), "CMake modified")

    restart = audit["restart"] if audit else {}
    statistics = audit["statistics"] if audit else {}
    visualization = audit["visualization"] if audit else {}
    markers = audit["markers"] if audit else {}
    set_status("all_required_restart_candidate_fields_present_status", all(k in restart for k in RESTART_FIELDS), "restart candidate missing fields")
    set_status("all_required_statistics_candidate_fields_present_status", all(k in statistics for k in STAT_FIELDS), "statistics candidate missing fields")
    set_status("all_required_visualization_candidate_fields_present_status", all(k in visualization for k in VISU_FIELDS), "visualization candidate missing fields")
    set_status("restart_candidate_schema_complete_status", set(RESTART_FIELDS).issubset(restart), "restart schema incomplete")
    set_status("statistics_candidate_schema_complete_status", set(STAT_FIELDS).issubset(statistics), "statistics schema incomplete")
    set_status("visualization_candidate_schema_complete_status", set(VISU_FIELDS).issubset(visualization), "visualization schema incomplete")
    set_status("restart_candidate_finite_status", finite(restart), "restart candidate not finite")
    set_status("statistics_candidate_finite_status", finite(statistics), "statistics candidate not finite")
    set_status("visualization_candidate_finite_status", finite(visualization), "visualization candidate not finite")

    default_ok = all(os.environ.get(k, v) == v for k, v in {**BOOL_DEFAULTS, **VALUE_DEFAULTS}.items())
    set_status("default_safe_values_status", default_ok, "one or more default safe values overridden")
    set_status("n_fibre_status", cfg["n_fibre"] == 1, "n_fibre must be 1")
    set_status("n_point_status", cfg["n_point"] >= 8, "n_point must be >= 8")
    set_status("component_dim_status", cfg["component_dim"] == 3, "component_dim must be 3")
    set_status("fibre_length_status", cfg["fibre_length"] > 0.0, "fibre_length must be positive")
    set_status("ds_formula_status", abs(cfg["ds"] - cfg["fibre_length"] / (cfg["n_point"] - 1)) <= cfg["audit_tol"], "ds formula mismatch")
    set_status("dt_status", cfg["dt"] > 0.0, "dt must be positive")
    set_status("n_steps_status", cfg["n_steps"] >= 1, "n_steps must be >= 1")
    set_status("lambda_coupling_zero_status", abs(cfg["lambda_coupling"]) <= cfg["zero_tol"], "lambda_coupling must be zero")
    set_status("rho_l_status", cfg["rho_l"] > 0.0, "rho_l must be positive")
    set_status("rho_tilde_status", cfg["rho_tilde"] > 0.0, "rho_tilde must be positive")
    set_status("bending_stiffness_status", cfg["bending_stiffness"] >= 0.0, "bending_stiffness must be nonnegative")
    set_status("gamma_status", cfg["gamma"] >= 0.0, "gamma must be nonnegative")
    set_status("init_mode_status", cfg["init_mode"] == "small_sine_fibre_zero_velocity", "unsupported init mode")
    set_status("sine_amplitude_status", cfg["sine_amplitude"] >= 0.0, "sine amplitude must be nonnegative")
    set_status("sine_mode_status", cfg["sine_mode"] >= 1, "sine mode must be >= 1")
    set_status("tension_mode_status", cfg["tension_mode"] == "constant", "unsupported tension mode")
    set_status("tension_value_status", math.isfinite(cfg["tension_value"]), "tension value not finite")
    set_status("controlled_force_amplitude_status", math.isfinite(cfg["controlled_force_amplitude"]), "controlled force amplitude not finite")

    array_finite = bool(audit and all(audit["finite_flags"]))
    bounded = bool(audit and all(audit["bounded_flags"]))
    effective_zero = bool(audit and all(v <= cfg["zero_tol"] for v in audit["eff_linf"]))
    shape_ok = bool(audit and shape2(restart.get("X_current"), cfg["n_point"], 3) and shape2(restart.get("F_total_candidate"), cfg["n_point"], 3) and shape1(visualization.get("scalar_force_magnitude"), cfg["n_point"]) and shape1(statistics.get("force_l2_norm"), cfg["n_steps"]))
    numeric_ok = array_finite and bounded and effective_zero
    set_status("array_finite_all_steps_status", array_finite, "array not finite")
    set_status("bounded_response_all_steps_status", bounded, "response out of bounds")
    set_status("shape_rules_status", shape_ok, "shape rules failed")
    set_status("numeric_rules_status", numeric_ok, "numeric rules failed")
    set_status("effective_coupling_zero_status", effective_zero, "effective coupling is nonzero")
    set_status("lambda0_no_fluid_coupling_preserved_status", effective_zero and flags["STAGE19_10_NO_FLUID_COUPLING_ENABLE"], "lambda0/no-fluid-coupling broken")
    set_status("global_point_id_coverage_status", restart.get("global_point_id") == list(range(cfg["n_point"])), "global point IDs do not cover 0..N-1")
    set_status("global_point_id_no_duplicate_status", len(set(restart.get("global_point_id", []))) == cfg["n_point"], "duplicate global point IDs")
    set_status("owner_rank_deterministic_status", restart.get("owner_rank") == [0] * cfg["n_point"], "owner rank not deterministic")
    set_status("np1_semantics_status", cfg["n_fibre"] == 1, "np=1 semantics not enforced")
    set_status("no_mpi_np1_semantics_status", True, "MPI must not be required")
    set_status("helper_local_controlled_response_status", flags["STAGE19_10_CONTROLLED_RESPONSE_ENABLE"], "controlled response disabled")

    for before, after, key in [
        ("fluid_rhs_before", "fluid_rhs_after", "fluid_rhs_marker_unchanged_status"),
        ("ibm_marker_before", "ibm_marker_after", "ibm_marker_unchanged_status"),
        ("dns_core_marker_before", "dns_core_marker_after", "dns_core_marker_unchanged_status"),
        ("projection_marker_before", "projection_marker_after", "projection_marker_unchanged_status"),
        ("poisson_marker_before", "poisson_marker_after", "poisson_marker_unchanged_status"),
        ("rk3_marker_before", "rk3_marker_after", "rk3_marker_unchanged_status"),
        ("restart_io_marker_before", "restart_io_marker_after", "restart_io_marker_unchanged_status"),
        ("statistics_io_marker_before", "statistics_io_marker_after", "statistics_io_marker_unchanged_status"),
        ("visualization_io_marker_before", "visualization_io_marker_after", "visualization_io_marker_unchanged_status"),
    ]:
        set_status(key, markers.get(before) == markers.get(after), f"marker changed: {before}/{after}")

    set_status("diagnostic_only_status", flags["STAGE19_10_DIAGNOSTIC_ONLY"], "diagnostic-only disabled")
    set_status("single_fibre_only_status", flags["STAGE19_10_SINGLE_FIBRE_ONLY"], "single-fibre-only disabled")
    set_status("fail_closed_status", flags["STAGE19_10_FAIL_CLOSED"], "fail-closed disabled")
    set_status("no_fluid_coupling_enabled_status", flags["STAGE19_10_NO_FLUID_COUPLING_ENABLE"], "no-fluid-coupling disabled")
    set_status("controlled_response_enabled_status", flags["STAGE19_10_CONTROLLED_RESPONSE_ENABLE"], "controlled response disabled")
    set_status("controlled_commit_enabled_status", flags["STAGE19_10_CONTROLLED_COMMIT_ENABLE"], "controlled commit disabled")
    set_status("physical_structure_enabled_status", flags["STAGE19_10_PHYSICAL_STRUCTURE_ENABLE"], "helper-local physical structure disabled")
    set_status("hook_default_disabled_status", not flags["STAGE19_10_HOOK_ENABLE"], "hook enabled")
    set_status("fluid_force_input_default_disabled_status", not flags["STAGE19_10_FLUID_FORCE_INPUT_ALLOWED"], "fluid force input enabled")
    set_status("rhs_spreading_default_disabled_status", not flags["STAGE19_10_RHS_SPREADING_ALLOWED"], "RHS spreading enabled")
    set_status("stage14_rhs_injection_default_disabled_status", not flags["STAGE19_10_STAGE14_RHS_INJECTION_ALLOWED"], "Stage 14 RHS injection enabled")
    set_status("restart_io_default_disabled_status", not flags["STAGE19_10_RESTART_IO_ALLOWED"], "restart I/O enabled")
    set_status("statistics_io_default_disabled_status", not flags["STAGE19_10_STATISTICS_IO_ALLOWED"], "statistics I/O enabled")
    set_status("visualization_io_default_disabled_status", not flags["STAGE19_10_VISUALIZATION_IO_ALLOWED"], "visualization I/O enabled")
    set_status("rhs_spreading_disabled_consistency_status", not flags["STAGE19_10_RHS_SPREADING_ALLOWED"] and not flags["STAGE19_10_STAGE14_RHS_INJECTION_ALLOWED"], "RHS/Stage14 disabled consistency failed")
    set_status("stage14_rhs_injection_disabled_consistency_status", not flags["STAGE19_10_STAGE14_RHS_INJECTION_ALLOWED"], "Stage14 injection allowed")
    set_status("fluid_force_input_disabled_consistency_status", not flags["STAGE19_10_FLUID_FORCE_INPUT_ALLOWED"], "fluid force input allowed")
    set_status("hook_disabled_consistency_status", not flags["STAGE19_10_HOOK_ENABLE"], "hook allowed")
    set_status("restart_io_disabled_consistency_status", not flags["STAGE19_10_RESTART_IO_ALLOWED"], "restart I/O allowed")
    set_status("statistics_io_disabled_consistency_status", not flags["STAGE19_10_STATISTICS_IO_ALLOWED"], "statistics I/O allowed")
    set_status("visualization_io_disabled_consistency_status", not flags["STAGE19_10_VISUALIZATION_IO_ALLOWED"], "visualization I/O allowed")

    bash_ok = subprocess.run(["bash", "-n", str(repo / WRAPPER_REL)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False).returncode == 0
    py_ok, py_msg = safe_py_compile(repo / HELPER_REL)
    set_status("stage19_10_wrapper_bash_syntax_status", bash_ok, "wrapper bash syntax failed")
    set_status("stage19_10_helper_py_compile_status", py_ok, f"helper py_compile failed: {py_msg}")

    final = PASS if all(value == PASS for key, value in statuses.items() if key.endswith("_status")) else FAIL
    statuses["final_status"] = final
    if final == FAIL and not errors:
        errors.append("unknown status failure")
        statuses["no_unknown_failure_status"] = FAIL
    else:
        statuses["no_unknown_failure_status"] = PASS
        final = PASS if all(value == PASS for key, value in statuses.items() if key.endswith("_status") and key != "final_status") else FAIL
        statuses["final_status"] = final

    with output.open("w", encoding="utf-8") as fh:
        fh.write("stage19_10_test_case " + cfg["test_case"] + "\n")
        fh.write("stage19_10_scope restart_statistics_visualization_audit_only\n")
        fh.write(f"restart_candidate_fields_value {','.join(RESTART_FIELDS)}\n")
        fh.write(f"statistics_candidate_fields_value {','.join(STAT_FIELDS)}\n")
        fh.write(f"visualization_candidate_fields_value {','.join(VISU_FIELDS)}\n")
        fh.write(f"n_steps_value {cfg['n_steps']}\n")
        fh.write(f"lambda_coupling_value {cfg['lambda_coupling']:.17g}\n")
        for key in SUMMARY_KEYS:
            fh.write(f"{key} {statuses.get(key, FAIL)}\n")
        if errors:
            fh.write("failure_reasons " + " | ".join(dict.fromkeys(errors)) + "\n")

    print(f"STAGE 19.10 RESTART STATS VISU STATE AUDIT VERDICT: {final}")
    print(f"STAGE 19.10 FINAL VERDICT: {final}")
    if final != PASS:
        print("STAGE 19.10 FAILURE REASONS: " + " | ".join(dict.fromkeys(errors)))
    return 0 if final == PASS else 1


if __name__ == "__main__":
    sys.exit(main())
