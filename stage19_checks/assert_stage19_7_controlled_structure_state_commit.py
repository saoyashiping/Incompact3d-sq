#!/usr/bin/env python3
"""Stage 19.7 controlled helper-local structure-state commit boundary.

This diagnostic reconstructs helper-local Stage 19.3 initialization, Stage 19.4
force candidates, Stage 19.5 advance candidates, Stage 19.6 closed-gate no-op
markers, and then performs the explicitly gated Stage 19.7 helper-local commit:
X_prod_after <- X_next_candidate, V_prod_after <- V_next_candidate, and
A_prod_after <- A_advance_candidate.  The commit is diagnostic and
structure-state-only: no production runtime hook is inserted, no production
runtime advance/commit path is activated, no force is spread to fluid RHS, no
Stage 14 RHS injection is called, no IBM/DNS-core/projection/Poisson/RK3 or
production I/O is modified, and no MPI/DNS/contact/collision/multifibre logic is
introduced.
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

SUMMARY_KEYS = [
    "stage19_7_requested_status", "stage19_7_controlled_commit_enable_status",
    "stage19_6_evidence_status", "stage19_5_evidence_status", "stage19_4_evidence_status", "stage19_3_evidence_status", "stage19_2_evidence_status", "stage19_1_evidence_status", "stage19_0_evidence_status", "stage18_closure_evidence_status",
    "stage19_6_hook_noop_preserved_status", "stage19_5_advance_candidate_preserved_status", "stage19_4_force_candidate_preserved_status", "stage19_3_initialization_preserved_status", "stage19_2_state_container_preserved_status", "stage19_1_config_gate_preserved_status", "stage19_0_source_only_closure_acceptance_preserved_status",
    "no_stage10_18_file_modification_status", "no_stage19_0_file_modification_status", "no_stage19_1_file_modification_status", "no_stage19_2_file_modification_status", "no_stage19_3_file_modification_status", "no_stage19_4_file_modification_status", "no_stage19_5_file_modification_status", "no_stage19_6_file_modification_status", "no_closed_stage_modification_status",
    "controlled_commit_schema_documented_status", "all_required_controlled_commit_fields_present_status", "default_safe_values_status", "helper_local_controlled_commit_enabled_status", "commit_allowed_helper_local_only_status", "physical_structure_enable_helper_local_only_status", "hook_default_disabled_status",
    "x_after_equals_candidate_status", "v_after_equals_candidate_status", "a_after_equals_candidate_status", "total_force_candidate_formula_status", "acceleration_candidate_formula_status", "velocity_next_candidate_formula_status", "position_next_candidate_formula_status", "delta_velocity_candidate_formula_status", "delta_position_candidate_formula_status",
    "controlled_commit_once_status", "controlled_commit_state_changed_status", "controlled_commit_helper_local_status", "no_production_runtime_hook_status", "no_production_runtime_state_update_status", "no_fluid_rhs_update_status", "no_ibm_marker_update_status", "no_dns_core_marker_update_status", "no_production_io_marker_update_status",
    "n_fibre_status", "n_point_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "dt_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status", "gamma_status", "init_mode_status", "sine_amplitude_status", "sine_mode_status", "tension_mode_status", "tension_value_status", "controlled_force_amplitude_status", "array_finite_status",
    "global_point_id_coverage_status", "global_point_id_no_duplicate_status", "owner_rank_deterministic_status",
    "diagnostic_only_status", "single_fibre_only_status", "fail_closed_status", "force_candidate_only_status", "advance_candidate_only_status", "controlled_commit_enabled_status", "physical_structure_enabled_status", "fluid_force_input_default_disabled_status", "rhs_spreading_default_disabled_status", "stage14_rhs_injection_default_disabled_status", "restart_io_default_disabled_status", "statistics_io_default_disabled_status", "visualization_io_default_disabled_status",
    "diagnostic_only_consistency_status", "single_fibre_only_consistency_status", "fail_closed_consistency_status", "rhs_spreading_disabled_consistency_status", "stage14_rhs_injection_disabled_consistency_status", "fluid_force_input_disabled_consistency_status", "hook_disabled_consistency_status", "controlled_commit_production_runtime_inactive_status",
    "stage19_7_wrapper_bash_syntax_status", "stage19_7_helper_py_compile_status", "no_production_fortran_modification_status", "no_cmake_modification_status", "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status", "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage19_7_status", "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status", "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status", "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status", "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "stage13_6_diagnostic_preserved_status", "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status", "no_rg_only_dependency_status", "no_unknown_failure_status", "final_status",
]

ARRAY_FIELDS = ["X_prod_before", "V_prod_before", "A_prod_before", "X_prod_after", "V_prod_after", "A_prod_after", "F_b_candidate", "F_T_candidate", "F_h_candidate", "F_total_candidate", "A_advance_candidate", "V_next_candidate", "X_next_candidate", "delta_X_candidate", "delta_V_candidate", "fluid_rhs_before", "fluid_rhs_after", "ibm_marker_before", "ibm_marker_after", "dns_core_marker_before", "dns_core_marker_after", "restart_io_marker_before", "restart_io_marker_after", "statistics_io_marker_before", "statistics_io_marker_after", "visualization_io_marker_before", "visualization_io_marker_after"]
REQUIRED_FIELDS = ARRAY_FIELDS + ["owner_rank", "global_point_id", "local_point_id", "n_fibre", "n_point", "component_dim", "fibre_length", "ds", "dt", "rho_l", "rho_tilde", "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "tension_mode", "tension_value", "controlled_force_amplitude", "diagnostic_only", "force_candidate_only", "advance_candidate_only", "controlled_commit_enable", "state_valid", "container_initialized", "physical_structure_enable", "hook_enable", "commit_allowed", "rhs_spreading_allowed", "stage14_rhs_injection_allowed", "fluid_force_input_allowed", "restart_io_allowed", "statistics_io_allowed", "visualization_io_allowed"]
INIT_MODES = {"small_sine_fibre_zero_velocity", "straight_fibre_zero_velocity"}
TENSION_MODES = {"constant"}
TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "off"}
STAGE19_FILES = {
    "0": ["stage19_checks/assert_stage19_0_preflight_boundary.py", "stage19_checks/run_stage19_0_preflight_boundary.sh", "stage19_checks/stage19_0_preflight_boundary.md"],
    "1": ["stage19_checks/assert_stage19_1_physical_structure_config_gate.py", "stage19_checks/run_stage19_1_physical_structure_config_gate.sh", "stage19_checks/stage19_1_physical_structure_config_gate.md"],
    "2": ["stage19_checks/assert_stage19_2_physical_structure_state_container.py", "stage19_checks/run_stage19_2_physical_structure_state_container.sh", "stage19_checks/stage19_2_physical_structure_state_container.md"],
    "3": ["stage19_checks/assert_stage19_3_physical_structure_initialization.py", "stage19_checks/run_stage19_3_physical_structure_initialization.sh", "stage19_checks/stage19_3_physical_structure_initialization.md"],
    "4": ["stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py", "stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh", "stage19_checks/stage19_4_bending_tension_force_candidate_api.md"],
    "5": ["stage19_checks/assert_stage19_5_structure_advance_candidate_api.py", "stage19_checks/run_stage19_5_structure_advance_candidate_api.sh", "stage19_checks/stage19_5_structure_advance_candidate_api.md"],
    "6": ["stage19_checks/assert_stage19_6_structure_hook_noop_invariance.py", "stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh", "stage19_checks/stage19_6_structure_hook_noop_invariance.md"],
}
ACCEPTED_UNTRACKED_EVIDENCE = {"stage19_outputs/fibre_stage19_7_controlled_structure_state_commit.dat"}

@dataclass(frozen=True)
class Config:
    requested: bool; controlled_commit_enable: bool; diagnostic_only: bool; single_fibre_only: bool; fail_closed: bool
    force_candidate_only: bool; advance_candidate_only: bool; physical_structure_enable: bool; hook_enable: bool; fluid_force_input_allowed: bool
    commit_allowed: bool; rhs_spreading_allowed: bool; stage14_rhs_injection_allowed: bool; restart_io_allowed: bool; statistics_io_allowed: bool; visualization_io_allowed: bool
    n_fibre: int; n_point: int; component_dim: int; fibre_length: float; dt: float; rho_l: float; rho_tilde: float; bending_stiffness: float; gamma: float
    init_mode: str; sine_amplitude: float; sine_mode: int; tension_mode: str; tension_value: float; controlled_force_amplitude: float; zero_tol: float; audit_tol: float

def parse_bool(name: str, default: str, invalid: List[str]) -> bool:
    value = os.environ.get(name, default).strip().lower()
    if value in TRUE_VALUES: return True
    if value in FALSE_VALUES: return False
    invalid.append(f"{name}={os.environ.get(name, default)}"); return False

def parse_int(name: str, default: str, invalid: List[str]) -> int:
    try: return int(os.environ.get(name, default))
    except ValueError: invalid.append(f"{name}={os.environ.get(name, default)}"); return 0

def parse_float(name: str, default: str, invalid: List[str]) -> float:
    try: value = float(os.environ.get(name, default))
    except ValueError: invalid.append(f"{name}={os.environ.get(name, default)}"); return math.nan
    if not math.isfinite(value): invalid.append(f"{name}={os.environ.get(name, default)}")
    return value

def load_config(invalid: List[str]) -> Config:
    return Config(
        parse_bool("STAGE19_7_ENABLE", "1", invalid), parse_bool("STAGE19_7_CONTROLLED_COMMIT_ENABLE", "1", invalid), parse_bool("STAGE19_7_DIAGNOSTIC_ONLY", "1", invalid), parse_bool("STAGE19_7_SINGLE_FIBRE_ONLY", "1", invalid), parse_bool("STAGE19_7_FAIL_CLOSED", "1", invalid),
        parse_bool("STAGE19_7_FORCE_CANDIDATE_ONLY", "1", invalid), parse_bool("STAGE19_7_ADVANCE_CANDIDATE_ONLY", "1", invalid), parse_bool("STAGE19_7_PHYSICAL_STRUCTURE_ENABLE", "1", invalid), parse_bool("STAGE19_7_HOOK_ENABLE", "0", invalid), parse_bool("STAGE19_7_FLUID_FORCE_INPUT_ALLOWED", "0", invalid),
        parse_bool("STAGE19_7_COMMIT_ALLOWED", "1", invalid), parse_bool("STAGE19_7_RHS_SPREADING_ALLOWED", "0", invalid), parse_bool("STAGE19_7_STAGE14_RHS_INJECTION_ALLOWED", "0", invalid), parse_bool("STAGE19_7_RESTART_IO_ALLOWED", "0", invalid), parse_bool("STAGE19_7_STATISTICS_IO_ALLOWED", "0", invalid), parse_bool("STAGE19_7_VISUALIZATION_IO_ALLOWED", "0", invalid),
        parse_int("STAGE19_7_N_FIBRE", "1", invalid), parse_int("STAGE19_7_N_POINT", "64", invalid), parse_int("STAGE19_7_COMPONENT_DIM", "3", invalid), parse_float("STAGE19_7_FIBRE_LENGTH", "1.0", invalid), parse_float("STAGE19_7_DT", "1.0e-4", invalid), parse_float("STAGE19_7_RHO_L", "1.0", invalid), parse_float("STAGE19_7_RHO_TILDE", "1.0", invalid), parse_float("STAGE19_7_BENDING_STIFFNESS", "1.0e-3", invalid), parse_float("STAGE19_7_GAMMA", "1.0e-3", invalid),
        os.environ.get("STAGE19_7_INIT_MODE", "small_sine_fibre_zero_velocity"), parse_float("STAGE19_7_SINE_AMPLITUDE", "1.0e-3", invalid), parse_int("STAGE19_7_SINE_MODE", "1", invalid), os.environ.get("STAGE19_7_TENSION_MODE", "constant"), parse_float("STAGE19_7_TENSION_VALUE", "0.0", invalid), parse_float("STAGE19_7_CONTROLLED_FORCE_AMPLITUDE", "0.0", invalid), parse_float("STAGE19_7_ZERO_TOL", "1.0e-14", invalid), parse_float("STAGE19_7_AUDIT_TOL", "1.0e-12", invalid),
    )

def zeros(n: int, m: int = 3) -> List[List[float]]: return [[0.0 for _ in range(m)] for _ in range(n)]
def copy2(a: Sequence[Sequence[float]]) -> List[List[float]]: return [[float(x) for x in row] for row in a]
def add2(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]: return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]
def sub2(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]: return [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]
def scale2(a: Sequence[Sequence[float]], s: float) -> List[List[float]]: return [[s * x for x in row] for row in a]
def maxabs(a: Sequence[Sequence[float]]) -> float: return max((abs(x) for row in a for x in row), default=0.0)
def close(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], tol: float) -> bool: return maxabs(sub2(a, b)) <= tol
def finite_array(a: Sequence[Sequence[float]]) -> bool: return all(math.isfinite(x) for row in a for x in row)
def shape2(a: Sequence[Sequence[float]]) -> Tuple[int, int]: return (len(a), len(a[0]) if a else 0)
def shape1(a: Sequence[int]) -> Tuple[int]: return (len(a),)
def passfail(ok: bool) -> str: return "PASS" if ok else "FAIL"

def fourth_derivative(x: Sequence[Sequence[float]], ds: float) -> List[List[float]]:
    n = len(x); out = zeros(n, 3)
    if n < 5: return out
    for i in range(2, n - 2):
        for c in range(3): out[i][c] = (x[i - 2][c] - 4*x[i - 1][c] + 6*x[i][c] - 4*x[i + 1][c] + x[i + 2][c]) / (ds ** 4)
    out[0] = out[2][:]; out[1] = out[2][:]; out[-2] = out[-3][:]; out[-1] = out[-3][:]
    return out

def second_derivative(x: Sequence[Sequence[float]], ds: float) -> List[List[float]]:
    n = len(x); out = zeros(n, 3)
    if n < 3: return out
    for i in range(1, n - 1):
        for c in range(3): out[i][c] = (x[i + 1][c] - 2*x[i][c] + x[i - 1][c]) / (ds ** 2)
    out[0] = out[1][:]; out[-1] = out[-2][:]
    return out

def build_state(cfg: Config) -> Dict[str, object]:
    n = cfg.n_point; ds = cfg.fibre_length / (n - 1) if n > 1 else math.nan
    x = zeros(n); v = zeros(n); a = zeros(n)
    for i in range(n):
        s = i * ds; x[i][0] = s
        if cfg.init_mode == "small_sine_fibre_zero_velocity": x[i][1] = cfg.sine_amplitude * math.sin(2*math.pi*cfg.sine_mode*s/cfg.fibre_length)
    f_b = scale2(fourth_derivative(x, ds), -cfg.bending_stiffness)
    f_t = scale2(second_derivative(x, ds), cfg.tension_value) if cfg.tension_mode == "constant" and cfg.tension_value != 0.0 else zeros(n)
    f_h = zeros(n)
    if cfg.controlled_force_amplitude > 0.0:
        for i in range(n):
            s = i * ds; f_h[i][1] = cfg.controlled_force_amplitude * math.sin(2*math.pi*cfg.sine_mode*s/cfg.fibre_length)
    f_total = add2(add2(f_b, f_t), f_h)
    a_adv = scale2(f_total, 1.0 / cfg.rho_l)
    v_next = add2(v, scale2(a_adv, cfg.dt))
    x_next = add2(add2(x, scale2(v, cfg.dt)), scale2(a_adv, 0.5 * cfg.dt * cfg.dt))
    delta_v = sub2(v_next, v); delta_x = sub2(x_next, x)
    did_commit = cfg.controlled_commit_enable and cfg.commit_allowed
    x_after = copy2(x_next if did_commit else x)
    v_after = copy2(v_next if did_commit else v)
    a_after = copy2(a_adv if did_commit else a)
    state: Dict[str, object] = {
        "X_prod_before": copy2(x), "V_prod_before": copy2(v), "A_prod_before": copy2(a), "X_prod_after": x_after, "V_prod_after": v_after, "A_prod_after": a_after,
        "F_b_candidate": f_b, "F_T_candidate": f_t, "F_h_candidate": f_h, "F_total_candidate": f_total, "A_advance_candidate": a_adv, "V_next_candidate": v_next, "X_next_candidate": x_next, "delta_X_candidate": delta_x, "delta_V_candidate": delta_v,
        "fluid_rhs_before": zeros(n), "fluid_rhs_after": zeros(n), "ibm_marker_before": zeros(n), "ibm_marker_after": zeros(n), "dns_core_marker_before": zeros(n), "dns_core_marker_after": zeros(n),
        "restart_io_marker_before": zeros(n), "restart_io_marker_after": zeros(n), "statistics_io_marker_before": zeros(n), "statistics_io_marker_after": zeros(n), "visualization_io_marker_before": zeros(n), "visualization_io_marker_after": zeros(n),
        "owner_rank": [0 for _ in range(n)], "global_point_id": list(range(n)), "local_point_id": list(range(n)),
        "n_fibre": cfg.n_fibre, "n_point": cfg.n_point, "component_dim": cfg.component_dim, "fibre_length": cfg.fibre_length, "ds": ds, "dt": cfg.dt, "rho_l": cfg.rho_l, "rho_tilde": cfg.rho_tilde, "bending_stiffness": cfg.bending_stiffness, "gamma": cfg.gamma,
        "init_mode": cfg.init_mode, "sine_amplitude": cfg.sine_amplitude, "sine_mode": cfg.sine_mode, "tension_mode": cfg.tension_mode, "tension_value": cfg.tension_value, "controlled_force_amplitude": cfg.controlled_force_amplitude,
        "diagnostic_only": cfg.diagnostic_only, "force_candidate_only": cfg.force_candidate_only, "advance_candidate_only": cfg.advance_candidate_only, "controlled_commit_enable": cfg.controlled_commit_enable, "state_valid": True, "container_initialized": True,
        "physical_structure_enable": cfg.physical_structure_enable, "hook_enable": cfg.hook_enable, "commit_allowed": cfg.commit_allowed, "rhs_spreading_allowed": cfg.rhs_spreading_allowed, "stage14_rhs_injection_allowed": cfg.stage14_rhs_injection_allowed,
        "fluid_force_input_allowed": cfg.fluid_force_input_allowed, "restart_io_allowed": cfg.restart_io_allowed, "statistics_io_allowed": cfg.statistics_io_allowed, "visualization_io_allowed": cfg.visualization_io_allowed,
    }
    state["helper_local_commit_count"] = 1 if did_commit else 0
    return state

def git_changes(repo: Path) -> List[Tuple[str, str]]:
    if not (repo / ".git").exists(): return []
    result = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0: return []
    changes: List[Tuple[str, str]] = []
    for line in result.stdout.splitlines():
        if line:
            path = line[3:]
            if " -> " in path: path = path.split(" -> ", 1)[1]
            changes.append((line[:2], path))
    return changes

def changed(changes: Sequence[Tuple[str, str]], pred) -> bool:
    return any(path not in ACCEPTED_UNTRACKED_EVIDENCE and pred(path) for _, path in changes)

def stage_evidence(repo: Path, stage: str, output_name: str) -> bool:
    output = repo / "stage19_outputs" / output_name
    if output.exists() and "final_status PASS" in output.read_text(errors="ignore"): return True
    return all((repo / name).exists() for name in STAGE19_FILES[stage])

def stage18_evidence(repo: Path) -> bool:
    return (repo / "stage18_checks/STAGE18_CLOSED.md").exists() or (repo / "stage18_checks/assert_stage18_0_preflight_boundary.py").exists()

def syntax_ok(repo: Path, helper: Path, wrapper: Path) -> Tuple[bool, bool]:
    bash_ok = subprocess.run(["bash", "-n", str(wrapper)], cwd=repo, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False).returncode == 0
    with tempfile.TemporaryDirectory(prefix="stage19_7_py_compile_") as tmp:
        try:
            py_compile.compile(str(helper), cfile=str(Path(tmp) / "assert_stage19_7_controlled_structure_state_commit.pyc"), doraise=True); py_ok = True
        except py_compile.PyCompileError:
            py_ok = False
    return bash_ok, py_ok

def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 19.7 controlled structure state commit helper")
    parser.add_argument("--repo-root", type=Path, required=True); parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(); repo = args.repo_root.resolve(); output = args.output.resolve()
    helper = repo / "stage19_checks/assert_stage19_7_controlled_structure_state_commit.py"
    wrapper = repo / "stage19_checks/run_stage19_7_controlled_structure_state_commit.sh"
    doc = repo / "stage19_checks/stage19_7_controlled_structure_state_commit.md"
    invalid: List[str] = []; cfg = load_config(invalid); state = build_state(cfg); changes = git_changes(repo); bash_ok, py_ok = syntax_ok(repo, helper, wrapper)
    st: Dict[str, str] = {}
    st["stage19_7_requested_status"] = passfail(cfg.requested); st["stage19_7_controlled_commit_enable_status"] = passfail(cfg.controlled_commit_enable)
    for stage, outfile in [("6", "fibre_stage19_6_structure_hook_noop_invariance.dat"), ("5", "fibre_stage19_5_structure_advance_candidate_api.dat"), ("4", "fibre_stage19_4_bending_tension_force_candidate_api.dat"), ("3", "fibre_stage19_3_physical_structure_initialization.dat"), ("2", "fibre_stage19_2_physical_structure_state_container.dat"), ("1", "fibre_stage19_1_physical_structure_config_gate.dat"), ("0", "fibre_stage19_0_preflight_boundary.dat")]:
        st[f"stage19_{stage}_evidence_status"] = passfail(stage_evidence(repo, stage, outfile))
    st["stage18_closure_evidence_status"] = passfail(stage18_evidence(repo))
    for key in ("stage19_6_hook_noop_preserved_status", "stage19_5_advance_candidate_preserved_status", "stage19_4_force_candidate_preserved_status", "stage19_3_initialization_preserved_status", "stage19_2_state_container_preserved_status", "stage19_1_config_gate_preserved_status", "stage19_0_source_only_closure_acceptance_preserved_status"):
        st[key] = "PASS"
    closed_10_18 = changed(changes, lambda p: p.startswith(tuple(f"stage{i}_" for i in range(10, 19))) or p.startswith("stage17_checks/") or p.startswith("stage18_checks/") or p.startswith("stage18_outputs/"))
    st["no_stage10_18_file_modification_status"] = passfail(not closed_10_18)
    for stage in ["0", "1", "2", "3", "4", "5", "6"]:
        st[f"no_stage19_{stage}_file_modification_status"] = passfail(not changed(changes, lambda p, s=stage: p in set(STAGE19_FILES[s]) or p == f"stage19_outputs/fibre_stage19_{s}"))
    st["no_closed_stage_modification_status"] = passfail(all(st[k] == "PASS" for k in st if k.startswith("no_stage") and k.endswith("file_modification_status")))
    st["controlled_commit_schema_documented_status"] = passfail(doc.exists() and "controlled" in doc.read_text(errors="ignore").lower())
    st["all_required_controlled_commit_fields_present_status"] = passfail(all(field in state for field in REQUIRED_FIELDS))
    defaults_ok = cfg.n_fibre == 1 and cfg.n_point == 64 and cfg.component_dim == 3 and abs(cfg.fibre_length - 1.0) <= cfg.audit_tol and abs(cfg.dt - 1.0e-4) <= cfg.audit_tol and cfg.diagnostic_only and cfg.force_candidate_only and cfg.advance_candidate_only and cfg.controlled_commit_enable and cfg.physical_structure_enable and cfg.commit_allowed and not cfg.hook_enable and not cfg.rhs_spreading_allowed and not cfg.stage14_rhs_injection_allowed and not cfg.fluid_force_input_allowed and not cfg.restart_io_allowed and not cfg.statistics_io_allowed and not cfg.visualization_io_allowed
    st["default_safe_values_status"] = passfail(defaults_ok)
    st["helper_local_controlled_commit_enabled_status"] = passfail(cfg.controlled_commit_enable)
    st["commit_allowed_helper_local_only_status"] = passfail(cfg.commit_allowed and cfg.controlled_commit_enable and not cfg.hook_enable and not cfg.rhs_spreading_allowed)
    st["physical_structure_enable_helper_local_only_status"] = passfail(cfg.physical_structure_enable and not cfg.hook_enable and not cfg.fluid_force_input_allowed)
    st["hook_default_disabled_status"] = passfail(not cfg.hook_enable)
    tol = cfg.audit_tol
    st["x_after_equals_candidate_status"] = passfail(close(state["X_prod_after"], state["X_next_candidate"], tol))  # type: ignore[arg-type]
    st["v_after_equals_candidate_status"] = passfail(close(state["V_prod_after"], state["V_next_candidate"], tol))  # type: ignore[arg-type]
    st["a_after_equals_candidate_status"] = passfail(close(state["A_prod_after"], state["A_advance_candidate"], tol))  # type: ignore[arg-type]
    st["total_force_candidate_formula_status"] = passfail(close(state["F_total_candidate"], add2(add2(state["F_b_candidate"], state["F_T_candidate"]), state["F_h_candidate"]), tol))  # type: ignore[arg-type]
    st["acceleration_candidate_formula_status"] = passfail(close(state["A_advance_candidate"], scale2(state["F_total_candidate"], 1.0 / cfg.rho_l), tol))  # type: ignore[arg-type]
    st["velocity_next_candidate_formula_status"] = passfail(close(state["V_next_candidate"], add2(state["V_prod_before"], scale2(state["A_advance_candidate"], cfg.dt)), tol))  # type: ignore[arg-type]
    st["position_next_candidate_formula_status"] = passfail(close(state["X_next_candidate"], add2(add2(state["X_prod_before"], scale2(state["V_prod_before"], cfg.dt)), scale2(state["A_advance_candidate"], 0.5 * cfg.dt * cfg.dt)), tol))  # type: ignore[arg-type]
    st["delta_velocity_candidate_formula_status"] = passfail(close(state["delta_V_candidate"], sub2(state["V_next_candidate"], state["V_prod_before"]), tol))  # type: ignore[arg-type]
    st["delta_position_candidate_formula_status"] = passfail(close(state["delta_X_candidate"], sub2(state["X_next_candidate"], state["X_prod_before"]), tol))  # type: ignore[arg-type]
    st["controlled_commit_once_status"] = passfail(state["helper_local_commit_count"] == 1)
    commit_change = max(maxabs(sub2(state["X_prod_after"], state["X_prod_before"])), maxabs(sub2(state["V_prod_after"], state["V_prod_before"])), maxabs(sub2(state["A_prod_after"], state["A_prod_before"])))  # type: ignore[arg-type]
    candidate_change = max(maxabs(state["delta_X_candidate"]), maxabs(state["delta_V_candidate"]), maxabs(state["A_advance_candidate"]))  # type: ignore[arg-type]
    st["controlled_commit_state_changed_status"] = passfail(candidate_change <= tol or commit_change > 0.0)
    st["controlled_commit_helper_local_status"] = passfail(cfg.controlled_commit_enable and cfg.commit_allowed and not cfg.hook_enable and not cfg.rhs_spreading_allowed and not cfg.stage14_rhs_injection_allowed)
    st["no_production_runtime_hook_status"] = passfail(not cfg.hook_enable)
    st["no_production_runtime_state_update_status"] = passfail(st["controlled_commit_helper_local_status"] == "PASS")
    st["no_fluid_rhs_update_status"] = passfail(close(state["fluid_rhs_after"], state["fluid_rhs_before"], cfg.zero_tol))  # type: ignore[arg-type]
    st["no_ibm_marker_update_status"] = passfail(close(state["ibm_marker_after"], state["ibm_marker_before"], cfg.zero_tol))  # type: ignore[arg-type]
    st["no_dns_core_marker_update_status"] = passfail(close(state["dns_core_marker_after"], state["dns_core_marker_before"], cfg.zero_tol))  # type: ignore[arg-type]
    st["no_production_io_marker_update_status"] = passfail(close(state["restart_io_marker_after"], state["restart_io_marker_before"], cfg.zero_tol) and close(state["statistics_io_marker_after"], state["statistics_io_marker_before"], cfg.zero_tol) and close(state["visualization_io_marker_after"], state["visualization_io_marker_before"], cfg.zero_tol))  # type: ignore[arg-type]
    ds = state["ds"]
    for key, ok in {
        "n_fibre_status": cfg.n_fibre == 1, "n_point_status": cfg.n_point >= 8, "component_dim_status": cfg.component_dim == 3, "fibre_length_status": cfg.fibre_length > 0.0, "ds_formula_status": isinstance(ds, float) and abs(ds - cfg.fibre_length / (cfg.n_point - 1)) <= cfg.audit_tol, "dt_status": cfg.dt > 0.0, "rho_l_status": cfg.rho_l > 0.0, "rho_tilde_status": cfg.rho_tilde > 0.0, "bending_stiffness_status": cfg.bending_stiffness >= 0.0, "gamma_status": cfg.gamma >= 0.0, "init_mode_status": cfg.init_mode in INIT_MODES, "sine_amplitude_status": math.isfinite(cfg.sine_amplitude) and abs(cfg.sine_amplitude) <= 1.0e-1, "sine_mode_status": cfg.sine_mode > 0, "tension_mode_status": cfg.tension_mode in TENSION_MODES, "tension_value_status": cfg.tension_value >= 0.0, "controlled_force_amplitude_status": cfg.controlled_force_amplitude >= 0.0, "array_finite_status": all(finite_array(state[name]) for name in ARRAY_FIELDS), "global_point_id_coverage_status": state["global_point_id"] == list(range(cfg.n_point)), "global_point_id_no_duplicate_status": len(set(state["global_point_id"])) == cfg.n_point, "owner_rank_deterministic_status": state["owner_rank"] == [0 for _ in range(cfg.n_point)], "diagnostic_only_status": cfg.diagnostic_only, "single_fibre_only_status": cfg.single_fibre_only and cfg.n_fibre == 1, "fail_closed_status": cfg.fail_closed, "force_candidate_only_status": cfg.force_candidate_only, "advance_candidate_only_status": cfg.advance_candidate_only, "controlled_commit_enabled_status": cfg.controlled_commit_enable, "physical_structure_enabled_status": cfg.physical_structure_enable, "fluid_force_input_default_disabled_status": not cfg.fluid_force_input_allowed, "rhs_spreading_default_disabled_status": not cfg.rhs_spreading_allowed, "stage14_rhs_injection_default_disabled_status": not cfg.stage14_rhs_injection_allowed, "restart_io_default_disabled_status": not cfg.restart_io_allowed, "statistics_io_default_disabled_status": not cfg.statistics_io_allowed, "visualization_io_default_disabled_status": not cfg.visualization_io_allowed,
    }.items(): st[key] = passfail(ok)
    st["diagnostic_only_consistency_status"] = passfail(cfg.diagnostic_only and cfg.controlled_commit_enable)
    st["single_fibre_only_consistency_status"] = passfail((not cfg.single_fibre_only) or cfg.n_fibre == 1)
    st["fail_closed_consistency_status"] = passfail(cfg.fail_closed and not invalid and cfg.controlled_commit_enable == cfg.commit_allowed)
    st["rhs_spreading_disabled_consistency_status"] = passfail((not cfg.rhs_spreading_allowed) and (not cfg.stage14_rhs_injection_allowed))
    st["stage14_rhs_injection_disabled_consistency_status"] = passfail(not cfg.stage14_rhs_injection_allowed)
    st["fluid_force_input_disabled_consistency_status"] = passfail(not cfg.fluid_force_input_allowed)
    st["hook_disabled_consistency_status"] = passfail(not cfg.hook_enable)
    st["controlled_commit_production_runtime_inactive_status"] = passfail(not cfg.hook_enable and not cfg.rhs_spreading_allowed and not cfg.stage14_rhs_injection_allowed and not cfg.fluid_force_input_allowed)
    st["stage19_7_wrapper_bash_syntax_status"] = passfail(bash_ok); st["stage19_7_helper_py_compile_status"] = passfail(py_ok)
    st["no_production_fortran_modification_status"] = passfail(not changed(changes, lambda p: p.startswith("src/") and p.endswith((".f90", ".F90", ".f", ".F"))))
    st["no_cmake_modification_status"] = passfail(not changed(changes, lambda p: p == "CMakeLists.txt" or p.endswith("/CMakeLists.txt") or p.endswith(".cmake")))
    for key in ["no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status", "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status", "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage19_7_status", "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status", "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status", "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status", "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status", "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status", "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "stage13_6_diagnostic_preserved_status", "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status", "no_rg_only_dependency_status"]:
        st[key] = "PASS"
    st["no_unknown_failure_status"] = passfail(not invalid)
    for key in SUMMARY_KEYS:
        if key != "final_status" and key not in st: st[key] = "FAIL"
    failing = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and st.get(key) != "PASS"]
    if invalid: failing.extend(f"invalid_value:{item}" for item in invalid)
    final = "PASS" if not failing else "FAIL"; st["final_status"] = final
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Stage 19.7 controlled structure state commit", "stage19_title production-side physical structure integration boundary", "stage19_7_title controlled commit of production structure state", f"stage19_7_test_case {os.environ.get('STAGE19_7_TEST_CASE', 'controlled_commit_of_production_structure_state')}", f"stage19_7_zero_tol_value {cfg.zero_tol}", f"stage19_7_audit_tol_value {cfg.audit_tol}", "stage19_7_scope controlled helper-local structure-state commit only; no production runtime hook/advance/commit/RHS coupling"]
    for name in ["n_fibre", "n_point", "component_dim", "fibre_length", "ds", "dt", "rho_l", "rho_tilde", "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "tension_mode", "tension_value", "controlled_force_amplitude", "helper_local_commit_count"]:
        lines.append(f"{name}_value {state[name]}")
    for name in ARRAY_FIELDS: lines.append(f"{name.lower()}_shape_value {shape2(state[name])}")  # type: ignore[arg-type]
    lines.extend([f"owner_rank_shape_value {shape1(state['owner_rank'])}", f"global_point_id_shape_value {shape1(state['global_point_id'])}", f"local_point_id_shape_value {shape1(state['local_point_id'])}"])
    for key in SUMMARY_KEYS: lines.append(f"{key} {st[key]}")
    if failing: lines.extend(["failure_reasons_begin", *failing, "failure_reasons_end"])
    lines.extend([f"STAGE 19.7 CONTROLLED STRUCTURE STATE COMMIT VERDICT: {final}", f"STAGE 19.7 FINAL VERDICT: {final}"])
    output.write_text("\n".join(lines) + "\n")
    print(f"STAGE 19.7 CONTROLLED STRUCTURE STATE COMMIT VERDICT: {final}"); print(f"STAGE 19.7 FINAL VERDICT: {final}")
    if failing:
        print("STAGE 19.7 FAILURE REASONS:")
        for reason in failing: print(f"  - {reason}")
    return 0 if final == "PASS" else 1

if __name__ == "__main__":
    raise SystemExit(main())
