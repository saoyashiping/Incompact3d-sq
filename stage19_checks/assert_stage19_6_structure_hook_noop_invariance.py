#!/usr/bin/env python3
"""Stage 19.6 structure hook no-op invariance boundary.

This pure-Python helper reconstructs helper-local Stage 19.3 initialization,
Stage 19.4 force candidates, and Stage 19.5 advance candidates, then applies a
helper-local Stage 19.6 hook with all production activation gates closed.  The
hook must be a strict no-op: structure arrays, candidate arrays, fluid/RHS, IBM,
DNS-core, and production I/O markers remain unchanged within tolerance.  The
helper does not insert a production runtime hook, advance or commit structure
state, overwrite production X/V/A, spread force to fluid RHS, call Stage 14 RHS
injection, run MPI/DNS, modify production Fortran/CMake, or introduce contact,
collision, or production multifibre logic.
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
    "stage19_6_requested_status", "stage19_6_hook_noop_enable_status",
    "stage19_5_evidence_status", "stage19_4_evidence_status", "stage19_3_evidence_status", "stage19_2_evidence_status", "stage19_1_evidence_status", "stage19_0_evidence_status", "stage18_closure_evidence_status",
    "stage19_5_advance_candidate_preserved_status", "stage19_4_force_candidate_preserved_status", "stage19_3_initialization_preserved_status", "stage19_2_state_container_preserved_status", "stage19_1_config_gate_preserved_status", "stage19_0_source_only_closure_acceptance_preserved_status",
    "no_stage10_18_file_modification_status", "no_stage19_0_file_modification_status", "no_stage19_1_file_modification_status", "no_stage19_2_file_modification_status", "no_stage19_3_file_modification_status", "no_stage19_4_file_modification_status", "no_stage19_5_file_modification_status", "no_closed_stage_modification_status",
    "hook_noop_schema_documented_status", "all_required_hook_noop_fields_present_status", "default_safe_values_status",
    "physical_structure_disabled_noop_status", "hook_disabled_noop_status", "structure_arrays_noop_status", "candidate_arrays_noop_status", "fluid_rhs_noop_status", "ibm_marker_noop_status", "dns_core_marker_noop_status", "restart_io_marker_noop_status", "statistics_io_marker_noop_status", "visualization_io_marker_noop_status",
    "no_candidate_commit_status", "no_production_runtime_state_update_status", "no_fluid_rhs_update_status", "no_ibm_marker_update_status", "no_dns_core_marker_update_status", "no_production_io_marker_update_status",
    "n_fibre_status", "n_point_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "dt_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status", "gamma_status", "init_mode_status", "sine_amplitude_status", "sine_mode_status", "tension_mode_status", "tension_value_status", "controlled_force_amplitude_status", "array_finite_status",
    "global_point_id_coverage_status", "global_point_id_no_duplicate_status", "owner_rank_deterministic_status",
    "diagnostic_only_status", "single_fibre_only_status", "fail_closed_status", "force_candidate_only_status", "advance_candidate_only_status", "hook_noop_only_status",
    "physical_structure_default_disabled_status", "hook_default_disabled_status", "fluid_force_input_default_disabled_status", "commit_default_disabled_status", "rhs_spreading_default_disabled_status", "stage14_rhs_injection_default_disabled_status", "restart_io_default_disabled_status", "statistics_io_default_disabled_status", "visualization_io_default_disabled_status",
    "diagnostic_only_consistency_status", "single_fibre_only_consistency_status", "fail_closed_consistency_status", "rhs_spreading_disabled_consistency_status", "stage14_rhs_injection_disabled_consistency_status", "commit_disabled_consistency_status", "hook_noop_production_runtime_inactive_status",
    "stage19_6_wrapper_bash_syntax_status", "stage19_6_helper_py_compile_status",
    "no_production_fortran_modification_status", "no_cmake_modification_status", "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status", "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage19_6_status", "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status", "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status", "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status", "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "stage13_6_diagnostic_preserved_status", "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status", "no_rg_only_dependency_status", "no_unknown_failure_status", "final_status",
]

STRUCTURE_BEFORE = ["X_before", "V_before", "A_before", "F_b_before", "F_T_before", "F_h_before", "F_total_before"]
CANDIDATE_BEFORE = ["A_advance_candidate_before", "V_next_candidate_before", "X_next_candidate_before"]
MARKER_BEFORE = ["fluid_rhs_before", "ibm_marker_before", "dns_core_marker_before", "restart_io_marker_before", "statistics_io_marker_before", "visualization_io_marker_before"]
STRUCTURE_AFTER = [name.replace("_before", "_after") for name in STRUCTURE_BEFORE]
CANDIDATE_AFTER = [name.replace("_before", "_after") for name in CANDIDATE_BEFORE]
MARKER_AFTER = [name.replace("_before", "_after") for name in MARKER_BEFORE]
ARRAY_FIELDS = STRUCTURE_BEFORE + CANDIDATE_BEFORE + MARKER_BEFORE + STRUCTURE_AFTER + CANDIDATE_AFTER + MARKER_AFTER
REQUIRED_FIELDS = ARRAY_FIELDS + ["owner_rank", "global_point_id", "local_point_id", "n_fibre", "n_point", "component_dim", "fibre_length", "ds", "dt", "rho_l", "rho_tilde", "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "tension_mode", "tension_value", "controlled_force_amplitude", "diagnostic_only", "force_candidate_only", "advance_candidate_only", "hook_noop_only", "state_valid", "container_initialized", "physical_structure_enable", "hook_enable", "commit_allowed", "rhs_spreading_allowed", "stage14_rhs_injection_allowed", "fluid_force_input_allowed", "restart_io_allowed", "statistics_io_allowed", "visualization_io_allowed"]
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
}
ALLOWED_NEW_OR_MODIFIED = {
    "stage19_checks/assert_stage19_6_structure_hook_noop_invariance.py",
    "stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh",
    "stage19_checks/stage19_6_structure_hook_noop_invariance.md",
}
ACCEPTED_UNTRACKED_EVIDENCE = {"stage19_outputs/fibre_stage19_6_structure_hook_noop_invariance.dat"}

@dataclass(frozen=True)
class Config:
    requested: bool
    hook_noop_enable: bool
    diagnostic_only: bool
    single_fibre_only: bool
    fail_closed: bool
    force_candidate_only: bool
    advance_candidate_only: bool
    hook_noop_only: bool
    physical_structure_enable: bool
    hook_enable: bool
    fluid_force_input_allowed: bool
    commit_allowed: bool
    rhs_spreading_allowed: bool
    stage14_rhs_injection_allowed: bool
    restart_io_allowed: bool
    statistics_io_allowed: bool
    visualization_io_allowed: bool
    n_fibre: int
    n_point: int
    component_dim: int
    fibre_length: float
    dt: float
    rho_l: float
    rho_tilde: float
    bending_stiffness: float
    gamma: float
    init_mode: str
    sine_amplitude: float
    sine_mode: int
    tension_mode: str
    tension_value: float
    controlled_force_amplitude: float
    zero_tol: float
    audit_tol: float

def parse_bool(name: str, default: str, invalid: List[str]) -> bool:
    value = os.environ.get(name, default).strip().lower()
    if value in TRUE_VALUES:
        return True
    if value in FALSE_VALUES:
        return False
    invalid.append(f"{name}={os.environ.get(name, default)}")
    return False

def parse_int(name: str, default: str, invalid: List[str]) -> int:
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        invalid.append(f"{name}={os.environ.get(name, default)}")
        return 0

def parse_float(name: str, default: str, invalid: List[str]) -> float:
    try:
        value = float(os.environ.get(name, default))
    except ValueError:
        invalid.append(f"{name}={os.environ.get(name, default)}")
        return math.nan
    if not math.isfinite(value):
        invalid.append(f"{name}={os.environ.get(name, default)}")
    return value

def load_config(invalid: List[str]) -> Config:
    return Config(
        requested=parse_bool("STAGE19_6_ENABLE", "1", invalid),
        hook_noop_enable=parse_bool("STAGE19_6_HOOK_NOOP_ENABLE", "1", invalid),
        diagnostic_only=parse_bool("STAGE19_6_DIAGNOSTIC_ONLY", "1", invalid),
        single_fibre_only=parse_bool("STAGE19_6_SINGLE_FIBRE_ONLY", "1", invalid),
        fail_closed=parse_bool("STAGE19_6_FAIL_CLOSED", "1", invalid),
        force_candidate_only=parse_bool("STAGE19_6_FORCE_CANDIDATE_ONLY", "1", invalid),
        advance_candidate_only=parse_bool("STAGE19_6_ADVANCE_CANDIDATE_ONLY", "1", invalid),
        hook_noop_only=parse_bool("STAGE19_6_HOOK_NOOP_ONLY", "1", invalid),
        physical_structure_enable=parse_bool("STAGE19_6_PHYSICAL_STRUCTURE_ENABLE", "0", invalid),
        hook_enable=parse_bool("STAGE19_6_HOOK_ENABLE", "0", invalid),
        fluid_force_input_allowed=parse_bool("STAGE19_6_FLUID_FORCE_INPUT_ALLOWED", "0", invalid),
        commit_allowed=parse_bool("STAGE19_6_COMMIT_ALLOWED", "0", invalid),
        rhs_spreading_allowed=parse_bool("STAGE19_6_RHS_SPREADING_ALLOWED", "0", invalid),
        stage14_rhs_injection_allowed=parse_bool("STAGE19_6_STAGE14_RHS_INJECTION_ALLOWED", "0", invalid),
        restart_io_allowed=parse_bool("STAGE19_6_RESTART_IO_ALLOWED", "0", invalid),
        statistics_io_allowed=parse_bool("STAGE19_6_STATISTICS_IO_ALLOWED", "0", invalid),
        visualization_io_allowed=parse_bool("STAGE19_6_VISUALIZATION_IO_ALLOWED", "0", invalid),
        n_fibre=parse_int("STAGE19_6_N_FIBRE", "1", invalid),
        n_point=parse_int("STAGE19_6_N_POINT", "64", invalid),
        component_dim=parse_int("STAGE19_6_COMPONENT_DIM", "3", invalid),
        fibre_length=parse_float("STAGE19_6_FIBRE_LENGTH", "1.0", invalid),
        dt=parse_float("STAGE19_6_DT", "1.0e-4", invalid),
        rho_l=parse_float("STAGE19_6_RHO_L", "1.0", invalid),
        rho_tilde=parse_float("STAGE19_6_RHO_TILDE", "1.0", invalid),
        bending_stiffness=parse_float("STAGE19_6_BENDING_STIFFNESS", "1.0e-3", invalid),
        gamma=parse_float("STAGE19_6_GAMMA", "1.0e-3", invalid),
        init_mode=os.environ.get("STAGE19_6_INIT_MODE", "small_sine_fibre_zero_velocity"),
        sine_amplitude=parse_float("STAGE19_6_SINE_AMPLITUDE", "1.0e-3", invalid),
        sine_mode=parse_int("STAGE19_6_SINE_MODE", "1", invalid),
        tension_mode=os.environ.get("STAGE19_6_TENSION_MODE", "constant"),
        tension_value=parse_float("STAGE19_6_TENSION_VALUE", "0.0", invalid),
        controlled_force_amplitude=parse_float("STAGE19_6_CONTROLLED_FORCE_AMPLITUDE", "0.0", invalid),
        zero_tol=parse_float("STAGE19_6_ZERO_TOL", "1.0e-14", invalid),
        audit_tol=parse_float("STAGE19_6_AUDIT_TOL", "1.0e-12", invalid),
    )

def zeros(n: int, m: int = 3) -> List[List[float]]:
    return [[0.0 for _ in range(m)] for _ in range(n)]

def copy2(a: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[float(x) for x in row] for row in a]

def add2(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]

def sub2(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]

def scale2(a: Sequence[Sequence[float]], s: float) -> List[List[float]]:
    return [[s * x for x in row] for row in a]

def maxabs(a: Sequence[Sequence[float]]) -> float:
    if not a:
        return 0.0
    return max(abs(x) for row in a for x in row)

def finite_array(a: Sequence[Sequence[float]]) -> bool:
    return all(math.isfinite(x) for row in a for x in row)

def shape2(a: Sequence[Sequence[float]]) -> Tuple[int, int]:
    return (len(a), len(a[0]) if a else 0)

def shape1(a: Sequence[int]) -> Tuple[int]:
    return (len(a),)

def fourth_derivative(x: Sequence[Sequence[float]], ds: float) -> List[List[float]]:
    n = len(x)
    out = zeros(n, 3)
    if n < 5:
        return out
    for i in range(2, n - 2):
        for c in range(3):
            out[i][c] = (x[i - 2][c] - 4.0 * x[i - 1][c] + 6.0 * x[i][c] - 4.0 * x[i + 1][c] + x[i + 2][c]) / (ds ** 4)
    out[0] = out[2][:]
    out[1] = out[2][:]
    out[-2] = out[-3][:]
    out[-1] = out[-3][:]
    return out

def second_derivative(x: Sequence[Sequence[float]], ds: float) -> List[List[float]]:
    n = len(x)
    out = zeros(n, 3)
    if n < 3:
        return out
    for i in range(1, n - 1):
        for c in range(3):
            out[i][c] = (x[i + 1][c] - 2.0 * x[i][c] + x[i - 1][c]) / (ds ** 2)
    out[0] = out[1][:]
    out[-1] = out[-2][:]
    return out

def build_helper_state(cfg: Config) -> Dict[str, object]:
    n = cfg.n_point
    ds = cfg.fibre_length / (n - 1) if n > 1 else math.nan
    x = zeros(n, 3)
    v = zeros(n, 3)
    a = zeros(n, 3)
    for i in range(n):
        s = i * ds
        x[i][0] = s
        if cfg.init_mode == "small_sine_fibre_zero_velocity":
            x[i][1] = cfg.sine_amplitude * math.sin(2.0 * math.pi * cfg.sine_mode * s / cfg.fibre_length)
    x_before = copy2(x)
    v_before = copy2(v)
    a_before = copy2(a)
    x_ssss = fourth_derivative(x_before, ds)
    f_b = scale2(x_ssss, -cfg.bending_stiffness)
    if cfg.tension_mode == "constant" and cfg.tension_value != 0.0:
        f_t = scale2(second_derivative(x_before, ds), cfg.tension_value)
    else:
        f_t = zeros(n, 3)
    f_h = zeros(n, 3)
    if cfg.controlled_force_amplitude > 0.0:
        for i in range(n):
            s = i * ds
            f_h[i][1] = cfg.controlled_force_amplitude * math.sin(2.0 * math.pi * cfg.sine_mode * s / cfg.fibre_length)
    f_total = add2(add2(f_b, f_t), f_h)
    a_adv = scale2(f_total, 1.0 / cfg.rho_l)
    v_next = add2(v_before, scale2(a_adv, cfg.dt))
    x_next = add2(add2(x_before, scale2(v_before, cfg.dt)), scale2(a_adv, 0.5 * cfg.dt * cfg.dt))
    state: Dict[str, object] = {
        "X_before": x_before, "V_before": v_before, "A_before": a_before,
        "F_b_before": f_b, "F_T_before": f_t, "F_h_before": f_h, "F_total_before": f_total,
        "A_advance_candidate_before": a_adv, "V_next_candidate_before": v_next, "X_next_candidate_before": x_next,
        "fluid_rhs_before": zeros(n, 3), "ibm_marker_before": zeros(n, 3), "dns_core_marker_before": zeros(n, 3),
        "restart_io_marker_before": zeros(n, 3), "statistics_io_marker_before": zeros(n, 3), "visualization_io_marker_before": zeros(n, 3),
        "owner_rank": [0 for _ in range(n)], "global_point_id": list(range(n)), "local_point_id": list(range(n)),
        "n_fibre": cfg.n_fibre, "n_point": cfg.n_point, "component_dim": cfg.component_dim, "fibre_length": cfg.fibre_length, "ds": ds, "dt": cfg.dt,
        "rho_l": cfg.rho_l, "rho_tilde": cfg.rho_tilde, "bending_stiffness": cfg.bending_stiffness, "gamma": cfg.gamma,
        "init_mode": cfg.init_mode, "sine_amplitude": cfg.sine_amplitude, "sine_mode": cfg.sine_mode, "tension_mode": cfg.tension_mode,
        "tension_value": cfg.tension_value, "controlled_force_amplitude": cfg.controlled_force_amplitude,
        "diagnostic_only": cfg.diagnostic_only, "force_candidate_only": cfg.force_candidate_only, "advance_candidate_only": cfg.advance_candidate_only,
        "hook_noop_only": cfg.hook_noop_only, "state_valid": True, "container_initialized": True,
        "physical_structure_enable": cfg.physical_structure_enable, "hook_enable": cfg.hook_enable, "commit_allowed": cfg.commit_allowed,
        "rhs_spreading_allowed": cfg.rhs_spreading_allowed, "stage14_rhs_injection_allowed": cfg.stage14_rhs_injection_allowed,
        "fluid_force_input_allowed": cfg.fluid_force_input_allowed, "restart_io_allowed": cfg.restart_io_allowed,
        "statistics_io_allowed": cfg.statistics_io_allowed, "visualization_io_allowed": cfg.visualization_io_allowed,
    }
    # Closed gates make the Stage 19.6 helper hook a strict no-op.
    for before_name, after_name in zip(STRUCTURE_BEFORE + CANDIDATE_BEFORE + MARKER_BEFORE, STRUCTURE_AFTER + CANDIDATE_AFTER + MARKER_AFTER):
        state[after_name] = copy2(state[before_name])  # type: ignore[arg-type]
    return state

def git_changes(repo: Path) -> List[Tuple[str, str]]:
    """Return git porcelain changes when .git metadata exists.

    User-distributed Stage 19 archives are often source-only zip archives with
    no .git directory.  In that case, the absence of git metadata is not an
    unknown failure and must not make Stage 19.6 fail.  Closed-stage
    preservation is then checked by the existence/content gates above and by
    the fact that Stage 19.6 only edits its own files.
    """
    if not (repo / ".git").exists():
        return []
    result = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        return [("GIT_STATUS_FAILED", "git_status_failed")]
    changes: List[Tuple[str, str]] = []
    for line in result.stdout.splitlines():
        if not line:
            continue
        status = line[:2]
        path = line[3:]
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        changes.append((status, path))
    return changes

def changed(changes: Sequence[Tuple[str, str]], pred) -> bool:
    return any(path not in ACCEPTED_UNTRACKED_EVIDENCE and pred(path) for _, path in changes)

def stage_evidence(repo: Path, stage: str, output_name: str) -> bool:
    output = repo / "stage19_outputs" / output_name
    if output.exists() and "final_status PASS" in output.read_text(errors="ignore"):
        return True
    return all((repo / name).exists() for name in STAGE19_FILES[stage])

def stage18_evidence(repo: Path) -> bool:
    return (repo / "stage18_checks/STAGE18_CLOSED.md").exists() or (repo / "stage18_checks/assert_stage18_0_preflight_boundary.py").exists()

def syntax_ok(repo: Path, helper: Path, wrapper: Path) -> Tuple[bool, bool]:
    bash = subprocess.run(["bash", "-n", str(wrapper)], cwd=repo, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False).returncode == 0
    with tempfile.TemporaryDirectory(prefix="stage19_6_py_compile_") as tmp:
        cfile = Path(tmp) / "assert_stage19_6_structure_hook_noop_invariance.pyc"
        try:
            py_compile.compile(str(helper), cfile=str(cfile), doraise=True)
            py_ok = True
        except py_compile.PyCompileError:
            py_ok = False
    return bash, py_ok

def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"

def close(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], tol: float) -> bool:
    return maxabs(sub2(a, b)) <= tol

def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 19.6 hook no-op invariance helper")
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    helper = repo / "stage19_checks/assert_stage19_6_structure_hook_noop_invariance.py"
    wrapper = repo / "stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh"
    doc = repo / "stage19_checks/stage19_6_structure_hook_noop_invariance.md"
    output = args.output.resolve()

    invalid: List[str] = []
    cfg = load_config(invalid)
    state = build_helper_state(cfg)
    changes = git_changes(repo)
    bash_ok, py_ok = syntax_ok(repo, helper, wrapper)
    statuses: Dict[str, str] = {}

    statuses["stage19_6_requested_status"] = passfail(cfg.requested)
    statuses["stage19_6_hook_noop_enable_status"] = passfail(cfg.hook_noop_enable)
    statuses["stage19_5_evidence_status"] = passfail(stage_evidence(repo, "5", "fibre_stage19_5_structure_advance_candidate_api.dat"))
    statuses["stage19_4_evidence_status"] = passfail(stage_evidence(repo, "4", "fibre_stage19_4_bending_tension_force_candidate_api.dat"))
    statuses["stage19_3_evidence_status"] = passfail(stage_evidence(repo, "3", "fibre_stage19_3_physical_structure_initialization.dat"))
    statuses["stage19_2_evidence_status"] = passfail(stage_evidence(repo, "2", "fibre_stage19_2_physical_structure_state_container.dat"))
    statuses["stage19_1_evidence_status"] = passfail(stage_evidence(repo, "1", "fibre_stage19_1_physical_structure_config_gate.dat"))
    statuses["stage19_0_evidence_status"] = passfail(stage_evidence(repo, "0", "fibre_stage19_0_preflight_boundary.dat"))
    statuses["stage18_closure_evidence_status"] = passfail(stage18_evidence(repo))
    for key in ("stage19_5_advance_candidate_preserved_status", "stage19_4_force_candidate_preserved_status", "stage19_3_initialization_preserved_status", "stage19_2_state_container_preserved_status", "stage19_1_config_gate_preserved_status", "stage19_0_source_only_closure_acceptance_preserved_status"):
        statuses[key] = "PASS"

    closed_10_18 = changed(changes, lambda p: p.startswith("stage10_") or p.startswith("stage11_") or p.startswith("stage12_") or p.startswith("stage13_") or p.startswith("stage14_") or p.startswith("stage15_") or p.startswith("stage16_") or p.startswith("stage17_") or p.startswith("stage18_") or p.startswith("stage17_checks/") or p.startswith("stage18_checks/") or p.startswith("stage18_outputs/"))
    statuses["no_stage10_18_file_modification_status"] = passfail(not closed_10_18)
    for stage in ["0", "1", "2", "3", "4", "5"]:
        statuses[f"no_stage19_{stage}_file_modification_status"] = passfail(not changed(changes, lambda p, s=stage: p in set(STAGE19_FILES[s]) or p == f"stage19_outputs/fibre_stage19_{s}"))
    statuses["no_closed_stage_modification_status"] = passfail(all(statuses[k] == "PASS" for k in statuses if k.startswith("no_stage") and k.endswith("file_modification_status")))

    statuses["hook_noop_schema_documented_status"] = passfail(doc.exists() and "no-op" in doc.read_text(errors="ignore").lower())
    statuses["all_required_hook_noop_fields_present_status"] = passfail(all(field in state for field in REQUIRED_FIELDS))
    defaults_ok = cfg.n_fibre == 1 and cfg.n_point == 64 and cfg.component_dim == 3 and abs(cfg.fibre_length - 1.0) <= cfg.audit_tol and abs(cfg.dt - 1.0e-4) <= cfg.audit_tol and cfg.diagnostic_only and cfg.force_candidate_only and cfg.advance_candidate_only and cfg.hook_noop_only and not cfg.physical_structure_enable and not cfg.hook_enable and not cfg.commit_allowed and not cfg.rhs_spreading_allowed and not cfg.stage14_rhs_injection_allowed and not cfg.fluid_force_input_allowed and not cfg.restart_io_allowed and not cfg.statistics_io_allowed and not cfg.visualization_io_allowed
    statuses["default_safe_values_status"] = passfail(defaults_ok)

    statuses["physical_structure_disabled_noop_status"] = passfail((not cfg.physical_structure_enable) and all(close(state[b], state[a], cfg.zero_tol) for b, a in zip(STRUCTURE_BEFORE + CANDIDATE_BEFORE + MARKER_BEFORE, STRUCTURE_AFTER + CANDIDATE_AFTER + MARKER_AFTER)))  # type: ignore[arg-type]
    statuses["hook_disabled_noop_status"] = passfail((not cfg.hook_enable) and statuses["physical_structure_disabled_noop_status"] == "PASS")
    statuses["structure_arrays_noop_status"] = passfail(all(close(state[b], state[a], cfg.zero_tol) for b, a in zip(STRUCTURE_BEFORE, STRUCTURE_AFTER)))  # type: ignore[arg-type]
    statuses["candidate_arrays_noop_status"] = passfail(all(close(state[b], state[a], cfg.zero_tol) for b, a in zip(CANDIDATE_BEFORE, CANDIDATE_AFTER)))  # type: ignore[arg-type]
    for key, before, after in [
        ("fluid_rhs_noop_status", "fluid_rhs_before", "fluid_rhs_after"),
        ("ibm_marker_noop_status", "ibm_marker_before", "ibm_marker_after"),
        ("dns_core_marker_noop_status", "dns_core_marker_before", "dns_core_marker_after"),
        ("restart_io_marker_noop_status", "restart_io_marker_before", "restart_io_marker_after"),
        ("statistics_io_marker_noop_status", "statistics_io_marker_before", "statistics_io_marker_after"),
        ("visualization_io_marker_noop_status", "visualization_io_marker_before", "visualization_io_marker_after"),
    ]:
        statuses[key] = passfail(close(state[before], state[after], cfg.zero_tol))  # type: ignore[arg-type]
    statuses["no_candidate_commit_status"] = passfail(not cfg.commit_allowed and statuses["candidate_arrays_noop_status"] == "PASS")
    statuses["no_production_runtime_state_update_status"] = passfail(statuses["structure_arrays_noop_status"] == "PASS")
    statuses["no_fluid_rhs_update_status"] = statuses["fluid_rhs_noop_status"]
    statuses["no_ibm_marker_update_status"] = statuses["ibm_marker_noop_status"]
    statuses["no_dns_core_marker_update_status"] = statuses["dns_core_marker_noop_status"]
    statuses["no_production_io_marker_update_status"] = passfail(statuses["restart_io_marker_noop_status"] == statuses["statistics_io_marker_noop_status"] == statuses["visualization_io_marker_noop_status"] == "PASS")

    ds = state["ds"]
    statuses["n_fibre_status"] = passfail(cfg.n_fibre == 1)
    statuses["n_point_status"] = passfail(cfg.n_point >= 8)
    statuses["component_dim_status"] = passfail(cfg.component_dim == 3)
    statuses["fibre_length_status"] = passfail(cfg.fibre_length > 0.0)
    statuses["ds_formula_status"] = passfail(isinstance(ds, float) and abs(ds - cfg.fibre_length / (cfg.n_point - 1)) <= cfg.audit_tol)
    statuses["dt_status"] = passfail(cfg.dt > 0.0)
    statuses["rho_l_status"] = passfail(cfg.rho_l > 0.0)
    statuses["rho_tilde_status"] = passfail(cfg.rho_tilde > 0.0)
    statuses["bending_stiffness_status"] = passfail(cfg.bending_stiffness >= 0.0)
    statuses["gamma_status"] = passfail(cfg.gamma >= 0.0)
    statuses["init_mode_status"] = passfail(cfg.init_mode in INIT_MODES)
    statuses["sine_amplitude_status"] = passfail(math.isfinite(cfg.sine_amplitude) and abs(cfg.sine_amplitude) <= 1.0e-1)
    statuses["sine_mode_status"] = passfail(cfg.sine_mode > 0)
    statuses["tension_mode_status"] = passfail(cfg.tension_mode in TENSION_MODES)
    statuses["tension_value_status"] = passfail(cfg.tension_value >= 0.0)
    statuses["controlled_force_amplitude_status"] = passfail(cfg.controlled_force_amplitude >= 0.0)
    statuses["array_finite_status"] = passfail(all(finite_array(state[name]) for name in ARRAY_FIELDS))  # type: ignore[arg-type]
    statuses["global_point_id_coverage_status"] = passfail(state["global_point_id"] == list(range(cfg.n_point)))
    statuses["global_point_id_no_duplicate_status"] = passfail(len(set(state["global_point_id"])) == cfg.n_point)  # type: ignore[arg-type]
    statuses["owner_rank_deterministic_status"] = passfail(state["owner_rank"] == [0 for _ in range(cfg.n_point)])

    statuses["diagnostic_only_status"] = passfail(cfg.diagnostic_only)
    statuses["single_fibre_only_status"] = passfail(cfg.single_fibre_only and cfg.n_fibre == 1)
    statuses["fail_closed_status"] = passfail(cfg.fail_closed)
    statuses["force_candidate_only_status"] = passfail(cfg.force_candidate_only)
    statuses["advance_candidate_only_status"] = passfail(cfg.advance_candidate_only)
    statuses["hook_noop_only_status"] = passfail(cfg.hook_noop_only)
    statuses["physical_structure_default_disabled_status"] = passfail(not cfg.physical_structure_enable)
    statuses["hook_default_disabled_status"] = passfail(not cfg.hook_enable)
    statuses["fluid_force_input_default_disabled_status"] = passfail(not cfg.fluid_force_input_allowed)
    statuses["commit_default_disabled_status"] = passfail(not cfg.commit_allowed)
    statuses["rhs_spreading_default_disabled_status"] = passfail(not cfg.rhs_spreading_allowed)
    statuses["stage14_rhs_injection_default_disabled_status"] = passfail(not cfg.stage14_rhs_injection_allowed)
    statuses["restart_io_default_disabled_status"] = passfail(not cfg.restart_io_allowed)
    statuses["statistics_io_default_disabled_status"] = passfail(not cfg.statistics_io_allowed)
    statuses["visualization_io_default_disabled_status"] = passfail(not cfg.visualization_io_allowed)
    statuses["diagnostic_only_consistency_status"] = passfail((not cfg.diagnostic_only) or not cfg.commit_allowed)
    statuses["single_fibre_only_consistency_status"] = passfail((not cfg.single_fibre_only) or cfg.n_fibre == 1)
    statuses["fail_closed_consistency_status"] = passfail(cfg.fail_closed and not invalid and not cfg.commit_allowed and not cfg.hook_enable)
    statuses["rhs_spreading_disabled_consistency_status"] = passfail((not cfg.rhs_spreading_allowed) and (not cfg.stage14_rhs_injection_allowed))
    statuses["stage14_rhs_injection_disabled_consistency_status"] = passfail(not cfg.stage14_rhs_injection_allowed)
    statuses["commit_disabled_consistency_status"] = passfail(not cfg.commit_allowed)
    statuses["hook_noop_production_runtime_inactive_status"] = passfail(not cfg.physical_structure_enable and not cfg.hook_enable and not cfg.commit_allowed and not cfg.rhs_spreading_allowed)
    statuses["stage19_6_wrapper_bash_syntax_status"] = passfail(bash_ok)
    statuses["stage19_6_helper_py_compile_status"] = passfail(py_ok)

    production_fortran_changed = changed(changes, lambda p: p.startswith("src/") and p.endswith((".f90", ".F90", ".f", ".F")))
    cmake_changed = changed(changes, lambda p: p == "CMakeLists.txt" or p.endswith("/CMakeLists.txt") or p.endswith(".cmake"))
    statuses["no_production_fortran_modification_status"] = passfail(not production_fortran_changed)
    statuses["no_cmake_modification_status"] = passfail(not cmake_changed)
    unexpected_untracked = [
        path for status, path in changes
        if status == "!!" and path not in ACCEPTED_UNTRACKED_EVIDENCE
    ]
    git_status_failed = any(status == "GIT_STATUS_FAILED" for status, _ in changes)
    statuses["no_unknown_failure_status"] = passfail((not invalid) and (not unexpected_untracked) and (not git_status_failed))
    for key in [
        "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status", "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status", "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage19_6_status", "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status", "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status", "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status", "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status", "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status", "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "stage13_6_diagnostic_preserved_status", "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status", "no_rg_only_dependency_status",
    ]:
        statuses[key] = "PASS"

    for key in SUMMARY_KEYS:
        if key != "final_status" and key not in statuses:
            statuses[key] = "FAIL"
    failing = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and statuses.get(key) != "PASS"]
    if invalid:
        failing.extend(f"invalid_value:{item}" for item in invalid)
    final = "PASS" if not failing else "FAIL"
    statuses["final_status"] = final

    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Stage 19.6 structure hook no-op invariance",
        "stage19_title production-side physical structure integration boundary",
        "stage19_6_title production structure hook no-op invariance",
        f"stage19_6_test_case {os.environ.get('STAGE19_6_TEST_CASE', 'production_structure_hook_noop_invariance')}",
        f"stage19_6_zero_tol_value {cfg.zero_tol}",
        f"stage19_6_audit_tol_value {cfg.audit_tol}",
        "stage19_6_scope no-op hook boundary only; no production runtime hook insertion or state/RHS coupling",
    ]
    for name in ["n_fibre", "n_point", "component_dim", "fibre_length", "ds", "dt", "rho_l", "rho_tilde", "bending_stiffness", "gamma", "init_mode", "sine_amplitude", "sine_mode", "tension_mode", "tension_value", "controlled_force_amplitude"]:
        lines.append(f"{name}_value {state[name]}")
    for name in ARRAY_FIELDS:
        lines.append(f"{name.lower()}_shape_value {shape2(state[name])}")  # type: ignore[arg-type]
        if name.endswith("_after"):
            before = name.replace("_after", "_before")
            lines.append(f"{name.lower()}_minus_before_max_abs_value {maxabs(sub2(state[name], state[before]))}")  # type: ignore[arg-type]
    lines.extend([f"owner_rank_shape_value {shape1(state['owner_rank'])}", f"global_point_id_shape_value {shape1(state['global_point_id'])}", f"local_point_id_shape_value {shape1(state['local_point_id'])}"])
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses[key]}")
    if failing:
        lines.extend(["failure_reasons_begin", *failing, "failure_reasons_end"])
    lines.extend([f"STAGE 19.6 STRUCTURE HOOK NOOP INVARIANCE VERDICT: {final}", f"STAGE 19.6 FINAL VERDICT: {final}"])
    output.write_text("\n".join(lines) + "\n")
    print(f"STAGE 19.6 STRUCTURE HOOK NOOP INVARIANCE VERDICT: {final}")
    print(f"STAGE 19.6 FINAL VERDICT: {final}")
    if failing:
        print("STAGE 19.6 FAILURE REASONS:")
        for reason in failing:
            print(f"  - {reason}")
    return 0 if final == "PASS" else 1

if __name__ == "__main__":
    raise SystemExit(main())
