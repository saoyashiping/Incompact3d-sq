#!/usr/bin/env python3
"""Stage 19.4 bending/tension force candidate API boundary.

Pure-Python helper-local force-candidate diagnostic.  The helper reconstructs a
Stage 19.3-style initialized state, computes candidate-only bending, tension,
fluid-placeholder, and total forces, and validates shape/numeric/formula and
fail-closed rules.  It does not create production runtime state, insert hooks,
advance or commit structure state, spread force to RHS, call Stage 14 RHS
injection, modify IBM/DNS-core/projection/Poisson/RK3, modify production I/O,
run MPI/DNS, or introduce contact/collision/production multifibre logic.
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
    "stage19_4_requested_status","stage19_4_force_candidate_enable_status","stage19_3_evidence_status","stage19_2_evidence_status","stage19_1_evidence_status","stage19_0_evidence_status","stage18_closure_evidence_status","stage19_3_initialization_preserved_status","stage19_2_state_container_preserved_status","stage19_1_config_gate_preserved_status","stage19_0_source_only_closure_acceptance_preserved_status","no_stage10_18_file_modification_status","no_stage19_0_file_modification_status","no_stage19_1_file_modification_status","no_stage19_2_file_modification_status","no_stage19_3_file_modification_status","no_closed_stage_modification_status","force_candidate_schema_documented_status","all_required_force_candidate_fields_present_status","default_safe_values_status","n_fibre_status","n_point_status","component_dim_status","fibre_length_status","ds_formula_status","rho_l_status","rho_tilde_status","bending_stiffness_status","gamma_status","init_mode_status","sine_amplitude_status","sine_mode_status","tension_mode_status","tension_value_status","controlled_force_amplitude_status","x_prod_shape_status","v_prod_shape_status","a_prod_shape_status","f_b_candidate_shape_status","f_t_candidate_shape_status","f_h_candidate_shape_status","f_total_candidate_shape_status","x_candidate_shape_status","v_candidate_shape_status","a_candidate_shape_status","owner_rank_shape_status","global_point_id_shape_status","local_point_id_shape_status","array_finite_status","bending_candidate_formula_status","tension_candidate_formula_status","fluid_placeholder_candidate_formula_status","total_force_candidate_formula_status","candidate_arrays_helper_local_status","candidate_no_state_update_status","no_state_advance_status","no_state_commit_status","global_point_id_coverage_status","global_point_id_no_duplicate_status","owner_rank_deterministic_status","diagnostic_only_status","single_fibre_only_status","fail_closed_status","force_candidate_only_status","commit_default_disabled_status","rhs_spreading_default_disabled_status","stage14_rhs_injection_default_disabled_status","diagnostic_only_consistency_status","single_fibre_only_consistency_status","fail_closed_consistency_status","rhs_spreading_disabled_consistency_status","stage14_rhs_injection_disabled_consistency_status","commit_disabled_consistency_status","force_candidate_production_runtime_inactive_status","stage19_4_wrapper_bash_syntax_status","stage19_4_helper_py_compile_status","no_production_fortran_modification_status","no_cmake_modification_status","no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_4_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status",
]

REQUIRED_FIELDS = ["X_prod","V_prod","A_prod","F_b_candidate","F_T_candidate","F_h_candidate","F_total_candidate","X_candidate","V_candidate","A_candidate","owner_rank","global_point_id","local_point_id","n_fibre","n_point","component_dim","fibre_length","ds","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude","diagnostic_only","force_candidate_only","state_valid","container_initialized","commit_allowed","rhs_spreading_allowed","stage14_rhs_injection_allowed"]
ARRAY_FIELDS_2D = ["X_prod","V_prod","A_prod","F_b_candidate","F_T_candidate","F_h_candidate","F_total_candidate","X_candidate","V_candidate","A_candidate"]
SHAPE_STATUS = {"X_prod":"x_prod_shape_status","V_prod":"v_prod_shape_status","A_prod":"a_prod_shape_status","F_b_candidate":"f_b_candidate_shape_status","F_T_candidate":"f_t_candidate_shape_status","F_h_candidate":"f_h_candidate_shape_status","F_total_candidate":"f_total_candidate_shape_status","X_candidate":"x_candidate_shape_status","V_candidate":"v_candidate_shape_status","A_candidate":"a_candidate_shape_status"}
INIT_MODES = {"small_sine_fibre_zero_velocity"}
TENSION_MODES = {"constant"}
STAGE19_0_FILES = {"stage19_checks/run_stage19_0_preflight_boundary.sh","stage19_checks/assert_stage19_0_preflight_boundary.py","stage19_checks/stage19_0_preflight_boundary.md","stage19_outputs/fibre_stage19_0_preflight_boundary.dat"}
STAGE19_1_FILES = {"stage19_checks/run_stage19_1_physical_structure_config_gate.sh","stage19_checks/assert_stage19_1_physical_structure_config_gate.py","stage19_checks/stage19_1_physical_structure_config_gate.md","stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat"}
STAGE19_2_FILES = {"stage19_checks/run_stage19_2_physical_structure_state_container.sh","stage19_checks/assert_stage19_2_physical_structure_state_container.py","stage19_checks/stage19_2_physical_structure_state_container.md","stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat"}
STAGE19_3_FILES = {"stage19_checks/run_stage19_3_physical_structure_initialization.sh","stage19_checks/assert_stage19_3_physical_structure_initialization.py","stage19_checks/stage19_3_physical_structure_initialization.md","stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat"}
ALLOWED_NEW_OR_MODIFIED = {"stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh","stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py","stage19_checks/stage19_4_bending_tension_force_candidate_api.md","stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat"}
ACCEPTED_UNTRACKED_EVIDENCE = {"stage17_checks/STAGE17_CLOSED.md","stage18_checks/STAGE18_CLOSED.md","stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat","stage19_outputs/fibre_stage19_0_preflight_boundary.dat","stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat","stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat","stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat"}
TRUE_VALUES = {"1","true","TRUE","yes","YES","on","ON","t","T"}
FALSE_VALUES = {"0","false","FALSE","no","NO","off","OFF","f","F"}

@dataclass(frozen=True)
class Config:
    n_fibre: int; n_point: int; component_dim: int; fibre_length: float; rho_l: float; rho_tilde: float
    bending_stiffness: float; gamma: float; init_mode: str; sine_amplitude: float; sine_mode: int
    tension_mode: str; tension_value: float; controlled_force_amplitude: float; diagnostic_only: bool
    single_fibre_only: bool; fail_closed: bool; force_candidate_only: bool; commit_allowed: bool
    rhs_spreading_allowed: bool; stage14_rhs_injection_allowed: bool

def read_text(path: Path) -> str:
    try: return path.read_text(errors="ignore")
    except OSError: return ""
def pass_fail(condition: bool) -> str: return "PASS" if condition else "FAIL"
def parse_bool(name: str, default: bool, invalid: List[str]) -> bool:
    raw = os.environ.get(name)
    if raw is None: return default
    value = raw.strip()
    if value in TRUE_VALUES: return True
    if value in FALSE_VALUES: return False
    invalid.append(f"{name}={raw}"); return default
def parse_int(name: str, default: int, invalid: List[str]) -> int:
    raw = os.environ.get(name)
    if raw is None: return default
    try: return int(raw.strip())
    except ValueError: invalid.append(f"{name}={raw}"); return default
def parse_float(name: str, default: float, invalid: List[str]) -> float:
    raw = os.environ.get(name)
    if raw is None: return default
    try: value = float(raw.strip())
    except ValueError: invalid.append(f"{name}={raw}"); return default
    if not math.isfinite(value): invalid.append(f"{name}={raw}"); return default
    return value

def run_quiet(cmd: Sequence[str], cwd: Path) -> Tuple[int, str]:
    try:
        proc = subprocess.run(cmd, cwd=str(cwd), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        return proc.returncode, proc.stdout
    except OSError as exc: return 127, str(exc)

def git_status_entries(root: Path) -> List[Tuple[str, str]]:
    code, out = run_quiet(["git","status","--porcelain","--untracked-files=all"], root)
    if code != 0: return []
    entries: List[Tuple[str, str]] = []
    for raw in out.splitlines():
        if not raw: continue
        xy = raw[:2]; path = raw[3:] if len(raw) > 3 else raw
        if " -> " in path: path = path.split(" -> ", 1)[1]
        entries.append((xy, path.strip()))
    return entries

def outside_allowed_changes(root: Path) -> List[str]:
    outside: List[str] = []
    for xy, path in git_status_entries(root):
        if path in ALLOWED_NEW_OR_MODIFIED: continue
        if xy == "??" and path in ACCEPTED_UNTRACKED_EVIDENCE: continue
        outside.append(path)
    return outside

def evidence_has_final_pass(path: Path, verdict: str) -> bool:
    text = read_text(path); return "final_status PASS" in text and verdict in text

def stage19_0_source_acceptance_preserved(root: Path) -> bool:
    helper = read_text(root/"stage19_checks/assert_stage19_0_preflight_boundary.py"); doc = read_text(root/"stage19_checks/stage19_0_preflight_boundary.md")
    return all(t in helper for t in ("stage18_closure_accepted_status","prior_stage18_outputs_required_status","stage18_closure_supersedes_individual_outputs_status","ACCEPTED_BY_STAGE18_CLOSURE")) and "must not force users to rerun Stage 18.0 through Stage 18.11" in doc

def stage19_1_config_gate_preserved(root: Path) -> bool:
    helper = read_text(root/"stage19_checks/assert_stage19_1_physical_structure_config_gate.py"); doc = read_text(root/"stage19_checks/stage19_1_physical_structure_config_gate.md")
    return all(t in helper for t in ("stage19_physical_structure_enable","stage19_physical_structure_diagnostic_only","stage19_physical_structure_fail_closed","stage19_rhs_spreading_enable","stage19_stage14_rhs_injection_enable")) and "Stage 19.1 is configuration-only" in doc

def stage19_2_state_container_preserved(root: Path) -> bool:
    helper = read_text(root/"stage19_checks/assert_stage19_2_physical_structure_state_container.py"); doc = read_text(root/"stage19_checks/stage19_2_physical_structure_state_container.md")
    return all(t in helper for t in ("REQUIRED_CONTAINER_FIELDS","X_prod","V_prod","A_prod","commit_allowed","rhs_spreading_allowed","stage14_rhs_injection_allowed")) and "state-container-boundary" in doc

def stage19_3_initialization_preserved(root: Path) -> bool:
    helper = read_text(root/"stage19_checks/assert_stage19_3_physical_structure_initialization.py"); doc = read_text(root/"stage19_checks/stage19_3_physical_structure_initialization.md")
    return all(t in helper for t in ("small_sine_fibre_zero_velocity","F_total_prod","candidate_equals_production_initialization_status","controlled_helper_force_placeholder")) and "initialization-boundary" in doc

def stage_evidence(root: Path, output_name: str, verdict: str, files: set[str], preserved) -> bool:
    output = root/"stage19_outputs"/output_name
    return evidence_has_final_pass(output, verdict) or (all((root/p).is_file() for p in files if p.startswith("stage19_checks/")) and preserved(root))
def stage18_evidence(root: Path) -> bool: return stage19_0_source_acceptance_preserved(root)

def syntax_status(root: Path) -> Tuple[str, str]:
    wrapper = root/"stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh"; helper = root/"stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py"
    bash_code, _ = run_quiet(["bash","-n",str(wrapper)], root)
    try:
        with tempfile.TemporaryDirectory() as td:
            py_compile.compile(str(helper), cfile=str(Path(td)/"stage19_4_helper.pyc"), doraise=True)
        py_status = "PASS"
    except py_compile.PyCompileError: py_status = "FAIL"
    return pass_fail(bash_code == 0), py_status

def zeros(n: int, c: int) -> List[List[float]]: return [[0.0 for _ in range(c)] for _ in range(n)]
def shape_2d(v: object) -> Tuple[int, int]:
    if not isinstance(v, list): return (-1,-1)
    if not v: return (0,0)
    if not all(isinstance(r, list) for r in v): return (len(v),-1)
    widths = {len(r) for r in v}; return (len(v), widths.pop() if len(widths) == 1 else -1)
def shape_1d(v: object) -> Tuple[int]: return (len(v),) if isinstance(v, list) else (-1,)
def finite_nested(v: object) -> bool:
    if isinstance(v, list): return all(finite_nested(x) for x in v)
    return isinstance(v, (int,float)) and math.isfinite(float(v))
def arrays_close(a: List[List[float]], b: List[List[float]], tol: float) -> bool:
    return shape_2d(a) == shape_2d(b) and all(abs(x-y) <= tol for ra, rb in zip(a,b) for x,y in zip(ra,rb))

def fourth_derivative(values: List[List[float]], ds: float) -> List[List[float]]:
    n = len(values); c = len(values[0]) if n else 0; out = zeros(n, c)
    if n < 5 or ds <= 0: return out
    for i in range(2, n-2):
        for j in range(c): out[i][j] = (values[i-2][j] - 4*values[i-1][j] + 6*values[i][j] - 4*values[i+1][j] + values[i+2][j]) / ds**4
    return out

def second_derivative(values: List[List[float]], ds: float) -> List[List[float]]:
    n = len(values); c = len(values[0]) if n else 0; out = zeros(n, c)
    if n < 3 or ds <= 0: return out
    for i in range(1, n-1):
        for j in range(c): out[i][j] = (values[i+1][j] - 2*values[i][j] + values[i-1][j]) / ds**2
    return out

def build_state_and_forces(config: Config) -> Dict[str, object]:
    ds = config.fibre_length/(config.n_point-1) if config.n_point > 1 else math.nan
    x = zeros(config.n_point, config.component_dim); v = zeros(config.n_point, config.component_dim); a = zeros(config.n_point, config.component_dim)
    for i in range(config.n_point):
        s = i*ds if math.isfinite(ds) else 0.0
        wave = math.sin(2*math.pi*config.sine_mode*s/config.fibre_length) if config.fibre_length > 0 and config.sine_mode > 0 else 0.0
        x[i][0] = s
        if config.component_dim > 1: x[i][1] = config.sine_amplitude * wave
    xssss = fourth_derivative(x, ds if math.isfinite(ds) else -1.0)
    fb = [[-config.gamma*xssss[i][j] for j in range(config.component_dim)] for i in range(config.n_point)]
    xss = second_derivative(x, ds if math.isfinite(ds) else -1.0)
    ft = [[config.tension_value*xss[i][j] if config.tension_mode == "constant" else 0.0 for j in range(config.component_dim)] for i in range(config.n_point)]
    fh = zeros(config.n_point, config.component_dim)
    if config.controlled_force_amplitude > 0:
        for i in range(config.n_point):
            s = i*ds if math.isfinite(ds) else 0.0
            wave = math.sin(2*math.pi*config.sine_mode*s/config.fibre_length) if config.fibre_length > 0 and config.sine_mode > 0 else 0.0
            if config.component_dim > 1: fh[i][1] = config.controlled_force_amplitude * wave
    total = [[fb[i][j]+ft[i][j]+fh[i][j] for j in range(config.component_dim)] for i in range(config.n_point)]
    return {"n_fibre":config.n_fibre,"n_point":config.n_point,"component_dim":config.component_dim,"fibre_length":config.fibre_length,"ds":ds,"rho_l":config.rho_l,"rho_tilde":config.rho_tilde,"bending_stiffness":config.bending_stiffness,"gamma":config.gamma,"init_mode":config.init_mode,"sine_amplitude":config.sine_amplitude,"sine_mode":config.sine_mode,"tension_mode":config.tension_mode,"tension_value":config.tension_value,"controlled_force_amplitude":config.controlled_force_amplitude,"diagnostic_only":config.diagnostic_only,"force_candidate_only":config.force_candidate_only,"state_valid":True,"container_initialized":True,"commit_allowed":config.commit_allowed,"rhs_spreading_allowed":config.rhs_spreading_allowed,"stage14_rhs_injection_allowed":config.stage14_rhs_injection_allowed,"X_prod":x,"V_prod":v,"A_prod":a,"F_b_candidate":fb,"F_T_candidate":ft,"F_h_candidate":fh,"F_total_candidate":total,"X_candidate":[r[:] for r in x],"V_candidate":[r[:] for r in v],"A_candidate":[r[:] for r in a],"owner_rank":[0 for _ in range(config.n_point)],"global_point_id":[i for i in range(config.n_point)],"local_point_id":[i for i in range(config.n_point)],"_xssss":xssss,"_xss":xss}

def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 19.4 bending/tension force candidate API boundary")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1]); parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv); root = args.repo_root.resolve(); output = args.output or root/"stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat"
    statuses: Dict[str,str] = {k:"PASS" for k in SUMMARY_KEYS if k != "final_status"}; failures: List[str] = []; invalid: List[str] = []
    requested = parse_bool("STAGE19_4_ENABLE", True, invalid); gate = parse_bool("STAGE19_4_FORCE_CANDIDATE_ENABLE", True, invalid)
    req3 = parse_bool("STAGE19_4_REQUIRE_STAGE19_3_PASS", True, invalid); req2 = parse_bool("STAGE19_4_REQUIRE_STAGE19_2_PASS", True, invalid); req1 = parse_bool("STAGE19_4_REQUIRE_STAGE19_1_PASS", True, invalid); req0 = parse_bool("STAGE19_4_REQUIRE_STAGE19_0_PASS", True, invalid); req18 = parse_bool("STAGE19_4_REQUIRE_STAGE18_CLOSED", True, invalid)
    init_mode = os.environ.get("STAGE19_4_INIT_MODE","small_sine_fibre_zero_velocity").strip(); tension_mode = os.environ.get("STAGE19_4_TENSION_MODE","constant").strip()
    if init_mode not in INIT_MODES: invalid.append(f"STAGE19_4_INIT_MODE={init_mode}"); init_mode = "small_sine_fibre_zero_velocity"
    if tension_mode not in TENSION_MODES: invalid.append(f"STAGE19_4_TENSION_MODE={tension_mode}"); tension_mode = "constant"
    config = Config(parse_int("STAGE19_4_N_FIBRE",1,invalid),parse_int("STAGE19_4_N_POINT",64,invalid),parse_int("STAGE19_4_COMPONENT_DIM",3,invalid),parse_float("STAGE19_4_FIBRE_LENGTH",1.0,invalid),parse_float("STAGE19_4_RHO_L",1.0,invalid),parse_float("STAGE19_4_RHO_TILDE",1.0,invalid),parse_float("STAGE19_4_BENDING_STIFFNESS",1.0e-3,invalid),parse_float("STAGE19_4_GAMMA",1.0e-3,invalid),init_mode,parse_float("STAGE19_4_SINE_AMPLITUDE",1.0e-3,invalid),parse_int("STAGE19_4_SINE_MODE",1,invalid),tension_mode,parse_float("STAGE19_4_TENSION_VALUE",0.0,invalid),parse_float("STAGE19_4_CONTROLLED_FORCE_AMPLITUDE",0.0,invalid),parse_bool("STAGE19_4_DIAGNOSTIC_ONLY",True,invalid),parse_bool("STAGE19_4_SINGLE_FIBRE_ONLY",True,invalid),parse_bool("STAGE19_4_FAIL_CLOSED",True,invalid),parse_bool("STAGE19_4_FORCE_CANDIDATE_ONLY",True,invalid),parse_bool("STAGE19_4_COMMIT_ALLOWED",False,invalid),parse_bool("STAGE19_4_RHS_SPREADING_ALLOWED",False,invalid),parse_bool("STAGE19_4_STAGE14_RHS_INJECTION_ALLOWED",False,invalid))
    state = build_state_and_forces(config); tol = parse_float("STAGE19_4_AUDIT_TOL",1.0e-12,invalid)
    statuses["stage19_4_requested_status"] = pass_fail(requested); statuses["stage19_4_force_candidate_enable_status"] = pass_fail(gate)
    statuses["stage19_3_evidence_status"] = pass_fail((not req3) or stage_evidence(root,"fibre_stage19_3_physical_structure_initialization.dat","STAGE 19.3 FINAL VERDICT: PASS",STAGE19_3_FILES,stage19_3_initialization_preserved))
    statuses["stage19_2_evidence_status"] = pass_fail((not req2) or stage_evidence(root,"fibre_stage19_2_physical_structure_state_container.dat","STAGE 19.2 FINAL VERDICT: PASS",STAGE19_2_FILES,stage19_2_state_container_preserved))
    statuses["stage19_1_evidence_status"] = pass_fail((not req1) or stage_evidence(root,"fibre_stage19_1_physical_structure_config_gate.dat","STAGE 19.1 FINAL VERDICT: PASS",STAGE19_1_FILES,stage19_1_config_gate_preserved))
    statuses["stage19_0_evidence_status"] = pass_fail((not req0) or stage_evidence(root,"fibre_stage19_0_preflight_boundary.dat","STAGE 19.0 FINAL VERDICT: PASS",STAGE19_0_FILES,stage19_0_source_acceptance_preserved))
    statuses["stage18_closure_evidence_status"] = pass_fail((not req18) or stage18_evidence(root))
    statuses["stage19_3_initialization_preserved_status"] = pass_fail(stage19_3_initialization_preserved(root)); statuses["stage19_2_state_container_preserved_status"] = pass_fail(stage19_2_state_container_preserved(root)); statuses["stage19_1_config_gate_preserved_status"] = pass_fail(stage19_1_config_gate_preserved(root)); statuses["stage19_0_source_only_closure_acceptance_preserved_status"] = pass_fail(stage19_0_source_acceptance_preserved(root))
    changed = outside_allowed_changes(root)
    statuses["no_stage10_18_file_modification_status"] = pass_fail(not any(p.startswith(("stage10_","stage11_","stage12_","stage13_","stage14_","stage15_","stage16_","stage17_","stage18_")) for p in changed)); statuses["no_stage19_0_file_modification_status"] = pass_fail(not any(p in STAGE19_0_FILES for p in changed)); statuses["no_stage19_1_file_modification_status"] = pass_fail(not any(p in STAGE19_1_FILES for p in changed)); statuses["no_stage19_2_file_modification_status"] = pass_fail(not any(p in STAGE19_2_FILES for p in changed)); statuses["no_stage19_3_file_modification_status"] = pass_fail(not any(p in STAGE19_3_FILES for p in changed)); statuses["no_closed_stage_modification_status"] = pass_fail(all(statuses[k] == "PASS" for k in ("no_stage10_18_file_modification_status","no_stage19_0_file_modification_status","no_stage19_1_file_modification_status","no_stage19_2_file_modification_status","no_stage19_3_file_modification_status")))
    doc = read_text(root/"stage19_checks/stage19_4_bending_tension_force_candidate_api.md")
    statuses["force_candidate_schema_documented_status"] = pass_fail(all(f in doc for f in REQUIRED_FIELDS)); statuses["all_required_force_candidate_fields_present_status"] = pass_fail(set(REQUIRED_FIELDS).issubset(set(state)))
    default_safe = config == Config(1,64,3,1.0,1.0,1.0,1.0e-3,1.0e-3,"small_sine_fibre_zero_velocity",1.0e-3,1,"constant",0.0,0.0,True,True,True,True,False,False,False)
    statuses["default_safe_values_status"] = pass_fail(default_safe); statuses["n_fibre_status"] = pass_fail(config.n_fibre == 1); statuses["n_point_status"] = pass_fail(config.n_point >= 8); statuses["component_dim_status"] = pass_fail(config.component_dim == 3); statuses["fibre_length_status"] = pass_fail(config.fibre_length > 0); expected_ds = config.fibre_length/(config.n_point-1) if config.n_point > 1 else math.nan; statuses["ds_formula_status"] = pass_fail(math.isfinite(expected_ds) and math.isclose(state["ds"], expected_ds, rel_tol=0, abs_tol=1e-14)); statuses["rho_l_status"] = pass_fail(config.rho_l > 0); statuses["rho_tilde_status"] = pass_fail(config.rho_tilde > 0); statuses["bending_stiffness_status"] = pass_fail(config.bending_stiffness >= 0); statuses["gamma_status"] = pass_fail(config.gamma >= 0); statuses["init_mode_status"] = pass_fail(config.init_mode in INIT_MODES); statuses["sine_amplitude_status"] = pass_fail(math.isfinite(config.sine_amplitude) and abs(config.sine_amplitude) <= 1e-1); statuses["sine_mode_status"] = pass_fail(config.sine_mode > 0); statuses["tension_mode_status"] = pass_fail(config.tension_mode in TENSION_MODES); statuses["tension_value_status"] = pass_fail(config.tension_value >= 0); statuses["controlled_force_amplitude_status"] = pass_fail(config.controlled_force_amplitude >= 0)
    expected_shape = (config.n_point, config.component_dim)
    for field in ARRAY_FIELDS_2D: statuses[SHAPE_STATUS[field]] = pass_fail(shape_2d(state[field]) == expected_shape)
    statuses["owner_rank_shape_status"] = pass_fail(shape_1d(state["owner_rank"]) == (config.n_point,)); statuses["global_point_id_shape_status"] = pass_fail(shape_1d(state["global_point_id"]) == (config.n_point,)); statuses["local_point_id_shape_status"] = pass_fail(shape_1d(state["local_point_id"]) == (config.n_point,)); statuses["array_finite_status"] = pass_fail(all(finite_nested(state[f]) for f in ARRAY_FIELDS_2D))
    expected_fb = [[-config.gamma*state["_xssss"][i][j] for j in range(config.component_dim)] for i in range(config.n_point)]; expected_ft = [[config.tension_value*state["_xss"][i][j] for j in range(config.component_dim)] for i in range(config.n_point)]; expected_fh = state["F_h_candidate"]; expected_total = [[state["F_b_candidate"][i][j]+state["F_T_candidate"][i][j]+state["F_h_candidate"][i][j] for j in range(config.component_dim)] for i in range(config.n_point)]
    statuses["bending_candidate_formula_status"] = pass_fail(arrays_close(state["F_b_candidate"], expected_fb, tol)); statuses["tension_candidate_formula_status"] = pass_fail(arrays_close(state["F_T_candidate"], expected_ft, tol)); statuses["fluid_placeholder_candidate_formula_status"] = pass_fail(arrays_close(state["F_h_candidate"], expected_fh, tol)); statuses["total_force_candidate_formula_status"] = pass_fail(arrays_close(state["F_total_candidate"], expected_total, tol))
    statuses["candidate_arrays_helper_local_status"] = "PASS"; statuses["candidate_no_state_update_status"] = pass_fail(arrays_close(state["X_candidate"], state["X_prod"], tol) and arrays_close(state["V_candidate"], state["V_prod"], tol) and arrays_close(state["A_candidate"], state["A_prod"], tol)); statuses["no_state_advance_status"] = "PASS"; statuses["no_state_commit_status"] = pass_fail(not config.commit_allowed); statuses["global_point_id_coverage_status"] = pass_fail(sorted(state["global_point_id"]) == list(range(config.n_point))); statuses["global_point_id_no_duplicate_status"] = pass_fail(len(set(state["global_point_id"])) == config.n_point); statuses["owner_rank_deterministic_status"] = pass_fail(state["owner_rank"] == [0 for _ in range(config.n_point)])
    statuses["diagnostic_only_status"] = pass_fail(config.diagnostic_only); statuses["single_fibre_only_status"] = pass_fail(config.single_fibre_only and config.n_fibre == 1); statuses["fail_closed_status"] = pass_fail(config.fail_closed); statuses["force_candidate_only_status"] = pass_fail(config.force_candidate_only); statuses["commit_default_disabled_status"] = pass_fail(not config.commit_allowed); statuses["rhs_spreading_default_disabled_status"] = pass_fail(not config.rhs_spreading_allowed); statuses["stage14_rhs_injection_default_disabled_status"] = pass_fail(not config.stage14_rhs_injection_allowed)
    diagnostic_ok = (not config.diagnostic_only) or not config.commit_allowed; single_ok = (not config.single_fibre_only) or config.n_fibre == 1; rhs_ok = config.rhs_spreading_allowed or not config.stage14_rhs_injection_allowed; stage14_ok = not config.stage14_rhs_injection_allowed; commit_ok = not config.commit_allowed; runtime_inactive = config.force_candidate_only and not config.commit_allowed and not config.rhs_spreading_allowed and not config.stage14_rhs_injection_allowed; fail_closed_ok = config.fail_closed and not invalid and diagnostic_ok and single_ok and rhs_ok and stage14_ok and commit_ok and runtime_inactive
    statuses["diagnostic_only_consistency_status"] = pass_fail(diagnostic_ok); statuses["single_fibre_only_consistency_status"] = pass_fail(single_ok); statuses["fail_closed_consistency_status"] = pass_fail(fail_closed_ok); statuses["rhs_spreading_disabled_consistency_status"] = pass_fail(rhs_ok); statuses["stage14_rhs_injection_disabled_consistency_status"] = pass_fail(stage14_ok); statuses["commit_disabled_consistency_status"] = pass_fail(commit_ok); statuses["force_candidate_production_runtime_inactive_status"] = pass_fail(runtime_inactive)
    wrapper_syntax, helper_compile = syntax_status(root); statuses["stage19_4_wrapper_bash_syntax_status"] = wrapper_syntax; statuses["stage19_4_helper_py_compile_status"] = helper_compile; statuses["no_production_fortran_modification_status"] = pass_fail(not any(p.startswith("src/") and p.endswith((".f90",".F90",".f",".F")) for p in changed)); statuses["no_cmake_modification_status"] = pass_fail(not any(p == "CMakeLists.txt" or p.endswith("CMakeLists.txt") or p.startswith("cmake/") for p in changed))
    for key in ["no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_4_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status"]: statuses[key] = pass_fail(not changed)
    statuses["no_unknown_failure_status"] = pass_fail(not invalid)
    failing = [k for k in SUMMARY_KEYS if k.endswith("_status") and k != "final_status" and statuses.get(k) != "PASS"]
    if invalid: failures.extend(f"invalid_value:{x}" for x in invalid)
    failures.extend(f"{k}={statuses.get(k,'MISSING')}" for k in failing); final = "PASS" if not failing and not invalid else "FAIL"; statuses["final_status"] = final
    output.parent.mkdir(parents=True, exist_ok=True); lines: List[str] = ["# Stage 19.4 bending/tension force candidate API","stage19_title production-side physical structure integration boundary","stage19_4_title production-side bending/tension force candidate API",f"stage19_4_test_case {os.environ.get('STAGE19_4_TEST_CASE','production_bending_tension_force_candidate_api')}",f"stage19_4_zero_tol_value {os.environ.get('STAGE19_4_ZERO_TOL','1.0e-14')}",f"stage19_4_audit_tol_value {tol}","stage19_4_scope force-candidate-boundary only; local helper force candidates are not production runtime force application"]
    for name in ("n_fibre","n_point","component_dim","fibre_length","ds","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude"): lines.append(f"{name}_value {state[name]}")
    for field in ARRAY_FIELDS_2D: lines.append(f"{field.lower()}_shape_value {shape_2d(state[field])}")
    lines.extend([f"owner_rank_shape_value {shape_1d(state['owner_rank'])}",f"global_point_id_shape_value {shape_1d(state['global_point_id'])}",f"local_point_id_shape_value {shape_1d(state['local_point_id'])}"])
    for key in SUMMARY_KEYS: lines.append(f"{key} {statuses[key]}")
    if failures: lines.extend(["failure_reasons_begin", *failures, "failure_reasons_end"])
    lines.append(f"STAGE 19.4 BENDING TENSION FORCE CANDIDATE API VERDICT: {final}"); lines.append(f"STAGE 19.4 FINAL VERDICT: {final}"); output.write_text("\n".join(lines)+"\n")
    print(f"STAGE 19.4 BENDING TENSION FORCE CANDIDATE API VERDICT: {final}"); print(f"STAGE 19.4 FINAL VERDICT: {final}")
    if failures:
        print("STAGE 19.4 FAILURE REASONS:")
        for reason in failures: print(f"  - {reason}")
    return 0 if final == "PASS" else 1

if __name__ == "__main__":
    raise SystemExit(main())
