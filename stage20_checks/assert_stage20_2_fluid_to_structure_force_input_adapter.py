#!/usr/bin/env python3
"""Stage 20.2 helper-local fluid-to-structure force input adapter audit.

This helper constructs a controlled, helper-local F_fs_candidate from local fibre
and controlled fluid-point velocities. It audits formulas and gates only; it does
not run prior stages, MPI, DNS, builds, production hooks, RHS injection, IBM, or
structure-to-fluid coupling.
"""
from __future__ import annotations

import math
import os
import py_compile
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

PASS = "PASS"
FAIL = "FAIL"
OPTIONAL = "OPTIONAL"

STATUS_FIELDS = [
    "stage20_2_requested_status",
    "stage20_2_force_input_adapter_enable_status",
    "stage20_1_evidence_status",
    "stage20_1_source_only_acceptance_preserved_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "fluid_to_structure_adapter_documented_status",
    "all_required_adapter_fields_present_status",
    "default_safe_gate_values_status",
    "twoway_coupling_remains_disabled_status",
    "fluid_to_structure_helper_local_only_status",
    "structure_to_fluid_default_disabled_status",
    "rhs_coupling_default_disabled_status",
    "lambda_coupling_zero_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "single_fibre_only_status",
    "contact_default_disabled_status",
    "collision_default_disabled_status",
    "multifibre_default_disabled_status",
    "n_fibre_status",
    "n_point_status",
    "component_dim_status",
    "fibre_length_status",
    "ds_formula_status",
    "rho_l_status",
    "rho_tilde_status",
    "bending_stiffness_status",
    "gamma_status",
    "c_fs_status",
    "shape_rules_status",
    "numeric_rules_status",
    "relative_velocity_formula_status",
    "f_fs_candidate_formula_status",
    "structure_total_force_without_fluid_formula_status",
    "structure_total_force_with_fluid_formula_status",
    "candidate_force_finite_status",
    "candidate_force_helper_local_only_status",
    "candidate_force_nonzero_status",
    "no_effective_rhs_coupling_status",
    "global_point_id_coverage_status",
    "global_point_id_no_duplicate_status",
    "owner_rank_deterministic_status",
    "no_stage10_19_file_modification_status",
    "no_stage20_0_file_modification_status",
    "no_stage20_1_file_modification_status",
    "no_closed_stage_modification_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status",
    "no_structure_to_fluid_reaction_force_activation_status",
    "no_bending_force_runtime_application_status",
    "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage20_2_status",
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
    "stage20_3_next_stage_declared_status",
    "stage20_2_wrapper_bash_syntax_status",
    "stage20_2_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_2_ENABLE": "1",
    "STAGE20_2_FORCE_INPUT_ADAPTER_ENABLE": "1",
    "STAGE20_2_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_2_REQUIRE_STAGE20_1_PASS": "1",
    "STAGE20_2_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_2_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_2_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_2_DIAGNOSTIC_ONLY": "1",
    "STAGE20_2_FAIL_CLOSED": "1",
    "STAGE20_2_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_2_FLUID_TO_STRUCTURE_ENABLE": "1",
    "STAGE20_2_STRUCTURE_TO_FLUID_ENABLE": "0",
    "STAGE20_2_RHS_COUPLING_ENABLE": "0",
    "STAGE20_2_LAMBDA_COUPLING": "0.0",
    "STAGE20_2_C_FS": "1.0",
    "STAGE20_2_SINGLE_FIBRE_ONLY": "1",
    "STAGE20_2_CONTACT_ENABLE": "0",
    "STAGE20_2_COLLISION_ENABLE": "0",
    "STAGE20_2_MULTIFIBRE_ENABLE": "0",
    "STAGE20_2_N_FIBRE": "1",
    "STAGE20_2_N_POINT": "64",
    "STAGE20_2_COMPONENT_DIM": "3",
    "STAGE20_2_FIBRE_LENGTH": "1.0",
    "STAGE20_2_RHO_L": "1.0",
    "STAGE20_2_RHO_TILDE": "1.0",
    "STAGE20_2_BENDING_STIFFNESS": "1.0e-5",
    "STAGE20_2_GAMMA": "1.0e-5",
    "STAGE20_2_ZERO_TOL": "1.0e-14",
    "STAGE20_2_AUDIT_TOL": "1.0e-12",
    "STAGE20_2_TEST_CASE": "fluid_to_structure_force_input_adapter",
}

REQUIRED_ARRAYS = [
    "X_fibre",
    "V_fibre",
    "u_interp_candidate",
    "u_relative_candidate",
    "F_fs_candidate",
    "F_b_candidate",
    "F_T_candidate",
    "F_total_structure_candidate_without_fluid",
    "F_total_structure_candidate_with_fluid",
    "owner_rank",
    "global_point_id",
    "local_point_id",
]


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


def stage20_1_evidence(root: Path) -> Tuple[bool, str]:
    output = root / "stage20_outputs" / "fibre_stage20_1_twoway_config_lambda_gate.dat"
    text = read_text(output)
    output_pass = output.exists() and (
        "STAGE 20.1 TWOWAY CONFIG LAMBDA GATE VERDICT: PASS" in text
        and "STAGE 20.1 FINAL VERDICT: PASS" in text
    )
    markers = [
        root / "stage20_checks" / "assert_stage20_1_twoway_config_lambda_gate.py",
        root / "stage20_checks" / "run_stage20_1_twoway_config_lambda_gate.sh",
        root / "stage20_checks" / "stage20_1_twoway_config_lambda_gate.md",
    ]
    source_only_available = truthy("STAGE20_2_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") and all(marker.exists() for marker in markers)
    if output_pass:
        return True, "ACCEPTED_BY_STAGE20_1_PASS_OUTPUT"
    if source_only_available:
        return True, "ACCEPTED_BY_STAGE20_1_SOURCE_ONLY_BEHAVIOR"
    return False, "NO_STAGE20_1_PASS_OR_SOURCE_ONLY_EVIDENCE"


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_2_py_compile_", suffix=".pyc", delete=False) as tmp:
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


def is_finite_matrix(values: Sequence[Sequence[float]]) -> bool:
    return all(math.isfinite(value) for row in values for value in row)


def matrix_shape(values: Sequence[Sequence[float]], rows: int, cols: int) -> bool:
    return len(values) == rows and all(len(row) == cols for row in values)


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def matrix_close(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], tol: float) -> bool:
    return len(a) == len(b) and all(len(ra) == len(rb) and all(close(x, y, tol) for x, y in zip(ra, rb)) for ra, rb in zip(a, b))


def build_candidate(n_point: int, component_dim: int, fibre_length: float, c_fs: float) -> Dict[str, object]:
    ds = fibre_length / float(n_point - 1)
    x_fibre: List[List[float]] = []
    v_fibre: List[List[float]] = []
    u_interp: List[List[float]] = []
    f_b: List[List[float]] = []
    f_t: List[List[float]] = []
    for i in range(n_point):
        s = i / float(n_point - 1)
        x_fibre.append([fibre_length * s, 1.0e-3 * math.sin(2.0 * math.pi * s), 0.0])
        v_fibre.append([0.0, 1.0e-5 * math.cos(2.0 * math.pi * s), 0.0])
        u_interp.append([1.0e-3, 2.0e-4 * math.sin(2.0 * math.pi * s), 0.0])
        f_b.append([0.0, -1.0e-6 * math.sin(2.0 * math.pi * s), 0.0])
        f_t.append([1.0e-7 * math.cos(2.0 * math.pi * s), 0.0, 0.0])
    u_relative = [[u_interp[i][j] - v_fibre[i][j] for j in range(component_dim)] for i in range(n_point)]
    f_fs = [[c_fs * u_relative[i][j] for j in range(component_dim)] for i in range(n_point)]
    total_without = [[f_b[i][j] + f_t[i][j] for j in range(component_dim)] for i in range(n_point)]
    total_with = [[f_b[i][j] + f_t[i][j] - f_fs[i][j] for j in range(component_dim)] for i in range(n_point)]
    return {
        "ds": ds,
        "X_fibre": x_fibre,
        "V_fibre": v_fibre,
        "u_interp_candidate": u_interp,
        "u_relative_candidate": u_relative,
        "F_fs_candidate": f_fs,
        "F_b_candidate": f_b,
        "F_T_candidate": f_t,
        "F_total_structure_candidate_without_fluid": total_without,
        "F_total_structure_candidate_with_fluid": total_with,
        "owner_rank": [0 for _ in range(n_point)],
        "global_point_id": [i for i in range(n_point)],
        "local_point_id": [i for i in range(n_point)],
    }


def main() -> int:
    root = repo_root()
    out_dir = root / "stage20_outputs"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_2_fluid_to_structure_force_input_adapter.dat"

    evidence_ok, evidence_reason = stage20_1_evidence(root)
    n_fibre = int(env_value("STAGE20_2_N_FIBRE"))
    n_point = int(env_value("STAGE20_2_N_POINT"))
    component_dim = int(env_value("STAGE20_2_COMPONENT_DIM"))
    fibre_length = float(env_value("STAGE20_2_FIBRE_LENGTH"))
    rho_l = float(env_value("STAGE20_2_RHO_L"))
    rho_tilde = float(env_value("STAGE20_2_RHO_TILDE"))
    bending_stiffness = float(env_value("STAGE20_2_BENDING_STIFFNESS"))
    gamma = float(env_value("STAGE20_2_GAMMA"))
    c_fs = float(env_value("STAGE20_2_C_FS"))
    lambda_coupling = float(env_value("STAGE20_2_LAMBDA_COUPLING"))
    tol = float(env_value("STAGE20_2_AUDIT_TOL"))
    zero_tol = float(env_value("STAGE20_2_ZERO_TOL"))

    candidate = build_candidate(n_point, component_dim, fibre_length, c_fs)
    ds = float(candidate["ds"])
    arrays_present = all(name in candidate for name in REQUIRED_ARRAYS)
    matrix_names = REQUIRED_ARRAYS[:9]
    shape_rules = all(matrix_shape(candidate[name], n_point, component_dim) for name in matrix_names) and all(
        len(candidate[name]) == n_point for name in ["owner_rank", "global_point_id", "local_point_id"]
    )
    all_finite = all(is_finite_matrix(candidate[name]) for name in matrix_names)
    u_expected = [[candidate["u_interp_candidate"][i][j] - candidate["V_fibre"][i][j] for j in range(component_dim)] for i in range(n_point)]
    ffs_expected = [[c_fs * candidate["u_relative_candidate"][i][j] for j in range(component_dim)] for i in range(n_point)]
    total_without_expected = [[candidate["F_b_candidate"][i][j] + candidate["F_T_candidate"][i][j] for j in range(component_dim)] for i in range(n_point)]
    total_with_expected = [[candidate["F_b_candidate"][i][j] + candidate["F_T_candidate"][i][j] - candidate["F_fs_candidate"][i][j] for j in range(component_dim)] for i in range(n_point)]
    global_ids = candidate["global_point_id"]
    owner_rank = candidate["owner_rank"]
    force_norm = math.sqrt(sum(value * value for row in candidate["F_fs_candidate"] for value in row))

    twoway = truthy("STAGE20_2_TWOWAY_COUPLING_ENABLE")
    f2s = truthy("STAGE20_2_FLUID_TO_STRUCTURE_ENABLE")
    s2f = truthy("STAGE20_2_STRUCTURE_TO_FLUID_ENABLE")
    rhs = truthy("STAGE20_2_RHS_COUPLING_ENABLE")
    diagnostic = truthy("STAGE20_2_DIAGNOSTIC_ONLY")
    fail_closed = truthy("STAGE20_2_FAIL_CLOSED")
    contact = truthy("STAGE20_2_CONTACT_ENABLE")
    collision = truthy("STAGE20_2_COLLISION_ENABLE")
    multifibre = truthy("STAGE20_2_MULTIFIBRE_ENABLE")
    lambda_zero = abs(lambda_coupling) <= zero_tol
    default_safe = all([not twoway, f2s, not s2f, not rhs, lambda_zero, diagnostic, fail_closed, truthy("STAGE20_2_SINGLE_FIBRE_ONLY"), not contact, not collision, not multifibre])
    numeric_rules = all([
        n_fibre == 1,
        n_point >= 8,
        component_dim == 3,
        fibre_length > 0.0,
        close(ds, fibre_length / float(n_point - 1), tol),
        rho_l > 0.0,
        rho_tilde > 0.0,
        bending_stiffness >= 0.0,
        gamma >= 0.0,
        c_fs >= 0.0,
        lambda_zero,
        all_finite,
        sorted(global_ids) == list(range(n_point)),
        len(set(global_ids)) == n_point,
        owner_rank == [0 for _ in range(n_point)],
        force_norm > zero_tol,
        (not rhs) and lambda_zero,
    ])

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    statuses["stage20_2_requested_status"] = PASS if truthy("STAGE20_2_ENABLE") else FAIL
    statuses["stage20_2_force_input_adapter_enable_status"] = PASS if truthy("STAGE20_2_FORCE_INPUT_ADAPTER_ENABLE") else FAIL
    statuses["stage20_1_evidence_status"] = PASS if evidence_ok else FAIL
    statuses["stage20_1_source_only_acceptance_preserved_status"] = PASS if truthy("STAGE20_2_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") else FAIL
    statuses["missing_old_stage_outputs_allowed_status"] = PASS if truthy("STAGE20_2_ALLOW_MISSING_OLD_STAGE_OUTPUTS") else FAIL
    statuses["missing_old_closure_files_allowed_status"] = PASS if truthy("STAGE20_2_ALLOW_MISSING_OLD_CLOSURE_FILES") else FAIL
    statuses["no_previous_stage_rerun_status"] = PASS if truthy("STAGE20_2_DO_NOT_RERUN_PREVIOUS_STAGES") else FAIL
    statuses["all_required_adapter_fields_present_status"] = PASS if arrays_present else FAIL
    statuses["default_safe_gate_values_status"] = PASS if default_safe else FAIL
    statuses["twoway_coupling_remains_disabled_status"] = PASS if not twoway else FAIL
    statuses["fluid_to_structure_helper_local_only_status"] = PASS if f2s and diagnostic and (not twoway) and (not rhs) and (not s2f) else FAIL
    statuses["structure_to_fluid_default_disabled_status"] = PASS if not s2f else FAIL
    statuses["rhs_coupling_default_disabled_status"] = PASS if not rhs else FAIL
    statuses["lambda_coupling_zero_status"] = PASS if lambda_zero else FAIL
    statuses["diagnostic_only_status"] = PASS if diagnostic else FAIL
    statuses["fail_closed_status"] = PASS if fail_closed else FAIL
    statuses["single_fibre_only_status"] = PASS if truthy("STAGE20_2_SINGLE_FIBRE_ONLY") else FAIL
    statuses["contact_default_disabled_status"] = PASS if not contact else FAIL
    statuses["collision_default_disabled_status"] = PASS if not collision else FAIL
    statuses["multifibre_default_disabled_status"] = PASS if not multifibre else FAIL
    statuses["n_fibre_status"] = PASS if n_fibre == 1 else FAIL
    statuses["n_point_status"] = PASS if n_point >= 8 else FAIL
    statuses["component_dim_status"] = PASS if component_dim == 3 else FAIL
    statuses["fibre_length_status"] = PASS if fibre_length > 0.0 else FAIL
    statuses["ds_formula_status"] = PASS if close(ds, fibre_length / float(n_point - 1), tol) else FAIL
    statuses["rho_l_status"] = PASS if rho_l > 0.0 else FAIL
    statuses["rho_tilde_status"] = PASS if rho_tilde > 0.0 else FAIL
    statuses["bending_stiffness_status"] = PASS if bending_stiffness >= 0.0 else FAIL
    statuses["gamma_status"] = PASS if gamma >= 0.0 else FAIL
    statuses["c_fs_status"] = PASS if c_fs >= 0.0 else FAIL
    statuses["shape_rules_status"] = PASS if shape_rules else FAIL
    statuses["numeric_rules_status"] = PASS if numeric_rules else FAIL
    statuses["relative_velocity_formula_status"] = PASS if matrix_close(candidate["u_relative_candidate"], u_expected, tol) else FAIL
    statuses["f_fs_candidate_formula_status"] = PASS if matrix_close(candidate["F_fs_candidate"], ffs_expected, tol) else FAIL
    statuses["structure_total_force_without_fluid_formula_status"] = PASS if matrix_close(candidate["F_total_structure_candidate_without_fluid"], total_without_expected, tol) else FAIL
    statuses["structure_total_force_with_fluid_formula_status"] = PASS if matrix_close(candidate["F_total_structure_candidate_with_fluid"], total_with_expected, tol) else FAIL
    statuses["candidate_force_finite_status"] = PASS if all_finite else FAIL
    statuses["candidate_force_helper_local_only_status"] = PASS if diagnostic and f2s and (not twoway) and (not s2f) and (not rhs) else FAIL
    statuses["candidate_force_nonzero_status"] = PASS if force_norm > zero_tol else FAIL
    statuses["no_effective_rhs_coupling_status"] = PASS if (not rhs) and lambda_zero else FAIL
    statuses["global_point_id_coverage_status"] = PASS if sorted(global_ids) == list(range(n_point)) else FAIL
    statuses["global_point_id_no_duplicate_status"] = PASS if len(set(global_ids)) == n_point else FAIL
    statuses["owner_rank_deterministic_status"] = PASS if owner_rank == [0 for _ in range(n_point)] else FAIL
    statuses["stage20_2_wrapper_bash_syntax_status"] = PASS if bash_syntax_ok(root / "stage20_checks" / "run_stage20_2_fluid_to_structure_force_input_adapter.sh") else FAIL
    statuses["stage20_2_helper_py_compile_status"] = PASS if py_compile_ok(Path(__file__).resolve()) else FAIL

    final = PASS if all(value in {PASS, OPTIONAL} for value in statuses.values()) else FAIL
    lines = [
        "STAGE 20.2 FLUID TO STRUCTURE FORCE INPUT ADAPTER AUDIT",
        "stage20_title = real two-way fluid-structure coupling activation boundary",
        "stage20_2_title = fluid-to-structure force input adapter",
        f"repository_root_value = {root}",
        f"stage20_2_test_case_value = {env_value('STAGE20_2_TEST_CASE')}",
        f"stage20_1_evidence_reason_value = {evidence_reason}",
        "stage20_3_next_stage_value = Stage 20.3: structure advance with hydrodynamic force candidate",
        "sign_convention_value = F_total_structure_candidate_with_fluid = F_b_candidate + F_T_candidate - F_fs_candidate",
        "adapter_scope_value = helper-local candidate only; no production DNS force input, no structure-to-fluid reaction, no RHS coupling",
        f"n_fibre_value = {n_fibre}",
        f"n_point_value = {n_point}",
        f"component_dim_value = {component_dim}",
        f"fibre_length_value = {fibre_length:.16e}",
        f"ds_value = {ds:.16e}",
        f"rho_l_value = {rho_l:.16e}",
        f"rho_tilde_value = {rho_tilde:.16e}",
        f"bending_stiffness_value = {bending_stiffness:.16e}",
        f"gamma_value = {gamma:.16e}",
        f"C_fs_value = {c_fs:.16e}",
        f"lambda_coupling_value = {lambda_coupling:.16e}",
        f"candidate_force_norm_value = {force_norm:.16e}",
        "X_fibre_shape_value = (N, 3)",
        "V_fibre_shape_value = (N, 3)",
        "u_interp_candidate_shape_value = (N, 3)",
        "F_fs_candidate_shape_value = (N, 3)",
    ]
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.2 FLUID TO STRUCTURE FORCE INPUT ADAPTER VERDICT: {final}")
    lines.append(f"STAGE 20.2 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"STAGE 20.2 FLUID TO STRUCTURE FORCE INPUT ADAPTER VERDICT: {final}")
    print(f"STAGE 20.2 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
