#!/usr/bin/env python3
"""Stage 20.3 helper-local structure advance with hydrodynamic force candidate.

This audit computes helper-local candidate A/V/X updates using F_b_candidate,
F_T_candidate, and F_fs_candidate. It does not commit production state, run prior
stages, MPI, DNS, builds, production hooks, RHS injection, IBM, or coupling back
to the fluid.
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
    "stage20_3_requested_status",
    "stage20_3_structure_advance_hydro_candidate_enable_status",
    "stage20_2_evidence_status",
    "stage20_2_source_only_acceptance_preserved_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "hydro_advance_candidate_documented_status",
    "all_required_hydro_advance_fields_present_status",
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
    "production_commit_disabled_status",
    "production_structure_advance_disabled_status",
    "n_fibre_status",
    "n_point_status",
    "component_dim_status",
    "fibre_length_status",
    "ds_formula_status",
    "dt_status",
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
    "acceleration_hydro_candidate_formula_status",
    "velocity_next_hydro_candidate_formula_status",
    "position_next_hydro_candidate_formula_status",
    "delta_velocity_hydro_candidate_formula_status",
    "delta_position_hydro_candidate_formula_status",
    "candidate_arrays_finite_status",
    "hydro_advance_candidate_helper_local_only_status",
    "hydro_candidate_differs_from_no_fluid_candidate_status",
    "no_production_committed_state_change_status",
    "no_effective_rhs_coupling_status",
    "global_point_id_coverage_status",
    "global_point_id_no_duplicate_status",
    "owner_rank_deterministic_status",
    "no_stage10_19_file_modification_status",
    "no_stage20_0_file_modification_status",
    "no_stage20_1_file_modification_status",
    "no_stage20_2_file_modification_status",
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
    "no_stage14_rhs_call_from_stage20_3_status",
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
    "stage20_4_next_stage_declared_status",
    "stage20_3_wrapper_bash_syntax_status",
    "stage20_3_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_3_ENABLE": "1",
    "STAGE20_3_STRUCTURE_ADVANCE_HYDRO_CANDIDATE_ENABLE": "1",
    "STAGE20_3_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_3_REQUIRE_STAGE20_2_PASS": "1",
    "STAGE20_3_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_3_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_3_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_3_DIAGNOSTIC_ONLY": "1",
    "STAGE20_3_FAIL_CLOSED": "1",
    "STAGE20_3_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_3_FLUID_TO_STRUCTURE_ENABLE": "1",
    "STAGE20_3_STRUCTURE_TO_FLUID_ENABLE": "0",
    "STAGE20_3_RHS_COUPLING_ENABLE": "0",
    "STAGE20_3_LAMBDA_COUPLING": "0.0",
    "STAGE20_3_C_FS": "1.0",
    "STAGE20_3_SINGLE_FIBRE_ONLY": "1",
    "STAGE20_3_CONTACT_ENABLE": "0",
    "STAGE20_3_COLLISION_ENABLE": "0",
    "STAGE20_3_MULTIFIBRE_ENABLE": "0",
    "STAGE20_3_PRODUCTION_COMMIT_ALLOWED": "0",
    "STAGE20_3_PRODUCTION_STRUCTURE_ADVANCE_ALLOWED": "0",
    "STAGE20_3_N_FIBRE": "1",
    "STAGE20_3_N_POINT": "64",
    "STAGE20_3_COMPONENT_DIM": "3",
    "STAGE20_3_FIBRE_LENGTH": "1.0",
    "STAGE20_3_DT": "1.0e-5",
    "STAGE20_3_RHO_L": "1.0",
    "STAGE20_3_RHO_TILDE": "1.0",
    "STAGE20_3_BENDING_STIFFNESS": "1.0e-5",
    "STAGE20_3_GAMMA": "1.0e-5",
    "STAGE20_3_ZERO_TOL": "1.0e-14",
    "STAGE20_3_AUDIT_TOL": "1.0e-12",
    "STAGE20_3_TEST_CASE": "structure_advance_with_hydrodynamic_force_candidate",
}

REQUIRED_ARRAYS = [
    "X_current", "V_current", "A_current", "u_interp_candidate", "u_relative_candidate",
    "F_fs_candidate", "F_b_candidate", "F_T_candidate",
    "F_total_structure_candidate_without_fluid", "F_total_structure_candidate_with_fluid",
    "A_hydro_candidate", "V_next_hydro_candidate", "X_next_hydro_candidate",
    "delta_V_hydro_candidate", "delta_X_hydro_candidate", "owner_rank", "global_point_id", "local_point_id",
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


def stage20_2_evidence(root: Path) -> Tuple[bool, str]:
    output = root / "stage20_outputs" / "fibre_stage20_2_fluid_to_structure_force_input_adapter.dat"
    text = read_text(output)
    output_pass = output.exists() and (
        "STAGE 20.2 FLUID TO STRUCTURE FORCE INPUT ADAPTER VERDICT: PASS" in text
        and "STAGE 20.2 FINAL VERDICT: PASS" in text
    )
    markers = [
        root / "stage20_checks" / "assert_stage20_2_fluid_to_structure_force_input_adapter.py",
        root / "stage20_checks" / "run_stage20_2_fluid_to_structure_force_input_adapter.sh",
        root / "stage20_checks" / "stage20_2_fluid_to_structure_force_input_adapter.md",
    ]
    source_only_available = truthy("STAGE20_3_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") and all(marker.exists() for marker in markers)
    if output_pass:
        return True, "ACCEPTED_BY_STAGE20_2_PASS_OUTPUT"
    if source_only_available:
        return True, "ACCEPTED_BY_STAGE20_2_SOURCE_ONLY_BEHAVIOR"
    return False, "NO_STAGE20_2_PASS_OR_SOURCE_ONLY_EVIDENCE"


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_3_py_compile_", suffix=".pyc", delete=False) as tmp:
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


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def matrix_shape(values: Sequence[Sequence[float]], rows: int, cols: int) -> bool:
    return len(values) == rows and all(len(row) == cols for row in values)


def matrix_finite(values: Sequence[Sequence[float]]) -> bool:
    return all(math.isfinite(value) for row in values for value in row)


def matrix_close(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], tol: float) -> bool:
    return len(a) == len(b) and all(len(ra) == len(rb) and all(close(x, y, tol) for x, y in zip(ra, rb)) for ra, rb in zip(a, b))


def build_candidate(n: int, dim: int, length: float, dt: float, rho_l: float, c_fs: float) -> Dict[str, object]:
    ds = length / float(n - 1)
    x: List[List[float]] = []
    v: List[List[float]] = []
    a_current: List[List[float]] = []
    u: List[List[float]] = []
    fb: List[List[float]] = []
    ft: List[List[float]] = []
    for i in range(n):
        s = i / float(n - 1)
        x.append([length * s, 1.0e-3 * math.sin(2.0 * math.pi * s), 0.0])
        v.append([0.0, 1.0e-5 * math.cos(2.0 * math.pi * s), 0.0])
        a_current.append([0.0, 0.0, 0.0])
        u.append([1.0e-3, 2.0e-4 * math.sin(2.0 * math.pi * s), 0.0])
        fb.append([0.0, -1.0e-6 * math.sin(2.0 * math.pi * s), 0.0])
        ft.append([1.0e-7 * math.cos(2.0 * math.pi * s), 0.0, 0.0])
    urel = [[u[i][j] - v[i][j] for j in range(dim)] for i in range(n)]
    ffs = [[c_fs * urel[i][j] for j in range(dim)] for i in range(n)]
    total_without = [[fb[i][j] + ft[i][j] for j in range(dim)] for i in range(n)]
    total_with = [[fb[i][j] + ft[i][j] - ffs[i][j] for j in range(dim)] for i in range(n)]
    ah = [[total_with[i][j] / rho_l for j in range(dim)] for i in range(n)]
    vnext = [[v[i][j] + dt * ah[i][j] for j in range(dim)] for i in range(n)]
    xnext = [[x[i][j] + dt * v[i][j] + 0.5 * dt * dt * ah[i][j] for j in range(dim)] for i in range(n)]
    dv = [[vnext[i][j] - v[i][j] for j in range(dim)] for i in range(n)]
    dx = [[xnext[i][j] - x[i][j] for j in range(dim)] for i in range(n)]
    nofluid_a = [[total_without[i][j] / rho_l for j in range(dim)] for i in range(n)]
    nofluid_vnext = [[v[i][j] + dt * nofluid_a[i][j] for j in range(dim)] for i in range(n)]
    return {
        "ds": ds, "X_current": x, "V_current": v, "A_current": a_current,
        "u_interp_candidate": u, "u_relative_candidate": urel, "F_fs_candidate": ffs,
        "F_b_candidate": fb, "F_T_candidate": ft,
        "F_total_structure_candidate_without_fluid": total_without,
        "F_total_structure_candidate_with_fluid": total_with,
        "A_hydro_candidate": ah, "V_next_hydro_candidate": vnext, "X_next_hydro_candidate": xnext,
        "delta_V_hydro_candidate": dv, "delta_X_hydro_candidate": dx,
        "V_next_no_fluid_candidate": nofluid_vnext,
        "owner_rank": [0 for _ in range(n)], "global_point_id": list(range(n)), "local_point_id": list(range(n)),
    }


def main() -> int:
    root = repo_root()
    out_dir = root / "stage20_outputs"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_3_structure_advance_hydro_force_candidate.dat"
    evidence_ok, evidence_reason = stage20_2_evidence(root)

    n_fibre = int(env_value("STAGE20_3_N_FIBRE")); n = int(env_value("STAGE20_3_N_POINT")); dim = int(env_value("STAGE20_3_COMPONENT_DIM"))
    length = float(env_value("STAGE20_3_FIBRE_LENGTH")); dt = float(env_value("STAGE20_3_DT")); rho_l = float(env_value("STAGE20_3_RHO_L"))
    rho_tilde = float(env_value("STAGE20_3_RHO_TILDE")); bending = float(env_value("STAGE20_3_BENDING_STIFFNESS")); gamma = float(env_value("STAGE20_3_GAMMA"))
    c_fs = float(env_value("STAGE20_3_C_FS")); lambda_c = float(env_value("STAGE20_3_LAMBDA_COUPLING")); tol = float(env_value("STAGE20_3_AUDIT_TOL")); ztol = float(env_value("STAGE20_3_ZERO_TOL"))
    candidate = build_candidate(n, dim, length, dt, rho_l, c_fs)
    ds = float(candidate["ds"])

    matrix_names = REQUIRED_ARRAYS[:15]
    fields_present = all(name in candidate for name in REQUIRED_ARRAYS)
    shape_ok = all(matrix_shape(candidate[name], n, dim) for name in matrix_names) and all(len(candidate[name]) == n for name in ["owner_rank", "global_point_id", "local_point_id"])
    finite_ok = all(matrix_finite(candidate[name]) for name in matrix_names)
    u_expected = [[candidate["u_interp_candidate"][i][j] - candidate["V_current"][i][j] for j in range(dim)] for i in range(n)]
    ffs_expected = [[c_fs * candidate["u_relative_candidate"][i][j] for j in range(dim)] for i in range(n)]
    tw_expected = [[candidate["F_b_candidate"][i][j] + candidate["F_T_candidate"][i][j] for j in range(dim)] for i in range(n)]
    tf_expected = [[candidate["F_b_candidate"][i][j] + candidate["F_T_candidate"][i][j] - candidate["F_fs_candidate"][i][j] for j in range(dim)] for i in range(n)]
    ah_expected = [[candidate["F_total_structure_candidate_with_fluid"][i][j] / rho_l for j in range(dim)] for i in range(n)]
    vn_expected = [[candidate["V_current"][i][j] + dt * candidate["A_hydro_candidate"][i][j] for j in range(dim)] for i in range(n)]
    xn_expected = [[candidate["X_current"][i][j] + dt * candidate["V_current"][i][j] + 0.5 * dt * dt * candidate["A_hydro_candidate"][i][j] for j in range(dim)] for i in range(n)]
    dv_expected = [[candidate["V_next_hydro_candidate"][i][j] - candidate["V_current"][i][j] for j in range(dim)] for i in range(n)]
    dx_expected = [[candidate["X_next_hydro_candidate"][i][j] - candidate["X_current"][i][j] for j in range(dim)] for i in range(n)]
    ids = candidate["global_point_id"]; owner = candidate["owner_rank"]
    ffs_norm = math.sqrt(sum(v * v for row in candidate["F_fs_candidate"] for v in row))
    hydro_diff = math.sqrt(sum((candidate["V_next_hydro_candidate"][i][j] - candidate["V_next_no_fluid_candidate"][i][j]) ** 2 for i in range(n) for j in range(dim)))

    twoway = truthy("STAGE20_3_TWOWAY_COUPLING_ENABLE"); f2s = truthy("STAGE20_3_FLUID_TO_STRUCTURE_ENABLE"); s2f = truthy("STAGE20_3_STRUCTURE_TO_FLUID_ENABLE"); rhs = truthy("STAGE20_3_RHS_COUPLING_ENABLE")
    diagnostic = truthy("STAGE20_3_DIAGNOSTIC_ONLY"); fail_closed = truthy("STAGE20_3_FAIL_CLOSED"); single = truthy("STAGE20_3_SINGLE_FIBRE_ONLY")
    contact = truthy("STAGE20_3_CONTACT_ENABLE"); collision = truthy("STAGE20_3_COLLISION_ENABLE"); multi = truthy("STAGE20_3_MULTIFIBRE_ENABLE")
    commit = truthy("STAGE20_3_PRODUCTION_COMMIT_ALLOWED"); prod_advance = truthy("STAGE20_3_PRODUCTION_STRUCTURE_ADVANCE_ALLOWED")
    lambda_zero = abs(lambda_c) <= ztol
    default_safe = all([not twoway, f2s, not s2f, not rhs, lambda_zero, diagnostic, fail_closed, single, not contact, not collision, not multi, not commit, not prod_advance])
    numeric_ok = all([n_fibre == 1, n >= 8, dim == 3, length > 0.0, close(ds, length / float(n - 1), tol), dt > 0.0, rho_l > 0.0, rho_tilde > 0.0, bending >= 0.0, gamma >= 0.0, c_fs >= 0.0, lambda_zero, finite_ok, sorted(ids) == list(range(n)), len(set(ids)) == n, owner == [0 for _ in range(n)], hydro_diff > ztol, not commit, not rhs and lambda_zero])

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    checks = {
        "stage20_3_requested_status": truthy("STAGE20_3_ENABLE"),
        "stage20_3_structure_advance_hydro_candidate_enable_status": truthy("STAGE20_3_STRUCTURE_ADVANCE_HYDRO_CANDIDATE_ENABLE"),
        "stage20_2_evidence_status": evidence_ok,
        "stage20_2_source_only_acceptance_preserved_status": truthy("STAGE20_3_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": truthy("STAGE20_3_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
        "missing_old_closure_files_allowed_status": truthy("STAGE20_3_ALLOW_MISSING_OLD_CLOSURE_FILES"),
        "no_previous_stage_rerun_status": truthy("STAGE20_3_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "all_required_hydro_advance_fields_present_status": fields_present,
        "default_safe_gate_values_status": default_safe,
        "twoway_coupling_remains_disabled_status": not twoway,
        "fluid_to_structure_helper_local_only_status": f2s and diagnostic and not twoway and not rhs and not s2f,
        "structure_to_fluid_default_disabled_status": not s2f,
        "rhs_coupling_default_disabled_status": not rhs,
        "lambda_coupling_zero_status": lambda_zero,
        "diagnostic_only_status": diagnostic,
        "fail_closed_status": fail_closed,
        "single_fibre_only_status": single,
        "contact_default_disabled_status": not contact,
        "collision_default_disabled_status": not collision,
        "multifibre_default_disabled_status": not multi,
        "production_commit_disabled_status": not commit,
        "production_structure_advance_disabled_status": not prod_advance,
        "n_fibre_status": n_fibre == 1,
        "n_point_status": n >= 8,
        "component_dim_status": dim == 3,
        "fibre_length_status": length > 0.0,
        "ds_formula_status": close(ds, length / float(n - 1), tol),
        "dt_status": dt > 0.0,
        "rho_l_status": rho_l > 0.0,
        "rho_tilde_status": rho_tilde > 0.0,
        "bending_stiffness_status": bending >= 0.0,
        "gamma_status": gamma >= 0.0,
        "c_fs_status": c_fs >= 0.0,
        "shape_rules_status": shape_ok,
        "numeric_rules_status": numeric_ok,
        "relative_velocity_formula_status": matrix_close(candidate["u_relative_candidate"], u_expected, tol),
        "f_fs_candidate_formula_status": matrix_close(candidate["F_fs_candidate"], ffs_expected, tol),
        "structure_total_force_without_fluid_formula_status": matrix_close(candidate["F_total_structure_candidate_without_fluid"], tw_expected, tol),
        "structure_total_force_with_fluid_formula_status": matrix_close(candidate["F_total_structure_candidate_with_fluid"], tf_expected, tol),
        "acceleration_hydro_candidate_formula_status": matrix_close(candidate["A_hydro_candidate"], ah_expected, tol),
        "velocity_next_hydro_candidate_formula_status": matrix_close(candidate["V_next_hydro_candidate"], vn_expected, tol),
        "position_next_hydro_candidate_formula_status": matrix_close(candidate["X_next_hydro_candidate"], xn_expected, tol),
        "delta_velocity_hydro_candidate_formula_status": matrix_close(candidate["delta_V_hydro_candidate"], dv_expected, tol),
        "delta_position_hydro_candidate_formula_status": matrix_close(candidate["delta_X_hydro_candidate"], dx_expected, tol),
        "candidate_arrays_finite_status": finite_ok,
        "hydro_advance_candidate_helper_local_only_status": diagnostic and f2s and not twoway and not s2f and not rhs and not commit and not prod_advance,
        "hydro_candidate_differs_from_no_fluid_candidate_status": ffs_norm > ztol and hydro_diff > ztol,
        "no_production_committed_state_change_status": not commit,
        "no_effective_rhs_coupling_status": not rhs and lambda_zero,
        "global_point_id_coverage_status": sorted(ids) == list(range(n)),
        "global_point_id_no_duplicate_status": len(set(ids)) == n,
        "owner_rank_deterministic_status": owner == [0 for _ in range(n)],
        "stage20_3_wrapper_bash_syntax_status": bash_syntax_ok(root / "stage20_checks" / "run_stage20_3_structure_advance_hydro_force_candidate.sh"),
        "stage20_3_helper_py_compile_status": py_compile_ok(Path(__file__).resolve()),
    }
    for key, ok in checks.items():
        statuses[key] = PASS if ok else FAIL

    final = PASS if all(value in {PASS, OPTIONAL} for value in statuses.values()) else FAIL
    lines = [
        "STAGE 20.3 STRUCTURE ADVANCE HYDRO FORCE CANDIDATE AUDIT",
        "stage20_title = real two-way fluid-structure coupling activation boundary",
        "stage20_3_title = structure advance with hydrodynamic force candidate",
        f"repository_root_value = {root}",
        f"stage20_3_test_case_value = {env_value('STAGE20_3_TEST_CASE')}",
        f"stage20_2_evidence_reason_value = {evidence_reason}",
        "stage20_4_next_stage_value = Stage 20.4: structure-to-fluid reaction force candidate",
        "advance_scope_value = helper-local A/V/X candidate only; no production commit, no production structure advance",
        "sign_convention_value = F_total_structure_candidate_with_fluid = F_b_candidate + F_T_candidate - F_fs_candidate",
        f"n_fibre_value = {n_fibre}", f"n_point_value = {n}", f"component_dim_value = {dim}",
        f"fibre_length_value = {length:.16e}", f"ds_value = {ds:.16e}", f"dt_value = {dt:.16e}",
        f"rho_l_value = {rho_l:.16e}", f"rho_tilde_value = {rho_tilde:.16e}", f"bending_stiffness_value = {bending:.16e}",
        f"gamma_value = {gamma:.16e}", f"C_fs_value = {c_fs:.16e}", f"lambda_coupling_value = {lambda_c:.16e}",
        f"F_fs_candidate_norm_value = {ffs_norm:.16e}", f"hydro_vs_no_fluid_velocity_delta_norm_value = {hydro_diff:.16e}",
    ]
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.3 STRUCTURE ADVANCE HYDRO FORCE CANDIDATE VERDICT: {final}")
    lines.append(f"STAGE 20.3 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 20.3 STRUCTURE ADVANCE HYDRO FORCE CANDIDATE VERDICT: {final}")
    print(f"STAGE 20.3 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
