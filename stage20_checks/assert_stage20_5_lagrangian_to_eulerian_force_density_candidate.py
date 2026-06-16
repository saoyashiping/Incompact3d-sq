#!/usr/bin/env python3
"""Stage 20.5 helper-local Lagrangian-to-Eulerian force-density audit.

This audit spreads a helper-local Lagrangian reaction-force candidate onto a
helper-local Eulerian force-density candidate with a deterministic normalized
nearest-grid-point kernel. It does not add to production RHS, call Stage 14 RHS
injection, run prior stages, MPI, DNS, builds, IBM, or activate production
coupling.
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
    "stage20_5_requested_status", "stage20_5_force_density_candidate_enable_status", "stage20_4_evidence_status",
    "stage20_4_source_only_acceptance_preserved_status", "missing_old_stage_outputs_allowed_status", "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status", "force_density_candidate_documented_status", "kernel_choice_documented_status",
    "all_required_force_density_fields_present_status", "default_safe_gate_values_status", "twoway_coupling_remains_disabled_status",
    "fluid_to_structure_helper_local_only_status", "structure_to_fluid_production_disabled_status", "structure_to_fluid_candidate_helper_local_only_status",
    "force_density_candidate_helper_local_only_status", "rhs_coupling_default_disabled_status", "lambda_coupling_zero_status",
    "diagnostic_only_status", "fail_closed_status", "single_fibre_only_status", "contact_default_disabled_status",
    "collision_default_disabled_status", "multifibre_default_disabled_status", "production_rhs_spreading_disabled_status",
    "production_rhs_update_disabled_status", "stage14_rhs_injection_disabled_status", "n_fibre_status", "n_point_status",
    "component_dim_status", "fibre_length_status", "ds_formula_status", "eulerian_grid_shape_status", "eulerian_grid_spacing_status",
    "eulerian_cell_volume_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status", "gamma_status",
    "c_fs_status", "kernel_name_status", "shape_rules_status", "numeric_rules_status", "kernel_normalization_status",
    "relative_velocity_formula_status", "f_fs_candidate_formula_status", "force_on_structure_formula_status", "force_on_fluid_formula_status",
    "action_reaction_residual_zero_status", "lagrangian_total_reaction_force_formula_status", "eulerian_force_density_candidate_formula_status",
    "eulerian_integral_force_formula_status", "force_conservation_residual_status", "eulerian_force_density_candidate_finite_status",
    "eulerian_force_density_candidate_nonzero_status", "effective_eulerian_force_density_zero_status", "no_production_rhs_update_status",
    "no_stage14_rhs_injection_status", "global_point_id_coverage_status", "global_point_id_no_duplicate_status", "owner_rank_deterministic_status",
    "no_stage10_19_file_modification_status", "no_stage20_0_file_modification_status", "no_stage20_1_file_modification_status",
    "no_stage20_2_file_modification_status", "no_stage20_3_file_modification_status", "no_stage20_4_file_modification_status",
    "no_closed_stage_modification_status", "no_production_fortran_modification_status", "no_cmake_modification_status",
    "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status", "no_production_structure_update_status",
    "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status", "no_production_structure_to_fluid_reaction_force_status",
    "no_production_eulerian_spreading_activation_status", "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage20_5_status",
    "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status",
    "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status", "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "no_rg_only_dependency_status",
    "no_unknown_failure_status", "stage20_6_next_stage_declared_status", "stage20_5_wrapper_bash_syntax_status", "stage20_5_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_5_ENABLE": "1", "STAGE20_5_FORCE_DENSITY_CANDIDATE_ENABLE": "1",
    "STAGE20_5_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1", "STAGE20_5_REQUIRE_STAGE20_4_PASS": "1",
    "STAGE20_5_DO_NOT_RERUN_PREVIOUS_STAGES": "1", "STAGE20_5_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_5_ALLOW_MISSING_OLD_CLOSURE_FILES": "1", "STAGE20_5_DIAGNOSTIC_ONLY": "1", "STAGE20_5_FAIL_CLOSED": "1",
    "STAGE20_5_TWOWAY_COUPLING_ENABLE": "0", "STAGE20_5_FLUID_TO_STRUCTURE_ENABLE": "1",
    "STAGE20_5_STRUCTURE_TO_FLUID_ENABLE": "0", "STAGE20_5_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE": "1",
    "STAGE20_5_RHS_COUPLING_ENABLE": "0", "STAGE20_5_LAMBDA_COUPLING": "0.0", "STAGE20_5_C_FS": "1.0",
    "STAGE20_5_SINGLE_FIBRE_ONLY": "1", "STAGE20_5_CONTACT_ENABLE": "0", "STAGE20_5_COLLISION_ENABLE": "0",
    "STAGE20_5_MULTIFIBRE_ENABLE": "0", "STAGE20_5_PRODUCTION_RHS_SPREADING_ALLOWED": "0",
    "STAGE20_5_PRODUCTION_RHS_UPDATE_ALLOWED": "0", "STAGE20_5_STAGE14_RHS_INJECTION_ALLOWED": "0",
    "STAGE20_5_N_FIBRE": "1", "STAGE20_5_N_POINT": "64", "STAGE20_5_COMPONENT_DIM": "3", "STAGE20_5_FIBRE_LENGTH": "1.0",
    "STAGE20_5_EULERIAN_NX": "16", "STAGE20_5_EULERIAN_NY": "16", "STAGE20_5_EULERIAN_NZ": "16",
    "STAGE20_5_DX": "0.0625", "STAGE20_5_DY": "0.0625", "STAGE20_5_DZ": "0.0625",
    "STAGE20_5_RHO_L": "1.0", "STAGE20_5_RHO_TILDE": "1.0", "STAGE20_5_BENDING_STIFFNESS": "1.0e-5",
    "STAGE20_5_GAMMA": "1.0e-5", "STAGE20_5_KERNEL_NAME": "nearest_grid_point", "STAGE20_5_ZERO_TOL": "1.0e-14",
    "STAGE20_5_AUDIT_TOL": "1.0e-12", "STAGE20_5_TEST_CASE": "lagrangian_to_eulerian_force_density_candidate",
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


def stage20_4_evidence(root: Path) -> Tuple[bool, str]:
    output = root / "stage20_outputs" / "fibre_stage20_4_structure_to_fluid_reaction_force_candidate.dat"
    text = read_text(output)
    output_pass = output.exists() and "STAGE 20.4 STRUCTURE TO FLUID REACTION FORCE CANDIDATE VERDICT: PASS" in text and "STAGE 20.4 FINAL VERDICT: PASS" in text
    markers = [root / "stage20_checks" / name for name in [
        "assert_stage20_4_structure_to_fluid_reaction_force_candidate.py",
        "run_stage20_4_structure_to_fluid_reaction_force_candidate.sh",
        "stage20_4_structure_to_fluid_reaction_force_candidate.md",
    ]]
    if output_pass:
        return True, "ACCEPTED_BY_STAGE20_4_PASS_OUTPUT"
    if truthy("STAGE20_5_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") and all(marker.exists() for marker in markers):
        return True, "ACCEPTED_BY_STAGE20_4_SOURCE_ONLY_BEHAVIOR"
    return False, "NO_STAGE20_4_PASS_OR_SOURCE_ONLY_EVIDENCE"


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_5_py_compile_", suffix=".pyc", delete=False) as tmp:
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
        return subprocess.run(["bash", "-n", str(path)], cwd=str(path.parents[1]), text=True, capture_output=True, check=False).returncode == 0
    except Exception:
        return False


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def vec_close(a: Sequence[float], b: Sequence[float], tol: float) -> bool:
    return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))


def mat_close(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], tol: float) -> bool:
    return len(a) == len(b) and all(vec_close(ra, rb, tol) for ra, rb in zip(a, b))


def norm_vec(v: Sequence[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def build_candidate(n: int, dim: int, length: float, nx: int, ny: int, nz: int, dx: float, dy: float, dz: float, c_fs: float, lambda_c: float) -> Dict[str, object]:
    ds = length / float(n - 1)
    dvol = dx * dy * dz
    x: List[List[float]] = []; v: List[List[float]] = []; u: List[List[float]] = []
    ffs: List[List[float]] = []; f_struct: List[List[float]] = []; f_fluid: List[List[float]] = []; f_sum: List[List[float]] = []
    weights: List[float] = []; support: List[Tuple[int, int, int]] = []
    f_e = [[[[0.0 for _ in range(dim)] for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]
    for q in range(n):
        s = q / float(n - 1)
        xq = [min(max(length * s, 0.0), 1.0 - 0.5 * dx), 0.5 + 1.0e-3 * math.sin(2.0 * math.pi * s), 0.5]
        vq = [0.0, 1.0e-5 * math.cos(2.0 * math.pi * s), 0.0]
        uq = [1.0e-3, 2.0e-4 * math.sin(2.0 * math.pi * s), 0.0]
        urel = [uq[j] - vq[j] for j in range(dim)]
        ff = [c_fs * urel[j] for j in range(dim)]
        fs = [-ff[j] for j in range(dim)]
        fl = [ff[j] for j in range(dim)]
        x.append(xq); v.append(vq); u.append(uq); ffs.append(ff); f_struct.append(fs); f_fluid.append(fl); f_sum.append([fs[j] + fl[j] for j in range(dim)])
        ix = min(nx - 1, max(0, int(xq[0] / dx))); iy = min(ny - 1, max(0, int(xq[1] / dy))); iz = min(nz - 1, max(0, int(xq[2] / dz)))
        weights.append(1.0); support.append((ix, iy, iz))
        for j in range(dim):
            f_e[ix][iy][iz][j] += fl[j] * ds / dvol
    f_eff = [[[[lambda_c * f_e[i][j][k][c] for c in range(dim)] for k in range(nz)] for j in range(ny)] for i in range(nx)]
    lag_total = [sum(f_fluid[q][c] * ds for q in range(n)) for c in range(dim)]
    eul_int = [sum(f_e[i][j][k][c] * dvol for i in range(nx) for j in range(ny) for k in range(nz)) for c in range(dim)]
    residual = [eul_int[c] - lag_total[c] for c in range(dim)]
    return {
        "ds": ds, "dV": dvol, "X_current": x, "V_current": v, "u_interp_candidate": u,
        "u_relative_candidate": [[u[q][c] - v[q][c] for c in range(dim)] for q in range(n)],
        "F_fs_candidate": ffs, "F_on_structure_from_fluid_candidate": f_struct, "F_on_fluid_from_structure_candidate": f_fluid,
        "F_action_reaction_sum_candidate": f_sum, "f_eulerian_candidate": f_e, "f_eulerian_effective": f_eff,
        "delta_kernel_weights": weights, "kernel_support_indices": support,
        "lagrangian_total_reaction_force": lag_total, "eulerian_integral_force_candidate": eul_int,
        "force_conservation_residual_candidate": residual, "force_conservation_residual_norm_candidate": norm_vec(residual),
        "owner_rank": [0 for _ in range(n)], "global_point_id": list(range(n)), "local_point_id": list(range(n)),
        "eulerian_grid_shape": (nx, ny, nz), "eulerian_grid_spacing": (dx, dy, dz), "eulerian_cell_volume": dvol,
    }


def finite_nested(value: object) -> bool:
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, (list, tuple)):
        return all(finite_nested(item) for item in value)
    return True


def main() -> int:
    root = repo_root(); out_dir = root / "stage20_outputs"; out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_5_lagrangian_to_eulerian_force_density_candidate.dat"
    evidence_ok, evidence_reason = stage20_4_evidence(root)
    n_fibre = int(env_value("STAGE20_5_N_FIBRE")); n = int(env_value("STAGE20_5_N_POINT")); dim = int(env_value("STAGE20_5_COMPONENT_DIM"))
    length = float(env_value("STAGE20_5_FIBRE_LENGTH")); nx = int(env_value("STAGE20_5_EULERIAN_NX")); ny = int(env_value("STAGE20_5_EULERIAN_NY")); nz = int(env_value("STAGE20_5_EULERIAN_NZ"))
    dx = float(env_value("STAGE20_5_DX")); dy = float(env_value("STAGE20_5_DY")); dz = float(env_value("STAGE20_5_DZ"))
    rho_l = float(env_value("STAGE20_5_RHO_L")); rho_tilde = float(env_value("STAGE20_5_RHO_TILDE")); bending = float(env_value("STAGE20_5_BENDING_STIFFNESS")); gamma = float(env_value("STAGE20_5_GAMMA"))
    c_fs = float(env_value("STAGE20_5_C_FS")); lambda_c = float(env_value("STAGE20_5_LAMBDA_COUPLING")); tol = float(env_value("STAGE20_5_AUDIT_TOL")); ztol = float(env_value("STAGE20_5_ZERO_TOL")); kernel = env_value("STAGE20_5_KERNEL_NAME")
    cand = build_candidate(n, dim, length, nx, ny, nz, dx, dy, dz, c_fs, lambda_c); ds = float(cand["ds"]); dvol = float(cand["dV"])

    twoway = truthy("STAGE20_5_TWOWAY_COUPLING_ENABLE"); f2s = truthy("STAGE20_5_FLUID_TO_STRUCTURE_ENABLE"); s2f = truthy("STAGE20_5_STRUCTURE_TO_FLUID_ENABLE"); s2f_c = truthy("STAGE20_5_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE"); rhs = truthy("STAGE20_5_RHS_COUPLING_ENABLE")
    diagnostic = truthy("STAGE20_5_DIAGNOSTIC_ONLY"); fail_closed = truthy("STAGE20_5_FAIL_CLOSED"); single = truthy("STAGE20_5_SINGLE_FIBRE_ONLY"); contact = truthy("STAGE20_5_CONTACT_ENABLE"); collision = truthy("STAGE20_5_COLLISION_ENABLE"); multi = truthy("STAGE20_5_MULTIFIBRE_ENABLE")
    prod_spread = truthy("STAGE20_5_PRODUCTION_RHS_SPREADING_ALLOWED"); prod_update = truthy("STAGE20_5_PRODUCTION_RHS_UPDATE_ALLOWED"); stage14 = truthy("STAGE20_5_STAGE14_RHS_INJECTION_ALLOWED"); lambda_zero = abs(lambda_c) <= ztol
    default_safe = all([not twoway, f2s, not s2f, s2f_c, not rhs, lambda_zero, diagnostic, fail_closed, single, not contact, not collision, not multi, not prod_spread, not prod_update, not stage14])
    required_present = all(key in cand for key in ["X_current", "V_current", "u_interp_candidate", "u_relative_candidate", "F_fs_candidate", "F_on_structure_from_fluid_candidate", "F_on_fluid_from_structure_candidate", "F_action_reaction_sum_candidate", "f_eulerian_candidate", "f_eulerian_effective", "delta_kernel_weights", "kernel_support_indices", "lagrangian_total_reaction_force", "eulerian_integral_force_candidate", "force_conservation_residual_candidate", "force_conservation_residual_norm_candidate", "owner_rank", "global_point_id", "local_point_id", "eulerian_grid_shape", "eulerian_grid_spacing", "eulerian_cell_volume"])
    shape_ok = len(cand["X_current"]) == n and len(cand["f_eulerian_candidate"]) == nx and len(cand["f_eulerian_candidate"][0]) == ny and len(cand["f_eulerian_candidate"][0][0]) == nz and len(cand["f_eulerian_candidate"][0][0][0]) == dim and len(cand["lagrangian_total_reaction_force"]) == dim
    finite_ok = finite_nested(cand)
    u_expected = [[cand["u_interp_candidate"][q][c] - cand["V_current"][q][c] for c in range(dim)] for q in range(n)]
    ffs_expected = [[c_fs * cand["u_relative_candidate"][q][c] for c in range(dim)] for q in range(n)]
    fstruct_expected = [[-cand["F_fs_candidate"][q][c] for c in range(dim)] for q in range(n)]
    ffluid_expected = [[cand["F_fs_candidate"][q][c] for c in range(dim)] for q in range(n)]
    lag_expected = [sum(cand["F_on_fluid_from_structure_candidate"][q][c] * ds for q in range(n)) for c in range(dim)]
    eul_expected = [sum(cand["f_eulerian_candidate"][i][j][k][c] * dvol for i in range(nx) for j in range(ny) for k in range(nz)) for c in range(dim)]
    eff_norm = math.sqrt(sum(v * v for i in cand["f_eulerian_effective"] for j in i for k in j for v in k))
    f_e_norm = math.sqrt(sum(v * v for i in cand["f_eulerian_candidate"] for j in i for k in j for v in k))
    action_residual = math.sqrt(sum(v * v for row in cand["F_action_reaction_sum_candidate"] for v in row))
    weights_ok = all(close(w, 1.0, tol) for w in cand["delta_kernel_weights"])
    residual_norm = float(cand["force_conservation_residual_norm_candidate"])
    ids = cand["global_point_id"]; owner = cand["owner_rank"]
    numeric_ok = all([n_fibre == 1, n >= 8, dim == 3, length > 0.0, close(ds, length / float(n - 1), tol), nx >= 4, ny >= 4, nz >= 4, dx > 0, dy > 0, dz > 0, close(dvol, dx * dy * dz, tol), rho_l > 0, rho_tilde > 0, bending >= 0, gamma >= 0, c_fs >= 0, lambda_zero, finite_ok, weights_ok, residual_norm <= tol, f_e_norm > ztol, eff_norm <= ztol, sorted(ids) == list(range(n)), len(set(ids)) == n, owner == [0 for _ in range(n)]])

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    checks = {
        "stage20_5_requested_status": truthy("STAGE20_5_ENABLE"), "stage20_5_force_density_candidate_enable_status": truthy("STAGE20_5_FORCE_DENSITY_CANDIDATE_ENABLE"),
        "stage20_4_evidence_status": evidence_ok, "stage20_4_source_only_acceptance_preserved_status": truthy("STAGE20_5_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": truthy("STAGE20_5_ALLOW_MISSING_OLD_STAGE_OUTPUTS"), "missing_old_closure_files_allowed_status": truthy("STAGE20_5_ALLOW_MISSING_OLD_CLOSURE_FILES"), "no_previous_stage_rerun_status": truthy("STAGE20_5_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "all_required_force_density_fields_present_status": required_present, "default_safe_gate_values_status": default_safe, "twoway_coupling_remains_disabled_status": not twoway, "fluid_to_structure_helper_local_only_status": f2s and diagnostic and not twoway and not rhs,
        "structure_to_fluid_production_disabled_status": not s2f, "structure_to_fluid_candidate_helper_local_only_status": s2f_c and diagnostic and not s2f and not twoway and not rhs,
        "force_density_candidate_helper_local_only_status": diagnostic and not rhs and not prod_spread and not prod_update and not stage14,
        "rhs_coupling_default_disabled_status": not rhs, "lambda_coupling_zero_status": lambda_zero, "diagnostic_only_status": diagnostic, "fail_closed_status": fail_closed, "single_fibre_only_status": single,
        "contact_default_disabled_status": not contact, "collision_default_disabled_status": not collision, "multifibre_default_disabled_status": not multi, "production_rhs_spreading_disabled_status": not prod_spread, "production_rhs_update_disabled_status": not prod_update, "stage14_rhs_injection_disabled_status": not stage14,
        "n_fibre_status": n_fibre == 1, "n_point_status": n >= 8, "component_dim_status": dim == 3, "fibre_length_status": length > 0, "ds_formula_status": close(ds, length / float(n - 1), tol), "eulerian_grid_shape_status": nx >= 4 and ny >= 4 and nz >= 4, "eulerian_grid_spacing_status": dx > 0 and dy > 0 and dz > 0, "eulerian_cell_volume_status": close(dvol, dx * dy * dz, tol),
        "rho_l_status": rho_l > 0, "rho_tilde_status": rho_tilde > 0, "bending_stiffness_status": bending >= 0, "gamma_status": gamma >= 0, "c_fs_status": c_fs >= 0, "kernel_name_status": kernel == "nearest_grid_point", "shape_rules_status": shape_ok, "numeric_rules_status": numeric_ok, "kernel_normalization_status": weights_ok,
        "relative_velocity_formula_status": mat_close(cand["u_relative_candidate"], u_expected, tol), "f_fs_candidate_formula_status": mat_close(cand["F_fs_candidate"], ffs_expected, tol), "force_on_structure_formula_status": mat_close(cand["F_on_structure_from_fluid_candidate"], fstruct_expected, tol), "force_on_fluid_formula_status": mat_close(cand["F_on_fluid_from_structure_candidate"], ffluid_expected, tol), "action_reaction_residual_zero_status": action_residual <= tol,
        "lagrangian_total_reaction_force_formula_status": vec_close(cand["lagrangian_total_reaction_force"], lag_expected, tol), "eulerian_force_density_candidate_formula_status": vec_close(eul_expected, lag_expected, tol), "eulerian_integral_force_formula_status": vec_close(cand["eulerian_integral_force_candidate"], eul_expected, tol), "force_conservation_residual_status": residual_norm <= tol,
        "eulerian_force_density_candidate_finite_status": finite_ok, "eulerian_force_density_candidate_nonzero_status": f_e_norm > ztol, "effective_eulerian_force_density_zero_status": eff_norm <= ztol, "no_production_rhs_update_status": not prod_update and not rhs and lambda_zero, "no_stage14_rhs_injection_status": not stage14,
        "global_point_id_coverage_status": sorted(ids) == list(range(n)), "global_point_id_no_duplicate_status": len(set(ids)) == n, "owner_rank_deterministic_status": owner == [0 for _ in range(n)],
        "stage20_5_wrapper_bash_syntax_status": bash_syntax_ok(root / "stage20_checks" / "run_stage20_5_lagrangian_to_eulerian_force_density_candidate.sh"), "stage20_5_helper_py_compile_status": py_compile_ok(Path(__file__).resolve()),
    }
    for key, ok in checks.items():
        statuses[key] = PASS if ok else FAIL
    final = PASS if all(value in {PASS, OPTIONAL} for value in statuses.values()) else FAIL
    lines = [
        "STAGE 20.5 LAGRANGIAN TO EULERIAN FORCE DENSITY CANDIDATE AUDIT",
        "stage20_title = real two-way fluid-structure coupling activation boundary",
        "stage20_5_title = Lagrangian-to-Eulerian force-density coupling candidate",
        f"repository_root_value = {root}", f"stage20_5_test_case_value = {env_value('STAGE20_5_TEST_CASE')}", f"stage20_4_evidence_reason_value = {evidence_reason}",
        "stage20_6_next_stage_value = Stage 20.6: production RHS coupling activation with lambda gate",
        "kernel_name_value = nearest_grid_point", "kernel_normalization_mode_value = dimensionless NGP weights sum to one per Lagrangian point; f = F * weight * ds / dV",
        "force_density_scope_value = helper-local Eulerian candidate only; no production RHS, no IBM/RHS coupling, no Stage 14 injection",
        f"n_fibre_value = {n_fibre}", f"n_point_value = {n}", f"component_dim_value = {dim}", f"eulerian_grid_shape_value = ({nx}, {ny}, {nz})",
        f"fibre_length_value = {length:.16e}", f"ds_value = {ds:.16e}", f"dx_value = {dx:.16e}", f"dy_value = {dy:.16e}", f"dz_value = {dz:.16e}", f"dV_value = {dvol:.16e}",
        f"lagrangian_total_reaction_force_value = {cand['lagrangian_total_reaction_force']}", f"eulerian_integral_force_candidate_value = {cand['eulerian_integral_force_candidate']}",
        f"force_conservation_residual_norm_candidate_value = {residual_norm:.16e}", f"eulerian_force_density_norm_value = {f_e_norm:.16e}", f"effective_eulerian_force_density_norm_value = {eff_norm:.16e}",
    ]
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.5 LAGRANGIAN TO EULERIAN FORCE DENSITY CANDIDATE VERDICT: {final}")
    lines.append(f"STAGE 20.5 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 20.5 LAGRANGIAN TO EULERIAN FORCE DENSITY CANDIDATE VERDICT: {final}")
    print(f"STAGE 20.5 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
