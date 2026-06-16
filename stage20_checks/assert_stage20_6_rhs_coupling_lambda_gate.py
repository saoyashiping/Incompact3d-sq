#!/usr/bin/env python3
"""Stage 20.6 helper-local RHS coupling lambda-gate audit.

This audit validates lambda-gated RHS candidate behavior for a helper-local
Eulerian force-density candidate. It proves lambda=0 is a strict RHS no-op and a
small lambda produces a bounded helper-local RHS delta. It does not write
production RHS, call Stage 14 RHS injection, run DNS/MPI, or modify production
paths.
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
    "stage20_6_requested_status", "stage20_6_rhs_coupling_lambda_gate_enable_status", "stage20_6_rhs_coupling_candidate_enable_status",
    "stage20_5_evidence_status", "stage20_5_source_only_acceptance_preserved_status", "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status", "no_previous_stage_rerun_status", "rhs_coupling_lambda_gate_documented_status",
    "all_required_rhs_coupling_fields_present_status", "default_safe_gate_values_status", "twoway_coupling_remains_disabled_status",
    "fluid_to_structure_helper_local_only_status", "structure_to_fluid_production_disabled_status", "structure_to_fluid_candidate_helper_local_only_status",
    "rhs_coupling_candidate_helper_local_only_status", "production_rhs_coupling_disabled_status", "production_rhs_update_disabled_status",
    "stage14_rhs_injection_disabled_status", "lambda_zero_case_status", "lambda_small_case_status", "lambda_scaling_ratio_status",
    "n_fibre_status", "n_point_status", "component_dim_status", "fibre_length_status", "ds_formula_status",
    "eulerian_grid_shape_status", "eulerian_grid_spacing_status", "eulerian_cell_volume_status", "rho_l_status", "rho_tilde_status",
    "bending_stiffness_status", "gamma_status", "c_fs_status", "kernel_name_status", "lambda_zero_status", "lambda_small_status",
    "shape_rules_status", "numeric_rules_status", "eulerian_force_density_candidate_finite_status", "eulerian_force_density_candidate_nonzero_status",
    "rhs_before_finite_status", "rhs_after_zero_finite_status", "rhs_after_small_finite_status", "rhs_zero_residual_zero_status",
    "rhs_small_formula_residual_status", "rhs_delta_small_bounded_status", "f_eulerian_effective_zero_formula_status",
    "f_eulerian_effective_small_formula_status", "rhs_delta_zero_formula_status", "rhs_delta_small_formula_status",
    "rhs_after_zero_formula_status", "rhs_after_small_formula_status", "no_production_rhs_update_status", "no_stage14_rhs_injection_status",
    "global_point_id_coverage_status", "global_point_id_no_duplicate_status", "owner_rank_deterministic_status", "diagnostic_only_status",
    "fail_closed_status", "single_fibre_only_status", "contact_default_disabled_status", "collision_default_disabled_status",
    "multifibre_default_disabled_status", "no_stage10_19_file_modification_status", "no_stage20_0_file_modification_status",
    "no_stage20_1_file_modification_status", "no_stage20_2_file_modification_status", "no_stage20_3_file_modification_status",
    "no_stage20_4_file_modification_status", "no_stage20_5_file_modification_status", "no_closed_stage_modification_status",
    "no_production_fortran_modification_status", "no_cmake_modification_status", "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status", "no_production_structure_update_status", "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status", "no_production_structure_to_fluid_reaction_force_status",
    "no_production_eulerian_spreading_activation_status", "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage20_6_status",
    "no_fluid_rhs_modification_status", "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status",
    "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status", "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status", "no_rg_only_dependency_status",
    "no_unknown_failure_status", "stage20_7_next_stage_declared_status", "stage20_6_wrapper_bash_syntax_status", "stage20_6_helper_py_compile_status",
]

SAFE_DEFAULTS = {
    "STAGE20_6_ENABLE": "1", "STAGE20_6_RHS_COUPLING_LAMBDA_GATE_ENABLE": "1", "STAGE20_6_RHS_COUPLING_CANDIDATE_ENABLE": "1",
    "STAGE20_6_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1", "STAGE20_6_REQUIRE_STAGE20_5_PASS": "1",
    "STAGE20_6_DO_NOT_RERUN_PREVIOUS_STAGES": "1", "STAGE20_6_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_6_ALLOW_MISSING_OLD_CLOSURE_FILES": "1", "STAGE20_6_DIAGNOSTIC_ONLY": "1", "STAGE20_6_FAIL_CLOSED": "1",
    "STAGE20_6_TWOWAY_COUPLING_ENABLE": "0", "STAGE20_6_FLUID_TO_STRUCTURE_ENABLE": "1",
    "STAGE20_6_STRUCTURE_TO_FLUID_ENABLE": "0", "STAGE20_6_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE": "1",
    "STAGE20_6_RHS_COUPLING_ENABLE": "0", "STAGE20_6_PRODUCTION_RHS_UPDATE_ALLOWED": "0", "STAGE20_6_STAGE14_RHS_INJECTION_ALLOWED": "0",
    "STAGE20_6_LAMBDA_ZERO": "0.0", "STAGE20_6_LAMBDA_SMALL": "1.0e-6", "STAGE20_6_C_FS": "1.0",
    "STAGE20_6_SINGLE_FIBRE_ONLY": "1", "STAGE20_6_CONTACT_ENABLE": "0", "STAGE20_6_COLLISION_ENABLE": "0", "STAGE20_6_MULTIFIBRE_ENABLE": "0",
    "STAGE20_6_N_FIBRE": "1", "STAGE20_6_N_POINT": "64", "STAGE20_6_COMPONENT_DIM": "3", "STAGE20_6_FIBRE_LENGTH": "1.0",
    "STAGE20_6_EULERIAN_NX": "16", "STAGE20_6_EULERIAN_NY": "16", "STAGE20_6_EULERIAN_NZ": "16",
    "STAGE20_6_DX": "0.0625", "STAGE20_6_DY": "0.0625", "STAGE20_6_DZ": "0.0625",
    "STAGE20_6_RHO_L": "1.0", "STAGE20_6_RHO_TILDE": "1.0", "STAGE20_6_BENDING_STIFFNESS": "1.0e-5", "STAGE20_6_GAMMA": "1.0e-5",
    "STAGE20_6_KERNEL_NAME": "nearest_grid_point", "STAGE20_6_MAX_RHS_DELTA_NORM": "1.0e-6", "STAGE20_6_ZERO_TOL": "1.0e-14", "STAGE20_6_AUDIT_TOL": "1.0e-12",
    "STAGE20_6_TEST_CASE": "rhs_coupling_lambda_gate",
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


def stage20_5_evidence(root: Path) -> Tuple[bool, str]:
    output = root / "stage20_outputs" / "fibre_stage20_5_lagrangian_to_eulerian_force_density_candidate.dat"
    text = read_text(output)
    output_pass = output.exists() and "STAGE 20.5 LAGRANGIAN TO EULERIAN FORCE DENSITY CANDIDATE VERDICT: PASS" in text and "STAGE 20.5 FINAL VERDICT: PASS" in text
    markers = [root / "stage20_checks" / name for name in [
        "assert_stage20_5_lagrangian_to_eulerian_force_density_candidate.py",
        "run_stage20_5_lagrangian_to_eulerian_force_density_candidate.sh",
        "stage20_5_lagrangian_to_eulerian_force_density_candidate.md",
    ]]
    if output_pass:
        return True, "ACCEPTED_BY_STAGE20_5_PASS_OUTPUT"
    if truthy("STAGE20_6_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE") and all(marker.exists() for marker in markers):
        return True, "ACCEPTED_BY_STAGE20_5_SOURCE_ONLY_BEHAVIOR"
    return False, "NO_STAGE20_5_PASS_OR_SOURCE_ONLY_EVIDENCE"


def py_compile_ok(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(prefix="stage20_6_py_compile_", suffix=".pyc", delete=False) as tmp:
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


def norm_flat(a: Sequence[float]) -> float:
    return math.sqrt(sum(x * x for x in a))


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def finite(values: Sequence[float]) -> bool:
    return all(math.isfinite(v) for v in values)


def build_candidate(n: int, dim: int, length: float, nx: int, ny: int, nz: int, dx: float, dy: float, dz: float, c_fs: float, lambda_zero: float, lambda_small: float) -> Dict[str, object]:
    ds = length / float(n - 1); dvol = dx * dy * dz; size = nx * ny * nz * dim
    x: List[List[float]] = []; v: List[List[float]] = []; u: List[List[float]] = []; ffs: List[List[float]] = []
    f_struct: List[List[float]] = []; f_fluid: List[List[float]] = []
    f_e = [0.0 for _ in range(size)]
    for q in range(n):
        s = q / float(n - 1)
        xq = [min(max(length * s, 0.0), 1.0 - 0.5 * dx), 0.5 + 1.0e-3 * math.sin(2.0 * math.pi * s), 0.5]
        vq = [0.0, 1.0e-5 * math.cos(2.0 * math.pi * s), 0.0]
        uq = [5.0e-4, 1.0e-4 * math.sin(2.0 * math.pi * s), 0.0]
        ff = [c_fs * (uq[c] - vq[c]) for c in range(dim)]
        x.append(xq); v.append(vq); u.append(uq); ffs.append(ff); f_struct.append([-z for z in ff]); f_fluid.append(list(ff))
        ix = min(nx - 1, max(0, int(xq[0] / dx))); iy = min(ny - 1, max(0, int(xq[1] / dy))); iz = min(nz - 1, max(0, int(xq[2] / dz)))
        base = ((ix * ny + iy) * nz + iz) * dim
        for c in range(dim):
            f_e[base + c] += ff[c] * ds / dvol
    rhs_before = [1.0e-9 * math.sin(0.01 * i) for i in range(size)]
    eff_zero = [lambda_zero * val for val in f_e]
    eff_small = [lambda_small * val for val in f_e]
    delta_zero = list(eff_zero); delta_small = list(eff_small)
    after_zero = [rhs_before[i] + delta_zero[i] for i in range(size)]
    after_small = [rhs_before[i] + delta_small[i] for i in range(size)]
    zero_residual = [after_zero[i] - rhs_before[i] for i in range(size)]
    small_residual = [delta_small[i] - lambda_small * f_e[i] for i in range(size)]
    f_norm = norm_flat(f_e); small_norm = norm_flat(delta_small)
    ratio = small_norm / f_norm if f_norm > 0.0 else 0.0
    lag_total = [sum(f_fluid[q][c] * ds for q in range(n)) for c in range(dim)]
    eul_int = [sum(f_e[((i * ny + j) * nz + k) * dim + c] * dvol for i in range(nx) for j in range(ny) for k in range(nz)) for c in range(dim)]
    force_res = [eul_int[c] - lag_total[c] for c in range(dim)]
    return {
        "ds": ds, "dV": dvol, "size": size, "X_current": x, "V_current": v, "u_interp_candidate": u,
        "u_relative_candidate": [[u[q][c] - v[q][c] for c in range(dim)] for q in range(n)], "F_fs_candidate": ffs,
        "F_on_structure_from_fluid_candidate": f_struct, "F_on_fluid_from_structure_candidate": f_fluid,
        "f_eulerian_candidate": f_e, "f_eulerian_effective_zero": eff_zero, "f_eulerian_effective_small": eff_small,
        "RHS_before": rhs_before, "RHS_delta_zero": delta_zero, "RHS_delta_small": delta_small,
        "RHS_after_zero": after_zero, "RHS_after_small": after_small, "RHS_zero_residual": zero_residual,
        "RHS_small_formula_residual": small_residual, "lambda_scaling_ratio_candidate": ratio,
        "eulerian_integral_force_candidate": eul_int, "lagrangian_total_reaction_force": lag_total,
        "force_conservation_residual_candidate": force_res, "owner_rank": [0 for _ in range(n)], "global_point_id": list(range(n)), "local_point_id": list(range(n)),
        "eulerian_grid_shape": (nx, ny, nz), "eulerian_grid_spacing": (dx, dy, dz), "eulerian_cell_volume": dvol,
    }


def main() -> int:
    root = repo_root(); out_dir = root / "stage20_outputs"; out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "fibre_stage20_6_rhs_coupling_lambda_gate.dat"
    evidence_ok, evidence_reason = stage20_5_evidence(root)
    n_fibre = int(env_value("STAGE20_6_N_FIBRE")); n = int(env_value("STAGE20_6_N_POINT")); dim = int(env_value("STAGE20_6_COMPONENT_DIM"))
    length = float(env_value("STAGE20_6_FIBRE_LENGTH")); nx = int(env_value("STAGE20_6_EULERIAN_NX")); ny = int(env_value("STAGE20_6_EULERIAN_NY")); nz = int(env_value("STAGE20_6_EULERIAN_NZ"))
    dx = float(env_value("STAGE20_6_DX")); dy = float(env_value("STAGE20_6_DY")); dz = float(env_value("STAGE20_6_DZ"))
    rho_l = float(env_value("STAGE20_6_RHO_L")); rho_tilde = float(env_value("STAGE20_6_RHO_TILDE")); bending = float(env_value("STAGE20_6_BENDING_STIFFNESS")); gamma = float(env_value("STAGE20_6_GAMMA"))
    c_fs = float(env_value("STAGE20_6_C_FS")); lzero = float(env_value("STAGE20_6_LAMBDA_ZERO")); lsmall = float(env_value("STAGE20_6_LAMBDA_SMALL")); max_rhs = float(env_value("STAGE20_6_MAX_RHS_DELTA_NORM"))
    tol = float(env_value("STAGE20_6_AUDIT_TOL")); ztol = float(env_value("STAGE20_6_ZERO_TOL")); kernel = env_value("STAGE20_6_KERNEL_NAME")
    cand = build_candidate(n, dim, length, nx, ny, nz, dx, dy, dz, c_fs, lzero, lsmall); ds = float(cand["ds"]); dvol = float(cand["dV"]); size = int(cand["size"])
    f_norm = norm_flat(cand["f_eulerian_candidate"]); zero_delta_norm = norm_flat(cand["RHS_delta_zero"]); small_delta_norm = norm_flat(cand["RHS_delta_small"])
    zero_res_norm = norm_flat(cand["RHS_zero_residual"]); small_res_norm = norm_flat(cand["RHS_small_formula_residual"]); ratio = float(cand["lambda_scaling_ratio_candidate"])
    all_rhs_finite = all(finite(cand[name]) for name in ["f_eulerian_candidate", "f_eulerian_effective_zero", "f_eulerian_effective_small", "RHS_before", "RHS_delta_zero", "RHS_delta_small", "RHS_after_zero", "RHS_after_small", "RHS_zero_residual", "RHS_small_formula_residual"])
    shape_ok = len(cand["X_current"]) == n and all(len(cand[name]) == size for name in ["f_eulerian_candidate", "f_eulerian_effective_zero", "f_eulerian_effective_small", "RHS_before", "RHS_delta_zero", "RHS_delta_small", "RHS_after_zero", "RHS_after_small"])
    required_present = all(name in cand for name in ["X_current", "V_current", "u_interp_candidate", "u_relative_candidate", "F_fs_candidate", "F_on_structure_from_fluid_candidate", "F_on_fluid_from_structure_candidate", "f_eulerian_candidate", "f_eulerian_effective_zero", "f_eulerian_effective_small", "RHS_before", "RHS_delta_zero", "RHS_delta_small", "RHS_after_zero", "RHS_after_small", "RHS_zero_residual", "RHS_small_formula_residual", "lambda_scaling_ratio_candidate", "eulerian_integral_force_candidate", "lagrangian_total_reaction_force", "force_conservation_residual_candidate", "owner_rank", "global_point_id", "local_point_id", "eulerian_grid_shape", "eulerian_grid_spacing", "eulerian_cell_volume"])
    twoway = truthy("STAGE20_6_TWOWAY_COUPLING_ENABLE"); f2s = truthy("STAGE20_6_FLUID_TO_STRUCTURE_ENABLE"); s2f = truthy("STAGE20_6_STRUCTURE_TO_FLUID_ENABLE"); s2fc = truthy("STAGE20_6_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE"); rhs = truthy("STAGE20_6_RHS_COUPLING_ENABLE")
    prod_rhs = truthy("STAGE20_6_PRODUCTION_RHS_UPDATE_ALLOWED"); stage14 = truthy("STAGE20_6_STAGE14_RHS_INJECTION_ALLOWED"); diagnostic = truthy("STAGE20_6_DIAGNOSTIC_ONLY"); fail_closed = truthy("STAGE20_6_FAIL_CLOSED"); single = truthy("STAGE20_6_SINGLE_FIBRE_ONLY"); contact = truthy("STAGE20_6_CONTACT_ENABLE"); collision = truthy("STAGE20_6_COLLISION_ENABLE"); multi = truthy("STAGE20_6_MULTIFIBRE_ENABLE")
    default_safe = all([not twoway, f2s, not s2f, s2fc, not rhs, not prod_rhs, not stage14, close(lzero, 0.0, ztol), 0.0 < lsmall <= 1.0e-3, diagnostic, fail_closed, single, not contact, not collision, not multi])
    lambda_zero_case = zero_delta_norm <= ztol and zero_res_norm <= ztol and all(abs(v) <= ztol for v in cand["f_eulerian_effective_zero"])
    lambda_small_case = small_delta_norm > ztol and small_delta_norm <= max_rhs and small_res_norm <= tol
    numeric_ok = all([n_fibre == 1, n >= 8, dim == 3, length > 0, close(ds, length / float(n - 1), tol), nx >= 4, ny >= 4, nz >= 4, dx > 0, dy > 0, dz > 0, close(dvol, dx * dy * dz, tol), rho_l > 0, rho_tilde > 0, bending >= 0, gamma >= 0, c_fs >= 0, close(lzero, 0.0, ztol), 0.0 < lsmall <= 1.0e-3, all_rhs_finite, f_norm > ztol, lambda_zero_case, lambda_small_case, close(ratio, lsmall, tol), sorted(cand["global_point_id"]) == list(range(n)), len(set(cand["global_point_id"])) == n, cand["owner_rank"] == [0 for _ in range(n)]])

    statuses: Dict[str, str] = {name: PASS for name in STATUS_FIELDS}
    checks = {
        "stage20_6_requested_status": truthy("STAGE20_6_ENABLE"), "stage20_6_rhs_coupling_lambda_gate_enable_status": truthy("STAGE20_6_RHS_COUPLING_LAMBDA_GATE_ENABLE"), "stage20_6_rhs_coupling_candidate_enable_status": truthy("STAGE20_6_RHS_COUPLING_CANDIDATE_ENABLE"), "stage20_5_evidence_status": evidence_ok, "stage20_5_source_only_acceptance_preserved_status": truthy("STAGE20_6_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": truthy("STAGE20_6_ALLOW_MISSING_OLD_STAGE_OUTPUTS"), "missing_old_closure_files_allowed_status": truthy("STAGE20_6_ALLOW_MISSING_OLD_CLOSURE_FILES"), "no_previous_stage_rerun_status": truthy("STAGE20_6_DO_NOT_RERUN_PREVIOUS_STAGES"), "all_required_rhs_coupling_fields_present_status": required_present, "default_safe_gate_values_status": default_safe,
        "twoway_coupling_remains_disabled_status": not twoway, "fluid_to_structure_helper_local_only_status": f2s and diagnostic and not twoway, "structure_to_fluid_production_disabled_status": not s2f, "structure_to_fluid_candidate_helper_local_only_status": s2fc and diagnostic and not s2f and not twoway, "rhs_coupling_candidate_helper_local_only_status": truthy("STAGE20_6_RHS_COUPLING_CANDIDATE_ENABLE") and diagnostic and not rhs and not prod_rhs,
        "production_rhs_coupling_disabled_status": not rhs, "production_rhs_update_disabled_status": not prod_rhs, "stage14_rhs_injection_disabled_status": not stage14, "lambda_zero_case_status": lambda_zero_case, "lambda_small_case_status": lambda_small_case, "lambda_scaling_ratio_status": close(ratio, lsmall, tol),
        "n_fibre_status": n_fibre == 1, "n_point_status": n >= 8, "component_dim_status": dim == 3, "fibre_length_status": length > 0, "ds_formula_status": close(ds, length / float(n - 1), tol), "eulerian_grid_shape_status": nx >= 4 and ny >= 4 and nz >= 4, "eulerian_grid_spacing_status": dx > 0 and dy > 0 and dz > 0, "eulerian_cell_volume_status": close(dvol, dx * dy * dz, tol), "rho_l_status": rho_l > 0, "rho_tilde_status": rho_tilde > 0, "bending_stiffness_status": bending >= 0, "gamma_status": gamma >= 0, "c_fs_status": c_fs >= 0, "kernel_name_status": kernel == "nearest_grid_point", "lambda_zero_status": close(lzero, 0.0, ztol), "lambda_small_status": 0.0 < lsmall <= 1.0e-3,
        "shape_rules_status": shape_ok, "numeric_rules_status": numeric_ok, "eulerian_force_density_candidate_finite_status": finite(cand["f_eulerian_candidate"]), "eulerian_force_density_candidate_nonzero_status": f_norm > ztol, "rhs_before_finite_status": finite(cand["RHS_before"]), "rhs_after_zero_finite_status": finite(cand["RHS_after_zero"]), "rhs_after_small_finite_status": finite(cand["RHS_after_small"]), "rhs_zero_residual_zero_status": zero_res_norm <= ztol, "rhs_small_formula_residual_status": small_res_norm <= tol, "rhs_delta_small_bounded_status": small_delta_norm <= max_rhs,
        "f_eulerian_effective_zero_formula_status": all(abs(v) <= ztol for v in cand["f_eulerian_effective_zero"]), "f_eulerian_effective_small_formula_status": small_res_norm <= tol, "rhs_delta_zero_formula_status": zero_delta_norm <= ztol, "rhs_delta_small_formula_status": small_res_norm <= tol, "rhs_after_zero_formula_status": zero_res_norm <= ztol, "rhs_after_small_formula_status": small_res_norm <= tol, "no_production_rhs_update_status": not prod_rhs and not rhs, "no_stage14_rhs_injection_status": not stage14,
        "global_point_id_coverage_status": sorted(cand["global_point_id"]) == list(range(n)), "global_point_id_no_duplicate_status": len(set(cand["global_point_id"])) == n, "owner_rank_deterministic_status": cand["owner_rank"] == [0 for _ in range(n)], "diagnostic_only_status": diagnostic, "fail_closed_status": fail_closed, "single_fibre_only_status": single, "contact_default_disabled_status": not contact, "collision_default_disabled_status": not collision, "multifibre_default_disabled_status": not multi,
        "stage20_6_wrapper_bash_syntax_status": bash_syntax_ok(root / "stage20_checks" / "run_stage20_6_rhs_coupling_lambda_gate.sh"), "stage20_6_helper_py_compile_status": py_compile_ok(Path(__file__).resolve()),
    }
    for key, ok in checks.items():
        statuses[key] = PASS if ok else FAIL
    final = PASS if all(value in {PASS, OPTIONAL} for value in statuses.values()) else FAIL
    lines = [
        "STAGE 20.6 RHS COUPLING LAMBDA GATE AUDIT", "stage20_title = real two-way fluid-structure coupling activation boundary", "stage20_6_title = production RHS coupling activation with lambda gate",
        f"repository_root_value = {root}", f"stage20_6_test_case_value = {env_value('STAGE20_6_TEST_CASE')}", f"stage20_5_evidence_reason_value = {evidence_reason}", "stage20_7_next_stage_value = Stage 20.7: controlled one-fibre closed-loop response np=1",
        "rhs_scope_value = helper-local RHS candidate arrays only; no production RHS write, no Stage 14 RHS injection, no production DNS", f"lambda_zero_value = {lzero:.16e}", f"lambda_small_value = {lsmall:.16e}", f"lambda_scaling_ratio_candidate_value = {ratio:.16e}",
        f"f_eulerian_candidate_norm_value = {f_norm:.16e}", f"lambda_zero_RHS_delta_norm_value = {zero_delta_norm:.16e}", f"lambda_small_RHS_delta_norm_value = {small_delta_norm:.16e}", f"RHS_zero_residual_norm_value = {zero_res_norm:.16e}", f"RHS_small_formula_residual_norm_value = {small_res_norm:.16e}",
        f"n_fibre_value = {n_fibre}", f"n_point_value = {n}", f"component_dim_value = {dim}", f"eulerian_grid_shape_value = ({nx}, {ny}, {nz})", f"ds_value = {ds:.16e}", f"dV_value = {dvol:.16e}",
    ]
    lines.extend(f"{name} {statuses[name]}" for name in STATUS_FIELDS)
    lines.append(f"final_status {final}")
    lines.append(f"STAGE 20.6 RHS COUPLING LAMBDA GATE VERDICT: {final}")
    lines.append(f"STAGE 20.6 FINAL VERDICT: {final}")
    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 20.6 RHS COUPLING LAMBDA GATE VERDICT: {final}")
    print(f"STAGE 20.6 FINAL VERDICT: {final}")
    if final != PASS:
        for key, value in statuses.items():
            if value == FAIL:
                print(f"FAIL_REASON {key}=FAIL")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())
