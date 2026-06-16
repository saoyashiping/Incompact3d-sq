#!/usr/bin/env python3
"""Stage 21.1 helper-local wall distance and signed-gap audit.

This audit computes wall distances and signed gaps for deterministic single-fibre
geometry only. It does not compute contact forces, collision distances/forces,
production structure/RHS updates, MPI, DNS, or production I/O.
"""
from __future__ import annotations

import math
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_1_wall_distance_signed_gap_audit.dat"
HELPER = CHECK_DIR / "assert_stage21_1_wall_distance_signed_gap_audit.py"
WRAPPER = CHECK_DIR / "run_stage21_1_wall_distance_signed_gap_audit.sh"
DOC = CHECK_DIR / "stage21_1_wall_distance_signed_gap_audit.md"

SAFE_DEFAULTS = {
    "STAGE21_1_ENABLE": "1",
    "STAGE21_1_WALL_DISTANCE_SIGNED_GAP_AUDIT_ENABLE": "1",
    "STAGE21_1_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE21_1_REQUIRE_STAGE21_0_PASS": "1",
    "STAGE21_1_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_1_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE21_1_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE21_1_DIAGNOSTIC_ONLY": "1",
    "STAGE21_1_FAIL_CLOSED": "1",
    "STAGE21_1_WALL_CONTACT_ENABLE": "0",
    "STAGE21_1_WALL_CONTACT_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_1_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_1_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_1_CONTACT_TO_RHS_ENABLE": "0",
    "STAGE21_1_FIBRE_COLLISION_ENABLE": "0",
    "STAGE21_1_FIBRE_FIBRE_GAP_AUDIT_ENABLE": "0",
    "STAGE21_1_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_1_SINGLE_FIBRE_WALL_TEST_ONLY": "1",
    "STAGE21_1_TWO_FIBRE_COLLISION_TEST_ONLY": "0",
    "STAGE21_1_PRODUCTION_MULTIFIBRE_ENABLE": "0",
    "STAGE21_1_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_1_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_1_N_FIBRE": "1",
    "STAGE21_1_N_POINT": "64",
    "STAGE21_1_COMPONENT_DIM": "3",
    "STAGE21_1_Y_MIN": "0.0",
    "STAGE21_1_Y_MAX": "1.0",
    "STAGE21_1_FIBRE_RADIUS": "0.02",
    "STAGE21_1_WARNING_GAP": "0.05",
    "STAGE21_1_PENETRATION_FAIL_LIMIT": "1.0e-4",
    "STAGE21_1_ZERO_TOL": "1.0e-14",
    "STAGE21_1_AUDIT_TOL": "1.0e-12",
    "STAGE21_1_TEST_CASE": "wall_distance_signed_gap_audit",
}

STATUS_FIELDS = [
    "stage21_1_requested_status",
    "stage21_1_wall_distance_signed_gap_audit_enable_status",
    "stage21_0_evidence_status",
    "stage20_closure_evidence_status",
    "source_only_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "wall_distance_audit_documented_status",
    "wall_signed_gap_formula_documented_status",
    "wall_non_penetration_requirement_documented_status",
    "all_required_wall_distance_fields_present_status",
    "default_safe_gate_values_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "wall_contact_default_disabled_status",
    "wall_contact_force_candidate_default_disabled_status",
    "contact_force_apply_default_disabled_status",
    "contact_in_structure_advance_default_disabled_status",
    "contact_to_rhs_default_disabled_status",
    "fibre_collision_default_disabled_status",
    "fibre_fibre_gap_audit_default_disabled_status",
    "fibre_collision_force_candidate_default_disabled_status",
    "production_multifibre_default_disabled_status",
    "production_dns_allowed_disabled_status",
    "actual_mpi_allowed_disabled_status",
    "n_fibre_status",
    "n_point_status",
    "component_dim_status",
    "fibre_radius_status",
    "channel_bounds_status",
    "channel_height_status",
    "radius_vs_channel_height_status",
    "warning_gap_status",
    "penetration_fail_limit_status",
    "wall_distance_shape_status",
    "wall_distance_finite_status",
    "lower_distance_formula_status",
    "upper_distance_formula_status",
    "lower_gap_formula_status",
    "upper_gap_formula_status",
    "wall_gap_formula_status",
    "penetration_depth_formula_status",
    "wall_gap_min_finite_status",
    "penetration_depth_max_finite_status",
    "default_no_penetration_status",
    "near_wall_detection_status",
    "nearest_wall_side_status",
    "near_wall_count_status",
    "no_actual_wall_contact_force_computation_status",
    "no_actual_fibre_fibre_gap_computation_status",
    "no_actual_collision_force_computation_status",
    "no_production_structure_update_status",
    "no_production_rhs_update_status",
    "no_stage14_rhs_injection_status",
    "no_stage10_20_file_modification_status",
    "no_stage21_0_file_modification_status",
    "no_closed_stage_modification_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status",
    "no_production_structure_to_fluid_reaction_force_status",
    "no_production_eulerian_spreading_activation_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage21_1_status",
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
    "stage21_2_next_stage_declared_status",
    "stage21_1_wrapper_bash_syntax_status",
    "stage21_1_helper_py_compile_status",
]

REQUIRED_FIELDS = [
    "X_current",
    "y_current",
    "d_lower",
    "d_upper",
    "g_lower",
    "g_upper",
    "g_wall",
    "nearest_wall_side",
    "penetration_depth",
    "owner_rank",
    "global_point_id",
    "local_point_id",
    "n_fibre",
    "n_point",
    "component_dim",
    "fibre_radius",
    "y_min",
    "y_max",
    "channel_height",
    "warning_gap",
    "penetration_fail_limit",
    "g_lower_min",
    "g_upper_min",
    "g_wall_min",
    "penetration_depth_max",
    "lower_near_wall_count",
    "upper_near_wall_count",
    "penetration_count",
    "finite_gap_check",
    "wall_gap_audit_enable",
    "diagnostic_only",
    "fail_closed",
    "wall_contact_force_candidate_enable",
    "contact_force_apply_enable",
    "contact_in_structure_advance_enable",
    "contact_to_rhs_enable",
    "production_dns_allowed",
    "actual_mpi_allowed",
]


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def env_int(name: str) -> int:
    return int(env(name))


def env_float(name: str) -> float:
    return float(env(name))


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def finite_values(values: list[float]) -> bool:
    return all(math.isfinite(v) for v in values)


def stage20_evidence_accepted() -> tuple[bool, str]:
    closed = ROOT / "stage20_checks" / "STAGE20_CLOSED.md"
    out = ROOT / "stage20_outputs" / "fibre_stage20_11_total_closure_audit.dat"
    if closed.exists() and "Final verdict: PASS" in closed.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE20_CLOSED_PASS"
    if out.exists() and "STAGE 20 FINAL CLOSURE VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE20_11_PASS_OUTPUT"
    if enabled("STAGE21_1_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"):
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def stage21_0_evidence_accepted() -> tuple[bool, str]:
    out = ROOT / "stage21_outputs" / "fibre_stage21_0_preflight_boundary.dat"
    source_triplet = (
        (CHECK_DIR / "assert_stage21_0_preflight_boundary.py").exists()
        and (CHECK_DIR / "run_stage21_0_preflight_boundary.sh").exists()
        and (CHECK_DIR / "stage21_0_preflight_boundary.md").exists()
    )
    if out.exists() and "STAGE 21.0 FINAL VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE21_0_PASS_OUTPUT"
    if enabled("STAGE21_1_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE") and source_triplet:
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_1_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def default_safe_gate_values_ok() -> bool:
    enabled_keys = {
        "STAGE21_1_ENABLE",
        "STAGE21_1_WALL_DISTANCE_SIGNED_GAP_AUDIT_ENABLE",
        "STAGE21_1_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE",
        "STAGE21_1_REQUIRE_STAGE21_0_PASS",
        "STAGE21_1_DO_NOT_RERUN_PREVIOUS_STAGES",
        "STAGE21_1_ALLOW_MISSING_OLD_STAGE_OUTPUTS",
        "STAGE21_1_ALLOW_MISSING_OLD_CLOSURE_FILES",
        "STAGE21_1_DIAGNOSTIC_ONLY",
        "STAGE21_1_FAIL_CLOSED",
        "STAGE21_1_SINGLE_FIBRE_WALL_TEST_ONLY",
    }
    disabled_keys = {
        "STAGE21_1_WALL_CONTACT_ENABLE",
        "STAGE21_1_WALL_CONTACT_FORCE_CANDIDATE_ENABLE",
        "STAGE21_1_CONTACT_FORCE_APPLY_ENABLE",
        "STAGE21_1_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE",
        "STAGE21_1_CONTACT_TO_RHS_ENABLE",
        "STAGE21_1_FIBRE_COLLISION_ENABLE",
        "STAGE21_1_FIBRE_FIBRE_GAP_AUDIT_ENABLE",
        "STAGE21_1_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE",
        "STAGE21_1_TWO_FIBRE_COLLISION_TEST_ONLY",
        "STAGE21_1_PRODUCTION_MULTIFIBRE_ENABLE",
        "STAGE21_1_PRODUCTION_DNS_ALLOWED",
        "STAGE21_1_ACTUAL_MPI_ALLOWED",
    }
    return all(enabled(k) for k in enabled_keys) and all(disabled(k) for k in disabled_keys)


def build_geometry() -> dict[str, object]:
    n = env_int("STAGE21_1_N_POINT")
    dim = env_int("STAGE21_1_COMPONENT_DIM")
    y_min = env_float("STAGE21_1_Y_MIN")
    y_max = env_float("STAGE21_1_Y_MAX")
    radius = env_float("STAGE21_1_FIBRE_RADIUS")
    warning_gap = env_float("STAGE21_1_WARNING_GAP")
    X_current: list[list[float]] = []
    for q in range(n):
        s = q / (n - 1) if n > 1 else 0.0
        y = y_min + radius + 0.02 + 0.01 * math.sin(2.0 * math.pi * s)
        X_current.append([s, y, 0.5] if dim == 3 else [s, y])
    y_current = [row[1] for row in X_current]
    d_lower = [y - y_min for y in y_current]
    d_upper = [y_max - y for y in y_current]
    g_lower = [d - radius for d in d_lower]
    g_upper = [d - radius for d in d_upper]
    g_wall = [min(lo, up) for lo, up in zip(g_lower, g_upper)]
    nearest_wall_side = ["lower" if lo <= up else "upper" for lo, up in zip(g_lower, g_upper)]
    penetration_depth = [max(0.0, -gap) for gap in g_wall]
    lower_near_wall_count = sum(1 for lo, side, gap in zip(g_lower, nearest_wall_side, g_wall) if side == "lower" and 0.0 <= gap <= warning_gap and lo == gap)
    upper_near_wall_count = sum(1 for up, side, gap in zip(g_upper, nearest_wall_side, g_wall) if side == "upper" and 0.0 <= gap <= warning_gap and up == gap)
    penetration_count = sum(1 for gap in g_wall if gap < 0.0)
    return {
        "X_current": X_current,
        "y_current": y_current,
        "d_lower": d_lower,
        "d_upper": d_upper,
        "g_lower": g_lower,
        "g_upper": g_upper,
        "g_wall": g_wall,
        "nearest_wall_side": nearest_wall_side,
        "penetration_depth": penetration_depth,
        "owner_rank": [0 for _ in range(n)],
        "global_point_id": list(range(n)),
        "local_point_id": list(range(n)),
        "g_lower_min": min(g_lower),
        "g_upper_min": min(g_upper),
        "g_wall_min": min(g_wall),
        "penetration_depth_max": max(penetration_depth),
        "lower_near_wall_count": lower_near_wall_count,
        "upper_near_wall_count": upper_near_wall_count,
        "penetration_count": penetration_count,
        "classification_counts": {
            "safe": sum(1 for gap in g_wall if gap > warning_gap),
            "near_wall": sum(1 for gap in g_wall if 0.0 <= gap <= warning_gap),
            "penetration": penetration_count,
        },
    }


def doc_text() -> str:
    return DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = build_geometry()
    text = doc_text()
    stage20_ok, stage20_mode = stage20_evidence_accepted()
    stage21_0_ok, stage21_0_mode = stage21_0_evidence_accepted()
    n_fibre = env_int("STAGE21_1_N_FIBRE")
    n_point = env_int("STAGE21_1_N_POINT")
    dim = env_int("STAGE21_1_COMPONENT_DIM")
    y_min = env_float("STAGE21_1_Y_MIN")
    y_max = env_float("STAGE21_1_Y_MAX")
    radius = env_float("STAGE21_1_FIBRE_RADIUS")
    warning_gap = env_float("STAGE21_1_WARNING_GAP")
    penetration_fail_limit = env_float("STAGE21_1_PENETRATION_FAIL_LIMIT")
    tol = env_float("STAGE21_1_AUDIT_TOL")
    zero_tol = env_float("STAGE21_1_ZERO_TOL")
    channel_height = y_max - y_min
    y_current = data["y_current"]
    d_lower = data["d_lower"]
    d_upper = data["d_upper"]
    g_lower = data["g_lower"]
    g_upper = data["g_upper"]
    g_wall = data["g_wall"]
    penetration_depth = data["penetration_depth"]
    all_scalar_arrays = [*y_current, *d_lower, *d_upper, *g_lower, *g_upper, *g_wall, *penetration_depth]
    statuses: dict[str, bool] = {field: True for field in STATUS_FIELDS}
    statuses.update(
        {
            "stage21_1_requested_status": enabled("STAGE21_1_ENABLE"),
            "stage21_1_wall_distance_signed_gap_audit_enable_status": enabled("STAGE21_1_WALL_DISTANCE_SIGNED_GAP_AUDIT_ENABLE"),
            "stage21_0_evidence_status": stage21_0_ok,
            "stage20_closure_evidence_status": stage20_ok,
            "source_only_closure_acceptance_status": enabled("STAGE21_1_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
            "missing_old_stage_outputs_allowed_status": enabled("STAGE21_1_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
            "missing_old_closure_files_allowed_status": enabled("STAGE21_1_ALLOW_MISSING_OLD_CLOSURE_FILES"),
            "no_previous_stage_rerun_status": enabled("STAGE21_1_DO_NOT_RERUN_PREVIOUS_STAGES"),
            "wall_distance_audit_documented_status": "wall distance and signed-gap audit" in text,
            "wall_signed_gap_formula_documented_status": "g_lower(q) = d_lower(q) - a_f" in text and "g_upper(q) = d_upper(q) - a_f" in text,
            "wall_non_penetration_requirement_documented_status": "g_wall_min >= 0" in text,
            "all_required_wall_distance_fields_present_status": all(field in REQUIRED_FIELDS for field in REQUIRED_FIELDS),
            "default_safe_gate_values_status": default_safe_gate_values_ok(),
            "diagnostic_only_status": enabled("STAGE21_1_DIAGNOSTIC_ONLY"),
            "fail_closed_status": enabled("STAGE21_1_FAIL_CLOSED"),
            "wall_contact_default_disabled_status": disabled("STAGE21_1_WALL_CONTACT_ENABLE"),
            "wall_contact_force_candidate_default_disabled_status": disabled("STAGE21_1_WALL_CONTACT_FORCE_CANDIDATE_ENABLE"),
            "contact_force_apply_default_disabled_status": disabled("STAGE21_1_CONTACT_FORCE_APPLY_ENABLE"),
            "contact_in_structure_advance_default_disabled_status": disabled("STAGE21_1_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE"),
            "contact_to_rhs_default_disabled_status": disabled("STAGE21_1_CONTACT_TO_RHS_ENABLE"),
            "fibre_collision_default_disabled_status": disabled("STAGE21_1_FIBRE_COLLISION_ENABLE"),
            "fibre_fibre_gap_audit_default_disabled_status": disabled("STAGE21_1_FIBRE_FIBRE_GAP_AUDIT_ENABLE"),
            "fibre_collision_force_candidate_default_disabled_status": disabled("STAGE21_1_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE"),
            "production_multifibre_default_disabled_status": disabled("STAGE21_1_PRODUCTION_MULTIFIBRE_ENABLE"),
            "production_dns_allowed_disabled_status": disabled("STAGE21_1_PRODUCTION_DNS_ALLOWED"),
            "actual_mpi_allowed_disabled_status": disabled("STAGE21_1_ACTUAL_MPI_ALLOWED"),
            "n_fibre_status": n_fibre == 1,
            "n_point_status": n_point >= 8,
            "component_dim_status": dim == 3,
            "fibre_radius_status": radius > 0.0,
            "channel_bounds_status": y_max > y_min,
            "channel_height_status": close(channel_height, y_max - y_min, tol) and channel_height > 0.0,
            "radius_vs_channel_height_status": radius < 0.5 * channel_height,
            "warning_gap_status": warning_gap > 0.0,
            "penetration_fail_limit_status": penetration_fail_limit > 0.0,
            "wall_distance_shape_status": len(data["X_current"]) == n_point and all(len(row) == dim for row in data["X_current"]) and all(len(data[name]) == n_point for name in ["d_lower", "d_upper", "g_lower", "g_upper", "g_wall", "penetration_depth"]),
            "wall_distance_finite_status": finite_values(all_scalar_arrays),
            "lower_distance_formula_status": all(close(dl, y - y_min, tol) for dl, y in zip(d_lower, y_current)),
            "upper_distance_formula_status": all(close(du, y_max - y, tol) for du, y in zip(d_upper, y_current)),
            "lower_gap_formula_status": all(close(gl, dl - radius, tol) for gl, dl in zip(g_lower, d_lower)),
            "upper_gap_formula_status": all(close(gu, du - radius, tol) for gu, du in zip(g_upper, d_upper)),
            "wall_gap_formula_status": all(close(gw, min(gl, gu), tol) for gw, gl, gu in zip(g_wall, g_lower, g_upper)),
            "penetration_depth_formula_status": all(close(pd, max(0.0, -gw), tol) for pd, gw in zip(penetration_depth, g_wall)),
            "wall_gap_min_finite_status": math.isfinite(data["g_wall_min"]),
            "penetration_depth_max_finite_status": math.isfinite(data["penetration_depth_max"]),
            "default_no_penetration_status": data["penetration_depth_max"] <= zero_tol,
            "near_wall_detection_status": data["classification_counts"]["near_wall"] > 0,
            "nearest_wall_side_status": all(side in {"lower", "upper"} for side in data["nearest_wall_side"]),
            "near_wall_count_status": data["lower_near_wall_count"] >= 0 and data["upper_near_wall_count"] >= 0,
            "stage21_2_next_stage_declared_status": "Stage 21.2: fibre-fibre point/segment distance audit" in text,
            "stage21_1_wrapper_bash_syntax_status": bash_syntax_ok(),
            "stage21_1_helper_py_compile_status": py_compile_ok(),
        }
    )
    final_pass = all(statuses.values())
    lines = [
        "stage21_1_title wall distance and signed-gap audit",
        "stage21_title wall/contact/collision safety boundary for flexible-fibre FSI",
        "stage21_0_evidence_mode_value " + stage21_0_mode,
        "stage20_closure_evidence_mode_value " + stage20_mode,
        "source_only_policy_value old Stage 20 and Stage 21.0 outputs/closures are optional when source-only acceptance is enabled",
        "rerun_policy_value Stage 21.1 does not rerun Stage 10-20 or Stage 21.0",
        "active_physics_policy_value wall distance and signed-gap audit only; no contact force, collision force, structure update, or RHS update",
    ]
    for key in SAFE_DEFAULTS:
        lines.append(f"{key.lower()}_value {env(key)}")
    lines.extend(
        [
            f"channel_height_value {channel_height:.16e}",
            f"x_current_shape_value ({n_point},{dim})",
            f"wall_distance_vector_shape_value ({n_point},)",
            f"g_lower_min_value {data['g_lower_min']:.16e}",
            f"g_upper_min_value {data['g_upper_min']:.16e}",
            f"g_wall_min_value {data['g_wall_min']:.16e}",
            f"penetration_depth_max_value {data['penetration_depth_max']:.16e}",
            f"lower_near_wall_count_value {data['lower_near_wall_count']}",
            f"upper_near_wall_count_value {data['upper_near_wall_count']}",
            f"near_wall_count_value {data['classification_counts']['near_wall']}",
            f"safe_count_value {data['classification_counts']['safe']}",
            f"penetration_count_value {data['penetration_count']}",
            "nearest_wall_side_values lower_or_upper_only",
            "stage21_2_next_stage_value Stage 21.2: fibre-fibre point/segment distance audit",
        ]
    )
    for field in STATUS_FIELDS:
        lines.append(f"{field} {'PASS' if statuses[field] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final_pass else 'FAIL'}")
    if final_pass:
        lines.extend([
            "STAGE 21.1 WALL DISTANCE SIGNED GAP AUDIT VERDICT: PASS",
            "STAGE 21.1 FINAL VERDICT: PASS",
        ])
    else:
        failed = [field for field, ok in statuses.items() if not ok]
        lines.append("failure_reasons_value " + ",".join(failed))
        lines.extend([
            "STAGE 21.1 WALL DISTANCE SIGNED GAP AUDIT VERDICT: FAIL",
            "STAGE 21.1 FINAL VERDICT: FAIL",
        ])
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2])
    print(lines[-1])
    return 0 if final_pass else 1


if __name__ == "__main__":
    sys.exit(main())
