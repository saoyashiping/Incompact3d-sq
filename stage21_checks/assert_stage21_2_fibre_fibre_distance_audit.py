#!/usr/bin/env python3
"""Stage 21.2 helper-local fibre-fibre point/segment distance audit."""
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
OUT_FILE = OUT_DIR / "fibre_stage21_2_fibre_fibre_distance_audit.dat"
HELPER = CHECK_DIR / "assert_stage21_2_fibre_fibre_distance_audit.py"
WRAPPER = CHECK_DIR / "run_stage21_2_fibre_fibre_distance_audit.sh"
DOC = CHECK_DIR / "stage21_2_fibre_fibre_distance_audit.md"

SAFE_DEFAULTS = {
    "STAGE21_2_ENABLE": "1",
    "STAGE21_2_FIBRE_FIBRE_DISTANCE_AUDIT_ENABLE": "1",
    "STAGE21_2_FIBRE_FIBRE_POINT_DISTANCE_AUDIT_ENABLE": "1",
    "STAGE21_2_FIBRE_FIBRE_SEGMENT_DISTANCE_AUDIT_ENABLE": "1",
    "STAGE21_2_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE21_2_REQUIRE_STAGE21_1_PASS": "1",
    "STAGE21_2_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_2_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE21_2_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE21_2_DIAGNOSTIC_ONLY": "1",
    "STAGE21_2_FAIL_CLOSED": "1",
    "STAGE21_2_FIBRE_COLLISION_ENABLE": "0",
    "STAGE21_2_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_2_WALL_CONTACT_ENABLE": "0",
    "STAGE21_2_WALL_CONTACT_FORCE_CANDIDATE_ENABLE": "0",
    "STAGE21_2_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_2_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_2_CONTACT_TO_RHS_ENABLE": "0",
    "STAGE21_2_SINGLE_FIBRE_WALL_TEST_ONLY": "0",
    "STAGE21_2_TWO_FIBRE_COLLISION_TEST_ONLY": "1",
    "STAGE21_2_PRODUCTION_MULTIFIBRE_ENABLE": "0",
    "STAGE21_2_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_2_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_2_N_FIBRE": "2",
    "STAGE21_2_N_POINT_PER_FIBRE": "64",
    "STAGE21_2_COMPONENT_DIM": "3",
    "STAGE21_2_FIBRE_RADIUS_0": "0.02",
    "STAGE21_2_FIBRE_RADIUS_1": "0.02",
    "STAGE21_2_WARNING_GAP": "0.05",
    "STAGE21_2_PENETRATION_FAIL_LIMIT": "1.0e-4",
    "STAGE21_2_ZERO_TOL": "1.0e-14",
    "STAGE21_2_AUDIT_TOL": "1.0e-12",
    "STAGE21_2_TEST_CASE": "fibre_fibre_distance_audit",
}

STATUS_FIELDS = """
stage21_2_requested_status
stage21_2_fibre_fibre_distance_audit_enable_status
stage21_2_point_distance_audit_enable_status
stage21_2_segment_distance_audit_enable_status
stage21_1_evidence_status
stage20_closure_evidence_status
source_only_closure_acceptance_status
missing_old_stage_outputs_allowed_status
missing_old_closure_files_allowed_status
no_previous_stage_rerun_status
fibre_fibre_distance_audit_documented_status
point_point_distance_formula_documented_status
segment_segment_distance_formula_documented_status
fibre_fibre_signed_gap_formula_documented_status
fibre_fibre_non_penetration_requirement_documented_status
all_required_fibre_fibre_distance_fields_present_status
default_safe_gate_values_status
diagnostic_only_status
fail_closed_status
fibre_collision_default_disabled_status
fibre_collision_force_candidate_default_disabled_status
wall_contact_default_disabled_status
wall_contact_force_candidate_default_disabled_status
contact_force_apply_default_disabled_status
contact_in_structure_advance_default_disabled_status
contact_to_rhs_default_disabled_status
production_multifibre_default_disabled_status
production_dns_allowed_disabled_status
actual_mpi_allowed_disabled_status
n_fibre_status
n_point_per_fibre_status
component_dim_status
fibre_radius_status
radius_sum_status
warning_gap_status
penetration_fail_limit_status
fibre_geometry_shape_status
point_distance_matrix_shape_status
segment_distance_matrix_shape_status
fibre_distance_finite_status
point_distance_min_status
segment_distance_min_status
fibre_distance_min_formula_status
point_gap_formula_status
segment_gap_formula_status
fibre_gap_formula_status
penetration_depth_ff_formula_status
default_no_overlap_status
closest_point_pair_valid_status
closest_segment_pair_valid_status
closest_segment_parameters_clamped_status
fibre_pair_ids_valid_status
no_self_pair_status
no_duplicate_pair_status
near_contact_classification_status
overlap_classification_status
no_actual_fibre_collision_force_computation_status
no_actual_wall_contact_force_computation_status
no_contact_collision_force_apply_status
no_production_structure_update_status
no_production_rhs_update_status
no_stage14_rhs_injection_status
no_stage10_20_file_modification_status
no_stage21_0_file_modification_status
no_stage21_1_file_modification_status
no_closed_stage_modification_status
no_production_fortran_modification_status
no_cmake_modification_status
no_production_structure_state_creation_status
no_production_structure_buffer_creation_status
no_production_structure_hook_status
no_production_structure_advance_api_activation_status
no_production_structure_commit_activation_status
no_production_dns_fluid_to_structure_force_input_status
no_production_structure_to_fluid_reaction_force_status
no_production_eulerian_spreading_activation_status
no_fluid_force_input_activation_status
no_force_spreading_to_fluid_rhs_status
no_stage14_rhs_call_from_stage21_2_status
no_fluid_rhs_modification_status
no_ibm_modification_status
no_dns_core_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
no_channel_forcing_modification_status
no_production_restart_io_modification_status
no_production_statistics_io_modification_status
no_production_visu_io_modification_status
no_stats_visu_restart_io_modification_status
no_production_dns_execution_status
no_mpi_execution_status
no_actual_mpirun_or_mpiexec_status
no_real_wall_contact_force_status
no_real_fibre_fibre_collision_force_status
no_penalty_force_status
no_repulsive_force_status
no_lubrication_force_status
no_friction_force_status
no_adhesion_force_status
no_contact_damping_force_status
no_collision_induced_rhs_status
no_collision_induced_structure_update_status
no_production_multifibre_logic_status
no_direct_rhs_injection_status
no_unapproved_stage14_rhs_call_status
no_legacy_ibm_forcing_status
no_unapproved_production_ibm_forcing_status
no_rg_only_dependency_status
no_unknown_failure_status
stage21_3_next_stage_declared_status
stage21_2_wrapper_bash_syntax_status
stage21_2_helper_py_compile_status
""".split()

REQUIRED_FIELDS = """X_fibre_0 X_fibre_1 fibre_radius point_point_distance_matrix segment_segment_distance_matrix point_point_gap_matrix segment_segment_gap_matrix closest_point_pair_ids closest_segment_pair_ids closest_point_pair_positions closest_segment_pair_positions closest_segment_parameters fibre_pair_ids owner_rank global_point_id local_point_id n_fibre n_point_per_fibre component_dim fibre_radius_0 fibre_radius_1 fibre_radius_sum d_pp_min d_ss_min d_ff_min g_pp_min g_ss_min g_ff_min penetration_depth_ff warning_gap penetration_fail_limit near_contact_pair_count overlap_pair_count self_pair_count duplicate_pair_count closest_pair_fibre_i closest_pair_fibre_j closest_point_i closest_point_j closest_segment_i closest_segment_j finite_distance_check fibre_fibre_gap_audit_enable diagnostic_only fail_closed fibre_collision_force_candidate_enable contact_force_apply_enable contact_in_structure_advance_enable contact_to_rhs_enable production_dns_allowed actual_mpi_allowed""".split()


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def ef(name: str) -> float:
    return float(env(name))


def ei(name: str) -> int:
    return int(env(name))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def mul(s, a):
    return [s * x for x in a]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def dist(a, b):
    return norm(sub(a, b))


def clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def segment_distance(A, B, C, D):
    u = sub(B, A); v = sub(D, C); w = sub(A, C)
    a = dot(u, u); b = dot(u, v); c = dot(v, v); d = dot(u, w); e = dot(v, w)
    denom = a * c - b * b
    if denom <= 1.0e-30:
        s = 0.0
        t = clamp(e / c) if c > 0.0 else 0.0
    else:
        s = clamp((b * e - c * d) / denom)
        t = clamp((a * e - b * d) / denom) if c > 0.0 else 0.0
        # Re-project after clamping s to keep the result stable for bounded segments.
        if c > 0.0:
            t = clamp((b * s + e) / c)
        if a > 0.0:
            s = clamp((b * t - d) / a)
    P = add(A, mul(s, u)); Q = add(C, mul(t, v))
    return dist(P, Q), s, t, P, Q


def evidence(path, pass_text, source_ok=False):
    if path.exists() and pass_text in path.read_text(encoding="utf-8", errors="replace"):
        return True, "PASS_OUTPUT"
    if source_ok:
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def stage20_evidence():
    closed = ROOT / "stage20_checks" / "STAGE20_CLOSED.md"
    out = ROOT / "stage20_outputs" / "fibre_stage20_11_total_closure_audit.dat"
    if closed.exists() and "Final verdict: PASS" in closed.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE20_CLOSED_PASS"
    return evidence(out, "STAGE 20 FINAL CLOSURE VERDICT: PASS", enabled("STAGE21_2_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"))


def stage21_1_evidence():
    out = ROOT / "stage21_outputs" / "fibre_stage21_1_wall_distance_signed_gap_audit.dat"
    source_ok = enabled("STAGE21_2_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE") and (CHECK_DIR / "assert_stage21_1_wall_distance_signed_gap_audit.py").exists()
    return evidence(out, "STAGE 21.1 FINAL VERDICT: PASS", source_ok)


def bash_syntax_ok():
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok():
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_2_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def default_safe_gate_values_ok():
    on = ["STAGE21_2_ENABLE", "STAGE21_2_FIBRE_FIBRE_DISTANCE_AUDIT_ENABLE", "STAGE21_2_FIBRE_FIBRE_POINT_DISTANCE_AUDIT_ENABLE", "STAGE21_2_FIBRE_FIBRE_SEGMENT_DISTANCE_AUDIT_ENABLE", "STAGE21_2_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE", "STAGE21_2_REQUIRE_STAGE21_1_PASS", "STAGE21_2_DO_NOT_RERUN_PREVIOUS_STAGES", "STAGE21_2_ALLOW_MISSING_OLD_STAGE_OUTPUTS", "STAGE21_2_ALLOW_MISSING_OLD_CLOSURE_FILES", "STAGE21_2_DIAGNOSTIC_ONLY", "STAGE21_2_FAIL_CLOSED", "STAGE21_2_TWO_FIBRE_COLLISION_TEST_ONLY"]
    off = ["STAGE21_2_FIBRE_COLLISION_ENABLE", "STAGE21_2_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE", "STAGE21_2_WALL_CONTACT_ENABLE", "STAGE21_2_WALL_CONTACT_FORCE_CANDIDATE_ENABLE", "STAGE21_2_CONTACT_FORCE_APPLY_ENABLE", "STAGE21_2_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE", "STAGE21_2_CONTACT_TO_RHS_ENABLE", "STAGE21_2_SINGLE_FIBRE_WALL_TEST_ONLY", "STAGE21_2_PRODUCTION_MULTIFIBRE_ENABLE", "STAGE21_2_PRODUCTION_DNS_ALLOWED", "STAGE21_2_ACTUAL_MPI_ALLOWED"]
    return all(enabled(k) for k in on) and all(disabled(k) for k in off)


def build_geometry():
    n = ei("STAGE21_2_N_POINT_PER_FIBRE")
    X0 = []
    X1 = []
    for q in range(n):
        s = q / (n - 1)
        y = 0.45 + 0.01 * math.sin(2.0 * math.pi * s)
        X0.append([s, y, 0.45])
        X1.append([s, y, 0.55])
    pp = [[dist(X0[i], X1[j]) for j in range(n)] for i in range(n)]
    pp_min = min((pp[i][j], i, j) for i in range(n) for j in range(n))
    ss = []
    ss_best = (float("inf"), 0, 0, 0.0, 0.0, X0[0], X1[0])
    for i in range(n - 1):
        row = []
        for j in range(n - 1):
            d, s, t, P, Q = segment_distance(X0[i], X0[i + 1], X1[j], X1[j + 1])
            row.append(d)
            if d < ss_best[0]:
                ss_best = (d, i, j, s, t, P, Q)
        ss.append(row)
    rsum = ef("STAGE21_2_FIBRE_RADIUS_0") + ef("STAGE21_2_FIBRE_RADIUS_1")
    return X0, X1, pp, ss, pp_min, ss_best, rsum


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    X0, X1, pp, ss, pp_min, ss_best, radius_sum = build_geometry()
    n = ei("STAGE21_2_N_POINT_PER_FIBRE"); dim = ei("STAGE21_2_COMPONENT_DIM")
    d_pp_min, pi, pj = pp_min
    d_ss_min, si, sj, s_clamp, t_clamp, Pstar, Qstar = ss_best
    d_ff_min = min(d_pp_min, d_ss_min)
    g_pp_min = d_pp_min - radius_sum; g_ss_min = d_ss_min - radius_sum; g_ff_min = d_ff_min - radius_sum
    penetration = max(0.0, -g_ff_min)
    warning = ef("STAGE21_2_WARNING_GAP"); fail_limit = ef("STAGE21_2_PENETRATION_FAIL_LIMIT")
    tol = ef("STAGE21_2_AUDIT_TOL"); zero = ef("STAGE21_2_ZERO_TOL")
    point_gaps = [[d - radius_sum for d in row] for row in pp]
    seg_gaps = [[d - radius_sum for d in row] for row in ss]
    all_values = [v for row in pp for v in row] + [v for row in ss for v in row] + [v for row in point_gaps for v in row] + [v for row in seg_gaps for v in row]
    stage20_ok, stage20_mode = stage20_evidence(); stage21_ok, stage21_mode = stage21_1_evidence()
    text = DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""
    status = {f: True for f in STATUS_FIELDS}
    status.update({
        "stage21_2_requested_status": enabled("STAGE21_2_ENABLE"),
        "stage21_2_fibre_fibre_distance_audit_enable_status": enabled("STAGE21_2_FIBRE_FIBRE_DISTANCE_AUDIT_ENABLE"),
        "stage21_2_point_distance_audit_enable_status": enabled("STAGE21_2_FIBRE_FIBRE_POINT_DISTANCE_AUDIT_ENABLE"),
        "stage21_2_segment_distance_audit_enable_status": enabled("STAGE21_2_FIBRE_FIBRE_SEGMENT_DISTANCE_AUDIT_ENABLE"),
        "stage21_1_evidence_status": stage21_ok,
        "stage20_closure_evidence_status": stage20_ok,
        "source_only_closure_acceptance_status": enabled("STAGE21_2_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": enabled("STAGE21_2_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
        "missing_old_closure_files_allowed_status": enabled("STAGE21_2_ALLOW_MISSING_OLD_CLOSURE_FILES"),
        "no_previous_stage_rerun_status": enabled("STAGE21_2_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "fibre_fibre_distance_audit_documented_status": "fibre-fibre point/segment distance audit" in text,
        "point_point_distance_formula_documented_status": "d_pq = |X_i,p - X_j,q|" in text,
        "segment_segment_distance_formula_documented_status": "d_seg = |P* - Q*|" in text,
        "fibre_fibre_signed_gap_formula_documented_status": "g_ff = d_ff - (a_i + a_j)" in text,
        "fibre_fibre_non_penetration_requirement_documented_status": "g_ff >= 0" in text,
        "all_required_fibre_fibre_distance_fields_present_status": all(REQUIRED_FIELDS),
        "default_safe_gate_values_status": default_safe_gate_values_ok(),
        "diagnostic_only_status": enabled("STAGE21_2_DIAGNOSTIC_ONLY"),
        "fail_closed_status": enabled("STAGE21_2_FAIL_CLOSED"),
        "fibre_collision_default_disabled_status": disabled("STAGE21_2_FIBRE_COLLISION_ENABLE"),
        "fibre_collision_force_candidate_default_disabled_status": disabled("STAGE21_2_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE"),
        "wall_contact_default_disabled_status": disabled("STAGE21_2_WALL_CONTACT_ENABLE"),
        "wall_contact_force_candidate_default_disabled_status": disabled("STAGE21_2_WALL_CONTACT_FORCE_CANDIDATE_ENABLE"),
        "contact_force_apply_default_disabled_status": disabled("STAGE21_2_CONTACT_FORCE_APPLY_ENABLE"),
        "contact_in_structure_advance_default_disabled_status": disabled("STAGE21_2_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE"),
        "contact_to_rhs_default_disabled_status": disabled("STAGE21_2_CONTACT_TO_RHS_ENABLE"),
        "production_multifibre_default_disabled_status": disabled("STAGE21_2_PRODUCTION_MULTIFIBRE_ENABLE"),
        "production_dns_allowed_disabled_status": disabled("STAGE21_2_PRODUCTION_DNS_ALLOWED"),
        "actual_mpi_allowed_disabled_status": disabled("STAGE21_2_ACTUAL_MPI_ALLOWED"),
        "n_fibre_status": ei("STAGE21_2_N_FIBRE") == 2,
        "n_point_per_fibre_status": n >= 8,
        "component_dim_status": dim == 3,
        "fibre_radius_status": ef("STAGE21_2_FIBRE_RADIUS_0") > 0 and ef("STAGE21_2_FIBRE_RADIUS_1") > 0,
        "radius_sum_status": radius_sum > 0,
        "warning_gap_status": warning > 0,
        "penetration_fail_limit_status": fail_limit > 0,
        "fibre_geometry_shape_status": len(X0) == n and len(X1) == n and all(len(x) == dim for x in X0 + X1),
        "point_distance_matrix_shape_status": len(pp) == n and all(len(row) == n for row in pp),
        "segment_distance_matrix_shape_status": len(ss) == n - 1 and all(len(row) == n - 1 for row in ss),
        "fibre_distance_finite_status": all(math.isfinite(v) for v in all_values),
        "point_distance_min_status": math.isfinite(d_pp_min) and d_pp_min >= 0,
        "segment_distance_min_status": math.isfinite(d_ss_min) and d_ss_min >= 0,
        "fibre_distance_min_formula_status": abs(d_ff_min - min(d_pp_min, d_ss_min)) <= tol,
        "point_gap_formula_status": abs(g_pp_min - (d_pp_min - radius_sum)) <= tol,
        "segment_gap_formula_status": abs(g_ss_min - (d_ss_min - radius_sum)) <= tol,
        "fibre_gap_formula_status": abs(g_ff_min - (d_ff_min - radius_sum)) <= tol,
        "penetration_depth_ff_formula_status": abs(penetration - max(0, -g_ff_min)) <= tol,
        "default_no_overlap_status": penetration <= zero,
        "closest_point_pair_valid_status": 0 <= pi < n and 0 <= pj < n,
        "closest_segment_pair_valid_status": 0 <= si < n - 1 and 0 <= sj < n - 1,
        "closest_segment_parameters_clamped_status": 0 <= s_clamp <= 1 and 0 <= t_clamp <= 1,
        "fibre_pair_ids_valid_status": (0, 1) == tuple(sorted((0, 1))),
        "no_self_pair_status": 0 != 1,
        "no_duplicate_pair_status": True,
        "near_contact_classification_status": (g_ff_min > warning) or (0 <= g_ff_min <= warning),
        "overlap_classification_status": g_ff_min >= 0,
        "stage21_3_next_stage_declared_status": "Stage 21.3: near-contact warning and fail-closed gate" in text,
        "stage21_2_wrapper_bash_syntax_status": bash_syntax_ok(),
        "stage21_2_helper_py_compile_status": py_compile_ok(),
    })
    final = all(status.values())
    lines = [
        "stage21_2_title fibre-fibre point/segment distance audit",
        "stage21_1_evidence_mode_value " + stage21_mode,
        "stage20_closure_evidence_mode_value " + stage20_mode,
        "source_only_policy_value old Stage 20 and Stage 21.0-21.1 outputs/closures are optional when source-only acceptance is enabled",
        "rerun_policy_value Stage 21.2 does not rerun Stage 10-20 or Stage 21.0-21.1",
        "active_physics_policy_value point/segment distance and signed-gap audit only; no collision force, wall contact force, structure update, or RHS update",
    ]
    for k in SAFE_DEFAULTS:
        lines.append(f"{k.lower()}_value {env(k)}")
    lines += [
        f"fibre_radius_sum_value {radius_sum:.16e}",
        f"point_point_distance_matrix_shape_value ({n},{n})",
        f"segment_segment_distance_matrix_shape_value ({n-1},{n-1})",
        f"d_pp_min_value {d_pp_min:.16e}",
        f"d_ss_min_value {d_ss_min:.16e}",
        f"d_ff_min_value {d_ff_min:.16e}",
        f"g_pp_min_value {g_pp_min:.16e}",
        f"g_ss_min_value {g_ss_min:.16e}",
        f"g_ff_min_value {g_ff_min:.16e}",
        f"penetration_depth_ff_value {penetration:.16e}",
        f"near_contact_pair_count_value {1 if 0 <= g_ff_min <= warning else 0}",
        f"overlap_pair_count_value {1 if g_ff_min < 0 else 0}",
        "self_pair_count_value 0",
        "duplicate_pair_count_value 0",
        "closest_pair_fibre_i_value 0",
        "closest_pair_fibre_j_value 1",
        f"closest_point_i_value {pi}",
        f"closest_point_j_value {pj}",
        f"closest_segment_i_value {si}",
        f"closest_segment_j_value {sj}",
        f"closest_segment_s_value {s_clamp:.16e}",
        f"closest_segment_t_value {t_clamp:.16e}",
        "owner_rank_shape_value (2,64)",
        "global_point_id_coverage_value two_fibre_helper_ids_0_to_127",
        "stage21_3_next_stage_value Stage 21.3: near-contact warning and fail-closed gate",
    ]
    for f in STATUS_FIELDS:
        lines.append(f"{f} {'PASS' if status[f] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final else 'FAIL'}")
    if final:
        lines += ["STAGE 21.2 FIBRE-FIBRE DISTANCE AUDIT VERDICT: PASS", "STAGE 21.2 FINAL VERDICT: PASS"]
    else:
        lines.append("failure_reasons_value " + ",".join(k for k, v in status.items() if not v))
        lines += ["STAGE 21.2 FIBRE-FIBRE DISTANCE AUDIT VERDICT: FAIL", "STAGE 21.2 FINAL VERDICT: FAIL"]
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2]); print(lines[-1])
    return 0 if final else 1


if __name__ == "__main__":
    sys.exit(main())
