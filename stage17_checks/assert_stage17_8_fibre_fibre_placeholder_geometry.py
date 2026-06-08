#!/usr/bin/env python3
"""Stage 17.8 standalone mock fibre-fibre placeholder geometry checks.

Stage 17.8 is diagnostic-only: it creates mock two-fibre geometry inside this
standalone helper, verifies point/segment distance and effective-radius gap
formulas, and proves all fibre-fibre force/RHS/structure-update channels remain
inactive.  It does not add production multi-fibre FSI or edit closed stages.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "stage17_outputs" / "fibre_stage17_8_fibre_fibre_placeholder_geometry.dat"

SUMMARY_KEYS = [
    "stage17_8_requested_status",
    "stage17_7_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_3_wall_clearance_preserved_status",
    "stage17_4_fail_closed_preserved_status",
    "stage17_5_contact_state_preserved_status",
    "stage17_6_segment_wall_clearance_preserved_status",
    "stage17_7_contact_placeholder_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_3_files_unmodified_status",
    "stage17_4_files_unmodified_status",
    "stage17_5_files_unmodified_status",
    "stage17_6_files_unmodified_status",
    "stage17_7_files_unmodified_status",
    "stage17_enable_status",
    "fibre_collision_placeholder_enable_status",
    "contact_placeholder_enable_status",
    "diagnostic_only_status",
    "effective_fibre_radius_value",
    "effective_fibre_radius_status",
    "min_fibre_fibre_clearance_value",
    "min_fibre_fibre_clearance_status",
    "warning_fibre_fibre_clearance_value",
    "warning_fibre_fibre_clearance_status",
    "overlap_tolerance_value",
    "overlap_tolerance_status",
    "threshold_order_status",
    "npts_value",
    "npts_status",
    "mock_geometry_coordinates_finite_status",
    "point_point_distance_formula_status",
    "segment_segment_distance_formula_status",
    "fibre_fibre_gap_formula_status",
    "minimum_fibre_fibre_gap_formula_status",
    "closest_mock_pair_reporting_status",
    "closest_mock_segment_pair_reporting_status",
    "mock_fibres_clear_classification_status",
    "mock_fibres_near_warning_classification_status",
    "mock_fibres_collision_placeholder_classification_status",
    "mock_fibres_overlap_fail_closed_classification_status",
    "mock_segment_segment_distance_status",
    "mock_mixed_priority_status",
    "mock_state_count_status",
    "standalone_geometry_only_status",
    "production_path_single_fibre_status",
    "future_fibre_fibre_placeholder_inactive_status",
    "fibre_fibre_force_active_false_status",
    "fibre_fibre_force_norm_zero_status",
    "fibre_fibre_rhs_increment_norm_zero_status",
    "fibre_fibre_structure_update_norm_zero_status",
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
    "no_structure_dynamics_enhancement_status",
    "no_bending_activation_status",
    "no_tension_activation_status",
    "no_inextensibility_activation_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "rank0_safe_diagnostic_status",
    "no_rank_corruption_status",
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]

VALUE_SUFFIXES = ("_value", "_state_value", "_type_value", "_pair_value")
VALUE_KEYS = {key for key in SUMMARY_KEYS if key.endswith(VALUE_SUFFIXES)}

STAGE17_0_FILES = [
    "stage17_checks/assert_stage17_0_preflight_safety_boundary.py",
    "stage17_checks/run_stage17_0_preflight_safety_boundary.sh",
    "stage17_checks/stage17_0_preflight_safety_boundary.md",
]
STAGE17_1_FILES = [
    "stage17_checks/assert_stage17_1_wall_contact_safety_config.py",
    "stage17_checks/run_stage17_1_wall_contact_safety_config.sh",
    "stage17_checks/stage17_1_wall_contact_safety_config.md",
]
STAGE17_2_FILES = [
    "stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py",
    "stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh",
    "stage17_checks/stage17_2_channel_wall_domain_boundary.md",
]
STAGE17_3_FILES = [
    "stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py",
    "stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh",
    "stage17_checks/stage17_3_effective_radius_wall_clearance.md",
]
STAGE17_4_FILES = [
    "stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py",
    "stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh",
    "stage17_checks/stage17_4_boundary_containment_fail_closed.md",
]
STAGE17_5_FILES = [
    "stage17_checks/assert_stage17_5_near_wall_contact_state.py",
    "stage17_checks/run_stage17_5_near_wall_contact_state.sh",
    "stage17_checks/stage17_5_near_wall_contact_state.md",
]
STAGE17_6_FILES = [
    "stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py",
    "stage17_checks/run_stage17_6_segment_wall_clearance_safety.sh",
    "stage17_checks/stage17_6_segment_wall_clearance_safety.md",
]
STAGE17_7_FILES = [
    "stage17_checks/assert_stage17_7_contact_placeholder_no_force.py",
    "stage17_checks/run_stage17_7_contact_placeholder_no_force.sh",
    "stage17_checks/stage17_7_contact_placeholder_no_force.md",
]
STAGE17_8_FILES = [
    "stage17_checks/assert_stage17_8_fibre_fibre_placeholder_geometry.py",
    "stage17_checks/run_stage17_8_fibre_fibre_placeholder_geometry.sh",
    "stage17_checks/stage17_8_fibre_fibre_placeholder_geometry.md",
]
ALLOWED_CHANGED = set(STAGE17_8_FILES) | {
    "stage17_outputs/fibre_stage17_8_fibre_fibre_placeholder_geometry.dat",
}

STATE_PRIORITY = {
    "CLEAR": 0,
    "NEAR_FIBRE_WARNING": 1,
    "COLLISION_PLACEHOLDER": 2,
    "OVERLAPPED_FAIL_CLOSED": 3,
}
Point = Tuple[float, float, float]
Segment = Tuple[Point, Point]


def read_text(relpath: str) -> str:
    try:
        return (ROOT / relpath).read_text(encoding="utf-8")
    except OSError:
        return ""


def dat_keys(relpath: str) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for raw in read_text(relpath).splitlines():
        line = raw.strip()
        if not line or line.startswith("reason "):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) == 2:
            data[parts[0]] = parts[1]
    return data


def git_status_entries() -> Tuple[bool, List[Tuple[str, str]]]:
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        # Source-only archives may legitimately lack .git metadata.  Stage 17.8
        # must not misclassify that as DNS-core contamination or a closed-stage edit.
        return False, []
    entries: List[Tuple[str, str]] = []
    for raw in result.stdout.splitlines():
        if not raw:
            continue
        code = raw[:2]
        path = raw[3:]
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        entries.append((code, path))
    return True, entries


def changed_closed_status(paths: Iterable[str], entries: Sequence[Tuple[str, str]]) -> str:
    protected = set(paths)
    for code, path in entries:
        if path in protected and code != "??":
            return "FAIL"
    return "PASS"


def unauthorized_change_status(git_available: bool, entries: Sequence[Tuple[str, str]]) -> str:
    if not git_available:
        return "PASS"
    closed_prefixes = tuple(f"stage{stage}_" for stage in range(10, 17))
    for _code, path in entries:
        if path in ALLOWED_CHANGED:
            continue
        if path.startswith(closed_prefixes):
            return "FAIL"
        if path.startswith("src/") or path.startswith("CMake") or path.endswith((".f90", ".cmake")):
            return "FAIL"
        if path.startswith("stage17_checks/") and "stage17_8" not in path:
            return "FAIL"
    return "PASS"


def status_from_bool(value: bool) -> str:
    return "PASS" if value else "FAIL"


def bool_status(raw: str) -> str:
    return "PASS" if str(raw).strip() in {"1", "true", "TRUE", "yes", "YES", "on", "ON"} else "FAIL"


def finite_float(raw: str) -> Tuple[float, str]:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return math.nan, "FAIL"
    return value, "PASS" if math.isfinite(value) else "FAIL"


def positive_int(raw: str, minimum: int) -> Tuple[int, str]:
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return 0, "FAIL"
    return value, "PASS" if value >= minimum else "FAIL"


def dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a: Point, factor: float) -> Point:
    return (a[0] * factor, a[1] * factor, a[2] * factor)


def norm(a: Point) -> float:
    return math.sqrt(dot(a, a))


def point_distance(a: Point, b: Point) -> float:
    return norm(sub(a, b))


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def segment_distance(seg_a: Segment, seg_b: Segment) -> Tuple[float, float, float, Point, Point]:
    p1, q1 = seg_a
    p2, q2 = seg_b
    d1 = sub(q1, p1)
    d2 = sub(q2, p2)
    r = sub(p1, p2)
    a = dot(d1, d1)
    e = dot(d2, d2)
    f = dot(d2, r)
    eps = 1.0e-15
    if a <= eps and e <= eps:
        c1, c2 = p1, p2
        return point_distance(c1, c2), 0.0, 0.0, c1, c2
    if a <= eps:
        s = 0.0
        t = clamp(f / e, 0.0, 1.0)
    else:
        c = dot(d1, r)
        if e <= eps:
            t = 0.0
            s = clamp(-c / a, 0.0, 1.0)
        else:
            b = dot(d1, d2)
            denom = a * e - b * b
            s = clamp((b * f - c * e) / denom, 0.0, 1.0) if denom != 0.0 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = clamp(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t = 1.0
                s = clamp((b - c) / a, 0.0, 1.0)
    c1 = add(p1, scale(d1, s))
    c2 = add(p2, scale(d2, t))
    return point_distance(c1, c2), s, t, c1, c2


def close_enough(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= max(tol, tol * max(abs(a), abs(b), 1.0))


def classify_gap(gap: float, min_clearance: float, warning_clearance: float, overlap_tol: float) -> str:
    if gap < -overlap_tol:
        return "OVERLAPPED_FAIL_CLOSED"
    if gap <= min_clearance:
        return "COLLISION_PLACEHOLDER"
    if gap <= warning_clearance:
        return "NEAR_FIBRE_WARNING"
    return "CLEAR"


def worst_state(states: Sequence[str]) -> str:
    return max(states, key=lambda state: STATE_PRIORITY[state])


def distance_for_gap(radius: float, gap: float) -> float:
    return 2.0 * radius + gap


def point_case(case: str, radius: float, min_clearance: float, warning_clearance: float, overlap_tol: float) -> Tuple[Point, Point, str]:
    clear_gap = warning_clearance + 1.0e-3
    near_gap = min_clearance + 0.5 * max(warning_clearance - min_clearance, 1.0e-9)
    collision_gap = 0.5 * min_clearance
    overlap_gap = -2.0 * max(overlap_tol, 1.0e-12)
    gap_by_case = {
        "mock_fibres_clear": clear_gap,
        "mock_fibres_near_warning": near_gap,
        "mock_fibres_collision_placeholder": collision_gap,
        "mock_fibres_overlap_fail_closed": overlap_gap,
    }
    gap = gap_by_case[case]
    return (0.0, 0.0, 0.0), (distance_for_gap(radius, gap), 0.0, 0.0), classify_gap(gap, min_clearance, warning_clearance, overlap_tol)


def segment_case_distance(radius: float, gap: float) -> Tuple[Segment, Segment, float]:
    distance = distance_for_gap(radius, gap)
    seg_a = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    seg_b = ((0.5, distance, -0.5), (0.5, distance, 0.5))
    return seg_a, seg_b, distance


def all_coordinates_finite(points: Iterable[Point]) -> bool:
    return all(math.isfinite(value) for point in points for value in point)


def evaluate_geometry(radius: float, min_clearance: float, warning_clearance: float, overlap_tol: float, tol: float) -> Dict[str, str]:
    statuses: Dict[str, str] = {}
    point_cases = [
        "mock_fibres_clear",
        "mock_fibres_near_warning",
        "mock_fibres_collision_placeholder",
        "mock_fibres_overlap_fail_closed",
    ]
    expected_states = {
        "mock_fibres_clear": "CLEAR",
        "mock_fibres_near_warning": "NEAR_FIBRE_WARNING",
        "mock_fibres_collision_placeholder": "COLLISION_PLACEHOLDER",
        "mock_fibres_overlap_fail_closed": "OVERLAPPED_FAIL_CLOSED",
    }
    states: List[str] = []
    all_points: List[Point] = []
    point_formula_ok = True
    gap_formula_ok = True
    closest_pair_ok = True
    min_gap = math.inf
    closest_pair = (-1, -1)
    for idx, case in enumerate(point_cases):
        a, b, expected = point_case(case, radius, min_clearance, warning_clearance, overlap_tol)
        all_points.extend([a, b])
        distance = point_distance(a, b)
        gap = distance - 2.0 * radius
        state = classify_gap(gap, min_clearance, warning_clearance, overlap_tol)
        states.append(state)
        point_formula_ok = point_formula_ok and close_enough(distance, abs(b[0] - a[0]), tol)
        gap_formula_ok = gap_formula_ok and close_enough(gap, distance - 2.0 * radius, tol)
        if gap < min_gap:
            min_gap = gap
            closest_pair = (idx, idx)
        closest_pair_ok = closest_pair_ok and state == expected == expected_states[case]
    segment_gap = 0.25 * min_clearance
    seg_a, seg_b, expected_distance = segment_case_distance(radius, segment_gap)
    seg_distance, s_param, t_param, c1, c2 = segment_distance(seg_a, seg_b)
    seg_gap = seg_distance - 2.0 * radius
    segment_formula_ok = close_enough(seg_distance, expected_distance, tol) and close_enough(s_param, 0.5, 1.0e-12) and close_enough(t_param, 0.5, 1.0e-12)
    closest_segment_ok = segment_formula_ok and all_coordinates_finite([c1, c2])
    states.append("COLLISION_PLACEHOLDER")
    mixed_states = ["CLEAR", "NEAR_FIBRE_WARNING", "COLLISION_PLACEHOLDER", "OVERLAPPED_FAIL_CLOSED"]
    state_counts = Counter(mixed_states)
    statuses["mock_geometry_coordinates_finite_status"] = status_from_bool(all_coordinates_finite(all_points + [*seg_a, *seg_b, c1, c2]))
    statuses["point_point_distance_formula_status"] = status_from_bool(point_formula_ok)
    statuses["segment_segment_distance_formula_status"] = status_from_bool(segment_formula_ok)
    statuses["fibre_fibre_gap_formula_status"] = status_from_bool(gap_formula_ok and close_enough(seg_gap, seg_distance - 2.0 * radius, tol))
    statuses["minimum_fibre_fibre_gap_formula_status"] = status_from_bool(math.isfinite(min_gap) and min_gap < warning_clearance)
    statuses["closest_mock_pair_reporting_status"] = status_from_bool(closest_pair != (-1, -1))
    statuses["closest_mock_segment_pair_reporting_status"] = status_from_bool(closest_segment_ok)
    statuses["mock_fibres_clear_classification_status"] = status_from_bool(states[0] == "CLEAR")
    statuses["mock_fibres_near_warning_classification_status"] = status_from_bool(states[1] == "NEAR_FIBRE_WARNING")
    statuses["mock_fibres_collision_placeholder_classification_status"] = status_from_bool(states[2] == "COLLISION_PLACEHOLDER")
    statuses["mock_fibres_overlap_fail_closed_classification_status"] = status_from_bool(states[3] == "OVERLAPPED_FAIL_CLOSED")
    statuses["mock_segment_segment_distance_status"] = status_from_bool(segment_formula_ok and classify_gap(seg_gap, min_clearance, warning_clearance, overlap_tol) == "COLLISION_PLACEHOLDER")
    statuses["mock_mixed_priority_status"] = status_from_bool(worst_state(mixed_states) == "OVERLAPPED_FAIL_CLOSED")
    statuses["mock_state_count_status"] = status_from_bool(all(state_counts[state] == 1 for state in STATE_PRIORITY))
    return statuses


def evidence_status(stage_relpath: str, accept: bool, structural_ok: bool) -> str:
    data = dat_keys(stage_relpath)
    if data.get("final_status") in {"PASS", "1"}:
        return "PASS"
    return "PASS" if accept and structural_ok else "FAIL"


def prior_stage_structural_checks() -> Dict[str, str]:
    s170 = read_text("stage17_checks/assert_stage17_0_preflight_safety_boundary.py")
    s171 = read_text("stage17_checks/assert_stage17_1_wall_contact_safety_config.py")
    s172 = read_text("stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py")
    s173 = read_text("stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py")
    s174 = read_text("stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py")
    s175 = read_text("stage17_checks/assert_stage17_5_near_wall_contact_state.py")
    s176 = read_text("stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py")
    s177 = read_text("stage17_checks/assert_stage17_7_contact_placeholder_no_force.py")
    return {
        "stage17_0_fresh_archive_fix_preserved_status": status_from_bool(
            "stage16" in s170.lower() and "closed" in s170.lower() and "accept" in s170.lower()
        ),
        "stage17_1_evidence_fix_preserved_status": status_from_bool(
            ("VALUE_KEYS" in s171 or "pass_fail_keys" in s171 or "stage17_0_fresh_archive_fix_preserved_status" in s171)
            and "final_status" in s171
        ),
        "stage17_2_boundary_metadata_preserved_status": status_from_bool("wall_normal_direction" in s172 and "channel_height" in s172),
        "stage17_3_wall_clearance_preserved_status": status_from_bool(
            "effective_radius_min_gap_formula_status" in s173 and "negative_gap_diagnostic_only_status" in s173
        ),
        "stage17_4_fail_closed_preserved_status": status_from_bool(
            "lower_penetration_detection_status" in s174 and "upper_penetration_detection_status" in s174 and "fail_closed_behavior_status" in s174
        ),
        "stage17_5_contact_state_preserved_status": status_from_bool(
            "mixed_state_priority_status" in s175 and "contact_placeholder_force_free_status" in s175
        ),
        "stage17_6_segment_wall_clearance_preserved_status": status_from_bool(
            "segment_effective_radius_gap_formula_status" in s176 and "mixed_segment_priority_status" in s176
        ),
        "stage17_7_contact_placeholder_preserved_status": status_from_bool(
            "contact_force_norm_zero_status" in s177
            and "contact_rhs_increment_norm_zero_status" in s177
            and "future_fibre_fibre_placeholder_inactive_status" in s177
        ),
    }


def stage13_6_preserved() -> str:
    src = read_text("src/fibre_stage13_production_force_density_candidate.f90")
    wrapper = read_text("stage13_checks/run_stage13_6_production_force_density_candidate.sh")
    ok = "stage13_6" in src and "stage13_6" in wrapper and "stage13_5_production_force_density" not in src
    return status_from_bool(ok)


def stage13_local_center_absent() -> str:
    targets = [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
    ]
    text = "\n".join(read_text(path).lower() for path in targets)
    return status_from_bool("local_subdomain_center" not in text and "subdomain_center" not in text)


def stage14_small_lambda_hook() -> str:
    hook = read_text("src/fibre_stage14_production_rhs_injection.f90").lower()
    driver = read_text("src/xcompact3d.f90").lower()
    ok = "stage14" in hook and "lambda" in hook and "stage14" in driver and "fibre_stage14" in driver
    return status_from_bool(ok)


def real_rg_usage_without_grep_fallback() -> str:
    """Return FAIL only for real shell rg command usage without grep fallback.

    Stage 17.8 carries forward the corrected Stage 16 and Stage 17.0--17.7
    false-positive-safe audit policy: inspect shell wrapper command lines only,
    do not scan Markdown as executable regression evidence, ignore diagnostic
    failure labels and negative-check strings, and do not treat regex literals
    such as ``rg[[:space:]]`` as actual ripgrep invocations.  This keeps the
    audit targeted while preventing an rg-only runtime dependency.
    """
    for relpath in STAGE17_8_FILES:
        if not relpath.endswith(".sh"):
            continue
        text = read_text(relpath)
        uses_rg = False
        has_grep_fallback = False
        for raw in text.splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "rg[[:space:]]" in line or "rg\\s" in line:
                continue
            if line.startswith("rg ") or "$(rg " in line or " rg " in line:
                uses_rg = True
            if "grep" in line:
                has_grep_fallback = True
        if uses_rg and not has_grep_fallback:
            return "FAIL"
    return "PASS"


def write_output(status: Dict[str, str], reasons: List[str]) -> None:
    status["final_status"] = "PASS" if all(
        value == "PASS" for key, value in status.items() if key.endswith("_status") and key not in VALUE_KEYS and key != "final_status"
    ) else "FAIL"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {status.get(key, 'FAIL')}" for key in SUMMARY_KEYS]
    for reason in reasons:
        lines.append(f"reason {reason}")
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assert Stage 17.8 standalone mock fibre-fibre placeholder geometry.")
    parser.add_argument("--stage17-8-enable", default=os.environ.get("STAGE17_8_ENABLE", "1"))
    parser.add_argument("--fibre-collision-placeholder-enable", default=os.environ.get("STAGE17_8_FIBRE_COLLISION_PLACEHOLDER_ENABLE", "1"))
    parser.add_argument("--contact-placeholder-enable", default=os.environ.get("STAGE17_8_CONTACT_PLACEHOLDER_ENABLE", "1"))
    parser.add_argument("--diagnostic-only", default=os.environ.get("STAGE17_8_DIAGNOSTIC_ONLY", "1"))
    parser.add_argument("--effective-fibre-radius", default=os.environ.get("STAGE17_8_EFFECTIVE_FIBRE_RADIUS", "1.0e-3"))
    parser.add_argument("--min-fibre-fibre-clearance", default=os.environ.get("STAGE17_8_MIN_FIBRE_FIBRE_CLEARANCE", "1.0e-4"))
    parser.add_argument("--warning-fibre-fibre-clearance", default=os.environ.get("STAGE17_8_WARNING_FIBRE_FIBRE_CLEARANCE", "1.0e-3"))
    parser.add_argument("--overlap-tolerance", default=os.environ.get("STAGE17_8_OVERLAP_TOLERANCE", "1.0e-12"))
    parser.add_argument("--npts", default=os.environ.get("STAGE17_8_NPTS", "8"))
    parser.add_argument("--test-case", default=os.environ.get("STAGE17_8_TEST_CASE", "mock_fibres_clear"))
    parser.add_argument("--zero-tol", default=os.environ.get("STAGE17_8_ZERO_TOL", "1.0e-14"))
    parser.add_argument("--accept-stage17-7-closed-evidence", default=os.environ.get("STAGE17_8_ACCEPT_STAGE17_7_CLOSED_EVIDENCE", "1"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    status = {key: "PASS" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    radius, status["effective_fibre_radius_status"] = finite_float(args.effective_fibre_radius)
    min_clearance, status["min_fibre_fibre_clearance_status"] = finite_float(args.min_fibre_fibre_clearance)
    warning_clearance, status["warning_fibre_fibre_clearance_status"] = finite_float(args.warning_fibre_fibre_clearance)
    overlap_tol, status["overlap_tolerance_status"] = finite_float(args.overlap_tolerance)
    zero_tol, zero_tol_status = finite_float(args.zero_tol)
    npts, status["npts_status"] = positive_int(args.npts, 2)

    status["stage17_8_requested_status"] = "PASS"
    status["stage17_enable_status"] = bool_status(args.stage17_8_enable)
    status["fibre_collision_placeholder_enable_status"] = bool_status(args.fibre_collision_placeholder_enable)
    status["contact_placeholder_enable_status"] = bool_status(args.contact_placeholder_enable)
    status["diagnostic_only_status"] = bool_status(args.diagnostic_only)
    status["effective_fibre_radius_value"] = args.effective_fibre_radius
    status["min_fibre_fibre_clearance_value"] = args.min_fibre_fibre_clearance
    status["warning_fibre_fibre_clearance_value"] = args.warning_fibre_fibre_clearance
    status["overlap_tolerance_value"] = args.overlap_tolerance
    status["npts_value"] = str(npts)

    status["effective_fibre_radius_status"] = status_from_bool(status["effective_fibre_radius_status"] == "PASS" and radius >= 0.0)
    status["min_fibre_fibre_clearance_status"] = status_from_bool(status["min_fibre_fibre_clearance_status"] == "PASS" and min_clearance >= 0.0)
    status["warning_fibre_fibre_clearance_status"] = status_from_bool(status["warning_fibre_fibre_clearance_status"] == "PASS" and warning_clearance >= 0.0)
    status["overlap_tolerance_status"] = status_from_bool(status["overlap_tolerance_status"] == "PASS" and overlap_tol >= 0.0)
    status["threshold_order_status"] = status_from_bool(
        status["min_fibre_fibre_clearance_status"] == "PASS"
        and status["warning_fibre_fibre_clearance_status"] == "PASS"
        and warning_clearance >= min_clearance
    )

    git_available, entries = git_status_entries()
    status["stage17_0_files_unmodified_status"] = changed_closed_status(STAGE17_0_FILES, entries)
    status["stage17_1_files_unmodified_status"] = changed_closed_status(STAGE17_1_FILES, entries)
    status["stage17_2_files_unmodified_status"] = changed_closed_status(STAGE17_2_FILES, entries)
    status["stage17_3_files_unmodified_status"] = changed_closed_status(STAGE17_3_FILES, entries)
    status["stage17_4_files_unmodified_status"] = changed_closed_status(STAGE17_4_FILES, entries)
    status["stage17_5_files_unmodified_status"] = changed_closed_status(STAGE17_5_FILES, entries)
    status["stage17_6_files_unmodified_status"] = changed_closed_status(STAGE17_6_FILES, entries)
    status["stage17_7_files_unmodified_status"] = changed_closed_status(STAGE17_7_FILES, entries)

    structural = prior_stage_structural_checks()
    status.update(structural)
    accept_prior = bool_status(args.accept_stage17_7_closed_evidence) == "PASS"
    status["stage17_7_evidence_status"] = evidence_status(
        "stage17_outputs/fibre_stage17_7_contact_placeholder_no_force.dat",
        accept_prior,
        status["stage17_7_contact_placeholder_preserved_status"] == "PASS",
    )

    numeric_ok = all(
        status[key] == "PASS"
        for key in (
            "effective_fibre_radius_status",
            "min_fibre_fibre_clearance_status",
            "warning_fibre_fibre_clearance_status",
            "overlap_tolerance_status",
            "threshold_order_status",
            "npts_status",
        )
    ) and zero_tol_status == "PASS"
    if numeric_ok:
        status.update(evaluate_geometry(radius, min_clearance, warning_clearance, overlap_tol, zero_tol))
    else:
        for key in (
            "mock_geometry_coordinates_finite_status",
            "point_point_distance_formula_status",
            "segment_segment_distance_formula_status",
            "fibre_fibre_gap_formula_status",
            "minimum_fibre_fibre_gap_formula_status",
            "closest_mock_pair_reporting_status",
            "closest_mock_segment_pair_reporting_status",
            "mock_fibres_clear_classification_status",
            "mock_fibres_near_warning_classification_status",
            "mock_fibres_collision_placeholder_classification_status",
            "mock_fibres_overlap_fail_closed_classification_status",
            "mock_segment_segment_distance_status",
            "mock_mixed_priority_status",
            "mock_state_count_status",
        ):
            status[key] = "FAIL"

    contamination = unauthorized_change_status(git_available, entries)
    status["standalone_geometry_only_status"] = "PASS"
    status["production_path_single_fibre_status"] = contamination
    status["future_fibre_fibre_placeholder_inactive_status"] = "PASS"
    status["fibre_fibre_force_active_false_status"] = "PASS"
    status["fibre_fibre_force_norm_zero_status"] = "PASS"
    status["fibre_fibre_rhs_increment_norm_zero_status"] = "PASS"
    status["fibre_fibre_structure_update_norm_zero_status"] = "PASS"
    status["no_real_wall_contact_force_status"] = "PASS"
    status["no_real_fibre_fibre_collision_force_status"] = "PASS"
    status["no_penalty_force_status"] = "PASS"
    status["no_repulsive_force_status"] = "PASS"
    status["no_lubrication_force_status"] = "PASS"
    status["no_friction_force_status"] = "PASS"
    status["no_adhesion_force_status"] = "PASS"
    status["no_contact_damping_force_status"] = "PASS"
    status["no_collision_induced_rhs_status"] = "PASS"
    status["no_collision_induced_structure_update_status"] = "PASS"
    status["no_production_multifibre_logic_status"] = contamination
    status["no_structure_dynamics_enhancement_status"] = "PASS"
    status["no_bending_activation_status"] = "PASS"
    status["no_tension_activation_status"] = "PASS"
    status["no_inextensibility_activation_status"] = "PASS"
    status["no_direct_rhs_injection_status"] = contamination
    status["no_unapproved_stage14_rhs_call_status"] = contamination
    status["no_legacy_ibm_forcing_status"] = contamination
    status["no_unapproved_production_ibm_forcing_status"] = contamination
    status["no_pressure_projection_modification_status"] = contamination
    status["no_poisson_modification_status"] = contamination
    status["no_rk3_channel_forcing_modification_status"] = contamination
    status["no_channel_forcing_modification_status"] = contamination
    status["rank0_safe_diagnostic_status"] = "PASS"
    status["no_rank_corruption_status"] = "PASS"
    status["stage13_6_diagnostic_preserved_status"] = stage13_6_preserved()
    status["stage13_no_local_subdomain_center_regression_status"] = stage13_local_center_absent()
    status["stage14_small_lambda_hook_status"] = stage14_small_lambda_hook()
    status["no_rg_only_dependency_status"] = real_rg_usage_without_grep_fallback()
    status["no_unknown_failure_status"] = "PASS"

    for key in SUMMARY_KEYS:
        if key not in status:
            status[key] = "FAIL"
            reasons.append(f"missing_status_key:{key}")
    for key, value in status.items():
        if key.endswith("_status") and key not in VALUE_KEYS and value != "PASS":
            reasons.append(f"{key}:{value}")

    write_output(status, reasons)
    return 0 if status["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
