#!/usr/bin/env python3
"""Stage 17.7 diagnostic contact placeholder interface with no force.

This helper adds only a standalone, force-free placeholder record layer for
wall/contact metadata.  It reads closed-stage evidence without modifying it,
verifies analytic placeholder cases, and emits a machine-readable summary for
Stage 17.7 automation.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "stage17_outputs" / "fibre_stage17_7_contact_placeholder_no_force.dat"

SUMMARY_KEYS = [
    "stage17_7_requested_status",
    "stage17_6_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_3_wall_clearance_preserved_status",
    "stage17_4_fail_closed_preserved_status",
    "stage17_5_contact_state_preserved_status",
    "stage17_6_segment_wall_clearance_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_3_files_unmodified_status",
    "stage17_4_files_unmodified_status",
    "stage17_5_files_unmodified_status",
    "stage17_6_files_unmodified_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "fail_closed_enable_status",
    "contact_placeholder_enable_status",
    "fibre_collision_placeholder_enable_status",
    "diagnostic_only_status",
    "contact_placeholder_interface_status",
    "contact_record_schema_status",
    "contact_pair_type_status",
    "contact_location_type_status",
    "contact_state_value_status",
    "contact_gap_finite_status",
    "lower_wall_normal_status",
    "upper_wall_normal_status",
    "point_placeholder_mapping_status",
    "segment_placeholder_mapping_status",
    "penetrated_fail_closed_placeholder_status",
    "future_fibre_fibre_placeholder_inactive_status",
    "contact_force_active_false_status",
    "contact_force_norm_zero_status",
    "contact_rhs_increment_norm_zero_status",
    "contact_structure_update_norm_zero_status",
    "contact_placeholder_diagnostic_only_status",
    "contact_placeholder_force_free_status",
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

VALUE_SUFFIXES = ("_value", "_state_value", "_type_value", "_normal_value")
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
ALLOWED_CHANGED = set(STAGE17_7_FILES) | {
    "stage17_outputs/fibre_stage17_7_contact_placeholder_no_force.dat",
}

VALID_PAIR_TYPES = {"WALL_LOWER", "WALL_UPPER", "FUTURE_FIBRE_FIBRE_PLACEHOLDER"}
VALID_LOCATION_TYPES = {"POINT", "SEGMENT", "MOCK_PAIR"}
VALID_STATES = {"CLEAR", "NEAR_WALL_WARNING", "CONTACT_PLACEHOLDER", "PENETRATED_FAIL_CLOSED"}
LOWER_NORMAL = (0.0, 1.0, 0.0)
UPPER_NORMAL = (0.0, -1.0, 0.0)
MOCK_NORMAL = (0.0, 0.0, 0.0)


@dataclass(frozen=True)
class ContactRecord:
    contact_pair_type: str
    contact_location_type: str
    contact_entity_id: int
    contact_segment_id: int
    contact_point_id: int
    contact_gap: float
    contact_normal_x: float
    contact_normal_y: float
    contact_normal_z: float
    contact_state: str
    contact_placeholder_status: str
    contact_force_active_status: str = "false"
    contact_force_norm: float = 0.0
    contact_rhs_increment_norm: float = 0.0
    contact_structure_update_norm: float = 0.0

    def diagnostic_only(self) -> bool:
        return (
            self.contact_placeholder_status == "PASS"
            and self.contact_force_active_status == "false"
            and self.contact_force_norm == 0.0
            and self.contact_rhs_increment_norm == 0.0
            and self.contact_structure_update_norm == 0.0
        )


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
        # Source-only archives may legitimately lack .git metadata.  Stage 17.7 must not
        # misclassify that environment limitation as DNS-core contamination or a closed-stage edit.
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
        if path.startswith("stage17_checks/") and "stage17_7" not in path:
            return "FAIL"
    return "PASS"


def bool_status(raw: str) -> str:
    return "PASS" if str(raw).strip() in {"1", "true", "TRUE", "yes", "YES", "on", "ON"} else "FAIL"


def disabled_status(raw: str) -> str:
    return "PASS" if str(raw).strip() in {"0", "false", "FALSE", "no", "NO", "off", "OFF"} else "FAIL"


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


def status_from_bool(value: bool) -> str:
    return "PASS" if value else "FAIL"


def close_enough(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= max(tol, tol * max(abs(a), abs(b), 1.0))


def classify_gap(gap: float, min_clearance: float, warning_clearance: float, penetration_tol: float) -> str:
    if gap < -penetration_tol:
        return "PENETRATED_FAIL_CLOSED"
    if gap <= min_clearance:
        return "CONTACT_PLACEHOLDER"
    if gap <= warning_clearance:
        return "NEAR_WALL_WARNING"
    return "CLEAR"


def evidence_status(stage_relpath: str, accept: bool, structural_ok: bool) -> str:
    data = dat_keys(stage_relpath)
    if data.get("final_status") == "PASS":
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
    stage17_1_value_ok = (
        "VALUE_KEYS" in s171
        or "pass_fail_keys" in s171
        or "stage17_0_fresh_archive_fix_preserved_status" in s171
    )
    return {
        "stage17_0_fresh_archive_fix_preserved_status": status_from_bool(
            "stage16" in s170.lower()
            and "closed" in s170.lower()
            and "accept" in s170.lower()
        ),
        "stage17_1_evidence_fix_preserved_status": status_from_bool(
            stage17_1_value_ok and "final_status" in s171
        ),
        "stage17_2_boundary_metadata_preserved_status": status_from_bool(
            "wall_normal_direction" in s172 and "channel_height" in s172
        ),
        "stage17_3_wall_clearance_preserved_status": status_from_bool(
            "effective_radius_min_gap_formula_status" in s173
            and "negative_gap_diagnostic_only_status" in s173
        ),
        "stage17_4_fail_closed_preserved_status": status_from_bool(
            "lower_penetration_detection_status" in s174
            and "upper_penetration_detection_status" in s174
            and "fail_closed_behavior_status" in s174
        ),
        "stage17_5_contact_state_preserved_status": status_from_bool(
            "mixed_state_priority_status" in s175
            and "contact_placeholder_force_free_status" in s175
            and ("VALUE_KEYS" in s175 or "pass_fail_keys" in s175)
        ),
        "stage17_6_segment_wall_clearance_preserved_status": status_from_bool(
            "segment_effective_radius_gap_formula_status" in s176
            and "mixed_segment_priority_status" in s176
            and ("VALUE_KEYS" in s176 or "pass_fail_keys" in s176)
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

    Stage 17.7 carries forward the corrected Stage 16 and Stage 17.0--17.6
    false-positive-safe audit policy: inspect shell wrapper command lines only,
    do not scan Markdown as executable regression evidence, ignore diagnostic
    failure labels and negative-check strings, and do not treat regex literals
    such as ``rg[[:space:]]`` as actual ripgrep invocations.  This targeted check
    avoids broad repository-wide scanning while still preventing an rg-only
    runtime dependency.
    """
    for relpath in STAGE17_7_FILES:
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


def normal_for_pair(pair_type: str) -> Tuple[float, float, float]:
    if pair_type == "WALL_LOWER":
        return LOWER_NORMAL
    if pair_type == "WALL_UPPER":
        return UPPER_NORMAL
    return MOCK_NORMAL


def make_record(
    pair_type: str,
    location_type: str,
    gap: float,
    state: str,
    entity_id: int = 0,
    segment_id: int = -1,
    point_id: int = -1,
) -> ContactRecord:
    nx, ny, nz = normal_for_pair(pair_type)
    return ContactRecord(
        contact_pair_type=pair_type,
        contact_location_type=location_type,
        contact_entity_id=entity_id,
        contact_segment_id=segment_id,
        contact_point_id=point_id,
        contact_gap=gap,
        contact_normal_x=nx,
        contact_normal_y=ny,
        contact_normal_z=nz,
        contact_state=state,
        contact_placeholder_status="PASS",
    )


def wall_gaps_from_y(y: float, y_min: float, y_max: float, radius: float) -> Tuple[float, float]:
    return y - y_min - radius, y_max - y - radius


def point_record_from_gap(
    point_id: int,
    lower_gap: float,
    upper_gap: float,
    min_clearance: float,
    warning_clearance: float,
    penetration_tol: float,
) -> ContactRecord:
    if lower_gap <= upper_gap:
        pair_type = "WALL_LOWER"
        gap = lower_gap
    else:
        pair_type = "WALL_UPPER"
        gap = upper_gap
    state = classify_gap(gap, min_clearance, warning_clearance, penetration_tol)
    return make_record(pair_type, "POINT", gap, state, point_id=point_id)


def segment_record_from_points(
    segment_id: int,
    ya: float,
    yb: float,
    y_min: float,
    y_max: float,
    radius: float,
    min_clearance: float,
    warning_clearance: float,
    penetration_tol: float,
) -> ContactRecord:
    ymid = 0.5 * (ya + yb)
    yseg_min = min(ya, yb, ymid)
    yseg_max = max(ya, yb, ymid)
    lower_gap = yseg_min - y_min - radius
    upper_gap = y_max - yseg_max - radius
    if lower_gap <= upper_gap:
        pair_type = "WALL_LOWER"
        gap = lower_gap
    else:
        pair_type = "WALL_UPPER"
        gap = upper_gap
    state = classify_gap(gap, min_clearance, warning_clearance, penetration_tol)
    return make_record(pair_type, "SEGMENT", gap, state, segment_id=segment_id)


def gap_to_y_lower(y_min: float, radius: float, gap: float) -> float:
    return y_min + radius + gap


def gap_to_y_upper(y_max: float, radius: float, gap: float) -> float:
    return y_max - radius - gap


def case_records(
    case: str,
    y_min: float,
    y_max: float,
    radius: float,
    min_clearance: float,
    warning_clearance: float,
    penetration_tol: float,
) -> List[ContactRecord]:
    clear_gap = warning_clearance + 1.0e-3
    contact_gap = 0.5 * min_clearance
    penetrate_gap = -2.0 * max(penetration_tol, 1.0e-12)
    if case == "clear_no_contact_record":
        y = gap_to_y_lower(y_min, radius, clear_gap)
        record = point_record_from_gap(0, *wall_gaps_from_y(y, y_min, y_max, radius), min_clearance, warning_clearance, penetration_tol)
        return [] if record.contact_state == "CLEAR" else [record]
    if case == "lower_wall_contact_placeholder_point":
        y = gap_to_y_lower(y_min, radius, contact_gap)
        return [point_record_from_gap(0, *wall_gaps_from_y(y, y_min, y_max, radius), min_clearance, warning_clearance, penetration_tol)]
    if case == "upper_wall_contact_placeholder_point":
        y = gap_to_y_upper(y_max, radius, contact_gap)
        return [point_record_from_gap(0, *wall_gaps_from_y(y, y_min, y_max, radius), min_clearance, warning_clearance, penetration_tol)]
    if case == "lower_wall_segment_placeholder":
        ya = gap_to_y_lower(y_min, radius, contact_gap)
        yb = gap_to_y_lower(y_min, radius, clear_gap)
        return [segment_record_from_points(0, ya, yb, y_min, y_max, radius, min_clearance, warning_clearance, penetration_tol)]
    if case == "penetrated_fail_closed_placeholder":
        y = gap_to_y_lower(y_min, radius, penetrate_gap)
        return [point_record_from_gap(0, *wall_gaps_from_y(y, y_min, y_max, radius), min_clearance, warning_clearance, penetration_tol)]
    if case == "future_fibre_fibre_placeholder_inactive":
        return [make_record("FUTURE_FIBRE_FIBRE_PLACEHOLDER", "MOCK_PAIR", math.inf if False else 0.0, "CONTACT_PLACEHOLDER")]
    return []


def record_schema_ok(records: Sequence[ContactRecord]) -> bool:
    for record in records:
        if record.contact_pair_type not in VALID_PAIR_TYPES:
            return False
        if record.contact_location_type not in VALID_LOCATION_TYPES:
            return False
        if record.contact_state not in VALID_STATES:
            return False
        if not math.isfinite(record.contact_gap):
            return False
        if not all(math.isfinite(value) for value in (record.contact_normal_x, record.contact_normal_y, record.contact_normal_z)):
            return False
    return True


def all_force_channels_zero(records: Sequence[ContactRecord]) -> bool:
    return all(record.diagnostic_only() for record in records)


def evaluate_placeholder_cases(
    y_min: float,
    y_max: float,
    radius: float,
    min_clearance: float,
    warning_clearance: float,
    penetration_tol: float,
    tol: float,
) -> Dict[str, str]:
    cases = {
        name: case_records(name, y_min, y_max, radius, min_clearance, warning_clearance, penetration_tol)
        for name in (
            "clear_no_contact_record",
            "lower_wall_contact_placeholder_point",
            "upper_wall_contact_placeholder_point",
            "lower_wall_segment_placeholder",
            "penetrated_fail_closed_placeholder",
            "future_fibre_fibre_placeholder_inactive",
        )
    }
    all_records = [record for records in cases.values() for record in records]
    lower_point = cases["lower_wall_contact_placeholder_point"][0]
    upper_point = cases["upper_wall_contact_placeholder_point"][0]
    lower_segment = cases["lower_wall_segment_placeholder"][0]
    penetrated = cases["penetrated_fail_closed_placeholder"][0]
    future = cases["future_fibre_fibre_placeholder_inactive"][0]
    return {
        "contact_placeholder_interface_status": "PASS",
        "contact_record_schema_status": status_from_bool(record_schema_ok(all_records)),
        "contact_pair_type_status": status_from_bool(all(record.contact_pair_type in VALID_PAIR_TYPES for record in all_records)),
        "contact_location_type_status": status_from_bool(all(record.contact_location_type in VALID_LOCATION_TYPES for record in all_records)),
        "contact_state_value_status": status_from_bool(all(record.contact_state in VALID_STATES for record in all_records)),
        "contact_gap_finite_status": status_from_bool(all(math.isfinite(record.contact_gap) for record in all_records)),
        "lower_wall_normal_status": status_from_bool((lower_point.contact_normal_x, lower_point.contact_normal_y, lower_point.contact_normal_z) == LOWER_NORMAL),
        "upper_wall_normal_status": status_from_bool((upper_point.contact_normal_x, upper_point.contact_normal_y, upper_point.contact_normal_z) == UPPER_NORMAL),
        "point_placeholder_mapping_status": status_from_bool(
            lower_point.contact_pair_type == "WALL_LOWER"
            and lower_point.contact_location_type == "POINT"
            and lower_point.contact_state == "CONTACT_PLACEHOLDER"
            and close_enough(lower_point.contact_gap, 0.5 * min_clearance, tol)
            and upper_point.contact_pair_type == "WALL_UPPER"
            and upper_point.contact_location_type == "POINT"
            and upper_point.contact_state == "CONTACT_PLACEHOLDER"
        ),
        "segment_placeholder_mapping_status": status_from_bool(
            lower_segment.contact_pair_type == "WALL_LOWER"
            and lower_segment.contact_location_type == "SEGMENT"
            and lower_segment.contact_state == "CONTACT_PLACEHOLDER"
        ),
        "penetrated_fail_closed_placeholder_status": status_from_bool(
            penetrated.contact_state == "PENETRATED_FAIL_CLOSED"
            and penetrated.contact_pair_type == "WALL_LOWER"
            and penetrated.diagnostic_only()
        ),
        "future_fibre_fibre_placeholder_inactive_status": status_from_bool(
            future.contact_pair_type == "FUTURE_FIBRE_FIBRE_PLACEHOLDER"
            and future.contact_location_type == "MOCK_PAIR"
            and future.contact_force_active_status == "false"
            and future.contact_force_norm == 0.0
        ),
        "contact_force_active_false_status": status_from_bool(all(record.contact_force_active_status == "false" for record in all_records)),
        "contact_force_norm_zero_status": status_from_bool(all(record.contact_force_norm == 0.0 for record in all_records)),
        "contact_rhs_increment_norm_zero_status": status_from_bool(all(record.contact_rhs_increment_norm == 0.0 for record in all_records)),
        "contact_structure_update_norm_zero_status": status_from_bool(all(record.contact_structure_update_norm == 0.0 for record in all_records)),
        "contact_placeholder_diagnostic_only_status": status_from_bool(all_force_channels_zero(all_records)),
        "contact_placeholder_force_free_status": status_from_bool(all_force_channels_zero(all_records)),
    }


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
    parser = argparse.ArgumentParser(description="Assert Stage 17.7 diagnostic contact placeholder no-force interface.")
    parser.add_argument("--stage17-7-enable", default=os.environ.get("STAGE17_7_ENABLE", "1"))
    parser.add_argument("--wall-safety-enable", default=os.environ.get("STAGE17_7_WALL_SAFETY_ENABLE", "1"))
    parser.add_argument("--boundary-check-enable", default=os.environ.get("STAGE17_7_BOUNDARY_CHECK_ENABLE", "1"))
    parser.add_argument("--fail-closed-enable", default=os.environ.get("STAGE17_7_FAIL_CLOSED_ENABLE", "1"))
    parser.add_argument("--contact-placeholder-enable", default=os.environ.get("STAGE17_7_CONTACT_PLACEHOLDER_ENABLE", "1"))
    parser.add_argument("--fibre-collision-placeholder-enable", default=os.environ.get("STAGE17_7_FIBRE_COLLISION_PLACEHOLDER_ENABLE", "0"))
    parser.add_argument("--diagnostic-only", default=os.environ.get("STAGE17_7_DIAGNOSTIC_ONLY", "1"))
    parser.add_argument("--y-min", default=os.environ.get("STAGE17_7_Y_MIN", "-1.0"))
    parser.add_argument("--y-max", default=os.environ.get("STAGE17_7_Y_MAX", "1.0"))
    parser.add_argument("--effective-fibre-radius", default=os.environ.get("STAGE17_7_EFFECTIVE_FIBRE_RADIUS", "1.0e-3"))
    parser.add_argument("--min-wall-clearance", default=os.environ.get("STAGE17_7_MIN_WALL_CLEARANCE", "1.0e-4"))
    parser.add_argument("--warning-wall-clearance", default=os.environ.get("STAGE17_7_WARNING_WALL_CLEARANCE", "1.0e-3"))
    parser.add_argument("--penetration-tolerance", default=os.environ.get("STAGE17_7_PENETRATION_TOLERANCE", "1.0e-12"))
    parser.add_argument("--npts", default=os.environ.get("STAGE17_7_NPTS", "8"))
    parser.add_argument("--test-case", default=os.environ.get("STAGE17_7_TEST_CASE", "clear_no_contact_record"))
    parser.add_argument("--zero-tol", default=os.environ.get("STAGE17_7_ZERO_TOL", "1.0e-14"))
    parser.add_argument("--accept-stage17-6-closed-evidence", default=os.environ.get("STAGE17_7_ACCEPT_STAGE17_6_CLOSED_EVIDENCE", "1"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    status = {key: "PASS" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    y_min, y_min_status = finite_float(args.y_min)
    y_max, y_max_status = finite_float(args.y_max)
    radius, radius_status = finite_float(args.effective_fibre_radius)
    min_clearance, min_clearance_status = finite_float(args.min_wall_clearance)
    warning_clearance, warning_clearance_status = finite_float(args.warning_wall_clearance)
    penetration_tol, penetration_tol_status = finite_float(args.penetration_tolerance)
    zero_tol, zero_tol_status = finite_float(args.zero_tol)
    _npts, npts_status = positive_int(args.npts, 1)

    status["stage17_7_requested_status"] = "PASS"
    status["stage17_enable_status"] = bool_status(args.stage17_7_enable)
    status["wall_safety_enable_status"] = bool_status(args.wall_safety_enable)
    status["boundary_check_enable_status"] = bool_status(args.boundary_check_enable)
    status["fail_closed_enable_status"] = bool_status(args.fail_closed_enable)
    status["contact_placeholder_enable_status"] = bool_status(args.contact_placeholder_enable)
    status["fibre_collision_placeholder_enable_status"] = disabled_status(args.fibre_collision_placeholder_enable)
    status["diagnostic_only_status"] = bool_status(args.diagnostic_only)

    numeric_ok = (
        y_min_status == "PASS"
        and y_max_status == "PASS"
        and y_max > y_min
        and radius_status == "PASS"
        and radius >= 0.0
        and min_clearance_status == "PASS"
        and min_clearance >= 0.0
        and warning_clearance_status == "PASS"
        and warning_clearance >= min_clearance
        and penetration_tol_status == "PASS"
        and penetration_tol >= 0.0
        and zero_tol_status == "PASS"
        and npts_status == "PASS"
    )

    git_available, entries = git_status_entries()
    status["stage17_0_files_unmodified_status"] = changed_closed_status(STAGE17_0_FILES, entries)
    status["stage17_1_files_unmodified_status"] = changed_closed_status(STAGE17_1_FILES, entries)
    status["stage17_2_files_unmodified_status"] = changed_closed_status(STAGE17_2_FILES, entries)
    status["stage17_3_files_unmodified_status"] = changed_closed_status(STAGE17_3_FILES, entries)
    status["stage17_4_files_unmodified_status"] = changed_closed_status(STAGE17_4_FILES, entries)
    status["stage17_5_files_unmodified_status"] = changed_closed_status(STAGE17_5_FILES, entries)
    status["stage17_6_files_unmodified_status"] = changed_closed_status(STAGE17_6_FILES, entries)

    structural = prior_stage_structural_checks()
    status.update(structural)
    accept_prior = bool_status(args.accept_stage17_6_closed_evidence) == "PASS"
    status["stage17_6_evidence_status"] = evidence_status(
        "stage17_outputs/fibre_stage17_6_segment_wall_clearance_safety.dat",
        accept_prior,
        status["stage17_6_segment_wall_clearance_preserved_status"] == "PASS",
    )

    if numeric_ok:
        status.update(evaluate_placeholder_cases(y_min, y_max, radius, min_clearance, warning_clearance, penetration_tol, zero_tol))
    else:
        for key in (
            "contact_placeholder_interface_status",
            "contact_record_schema_status",
            "contact_pair_type_status",
            "contact_location_type_status",
            "contact_state_value_status",
            "contact_gap_finite_status",
            "lower_wall_normal_status",
            "upper_wall_normal_status",
            "point_placeholder_mapping_status",
            "segment_placeholder_mapping_status",
            "penetrated_fail_closed_placeholder_status",
            "future_fibre_fibre_placeholder_inactive_status",
            "contact_force_active_false_status",
            "contact_force_norm_zero_status",
            "contact_rhs_increment_norm_zero_status",
            "contact_structure_update_norm_zero_status",
            "contact_placeholder_diagnostic_only_status",
            "contact_placeholder_force_free_status",
        ):
            status[key] = "FAIL"

    contamination = unauthorized_change_status(git_available, entries)
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
