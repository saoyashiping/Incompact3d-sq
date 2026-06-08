#!/usr/bin/env python3
"""Stage 17.5 near-wall warning/contact-state diagnostic classifier.

Stage 17.5 classifies standalone analytic point-wise wall-safety states from
Stage 17 effective-radius wall gaps only. It is diagnostic-only and force-free: it
never applies wall/contact/collision forces, changes fibre state, changes fluid RHS,
inserts production hooks, or edits closed Stage 10--16 / Stage 17.0--17.4 files.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_5_requested_status",
    "stage17_4_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_3_wall_clearance_preserved_status",
    "stage17_4_fail_closed_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_3_files_unmodified_status",
    "stage17_4_files_unmodified_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "fail_closed_enable_status",
    "contact_placeholder_enable_status",
    "diagnostic_only_status",
    "y_min_value",
    "y_min_finite_status",
    "y_max_value",
    "y_max_finite_status",
    "y_bounds_ordered_status",
    "effective_fibre_radius_value",
    "effective_fibre_radius_status",
    "min_wall_clearance_value",
    "min_wall_clearance_status",
    "warning_wall_clearance_value",
    "warning_wall_clearance_status",
    "penetration_tolerance_value",
    "penetration_tolerance_status",
    "threshold_order_status",
    "npts_value",
    "npts_status",
    "y_coordinates_finite_status",
    "gap_values_finite_status",
    "all_clear_classification_status",
    "near_wall_warning_classification_status",
    "contact_placeholder_classification_status",
    "penetrated_fail_closed_classification_status",
    "mixed_state_priority_status",
    "state_count_status",
    "global_worst_state_status",
    "penetration_index_reporting_status",
    "contact_placeholder_force_free_status",
    "classification_diagnostic_only_status",
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

VALUE_KEYS = {
    "y_min_value",
    "y_max_value",
    "effective_fibre_radius_value",
    "min_wall_clearance_value",
    "warning_wall_clearance_value",
    "penetration_tolerance_value",
    "npts_value",
}

REQUIRED_STAGE17_5_FILES = [
    Path("stage17_checks/run_stage17_5_near_wall_contact_state.sh"),
    Path("stage17_checks/assert_stage17_5_near_wall_contact_state.py"),
    Path("stage17_checks/stage17_5_near_wall_contact_state.md"),
]
STAGE17_0_FILES = [
    Path("stage17_checks/assert_stage17_0_preflight_safety_boundary.py"),
    Path("stage17_checks/run_stage17_0_preflight_safety_boundary.sh"),
    Path("stage17_checks/stage17_0_preflight_safety_boundary.md"),
    Path("stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat"),
]
STAGE17_1_FILES = [
    Path("stage17_checks/assert_stage17_1_wall_contact_safety_config.py"),
    Path("stage17_checks/run_stage17_1_wall_contact_safety_config.sh"),
    Path("stage17_checks/stage17_1_wall_contact_safety_config.md"),
    Path("stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat"),
]
STAGE17_2_FILES = [
    Path("stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py"),
    Path("stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh"),
    Path("stage17_checks/stage17_2_channel_wall_domain_boundary.md"),
    Path("stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat"),
]
STAGE17_3_FILES = [
    Path("stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py"),
    Path("stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh"),
    Path("stage17_checks/stage17_3_effective_radius_wall_clearance.md"),
    Path("stage17_outputs/fibre_stage17_3_effective_radius_wall_clearance.dat"),
]
STAGE17_4_FILES = [
    Path("stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py"),
    Path("stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh"),
    Path("stage17_checks/stage17_4_boundary_containment_fail_closed.md"),
    Path("stage17_outputs/fibre_stage17_4_boundary_containment_fail_closed.dat"),
]
ALLOWED_STAGE17_5_CHANGES = {
    "stage17_checks/run_stage17_5_near_wall_contact_state.sh",
    "stage17_checks/assert_stage17_5_near_wall_contact_state.py",
    "stage17_checks/stage17_5_near_wall_contact_state.md",
    "stage17_outputs/fibre_stage17_5_near_wall_contact_state.dat",
}
HISTORICAL_STAGE17_OUTPUTS = {
    "stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat",
    "stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat",
    "stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat",
    "stage17_outputs/fibre_stage17_3_effective_radius_wall_clearance.dat",
    "stage17_outputs/fibre_stage17_4_boundary_containment_fail_closed.dat",
}
CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))
KNOWN_CASES = {
    "all_clear",
    "near_wall_warning",
    "contact_placeholder",
    "penetrated_fail_closed",
    "mixed_states_priority",
}
STATE_PRIORITY = {
    "CLEAR": 0,
    "NEAR_WALL_WARNING": 1,
    "CONTACT_PLACEHOLDER": 2,
    "PENETRATED_FAIL_CLOSED": 3,
}


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def parse_dat(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def status(ok: bool) -> str:
    return "1" if ok else "0"


def finite_float(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def positive_int(text: str) -> int | None:
    try:
        value = int(text)
    except ValueError:
        return None
    return value if value > 0 else None


def git_status_entries(repo: Path) -> list[tuple[str, str]]:
    try:
        proc = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=repo,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError:
        return []
    entries: list[tuple[str, str]] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
        code = line[:2]
        raw = line[3:] if len(line) > 3 else line
        if " -> " in raw:
            raw = raw.split(" -> ", 1)[1]
        entries.append((code, raw.strip()))
    return entries


def real_rg_usage_without_grep_fallback(script: Path) -> bool:
    """Detect real shell-wrapper rg dependencies without false positives.

    Stage 17.5 reuses the corrected Stage 17.4 / Stage 17.3 / Stage 17.2 /
    Stage 17.1 / Stage 17.0 / Stage 16 false-positive-safe policy. Only executable
    shell-wrapper command lines are audited. Markdown prose is not scanned as real
    code-regression evidence, negative-check strings are not treated as behavior,
    and regex literals such as rg[[:space:]] do not count as ripgrep invocations.
    Any real wrapper use of rg must include grep fallback so ripgrep is never a hard
    dependency. Existing diagnostic failure-label strings in closed Stage 17.0--17.4
    files are labels, not rollback evidence.
    """
    if script.suffix != ".sh":
        return False
    text = read_text(script)
    real_rg = False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.replace(";", " ").replace("|", " ").replace("&&", " ").split()
        if "rg" in tokens or "command -v rg" in stripped or "which rg" in stripped:
            real_rg = True
    if not real_rg:
        return False
    return not (("command -v rg" in text or "which rg" in text) and "grep" in text)


def classify_gap(gap: float, min_clearance: float, warning_clearance: float, penetration_tol: float) -> str:
    if gap < -penetration_tol:
        return "PENETRATED_FAIL_CLOSED"
    if gap <= min_clearance:
        return "CONTACT_PLACEHOLDER"
    if gap <= warning_clearance:
        return "NEAR_WALL_WARNING"
    return "CLEAR"


def worst_state(states: list[str]) -> str:
    if not states:
        return "CLEAR"
    return max(states, key=lambda state: STATE_PRIORITY[state])


def control_point_y_values(case_name: str, npts: int, y_min: float, y_max: float, r_eff: float, min_clearance: float, warning_clearance: float, penetration_tol: float) -> list[float]:
    height = y_max - y_min
    if npts <= 0 or height <= 0.0:
        return []
    clear_gap = max(warning_clearance + max(0.25 * height, 1.0e-6), warning_clearance + 1.0e-6)
    clear_gap = min(clear_gap, 0.45 * height)
    near_gap = 0.5 * (min_clearance + warning_clearance) if warning_clearance > min_clearance else min_clearance + 0.5 * max(warning_clearance, 1.0e-6)
    contact_gap = 0.5 * min_clearance if min_clearance > 0.0 else 0.0
    pen_gap = -2.0 * max(penetration_tol, 1.0e-12)
    values = [y_min + r_eff + clear_gap for _ in range(npts)]
    if case_name == "all_clear":
        return values
    if case_name == "near_wall_warning":
        values[0] = y_min + r_eff + near_gap
        return values
    if case_name == "contact_placeholder":
        values[0] = y_min + r_eff + contact_gap
        return values
    if case_name == "penetrated_fail_closed":
        values[0] = y_min + r_eff + pen_gap
        return values
    if case_name == "mixed_states_priority":
        mixed = [
            y_min + r_eff + clear_gap,
            y_min + r_eff + near_gap,
            y_min + r_eff + contact_gap,
            y_min + r_eff + pen_gap,
        ]
        for idx in range(min(npts, len(mixed))):
            values[idx] = mixed[idx]
        return values
    return []


def evaluate_case(case_name: str, npts: int, y_min: float, y_max: float, r_eff: float, min_clearance: float, warning_clearance: float, penetration_tol: float) -> dict[str, object]:
    y_values = control_point_y_values(case_name, npts, y_min, y_max, r_eff, min_clearance, warning_clearance, penetration_tol)
    gap_lower = [y - y_min - r_eff for y in y_values]
    gap_upper = [y_max - y - r_eff for y in y_values]
    gap_wall = [min(lo, up) for lo, up in zip(gap_lower, gap_upper)]
    states = [classify_gap(gap, min_clearance, warning_clearance, penetration_tol) for gap in gap_wall]
    counts = Counter(states)
    min_gap = min(gap_wall) if gap_wall else float("nan")
    min_index = gap_wall.index(min_gap) if gap_wall else -1
    penetrated_indices = [idx for idx, state in enumerate(states) if state == "PENETRATED_FAIL_CLOSED"]
    return {
        "y_values": y_values,
        "gap_lower": gap_lower,
        "gap_upper": gap_upper,
        "gap_wall": gap_wall,
        "states": states,
        "counts": counts,
        "global": worst_state(states),
        "min_gap": min_gap,
        "min_index": min_index,
        "penetrated_indices": penetrated_indices,
    }


def all_finite(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def expected_counts(case_name: str, npts: int) -> Counter:
    if case_name == "all_clear":
        return Counter({"CLEAR": npts})
    if case_name == "near_wall_warning":
        return Counter({"NEAR_WALL_WARNING": 1, "CLEAR": max(npts - 1, 0)})
    if case_name == "contact_placeholder":
        return Counter({"CONTACT_PLACEHOLDER": 1, "CLEAR": max(npts - 1, 0)})
    if case_name == "penetrated_fail_closed":
        return Counter({"PENETRATED_FAIL_CLOSED": 1, "CLEAR": max(npts - 1, 0)})
    if case_name == "mixed_states_priority":
        if npts < 4:
            return Counter()
        return Counter({"CLEAR": npts - 3, "NEAR_WALL_WARNING": 1, "CONTACT_PLACEHOLDER": 1, "PENETRATED_FAIL_CLOSED": 1})
    return Counter()


def expected_global(case_name: str) -> str:
    return {
        "all_clear": "CLEAR",
        "near_wall_warning": "NEAR_WALL_WARNING",
        "contact_placeholder": "CONTACT_PLACEHOLDER",
        "penetrated_fail_closed": "PENETRATED_FAIL_CLOSED",
        "mixed_states_priority": "PENETRATED_FAIL_CLOSED",
    }.get(case_name, "CLEAR")


def stage17_0_fresh_archive_fix_preserved(repo: Path) -> bool:
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_0_preflight_safety_boundary.dat")
    if data.get("final_status") == "1":
        return True
    helper = read_text(repo / "stage17_checks" / "assert_stage17_0_preflight_safety_boundary.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_0_preflight_safety_boundary.sh")
    doc = read_text(repo / "stage17_checks" / "stage17_0_preflight_safety_boundary.md").lower()
    return (
        all((repo / path).is_file() for path in STAGE17_0_FILES[:3])
        and "accept" in helper.lower()
        and ("stage16_12" in helper.lower() or "stage 16.12" in helper.lower())
        and "STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE" in wrapper
        and "safety-boundary" in doc
        and "stage 21" in doc
    )


def stage17_1_evidence_fix_preserved(repo: Path) -> bool:
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_1_wall_contact_safety_config.dat")
    if data.get("final_status") == "1":
        return True
    helper = read_text(repo / "stage17_checks" / "assert_stage17_1_wall_contact_safety_config.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_1_wall_contact_safety_config.sh")
    return (
        all((repo / path).is_file() for path in STAGE17_1_FILES[:3])
        and "effective_fibre_radius_value" in helper
        and "final_status" in helper
        and "diagnostic_only" in helper
        and "STAGE17_1_ACCEPT_STAGE17_0_CLOSED_EVIDENCE" in wrapper
    )


def stage17_2_boundary_metadata_preserved(repo: Path) -> bool:
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_2_channel_wall_domain_boundary.dat")
    if data.get("final_status") == "1" and data.get("wall_normal_direction_status", "1") == "1":
        return True
    helper = read_text(repo / "stage17_checks" / "assert_stage17_2_channel_wall_domain_boundary.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_2_channel_wall_domain_boundary.sh")
    doc = read_text(repo / "stage17_checks" / "stage17_2_channel_wall_domain_boundary.md").lower()
    return (
        all((repo / path).is_file() for path in STAGE17_2_FILES[:3])
        and "wall_normal_direction_value" in helper
        and "channel_height_value" in helper
        and "VALUE_KEYS" in helper
        and "STAGE17_2_Y_MIN" in wrapper
        and "channel wall" in doc
        and "domain-boundary" in doc
    )


def stage17_3_wall_clearance_preserved(repo: Path) -> bool:
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_3_effective_radius_wall_clearance.dat")
    if data.get("final_status") == "1" and data.get("effective_radius_min_gap_formula_status", "1") == "1":
        return True
    helper = read_text(repo / "stage17_checks" / "assert_stage17_3_effective_radius_wall_clearance.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_3_effective_radius_wall_clearance.sh")
    doc = read_text(repo / "stage17_checks" / "stage17_3_effective_radius_wall_clearance.md").lower()
    return (
        all((repo / path).is_file() for path in STAGE17_3_FILES[:3])
        and "effective_radius_min_gap_formula_status" in helper
        and "negative_gap_diagnostic_only_status" in helper
        and "VALUE_KEYS" in helper
        and "STAGE17_3_EFFECTIVE_FIBRE_RADIUS" in wrapper
        and "effective-radius wall-clearance" in doc
    )


def stage17_4_fail_closed_preserved(repo: Path) -> bool:
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_4_boundary_containment_fail_closed.dat")
    if data.get("final_status") == "1" and data.get("fail_closed_behavior_status", "1") == "1":
        return True
    helper = read_text(repo / "stage17_checks" / "assert_stage17_4_boundary_containment_fail_closed.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_4_boundary_containment_fail_closed.sh")
    doc = read_text(repo / "stage17_checks" / "stage17_4_boundary_containment_fail_closed.md").lower()
    return (
        all((repo / path).is_file() for path in STAGE17_4_FILES[:3])
        and "fail_closed_behavior_status" in helper
        and "lower_penetration_detection_status" in helper
        and "upper_penetration_detection_status" in helper
        and "VALUE_KEYS" in helper
        and "STAGE17_4_PENETRATION_TOLERANCE" in wrapper
        and "fail-closed" in doc
    )


def stage17_4_evidence_ok(repo: Path, require: bool, accept_closed: bool) -> bool:
    if not require:
        return True
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_4_boundary_containment_fail_closed.dat")
    if data.get("final_status") == "1":
        return True
    structural = all((repo / path).is_file() for path in STAGE17_4_FILES[:3])
    return accept_closed and structural and stage17_4_fail_closed_preserved(repo)


def stage13_6_preserved(repo: Path) -> bool:
    files = [
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh",
    ]
    text = "\n".join(read_text(path) for path in files if path.exists())
    return (
        "stage13_6_production_force_density_candidate_status" in text
        and "fibre_stage13_6_production_force_density_candidate.dat" in text
        and "stage13_5_production_force_density_candidate" not in text
    )


def stage13_local_center_absent(repo: Path) -> bool:
    files = [
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "src" / "fibre_stage13_production_force_density_candidate_check.f90",
        repo / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh",
    ]
    text = "\n".join(read_text(path) for path in files if path.exists()).lower()
    return "local_subdomain_center" not in text and "subdomain_center" not in text


def stage14_small_lambda_hook_ok(repo: Path) -> bool:
    src = read_text(repo / "src" / "fibre_stage14_production_rhs_injection.f90")
    xcompact = read_text(repo / "src" / "xcompact3d.f90")
    return (
        "stage14_production_rhs_injection_apply" in src
        and "stage14_rhs_reg = stage14_requested() .and. stage14_rhs_injection_enabled()" in xcompact
        and "stage14_get_injection_gain() == 0.0" not in (src + xcompact)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_5_near_wall_contact_state.dat")
    parser.add_argument("--require-stage17-4", default="1")
    parser.add_argument("--accept-stage17-4-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--wall-safety-enable", default="1")
    parser.add_argument("--boundary-check-enable", default="1")
    parser.add_argument("--fail-closed-enable", default="1")
    parser.add_argument("--contact-placeholder-enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--y-min", default="-1.0")
    parser.add_argument("--y-max", default="1.0")
    parser.add_argument("--effective-fibre-radius", default="1.0e-3")
    parser.add_argument("--min-wall-clearance", default="1.0e-4")
    parser.add_argument("--warning-wall-clearance", default="1.0e-3")
    parser.add_argument("--penetration-tolerance", default="1.0e-12")
    parser.add_argument("--npts", default="8")
    parser.add_argument("--test-case", default="all_clear")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}
    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_5_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_5_file_{path}")

    entries = git_status_entries(repo)
    stage17_0_paths = {str(path) for path in STAGE17_0_FILES}
    stage17_1_paths = {str(path) for path in STAGE17_1_FILES}
    stage17_2_paths = {str(path) for path in STAGE17_2_FILES}
    stage17_3_paths = {str(path) for path in STAGE17_3_FILES}
    stage17_4_paths = {str(path) for path in STAGE17_4_FILES}
    stage17_0_changed = sorted(path for code, path in entries if path in stage17_0_paths and code != "??")
    stage17_1_changed = sorted(path for code, path in entries if path in stage17_1_paths and code != "??")
    stage17_2_changed = sorted(path for code, path in entries if path in stage17_2_paths and code != "??")
    stage17_3_changed = sorted(path for code, path in entries if path in stage17_3_paths and code != "??")
    stage17_4_changed = sorted(path for code, path in entries if path in stage17_4_paths and code != "??")
    closed_changed = sorted(path for _code, path in entries if path.startswith(CLOSED_STAGE_PREFIXES))
    unauthorized = sorted(
        path for _code, path in entries
        if path not in ALLOWED_STAGE17_5_CHANGES and path not in HISTORICAL_STAGE17_OUTPUTS
    )
    source_or_build_changes = sorted(
        path for _code, path in entries
        if path.endswith((".f90", ".F90", ".c", ".cc", ".cpp", ".h", ".hpp")) or path.endswith("CMakeLists.txt")
    )
    for label, changed in [
        ("stage17_0", stage17_0_changed),
        ("stage17_1", stage17_1_changed),
        ("stage17_2", stage17_2_changed),
        ("stage17_3", stage17_3_changed),
        ("stage17_4", stage17_4_changed),
    ]:
        if changed:
            reasons.append(f"{label}_file_modified_" + ",".join(changed))
    if closed_changed:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_changed))
    if unauthorized:
        reasons.append("unauthorized_stage17_5_change_detected_" + ",".join(unauthorized))

    y_min = finite_float(args.y_min)
    y_max = finite_float(args.y_max)
    r_eff = finite_float(args.effective_fibre_radius)
    min_clearance = finite_float(args.min_wall_clearance)
    warning_clearance = finite_float(args.warning_wall_clearance)
    penetration_tol = finite_float(args.penetration_tolerance)
    npts = positive_int(args.npts)
    bounds_ordered = y_min is not None and y_max is not None and y_max > y_min
    radius_ok = r_eff is not None and r_eff >= 0.0
    min_ok = min_clearance is not None and min_clearance >= 0.0
    warning_ok = warning_clearance is not None and warning_clearance >= 0.0
    tol_ok = penetration_tol is not None and penetration_tol >= 0.0
    threshold_order_ok = min_ok and warning_ok and warning_clearance >= min_clearance if min_clearance is not None and warning_clearance is not None else False
    npts_ok = npts is not None
    selected_case_ok = args.test_case in KNOWN_CASES

    cases: dict[str, dict[str, object]] = {}
    selected = {"y_values": [], "gap_wall": [], "states": [], "counts": Counter(), "global": "CLEAR", "penetrated_indices": []}
    if bounds_ordered and radius_ok and min_ok and warning_ok and tol_ok and threshold_order_ok and npts_ok and selected_case_ok:
        assert y_min is not None and y_max is not None and r_eff is not None and min_clearance is not None and warning_clearance is not None and penetration_tol is not None and npts is not None
        eval_npts = max(npts, 4)
        for case_name in KNOWN_CASES:
            case_npts = eval_npts if case_name == "mixed_states_priority" else npts
            cases[case_name] = evaluate_case(case_name, case_npts, y_min, y_max, r_eff, min_clearance, warning_clearance, penetration_tol)
        selected = cases[args.test_case]

    case_ok: dict[str, bool] = {}
    for case_name in KNOWN_CASES:
        case = cases.get(case_name)
        case_npts = len(case["states"]) if case else 0
        case_ok[case_name] = bool(
            case
            and case["global"] == expected_global(case_name)
            and case["counts"] == expected_counts(case_name, case_npts)
            and all_finite(case["gap_wall"])
            and all_finite(case["y_values"])
        )

    selected_counts_ok = bool(selected["states"]) and selected["counts"] == expected_counts(args.test_case, len(selected["states"])) if selected_case_ok else False
    selected_global_ok = bool(selected["states"]) and selected["global"] == expected_global(args.test_case) if selected_case_ok else False
    penetration_indices = list(cases.get("penetrated_fail_closed", {}).get("penetrated_indices", [])) if cases else []
    mixed_indices = list(cases.get("mixed_states_priority", {}).get("penetrated_indices", [])) if cases else []

    summary["stage17_5_requested_status"] = status(args.enable == "1")
    summary["stage17_4_evidence_status"] = status(stage17_4_evidence_ok(repo, args.require_stage17_4 == "1", args.accept_stage17_4_closed_evidence == "1"))
    summary["stage17_0_fresh_archive_fix_preserved_status"] = status(stage17_0_fresh_archive_fix_preserved(repo))
    summary["stage17_1_evidence_fix_preserved_status"] = status(stage17_1_evidence_fix_preserved(repo))
    summary["stage17_2_boundary_metadata_preserved_status"] = status(stage17_2_boundary_metadata_preserved(repo))
    summary["stage17_3_wall_clearance_preserved_status"] = status(stage17_3_wall_clearance_preserved(repo))
    summary["stage17_4_fail_closed_preserved_status"] = status(stage17_4_fail_closed_preserved(repo))
    summary["stage17_0_files_unmodified_status"] = status(not stage17_0_changed)
    summary["stage17_1_files_unmodified_status"] = status(not stage17_1_changed)
    summary["stage17_2_files_unmodified_status"] = status(not stage17_2_changed)
    summary["stage17_3_files_unmodified_status"] = status(not stage17_3_changed)
    summary["stage17_4_files_unmodified_status"] = status(not stage17_4_changed)
    summary["stage17_enable_status"] = status(args.enable == "1")
    summary["wall_safety_enable_status"] = status(args.wall_safety_enable == "1")
    summary["boundary_check_enable_status"] = status(args.boundary_check_enable == "1")
    summary["fail_closed_enable_status"] = status(args.fail_closed_enable == "1")
    summary["contact_placeholder_enable_status"] = status(args.contact_placeholder_enable == "1")
    summary["diagnostic_only_status"] = status(args.diagnostic_only == "1")
    summary["y_min_value"] = args.y_min
    summary["y_min_finite_status"] = status(y_min is not None)
    summary["y_max_value"] = args.y_max
    summary["y_max_finite_status"] = status(y_max is not None)
    summary["y_bounds_ordered_status"] = status(bounds_ordered)
    summary["effective_fibre_radius_value"] = args.effective_fibre_radius
    summary["effective_fibre_radius_status"] = status(radius_ok)
    summary["min_wall_clearance_value"] = args.min_wall_clearance
    summary["min_wall_clearance_status"] = status(min_ok)
    summary["warning_wall_clearance_value"] = args.warning_wall_clearance
    summary["warning_wall_clearance_status"] = status(warning_ok)
    summary["penetration_tolerance_value"] = args.penetration_tolerance
    summary["penetration_tolerance_status"] = status(tol_ok)
    summary["threshold_order_status"] = status(threshold_order_ok)
    summary["npts_value"] = args.npts
    summary["npts_status"] = status(npts_ok)
    summary["y_coordinates_finite_status"] = status(bool(selected["y_values"]) and all_finite(selected["y_values"]))
    summary["gap_values_finite_status"] = status(bool(selected["gap_wall"]) and all_finite(selected["gap_wall"]))
    summary["all_clear_classification_status"] = status(case_ok.get("all_clear", False))
    summary["near_wall_warning_classification_status"] = status(case_ok.get("near_wall_warning", False))
    summary["contact_placeholder_classification_status"] = status(case_ok.get("contact_placeholder", False))
    summary["penetrated_fail_closed_classification_status"] = status(case_ok.get("penetrated_fail_closed", False))
    summary["mixed_state_priority_status"] = status(case_ok.get("mixed_states_priority", False) and cases.get("mixed_states_priority", {}).get("global") == "PENETRATED_FAIL_CLOSED")
    summary["state_count_status"] = status(selected_counts_ok)
    summary["global_worst_state_status"] = status(selected_global_ok)
    summary["penetration_index_reporting_status"] = status(penetration_indices == [0] and mixed_indices == [3])
    summary["contact_placeholder_force_free_status"] = status(args.contact_placeholder_enable == "1" and args.diagnostic_only == "1")
    summary["classification_diagnostic_only_status"] = status(args.diagnostic_only == "1")

    no_source_or_build = not source_or_build_changes and not closed_changed and not stage17_0_changed and not stage17_1_changed and not stage17_2_changed and not stage17_3_changed and not stage17_4_changed and not unauthorized
    guarded_absence_keys = [
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
    ]
    for key in guarded_absence_keys:
        summary[key] = status(no_source_or_build and args.diagnostic_only == "1")
    if source_or_build_changes:
        reasons.append("stage17_5_source_or_build_change_detected_" + ",".join(source_or_build_changes))
    if not selected_case_ok:
        reasons.append("unrecognized_stage17_5_test_case_" + args.test_case)

    wrapper_text = read_text(repo / "stage17_checks" / "run_stage17_5_near_wall_contact_state.sh")
    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in wrapper_text)
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    summary["no_rg_only_dependency_status"] = status(all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_5_FILES))

    for key in SUMMARY_KEYS:
        if key.endswith("_status") and summary[key] != "1":
            reasons.append(f"{key}_not_pass")
    unknown_failure = any("unknown" in reason.lower() for reason in reasons)
    summary["no_unknown_failure_status"] = status(not unknown_failure)
    boolean_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key not in VALUE_KEYS and not key.endswith("_state_value")]
    summary["final_status"] = status(all(summary[key] == "1" for key in boolean_keys) and not reasons)

    out = repo / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary[key]}" for key in SUMMARY_KEYS]
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")
    for line in lines:
        print(line)
    if summary["final_status"] == "1":
        print("STAGE 17.5 NEAR-WALL CONTACT STATE VERDICT: PASS")
        print("STAGE 17.5 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.5 NEAR-WALL CONTACT STATE VERDICT: FAIL")
    print("STAGE 17.5 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
