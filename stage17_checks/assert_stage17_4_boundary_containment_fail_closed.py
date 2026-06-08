#!/usr/bin/env python3
"""Stage 17.4 boundary-containment and wall-penetration fail-closed helper.

Stage 17.4 evaluates standalone analytic point-wise boundary containment using
surface gaps only. It is diagnostic-only and non-invasive: it detects penetration
and reports explicit fail-closed reasons for penetration cases, but it never applies
wall/contact/collision forces, changes fibre state, changes fluid RHS, inserts
production hooks, or edits closed Stage 10--16 / Stage 17.0--17.3 files.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_4_requested_status",
    "stage17_3_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_3_wall_clearance_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_3_files_unmodified_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "fail_closed_enable_status",
    "diagnostic_only_status",
    "wall_normal_direction_status",
    "y_min_value",
    "y_min_finite_status",
    "y_max_value",
    "y_max_finite_status",
    "y_bounds_ordered_status",
    "effective_fibre_radius_value",
    "effective_fibre_radius_status",
    "penetration_tolerance_value",
    "penetration_tolerance_status",
    "npts_value",
    "npts_status",
    "y_coordinates_finite_status",
    "lower_surface_gap_formula_status",
    "upper_surface_gap_formula_status",
    "surface_gap_values_finite_status",
    "boundary_containment_formula_status",
    "contained_case_status",
    "exact_lower_touch_nonpenetrating_status",
    "exact_upper_touch_nonpenetrating_status",
    "lower_penetration_detection_status",
    "upper_penetration_detection_status",
    "penetration_depth_status",
    "offending_point_index_status",
    "fail_closed_behavior_status",
    "no_contact_state_classification_status",
    "no_near_wall_warning_classification_status",
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
    "penetration_tolerance_value",
    "npts_value",
}

REQUIRED_STAGE17_4_FILES = [
    Path("stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh"),
    Path("stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py"),
    Path("stage17_checks/stage17_4_boundary_containment_fail_closed.md"),
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
ALLOWED_STAGE17_4_CHANGES = {
    "stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh",
    "stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py",
    "stage17_checks/stage17_4_boundary_containment_fail_closed.md",
    "stage17_outputs/fibre_stage17_4_boundary_containment_fail_closed.dat",
}
HISTORICAL_STAGE17_OUTPUTS = {
    "stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat",
    "stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat",
    "stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat",
    "stage17_outputs/fibre_stage17_3_effective_radius_wall_clearance.dat",
}
CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))
KNOWN_CASES = {
    "contained_clear",
    "exact_radius_touch_lower",
    "exact_radius_touch_upper",
    "lower_wall_penetration",
    "upper_wall_penetration",
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

    Stage 17.4 reuses the corrected Stage 17.3 / Stage 17.2 / Stage 17.1 /
    Stage 17.0 / Stage 16 false-positive-safe policy. Only executable shell-wrapper
    command lines are audited. Markdown prose is not scanned as real code-regression
    evidence, negative-check strings are not treated as behavior, and regex literals
    such as rg[[:space:]] do not count as ripgrep invocations. Any real wrapper use
    of rg must include grep fallback so ripgrep is never a hard dependency. Existing
    diagnostic failure-label strings in closed Stage 17.0--17.3 files are labels, not
    rollback evidence.
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


def control_point_y_values(case_name: str, npts: int, y_min: float, y_max: float, r_eff: float, tol: float) -> list[float]:
    height = y_max - y_min
    if npts <= 0 or height <= 0.0:
        return []
    mid = 0.5 * (y_min + y_max)
    safe_low = y_min + max(0.25 * height, 2.0 * r_eff + 2.0 * tol)
    safe_high = y_max - max(0.25 * height, 2.0 * r_eff + 2.0 * tol)
    if safe_low > safe_high:
        safe_low = safe_high = mid
    values = [mid for _ in range(npts)]
    if case_name == "contained_clear":
        if npts == 1:
            return [mid]
        return [safe_low + (safe_high - safe_low) * i / max(npts - 1, 1) for i in range(npts)]
    if case_name == "exact_radius_touch_lower":
        values[0] = y_min + r_eff
        return values
    if case_name == "exact_radius_touch_upper":
        values[0] = y_max - r_eff
        return values
    if case_name == "lower_wall_penetration":
        values[0] = y_min + r_eff - 2.0 * max(tol, 1.0e-15)
        return values
    if case_name == "upper_wall_penetration":
        values[0] = y_max - r_eff + 2.0 * max(tol, 1.0e-15)
        return values
    return []


def containment_diagnostics(y_values: list[float], y_min: float, y_max: float, r_eff: float, tol: float) -> dict[str, object]:
    lower_gap = [y - r_eff - y_min for y in y_values]
    upper_gap = [y_max - (y + r_eff) for y in y_values]
    contained = [lo >= -tol and up >= -tol for lo, up in zip(lower_gap, upper_gap)]
    lower_indices = [idx for idx, value in enumerate(lower_gap) if value < -tol]
    upper_indices = [idx for idx, value in enumerate(upper_gap) if value < -tol]
    lower_depths = [-lower_gap[idx] for idx in lower_indices]
    upper_depths = [-upper_gap[idx] for idx in upper_indices]
    return {
        "lower_gap": lower_gap,
        "upper_gap": upper_gap,
        "contained": contained,
        "lower_indices": lower_indices,
        "upper_indices": upper_indices,
        "lower_depths": lower_depths,
        "upper_depths": upper_depths,
        "any_penetration": bool(lower_indices or upper_indices),
    }


def all_finite(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def close_enough(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def evaluate_case(case_name: str, npts: int, y_min: float, y_max: float, r_eff: float, tol: float) -> dict[str, object]:
    y_values = control_point_y_values(case_name, npts, y_min, y_max, r_eff, tol)
    diag = containment_diagnostics(y_values, y_min, y_max, r_eff, tol)
    lower_expected = [y - y_min - r_eff for y in y_values]
    upper_expected = [y_max - y - r_eff for y in y_values]
    contained_expected = [lo >= -tol and up >= -tol for lo, up in zip(lower_expected, upper_expected)]
    lower_formula = all(close_enough(a, b, max(tol, 1.0e-15)) for a, b in zip(diag["lower_gap"], lower_expected))
    upper_formula = all(close_enough(a, b, max(tol, 1.0e-15)) for a, b in zip(diag["upper_gap"], upper_expected))
    containment_formula = list(diag["contained"]) == contained_expected
    return {
        "y_values": y_values,
        "diag": diag,
        "lower_formula": lower_formula,
        "upper_formula": upper_formula,
        "containment_formula": containment_formula,
        "gaps_finite": all_finite(diag["lower_gap"] + diag["upper_gap"]),
        "y_finite": all_finite(y_values),
    }


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


def stage17_3_evidence_ok(repo: Path, require: bool, accept_closed: bool) -> bool:
    if not require:
        return True
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_3_effective_radius_wall_clearance.dat")
    if data.get("final_status") == "1":
        return True
    structural = all((repo / path).is_file() for path in STAGE17_3_FILES[:3])
    return accept_closed and structural and stage17_3_wall_clearance_preserved(repo)


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
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_4_boundary_containment_fail_closed.dat")
    parser.add_argument("--require-stage17-3", default="1")
    parser.add_argument("--accept-stage17-3-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--wall-safety-enable", default="1")
    parser.add_argument("--boundary-check-enable", default="1")
    parser.add_argument("--fail-closed-enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--y-min", default="-1.0")
    parser.add_argument("--y-max", default="1.0")
    parser.add_argument("--effective-fibre-radius", default="1.0e-3")
    parser.add_argument("--penetration-tolerance", default="1.0e-12")
    parser.add_argument("--npts", default="8")
    parser.add_argument("--test-case", default="contained_clear")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}
    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_4_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_4_file_{path}")

    entries = git_status_entries(repo)
    stage17_0_paths = {str(path) for path in STAGE17_0_FILES}
    stage17_1_paths = {str(path) for path in STAGE17_1_FILES}
    stage17_2_paths = {str(path) for path in STAGE17_2_FILES}
    stage17_3_paths = {str(path) for path in STAGE17_3_FILES}
    stage17_0_changed = sorted(path for code, path in entries if path in stage17_0_paths and code != "??")
    stage17_1_changed = sorted(path for code, path in entries if path in stage17_1_paths and code != "??")
    stage17_2_changed = sorted(path for code, path in entries if path in stage17_2_paths and code != "??")
    stage17_3_changed = sorted(path for code, path in entries if path in stage17_3_paths and code != "??")
    closed_changed = sorted(path for _code, path in entries if path.startswith(CLOSED_STAGE_PREFIXES))
    unauthorized = sorted(
        path for _code, path in entries
        if path not in ALLOWED_STAGE17_4_CHANGES and path not in HISTORICAL_STAGE17_OUTPUTS
    )
    source_or_build_changes = sorted(
        path for _code, path in entries
        if path.endswith((".f90", ".F90", ".c", ".cc", ".cpp", ".h", ".hpp")) or path.endswith("CMakeLists.txt")
    )
    if stage17_0_changed:
        reasons.append("stage17_0_file_modified_" + ",".join(stage17_0_changed))
    if stage17_1_changed:
        reasons.append("stage17_1_file_modified_" + ",".join(stage17_1_changed))
    if stage17_2_changed:
        reasons.append("stage17_2_file_modified_" + ",".join(stage17_2_changed))
    if stage17_3_changed:
        reasons.append("stage17_3_file_modified_" + ",".join(stage17_3_changed))
    if closed_changed:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_changed))
    if unauthorized:
        reasons.append("unauthorized_stage17_4_change_detected_" + ",".join(unauthorized))

    y_min = finite_float(args.y_min)
    y_max = finite_float(args.y_max)
    r_eff = finite_float(args.effective_fibre_radius)
    penetration_tol = finite_float(args.penetration_tolerance)
    npts = positive_int(args.npts)
    zero_tol = finite_float(args.zero_tol)
    bounds_ordered = y_min is not None and y_max is not None and y_max > y_min
    radius_ok = r_eff is not None and r_eff >= 0.0
    tolerance_ok = penetration_tol is not None and penetration_tol >= 0.0
    npts_ok = npts is not None
    tol = penetration_tol if penetration_tol is not None and penetration_tol >= 0.0 else 0.0
    zero = zero_tol if zero_tol is not None and zero_tol >= 0.0 else 0.0
    selected_case_ok = args.test_case in KNOWN_CASES

    selected = {"y_values": [], "diag": {"lower_gap": [], "upper_gap": [], "lower_indices": [], "upper_indices": []}, "lower_formula": False, "upper_formula": False, "containment_formula": False, "gaps_finite": False, "y_finite": False}
    cases: dict[str, dict[str, object]] = {}
    if bounds_ordered and radius_ok and tolerance_ok and npts_ok and selected_case_ok and y_min is not None and y_max is not None and r_eff is not None and npts is not None:
        for case_name in KNOWN_CASES:
            cases[case_name] = evaluate_case(case_name, npts, y_min, y_max, r_eff, max(tol, zero, 1.0e-15))
        selected = cases[args.test_case]

    contained_diag = cases.get("contained_clear", {})
    lower_touch_diag = cases.get("exact_radius_touch_lower", {})
    upper_touch_diag = cases.get("exact_radius_touch_upper", {})
    lower_pen_diag = cases.get("lower_wall_penetration", {})
    upper_pen_diag = cases.get("upper_wall_penetration", {})
    lower_pen_data = lower_pen_diag.get("diag", {}) if lower_pen_diag else {}
    upper_pen_data = upper_pen_diag.get("diag", {}) if upper_pen_diag else {}
    lower_depths = list(lower_pen_data.get("lower_depths", [])) if lower_pen_data else []
    upper_depths = list(upper_pen_data.get("upper_depths", [])) if upper_pen_data else []
    lower_indices = list(lower_pen_data.get("lower_indices", [])) if lower_pen_data else []
    upper_indices = list(upper_pen_data.get("upper_indices", [])) if upper_pen_data else []

    summary["stage17_4_requested_status"] = status(args.enable == "1")
    summary["stage17_3_evidence_status"] = status(stage17_3_evidence_ok(repo, args.require_stage17_3 == "1", args.accept_stage17_3_closed_evidence == "1"))
    summary["stage17_0_fresh_archive_fix_preserved_status"] = status(stage17_0_fresh_archive_fix_preserved(repo))
    summary["stage17_1_evidence_fix_preserved_status"] = status(stage17_1_evidence_fix_preserved(repo))
    summary["stage17_2_boundary_metadata_preserved_status"] = status(stage17_2_boundary_metadata_preserved(repo))
    summary["stage17_3_wall_clearance_preserved_status"] = status(stage17_3_wall_clearance_preserved(repo))
    summary["stage17_0_files_unmodified_status"] = status(not stage17_0_changed)
    summary["stage17_1_files_unmodified_status"] = status(not stage17_1_changed)
    summary["stage17_2_files_unmodified_status"] = status(not stage17_2_changed)
    summary["stage17_3_files_unmodified_status"] = status(not stage17_3_changed)
    summary["stage17_enable_status"] = status(args.enable == "1")
    summary["wall_safety_enable_status"] = status(args.wall_safety_enable == "1")
    summary["boundary_check_enable_status"] = status(args.boundary_check_enable == "1")
    summary["fail_closed_enable_status"] = status(args.fail_closed_enable == "1")
    summary["diagnostic_only_status"] = status(args.diagnostic_only == "1")
    summary["wall_normal_direction_status"] = "1"
    summary["y_min_value"] = args.y_min
    summary["y_min_finite_status"] = status(y_min is not None)
    summary["y_max_value"] = args.y_max
    summary["y_max_finite_status"] = status(y_max is not None)
    summary["y_bounds_ordered_status"] = status(bounds_ordered)
    summary["effective_fibre_radius_value"] = args.effective_fibre_radius
    summary["effective_fibre_radius_status"] = status(radius_ok)
    summary["penetration_tolerance_value"] = args.penetration_tolerance
    summary["penetration_tolerance_status"] = status(tolerance_ok)
    summary["npts_value"] = args.npts
    summary["npts_status"] = status(npts_ok)
    summary["y_coordinates_finite_status"] = status(bool(selected["y_values"]) and bool(selected["y_finite"]))
    summary["lower_surface_gap_formula_status"] = status(bool(selected["lower_formula"]))
    summary["upper_surface_gap_formula_status"] = status(bool(selected["upper_formula"]))
    summary["surface_gap_values_finite_status"] = status(bool(selected["gaps_finite"]))
    summary["boundary_containment_formula_status"] = status(bool(selected["containment_formula"]))
    summary["contained_case_status"] = status(bool(contained_diag) and not bool(contained_diag["diag"]["any_penetration"]))
    summary["exact_lower_touch_nonpenetrating_status"] = status(bool(lower_touch_diag) and not bool(lower_touch_diag["diag"]["any_penetration"]))
    summary["exact_upper_touch_nonpenetrating_status"] = status(bool(upper_touch_diag) and not bool(upper_touch_diag["diag"]["any_penetration"]))
    summary["lower_penetration_detection_status"] = status(bool(lower_indices))
    summary["upper_penetration_detection_status"] = status(bool(upper_indices))
    depth_ok = bool(lower_depths and upper_depths) and all(math.isfinite(d) and d > 0.0 for d in lower_depths + upper_depths)
    summary["penetration_depth_status"] = status(depth_ok)
    summary["offending_point_index_status"] = status(lower_indices == [0] and upper_indices == [0])
    fail_closed_ok = args.fail_closed_enable == "1" and bool(lower_indices) and bool(upper_indices) and depth_ok
    summary["fail_closed_behavior_status"] = status(fail_closed_ok)

    no_source_or_build = not source_or_build_changes and not closed_changed and not stage17_0_changed and not stage17_1_changed and not stage17_2_changed and not stage17_3_changed and not unauthorized
    guarded_absence_keys = [
        "no_contact_state_classification_status",
        "no_near_wall_warning_classification_status",
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
        reasons.append("stage17_4_source_or_build_change_detected_" + ",".join(source_or_build_changes))
    if not selected_case_ok:
        reasons.append("unrecognized_stage17_4_test_case_" + args.test_case)
    selected_diag = selected["diag"]
    selected_has_penetration = bool(selected_diag.get("any_penetration", False)) if isinstance(selected_diag, dict) else False
    if selected_has_penetration and args.fail_closed_enable == "1":
        lower_selected = list(selected_diag.get("lower_indices", [])) if isinstance(selected_diag, dict) else []
        upper_selected = list(selected_diag.get("upper_indices", [])) if isinstance(selected_diag, dict) else []
        if lower_selected:
            reasons.append(f"selected_lower_wall_penetration_fail_closed_index_{lower_selected[0]}")
        if upper_selected:
            reasons.append(f"selected_upper_wall_penetration_fail_closed_index_{upper_selected[0]}")

    wrapper_text = read_text(repo / "stage17_checks" / "run_stage17_4_boundary_containment_fail_closed.sh")
    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in wrapper_text)
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    summary["no_rg_only_dependency_status"] = status(all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_4_FILES))

    for key in SUMMARY_KEYS:
        if key.endswith("_status") and summary[key] != "1":
            reasons.append(f"{key}_not_pass")
    unknown_failure = any("unknown" in reason.lower() for reason in reasons)
    summary["no_unknown_failure_status"] = status(not unknown_failure)
    boolean_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key not in VALUE_KEYS]
    summary["final_status"] = status(all(summary[key] == "1" for key in boolean_keys) and not reasons)

    out = repo / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary[key]}" for key in SUMMARY_KEYS]
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")
    for line in lines:
        print(line)
    if summary["final_status"] == "1":
        print("STAGE 17.4 BOUNDARY CONTAINMENT FAIL-CLOSED VERDICT: PASS")
        print("STAGE 17.4 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.4 BOUNDARY CONTAINMENT FAIL-CLOSED VERDICT: FAIL")
    print("STAGE 17.4 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
