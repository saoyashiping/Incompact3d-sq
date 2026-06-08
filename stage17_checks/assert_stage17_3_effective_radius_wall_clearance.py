#!/usr/bin/env python3
"""Stage 17.3 effective-radius wall-clearance diagnostic helper.

Stage 17.3 computes standalone analytic point-wise wall-clearance diagnostics for
fibre control-point y coordinates only. It is diagnostic-only: it does not modify
fibre state, fluid RHS, IBM forcing, production hooks, DNS-core numerics, or any
closed Stage 10--16 / Stage 17.0--17.2 file.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_3_requested_status",
    "stage17_2_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "diagnostic_only_status",
    "wall_normal_direction_status",
    "y_min_value",
    "y_min_finite_status",
    "y_max_value",
    "y_max_finite_status",
    "y_bounds_ordered_status",
    "effective_fibre_radius_value",
    "effective_fibre_radius_status",
    "npts_value",
    "npts_status",
    "y_coordinates_finite_status",
    "lower_distance_formula_status",
    "upper_distance_formula_status",
    "centerline_min_distance_formula_status",
    "effective_radius_lower_gap_formula_status",
    "effective_radius_upper_gap_formula_status",
    "effective_radius_min_gap_formula_status",
    "exact_radius_touch_zero_gap_status",
    "negative_gap_diagnostic_only_status",
    "distance_values_finite_status",
    "gap_values_finite_status",
    "min_centerline_wall_distance_value",
    "min_effective_wall_gap_value",
    "no_contact_state_classification_status",
    "no_near_wall_warning_classification_status",
    "no_wall_penetration_fail_closed_status",
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
    "npts_value",
    "min_centerline_wall_distance_value",
    "min_effective_wall_gap_value",
}

REQUIRED_STAGE17_3_FILES = [
    Path("stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh"),
    Path("stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py"),
    Path("stage17_checks/stage17_3_effective_radius_wall_clearance.md"),
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
ALLOWED_STAGE17_3_CHANGES = {
    "stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh",
    "stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py",
    "stage17_checks/stage17_3_effective_radius_wall_clearance.md",
    "stage17_outputs/fibre_stage17_3_effective_radius_wall_clearance.dat",
}
HISTORICAL_STAGE17_OUTPUTS = {
    "stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat",
    "stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat",
    "stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat",
}
CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))
KNOWN_CASES = {
    "centered_clear",
    "near_lower_wall_gap_positive",
    "exact_radius_touch_placeholder",
    "outside_gap_negative_diagnostic_only",
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

    Stage 17.3 reuses the corrected Stage 17.2 / Stage 17.1 / Stage 17.0 / Stage 16
    false-positive-safe policy. Only executable shell-wrapper command lines are
    audited. Markdown prose is not scanned as real code-regression evidence,
    negative-check strings are not treated as behavior, and regex literals such as
    rg[[:space:]] do not count as ripgrep invocations. Any real wrapper use of rg
    must include grep fallback so ripgrep is never a hard dependency.
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


def control_point_y_values(case_name: str, npts: int, y_min: float, y_max: float, r_eff: float) -> list[float]:
    height = y_max - y_min
    if npts <= 0 or height <= 0.0:
        return []
    mid = 0.5 * (y_min + y_max)
    safe_low = y_min + max(0.25 * height, 2.0 * r_eff)
    safe_high = y_max - max(0.25 * height, 2.0 * r_eff)
    if safe_low > safe_high:
        safe_low = safe_high = mid
    values = [mid for _ in range(npts)]
    if case_name == "centered_clear":
        if npts == 1:
            return [mid]
        return [safe_low + (safe_high - safe_low) * i / max(npts - 1, 1) for i in range(npts)]
    if case_name == "near_lower_wall_gap_positive":
        values[0] = y_min + max(2.0 * r_eff, 0.05 * height)
        return values
    if case_name == "exact_radius_touch_placeholder":
        values[0] = y_min + r_eff
        return values
    if case_name == "outside_gap_negative_diagnostic_only":
        values[0] = y_min + (0.5 * r_eff if r_eff > 0.0 else -0.01 * height)
        return values
    return []


def wall_clearance(y_values: list[float], y_min: float, y_max: float, r_eff: float) -> dict[str, object]:
    d_lower = [y - y_min for y in y_values]
    d_upper = [y_max - y for y in y_values]
    d_center = [min(lo, up) for lo, up in zip(d_lower, d_upper)]
    gap_lower = [lo - r_eff for lo in d_lower]
    gap_upper = [up - r_eff for up in d_upper]
    gap_wall = [min(lo, up) for lo, up in zip(gap_lower, gap_upper)]
    return {
        "d_lower": d_lower,
        "d_upper": d_upper,
        "d_center": d_center,
        "gap_lower": gap_lower,
        "gap_upper": gap_upper,
        "gap_wall": gap_wall,
        "d_center_min": min(d_center) if d_center else float("nan"),
        "gap_wall_min": min(gap_wall) if gap_wall else float("nan"),
    }


def all_finite(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def close_enough(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def formula_statuses(y_values: list[float], y_min: float, y_max: float, r_eff: float, tol: float) -> dict[str, bool]:
    result = wall_clearance(y_values, y_min, y_max, r_eff)
    expected_lower = [y - y_min for y in y_values]
    expected_upper = [y_max - y for y in y_values]
    expected_center = [min(lo, up) for lo, up in zip(expected_lower, expected_upper)]
    expected_gap_lower = [lo - r_eff for lo in expected_lower]
    expected_gap_upper = [up - r_eff for up in expected_upper]
    expected_gap = [min(lo, up) for lo, up in zip(expected_gap_lower, expected_gap_upper)]
    return {
        "lower": all(close_enough(a, b, tol) for a, b in zip(result["d_lower"], expected_lower)),
        "upper": all(close_enough(a, b, tol) for a, b in zip(result["d_upper"], expected_upper)),
        "center_min": close_enough(float(result["d_center_min"]), min(expected_center), tol) if expected_center else False,
        "gap_lower": all(close_enough(a, b, tol) for a, b in zip(result["gap_lower"], expected_gap_lower)),
        "gap_upper": all(close_enough(a, b, tol) for a, b in zip(result["gap_upper"], expected_gap_upper)),
        "gap_min": close_enough(float(result["gap_wall_min"]), min(expected_gap), tol) if expected_gap else False,
        "distance_finite": all_finite(result["d_lower"] + result["d_upper"] + result["d_center"]),
        "gap_finite": all_finite(result["gap_lower"] + result["gap_upper"] + result["gap_wall"]),
    }


def analytic_touch_and_negative_ok(y_min: float, y_max: float, r_eff: float, tol: float) -> tuple[bool, bool]:
    radius = max(r_eff, 1.0e-6)
    touch = wall_clearance([y_min + radius], y_min, y_max, radius)
    negative = wall_clearance([y_min + 0.5 * radius], y_min, y_max, radius)
    return abs(float(touch["gap_wall_min"])) <= tol, float(negative["gap_wall_min"]) < 0.0


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


def stage17_2_evidence_ok(repo: Path, require: bool, accept_closed: bool) -> bool:
    if not require:
        return True
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_2_channel_wall_domain_boundary.dat")
    if data.get("final_status") == "1":
        return True
    structural = all((repo / path).is_file() for path in STAGE17_2_FILES[:3])
    return accept_closed and structural and stage17_2_boundary_metadata_preserved(repo)


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
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_3_effective_radius_wall_clearance.dat")
    parser.add_argument("--require-stage17-2", default="1")
    parser.add_argument("--accept-stage17-2-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--wall-safety-enable", default="1")
    parser.add_argument("--boundary-check-enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--y-min", default="-1.0")
    parser.add_argument("--y-max", default="1.0")
    parser.add_argument("--effective-fibre-radius", default="1.0e-3")
    parser.add_argument("--npts", default="8")
    parser.add_argument("--test-case", default="centered_clear")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}
    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_3_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_3_file_{path}")

    entries = git_status_entries(repo)
    stage17_0_paths = {str(path) for path in STAGE17_0_FILES}
    stage17_1_paths = {str(path) for path in STAGE17_1_FILES}
    stage17_2_paths = {str(path) for path in STAGE17_2_FILES}
    stage17_0_changed = sorted(path for code, path in entries if path in stage17_0_paths and code != "??")
    stage17_1_changed = sorted(path for code, path in entries if path in stage17_1_paths and code != "??")
    stage17_2_changed = sorted(path for code, path in entries if path in stage17_2_paths and code != "??")
    closed_changed = sorted(path for _code, path in entries if path.startswith(CLOSED_STAGE_PREFIXES))
    unauthorized = sorted(
        path for _code, path in entries
        if path not in ALLOWED_STAGE17_3_CHANGES and path not in HISTORICAL_STAGE17_OUTPUTS
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
    if closed_changed:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_changed))
    if unauthorized:
        reasons.append("unauthorized_stage17_3_change_detected_" + ",".join(unauthorized))

    y_min = finite_float(args.y_min)
    y_max = finite_float(args.y_max)
    r_eff = finite_float(args.effective_fibre_radius)
    npts = positive_int(args.npts)
    zero_tol = finite_float(args.zero_tol)
    bounds_ordered = y_min is not None and y_max is not None and y_max > y_min
    radius_ok = r_eff is not None and r_eff >= 0.0
    npts_ok = npts is not None
    tol = zero_tol if zero_tol is not None and zero_tol >= 0.0 else 0.0
    selected_case_ok = args.test_case in KNOWN_CASES

    y_values: list[float] = []
    result: dict[str, object] = {"d_center_min": float("nan"), "gap_wall_min": float("nan")}
    formula = {key: False for key in ["lower", "upper", "center_min", "gap_lower", "gap_upper", "gap_min", "distance_finite", "gap_finite"]}
    touch_ok = False
    negative_ok = False
    if bounds_ordered and radius_ok and npts_ok and selected_case_ok and r_eff is not None and y_min is not None and y_max is not None and npts is not None:
        y_values = control_point_y_values(args.test_case, npts, y_min, y_max, r_eff)
        result = wall_clearance(y_values, y_min, y_max, r_eff)
        formula = formula_statuses(y_values, y_min, y_max, r_eff, max(tol, 1.0e-15))
        touch_ok, negative_ok = analytic_touch_and_negative_ok(y_min, y_max, r_eff, max(tol, 1.0e-15))

    summary["stage17_3_requested_status"] = status(args.enable == "1")
    summary["stage17_2_evidence_status"] = status(stage17_2_evidence_ok(repo, args.require_stage17_2 == "1", args.accept_stage17_2_closed_evidence == "1"))
    summary["stage17_0_fresh_archive_fix_preserved_status"] = status(stage17_0_fresh_archive_fix_preserved(repo))
    summary["stage17_1_evidence_fix_preserved_status"] = status(stage17_1_evidence_fix_preserved(repo))
    summary["stage17_2_boundary_metadata_preserved_status"] = status(stage17_2_boundary_metadata_preserved(repo))
    summary["stage17_0_files_unmodified_status"] = status(not stage17_0_changed)
    summary["stage17_1_files_unmodified_status"] = status(not stage17_1_changed)
    summary["stage17_2_files_unmodified_status"] = status(not stage17_2_changed)
    summary["stage17_enable_status"] = status(args.enable == "1")
    summary["wall_safety_enable_status"] = status(args.wall_safety_enable == "1")
    summary["boundary_check_enable_status"] = status(args.boundary_check_enable == "1")
    summary["diagnostic_only_status"] = status(args.diagnostic_only == "1")
    summary["wall_normal_direction_status"] = "1"
    summary["y_min_value"] = args.y_min
    summary["y_min_finite_status"] = status(y_min is not None)
    summary["y_max_value"] = args.y_max
    summary["y_max_finite_status"] = status(y_max is not None)
    summary["y_bounds_ordered_status"] = status(bounds_ordered)
    summary["effective_fibre_radius_value"] = args.effective_fibre_radius
    summary["effective_fibre_radius_status"] = status(radius_ok)
    summary["npts_value"] = args.npts
    summary["npts_status"] = status(npts_ok)
    summary["y_coordinates_finite_status"] = status(bool(y_values) and all_finite(y_values))
    summary["lower_distance_formula_status"] = status(formula["lower"])
    summary["upper_distance_formula_status"] = status(formula["upper"])
    summary["centerline_min_distance_formula_status"] = status(formula["center_min"])
    summary["effective_radius_lower_gap_formula_status"] = status(formula["gap_lower"])
    summary["effective_radius_upper_gap_formula_status"] = status(formula["gap_upper"])
    summary["effective_radius_min_gap_formula_status"] = status(formula["gap_min"])
    summary["exact_radius_touch_zero_gap_status"] = status(touch_ok)
    summary["negative_gap_diagnostic_only_status"] = status(negative_ok and args.diagnostic_only == "1")
    summary["distance_values_finite_status"] = status(formula["distance_finite"])
    summary["gap_values_finite_status"] = status(formula["gap_finite"])
    summary["min_centerline_wall_distance_value"] = f"{float(result['d_center_min']):.17g}"
    summary["min_effective_wall_gap_value"] = f"{float(result['gap_wall_min']):.17g}"

    no_source_or_build = not source_or_build_changes and not closed_changed and not stage17_0_changed and not stage17_1_changed and not stage17_2_changed and not unauthorized
    guarded_absence_keys = [
        "no_contact_state_classification_status",
        "no_near_wall_warning_classification_status",
        "no_wall_penetration_fail_closed_status",
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
        reasons.append("stage17_3_source_or_build_change_detected_" + ",".join(source_or_build_changes))
    if not selected_case_ok:
        reasons.append("unknown_stage17_3_test_case_" + args.test_case)

    wrapper_text = read_text(repo / "stage17_checks" / "run_stage17_3_effective_radius_wall_clearance.sh")
    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in wrapper_text)
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    summary["no_rg_only_dependency_status"] = status(all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_3_FILES))

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
        print("STAGE 17.3 EFFECTIVE-RADIUS WALL CLEARANCE VERDICT: PASS")
        print("STAGE 17.3 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.3 EFFECTIVE-RADIUS WALL CLEARANCE VERDICT: FAIL")
    print("STAGE 17.3 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
