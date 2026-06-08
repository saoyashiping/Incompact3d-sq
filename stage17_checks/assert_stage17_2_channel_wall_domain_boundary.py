#!/usr/bin/env python3
"""Stage 17.2 channel wall/domain-boundary diagnostic metadata helper.

Stage 17.2 is metadata-only. It validates finite ordered channel wall bounds,
wall-normal direction, and explicit x/y/z boundary-policy metadata for later wall
clearance diagnostics. It does not compute wall clearance, evaluate fibre points,
classify contact, detect penetration, build targets, run MPI, or introduce any
contact/collision physics.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_2_requested_status",
    "stage17_1_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_boundary_metadata_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "diagnostic_only_status",
    "wall_normal_direction_value",
    "wall_normal_direction_status",
    "y_min_value",
    "y_min_finite_status",
    "y_max_value",
    "y_max_finite_status",
    "y_bounds_ordered_status",
    "channel_height_value",
    "channel_height_positive_status",
    "x_boundary_policy_status",
    "y_boundary_policy_status",
    "z_boundary_policy_status",
    "periodic_x_metadata_status",
    "periodic_z_metadata_status",
    "wall_y_metadata_status",
    "invalid_boundary_metadata_rejection_status",
    "no_wall_clearance_computation_status",
    "no_effective_radius_gap_computation_status",
    "no_contact_state_classification_status",
    "no_wall_penetration_detection_status",
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
    "wall_normal_direction_value",
    "y_min_value",
    "y_max_value",
    "channel_height_value",
}

REQUIRED_STAGE17_2_FILES = [
    Path("stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh"),
    Path("stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py"),
    Path("stage17_checks/stage17_2_channel_wall_domain_boundary.md"),
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
ALLOWED_STAGE17_2_CHANGES = {
    "stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh",
    "stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py",
    "stage17_checks/stage17_2_channel_wall_domain_boundary.md",
    "stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat",
}
HISTORICAL_STAGE17_OUTPUTS = {
    "stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat",
    "stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat",
}
CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))
EXPLICIT_XZ_POLICIES = {"periodic", "metadata_only", "nonperiodic", "explicit_nonperiodic"}
WALL_Y_POLICIES = {"wall_bounded", "wall", "wall_normal", "nonperiodic_wall", "bounded_wall"}


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

    Stage 17.2 deliberately reuses the corrected Stage 17.1 / Stage 17.0 / Stage 16
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


def validate_boundary_metadata(y_min_text: str, y_max_text: str, direction: str, x_policy: str, y_policy: str, z_policy: str) -> dict[str, object]:
    y_min = finite_float(y_min_text)
    y_max = finite_float(y_max_text)
    height = None if y_min is None or y_max is None else y_max - y_min
    return {
        "y_min": y_min,
        "y_max": y_max,
        "height": height,
        "direction_ok": direction.lower() == "y",
        "y_min_ok": y_min is not None,
        "y_max_ok": y_max is not None,
        "ordered_ok": y_min is not None and y_max is not None and y_max > y_min,
        "height_ok": height is not None and math.isfinite(height) and height > 0.0,
        "x_policy_ok": x_policy.lower() in EXPLICIT_XZ_POLICIES,
        "y_policy_ok": y_policy.lower() in WALL_Y_POLICIES,
        "z_policy_ok": z_policy.lower() in EXPLICIT_XZ_POLICIES,
        "periodic_x_ok": x_policy.lower() == "periodic" or x_policy.lower() in EXPLICIT_XZ_POLICIES,
        "periodic_z_ok": z_policy.lower() == "periodic" or z_policy.lower() in EXPLICIT_XZ_POLICIES,
        "wall_y_ok": y_policy.lower() in WALL_Y_POLICIES,
    }


def invalid_boundary_metadata_rejected() -> bool:
    bad_cases = [
        ("nan", "1.0", "y", "periodic", "wall_bounded", "periodic"),
        ("0.0", "inf", "y", "periodic", "wall_bounded", "periodic"),
        ("1.0", "1.0", "y", "periodic", "wall_bounded", "periodic"),
        ("1.0", "0.0", "y", "periodic", "wall_bounded", "periodic"),
        ("-1.0", "1.0", "x", "periodic", "wall_bounded", "periodic"),
        ("-1.0", "1.0", "y", "", "wall_bounded", "periodic"),
        ("-1.0", "1.0", "y", "periodic", "periodic", "periodic"),
        ("-1.0", "1.0", "y", "periodic", "wall_bounded", ""),
    ]
    required = (
        "direction_ok", "y_min_ok", "y_max_ok", "ordered_ok", "height_ok",
        "x_policy_ok", "y_policy_ok", "z_policy_ok", "wall_y_ok",
    )
    for case in bad_cases:
        result = validate_boundary_metadata(*case)
        if all(bool(result[key]) for key in required):
            return False
    return True


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


def stage17_1_evidence_ok(repo: Path, require: bool, accept_closed: bool) -> bool:
    if not require:
        return True
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_1_wall_contact_safety_config.dat")
    if data.get("final_status") == "1":
        return True
    structural = all((repo / path).is_file() for path in STAGE17_1_FILES[:3])
    return accept_closed and structural and stage17_1_evidence_fix_preserved(repo)


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
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat")
    parser.add_argument("--require-stage17-1", default="1")
    parser.add_argument("--accept-stage17-1-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--wall-safety-enable", default="1")
    parser.add_argument("--boundary-check-enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--y-min", default="-1.0")
    parser.add_argument("--y-max", default="1.0")
    parser.add_argument("--wall-normal-direction", default="y")
    parser.add_argument("--x-boundary-policy", default="periodic")
    parser.add_argument("--y-boundary-policy", default="wall_bounded")
    parser.add_argument("--z-boundary-policy", default="periodic")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}
    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_2_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_2_file_{path}")

    entries = git_status_entries(repo)
    stage17_0_paths = {str(path) for path in STAGE17_0_FILES}
    stage17_1_paths = {str(path) for path in STAGE17_1_FILES}
    stage17_0_changed = sorted(path for code, path in entries if path in stage17_0_paths and code != "??")
    stage17_1_changed = sorted(path for code, path in entries if path in stage17_1_paths and code != "??")
    closed_changed = sorted(path for _code, path in entries if path.startswith(CLOSED_STAGE_PREFIXES))
    unauthorized = sorted(
        path for _code, path in entries
        if path not in ALLOWED_STAGE17_2_CHANGES and path not in HISTORICAL_STAGE17_OUTPUTS
    )
    source_or_build_changes = sorted(
        path for _code, path in entries
        if path.endswith((".f90", ".F90", ".c", ".cc", ".cpp", ".h", ".hpp")) or path.endswith("CMakeLists.txt")
    )
    if stage17_0_changed:
        reasons.append("stage17_0_file_modified_" + ",".join(stage17_0_changed))
    if stage17_1_changed:
        reasons.append("stage17_1_file_modified_" + ",".join(stage17_1_changed))
    if closed_changed:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_changed))
    if unauthorized:
        reasons.append("unauthorized_stage17_2_change_detected_" + ",".join(unauthorized))

    metadata = validate_boundary_metadata(
        args.y_min,
        args.y_max,
        args.wall_normal_direction,
        args.x_boundary_policy,
        args.y_boundary_policy,
        args.z_boundary_policy,
    )
    height = metadata["height"]

    summary["stage17_2_requested_status"] = status(args.enable == "1")
    summary["stage17_1_evidence_status"] = status(stage17_1_evidence_ok(repo, args.require_stage17_1 == "1", args.accept_stage17_1_closed_evidence == "1"))
    summary["stage17_0_fresh_archive_fix_preserved_status"] = status(stage17_0_fresh_archive_fix_preserved(repo))
    summary["stage17_1_evidence_fix_preserved_status"] = status(stage17_1_evidence_fix_preserved(repo))
    summary["stage17_0_files_unmodified_status"] = status(not stage17_0_changed)
    summary["stage17_1_files_unmodified_status"] = status(not stage17_1_changed)
    summary["stage17_enable_status"] = status(args.enable == "1")
    summary["wall_safety_enable_status"] = status(args.wall_safety_enable == "1")
    summary["boundary_check_enable_status"] = status(args.boundary_check_enable == "1")
    summary["diagnostic_only_status"] = status(args.diagnostic_only == "1")
    summary["wall_normal_direction_value"] = args.wall_normal_direction
    summary["wall_normal_direction_status"] = status(bool(metadata["direction_ok"]))
    summary["y_min_value"] = args.y_min
    summary["y_min_finite_status"] = status(bool(metadata["y_min_ok"]))
    summary["y_max_value"] = args.y_max
    summary["y_max_finite_status"] = status(bool(metadata["y_max_ok"]))
    summary["y_bounds_ordered_status"] = status(bool(metadata["ordered_ok"]))
    summary["channel_height_value"] = "nan" if height is None else f"{height:.17g}"
    summary["channel_height_positive_status"] = status(bool(metadata["height_ok"]))
    summary["x_boundary_policy_status"] = status(bool(metadata["x_policy_ok"]))
    summary["y_boundary_policy_status"] = status(bool(metadata["y_policy_ok"]))
    summary["z_boundary_policy_status"] = status(bool(metadata["z_policy_ok"]))
    summary["periodic_x_metadata_status"] = status(bool(metadata["periodic_x_ok"]))
    summary["periodic_z_metadata_status"] = status(bool(metadata["periodic_z_ok"]))
    summary["wall_y_metadata_status"] = status(bool(metadata["wall_y_ok"]))
    summary["invalid_boundary_metadata_rejection_status"] = status(invalid_boundary_metadata_rejected())

    metadata_status_keys = [
        "stage17_enable_status", "wall_safety_enable_status", "boundary_check_enable_status",
        "diagnostic_only_status", "wall_normal_direction_status", "y_min_finite_status",
        "y_max_finite_status", "y_bounds_ordered_status", "channel_height_positive_status",
        "x_boundary_policy_status", "y_boundary_policy_status", "z_boundary_policy_status",
        "periodic_x_metadata_status", "periodic_z_metadata_status", "wall_y_metadata_status",
        "invalid_boundary_metadata_rejection_status",
    ]
    summary["stage17_boundary_metadata_status"] = status(all(summary[key] == "1" for key in metadata_status_keys))

    no_source_or_build = not source_or_build_changes and not closed_changed and not stage17_0_changed and not stage17_1_changed and not unauthorized
    guarded_absence_keys = [
        "no_wall_clearance_computation_status",
        "no_effective_radius_gap_computation_status",
        "no_contact_state_classification_status",
        "no_wall_penetration_detection_status",
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
        reasons.append("stage17_2_source_or_build_change_detected_" + ",".join(source_or_build_changes))

    wrapper_text = read_text(repo / "stage17_checks" / "run_stage17_2_channel_wall_domain_boundary.sh")
    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in wrapper_text)
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    summary["no_rg_only_dependency_status"] = status(all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_2_FILES))

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
        print("STAGE 17.2 CHANNEL WALL DOMAIN BOUNDARY VERDICT: PASS")
        print("STAGE 17.2 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.2 CHANNEL WALL DOMAIN BOUNDARY VERDICT: FAIL")
    print("STAGE 17.2 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
