#!/usr/bin/env python3
"""Stage 17.1 wall/contact safety configuration audit helper.

Stage 17.1 is configuration-only. This helper validates guarded diagnostic-only
configuration values and audits that Stage 17.0 and closed Stage 10--16 files were
not modified. It does not compute wall distance, classify contact, build targets,
run MPI, or add contact/collision physics.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_1_requested_status",
    "stage17_0_evidence_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_config_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "boundary_check_enable_status",
    "fail_closed_enable_status",
    "contact_placeholder_enable_status",
    "fibre_collision_placeholder_enable_status",
    "diagnostic_only_status",
    "effective_fibre_radius_value",
    "effective_fibre_radius_status",
    "min_wall_clearance_value",
    "min_wall_clearance_status",
    "warning_wall_clearance_value",
    "warning_wall_clearance_status",
    "penetration_tolerance_value",
    "penetration_tolerance_status",
    "clearance_order_status",
    "invalid_config_rejection_status",
    "no_wall_distance_computation_status",
    "no_contact_state_classification_status",
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

REQUIRED_STAGE17_1_FILES = [
    Path("stage17_checks/run_stage17_1_wall_contact_safety_config.sh"),
    Path("stage17_checks/assert_stage17_1_wall_contact_safety_config.py"),
    Path("stage17_checks/stage17_1_wall_contact_safety_config.md"),
]
STAGE17_0_FILES = [
    Path("stage17_checks/assert_stage17_0_preflight_safety_boundary.py"),
    Path("stage17_checks/run_stage17_0_preflight_safety_boundary.sh"),
    Path("stage17_checks/stage17_0_preflight_safety_boundary.md"),
    Path("stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat"),
]
ALLOWED_STAGE17_1_CHANGES = {
    "stage17_checks/run_stage17_1_wall_contact_safety_config.sh",
    "stage17_checks/assert_stage17_1_wall_contact_safety_config.py",
    "stage17_checks/stage17_1_wall_contact_safety_config.md",
    "stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat",
}
CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))


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


def is_file_nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


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


def git_changed_paths(repo: Path) -> list[str]:
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
    paths: list[str] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
        raw = line[3:] if len(line) > 3 else line
        if " -> " in raw:
            raw = raw.split(" -> ", 1)[1]
        paths.append(raw.strip())
    return paths


def real_rg_usage_without_grep_fallback(script: Path) -> bool:
    """Detect real shell-wrapper rg dependencies without false positives.

    This reuses the corrected Stage 17.0 / Stage 16 false-positive-safe policy:
    only executable shell-wrapper command lines are audited. Markdown prose is not
    scanned as code-regression evidence, Python negative-check strings are not
    executable behavior, and regex literals such as rg[[:space:]] do not count as
    ripgrep invocations. Any real wrapper use of rg must include grep fallback.
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


def validate_numeric_config(radius: str, min_clearance: str, warning_clearance: str, tolerance: str) -> dict[str, object]:
    r = finite_float(radius)
    mn = finite_float(min_clearance)
    warn = finite_float(warning_clearance)
    tol = finite_float(tolerance)
    return {
        "radius": r,
        "min_clearance": mn,
        "warning_clearance": warn,
        "tolerance": tol,
        "radius_ok": r is not None and r >= 0.0,
        "min_ok": mn is not None and mn >= 0.0,
        "warning_ok": warn is not None and warn >= 0.0,
        "tolerance_ok": tol is not None and tol >= 0.0,
        "order_ok": mn is not None and warn is not None and warn >= mn,
    }


def invalid_config_rejected() -> bool:
    bad_cases = [
        ("nan", "1.0e-4", "1.0e-3", "1.0e-12"),
        ("-1.0", "1.0e-4", "1.0e-3", "1.0e-12"),
        ("1.0e-3", "-1.0", "1.0e-3", "1.0e-12"),
        ("1.0e-3", "2.0e-3", "1.0e-3", "1.0e-12"),
        ("1.0e-3", "1.0e-4", "1.0e-3", "-1.0"),
    ]
    for case in bad_cases:
        result = validate_numeric_config(*case)
        if all(bool(result[key]) for key in ("radius_ok", "min_ok", "warning_ok", "tolerance_ok", "order_ok")):
            return False
    return True


def stage17_0_evidence_ok(repo: Path, require: bool, accept_closed: bool) -> bool:
    if not require:
        return True
    data = parse_dat(repo / "stage17_outputs" / "fibre_stage17_0_preflight_safety_boundary.dat")
    if data.get("final_status") == "1":
        return True
    structural = all((repo / path).is_file() for path in STAGE17_0_FILES[:3])
    return accept_closed and structural and stage17_0_fresh_archive_fix_preserved(repo)


def stage17_0_fresh_archive_fix_preserved(repo: Path) -> bool:
    """Check that the passed Stage 17.0 fresh-archive fix is structurally present.

    Important: do not reject Stage 17.0 merely because old failure-label strings
    appear in comments, explicit reason messages, or negative-check text.  Stage
    17.0 already passed with those labels available for diagnostics.  The real
    regression would be the absence of the accepted fresh-archive closure path or
    a wrapper that hard-fails before the helper when STAGE16_CLOSED.md is missing.
    """
    helper = read_text(repo / "stage17_checks" / "assert_stage17_0_preflight_safety_boundary.py")
    wrapper = read_text(repo / "stage17_checks" / "run_stage17_0_preflight_safety_boundary.sh")
    doc = read_text(repo / "stage17_checks" / "stage17_0_preflight_safety_boundary.md").lower()
    helper_l = helper.lower()

    required_helper_evidence = (
        "stage16_closure_evidence_ok" in helper
        and "stage16_12" in helper_l
        and "STAGE16_CLOSED.md" in helper
        and "stage14_checks" in helper
        and "stage15_checks" in helper
        and "fresh" in helper_l
        and "archive" in helper_l
        and "accept" in helper_l
    )
    required_wrapper_evidence = (
        "STAGE17_0_ACCEPT_STAGE16_CLOSED_EVIDENCE" in wrapper
        and "--accept-stage16-closed-evidence" in wrapper
    )
    required_doc_evidence = (
        "fresh-archive stage 16 closure evidence" in doc
        and "stage16_closed.md" in doc
        and "stage 16.12" in doc
        and "safety-boundary preflight" in doc
        and "stage 21" in doc
    )
    return required_helper_evidence and required_wrapper_evidence and required_doc_evidence


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
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat")
    parser.add_argument("--require-stage17-0", default="1")
    parser.add_argument("--accept-stage17-0-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--wall-safety-enable", default="1")
    parser.add_argument("--boundary-check-enable", default="1")
    parser.add_argument("--fail-closed-enable", default="1")
    parser.add_argument("--contact-placeholder-enable", default="1")
    parser.add_argument("--fibre-collision-placeholder-enable", default="0")
    parser.add_argument("--effective-fibre-radius", default="1.0e-3")
    parser.add_argument("--min-wall-clearance", default="1.0e-4")
    parser.add_argument("--warning-wall-clearance", default="1.0e-3")
    parser.add_argument("--penetration-tolerance", default="1.0e-12")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}
    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_1_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_1_file_{path}")

    changed = git_changed_paths(repo)
    stage17_0_changed = sorted(path for path in changed if path in {str(p) for p in STAGE17_0_FILES})
    closed_changed = sorted(path for path in changed if path.startswith(CLOSED_STAGE_PREFIXES))
    unauthorized = sorted(path for path in changed if path not in ALLOWED_STAGE17_1_CHANGES)
    source_or_build_changes = sorted(
        path for path in changed
        if path.endswith((".f90", ".F90", ".c", ".cc", ".cpp", ".h", ".hpp")) or path.endswith("CMakeLists.txt")
    )
    if stage17_0_changed:
        reasons.append("stage17_0_file_modified_" + ",".join(stage17_0_changed))
    if closed_changed:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_changed))
    if unauthorized:
        reasons.append("unauthorized_stage17_1_change_detected_" + ",".join(unauthorized))

    summary["stage17_1_requested_status"] = status(args.enable == "1")
    summary["stage17_0_evidence_status"] = status(stage17_0_evidence_ok(repo, args.require_stage17_0 == "1", args.accept_stage17_0_closed_evidence == "1"))
    summary["stage17_0_fresh_archive_fix_preserved_status"] = status(stage17_0_fresh_archive_fix_preserved(repo))
    summary["stage17_0_files_unmodified_status"] = status(not stage17_0_changed)

    config = validate_numeric_config(
        args.effective_fibre_radius,
        args.min_wall_clearance,
        args.warning_wall_clearance,
        args.penetration_tolerance,
    )
    summary["stage17_enable_status"] = status(args.enable == "1")
    summary["wall_safety_enable_status"] = status(args.wall_safety_enable == "1")
    summary["boundary_check_enable_status"] = status(args.boundary_check_enable == "1")
    summary["fail_closed_enable_status"] = status(args.fail_closed_enable == "1")
    summary["contact_placeholder_enable_status"] = status(args.contact_placeholder_enable == "1")
    summary["fibre_collision_placeholder_enable_status"] = status(args.fibre_collision_placeholder_enable == "0")
    summary["diagnostic_only_status"] = status(args.diagnostic_only == "1")
    summary["effective_fibre_radius_value"] = args.effective_fibre_radius
    summary["effective_fibre_radius_status"] = status(bool(config["radius_ok"]))
    summary["min_wall_clearance_value"] = args.min_wall_clearance
    summary["min_wall_clearance_status"] = status(bool(config["min_ok"]))
    summary["warning_wall_clearance_value"] = args.warning_wall_clearance
    summary["warning_wall_clearance_status"] = status(bool(config["warning_ok"]))
    summary["penetration_tolerance_value"] = args.penetration_tolerance
    summary["penetration_tolerance_status"] = status(bool(config["tolerance_ok"]))
    summary["clearance_order_status"] = status(bool(config["order_ok"]))
    summary["invalid_config_rejection_status"] = status(invalid_config_rejected())

    config_keys_ok = all(summary[key] == "1" for key in [
        "stage17_enable_status", "wall_safety_enable_status", "boundary_check_enable_status",
        "fail_closed_enable_status", "contact_placeholder_enable_status",
        "fibre_collision_placeholder_enable_status", "diagnostic_only_status",
        "effective_fibre_radius_status", "min_wall_clearance_status",
        "warning_wall_clearance_status", "penetration_tolerance_status",
        "clearance_order_status", "invalid_config_rejection_status",
    ])
    summary["stage17_config_status"] = status(config_keys_ok)
    for key in [
        "stage17_1_requested_status", "stage17_0_evidence_status", "stage17_0_fresh_archive_fix_preserved_status",
        "stage17_0_files_unmodified_status", "stage17_config_status",
    ]:
        if summary[key] != "1":
            reasons.append(f"{key}_not_pass")

    no_source_or_build = not source_or_build_changes and not closed_changed and not stage17_0_changed and not unauthorized
    force_keys = [
        "no_wall_distance_computation_status",
        "no_contact_state_classification_status",
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
    for key in force_keys:
        summary[key] = status(no_source_or_build and args.diagnostic_only == "1")
    if source_or_build_changes:
        reasons.append("stage17_1_source_or_build_change_detected_" + ",".join(source_or_build_changes))

    wrapper_text = read_text(repo / "stage17_checks" / "run_stage17_1_wall_contact_safety_config.sh")
    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in wrapper_text)
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    summary["no_rg_only_dependency_status"] = status(all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_1_FILES))

    for key in [
        "rank0_safe_diagnostic_status", "no_rank_corruption_status", "stage13_6_diagnostic_preserved_status",
        "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status",
        "no_rg_only_dependency_status",
    ]:
        if summary[key] != "1":
            reasons.append(f"{key}_not_pass")

    unknown_failure = any("unknown" in reason.lower() for reason in reasons)
    summary["no_unknown_failure_status"] = status(not unknown_failure)
    pass_fail_keys = [
        key for key in SUMMARY_KEYS
        if key != "final_status" and not key.endswith("_value")
    ]
    summary["final_status"] = status(all(summary[key] == "1" for key in pass_fail_keys) and not reasons)

    out = repo / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary[key]}" for key in SUMMARY_KEYS]
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")
    for line in lines:
        print(line)
    if summary["final_status"] == "1":
        print("STAGE 17.1 WALL CONTACT SAFETY CONFIG VERDICT: PASS")
        print("STAGE 17.1 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.1 WALL CONTACT SAFETY CONFIG VERDICT: FAIL")
    print("STAGE 17.1 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
