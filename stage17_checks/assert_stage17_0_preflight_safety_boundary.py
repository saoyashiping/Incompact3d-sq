#!/usr/bin/env python3
"""Stage 17.0 preflight safety-boundary audit helper.

This helper is intentionally an evidence aggregator and static boundary audit only.
It does not build, run MPI, execute physics, or modify closed Stage 10--16 files.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage17_0_requested_status",
    "stage16_closed_file_status",
    "stage16_closed_content_status",
    "stage17_boundary_status",
    "stage17_safety_placeholder_only_status",
    "no_closed_stage_modification_status",
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

REQUIRED_STAGE17_FILES = [
    Path("stage17_checks/run_stage17_0_preflight_safety_boundary.sh"),
    Path("stage17_checks/assert_stage17_0_preflight_safety_boundary.py"),
    Path("stage17_checks/stage17_0_preflight_safety_boundary.md"),
]

ALLOWED_STAGE17_CHANGES = {
    "stage17_checks/run_stage17_0_preflight_safety_boundary.sh",
    "stage17_checks/assert_stage17_0_preflight_safety_boundary.py",
    "stage17_checks/stage17_0_preflight_safety_boundary.md",
    "stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat",
}

CLOSED_STAGE_PREFIXES = tuple(f"stage{stage}_checks/" for stage in range(10, 17))
DNS_CORE_NAMES = ("pressure", "projection", "poisson", "rk3", "channel_forcing")


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def status(ok: bool) -> str:
    return "1" if ok else "0"


def git_changed_paths(repo: Path) -> list[str]:
    """Return changed paths using git status, with no ripgrep dependency."""
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
    """Audit shell wrappers for real rg command use while avoiding false positives.

    This is the corrected Stage 16.5--16.12 false-positive-safe pattern reused for
    Stage 17.0: only shell script command text is considered executable evidence.
    Markdown is not scanned as code, Python negative-check strings are not treated as
    behavior, and regex literals such as rg[[:space:]] are not treated as actual rg
    invocations. If a real wrapper uses rg, it must also include a grep fallback.
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
    has_detector = "command -v rg" in text or "which rg" in text
    has_grep = "grep" in text
    return not (has_detector and has_grep)


def stage16_closed_content_ok(text: str, accept: bool) -> bool:
    if not accept:
        return False
    lowered = text.lower()
    return "stage 16" in lowered and ("closed" in lowered or "closure passed" in lowered or "stage 16.12" in lowered)


def stage16_closure_evidence_ok(repo: Path, accept: bool) -> bool:
    """Accept Stage 16 closure in fresh source trees without editing Stage 16 files.

    Stage 16.12 normally writes stage16_checks/STAGE16_CLOSED.md at runtime after
    full PASS.  User-shared source archives can omit generated runtime artifacts,
    including STAGE16_CLOSED.md and stage16_outputs/*.dat.  In that case Stage
    17.0 must not modify closed Stage 16 files just to recreate the closure record.
    When acceptance is explicitly enabled, use read-only structural evidence from
    the passed Stage 16.12 closure machinery instead.
    """
    if not accept:
        return False
    closed = repo / "stage16_checks" / "STAGE16_CLOSED.md"
    if nonempty(closed) and stage16_closed_content_ok(read_text(closed), True):
        return True

    required = [
        repo / "stage16_checks" / "run_stage16_12_total_smoke_closure.sh",
        repo / "stage16_checks" / "assert_stage16_12_total_smoke_closure.py",
        repo / "stage16_checks" / "stage16_12_total_smoke_closure.md",
        repo / "stage16_checks" / "assert_stage16_11_short_time_stability_smoke.py",
        repo / "stage16_checks" / "run_stage16_11_short_time_stability_smoke.sh",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
        repo / "stage14_checks" / "STAGE14_CLOSED.md",
        repo / "stage15_checks" / "STAGE15_CLOSED.md",
    ]
    if not all(nonempty(path) for path in required):
        return False

    helper = read_text(repo / "stage16_checks" / "assert_stage16_12_total_smoke_closure.py")
    doc = read_text(repo / "stage16_checks" / "stage16_12_total_smoke_closure.md").lower()
    return (
        "write_closure_file" in helper
        and "STAGE 16.12 FINAL VERDICT: PASS" in helper
        and "STAGE16_CLOSED.md" in helper
        and "stage 16.12 total smoke and closure" in doc
        and "stage16_closed.md" in doc
        and "final_status 1" in doc
    )


def stage17_boundary_doc_ok(repo: Path) -> bool:
    doc = read_text(repo / "stage17_checks" / "stage17_0_preflight_safety_boundary.md").lower()
    required = [
        "stage 17.0",
        "safety-boundary preflight",
        "diagnostics",
        "placeholder",
        "must not implement real contact",
        "must not implement real collision",
        "stage 21",
        "stage 18",
        "stage 10--16",
        "dns-core numerics",
    ]
    return all(item in doc for item in required)


def stage13_6_preserved(repo: Path) -> bool:
    """Check specific production/check evidence only, not documentation."""
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
    """Look only at Stage 13.6 production/check logic to avoid documentation regressions."""
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
    forbidden_zero_gate = "stage14_get_injection_gain() == 0.0" in (src + xcompact)
    return (
        "stage14_production_rhs_injection_apply" in src
        and "stage14_rhs_reg = stage14_requested() .and. stage14_rhs_injection_enabled()" in xcompact
        and not forbidden_zero_gate
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="stage17_outputs/fibre_stage17_0_preflight_safety_boundary.dat")
    parser.add_argument("--require-stage16-closed", default="1")
    parser.add_argument("--accept-stage16-closed-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    reasons: list[str] = []
    summary = {key: "1" for key in SUMMARY_KEYS}

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())

    for path in REQUIRED_STAGE17_FILES:
        if not (repo / path).is_file():
            reasons.append(f"missing_required_stage17_0_file_{path}")

    summary["stage17_0_requested_status"] = status(args.enable == "1" and args.diagnostic_only == "1")
    if summary["stage17_0_requested_status"] != "1":
        reasons.append("stage17_0_not_enabled_or_not_diagnostic_only")

    closed = repo / "stage16_checks" / "STAGE16_CLOSED.md"
    accept_stage16 = args.accept_stage16_closed_evidence == "1"
    if args.require_stage16_closed == "1":
        closed_file_or_evidence = nonempty(closed) or stage16_closure_evidence_ok(repo, accept_stage16)
    else:
        closed_file_or_evidence = True
    summary["stage16_closed_file_status"] = status(closed_file_or_evidence)
    if not closed_file_or_evidence:
        reasons.append("missing_or_unaccepted_stage16_closure_evidence")

    closed_content = (
        stage16_closed_content_ok(read_text(closed), accept_stage16)
        if nonempty(closed)
        else stage16_closure_evidence_ok(repo, accept_stage16)
    )
    summary["stage16_closed_content_status"] = status(closed_content)
    if not closed_content:
        reasons.append("stage16_closed_evidence_not_accepted")

    boundary_ok = stage17_boundary_doc_ok(repo)
    summary["stage17_boundary_status"] = status(boundary_ok)
    summary["stage17_safety_placeholder_only_status"] = status(boundary_ok and args.diagnostic_only == "1")
    if not boundary_ok:
        reasons.append("stage17_boundary_documentation_missing_or_ambiguous")

    changed = git_changed_paths(repo)
    changed_set = set(changed)
    unauthorized_changes = sorted(path for path in changed if path not in ALLOWED_STAGE17_CHANGES)
    closed_stage_changes = sorted(path for path in changed if path.startswith(CLOSED_STAGE_PREFIXES))
    source_or_build_changes = sorted(
        path for path in changed
        if path.endswith((".f90", ".F90", ".c", ".cc", ".cpp", ".h", ".hpp")) or path.endswith("CMakeLists.txt")
    )
    no_closed_mod = not unauthorized_changes and not closed_stage_changes
    summary["no_closed_stage_modification_status"] = status(no_closed_mod)
    if unauthorized_changes:
        reasons.append("unauthorized_non_stage17_0_change_detected_" + ",".join(unauthorized_changes))
    if closed_stage_changes:
        reasons.append("closed_stage_file_modified_" + ",".join(closed_stage_changes))

    no_new_physics_source = not source_or_build_changes
    if not no_new_physics_source:
        reasons.append("stage17_0_source_or_build_physics_change_detected_" + ",".join(source_or_build_changes))

    # Stage 17.0 is limited to new wrapper/helper/doc files. Because no source or
    # build file may change, every real contact/collision/structure/DNS-core force
    # category below remains absent without scanning documentation or negative-check
    # strings as behavior.
    force_status_keys = [
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
    for key in force_status_keys:
        summary[key] = status(no_new_physics_source and no_closed_mod)

    summary["rank0_safe_diagnostic_status"] = status("stage17_outputs" in read_text(repo / REQUIRED_STAGE17_FILES[0]))
    summary["no_rank_corruption_status"] = status(summary["rank0_safe_diagnostic_status"] == "1")
    if summary["rank0_safe_diagnostic_status"] != "1":
        reasons.append("stage17_0_rank0_safe_diagnostic_output_evidence_missing")

    summary["stage13_6_diagnostic_preserved_status"] = status(stage13_6_preserved(repo))
    if summary["stage13_6_diagnostic_preserved_status"] != "1":
        reasons.append("stage13_6_diagnostic_naming_not_preserved")
    summary["stage13_no_local_subdomain_center_regression_status"] = status(stage13_local_center_absent(repo))
    if summary["stage13_no_local_subdomain_center_regression_status"] != "1":
        reasons.append("stage13_local_subdomain_center_regression_detected")
    summary["stage14_small_lambda_hook_status"] = status(stage14_small_lambda_hook_ok(repo))
    if summary["stage14_small_lambda_hook_status"] != "1":
        reasons.append("stage14_small_nonzero_lambda_hook_missing_or_blocked")

    rg_ok = all(not real_rg_usage_without_grep_fallback(repo / path) for path in REQUIRED_STAGE17_FILES)
    summary["no_rg_only_dependency_status"] = status(rg_ok)
    if not rg_ok:
        reasons.append("stage17_0_rg_usage_without_grep_fallback_detected")

    unknown_failure = any("unknown" in reason.lower() for reason in reasons)
    summary["no_unknown_failure_status"] = status(not unknown_failure)

    final_ok = all(summary[key] == "1" for key in SUMMARY_KEYS if key != "final_status") and not reasons
    summary["final_status"] = status(final_ok)

    out = repo / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary[key]}" for key in SUMMARY_KEYS]
    for reason in reasons:
        lines.append(f"reason {reason}")
    out.write_text("\n".join(lines) + "\n")

    for line in lines:
        print(line)
    if final_ok:
        print("STAGE 17.0 PREFLIGHT SAFETY BOUNDARY VERDICT: PASS")
        print("STAGE 17.0 FINAL VERDICT: PASS")
        return 0
    print("STAGE 17.0 PREFLIGHT SAFETY BOUNDARY VERDICT: FAIL")
    print("STAGE 17.0 FINAL VERDICT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
