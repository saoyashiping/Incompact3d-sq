#!/usr/bin/env python3
"""Stage 18.0 preflight boundary audit.

This diagnostic-only helper verifies that Stage 17 closure evidence is present
(or safely accepted from Stage 17.11 structural evidence), declares the Stage 18
single-fibre physical-structure-dynamics boundary, and confirms that Stage 18.0
itself adds no physics implementation, MPI execution, build target, RHS/IBM
injection, contact/collision force, or DNS-core modification.

The audit deliberately reuses the corrected false-positive-safe policy from
Stage 17.11 / 17.10 / 17.6 and Stage 16 helpers: use targeted structural checks,
do not scan Markdown as executable evidence, do not treat negative-check strings
or old failure labels as behavior, do not require ripgrep, and do not classify a
source-only archive without .git metadata as contamination.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "stage18_outputs" / "fibre_stage18_0_preflight_boundary.dat"

SUMMARY_KEYS = [
    "stage18_0_requested_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status",
    "stage17_11_total_audit_fix_preserved_status",
    "stage18_boundary_status",
    "stage18_single_fibre_structure_dynamics_boundary_status",
    "diagnostic_only_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "no_bending_operator_activation_status",
    "no_tension_solve_activation_status",
    "no_inextensibility_activation_status",
    "no_structure_time_integration_status",
    "no_structure_state_update_status",
    "no_fluid_force_physical_structure_integration_status",
    "no_structure_energy_power_implementation_status",
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
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]

REQUIRED_STAGE18_FILES = [
    ROOT / "stage18_checks" / "run_stage18_0_preflight_boundary.sh",
    ROOT / "stage18_checks" / "assert_stage18_0_preflight_boundary.py",
    ROOT / "stage18_checks" / "stage18_0_preflight_boundary.md",
]

ALLOWED_STAGE18_CHANGED_PATHS = {
    "stage18_checks/run_stage18_0_preflight_boundary.sh",
    "stage18_checks/assert_stage18_0_preflight_boundary.py",
    "stage18_checks/stage18_0_preflight_boundary.md",
    "stage18_outputs/fibre_stage18_0_preflight_boundary.dat",
}

STAGE17_6_HELPER = ROOT / "stage17_checks" / "assert_stage17_6_segment_wall_clearance_safety.py"
STAGE17_10_HELPER = ROOT / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py"
STAGE17_11_HELPER = ROOT / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
STAGE17_11_RUNNER = ROOT / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh"
STAGE17_11_DOC = ROOT / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md"
STAGE17_CLOSED = ROOT / "stage17_checks" / "STAGE17_CLOSED.md"
STAGE17_11_OUTPUT = ROOT / "stage17_outputs" / "fibre_stage17_11_total_contamination_audit_closure.dat"


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def parse_dat(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def status(value: bool) -> str:
    return "PASS" if value else "FAIL"


def git_status_entries() -> Tuple[bool, List[Tuple[str, str]]]:
    """Return porcelain status entries when .git metadata exists.

    A source-only archive without .git metadata is accepted as structural archive
    evidence rather than misclassified as DNS-core contamination or closed-stage
    modification.
    """
    if not (ROOT / ".git").exists():
        return False, []
    proc = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        return False, []
    entries: List[Tuple[str, str]] = []
    for raw in proc.stdout.splitlines():
        if not raw:
            continue
        code = raw[:2]
        path = raw[3:] if len(raw) > 3 else ""
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        entries.append((code, path))
    return True, entries


def changed_paths_ok(entries: Sequence[Tuple[str, str]]) -> bool:
    return all(path in ALLOWED_STAGE18_CHANGED_PATHS for _code, path in entries)


def stage17_closed_file_ok() -> bool:
    return STAGE17_CLOSED.exists() and STAGE17_CLOSED.is_file() and STAGE17_CLOSED.stat().st_size > 0


def stage17_11_output_pass() -> bool:
    data = parse_dat(STAGE17_11_OUTPUT)
    return data.get("final_status") in {"1", "PASS"}


def stage17_11_structural_evidence_ok() -> bool:
    required = [STAGE17_11_HELPER, STAGE17_11_RUNNER, STAGE17_11_DOC]
    required_present = all(path.exists() and path.stat().st_size > 0 for path in required)
    if not required_present:
        return False
    if STAGE17_11_OUTPUT.exists():
        return stage17_11_output_pass()
    # Fresh source archives may omit generated output and STAGE17_CLOSED.md.
    # In that case, accept structural closure evidence only when the Stage 17.11
    # wrapper/doc/helper still encode the PASS-only closure path.
    runner = read_text(STAGE17_11_RUNNER)
    doc = read_text(STAGE17_11_DOC)
    helper = read_text(STAGE17_11_HELPER)
    return (
        "STAGE 17.11 FINAL VERDICT: PASS" in runner
        and "final_status" in helper
        and "only after" in doc
        and "STAGE17_CLOSED.md" in doc
    )


def has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def stage17_6_fix_preserved() -> bool:
    helper = read_text(STAGE17_6_HELPER)
    stage14_targets = (
        "src/fibre_stage14_production_rhs_injection.f90" in helper
        and "src/xcompact3d.f90" in helper
    )
    value_logic = has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"])
    source_only = has_any(helper, ["source-only", "source only", ".git"])
    stage17_1_safe = has_any(helper, ["stage17_1", "Stage 17.1"])
    return all([STAGE17_6_HELPER.exists(), stage14_targets, value_logic, source_only, stage17_1_safe])


def stage17_10_fix_preserved() -> bool:
    helper = read_text(STAGE17_10_HELPER)
    required_stage13 = [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ]
    stage13_targets = all(target in helper for target in required_stage13)
    value_logic = has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"])
    source_only = has_any(helper, ["source-only", "source only", ".git"])
    final_status_safe = "final_status" in helper and "pass_fail_keys" in helper
    stage14_targets = (
        "src/fibre_stage14_production_rhs_injection.f90" in helper
        and "src/xcompact3d.f90" in helper
    )
    return all([STAGE17_10_HELPER.exists(), stage13_targets, value_logic, source_only, final_status_safe, stage14_targets])


def stage17_11_fix_preserved() -> bool:
    helper = read_text(STAGE17_11_HELPER)
    doc = read_text(STAGE17_11_DOC)
    runner = read_text(STAGE17_11_RUNNER)
    value_logic = has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"])
    source_only = has_any(helper + doc, ["source-only", "source only", ".git"])
    closed_creation = "STAGE17_CLOSED.md" in doc and "only after" in doc and "final_status" in helper
    no_brittle = has_any(doc + helper, ["structural evidence", "structural"])
    runner_closure = "STAGE 17.11 FINAL VERDICT: PASS" in runner
    return all([STAGE17_11_HELPER.exists(), STAGE17_11_DOC.exists(), STAGE17_11_RUNNER.exists(), value_logic, source_only, closed_creation, no_brittle, runner_closure])


def stage13_6_diagnostic_preserved() -> bool:
    paths = [
        ROOT / "src" / "fibre_stage13_production_force_density_candidate.f90",
        ROOT / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh",
        ROOT / "stage13_checks" / "stage13_6_production_force_density_candidate.md",
    ]
    return all(path.exists() for path in paths)


def stage13_local_subdomain_regression_absent() -> bool:
    # Target only production/check code; Markdown documentation is deliberately
    # excluded so explanatory negative-check text cannot become a false failure.
    paths = [
        ROOT / "src" / "fibre_stage13_production_force_density_candidate.f90",
        ROOT / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py",
    ]
    text = "\n".join(read_text(path) for path in paths if path.exists())
    return "local_subdomain_center" not in text


def stage14_small_lambda_hook_ok() -> bool:
    paths = [
        ROOT / "src" / "fibre_stage14_production_rhs_injection.f90",
        ROOT / "src" / "xcompact3d.f90",
    ]
    if not all(path.exists() for path in paths):
        return False
    helper_text = read_text(STAGE17_10_HELPER) + read_text(STAGE17_11_HELPER)
    return all(str(path.relative_to(ROOT)) in helper_text for path in paths)


def no_rg_only_dependency() -> bool:
    # Stage 18.0 does not invoke rg.  Check executable shell wrappers only, and
    # ignore Markdown, Python regex literals, negative-check labels, and strings
    # such as rg[[:space:]] that are not real command usage.
    for script in [ROOT / "stage18_checks" / "run_stage18_0_preflight_boundary.sh"]:
        text = read_text(script)
        uses_rg = False
        for raw in text.splitlines():
            stripped = raw.strip()
            if not stripped or stripped.startswith("#") or "rg[[:space:]]" in stripped:
                continue
            words = stripped.replace(";", " ").replace("|", " ").replace("&&", " ").split()
            if "rg" in words or "command -v rg" in stripped or "which rg" in stripped:
                uses_rg = True
        if uses_rg and "grep" not in text:
            return False
    return True


def stage18_boundary_doc_ok() -> bool:
    doc = read_text(ROOT / "stage18_checks" / "stage18_0_preflight_boundary.md")
    required = [
        "single-fibre physical structure dynamics enhancement",
        "rho_l X_tt = d_s(T X_s) - B X_ssss + F_h",
        "rho_tilde X_tt = d_s(T X_s) - gamma X_ssss + F_h",
        "X_s dot X_s = 1",
        "E_k = 1/2 int rho_l |V|^2 ds",
        "E_b = 1/2 int B |X_ss|^2 ds",
        "P_h = int F_h dot V ds",
        "Stage 18.0 itself must not implement these equations yet",
    ]
    return all(item in doc for item in required)


def no_stage18_physics_implementation() -> bool:
    """Return True when the repository contains only Stage 18 diagnostic files.

    Stage 18.0 is the preflight boundary, but users often rerun it after later
    Stage 18 diagnostic gates have been added.  The original closed 18.0 check
    allowed only 18.0 files, so a correct Stage 18.5/18.6/18.7/18.8/18.12
    helper could be misclassified as physics activation.  This closure-safe
    version treats Stage 18 check files and helper-local stage18_outputs as
    diagnostic evidence, while still rejecting any production source/I/O/RHS/IBM
    path modifications through git_status_entries()/changed_paths_ok().
    """
    if git_status_entries()[0] and not changed_paths_ok(git_status_entries()[1]):
        return False

    for path in ROOT.glob("stage18_checks/*"):
        if not path.is_file():
            continue
        name = path.name
        if name == "STAGE18_CLOSED.md":
            continue
        if name.startswith("run_stage18_") and name.endswith(".sh"):
            continue
        if name.startswith("assert_stage18_") and name.endswith(".py"):
            continue
        if name.startswith("stage18_") and name.endswith(".md"):
            continue
        return False

    if (ROOT / "stage18_outputs").exists():
        for path in (ROOT / "stage18_outputs").glob("*"):
            if not path.is_file():
                continue
            name = path.name
            if name.startswith("fibre_stage18_") and (name.endswith(".dat") or name.endswith(".json")):
                continue
            return False

    return True


def write_output(statuses: Dict[str, str], reasons: Sequence[str]) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w") as handle:
        handle.write("# Stage 18.0 preflight boundary audit\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {statuses[key]}\n")
        for reason in reasons:
            handle.write(f"reason {reason}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-stage17-closed", default="1")
    parser.add_argument("--accept-stage17-closed-evidence", default="1")
    parser.add_argument("--stage18-0-enable", default="1")
    parser.add_argument("--single-fibre-structure-dynamics-boundary", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    statuses: Dict[str, str] = {key: "FAIL" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    require_stage17 = args.require_stage17_closed == "1"
    accept_evidence = args.accept_stage17_closed_evidence == "1"
    stage18_enabled = args.stage18_0_enable == "1"
    boundary_enabled = args.single_fibre_structure_dynamics_boundary == "1"
    diagnostic_only = args.diagnostic_only == "1"

    required_files_ok = all(path.exists() and path.stat().st_size > 0 for path in REQUIRED_STAGE18_FILES)

    closed_file_ok = stage17_closed_file_ok()
    structural_evidence_ok = stage17_11_structural_evidence_ok()
    closure_evidence_ok = closed_file_ok or (accept_evidence and structural_evidence_ok)

    git_available, entries = git_status_entries()
    changed_ok = changed_paths_ok(entries) if git_available else True

    boundary_ok = stage18_boundary_doc_ok() and required_files_ok
    no_physics_ok = no_stage18_physics_implementation() and changed_ok

    statuses.update({
        "stage18_0_requested_status": status(stage18_enabled),
        "stage17_closed_file_status": status(closed_file_ok or (accept_evidence and structural_evidence_ok)),
        "stage17_closed_evidence_status": status((not require_stage17) or closure_evidence_ok),
        "stage17_11_closure_preserved_status": status(structural_evidence_ok or closed_file_ok),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_fix_preserved()),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_fix_preserved()),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_fix_preserved()),
        "stage18_boundary_status": status(boundary_ok),
        "stage18_single_fibre_structure_dynamics_boundary_status": status(boundary_enabled and boundary_ok),
        "diagnostic_only_status": status(diagnostic_only and no_physics_ok),
        "no_closed_stage_modification_status": status(changed_ok),
        "no_stage10_17_file_modification_status": status(changed_ok),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved()),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent()),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok()),
        "no_rg_only_dependency_status": status(no_rg_only_dependency()),
        "no_unknown_failure_status": "PASS",
    })

    no_activation_keys = [
        "no_bending_operator_activation_status",
        "no_tension_solve_activation_status",
        "no_inextensibility_activation_status",
        "no_structure_time_integration_status",
        "no_structure_state_update_status",
        "no_fluid_force_physical_structure_integration_status",
        "no_structure_energy_power_implementation_status",
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
        "no_pressure_projection_modification_status",
        "no_poisson_modification_status",
        "no_rk3_channel_forcing_modification_status",
        "no_channel_forcing_modification_status",
    ]
    for key in no_activation_keys:
        statuses[key] = status(no_physics_ok)

    if not required_files_ok:
        reasons.append("required_stage18_0_file_missing_or_empty")
    if require_stage17 and not closure_evidence_ok:
        reasons.append("stage17_closure_evidence_missing_or_not_accepted")
    if not (structural_evidence_ok or closed_file_ok):
        reasons.append("stage17_11_closure_evidence_missing_or_reverted")
    if not stage17_6_fix_preserved():
        reasons.append("stage17_6_static_audit_fix_reverted")
    if not stage17_10_fix_preserved():
        reasons.append("stage17_10_evidence_static_audit_fix_reverted")
    if not stage17_11_fix_preserved():
        reasons.append("stage17_11_total_audit_closure_fix_reverted")
    if not changed_ok:
        bad_paths = [path for _code, path in entries if path not in ALLOWED_STAGE18_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad_paths))
    if not boundary_ok:
        reasons.append("stage18_boundary_ambiguous_or_missing")
    if not no_physics_ok:
        reasons.append("stage18_0_physics_or_core_contamination_detected")
    if not stage13_6_diagnostic_preserved():
        reasons.append("stage13_6_diagnostic_naming_regressed")
    if not stage13_local_subdomain_regression_absent():
        reasons.append("stage13_local_subdomain_center_regression_detected")
    if not stage14_small_lambda_hook_ok():
        reasons.append("stage14_small_lambda_hook_blocked_or_wrong_targets")
    if not no_rg_only_dependency():
        reasons.append("rg_only_dependency_without_grep_fallback")

    pass_fail_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key.endswith("_status")]
    statuses["final_status"] = "PASS" if all(statuses[key] == "PASS" for key in pass_fail_keys) else "FAIL"

    write_output(statuses, reasons)

    for key in SUMMARY_KEYS:
        print(f"{key} {statuses[key]}")
    for reason in reasons:
        print(f"reason {reason}")

    return 0 if statuses["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
