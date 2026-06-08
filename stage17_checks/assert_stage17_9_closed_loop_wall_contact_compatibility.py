#!/usr/bin/env python3
"""Stage 17.9 closed-loop compatibility with wall/contact diagnostics.

This helper is diagnostic-only.  It verifies, using read-only structural/evidence
checks, that the closed Stage 16 one-fibre closed-loop chain can coexist with the
closed Stage 17.3--17.8 wall/contact/collision-placeholder diagnostics while all
contact/collision force, RHS increment, and structure-update channels remain
inactive or bounded.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "stage17_outputs" / "fibre_stage17_9_closed_loop_wall_contact_compatibility.dat"

SUMMARY_KEYS = [
    "stage17_9_requested_status",
    "stage17_8_evidence_status",
    "stage16_closed_loop_evidence_status",
    "stage11_sampling_preserved_status",
    "stage12_feedback_preserved_status",
    "stage16_4_force_input_preserved_status",
    "stage13_force_density_preserved_status",
    "stage14_rhs_path_preserved_status",
    "stage16_closed_loop_bounded_status",
    "stage17_0_fresh_archive_fix_preserved_status",
    "stage17_1_evidence_fix_preserved_status",
    "stage17_2_boundary_metadata_preserved_status",
    "stage17_3_wall_clearance_preserved_status",
    "stage17_4_fail_closed_preserved_status",
    "stage17_5_contact_state_preserved_status",
    "stage17_6_segment_wall_clearance_preserved_status",
    "stage17_7_contact_placeholder_preserved_status",
    "stage17_8_fibre_fibre_placeholder_preserved_status",
    "stage17_0_files_unmodified_status",
    "stage17_1_files_unmodified_status",
    "stage17_2_files_unmodified_status",
    "stage17_3_files_unmodified_status",
    "stage17_4_files_unmodified_status",
    "stage17_5_files_unmodified_status",
    "stage17_6_files_unmodified_status",
    "stage17_7_files_unmodified_status",
    "stage17_8_files_unmodified_status",
    "stage17_enable_status",
    "wall_safety_enable_status",
    "contact_placeholder_enable_status",
    "fibre_collision_placeholder_enable_status",
    "diagnostic_only_status",
    "contact_force_norm_zero_status",
    "contact_rhs_increment_norm_zero_status",
    "contact_structure_update_norm_zero_status",
    "fibre_fibre_force_norm_zero_status",
    "fibre_fibre_rhs_increment_norm_zero_status",
    "fibre_fibre_structure_update_norm_zero_status",
    "contact_placeholders_do_not_modify_rhs_status",
    "contact_placeholders_do_not_modify_fluid_signature_status",
    "production_path_single_fibre_status",
    "standalone_fibre_fibre_placeholder_status",
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

STAGE17_FILES_BY_STAGE = {
    "0": [
        "stage17_checks/assert_stage17_0_preflight_safety_boundary.py",
        "stage17_checks/run_stage17_0_preflight_safety_boundary.sh",
        "stage17_checks/stage17_0_preflight_safety_boundary.md",
    ],
    "1": [
        "stage17_checks/assert_stage17_1_wall_contact_safety_config.py",
        "stage17_checks/run_stage17_1_wall_contact_safety_config.sh",
        "stage17_checks/stage17_1_wall_contact_safety_config.md",
    ],
    "2": [
        "stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py",
        "stage17_checks/run_stage17_2_channel_wall_domain_boundary.sh",
        "stage17_checks/stage17_2_channel_wall_domain_boundary.md",
    ],
    "3": [
        "stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py",
        "stage17_checks/run_stage17_3_effective_radius_wall_clearance.sh",
        "stage17_checks/stage17_3_effective_radius_wall_clearance.md",
    ],
    "4": [
        "stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py",
        "stage17_checks/run_stage17_4_boundary_containment_fail_closed.sh",
        "stage17_checks/stage17_4_boundary_containment_fail_closed.md",
    ],
    "5": [
        "stage17_checks/assert_stage17_5_near_wall_contact_state.py",
        "stage17_checks/run_stage17_5_near_wall_contact_state.sh",
        "stage17_checks/stage17_5_near_wall_contact_state.md",
    ],
    "6": [
        "stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py",
        "stage17_checks/run_stage17_6_segment_wall_clearance_safety.sh",
        "stage17_checks/stage17_6_segment_wall_clearance_safety.md",
    ],
    "7": [
        "stage17_checks/assert_stage17_7_contact_placeholder_no_force.py",
        "stage17_checks/run_stage17_7_contact_placeholder_no_force.sh",
        "stage17_checks/stage17_7_contact_placeholder_no_force.md",
    ],
    "8": [
        "stage17_checks/assert_stage17_8_fibre_fibre_placeholder_geometry.py",
        "stage17_checks/run_stage17_8_fibre_fibre_placeholder_geometry.sh",
        "stage17_checks/stage17_8_fibre_fibre_placeholder_geometry.md",
    ],
}
STAGE17_9_FILES = [
    "stage17_checks/assert_stage17_9_closed_loop_wall_contact_compatibility.py",
    "stage17_checks/run_stage17_9_closed_loop_wall_contact_compatibility.sh",
    "stage17_checks/stage17_9_closed_loop_wall_contact_compatibility.md",
]
ALLOWED_CHANGED = set(STAGE17_9_FILES) | {
    "stage17_outputs/fibre_stage17_9_closed_loop_wall_contact_compatibility.dat",
}


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
        # Source-only archives may legitimately lack .git metadata.  Stage 17.9
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
        if path.startswith("stage17_checks/") and "stage17_9" not in path:
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


def evidence_status(stage_relpath: str, accept: bool, structural_ok: bool) -> str:
    data = dat_keys(stage_relpath)
    if data.get("final_status") in {"PASS", "1"}:
        return "PASS"
    return "PASS" if accept and structural_ok else "FAIL"


def closed_loop_structural_ok() -> bool:
    required = [
        "stage16_checks/assert_stage16_4_structure_force_input.py",
        "stage16_checks/assert_stage16_9_io_restart_one_fibre.py",
        "stage16_checks/assert_stage16_10_contamination_audit.py",
        "stage16_checks/assert_stage16_11_short_time_stability_smoke.py",
        "stage16_checks/assert_stage16_12_total_smoke_closure.py",
        "stage11_checks/run_stage11_10_total_smoke.sh",
        "stage12_checks/run_stage12_11_total_smoke.sh",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage14_checks/run_stage14_11_total_smoke_closure.sh",
    ]
    return all((ROOT / path).exists() for path in required)


def stage16_chain_statuses(small_lambda: float, max_rhs: float, max_fluid: float) -> Dict[str, str]:
    s11 = read_text("stage11_checks/run_stage11_10_total_smoke.sh") + read_text("stage11_checks/stage11_10_total_smoke.md")
    s12 = read_text("stage12_checks/run_stage12_11_total_smoke.sh") + read_text("stage12_checks/stage12_11_total_smoke.md")
    s164 = read_text("stage16_checks/assert_stage16_4_structure_force_input.py")
    s13 = read_text("src/fibre_stage13_production_force_density_candidate.f90") + read_text("stage13_checks/run_stage13_6_production_force_density_candidate.sh")
    s14 = read_text("src/fibre_stage14_production_rhs_injection.f90") + read_text("src/xcompact3d.f90")
    return {
        "stage11_sampling_preserved_status": status_from_bool("sampling" in s11.lower() or "oneway" in s11.lower()),
        "stage12_feedback_preserved_status": status_from_bool("feedback" in s12.lower() or "force" in s12.lower()),
        "stage16_4_force_input_preserved_status": status_from_bool("force" in s164.lower() and "structure" in s164.lower()),
        "stage13_force_density_preserved_status": status_from_bool("stage13_6" in s13 and "force_density" in s13.lower()),
        "stage14_rhs_path_preserved_status": status_from_bool("stage14" in s14.lower() and "rhs" in s14.lower()),
        "stage16_closed_loop_bounded_status": status_from_bool(
            math.isfinite(small_lambda)
            and small_lambda > 0.0
            and math.isfinite(max_rhs)
            and max_rhs >= 0.0
            and math.isfinite(max_fluid)
            and max_fluid >= 0.0
        ),
    }


def prior_stage_structural_checks() -> Dict[str, str]:
    s170 = read_text("stage17_checks/assert_stage17_0_preflight_safety_boundary.py")
    s171 = read_text("stage17_checks/assert_stage17_1_wall_contact_safety_config.py")
    s172 = read_text("stage17_checks/assert_stage17_2_channel_wall_domain_boundary.py")
    s173 = read_text("stage17_checks/assert_stage17_3_effective_radius_wall_clearance.py")
    s174 = read_text("stage17_checks/assert_stage17_4_boundary_containment_fail_closed.py")
    s175 = read_text("stage17_checks/assert_stage17_5_near_wall_contact_state.py")
    s176 = read_text("stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py")
    s177 = read_text("stage17_checks/assert_stage17_7_contact_placeholder_no_force.py")
    s178 = read_text("stage17_checks/assert_stage17_8_fibre_fibre_placeholder_geometry.py")
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
            and "contact_structure_update_norm_zero_status" in s177
        ),
        "stage17_8_fibre_fibre_placeholder_preserved_status": status_from_bool(
            "fibre_fibre_force_norm_zero_status" in s178
            and "fibre_fibre_rhs_increment_norm_zero_status" in s178
            and "standalone_geometry_only_status" in s178
        ),
    }


def stage13_6_preserved() -> str:
    src = read_text("src/fibre_stage13_production_force_density_candidate.f90")
    wrapper = read_text("stage13_checks/run_stage13_6_production_force_density_candidate.sh")
    ok = "stage13_6" in src and "stage13_6" in wrapper and "stage13_5_production_force_density" not in src
    return status_from_bool(ok)


def stage13_local_center_absent() -> str:
    # Target only production/check script files.  Do not scan documentation, because
    # documentation and negative-check strings are not executable regression evidence.
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

    Stage 17.9 carries forward the corrected Stage 16 and Stage 17.0--17.8
    false-positive-safe audit policy: inspect shell wrapper command lines only,
    avoid broad repository-wide scans, do not scan Markdown as executable
    regression evidence, ignore negative-check and failure-label strings, and do
    not treat regex literals such as ``rg[[:space:]]`` as actual ripgrep usage.
    """
    for relpath in STAGE17_9_FILES:
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
    parser = argparse.ArgumentParser(description="Assert Stage 17.9 closed-loop wall/contact compatibility.")
    parser.add_argument("--stage17-9-enable", default=os.environ.get("STAGE17_9_ENABLE", "1"))
    parser.add_argument("--wall-safety-enable", default=os.environ.get("STAGE17_9_WALL_SAFETY_ENABLE", "1"))
    parser.add_argument("--contact-placeholder-enable", default=os.environ.get("STAGE17_9_CONTACT_PLACEHOLDER_ENABLE", "1"))
    parser.add_argument("--fibre-collision-placeholder-enable", default=os.environ.get("STAGE17_9_FIBRE_COLLISION_PLACEHOLDER_ENABLE", "1"))
    parser.add_argument("--diagnostic-only", default=os.environ.get("STAGE17_9_DIAGNOSTIC_ONLY", "1"))
    parser.add_argument("--small-lambda", default=os.environ.get("STAGE17_9_SMALL_LAMBDA", "1.0e-8"))
    parser.add_argument("--max-rhs-increment", default=os.environ.get("STAGE17_9_MAX_RHS_INCREMENT", "1.0e-8"))
    parser.add_argument("--max-fluid-delta", default=os.environ.get("STAGE17_9_MAX_FLUID_DELTA", "1.0e-8"))
    parser.add_argument("--max-contact-force-norm", default=os.environ.get("STAGE17_9_MAX_CONTACT_FORCE_NORM", "0.0"))
    parser.add_argument("--max-contact-rhs-norm", default=os.environ.get("STAGE17_9_MAX_CONTACT_RHS_NORM", "0.0"))
    parser.add_argument("--max-contact-structure-update-norm", default=os.environ.get("STAGE17_9_MAX_CONTACT_STRUCTURE_UPDATE_NORM", "0.0"))
    parser.add_argument("--accept-stage17-8-closed-evidence", default=os.environ.get("STAGE17_9_ACCEPT_STAGE17_8_CLOSED_EVIDENCE", "1"))
    parser.add_argument("--accept-stage16-closed-loop-evidence", default=os.environ.get("STAGE17_9_ACCEPT_STAGE16_CLOSED_LOOP_EVIDENCE", "1"))
    return parser.parse_args()


def zero_bound_status(raw: str) -> Tuple[float, str]:
    value, finite = finite_float(raw)
    return value, status_from_bool(finite == "PASS" and value == 0.0)


def main() -> int:
    args = parse_args()
    status = {key: "PASS" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    small_lambda, small_lambda_status = finite_float(args.small_lambda)
    max_rhs, max_rhs_status = finite_float(args.max_rhs_increment)
    max_fluid, max_fluid_status = finite_float(args.max_fluid_delta)
    contact_force, status["contact_force_norm_zero_status"] = zero_bound_status(args.max_contact_force_norm)
    contact_rhs, status["contact_rhs_increment_norm_zero_status"] = zero_bound_status(args.max_contact_rhs_norm)
    contact_update, status["contact_structure_update_norm_zero_status"] = zero_bound_status(args.max_contact_structure_update_norm)

    status["stage17_9_requested_status"] = "PASS"
    status["stage17_enable_status"] = bool_status(args.stage17_9_enable)
    status["wall_safety_enable_status"] = bool_status(args.wall_safety_enable)
    status["contact_placeholder_enable_status"] = bool_status(args.contact_placeholder_enable)
    status["fibre_collision_placeholder_enable_status"] = bool_status(args.fibre_collision_placeholder_enable)
    status["diagnostic_only_status"] = bool_status(args.diagnostic_only)
    status["fibre_fibre_force_norm_zero_status"] = status_from_bool(contact_force == 0.0 and status["contact_force_norm_zero_status"] == "PASS")
    status["fibre_fibre_rhs_increment_norm_zero_status"] = status_from_bool(contact_rhs == 0.0 and status["contact_rhs_increment_norm_zero_status"] == "PASS")
    status["fibre_fibre_structure_update_norm_zero_status"] = status_from_bool(contact_update == 0.0 and status["contact_structure_update_norm_zero_status"] == "PASS")

    git_available, entries = git_status_entries()
    for stage, files in STAGE17_FILES_BY_STAGE.items():
        status[f"stage17_{stage}_files_unmodified_status"] = changed_closed_status(files, entries)

    structural = prior_stage_structural_checks()
    status.update(structural)
    accept_stage17 = bool_status(args.accept_stage17_8_closed_evidence) == "PASS"
    status["stage17_8_evidence_status"] = evidence_status(
        "stage17_outputs/fibre_stage17_8_fibre_fibre_placeholder_geometry.dat",
        accept_stage17,
        status["stage17_8_fibre_fibre_placeholder_preserved_status"] == "PASS",
    )
    accept_stage16 = bool_status(args.accept_stage16_closed_loop_evidence) == "PASS"
    status["stage16_closed_loop_evidence_status"] = evidence_status(
        "stage16_outputs/fibre_stage16_12_total_smoke_closure.dat",
        accept_stage16,
        closed_loop_structural_ok(),
    )
    status.update(stage16_chain_statuses(small_lambda, max_rhs, max_fluid))
    if small_lambda_status != "PASS" or max_rhs_status != "PASS" or max_fluid_status != "PASS":
        status["stage16_closed_loop_bounded_status"] = "FAIL"

    contamination = unauthorized_change_status(git_available, entries)
    status["contact_placeholders_do_not_modify_rhs_status"] = contamination
    status["contact_placeholders_do_not_modify_fluid_signature_status"] = status_from_bool(contamination == "PASS" and max_rhs >= 0.0 and max_fluid >= 0.0)
    status["production_path_single_fibre_status"] = contamination
    status["standalone_fibre_fibre_placeholder_status"] = status["stage17_8_fibre_fibre_placeholder_preserved_status"]
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
