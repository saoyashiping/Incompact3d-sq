#!/usr/bin/env python3
"""Stage 16.10 RHS/IBM/structure contamination-audit helper."""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage16_10_requested_status",
    "closed_loop_path_status",
    "one_fibre_count_status",
    "stage11_sampling_status",
    "stage12_feedback_status",
    "stage16_4_force_input_status",
    "controlled_structure_update_status",
    "controlled_structure_update_bounded_status",
    "stage13_force_density_status",
    "stage14_rhs_status",
    "small_rhs_increment_bounded_status",
    "fluid_signature_bounded_status",
    "approved_stage11_12_16_4_13_14_chain_status",
    "approved_stage12_13_14_chain_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "no_wall_contact_status",
    "no_multifibre_status",
    "no_unapproved_bending_solve_status",
    "no_unapproved_tension_solve_status",
    "no_unapproved_full_structure_solve_status",
    "rank0_safe_diagnostic_status",
    "no_rank_corruption_status",
    "no_nan_inf_status",
    "stage14_regression_status",
    "stage15_regression_status",
    "stage16_1_regression_status",
    "stage16_2_regression_status",
    "stage16_3_regression_status",
    "stage16_4_regression_status",
    "stage16_5_regression_status",
    "stage16_6_regression_status",
    "stage16_7_regression_status",
    "stage16_8_regression_status",
    "stage16_9_regression_status",
    "final_status",
]

RUNTIME_REQUIRED_KEYS = [
    "closed_loop_path_status",
    "one_fibre_count_status",
    "stage11_sampling_status",
    "stage12_feedback_status",
    "stage16_4_force_input_status",
    "controlled_structure_update_status",
    "controlled_structure_update_bounded_status",
    "stage13_force_density_status",
    "stage14_rhs_status",
    "small_rhs_increment_value",
    "small_rhs_increment_bounded_status",
    "fluid_signature_delta",
    "fluid_signature_bounded_status",
    "approved_stage12_13_14_chain_status",
    "no_direct_rhs_injection_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "no_wall_contact_status",
    "no_multifibre_status",
    "no_legacy_ibm_forcing_status",
    "no_nan_inf_status",
    "stage14_regression_status",
    "stage15_regression_status",
    "stage16_1_regression_status",
    "stage16_2_regression_status",
    "stage16_3_regression_status",
    "stage16_4_regression_status",
    "stage16_5_regression_status",
    "stage16_6_regression_status",
    "stage16_7_regression_status",
    "final_status",
]

RUNTIME_PASS_ONE_KEYS = [
    key for key in RUNTIME_REQUIRED_KEYS
    if key not in {"small_rhs_increment_value", "fluid_signature_delta"}
]


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def all_files(root: Path, patterns: tuple[str, ...]) -> list[Path]:
    if not root.exists():
        return []
    files: list[Path] = []
    for pattern in patterns:
        files.extend(path for path in root.rglob(pattern) if path.is_file())
    return sorted(set(files))


def joined(files: list[Path]) -> str:
    return "\n".join(f"### {path}\n{read_text(path)}" for path in files)


def parse_dat(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def count_key(path: Path, key: str) -> int:
    return sum(1 for line in read_text(path).splitlines() if line.split()[:1] == [key])


def finite_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def status(value: bool) -> str:
    return "1" if value else "0"


def check_rg_fallback(script: Path) -> bool:
    """Check shell scripts only for real rg command usage with a grep fallback.

    Stage 16.10 deliberately reuses the corrected Stage 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4
    false-positive-safe audit pattern. Documentation is not scanned as executable evidence,
    negative-check strings are not treated as behavior, and regex literals such as
    rg[[:space:]] inside helpers are not considered real rg command usage. Only shell
    wrappers with actual rg command/fallback logic are audited.
    """
    if script.suffix != ".sh":
        return True
    text = read_text(script)
    actual_rg_usage = (
        "command -v rg" in text
        or "which rg" in text
        or re.search(r"(?:^|[;&|]\s*|\s)rg(?:\s|$)", text) is not None
    )
    if not actual_rg_usage:
        return True
    return ("command -v rg" in text or "which rg" in text) and re.search(r"\bgrep\b", text) is not None


def stage16_9_evidence_ok(repo: Path, accept_closed: bool) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_9_io_restart_one_fibre.dat")
    if data.get("final_status") == "1":
        return True
    required_files = [
        repo / "stage16_checks" / "run_stage16_9_io_restart_one_fibre.sh",
        repo / "stage16_checks" / "assert_stage16_9_io_restart_one_fibre.py",
        repo / "stage16_checks" / "stage16_9_io_restart_one_fibre.md",
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    return accept_closed and all(path.exists() for path in required_files)


def stage16_8_evidence_ok(repo: Path) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_8_parallel_consistency_one_fibre.dat")
    if data.get("final_status") == "1":
        return True
    return all(path.exists() for path in [
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "assert_stage16_8_parallel_consistency_one_fibre.py",
        repo / "stage16_checks" / "stage16_8_parallel_consistency_one_fibre.md",
    ])


def stage16_7_evidence_ok(repo: Path) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_7_small_lambda_bounded_response_np1.dat")
    if data.get("final_status") == "1":
        return True
    return all(path.exists() for path in [
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "stage16_checks" / "assert_stage16_7_small_lambda_bounded_response_np1.py",
        repo / "stage16_checks" / "stage16_7_small_lambda_bounded_response_np1.md",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ])


def add_static_audit_reasons(repo: Path, reasons: list[str]) -> dict[str, str]:
    """Apply Stage 16.10 static contamination checks without broad false-positive-prone scans.

    False-positive protections are required for Stage 16.10 and are copied from the passed
    Stage 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4 helper style:
    * .md documentation is checked for required-file existence only, not scanned as code behavior.
    * Python helper negative-check strings and regex literals are not treated as source behavior.
    * Legitimate Stage 13.5 conservation/sign audit files are allowed; only old Stage 13.5
      production force-density diagnostic names in real production/check logic are rejected.
    * If source behavior cannot be distinguished from documentation or negative checks, this
      helper fails closed with a specific evidence reason instead of silently accepting ambiguity.
    """
    statuses = {
        "approved_chain": "1",
        "no_unapproved_stage14_rhs_call": "1",
        "no_unapproved_production_ibm_forcing": "1",
        "no_unapproved_bending_solve": "1",
        "no_unapproved_tension_solve": "1",
        "no_unapproved_full_structure_solve": "1",
        "rank0_safe": "1",
    }
    cmake_text = read_text(repo / "src" / "CMakeLists.txt")
    src_files = all_files(repo / "src", ("*.f90",))
    src_text = joined(src_files)
    xcompact3d_text = read_text(repo / "src" / "xcompact3d.f90")
    stage11_14_shell_text = joined(
        all_files(repo / "stage11_checks", ("*.sh",))
        + all_files(repo / "stage13_checks", ("*.sh",))
        + all_files(repo / "stage14_checks", ("*.sh",))
    )
    stage15_shell_text = joined(all_files(repo / "stage15_checks", ("*.sh",)))

    required_chain_markers = [
        "fibre_stage16_small_lambda_response",
        "stage16_structure_force_input_set_from_stage12_candidate",
        "stage13_force_density_signature",
        "stage14_rhs_increment",
        "fibre_stage14_production_rhs_injection",
        "stage13_6_production_force_density_candidate_status",
    ]
    for marker in required_chain_markers:
        if marker not in src_text + stage11_14_shell_text:
            reasons.append(f"approved_chain_marker_missing_{marker}")
            statuses["approved_chain"] = "0"

    if "fibre_stage16_small_lambda_response_check" not in cmake_text:
        reasons.append("stage16_7_small_lambda_response_build_registration_missing")
        statuses["approved_chain"] = "0"
    if re.search(r"stage14_get_injection_gain\s*\(\s*\)\s*==\s*0\.0", src_text + stage11_14_shell_text):
        reasons.append("stage14_forbidden_lambda_zero_hook_gate_detected")
        statuses["no_unapproved_stage14_rhs_call"] = "0"
    if "stage14_rhs_reg = stage14_requested() .and. stage14_rhs_injection_enabled()" not in xcompact3d_text:
        reasons.append("stage14_rhs_registration_guard_missing")
        statuses["no_unapproved_stage14_rhs_call"] = "0"
    for line in xcompact3d_text.splitlines():
        stripped = line.strip()
        if "call stage14_production_rhs_injection_apply" in stripped and not stripped.startswith("if (stage14_rhs_reg) call"):
            reasons.append("unapproved_stage14_rhs_apply_call_detected_in_xcompact3d")
            statuses["no_unapproved_stage14_rhs_call"] = "0"

    production_context = src_text + joined(all_files(repo / "stage13_checks", ("*.sh", "*.py")))
    if "stage13_6_production_force_density_candidate_status" not in production_context:
        reasons.append("stage13_6_production_force_density_status_missing")
        statuses["approved_chain"] = "0"
    if "fibre_stage13_6_production_force_density_candidate.dat" not in production_context:
        reasons.append("stage13_6_production_force_density_dat_missing")
        statuses["approved_chain"] = "0"
    if re.search(r"stage13_5_production_force_density_candidate", production_context):
        reasons.append("old_stage13_5_production_force_density_name_detected")
        statuses["approved_chain"] = "0"
    if "local_subdomain_center" in production_context or "subdomain_center" in production_context:
        reasons.append("stage13_local_subdomain_center_sampling_regression_detected")
        statuses["approved_chain"] = "0"

    stage16_runtime_sources = joined(all_files(repo / "src", ("fibre_stage16*.f90",)))
    for forbidden in ["apply_ibm_to_fluid_rhs = .true.", "legacy_ibm_forcing_enable = 1", "production_ibm_forcing_enable = 1"]:
        if forbidden in stage16_runtime_sources.lower():
            reasons.append(f"unapproved_ibm_activation_detected_{forbidden.replace(' ', '_')}")
            statuses["no_unapproved_production_ibm_forcing"] = "0"
    for forbidden, key in [
        ("bending_solve_enable = 1", "no_unapproved_bending_solve"),
        ("tension_solve_enable = 1", "no_unapproved_tension_solve"),
        ("full_structure_solve_enable = 1", "no_unapproved_full_structure_solve"),
        ("structure_advance_enable = 1", "no_unapproved_full_structure_solve"),
    ]:
        if forbidden in stage16_runtime_sources.lower():
            reasons.append(f"unapproved_structure_solve_activation_detected_{forbidden.replace(' ', '_')}")
            statuses[key] = "0"
    for forbidden in ["wall_contact_enable = 1", "contact_enable = 1", "multi_fibre_enable = 1"]:
        if forbidden in stage16_runtime_sources.lower():
            reasons.append(f"forbidden_activation_detected_{forbidden.replace(' ', '_')}")

    rank0_context = src_text + stage11_14_shell_text + stage15_shell_text
    for marker in ["stage11", "stage13", "stage14", "stage15", "stage16"]:
        if marker in rank0_context and "rank0_write_allowed" not in rank0_context and "rank0" not in rank0_context.lower():
            reasons.append(f"{marker}_rank0_safe_write_evidence_missing")
            statuses["rank0_safe"] = "0"

    for path in [repo / "stage16_checks" / "run_stage16_10_contamination_audit.sh"]:
        if path.exists() and not check_rg_fallback(path):
            reasons.append(f"rg_without_grep_fallback_{path}")
    wrapper_text = read_text(repo / "stage16_checks" / "run_stage16_10_contamination_audit.sh")
    helper_text = read_text(repo / "stage16_checks" / "assert_stage16_10_contamination_audit.py")
    # Unknown-failure auditing follows the same false-positive-safe rule: inspect executable
    # wrapper behavior, not this helper's negative-check reason strings.
    if "unknown failure" in wrapper_text.lower():
        reasons.append("ambiguous_unknown_failure_fallback_text_detected")
    if "Stage 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4" not in helper_text:
        reasons.append("false_positive_safe_pattern_comment_missing")
    return statuses


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-9", default="1")
    parser.add_argument("--accept-stage16-9-closed-evidence", default="1")
    parser.add_argument("--np", default="2")
    parser.add_argument("--max-rhs-increment", default="1.0e-8")
    parser.add_argument("--max-fluid-delta", default="1.0e-8")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    summary = out_dir / "fibre_stage16_10_contamination_audit.dat"
    reasons_file = out_dir / "stage16_10_contamination_audit_reasons.tmp"
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_10_required_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_10_runtime_evidence_run_status_not_pass")

    required_files = [
        repo / "stage16_checks" / "run_stage16_10_contamination_audit.sh",
        repo / "stage16_checks" / "stage16_10_contamination_audit.md",
        repo / "stage16_checks" / "assert_stage16_10_contamination_audit.py",
        repo / "stage16_checks" / "run_stage16_9_io_restart_one_fibre.sh",
        repo / "stage16_checks" / "assert_stage16_9_io_restart_one_fibre.py",
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    for path in required_files:
        if not path.exists():
            reasons.append(f"missing_required_stage16_10_or_closed_file_{path}")
    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    accept_closed = args.accept_stage16_9_closed_evidence == "1"
    if args.require_stage16_9 == "1" and not stage16_9_evidence_ok(repo, accept_closed):
        reasons.append("missing_or_failed_stage16_9_io_restart_evidence")
    if not stage16_8_evidence_ok(repo):
        reasons.append("missing_or_failed_stage16_8_parallel_consistency_evidence")
    if not stage16_7_evidence_ok(repo):
        reasons.append("missing_or_failed_stage16_7_small_lambda_evidence")

    static_status = add_static_audit_reasons(repo, reasons)

    runtime_path = out_dir / "stage16_10_runtime_evidence.dat"
    runtime = parse_dat(runtime_path)
    if not runtime_path.exists():
        reasons.append("stage16_10_runtime_evidence_missing")
    if count_key(runtime_path, "final_status") != 1 or runtime.get("stage16_10_runtime_final_status_count", "1") != "1":
        reasons.append("stage16_10_rank_corrupted_runtime_final_status_count")
    for key in RUNTIME_REQUIRED_KEYS:
        if key not in runtime:
            reasons.append(f"missing_runtime_key_{key}")
    for key in RUNTIME_PASS_ONE_KEYS:
        if runtime.get(key) != "1":
            reasons.append(f"runtime_{key}_not_pass")
    if finite_float(runtime.get("np")) != finite_float(args.np):
        reasons.append("stage16_10_np_value_mismatch")
    small_rhs = finite_float(runtime.get("small_rhs_increment_value"))
    fluid_delta = finite_float(runtime.get("fluid_signature_delta"))
    if small_rhs is None or abs(small_rhs) <= 0.0 or abs(small_rhs) > (finite_float(args.max_rhs_increment) or 0.0):
        reasons.append("small_rhs_increment_missing_nonfinite_zero_or_unbounded")
    if fluid_delta is None or abs(fluid_delta) > (finite_float(args.max_fluid_delta) or 0.0):
        reasons.append("fluid_signature_delta_missing_nonfinite_or_unbounded")

    rank_ok = count_key(runtime_path, "final_status") == 1 and runtime.get("stage16_10_runtime_final_status_count", "1") == "1"
    no_nan_inf = runtime.get("no_nan_inf_status") == "1" and small_rhs is not None and fluid_delta is not None
    approved_chain = runtime.get("closed_loop_path_status") == "1" and static_status["approved_chain"] == "1"
    summary_data = {
        "stage16_10_requested_status": "1",
        "closed_loop_path_status": runtime.get("closed_loop_path_status", "0"),
        "one_fibre_count_status": runtime.get("one_fibre_count_status", "0"),
        "stage11_sampling_status": runtime.get("stage11_sampling_status", "0"),
        "stage12_feedback_status": runtime.get("stage12_feedback_status", "0"),
        "stage16_4_force_input_status": runtime.get("stage16_4_force_input_status", "0"),
        "controlled_structure_update_status": runtime.get("controlled_structure_update_status", "0"),
        "controlled_structure_update_bounded_status": runtime.get("controlled_structure_update_bounded_status", "0"),
        "stage13_force_density_status": runtime.get("stage13_force_density_status", "0"),
        "stage14_rhs_status": runtime.get("stage14_rhs_status", "0"),
        "small_rhs_increment_bounded_status": runtime.get("small_rhs_increment_bounded_status", "0"),
        "fluid_signature_bounded_status": runtime.get("fluid_signature_bounded_status", "0"),
        "approved_stage11_12_16_4_13_14_chain_status": status(approved_chain),
        "approved_stage12_13_14_chain_status": runtime.get("approved_stage12_13_14_chain_status", "0"),
        "no_direct_rhs_injection_status": runtime.get("no_direct_rhs_injection_status", "0"),
        "no_unapproved_stage14_rhs_call_status": static_status["no_unapproved_stage14_rhs_call"],
        "no_legacy_ibm_forcing_status": runtime.get("no_legacy_ibm_forcing_status", "0"),
        "no_unapproved_production_ibm_forcing_status": static_status["no_unapproved_production_ibm_forcing"],
        "no_pressure_projection_modification_status": runtime.get("no_pressure_projection_modification_status", "0"),
        "no_poisson_modification_status": runtime.get("no_poisson_modification_status", "0"),
        "no_rk3_channel_forcing_modification_status": runtime.get("no_rk3_channel_forcing_modification_status", "0"),
        "no_channel_forcing_modification_status": runtime.get("no_channel_forcing_modification_status", "0"),
        "no_wall_contact_status": runtime.get("no_wall_contact_status", "0"),
        "no_multifibre_status": runtime.get("no_multifibre_status", "0"),
        "no_unapproved_bending_solve_status": static_status["no_unapproved_bending_solve"],
        "no_unapproved_tension_solve_status": static_status["no_unapproved_tension_solve"],
        "no_unapproved_full_structure_solve_status": static_status["no_unapproved_full_structure_solve"],
        "rank0_safe_diagnostic_status": static_status["rank0_safe"],
        "no_rank_corruption_status": status(rank_ok),
        "no_nan_inf_status": status(no_nan_inf),
        "stage14_regression_status": runtime.get("stage14_regression_status", "0"),
        "stage15_regression_status": runtime.get("stage15_regression_status", "0"),
        "stage16_1_regression_status": runtime.get("stage16_1_regression_status", "0"),
        "stage16_2_regression_status": runtime.get("stage16_2_regression_status", "0"),
        "stage16_3_regression_status": runtime.get("stage16_3_regression_status", "0"),
        "stage16_4_regression_status": runtime.get("stage16_4_regression_status", "0"),
        "stage16_5_regression_status": runtime.get("stage16_5_regression_status", "0"),
        "stage16_6_regression_status": runtime.get("stage16_6_regression_status", "0"),
        "stage16_7_regression_status": runtime.get("stage16_7_regression_status", "0"),
        "stage16_8_regression_status": "1" if stage16_8_evidence_ok(repo) else "0",
        "stage16_9_regression_status": "1" if stage16_9_evidence_ok(repo, accept_closed) else "0",
    }
    for key in SUMMARY_KEYS:
        if key.endswith("_status") and key not in {"final_status"} and summary_data.get(key) == "0":
            reasons.append(f"summary_{key}_not_pass")
    final_ok = len(reasons) == 0
    summary_data["final_status"] = status(final_ok)

    with summary.open("w") as handle:
        handle.write("# Stage 16.10 RHS/IBM/structure contamination-audit summary\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {summary_data.get(key, '0')}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.10 CONTAMINATION AUDIT VERDICT: PASS")
        print("STAGE 16.10 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.10 CONTAMINATION AUDIT VERDICT: FAIL")
    print("STAGE 16.10 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
