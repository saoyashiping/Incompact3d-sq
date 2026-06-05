#!/usr/bin/env python3
"""Stage 16.7 small-lambda bounded-response audit helper."""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

REQUIRED_KEYS = [
    "stage16_7_requested_status",
    "np",
    "zero_lambda_value",
    "small_lambda_value",
    "small_lambda_nonzero_status",
    "one_fibre_count_status",
    "closed_loop_path_status",
    "stage11_sampling_status",
    "stage12_feedback_status",
    "stage16_4_force_input_status",
    "structure_force_input_finite_status",
    "structure_force_input_bounded_status",
    "controlled_structure_update_status",
    "controlled_structure_update_bounded_status",
    "stage13_force_density_status",
    "stage14_rhs_status",
    "zero_rhs_increment_value",
    "zero_rhs_increment_status",
    "small_rhs_increment_value",
    "small_rhs_increment_nonzero_status",
    "small_rhs_increment_bounded_status",
    "rhs_scaling_status",
    "fluid_signature_delta",
    "fluid_signature_bounded_status",
    "approved_stage12_13_14_chain_status",
    "no_direct_rhs_injection_status",
    "no_production_hook_status",
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
    "final_status",
]
NUMERIC_KEYS = {
    "np",
    "zero_lambda_value",
    "small_lambda_value",
    "zero_rhs_increment_value",
    "small_rhs_increment_value",
    "fluid_signature_delta",
}
PASS_ONE_KEYS = [key for key in REQUIRED_KEYS if key not in NUMERIC_KEYS | {"final_status"}]

STAGE15_DOC_SCRIPT_STEMS = {
    1: "structure_state_buffer",
    2: "velocity_source_adapter",
    3: "structure_advance_formula",
    4: "production_structure_hook",
    5: "structure_noop_invariance",
    6: "controlled_structure_step_np1",
    7: "feedback_linkage",
    8: "controlled_structure_parallel_consistency",
    9: "io_restart_structure_state",
    10: "rhs_ibm_structure_contamination_audit",
    11: "total_smoke_closure",
}


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


def finite_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def check_rg_fallback(script: Path) -> bool:
    """Check shell scripts only for real rg command usage with a grep fallback.

    Stage 16.7 intentionally reuses the corrected Stage 16.6 / 16.5 / 16.4
    false-positive-safe audit pattern.  Documentation files are not scanned as executable
    evidence; negative-check strings are not counted as regressions; regex literals such as
    rg[[:space:]] inside helpers are not treated as command invocations.  Only shell wrappers
    with actual rg command/fallback logic are audited here.
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


def add_static_audit_reasons(repo: Path, reasons: list[str]) -> None:
    """Apply Stage 16.7 static checks without broad false-positive-prone scans.

    False-positive protections are required for Stage 16.7 and are copied from the passed
    Stage 16.6 / 16.5 / 16.4 helper style:
    * .md documentation is checked for required-file existence only, not scanned as code behavior.
    * Python helper negative-check strings and regex literals are not treated as source behavior.
    * Legitimate Stage 13.5 conservation/sign audit files are allowed; only old Stage 13.5
      production force-density diagnostic names in real production/check logic are rejected.
    * If source behavior cannot be distinguished from documentation or negative checks, this
      helper fails closed by reporting a specific missing/ambiguous evidence reason rather than
      accepting a permissive match.
    """
    cmake = repo / "src" / "CMakeLists.txt"
    src_files = all_files(repo / "src", ("*.f90",))
    src_text = joined(src_files)
    xcompact3d_text = read_text(repo / "src" / "xcompact3d.f90")
    stage11_14_shell_text = joined(
        all_files(repo / "stage11_checks", ("*.sh",))
        + all_files(repo / "stage13_checks", ("*.sh",))
        + all_files(repo / "stage14_checks", ("*.sh",))
    )
    stage15_shell_files = all_files(repo / "stage15_checks", ("*.sh",))
    stage15_text = joined(stage15_shell_files)

    if "fibre_stage16_small_lambda_response" in xcompact3d_text or "stage16_7" in xcompact3d_text.lower():
        reasons.append("stage16_7_production_hook_inserted_into_xcompact3d")
    if re.search(r"stage14_get_injection_gain\s*\(\s*\)\s*==\s*0\.0", src_text + stage11_14_shell_text):
        reasons.append("stage14_forbidden_lambda_zero_hook_gate_detected")
    if "fibre_stage14_production_rhs_injection" not in src_text:
        reasons.append("stage14_production_rhs_injection_source_missing")
    if "stage14_small_lambda" not in src_text + stage11_14_shell_text and "small_lambda" not in src_text + stage11_14_shell_text:
        reasons.append("stage14_small_lambda_diagnostic_evidence_missing")

    production_context = src_text + joined(all_files(repo / "stage13_checks", ("*.sh", "*.py")))
    if "stage13_6_production_force_density_candidate_status" not in production_context:
        reasons.append("stage13_6_production_force_density_status_missing")
    if "fibre_stage13_6_production_force_density_candidate.dat" not in production_context:
        reasons.append("stage13_6_production_force_density_dat_missing")
    if re.search(r"stage13_5_production_force_density_candidate", production_context):
        reasons.append("old_stage13_5_production_force_density_name_detected")
    if "local_subdomain_center" in production_context or "subdomain_center" in production_context:
        reasons.append("stage13_local_subdomain_center_sampling_regression_detected")

    rank0_context = src_text + stage11_14_shell_text + stage15_text
    for marker in ["stage11", "stage13", "stage14", "stage15", "stage16"]:
        if marker in rank0_context and "rank0_write_allowed" not in rank0_context and "rank0" not in rank0_context.lower():
            reasons.append(f"{marker}_rank0_safe_write_evidence_missing")

    cmake_text = read_text(cmake)
    for marker, label in [
        ("fibre_stage16_small_lambda_response_check", "stage16_7_build_registration_missing"),
        ("fibre_stage16_lambda0_no_contamination_check", "stage16_6_build_registration_missing"),
        ("fibre_stage16_closed_loop_dryrun_check", "stage16_5_build_registration_missing"),
        ("fibre_stage16_structure_force_input_check", "stage16_4_build_registration_missing"),
    ]:
        if marker not in cmake_text:
            reasons.append(label)

    for path in [repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh"]:
        if path.exists() and not check_rg_fallback(path):
            reasons.append(f"rg_without_grep_fallback_{path}")

    wrapper_text = read_text(repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh")
    helper_text = read_text(repo / "stage16_checks" / "assert_stage16_7_small_lambda_bounded_response_np1.py")
    # Unknown-failure auditing follows the same false-positive-safe rule: inspect executable
    # wrapper behavior, not this helper's negative-check reason strings.
    if "unknown failure" in wrapper_text.lower():
        reasons.append("ambiguous_unknown_failure_fallback_text_detected")
    if "Stage 16.6 / 16.5 / 16.4" not in helper_text:
        reasons.append("false_positive_safe_pattern_comment_missing")

    for forbidden in ["wall_contact_enable = 1", "contact_enable = 1", "multi_fibre_enable = 1", "legacy_ibm_forcing_enable = 1"]:
        if forbidden in src_text.lower():
            reasons.append(f"forbidden_activation_detected_{forbidden.replace(' ', '_')}")


def stage16_6_evidence_ok(repo: Path, accept_closed: bool) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_6_lambda0_no_fluid_contamination.dat")
    if data.get("final_status") == "1":
        return True
    required_files = [
        repo / "stage16_checks" / "run_stage16_6_lambda0_no_fluid_contamination.sh",
        repo / "stage16_checks" / "assert_stage16_6_lambda0_no_fluid_contamination.py",
        repo / "stage16_checks" / "stage16_6_lambda0_no_fluid_contamination.md",
        repo / "src" / "fibre_stage16_lambda0_no_contamination.f90",
        repo / "src" / "fibre_stage16_lambda0_no_contamination_check.f90",
    ]
    return accept_closed and all(path.exists() for path in required_files)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-6", default="1")
    parser.add_argument("--accept-stage16-6-closed-evidence", default="1")
    parser.add_argument("--max-zero-rhs-increment", default="1.0e-14")
    parser.add_argument("--max-rhs-increment", default="1.0e-8")
    parser.add_argument("--max-fluid-delta", default="1.0e-8")
    parser.add_argument("--min-rhs-increment", default="1.0e-20")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    dat = out_dir / "fibre_stage16_7_small_lambda_bounded_response_np1.dat"
    reasons_file = out_dir / "stage16_7_small_lambda_bounded_response_np1_reasons.tmp"
    data = parse_dat(dat)
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_7_small_lambda_response_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_7_small_lambda_response_target_run_status_not_pass")

    required_files = [
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "stage16_checks" / "stage16_7_small_lambda_bounded_response_np1.md",
        repo / "stage16_checks" / "assert_stage16_7_small_lambda_bounded_response_np1.py",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    for path in required_files:
        if not path.exists():
            reasons.append(f"missing_required_stage16_7_file_{path}")
    if "add_executable(fibre_stage16_small_lambda_response_check" not in read_text(repo / "src" / "CMakeLists.txt"):
        reasons.append("missing_fibre_stage16_small_lambda_response_check_build_registration")

    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    if args.require_stage16_6 == "1" and not stage16_6_evidence_ok(repo, args.accept_stage16_6_closed_evidence == "1"):
        reasons.append("missing_or_failed_stage16_6_lambda0_no_contamination_evidence")

    add_static_audit_reasons(repo, reasons)

    for key in REQUIRED_KEYS:
        if key not in data and key != "final_status":
            reasons.append(f"missing_required_output_key_{key}")
            data[key] = "0"
    for key in PASS_ONE_KEYS:
        if data.get(key) != "1":
            reasons.append(f"{key}_not_pass")

    max_zero = finite_float(args.max_zero_rhs_increment) or 0.0
    max_rhs = finite_float(args.max_rhs_increment) or 0.0
    max_fluid = finite_float(args.max_fluid_delta) or 0.0
    min_rhs = finite_float(args.min_rhs_increment) or 0.0

    np_value = finite_float(data.get("np", "nan"))
    zero_lambda = finite_float(data.get("zero_lambda_value", "nan"))
    small_lambda = finite_float(data.get("small_lambda_value", "nan"))
    zero_rhs = finite_float(data.get("zero_rhs_increment_value", "nan"))
    small_rhs = finite_float(data.get("small_rhs_increment_value", "nan"))
    fluid_delta = finite_float(data.get("fluid_signature_delta", "nan"))

    if np_value != 1.0:
        reasons.append("np_not_one")
    if zero_lambda != 0.0:
        reasons.append("zero_lambda_value_not_zero")
    if small_lambda is None or small_lambda <= 0.0:
        reasons.append("small_lambda_value_not_positive")
    if zero_rhs is None or abs(zero_rhs) > max_zero:
        reasons.append("zero_rhs_increment_exceeds_bound")
    if small_rhs is None or abs(small_rhs) < min_rhs:
        reasons.append("small_rhs_increment_below_minimum")
    if small_rhs is None or abs(small_rhs) > max_rhs:
        reasons.append("small_rhs_increment_exceeds_bound")
    if fluid_delta is None or abs(fluid_delta) > max_fluid:
        reasons.append("fluid_signature_delta_exceeds_bound")
    for key in NUMERIC_KEYS:
        if finite_float(data.get(key, "nan")) is None:
            reasons.append(f"{key}_not_finite_numeric")

    final_ok = len(reasons) == 0
    data["final_status"] = "1" if final_ok else "0"

    with dat.open("w") as handle:
        handle.write("# Stage 16.7 np=1 small-lambda bounded-response summary\n")
        for key in REQUIRED_KEYS:
            handle.write(f"{key} {data.get(key, '0')}\n")
        for key in ["numeric_parse_status", "numeric_bounds_status", "feedback_alpha", "stage13_force_density_signature"]:
            if key in data and key not in REQUIRED_KEYS:
                handle.write(f"{key} {data[key]}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.7 SMALL-LAMBDA BOUNDED RESPONSE NP1 VERDICT: PASS")
        print("STAGE 16.7 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.7 SMALL-LAMBDA BOUNDED RESPONSE NP1 VERDICT: FAIL")
    print("STAGE 16.7 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
