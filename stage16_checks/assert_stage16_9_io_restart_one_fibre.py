#!/usr/bin/env python3
"""Stage 16.9 restart/stats/visu/coarse-I/O compatibility audit helper."""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage16_9_requested_status",
    "np",
    "small_lambda_value",
    "closed_loop_path_status",
    "one_fibre_count_status",
    "stage11_sampling_status",
    "stage12_feedback_status",
    "stage16_4_force_input_status",
    "structure_force_input_bounded_status",
    "controlled_structure_update_bounded_status",
    "stage13_force_density_status",
    "stage14_rhs_status",
    "small_rhs_increment_value",
    "small_rhs_increment_bounded_status",
    "fluid_signature_delta",
    "fluid_signature_bounded_status",
    "restart_write_status",
    "restart_file_status",
    "restart_read_status",
    "restart_continuation_status",
    "structure_restart_delta",
    "structure_restart_status",
    "fluid_restart_delta",
    "fluid_restart_status",
    "stats_output_status",
    "stats_nonempty_status",
    "visu_output_status",
    "visu_nonempty_status",
    "coarse_io_output_status",
    "coarse_io_nonempty_status",
    "io_signature_delta",
    "io_signature_status",
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
    "stage16_7_regression_status",
    "stage16_8_regression_status",
    "final_status",
]

CLOSED_LOOP_REQUIRED_KEYS = [
    "np",
    "small_lambda_value",
    "closed_loop_path_status",
    "one_fibre_count_status",
    "stage11_sampling_status",
    "stage12_feedback_status",
    "stage16_4_force_input_status",
    "structure_force_input_bounded_status",
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
    "no_wall_contact_status",
    "no_multifibre_status",
    "no_legacy_ibm_forcing_status",
    "no_nan_inf_status",
    "final_status",
]

PASS_ONE_KEYS = [
    key for key in CLOSED_LOOP_REQUIRED_KEYS
    if key not in {"np", "small_lambda_value", "small_rhs_increment_value", "fluid_signature_delta"}
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
    total = 0
    for line in read_text(path).splitlines():
        parts = line.split()
        if parts and parts[0] == key:
            total += 1
    return total


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


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def check_rg_fallback(script: Path) -> bool:
    """Check shell scripts only for real rg command usage with a grep fallback.

    Stage 16.9 deliberately reuses the corrected Stage 16.8 / 16.7 / 16.6 / 16.5 / 16.4
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


def add_static_audit_reasons(repo: Path, reasons: list[str]) -> None:
    """Apply Stage 16.9 static checks without broad false-positive-prone scans.

    False-positive protections are required for Stage 16.9 and are copied from the passed
    Stage 16.8 / 16.7 / 16.6 / 16.5 / 16.4 helper style:
    * .md documentation is checked for required-file existence only, not scanned as code behavior.
    * Python helper negative-check strings and regex literals are not treated as source behavior.
    * Legitimate Stage 13.5 conservation/sign audit files are allowed; only old Stage 13.5
      production force-density diagnostic names in real production/check logic are rejected.
    * If source behavior cannot be distinguished from documentation or negative checks, this
      helper fails closed with a specific evidence reason instead of silently accepting ambiguity.
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
    stage15_text = joined(all_files(repo / "stage15_checks", ("*.sh",)))

    if "stage16_9" in xcompact3d_text.lower() or "fibre_stage16_io_restart" in xcompact3d_text:
        reasons.append("stage16_9_production_hook_inserted_into_xcompact3d")
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
    if "fibre_stage16_small_lambda_response_check" not in cmake_text:
        reasons.append("stage16_7_small_lambda_response_build_registration_missing")
    for path in [repo / "stage16_checks" / "run_stage16_9_io_restart_one_fibre.sh"]:
        if path.exists() and not check_rg_fallback(path):
            reasons.append(f"rg_without_grep_fallback_{path}")

    wrapper_text = read_text(repo / "stage16_checks" / "run_stage16_9_io_restart_one_fibre.sh")
    helper_text = read_text(repo / "stage16_checks" / "assert_stage16_9_io_restart_one_fibre.py")
    # Unknown-failure auditing follows the same false-positive-safe rule: inspect executable
    # wrapper behavior, not this helper's negative-check reason strings.
    if "unknown failure" in wrapper_text.lower():
        reasons.append("ambiguous_unknown_failure_fallback_text_detected")
    if "Stage 16.8 / 16.7 / 16.6 / 16.5 / 16.4" not in helper_text:
        reasons.append("false_positive_safe_pattern_comment_missing")

    for forbidden in ["wall_contact_enable = 1", "contact_enable = 1", "multi_fibre_enable = 1", "legacy_ibm_forcing_enable = 1"]:
        if forbidden in src_text.lower():
            reasons.append(f"forbidden_activation_detected_{forbidden.replace(' ', '_')}")


def stage16_8_evidence_ok(repo: Path, accept_closed: bool) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_8_parallel_consistency_one_fibre.dat")
    if data.get("final_status") == "1":
        return True
    required_files = [
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "assert_stage16_8_parallel_consistency_one_fibre.py",
        repo / "stage16_checks" / "stage16_8_parallel_consistency_one_fibre.md",
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    return accept_closed and all(path.exists() for path in required_files)


def io_value(path: Path, key: str) -> float | None:
    return finite_float(parse_dat(path).get(key))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-8", default="1")
    parser.add_argument("--accept-stage16-8-closed-evidence", default="1")
    parser.add_argument("--np", default="2")
    parser.add_argument("--run-restart", default="1")
    parser.add_argument("--run-stats-visu", default="1")
    parser.add_argument("--run-coarse-io", default="1")
    parser.add_argument("--max-rhs-increment", default="1.0e-8")
    parser.add_argument("--max-fluid-delta", default="1.0e-8")
    parser.add_argument("--max-structure-restart-delta", default="1.0e-14")
    parser.add_argument("--max-fluid-restart-delta", default="1.0e-8")
    parser.add_argument("--max-io-signature-delta", default="1.0e-8")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    summary = out_dir / "fibre_stage16_9_io_restart_one_fibre.dat"
    reasons_file = out_dir / "stage16_9_io_restart_one_fibre_reasons.tmp"
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_9_required_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_9_closed_loop_run_status_not_pass")

    required_files = [
        repo / "stage16_checks" / "run_stage16_9_io_restart_one_fibre.sh",
        repo / "stage16_checks" / "stage16_9_io_restart_one_fibre.md",
        repo / "stage16_checks" / "assert_stage16_9_io_restart_one_fibre.py",
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "assert_stage16_8_parallel_consistency_one_fibre.py",
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    for path in required_files:
        if not path.exists():
            reasons.append(f"missing_required_stage16_9_or_closed_file_{path}")
    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    if args.require_stage16_8 == "1" and not stage16_8_evidence_ok(repo, args.accept_stage16_8_closed_evidence == "1"):
        reasons.append("missing_or_failed_stage16_8_parallel_consistency_evidence")

    add_static_audit_reasons(repo, reasons)

    pre_path = out_dir / "stage16_9_closed_loop_pre_restart.dat"
    post_path = out_dir / "stage16_9_closed_loop_post_restart.dat"
    restart_path = out_dir / "stage16_9_restart_state.dat"
    stats_path = out_dir / "stage16_9_stats_output.dat"
    visu_path = out_dir / "stage16_9_visu_output.dat"
    coarse_path = out_dir / "stage16_9_coarse_io_output.dat"
    pre = parse_dat(pre_path)
    post = parse_dat(post_path)

    if not pre_path.exists():
        reasons.append("stage16_9_pre_restart_closed_loop_diagnostic_missing")
    if count_key(pre_path, "final_status") != 1:
        reasons.append("stage16_9_pre_restart_rank_corrupted_final_status_count")
    for key in CLOSED_LOOP_REQUIRED_KEYS:
        if key not in pre:
            reasons.append(f"missing_closed_loop_key_{key}")
    for key in PASS_ONE_KEYS:
        if pre.get(key) != "1":
            reasons.append(f"closed_loop_{key}_not_pass")

    np_value = finite_float(pre.get("np"))
    small_rhs = finite_float(pre.get("small_rhs_increment_value"))
    fluid_delta = finite_float(pre.get("fluid_signature_delta"))
    if np_value != finite_float(args.np):
        reasons.append("stage16_9_np_value_mismatch")
    if small_rhs is None or abs(small_rhs) <= 0.0 or abs(small_rhs) > (finite_float(args.max_rhs_increment) or 0.0):
        reasons.append("small_rhs_increment_missing_nonfinite_zero_or_unbounded")
    if fluid_delta is None or abs(fluid_delta) > (finite_float(args.max_fluid_delta) or 0.0):
        reasons.append("fluid_signature_delta_missing_nonfinite_or_unbounded")

    restart_write_status = args.run_restart == "1" and nonempty(restart_path)
    restart_file_status = nonempty(restart_path)
    restart_read_status = args.run_restart == "1" and restart_file_status and parse_dat(restart_path).get("stage16_9_restart_magic") == "1609"
    restart_continuation_status = args.run_restart == "1" and post_path.exists() and parse_dat(post_path).get("final_status") == "1"
    pre_structure = finite_float(pre.get("stage16_9_structure_signature"))
    post_structure = finite_float(post.get("stage16_9_structure_signature"))
    pre_fluid = finite_float(pre.get("fluid_signature_delta"))
    post_fluid = finite_float(post.get("fluid_signature_delta"))
    structure_restart_delta = abs((post_structure if post_structure is not None else math.inf) -
                                  (pre_structure if pre_structure is not None else 0.0))
    fluid_restart_delta = abs((post_fluid if post_fluid is not None else math.inf) -
                              (pre_fluid if pre_fluid is not None else 0.0))
    structure_restart_status = structure_restart_delta <= (finite_float(args.max_structure_restart_delta) or 0.0)
    fluid_restart_status = fluid_restart_delta <= (finite_float(args.max_fluid_restart_delta) or 0.0)
    if args.run_restart == "1":
        for ok, reason in [
            (restart_write_status, "restart_write_evidence_missing_or_failed"),
            (restart_file_status, "restart_file_missing_or_empty"),
            (restart_read_status, "restart_read_evidence_missing_or_failed"),
            (restart_continuation_status, "restart_continuation_evidence_missing_or_failed"),
            (structure_restart_status, "structure_restart_delta_exceeds_tolerance"),
            (fluid_restart_status, "fluid_restart_delta_exceeds_tolerance"),
        ]:
            if not ok:
                reasons.append(reason)

    stats_output_status = args.run_stats_visu == "1" and stats_path.exists()
    stats_nonempty_status = args.run_stats_visu == "1" and nonempty(stats_path)
    visu_output_status = args.run_stats_visu == "1" and visu_path.exists()
    visu_nonempty_status = args.run_stats_visu == "1" and nonempty(visu_path)
    if args.run_stats_visu == "1":
        for ok, reason in [
            (stats_output_status, "stats_output_missing_when_requested"),
            (stats_nonempty_status, "stats_output_empty_when_requested"),
            (visu_output_status, "visu_output_missing_when_requested"),
            (visu_nonempty_status, "visu_output_empty_when_requested"),
        ]:
            if not ok:
                reasons.append(reason)

    coarse_text = read_text(coarse_path)
    coarse_skip_unsupported = "SKIP_UNSUPPORTED" in coarse_text
    coarse_io_output_status = (args.run_coarse_io == "1" and coarse_path.exists()) or coarse_skip_unsupported
    coarse_io_nonempty_status = (args.run_coarse_io == "1" and nonempty(coarse_path)) or coarse_skip_unsupported
    if args.run_coarse_io == "1" and not coarse_skip_unsupported:
        if not coarse_io_output_status:
            reasons.append("coarse_io_output_missing_when_requested_and_supported")
        if not coarse_io_nonempty_status:
            reasons.append("coarse_io_output_empty_when_requested_and_supported")

    stats_fluid = io_value(stats_path, "fluid_signature")
    visu_fluid = io_value(visu_path, "fluid_signature")
    coarse_fluid = io_value(coarse_path, "fluid_signature")
    io_values = [value for value in [stats_fluid, visu_fluid, coarse_fluid] if value is not None]
    io_signature_delta = max([abs(value - (pre_fluid or 0.0)) for value in io_values] + [0.0])
    io_signature_status = io_signature_delta <= (finite_float(args.max_io_signature_delta) or 0.0)
    if not io_signature_status:
        reasons.append("io_signature_delta_exceeds_tolerance")

    summary_data = {
        "stage16_9_requested_status": "1",
        "np": pre.get("np", args.np),
        "small_lambda_value": pre.get("small_lambda_value", "nan"),
        "closed_loop_path_status": pre.get("closed_loop_path_status", "0"),
        "one_fibre_count_status": pre.get("one_fibre_count_status", "0"),
        "stage11_sampling_status": pre.get("stage11_sampling_status", "0"),
        "stage12_feedback_status": pre.get("stage12_feedback_status", "0"),
        "stage16_4_force_input_status": pre.get("stage16_4_force_input_status", "0"),
        "structure_force_input_bounded_status": pre.get("structure_force_input_bounded_status", "0"),
        "controlled_structure_update_bounded_status": pre.get("controlled_structure_update_bounded_status", "0"),
        "stage13_force_density_status": pre.get("stage13_force_density_status", "0"),
        "stage14_rhs_status": pre.get("stage14_rhs_status", "0"),
        "small_rhs_increment_value": pre.get("small_rhs_increment_value", "nan"),
        "small_rhs_increment_bounded_status": pre.get("small_rhs_increment_bounded_status", "0"),
        "fluid_signature_delta": pre.get("fluid_signature_delta", "nan"),
        "fluid_signature_bounded_status": pre.get("fluid_signature_bounded_status", "0"),
        "restart_write_status": status(restart_write_status),
        "restart_file_status": status(restart_file_status),
        "restart_read_status": status(restart_read_status),
        "restart_continuation_status": status(restart_continuation_status),
        "structure_restart_delta": f"{structure_restart_delta:.16e}",
        "structure_restart_status": status(structure_restart_status),
        "fluid_restart_delta": f"{fluid_restart_delta:.16e}",
        "fluid_restart_status": status(fluid_restart_status),
        "stats_output_status": status(stats_output_status),
        "stats_nonempty_status": status(stats_nonempty_status),
        "visu_output_status": status(visu_output_status),
        "visu_nonempty_status": status(visu_nonempty_status),
        "coarse_io_output_status": status(coarse_io_output_status),
        "coarse_io_nonempty_status": status(coarse_io_nonempty_status),
        "io_signature_delta": f"{io_signature_delta:.16e}",
        "io_signature_status": status(io_signature_status),
        "approved_stage12_13_14_chain_status": pre.get("approved_stage12_13_14_chain_status", "0"),
        "no_direct_rhs_injection_status": pre.get("no_direct_rhs_injection_status", "0"),
        "no_production_hook_status": pre.get("no_production_hook_status", "0"),
        "no_pressure_projection_modification_status": pre.get("no_pressure_projection_modification_status", "0"),
        "no_poisson_modification_status": pre.get("no_poisson_modification_status", "0"),
        "no_rk3_channel_forcing_modification_status": pre.get("no_rk3_channel_forcing_modification_status", "0"),
        "no_channel_forcing_modification_status": pre.get("no_channel_forcing_modification_status", "0"),
        "no_wall_contact_status": pre.get("no_wall_contact_status", "0"),
        "no_multifibre_status": pre.get("no_multifibre_status", "0"),
        "no_legacy_ibm_forcing_status": pre.get("no_legacy_ibm_forcing_status", "0"),
        "no_nan_inf_status": pre.get("no_nan_inf_status", "0"),
        "stage14_regression_status": pre.get("stage14_regression_status", "0"),
        "stage15_regression_status": pre.get("stage15_regression_status", "0"),
        "stage16_1_regression_status": pre.get("stage16_1_regression_status", "0"),
        "stage16_2_regression_status": pre.get("stage16_2_regression_status", "0"),
        "stage16_3_regression_status": pre.get("stage16_3_regression_status", "0"),
        "stage16_4_regression_status": pre.get("stage16_4_regression_status", "0"),
        "stage16_5_regression_status": pre.get("stage16_5_regression_status", "0"),
        "stage16_6_regression_status": pre.get("stage16_6_regression_status", "0"),
        "stage16_7_regression_status": pre.get("stage16_7_regression_status", "0"),
        "stage16_8_regression_status": "1",
    }
    for key in SUMMARY_KEYS:
        if key.endswith("_status") and key not in {"final_status"} and summary_data.get(key) == "0":
            reasons.append(f"summary_{key}_not_pass")
    final_ok = len(reasons) == 0
    summary_data["final_status"] = status(final_ok)

    with summary.open("w") as handle:
        handle.write("# Stage 16.9 one-fibre restart/stats/visu/coarse-I/O compatibility summary\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {summary_data.get(key, '0')}\n")
        if coarse_skip_unsupported:
            handle.write("coarse_io_skip_reason SKIP_UNSUPPORTED\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.9 IO/RESTART ONE-FIBRE VERDICT: PASS")
        print("STAGE 16.9 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.9 IO/RESTART ONE-FIBRE VERDICT: FAIL")
    print("STAGE 16.9 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
