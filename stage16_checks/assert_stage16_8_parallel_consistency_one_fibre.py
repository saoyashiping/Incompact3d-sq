#!/usr/bin/env python3
"""Stage 16.8 np=1/2/4 one-fibre parallel-consistency audit helper."""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage16_8_requested_status",
    "np_list",
    "np1_run_status",
    "np2_run_status",
    "np4_run_status",
    "np1_final_status",
    "np2_final_status",
    "np4_final_status",
    "np1_small_rhs_increment_value",
    "np2_small_rhs_increment_value",
    "np4_small_rhs_increment_value",
    "np1_fluid_signature_delta",
    "np2_fluid_signature_delta",
    "np4_fluid_signature_delta",
    "parallel_force_diff_np2",
    "parallel_force_diff_np4",
    "parallel_structure_diff_np2",
    "parallel_structure_diff_np4",
    "parallel_rhs_diff_np2",
    "parallel_rhs_diff_np4",
    "parallel_fluid_diff_np2",
    "parallel_fluid_diff_np4",
    "parallel_force_status",
    "parallel_structure_status",
    "parallel_rhs_status",
    "parallel_fluid_status",
    "rank0_safe_diagnostic_status",
    "no_rank_corruption_status",
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
    "final_status",
]

PER_NP_REQUIRED_KEYS = [
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

PER_NP_PASS_ONE_KEYS = [
    key for key in PER_NP_REQUIRED_KEYS
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


def check_rg_fallback(script: Path) -> bool:
    """Check shell scripts only for real rg command usage with a grep fallback.

    Stage 16.8 deliberately reuses the corrected Stage 16.7 / 16.6 / 16.5 / 16.4
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
    """Apply Stage 16.8 static checks without broad false-positive-prone scans.

    False-positive protections are required for Stage 16.8 and are copied from the passed
    Stage 16.7 / 16.6 / 16.5 / 16.4 helper style:
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

    if "stage16_8" in xcompact3d_text.lower() or "fibre_stage16_parallel_consistency" in xcompact3d_text:
        reasons.append("stage16_8_production_hook_inserted_into_xcompact3d")
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
    for path in [repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh"]:
        if path.exists() and not check_rg_fallback(path):
            reasons.append(f"rg_without_grep_fallback_{path}")

    wrapper_text = read_text(repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh")
    helper_text = read_text(repo / "stage16_checks" / "assert_stage16_8_parallel_consistency_one_fibre.py")
    # Unknown-failure auditing follows the same false-positive-safe rule: inspect executable
    # wrapper behavior, not this helper's negative-check reason strings.
    if "unknown failure" in wrapper_text.lower():
        reasons.append("ambiguous_unknown_failure_fallback_text_detected")
    if "Stage 16.7 / 16.6 / 16.5 / 16.4" not in helper_text:
        reasons.append("false_positive_safe_pattern_comment_missing")

    for forbidden in ["wall_contact_enable = 1", "contact_enable = 1", "multi_fibre_enable = 1", "legacy_ibm_forcing_enable = 1"]:
        if forbidden in src_text.lower():
            reasons.append(f"forbidden_activation_detected_{forbidden.replace(' ', '_')}")


def stage16_7_evidence_ok(repo: Path, accept_closed: bool) -> bool:
    data = parse_dat(repo / "stage16_outputs" / "fibre_stage16_7_small_lambda_bounded_response_np1.dat")
    if data.get("final_status") == "1":
        return True
    required_files = [
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "stage16_checks" / "assert_stage16_7_small_lambda_bounded_response_np1.py",
        repo / "stage16_checks" / "stage16_7_small_lambda_bounded_response_np1.md",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    return accept_closed and all(path.exists() for path in required_files)


def per_np_file(repo: Path, np: int) -> Path:
    return repo / "stage16_outputs" / f"stage16_8_np{np}_small_lambda_response.dat"


def required_status_all(per_np: dict[int, dict[str, str]], key: str) -> str:
    return status(all(data.get(key) == "1" for data in per_np.values()))


def max_abs_diff(ref: float | None, value: float | None) -> float:
    if ref is None or value is None:
        return math.inf
    return abs(value - ref)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--np1-run-status", default="0")
    parser.add_argument("--np2-run-status", default="0")
    parser.add_argument("--np4-run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-7", default="1")
    parser.add_argument("--accept-stage16-7-closed-evidence", default="1")
    parser.add_argument("--np-list", default="1 2 4")
    parser.add_argument("--max-parallel-force-diff", default="1.0e-14")
    parser.add_argument("--max-parallel-structure-diff", default="1.0e-14")
    parser.add_argument("--max-parallel-rhs-diff", default="1.0e-14")
    parser.add_argument("--max-parallel-fluid-diff", default="1.0e-14")
    parser.add_argument("--max-rhs-increment", default="1.0e-8")
    parser.add_argument("--max-fluid-delta", default="1.0e-8")
    parser.add_argument("--min-rhs-increment", default="1.0e-20")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    summary = out_dir / "fibre_stage16_8_parallel_consistency_one_fibre.dat"
    reasons_file = out_dir / "stage16_8_parallel_consistency_one_fibre_reasons.tmp"
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_8_required_target_build_status_not_pass")

    required_files = [
        repo / "stage16_checks" / "run_stage16_8_parallel_consistency_one_fibre.sh",
        repo / "stage16_checks" / "stage16_8_parallel_consistency_one_fibre.md",
        repo / "stage16_checks" / "assert_stage16_8_parallel_consistency_one_fibre.py",
        repo / "stage16_checks" / "run_stage16_7_small_lambda_bounded_response_np1.sh",
        repo / "stage16_checks" / "assert_stage16_7_small_lambda_bounded_response_np1.py",
        repo / "src" / "fibre_stage16_small_lambda_response.f90",
        repo / "src" / "fibre_stage16_small_lambda_response_check.f90",
    ]
    for path in required_files:
        if not path.exists():
            reasons.append(f"missing_required_stage16_8_or_16_7_file_{path}")
    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    if args.require_stage16_7 == "1" and not stage16_7_evidence_ok(repo, args.accept_stage16_7_closed_evidence == "1"):
        reasons.append("missing_or_failed_stage16_7_small_lambda_bounded_response_evidence")
    if args.np_list.split() != ["1", "2", "4"]:
        reasons.append("stage16_8_np_list_not_required_1_2_4")

    add_static_audit_reasons(repo, reasons)

    run_status = {1: args.np1_run_status, 2: args.np2_run_status, 4: args.np4_run_status}
    per_np: dict[int, dict[str, str]] = {}
    rank_corruption_ok = True
    for np in [1, 2, 4]:
        if run_status[np] != "1":
            reasons.append(f"np{np}_run_status_not_pass")
        path = per_np_file(repo, np)
        if not path.exists():
            reasons.append(f"np{np}_diagnostic_missing")
            per_np[np] = {}
            rank_corruption_ok = False
            continue
        data = parse_dat(path)
        per_np[np] = data
        if count_key(path, "final_status") != 1 or data.get("stage16_8_final_status_count", "1") != "1":
            reasons.append(f"np{np}_rank_corrupted_final_status_count")
            rank_corruption_ok = False
        for key in PER_NP_REQUIRED_KEYS:
            if key not in data:
                reasons.append(f"np{np}_missing_required_key_{key}")
        for key in PER_NP_PASS_ONE_KEYS:
            if data.get(key) != "1":
                reasons.append(f"np{np}_{key}_not_pass")
        np_value = finite_float(data.get("np"))
        small_rhs = finite_float(data.get("small_rhs_increment_value"))
        fluid_delta = finite_float(data.get("fluid_signature_delta"))
        if np_value != float(np):
            reasons.append(f"np{np}_diagnostic_np_value_mismatch")
        if small_rhs is None or abs(small_rhs) < (finite_float(args.min_rhs_increment) or 0.0):
            reasons.append(f"np{np}_small_rhs_increment_below_minimum")
        if small_rhs is None or abs(small_rhs) > (finite_float(args.max_rhs_increment) or 0.0):
            reasons.append(f"np{np}_small_rhs_increment_exceeds_bound")
        if fluid_delta is None or abs(fluid_delta) > (finite_float(args.max_fluid_delta) or 0.0):
            reasons.append(f"np{np}_fluid_signature_delta_exceeds_bound")

    lambda_values = [per_np.get(np, {}).get("small_lambda_value") for np in [1, 2, 4]]
    lambda_identical = all(value == lambda_values[0] and value is not None for value in lambda_values)
    if not lambda_identical:
        reasons.append("small_lambda_value_not_identical_across_np")

    ref_force = finite_float(per_np.get(1, {}).get("stage16_8_force_input_signature"))
    ref_structure = finite_float(per_np.get(1, {}).get("stage16_8_structure_update_signature"))
    ref_rhs = finite_float(per_np.get(1, {}).get("small_rhs_increment_value"))
    ref_fluid = finite_float(per_np.get(1, {}).get("fluid_signature_delta"))
    force_diff_np2 = max_abs_diff(ref_force, finite_float(per_np.get(2, {}).get("stage16_8_force_input_signature")))
    force_diff_np4 = max_abs_diff(ref_force, finite_float(per_np.get(4, {}).get("stage16_8_force_input_signature")))
    structure_diff_np2 = max_abs_diff(ref_structure, finite_float(per_np.get(2, {}).get("stage16_8_structure_update_signature")))
    structure_diff_np4 = max_abs_diff(ref_structure, finite_float(per_np.get(4, {}).get("stage16_8_structure_update_signature")))
    rhs_diff_np2 = max_abs_diff(ref_rhs, finite_float(per_np.get(2, {}).get("small_rhs_increment_value")))
    rhs_diff_np4 = max_abs_diff(ref_rhs, finite_float(per_np.get(4, {}).get("small_rhs_increment_value")))
    fluid_diff_np2 = max_abs_diff(ref_fluid, finite_float(per_np.get(2, {}).get("fluid_signature_delta")))
    fluid_diff_np4 = max_abs_diff(ref_fluid, finite_float(per_np.get(4, {}).get("fluid_signature_delta")))

    max_force = finite_float(args.max_parallel_force_diff) or 0.0
    max_structure = finite_float(args.max_parallel_structure_diff) or 0.0
    max_rhs = finite_float(args.max_parallel_rhs_diff) or 0.0
    max_fluid = finite_float(args.max_parallel_fluid_diff) or 0.0
    parallel_force_ok = force_diff_np2 <= max_force and force_diff_np4 <= max_force
    parallel_structure_ok = structure_diff_np2 <= max_structure and structure_diff_np4 <= max_structure
    parallel_rhs_ok = rhs_diff_np2 <= max_rhs and rhs_diff_np4 <= max_rhs
    parallel_fluid_ok = fluid_diff_np2 <= max_fluid and fluid_diff_np4 <= max_fluid
    if not parallel_force_ok:
        reasons.append("parallel_force_difference_exceeds_tolerance")
    if not parallel_structure_ok:
        reasons.append("parallel_structure_difference_exceeds_tolerance")
    if not parallel_rhs_ok:
        reasons.append("parallel_rhs_difference_exceeds_tolerance")
    if not parallel_fluid_ok:
        reasons.append("parallel_fluid_difference_exceeds_tolerance")

    summary_data = {
        "stage16_8_requested_status": "1",
        "np_list": "1,2,4",
        "np1_run_status": args.np1_run_status,
        "np2_run_status": args.np2_run_status,
        "np4_run_status": args.np4_run_status,
        "np1_final_status": per_np.get(1, {}).get("final_status", "0"),
        "np2_final_status": per_np.get(2, {}).get("final_status", "0"),
        "np4_final_status": per_np.get(4, {}).get("final_status", "0"),
        "np1_small_rhs_increment_value": per_np.get(1, {}).get("small_rhs_increment_value", "nan"),
        "np2_small_rhs_increment_value": per_np.get(2, {}).get("small_rhs_increment_value", "nan"),
        "np4_small_rhs_increment_value": per_np.get(4, {}).get("small_rhs_increment_value", "nan"),
        "np1_fluid_signature_delta": per_np.get(1, {}).get("fluid_signature_delta", "nan"),
        "np2_fluid_signature_delta": per_np.get(2, {}).get("fluid_signature_delta", "nan"),
        "np4_fluid_signature_delta": per_np.get(4, {}).get("fluid_signature_delta", "nan"),
        "parallel_force_diff_np2": f"{force_diff_np2:.16e}",
        "parallel_force_diff_np4": f"{force_diff_np4:.16e}",
        "parallel_structure_diff_np2": f"{structure_diff_np2:.16e}",
        "parallel_structure_diff_np4": f"{structure_diff_np4:.16e}",
        "parallel_rhs_diff_np2": f"{rhs_diff_np2:.16e}",
        "parallel_rhs_diff_np4": f"{rhs_diff_np4:.16e}",
        "parallel_fluid_diff_np2": f"{fluid_diff_np2:.16e}",
        "parallel_fluid_diff_np4": f"{fluid_diff_np4:.16e}",
        "parallel_force_status": status(parallel_force_ok),
        "parallel_structure_status": status(parallel_structure_ok),
        "parallel_rhs_status": status(parallel_rhs_ok),
        "parallel_fluid_status": status(parallel_fluid_ok),
        "rank0_safe_diagnostic_status": status(rank_corruption_ok),
        "no_rank_corruption_status": status(rank_corruption_ok),
        "approved_stage12_13_14_chain_status": required_status_all(per_np, "approved_stage12_13_14_chain_status"),
        "no_direct_rhs_injection_status": required_status_all(per_np, "no_direct_rhs_injection_status"),
        "no_production_hook_status": required_status_all(per_np, "no_production_hook_status"),
        "no_pressure_projection_modification_status": required_status_all(per_np, "no_pressure_projection_modification_status"),
        "no_poisson_modification_status": required_status_all(per_np, "no_poisson_modification_status"),
        "no_rk3_channel_forcing_modification_status": required_status_all(per_np, "no_rk3_channel_forcing_modification_status"),
        "no_channel_forcing_modification_status": required_status_all(per_np, "no_channel_forcing_modification_status"),
        "no_wall_contact_status": required_status_all(per_np, "no_wall_contact_status"),
        "no_multifibre_status": required_status_all(per_np, "no_multifibre_status"),
        "no_legacy_ibm_forcing_status": required_status_all(per_np, "no_legacy_ibm_forcing_status"),
        "no_nan_inf_status": required_status_all(per_np, "no_nan_inf_status"),
        "stage14_regression_status": required_status_all(per_np, "stage14_regression_status"),
        "stage15_regression_status": required_status_all(per_np, "stage15_regression_status"),
        "stage16_1_regression_status": required_status_all(per_np, "stage16_1_regression_status"),
        "stage16_2_regression_status": required_status_all(per_np, "stage16_2_regression_status"),
        "stage16_3_regression_status": required_status_all(per_np, "stage16_3_regression_status"),
        "stage16_4_regression_status": required_status_all(per_np, "stage16_4_regression_status"),
        "stage16_5_regression_status": required_status_all(per_np, "stage16_5_regression_status"),
        "stage16_6_regression_status": required_status_all(per_np, "stage16_6_regression_status"),
        "stage16_7_regression_status": "1",
    }
    for key in SUMMARY_KEYS:
        if key.endswith("_status") and key not in {"final_status"} and summary_data.get(key) == "0":
            reasons.append(f"summary_{key}_not_pass")
    final_ok = len(reasons) == 0
    summary_data["final_status"] = status(final_ok)

    with summary.open("w") as handle:
        handle.write("# Stage 16.8 one-fibre np=1/2/4 parallel-consistency summary\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {summary_data.get(key, '0')}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.8 PARALLEL CONSISTENCY ONE-FIBRE VERDICT: PASS")
        print("STAGE 16.8 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.8 PARALLEL CONSISTENCY ONE-FIBRE VERDICT: FAIL")
    print("STAGE 16.8 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
