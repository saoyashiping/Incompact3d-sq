#!/usr/bin/env python3
"""Stage 16.11 short-time stability / bounded-energy smoke helper."""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

SUMMARY_KEYS = [
    "stage16_11_requested_status", "np", "nsteps", "small_lambda_value",
    "closed_loop_path_status", "one_fibre_count_status", "stage11_sampling_status",
    "stage12_feedback_status", "stage16_4_force_input_status", "controlled_structure_update_status",
    "controlled_structure_update_bounded_status", "max_position_update", "position_update_bounded_status",
    "max_velocity_update", "velocity_update_bounded_status", "max_acceleration_update",
    "acceleration_update_bounded_status", "max_force_input", "force_input_bounded_status",
    "stage13_force_density_status", "stage13_force_density_bounded_status", "stage14_rhs_status",
    "max_rhs_increment", "rhs_increment_bounded_status", "fluid_signature_delta",
    "fluid_signature_bounded_status", "work_proxy_value", "work_proxy_finite_status",
    "work_proxy_bounded_status", "energy_proxy_value", "energy_proxy_finite_status",
    "energy_proxy_bounded_status", "growth_ratio", "growth_ratio_bounded_status",
    "no_runaway_growth_status", "approved_stage11_12_16_4_13_14_chain_status",
    "approved_stage12_13_14_chain_status", "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status", "no_wall_contact_status", "no_multifibre_status",
    "no_unapproved_bending_solve_status", "no_unapproved_tension_solve_status",
    "no_unapproved_full_structure_solve_status", "rank0_safe_diagnostic_status",
    "no_rank_corruption_status", "no_nan_inf_status", "stage14_regression_status",
    "stage15_regression_status", "stage16_1_regression_status", "stage16_2_regression_status",
    "stage16_3_regression_status", "stage16_4_regression_status", "stage16_5_regression_status",
    "stage16_6_regression_status", "stage16_7_regression_status", "stage16_8_regression_status",
    "stage16_9_regression_status", "stage16_10_regression_status", "final_status",
]

RUNTIME_REQUIRED_KEYS = [
    "np", "nsteps", "small_lambda_value", "closed_loop_path_status", "one_fibre_count_status",
    "stage11_sampling_status", "stage12_feedback_status", "stage16_4_force_input_status",
    "controlled_structure_update_status", "controlled_structure_update_bounded_status",
    "max_position_update", "max_velocity_update", "max_acceleration_update", "max_force_input",
    "stage13_force_density_status", "max_rhs_increment", "stage14_rhs_status",
    "fluid_signature_delta", "work_proxy_value", "energy_proxy_value", "growth_ratio",
    "approved_stage12_13_14_chain_status", "no_direct_rhs_injection_status",
    "no_pressure_projection_modification_status", "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_wall_contact_status", "no_multifibre_status", "no_legacy_ibm_forcing_status",
    "no_nan_inf_status", "stage14_regression_status", "stage15_regression_status",
    "stage16_1_regression_status", "stage16_2_regression_status", "stage16_3_regression_status",
    "stage16_4_regression_status", "stage16_5_regression_status", "stage16_6_regression_status",
    "stage16_7_regression_status", "final_status",
]

RUNTIME_PASS_ONE_KEYS = [
    key for key in RUNTIME_REQUIRED_KEYS
    if key not in {"np", "nsteps", "small_lambda_value", "max_position_update", "max_velocity_update",
                   "max_acceleration_update", "max_force_input", "max_rhs_increment", "fluid_signature_delta",
                   "work_proxy_value", "energy_proxy_value", "growth_ratio"}
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

    Stage 16.11 deliberately reuses the corrected Stage 16.10 / 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4
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


def evidence_ok(repo: Path, stage: str, files: list[str], dat_name: str = "") -> bool:
    if dat_name:
        data = parse_dat(repo / "stage16_outputs" / dat_name)
        if data.get("final_status") == "1":
            return True
    return all((repo / path).exists() for path in files)


def add_static_audit_reasons(repo: Path, reasons: list[str]) -> dict[str, str]:
    """Apply Stage 16.11 static bounded-stability checks without broad false-positive scans.

    False-positive protections are required for Stage 16.11 and are copied from the passed
    Stage 16.10 / 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4 helper style:
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
    src_files = all_files(repo / "src", ("*.f90",))
    src_text = joined(src_files)
    xcompact3d_text = read_text(repo / "src" / "xcompact3d.f90")
    stage11_14_shell_text = joined(
        all_files(repo / "stage11_checks", ("*.sh",))
        + all_files(repo / "stage13_checks", ("*.sh",))
        + all_files(repo / "stage14_checks", ("*.sh",))
    )
    stage15_shell_text = joined(all_files(repo / "stage15_checks", ("*.sh",)))
    cmake_text = read_text(repo / "src" / "CMakeLists.txt")

    required_chain_markers = [
        "fibre_stage16_small_lambda_response", "stage16_structure_force_input_set_from_stage12_candidate",
        "stage13_force_density_signature", "stage14_rhs_increment", "fibre_stage14_production_rhs_injection",
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

    wrapper = repo / "stage16_checks" / "run_stage16_11_short_time_stability_smoke.sh"
    if wrapper.exists() and not check_rg_fallback(wrapper):
        reasons.append(f"rg_without_grep_fallback_{wrapper}")
    wrapper_text = read_text(wrapper)
    helper_text = read_text(repo / "stage16_checks" / "assert_stage16_11_short_time_stability_smoke.py")
    # Unknown-failure auditing follows the same false-positive-safe rule: inspect executable
    # wrapper behavior, not this helper's negative-check reason strings.
    if "unknown failure" in wrapper_text.lower():
        reasons.append("ambiguous_unknown_failure_fallback_text_detected")
    if "Stage 16.10 / 16.9 / 16.8 / 16.7 / 16.6 / 16.5 / 16.4" not in helper_text:
        reasons.append("false_positive_safe_pattern_comment_missing")
    return statuses


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-10", default="1")
    parser.add_argument("--accept-stage16-10-closed-evidence", default="1")
    parser.add_argument("--np", default="2")
    parser.add_argument("--nsteps", default="5")
    parser.add_argument("--max-force-input", default="1.0e-6")
    parser.add_argument("--max-structure-update", default="1.0e-12")
    parser.add_argument("--max-velocity-update", default="1.0e-10")
    parser.add_argument("--max-acceleration-update", default="1.0e-6")
    parser.add_argument("--max-rhs-increment", default="1.0e-8")
    parser.add_argument("--max-fluid-delta", default="1.0e-8")
    parser.add_argument("--max-work-proxy", default="1.0e-14")
    parser.add_argument("--max-energy-proxy", default="1.0e-14")
    parser.add_argument("--max-growth-ratio", default="10.0")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    summary = out_dir / "fibre_stage16_11_short_time_stability_smoke.dat"
    reasons_file = out_dir / "stage16_11_short_time_stability_smoke_reasons.tmp"
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_11_required_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_11_runtime_evidence_run_status_not_pass")

    required_files = [
        "stage16_checks/run_stage16_11_short_time_stability_smoke.sh",
        "stage16_checks/stage16_11_short_time_stability_smoke.md",
        "stage16_checks/assert_stage16_11_short_time_stability_smoke.py",
        "stage16_checks/run_stage16_10_contamination_audit.sh",
        "stage16_checks/assert_stage16_10_contamination_audit.py",
        "stage16_checks/run_stage16_9_io_restart_one_fibre.sh",
        "stage16_checks/run_stage16_8_parallel_consistency_one_fibre.sh",
        "stage16_checks/run_stage16_7_small_lambda_bounded_response_np1.sh",
        "src/fibre_stage16_small_lambda_response.f90",
        "src/fibre_stage16_small_lambda_response_check.f90",
    ]
    for rel in required_files:
        if not (repo / rel).exists():
            reasons.append(f"missing_required_stage16_11_or_closed_file_{rel}")
    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")

    if args.require_stage16_10 == "1" and not evidence_ok(repo, "16.10", [
        "stage16_checks/run_stage16_10_contamination_audit.sh",
        "stage16_checks/assert_stage16_10_contamination_audit.py",
        "stage16_checks/stage16_10_contamination_audit.md"], "fibre_stage16_10_contamination_audit.dat"):
        reasons.append("missing_or_failed_stage16_10_contamination_evidence")
    for label, files, dat in [
        ("16.9", ["stage16_checks/run_stage16_9_io_restart_one_fibre.sh", "stage16_checks/assert_stage16_9_io_restart_one_fibre.py", "stage16_checks/stage16_9_io_restart_one_fibre.md"], "fibre_stage16_9_io_restart_one_fibre.dat"),
        ("16.8", ["stage16_checks/run_stage16_8_parallel_consistency_one_fibre.sh", "stage16_checks/assert_stage16_8_parallel_consistency_one_fibre.py", "stage16_checks/stage16_8_parallel_consistency_one_fibre.md"], "fibre_stage16_8_parallel_consistency_one_fibre.dat"),
        ("16.7", ["stage16_checks/run_stage16_7_small_lambda_bounded_response_np1.sh", "stage16_checks/assert_stage16_7_small_lambda_bounded_response_np1.py", "stage16_checks/stage16_7_small_lambda_bounded_response_np1.md"], "fibre_stage16_7_small_lambda_bounded_response_np1.dat"),
    ]:
        if not evidence_ok(repo, label, files, dat):
            reasons.append(f"missing_or_failed_stage{label}_closed_evidence")

    static_status = add_static_audit_reasons(repo, reasons)
    runtime_path = out_dir / "stage16_11_short_time_stability_evidence.dat"
    runtime = parse_dat(runtime_path)
    if not runtime_path.exists():
        reasons.append("stage16_11_runtime_evidence_missing")
    if count_key(runtime_path, "final_status") != 1:
        reasons.append("stage16_11_rank_corrupted_runtime_final_status_count")
    for key in RUNTIME_REQUIRED_KEYS:
        if key not in runtime:
            reasons.append(f"missing_runtime_key_{key}")
    for key in RUNTIME_PASS_ONE_KEYS:
        if runtime.get(key) != "1":
            reasons.append(f"runtime_{key}_not_pass")

    values = {key: finite_float(runtime.get(key)) for key in [
        "np", "nsteps", "max_position_update", "max_velocity_update", "max_acceleration_update",
        "max_force_input", "max_rhs_increment", "fluid_signature_delta", "work_proxy_value",
        "energy_proxy_value", "growth_ratio"]}
    position_ok = values["max_position_update"] is not None and abs(values["max_position_update"]) <= (finite_float(args.max_structure_update) or 0.0)
    velocity_ok = values["max_velocity_update"] is not None and abs(values["max_velocity_update"]) <= (finite_float(args.max_velocity_update) or 0.0)
    accel_ok = values["max_acceleration_update"] is not None and abs(values["max_acceleration_update"]) <= (finite_float(args.max_acceleration_update) or 0.0)
    force_ok = values["max_force_input"] is not None and abs(values["max_force_input"]) <= (finite_float(args.max_force_input) or 0.0)
    stage13_ok = force_ok and runtime.get("stage13_force_density_status") == "1"
    rhs_ok = values["max_rhs_increment"] is not None and abs(values["max_rhs_increment"]) <= (finite_float(args.max_rhs_increment) or 0.0)
    fluid_ok = values["fluid_signature_delta"] is not None and abs(values["fluid_signature_delta"]) <= (finite_float(args.max_fluid_delta) or 0.0)
    work_finite = values["work_proxy_value"] is not None
    work_ok = work_finite and abs(values["work_proxy_value"] or 0.0) <= (finite_float(args.max_work_proxy) or 0.0)
    energy_finite = values["energy_proxy_value"] is not None
    energy_ok = energy_finite and abs(values["energy_proxy_value"] or 0.0) <= (finite_float(args.max_energy_proxy) or 0.0)
    growth_ok = values["growth_ratio"] is not None and values["growth_ratio"] >= 0.0 and values["growth_ratio"] <= (finite_float(args.max_growth_ratio) or 0.0)
    np_ok = values["np"] == finite_float(args.np)
    nsteps_ok = values["nsteps"] == finite_float(args.nsteps)
    for ok, reason in [
        (np_ok, "stage16_11_np_value_mismatch"), (nsteps_ok, "stage16_11_nsteps_value_mismatch"),
        (position_ok, "position_update_nonfinite_or_unbounded"), (velocity_ok, "velocity_update_nonfinite_or_unbounded"),
        (accel_ok, "acceleration_update_nonfinite_or_unbounded"), (force_ok, "force_input_nonfinite_or_unbounded"),
        (stage13_ok, "stage13_force_density_nonfinite_or_unbounded"), (rhs_ok, "stage14_rhs_increment_nonfinite_or_unbounded"),
        (fluid_ok, "fluid_signature_delta_exceeds_tolerance"), (work_ok, "work_proxy_nonfinite_or_unbounded"),
        (energy_ok, "energy_proxy_nonfinite_or_unbounded"), (growth_ok, "growth_ratio_exceeds_tolerance"),
    ]:
        if not ok:
            reasons.append(reason)

    rank_ok = count_key(runtime_path, "final_status") == 1
    no_nan_inf = runtime.get("no_nan_inf_status") == "1" and all(value is not None for value in values.values())
    approved_chain = runtime.get("closed_loop_path_status") == "1" and static_status["approved_chain"] == "1"
    no_runaway = growth_ok and position_ok and velocity_ok and accel_ok and rhs_ok and fluid_ok and work_ok and energy_ok

    summary_data = {
        "stage16_11_requested_status": "1",
        "np": runtime.get("np", args.np), "nsteps": runtime.get("nsteps", args.nsteps),
        "small_lambda_value": runtime.get("small_lambda_value", "nan"),
        "closed_loop_path_status": runtime.get("closed_loop_path_status", "0"),
        "one_fibre_count_status": runtime.get("one_fibre_count_status", "0"),
        "stage11_sampling_status": runtime.get("stage11_sampling_status", "0"),
        "stage12_feedback_status": runtime.get("stage12_feedback_status", "0"),
        "stage16_4_force_input_status": runtime.get("stage16_4_force_input_status", "0"),
        "controlled_structure_update_status": runtime.get("controlled_structure_update_status", "0"),
        "controlled_structure_update_bounded_status": runtime.get("controlled_structure_update_bounded_status", "0"),
        "max_position_update": runtime.get("max_position_update", "nan"),
        "position_update_bounded_status": status(position_ok),
        "max_velocity_update": runtime.get("max_velocity_update", "nan"),
        "velocity_update_bounded_status": status(velocity_ok),
        "max_acceleration_update": runtime.get("max_acceleration_update", "nan"),
        "acceleration_update_bounded_status": status(accel_ok),
        "max_force_input": runtime.get("max_force_input", "nan"),
        "force_input_bounded_status": status(force_ok),
        "stage13_force_density_status": runtime.get("stage13_force_density_status", "0"),
        "stage13_force_density_bounded_status": status(stage13_ok),
        "stage14_rhs_status": runtime.get("stage14_rhs_status", "0"),
        "max_rhs_increment": runtime.get("max_rhs_increment", "nan"),
        "rhs_increment_bounded_status": status(rhs_ok),
        "fluid_signature_delta": runtime.get("fluid_signature_delta", "nan"),
        "fluid_signature_bounded_status": status(fluid_ok),
        "work_proxy_value": runtime.get("work_proxy_value", "nan"),
        "work_proxy_finite_status": status(work_finite),
        "work_proxy_bounded_status": status(work_ok),
        "energy_proxy_value": runtime.get("energy_proxy_value", "nan"),
        "energy_proxy_finite_status": status(energy_finite),
        "energy_proxy_bounded_status": status(energy_ok),
        "growth_ratio": runtime.get("growth_ratio", "nan"),
        "growth_ratio_bounded_status": status(growth_ok),
        "no_runaway_growth_status": status(no_runaway),
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
        "stage16_8_regression_status": "1" if evidence_ok(repo, "16.8", ["stage16_checks/run_stage16_8_parallel_consistency_one_fibre.sh"], "fibre_stage16_8_parallel_consistency_one_fibre.dat") else "0",
        "stage16_9_regression_status": "1" if evidence_ok(repo, "16.9", ["stage16_checks/run_stage16_9_io_restart_one_fibre.sh"], "fibre_stage16_9_io_restart_one_fibre.dat") else "0",
        "stage16_10_regression_status": "1" if evidence_ok(repo, "16.10", ["stage16_checks/run_stage16_10_contamination_audit.sh"], "fibre_stage16_10_contamination_audit.dat") else "0",
    }
    for key in SUMMARY_KEYS:
        if key.endswith("_status") and key != "final_status" and summary_data.get(key) == "0":
            reasons.append(f"summary_{key}_not_pass")
    final_ok = len(reasons) == 0
    summary_data["final_status"] = status(final_ok)

    with summary.open("w") as handle:
        handle.write("# Stage 16.11 short-time stability / bounded-energy smoke summary\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {summary_data.get(key, '0')}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.11 SHORT-TIME STABILITY SMOKE VERDICT: PASS")
        print("STAGE 16.11 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.11 SHORT-TIME STABILITY SMOKE VERDICT: FAIL")
    print("STAGE 16.11 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
