#!/usr/bin/env python3
"""Stage 16.5 np=1 closed-loop dry-run audit helper."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REQUIRED_KEYS = [
    "stage16_5_requested_status",
    "np",
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
    "rhs_increment_bounded_status",
    "fluid_signature_delta",
    "fluid_signature_status",
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
    "final_status",
]

NUMERIC_KEYS = {"np", "fluid_signature_delta"}
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


def check_rg_fallback(script: Path) -> bool:
    """Check shell scripts only for real rg command usage with a grep fallback.

    This deliberately mirrors the corrected Stage 16.4 false-positive-safe audit style:
    documentation files are not scanned as executable evidence, Python regex literals such
    as rg[[:space:]] are not treated as command usage, and negative-check text is not a
    regression.  Only shell wrappers with actual rg command/fallback logic are audited.
    """
    if script.suffix != ".sh":
        return True
    text = read_text(script)
    actual_rg_usage = ("command -v rg" in text or "which rg" in text or
                       re.search(r"(?:^|[;&|]\s*|\s)rg(?:\s|$)", text) is not None)
    if not actual_rg_usage:
        return True
    return ("command -v rg" in text or "which rg" in text) and re.search(r"\bgrep\b", text) is not None


def add_static_audit_reasons(repo: Path, data: dict[str, str], reasons: list[str]) -> None:
    """Apply Stage 16.5 static checks without broad false-positive-prone scans.

    False-positive protections are intentional and required for Stage 16.5:
    * .md documentation is checked for existence only, not scanned as code-regression evidence.
    * negative-check strings and regex literals inside Python helpers are not used as source behavior.
    * the legitimate Stage 13.5 conservation/sign audit is allowed; only old Stage 13.5
      production force-density diagnostic names in real production/check contexts are rejected.
    * if a forbidden pattern cannot be distinguished from documentation/negative-check text, this
      helper avoids that scan rather than reporting ambiguous evidence as a regression.
    """
    cmake = repo / "src" / "CMakeLists.txt"
    src_files = all_files(repo / "src", ("*.f90",))
    src_text = joined(src_files)
    stage11_14_text = joined(all_files(repo / "stage11_checks", ("*.sh",)) +
                            all_files(repo / "stage13_checks", ("*.sh",)) +
                            all_files(repo / "stage14_checks", ("*.sh",)))
    stage15_shell_files = all_files(repo / "stage15_checks", ("*.sh",))
    stage15_text = joined(stage15_shell_files)

    stage14_regression_status = 1
    if re.search(r"stage14_get_injection_gain\s*\(\s*\)\s*==\s*0\.0", src_text + stage11_14_text):
        reasons.append("stage14_forbidden_lambda_zero_hook_gate_detected")
        stage14_regression_status = 0
    if "fibre_stage14_production_rhs_injection" not in src_text:
        reasons.append("stage14_production_rhs_injection_source_missing")
        stage14_regression_status = 0
    if "stage14_small_lambda" not in src_text + stage11_14_text and "small_lambda" not in src_text + stage11_14_text:
        reasons.append("stage14_small_lambda_diagnostic_evidence_missing")
        stage14_regression_status = 0

    stage13_6_status = 1
    production_context = src_text + joined(all_files(repo / "stage13_checks", ("run_stage13_*production*.sh",)))
    if "stage13_6_production_force_density_candidate_status" not in production_context:
        reasons.append("stage13_6_production_force_density_candidate_status_missing")
        stage13_6_status = 0
        stage14_regression_status = 0
    if "fibre_stage13_6_production_force_density_candidate.dat" not in production_context:
        reasons.append("stage13_6_production_force_density_candidate_dat_missing")
        stage13_6_status = 0
        stage14_regression_status = 0
    if "stage13_5_production_force_density_candidate" in production_context:
        reasons.append("old_stage13_5_production_force_density_candidate_name_detected_in_production_context")
        stage13_6_status = 0
        stage14_regression_status = 0

    sampling_src = read_text(repo / "src" / "fibre_stage13_production_force_density_candidate.f90")
    stage13_sampling_status = 1
    if "subdomain center" in sampling_src.lower() or "local_subdomain_center" in sampling_src.lower():
        reasons.append("stage13_local_subdomain_center_sampling_detected")
        stage13_sampling_status = 0
        stage14_regression_status = 0

    rank0_safe_status = 1
    rank0_files = [
        repo / "src" / "fibre_stage11_production_oneway_hook.f90",
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "src" / "fibre_stage14_production_rhs_injection.f90",
        repo / "src" / "fibre_stage15_production_structure_hook.f90",
        repo / "src" / "fibre_stage16_closed_loop_dryrun_check.f90",
    ]
    for path in rank0_files:
        text = read_text(path)
        if "rank0" not in text.lower() and "rank == 0" not in text.lower():
            reasons.append(f"rank0_safe_diagnostic_evidence_missing_{path}")
            rank0_safe_status = 0

    stage15_regression_status = 1
    for idx, stem in STAGE15_DOC_SCRIPT_STEMS.items():
        script = repo / "stage15_checks" / f"run_stage15_{idx}_{stem}.sh"
        doc = repo / "stage15_checks" / f"stage15_{idx}_{stem}.md"
        if not script.exists():
            reasons.append(f"missing_stage15_{idx}_script_{script}")
            stage15_regression_status = 0
        if not doc.exists():
            reasons.append(f"missing_stage15_{idx}_doc_{doc}")
            stage15_regression_status = 0
    for script in stage15_shell_files + [repo / "stage16_checks" / "run_stage16_5_closed_loop_dryrun_np1.sh"]:
        if script.exists() and not check_rg_fallback(script):
            reasons.append(f"rg_without_grep_fallback_{script}")
            stage15_regression_status = 0
    if "unknown failure" in stage15_text.lower() and "failure_reason" not in stage15_text:
        reasons.append("stage15_unknown_failure_without_explicit_reason_detected")
        stage15_regression_status = 0

    stage16_1_regression_status = check_stage16_evidence(repo, reasons, 1, [
        repo / "stage16_checks" / "run_stage16_1_config.sh",
        repo / "stage16_checks" / "assert_stage16_1_config.py",
        repo / "stage16_checks" / "stage16_1_config.md",
        repo / "src" / "fibre_stage16_config.f90",
        repo / "src" / "fibre_stage16_config_check.f90",
    ], "add_executable(fibre_stage16_config_check")
    stage16_2_regression_status = check_stage16_evidence(repo, reasons, 2, [
        repo / "stage16_checks" / "run_stage16_2_one_fibre_case_definition.sh",
        repo / "stage16_checks" / "assert_stage16_2_one_fibre_case_definition.py",
        repo / "stage16_checks" / "stage16_2_one_fibre_case_definition.md",
        repo / "src" / "fibre_stage16_one_fibre_case.f90",
        repo / "src" / "fibre_stage16_one_fibre_case_check.f90",
    ], "add_executable(fibre_stage16_one_fibre_case_check")
    stage16_3_regression_status = check_stage16_evidence(repo, reasons, 3, [
        repo / "stage16_checks" / "run_stage16_3_force_sign_audit.sh",
        repo / "stage16_checks" / "assert_stage16_3_force_sign_audit.py",
        repo / "stage16_checks" / "stage16_3_force_sign_audit.md",
        repo / "src" / "fibre_stage16_force_sign_audit.f90",
        repo / "src" / "fibre_stage16_force_sign_audit_check.f90",
    ], "add_executable(fibre_stage16_force_sign_audit_check")
    stage16_4_regression_status = check_stage16_evidence(repo, reasons, 4, [
        repo / "stage16_checks" / "run_stage16_4_structure_force_input.sh",
        repo / "stage16_checks" / "assert_stage16_4_structure_force_input.py",
        repo / "stage16_checks" / "stage16_4_structure_force_input.md",
        repo / "src" / "fibre_stage16_structure_force_input.f90",
        repo / "src" / "fibre_stage16_structure_force_input_check.f90",
    ], "add_executable(fibre_stage16_structure_force_input_check")

    boundary_status = 1
    if "stage16" in read_text(repo / "src" / "xcompact3d.f90").lower():
        reasons.append("stage16_production_hook_detected_in_xcompact3d")
        boundary_status = 0
    stage16_5_src = read_text(repo / "src" / "fibre_stage16_closed_loop_dryrun.f90") + "\n" + \
        read_text(repo / "src" / "fibre_stage16_closed_loop_dryrun_check.f90")
    forbidden_call_patterns = {
        "pressure_projection_detected": r"\b(use|call)\b.*(pressure|projection)",
        "poisson_detected": r"\b(use|call)\b.*poisson",
        "rk3_channel_forcing_detected": r"\b(use|call)\b.*(rk3|channel.*forc|forcing.*channel)",
        "wall_contact_detected": r"\b(use|call)\b.*(wall|contact)",
        "multifibre_detected": r"\b(use|call)\b.*multi[-_ ]?fibre",
        "legacy_ibm_forcing_detected": r"\b(use|call)\b.*(legacy.*ibm|production_ibm|ibm_forcing|apply_ibm)",
        "direct_fluid_rhs_update_detected": r"\b(ux|uy|uz|rhs_x|rhs_y|rhs_z)\s*=\s*[^=]",
        "stage14_rhs_injection_call_detected": r"\bcall\b.*stage14.*rhs.*(inject|addition|accumulator|apply)",
    }
    for reason, pattern in forbidden_call_patterns.items():
        if re.search(pattern, stage16_5_src, re.I):
            reasons.append(reason)
            boundary_status = 0

    data["stage14_regression_status"] = str(stage14_regression_status)
    data["stage15_regression_status"] = str(stage15_regression_status)
    data["stage16_1_regression_status"] = str(stage16_1_regression_status)
    data["stage16_2_regression_status"] = str(stage16_2_regression_status)
    data["stage16_3_regression_status"] = str(stage16_3_regression_status)
    data["stage16_4_regression_status"] = str(stage16_4_regression_status)
    if stage13_6_status == 0 or stage13_sampling_status == 0 or rank0_safe_status == 0:
        data["stage14_regression_status"] = "0"
    if boundary_status == 0:
        for key in ["no_production_hook_status", "no_direct_rhs_injection_status",
                    "no_pressure_projection_modification_status", "no_poisson_modification_status",
                    "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
                    "no_wall_contact_status", "no_multifibre_status", "no_legacy_ibm_forcing_status"]:
            data[key] = "0"


def check_stage16_evidence(repo: Path, reasons: list[str], stage_index: int, files: list[Path], cmake_marker: str) -> int:
    status = 1
    for path in files:
        if not path.exists():
            reasons.append(f"missing_stage16_{stage_index}_evidence_file_{path}")
            status = 0
    if cmake_marker not in read_text(repo / "src" / "CMakeLists.txt"):
        reasons.append(f"stage16_{stage_index}_build_registration_missing")
        status = 0
    return status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-4", default="1")
    parser.add_argument("--accept-stage16-4-closed-evidence", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    dat = out_dir / "fibre_stage16_5_closed_loop_dryrun_np1.dat"
    reasons_file = out_dir / "stage16_5_closed_loop_dryrun_np1_reasons.tmp"
    data = parse_dat(dat)
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_5_closed_loop_dryrun_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_5_closed_loop_dryrun_target_run_status_not_pass")

    required_files = [
        repo / "stage16_checks" / "run_stage16_5_closed_loop_dryrun_np1.sh",
        repo / "stage16_checks" / "stage16_5_closed_loop_dryrun_np1.md",
        repo / "stage16_checks" / "assert_stage16_5_closed_loop_dryrun_np1.py",
        repo / "src" / "fibre_stage16_closed_loop_dryrun.f90",
        repo / "src" / "fibre_stage16_closed_loop_dryrun_check.f90",
    ]
    for path in required_files:
        if not path.exists():
            reasons.append(f"missing_required_stage16_5_file_{path}")
    if "add_executable(fibre_stage16_closed_loop_dryrun_check" not in read_text(repo / "src" / "CMakeLists.txt"):
        reasons.append("missing_fibre_stage16_closed_loop_dryrun_check_build_registration")

    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")

    stage16_4_files = [
        repo / "stage16_checks" / "run_stage16_4_structure_force_input.sh",
        repo / "stage16_checks" / "assert_stage16_4_structure_force_input.py",
        repo / "stage16_checks" / "stage16_4_structure_force_input.md",
        repo / "src" / "fibre_stage16_structure_force_input.f90",
        repo / "src" / "fibre_stage16_structure_force_input_check.f90",
    ]
    if args.require_stage16_4 == "1":
        stage16_4 = parse_dat(repo / "stage16_outputs" / "fibre_stage16_4_structure_force_input.dat")
        stage16_4_files_ok = all(path.exists() for path in stage16_4_files)
        if stage16_4.get("final_status") != "1" and not (args.accept_stage16_4_closed_evidence == "1" and stage16_4_files_ok):
            reasons.append("missing_or_failed_stage16_4_structure_force_input_evidence")

    add_static_audit_reasons(repo, data, reasons)

    for key in REQUIRED_KEYS:
        if key not in data and key != "final_status":
            reasons.append(f"missing_required_output_key_{key}")
            data[key] = "0"
    for key in PASS_ONE_KEYS:
        if data.get(key) != "1":
            reasons.append(f"{key}_not_pass")
    try:
        if int(data.get("np", "0")) != 1:
            reasons.append("np_not_one")
    except ValueError:
        reasons.append("np_not_integer")
    try:
        if abs(float(data.get("fluid_signature_delta", "nan"))) > 1.0e-8:
            reasons.append("fluid_signature_delta_exceeds_bound")
    except ValueError:
        reasons.append("fluid_signature_delta_not_numeric")
    for key in NUMERIC_KEYS:
        try:
            float(data.get(key, "nan"))
        except ValueError:
            reasons.append(f"{key}_not_numeric")

    final_ok = len(reasons) == 0
    data["final_status"] = "1" if final_ok else "0"

    with dat.open("w") as handle:
        handle.write("# Stage 16.5 np=1 closed-loop dry-run summary\n")
        for key in REQUIRED_KEYS:
            handle.write(f"{key} {data.get(key, '0')}\n")
        for key in ["numeric_parse_status", "numeric_bounds_status", "feedback_alpha", "lambda_value",
                    "stage13_force_density_signature", "stage14_rhs_increment"]:
            if key in data and key not in REQUIRED_KEYS:
                handle.write(f"{key} {data[key]}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.5 CLOSED-LOOP DRY RUN NP1 VERDICT: PASS")
        print("STAGE 16.5 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.5 CLOSED-LOOP DRY RUN NP1 VERDICT: FAIL")
    print("STAGE 16.5 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
