#!/usr/bin/env python3
"""Stage 16.1 configuration/global-switch audit helper."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REQUIRED_KEYS = [
    "stage16_1_requested_status",
    "default_off_status",
    "env_override_status",
    "master_enable_status",
    "one_fibre_fsi_enable_status",
    "structure_advance_enable_status",
    "two_way_rhs_enable_status",
    "diagnostic_only_status",
    "feedback_alpha",
    "lambda_value",
    "max_structure_update",
    "max_rhs_increment",
    "invalid_numeric_rejection_status",
    "invalid_flag_combination_rejection_status",
    "no_production_hook_status",
    "no_structure_advance_status",
    "no_rhs_modification_status",
    "no_bending_solve_status",
    "no_tension_solve_status",
    "no_wall_contact_status",
    "no_multifibre_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "stage14_regression_status",
    "stage15_regression_status",
    "final_status",
]

PASS_ONE_KEYS = [key for key in REQUIRED_KEYS if key not in {
    "feedback_alpha", "lambda_value", "max_structure_update", "max_rhs_increment",
    "master_enable_status", "one_fibre_fsi_enable_status",
    "structure_advance_enable_status", "two_way_rhs_enable_status", "final_status"
}]

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
    text = read_text(script)
    if not re.search(r"\brg\b", text):
        return True
    return ("command -v rg" in text or "which rg" in text) and re.search(r"\bgrep\b", text) is not None


def value_is_one(data: dict[str, str], key: str) -> bool:
    return data.get(key) == "1"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-0", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    dat = out_dir / "fibre_stage16_1_config.dat"
    reasons_file = out_dir / "stage16_1_config_reasons.tmp"
    reasons: list[str] = []
    data = parse_dat(dat)

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_1_config_check_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_1_config_check_target_run_status_not_pass")

    wrapper = repo / "stage16_checks" / "run_stage16_1_config.sh"
    doc = repo / "stage16_checks" / "stage16_1_config.md"
    helper = repo / "stage16_checks" / "assert_stage16_1_config.py"
    config_src = repo / "src" / "fibre_stage16_config.f90"
    check_src = repo / "src" / "fibre_stage16_config_check.f90"
    cmake = repo / "src" / "CMakeLists.txt"
    if not wrapper.exists() or not doc.exists() or not helper.exists() or not config_src.exists() or not check_src.exists():
        reasons.append("missing_stage16_1_wrapper_doc_helper_or_source")
    if "add_executable(fibre_stage16_config_check" not in read_text(cmake):
        reasons.append("missing_fibre_stage16_config_check_build_registration")

    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    if args.require_stage16_0 == "1":
        stage16_0_dat = repo / "stage16_outputs" / "fibre_stage16_0_preflight_closure_integrity.dat"
        stage16_0 = parse_dat(stage16_0_dat)
        if stage16_0.get("final_status") != "1":
            reasons.append("missing_or_failed_stage16_0_preflight_evidence")

    src_text = joined(all_files(repo / "src", ("*.f90",)))
    stage11_14_text = joined(all_files(repo / "stage11_checks", ("*.sh", "*.md")) +
                            all_files(repo / "stage13_checks", ("*.sh", "*.md")) +
                            all_files(repo / "stage14_checks", ("*.sh", "*.md")))
    stage15_text = joined(all_files(repo / "stage15_checks", ("*.sh", "*.md")))

    stage14_regression_status = 1
    if re.search(r"stage14_get_injection_gain\s*\(\s*\)\s*(?:==|\.eq\.)\s*0\.0", src_text + stage11_14_text, re.I):
        reasons.append("stage14_forbidden_lambda_zero_registration_gate_found")
        stage14_regression_status = 0
    for label, path, needle in [
        ("stage11_5_production_oneway_hook", repo / "stage11_checks" / "run_stage11_5_production_oneway_hook.sh", "stage11_5_production_oneway_hook"),
        ("stage13_6_production_force_density_candidate", repo / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh", "stage13_6_production_force_density_candidate"),
        ("stage14_5_production_rhs_hook", repo / "stage14_checks" / "run_stage14_5_production_rhs_hook.sh", "stage14_5_production_rhs_hook"),
    ]:
        if not path.exists() or needle not in src_text + stage11_14_text:
            reasons.append(f"missing_{label}_diagnostics")
            stage14_regression_status = 0
    if "stage13_5_production_force_density_candidate" in src_text + stage11_14_text + stage15_text:
        reasons.append("old_stage13_5_production_force_density_name_reappeared")
        stage14_regression_status = 0
    if "fibre_stage13_6_production_force_density_candidate.dat" not in src_text + stage11_14_text or "stage13_6_production_force_density_candidate_status" not in src_text + stage11_14_text:
        reasons.append("stage13_6_diagnostic_name_evidence_missing")
        stage14_regression_status = 0
    stage13_source = read_text(repo / "src" / "fibre_stage13_production_force_density_candidate.f90")
    center_pattern = re.compile(r"[ijk]0\s*=\s*\(\s*lbound\s*\(\s*ux\s*,\s*[123]\s*\)\s*\+\s*ubound\s*\(\s*ux\s*,\s*[123]\s*\)\s*\)\s*/\s*2", re.I)
    if center_pattern.search(stage13_source) or "lbound(ux, 1) + 2" not in stage13_source:
        reasons.append("stage13_force_density_sampling_repair_missing_or_local_center_detected")
        stage14_regression_status = 0
    for path in [
        repo / "src" / "fibre_stage11_production_oneway_hook.f90",
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "src" / "fibre_stage14_production_rhs_injection.f90",
        repo / "src" / "fibre_stage15_production_structure_hook.f90",
    ]:
        text = read_text(path)
        if "rank0_write_allowed" not in text and "MPI_COMM_RANK" not in text and "nrank" not in text:
            reasons.append(f"rank0_safe_diagnostic_writing_regressed_in_{path}")
            stage14_regression_status = 0

    stage15_regression_status = 1
    for stage, stem in STAGE15_DOC_SCRIPT_STEMS.items():
        if not (repo / "stage15_checks" / f"stage15_{stage}_{stem}.md").exists() or not (repo / "stage15_checks" / f"run_stage15_{stage}_{stem}.sh").exists():
            reasons.append(f"missing_stage15_{stage}_{stem}_script_or_doc")
            stage15_regression_status = 0
    if not all(check_rg_fallback(path) for path in all_files(repo / "stage15_checks", ("*.sh",)) + [wrapper]):
        reasons.append("script_has_rg_without_grep_fallback")
        stage15_regression_status = 0
    stage15_11 = read_text(repo / "stage15_checks" / "run_stage15_11_total_smoke_closure.sh")
    if "add_unmet_final_status_reasons" not in stage15_11 or not re.search(r"if \[ ! -s \"\$REASONS_FILE\" \].*add_unmet_final_status_reasons", stage15_11, re.S):
        reasons.append("stage15_11_unknown_failure_protection_missing")
        stage15_regression_status = 0

    stage16_boundary_status = 1
    xcompact = read_text(repo / "src" / "xcompact3d.f90")
    if "stage16" in xcompact.lower():
        reasons.append("stage16_production_hook_detected_in_xcompact3d")
        stage16_boundary_status = 0
    stage16_src = read_text(config_src) + "\n" + read_text(check_src)
    forbidden_call_patterns = {
        "production_structure_advance_detected": r"\b(use|call)\b.*(fibre_stage15_structure_advance|production_structure|controlled_structure_step)",
        "bending_solve_detected": r"\b(use|call)\b.*(bending|implicit_bending)",
        "tension_solve_detected": r"\b(use|call)\b.*tension",
        "wall_contact_detected": r"\b(use|call)\b.*(wall|contact)",
        "multifibre_detected": r"\b(use|call)\b.*multi[-_ ]?fibre",
        "rhs_modification_detected": r"\b(rhs|ux|uy|uz)\s*=|\b(use|call)\b.*stage14.*rhs.*(inject|addition|accumulator)",
        "pressure_projection_detected": r"\b(use|call)\b.*(pressure|projection)",
        "poisson_detected": r"\b(use|call)\b.*poisson",
        "rk3_channel_forcing_detected": r"\b(use|call)\b.*(rk3|channel.*forc|forcing.*channel)",
        "legacy_ibm_forcing_detected": r"\b(use|call)\b.*(legacy.*ibm|production_ibm|ibm_forcing|apply_ibm)",
    }
    for reason, pattern in forbidden_call_patterns.items():
        if re.search(pattern, stage16_src, re.I):
            reasons.append(reason)
            stage16_boundary_status = 0

    data["stage14_regression_status"] = str(stage14_regression_status)
    data["stage15_regression_status"] = str(stage15_regression_status)
    if stage16_boundary_status == 0:
        for key in ["no_production_hook_status", "no_structure_advance_status", "no_rhs_modification_status", "no_bending_solve_status",
                    "no_tension_solve_status", "no_wall_contact_status", "no_multifibre_status", "no_pressure_projection_modification_status",
                    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status"]:
            data[key] = "0"

    for key in REQUIRED_KEYS:
        if key not in data and key != "final_status":
            reasons.append(f"missing_required_output_key_{key}")
            data[key] = "0"
    for key in PASS_ONE_KEYS:
        if data.get(key) != "1":
            reasons.append(f"{key}_not_pass")
    for key in ["structure_advance_enable_status", "two_way_rhs_enable_status"]:
        if data.get(key) != "0":
            reasons.append(f"{key}_must_remain_zero_in_stage16_1")
    if data.get("no_legacy_ibm_forcing_status") not in (None, "1"):
        reasons.append("no_legacy_ibm_forcing_status_not_pass")
    for key in ["feedback_alpha", "lambda_value", "max_structure_update", "max_rhs_increment"]:
        try:
            float(data.get(key, "nan"))
        except ValueError:
            reasons.append(f"{key}_not_numeric")

    final_ok = len(reasons) == 0
    data["final_status"] = "1" if final_ok else "0"

    with dat.open("w") as handle:
        handle.write("# Stage 16.1 configuration/global-switch summary\n")
        for key in REQUIRED_KEYS:
            handle.write(f"{key} {data.get(key, '0')}\n")
        legacy = ["no_legacy_ibm_forcing_status", "config_status", "numeric_parse_status", "numeric_bounds_status"]
        for key in legacy:
            if key in data and key not in REQUIRED_KEYS:
                handle.write(f"{key} {data[key]}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.1 CONFIG VERDICT: PASS")
        print("STAGE 16.1 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.1 CONFIG VERDICT: FAIL")
    print("STAGE 16.1 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
