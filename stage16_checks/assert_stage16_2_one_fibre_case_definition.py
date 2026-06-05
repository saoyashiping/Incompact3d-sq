#!/usr/bin/env python3
"""Stage 16.2 one-fibre case-definition audit helper."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REQUIRED_KEYS = [
    "stage16_2_requested_status",
    "one_fibre_case_definition_status",
    "one_fibre_count_status",
    "npts",
    "npts_valid_status",
    "position_finite_status",
    "velocity_finite_status",
    "acceleration_finite_status",
    "initial_velocity_bound_status",
    "initial_acceleration_bound_status",
    "min_point_spacing",
    "point_spacing_status",
    "min_wall_clearance",
    "wall_clearance_status",
    "domain_containment_status",
    "invalid_npts_rejection_status",
    "invalid_geometry_rejection_status",
    "invalid_velocity_rejection_status",
    "invalid_acceleration_rejection_status",
    "wall_contact_rejection_status",
    "multifibre_rejection_status",
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
    "stage16_1_regression_status",
    "final_status",
]

NUMERIC_KEYS = {"npts", "min_point_spacing", "min_wall_clearance"}
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
    out: list[Path] = []
    for pattern in patterns:
        out.extend(path for path in root.rglob(pattern) if path.is_file())
    return sorted(set(out))


def joined(files: list[Path]) -> str:
    return "\n".join(f"### {path}\n{read_text(path)}" for path in files)


def parse_dat(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def active_shell_has_rg_command(script: Path) -> bool:
    """Return True only for executable rg command use, not quoted regex/doc text."""
    for raw in read_text(script).splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # Strip simple trailing comments; command detection should not treat quoted regex text as rg use.
        unquoted_prefix = line.split("#", 1)[0].strip()
        if re.match(r"^(?:command\s+-v\s+)?rg(?:\s|$)", unquoted_prefix):
            return True
        if re.search(r"(?:^|[;&|()])\s*rg\s+", unquoted_prefix):
            return True
    return False


def script_has_rg_only_dependency(script: Path) -> bool:
    if not active_shell_has_rg_command(script):
        return False
    text = read_text(script)
    has_fallback = ("command -v rg" in text or "which rg" in text) and re.search(r"\bgrep\b", text) is not None
    return not has_fallback


def stage16_1_closed_static_evidence_ok(repo: Path) -> bool:
    required = [
        repo / "stage16_checks" / "run_stage16_1_config.sh",
        repo / "stage16_checks" / "stage16_1_config.md",
        repo / "stage16_checks" / "assert_stage16_1_config.py",
        repo / "src" / "fibre_stage16_config.f90",
        repo / "src" / "fibre_stage16_config_check.f90",
        repo / "stage14_checks" / "STAGE14_CLOSED.md",
        repo / "stage15_checks" / "STAGE15_CLOSED.md",
    ]
    cmake_ok = "add_executable(fibre_stage16_config_check" in read_text(repo / "src" / "CMakeLists.txt")
    return cmake_ok and all(path.exists() and path.stat().st_size > 0 for path in required)


def stage14_lambda_zero_gate_found(repo: Path) -> bool:
    # Only executable/source evidence is scanned.  Markdown prohibition text must not be
    # misclassified as a real Stage 14 lambda==0 hook gate regression.
    files = all_files(repo / "src", ("*.f90",)) + all_files(repo / "stage14_checks", ("*.sh",))
    pattern = re.compile(r"stage14_get_injection_gain\s*\(\s*\)\s*(?:==|\.eq\.)\s*0\.0", re.I)
    return any(pattern.search(read_text(path)) for path in files)


def old_stage13_5_production_name_reappeared(repo: Path) -> bool:
    # Stage 13.5 conservation-sign audit is a legitimate closed substep.  Only the old
    # production force-density candidate name is forbidden, and negative audit text in
    # Stage 15/16 scripts/docs must not be treated as a regression.
    files = [
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "src" / "fibre_stage13_production_force_density_candidate_check.f90",
        repo / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh",
        repo / "stage13_checks" / "stage13_6_production_force_density_candidate.md",
    ]
    return any("stage13_5_production_force_density_candidate" in read_text(path) for path in files if path.exists())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-status", default="0")
    parser.add_argument("--run-status", default="0")
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--require-stage16-1", default="1")
    parser.add_argument("--accept-stage16-1-closed-evidence", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    dat = out_dir / "fibre_stage16_2_one_fibre_case_definition.dat"
    reasons_file = out_dir / "stage16_2_one_fibre_case_definition_reasons.tmp"
    data = parse_dat(dat)
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        reasons.extend(line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip())
    if args.build_status != "1":
        reasons.append("stage16_2_case_check_target_build_status_not_pass")
    if args.run_status != "1":
        reasons.append("stage16_2_case_check_target_run_status_not_pass")

    wrapper = repo / "stage16_checks" / "run_stage16_2_one_fibre_case_definition.sh"
    doc = repo / "stage16_checks" / "stage16_2_one_fibre_case_definition.md"
    helper = repo / "stage16_checks" / "assert_stage16_2_one_fibre_case_definition.py"
    case_src = repo / "src" / "fibre_stage16_one_fibre_case.f90"
    check_src = repo / "src" / "fibre_stage16_one_fibre_case_check.f90"
    cmake = repo / "src" / "CMakeLists.txt"
    for path in [wrapper, doc, helper, case_src, check_src]:
        if not path.exists():
            reasons.append(f"missing_required_stage16_2_file_{path}")
    if "add_executable(fibre_stage16_one_fibre_case_check" not in read_text(cmake):
        reasons.append("missing_fibre_stage16_one_fibre_case_check_build_registration")

    if args.require_stage14_closed == "1" and not (repo / "stage14_checks" / "STAGE14_CLOSED.md").exists():
        reasons.append("missing_stage14_closed_file")
    if args.require_stage15_closed == "1" and not (repo / "stage15_checks" / "STAGE15_CLOSED.md").exists():
        reasons.append("missing_stage15_closed_file")
    if args.require_stage16_1 == "1":
        stage16_1 = parse_dat(repo / "stage16_outputs" / "fibre_stage16_1_config.dat")
        if stage16_1.get("final_status") != "1":
            if args.accept_stage16_1_closed_evidence == "1" and stage16_1_closed_static_evidence_ok(repo):
                pass
            else:
                reasons.append("missing_or_failed_stage16_1_config_evidence")

    src_text = joined(all_files(repo / "src", ("*.f90",)))
    stage11_14_text = joined(all_files(repo / "stage11_checks", ("*.sh", "*.md")) +
                            all_files(repo / "stage13_checks", ("*.sh", "*.md")) +
                            all_files(repo / "stage14_checks", ("*.sh", "*.md")))
    stage15_text = joined(all_files(repo / "stage15_checks", ("*.sh", "*.md")))

    stage14_regression_status = 1
    if stage14_lambda_zero_gate_found(repo):
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
    if old_stage13_5_production_name_reappeared(repo):
        reasons.append("old_stage13_5_production_force_density_name_reappeared")
        stage14_regression_status = 0
    if "fibre_stage13_6_production_force_density_candidate.dat" not in src_text + stage11_14_text:
        reasons.append("stage13_6_force_density_output_name_missing")
        stage14_regression_status = 0
    if "stage13_6_production_force_density_candidate_status" not in src_text + stage11_14_text:
        reasons.append("stage13_6_force_density_status_name_missing")
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
        check_src,
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
    scripts_to_check = all_files(repo / "stage15_checks", ("*.sh",)) + all_files(repo / "stage16_checks", ("*.sh",))
    if any(script_has_rg_only_dependency(path) for path in scripts_to_check):
        reasons.append("script_has_rg_without_grep_fallback")
        stage15_regression_status = 0
    stage15_11 = read_text(repo / "stage15_checks" / "run_stage15_11_total_smoke_closure.sh")
    if "add_unmet_final_status_reasons" not in stage15_11 or not re.search(r"if \[ ! -s \"\$REASONS_FILE\" \].*add_unmet_final_status_reasons", stage15_11, re.S):
        reasons.append("stage15_11_unknown_failure_protection_missing")
        stage15_regression_status = 0

    stage16_1_regression_status = 1
    for path in [
        repo / "stage16_checks" / "run_stage16_0_preflight_closure_integrity.sh",
        repo / "stage16_checks" / "assert_stage16_0_preflight_closure_integrity.py",
        repo / "stage16_checks" / "stage16_0_preflight_closure_integrity.md",
        repo / "stage16_checks" / "run_stage16_1_config.sh",
        repo / "stage16_checks" / "assert_stage16_1_config.py",
        repo / "stage16_checks" / "stage16_1_config.md",
        repo / "src" / "fibre_stage16_config.f90",
        repo / "src" / "fibre_stage16_config_check.f90",
    ]:
        if not path.exists():
            reasons.append(f"missing_stage16_0_or_stage16_1_evidence_file_{path}")
            stage16_1_regression_status = 0
    if "add_executable(fibre_stage16_config_check" not in read_text(cmake):
        reasons.append("stage16_1_config_check_build_registration_missing")
        stage16_1_regression_status = 0

    stage16_boundary_status = 1
    if "stage16" in read_text(repo / "src" / "xcompact3d.f90").lower():
        reasons.append("stage16_production_hook_detected_in_xcompact3d")
        stage16_boundary_status = 0
    stage16_2_src = read_text(case_src) + "\n" + read_text(check_src)
    forbidden_call_patterns = {
        "production_structure_advance_detected": r"\b(use|call)\b.*(fibre_stage15_structure_advance|production_structure|controlled_structure_step)",
        "bending_solve_detected": r"\b(use|call)\b.*(bending|implicit_bending)",
        "tension_solve_detected": r"\b(use|call)\b.*tension",
        "rhs_modification_detected": r"\b(rhs|ux|uy|uz)\s*=|\b(use|call)\b.*stage14.*rhs.*(inject|addition|accumulator)",
        "pressure_projection_detected": r"\b(use|call)\b.*(pressure|projection)",
        "poisson_detected": r"\b(use|call)\b.*poisson",
        "rk3_channel_forcing_detected": r"\b(use|call)\b.*(rk3|channel.*forc|forcing.*channel)",
        "legacy_ibm_forcing_detected": r"\b(use|call)\b.*(legacy.*ibm|production_ibm|ibm_forcing|apply_ibm)",
    }
    for reason, pattern in forbidden_call_patterns.items():
        if re.search(pattern, stage16_2_src, re.I):
            reasons.append(reason)
            stage16_boundary_status = 0

    data["stage14_regression_status"] = str(stage14_regression_status)
    data["stage15_regression_status"] = str(stage15_regression_status)
    data["stage16_1_regression_status"] = str(stage16_1_regression_status)
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
    for key in NUMERIC_KEYS:
        try:
            float(data.get(key, "nan"))
        except ValueError:
            reasons.append(f"{key}_not_numeric")
    if data.get("no_legacy_ibm_forcing_status") not in (None, "1"):
        reasons.append("no_legacy_ibm_forcing_status_not_pass")

    final_ok = len(reasons) == 0
    data["final_status"] = "1" if final_ok else "0"

    with dat.open("w") as handle:
        handle.write("# Stage 16.2 one-fibre case-definition summary\n")
        for key in REQUIRED_KEYS:
            handle.write(f"{key} {data.get(key, '0')}\n")
        for key in ["no_legacy_ibm_forcing_status", "case_status", "stage16_1_config_status"]:
            if key in data and key not in REQUIRED_KEYS:
                handle.write(f"{key} {data[key]}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.2 ONE-FIBRE CASE DEFINITION VERDICT: PASS")
        print("STAGE 16.2 FINAL VERDICT: PASS")
        return 0
    print("STAGE 16.2 ONE-FIBRE CASE DEFINITION VERDICT: FAIL")
    print("STAGE 16.2 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
