#!/usr/bin/env python3
"""Stage 16.0 preflight closure-integrity audit.

This helper is intentionally static/diagnostic-only. It verifies that the
Stage 10-15 closure evidence and anti-regression guards remain present before
Stage 16 production FSI work begins, and it fails closed with explicit reasons.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

REQUIRED_STATUS_KEYS = [
    "stage16_0_requested_status",
    "stage14_closed_file_status",
    "stage15_closed_file_status",
    "stage15_closed_content_status",
    "stage14_regression_status",
    "stage15_regression_status",
    "stage13_6_diagnostic_name_status",
    "stage13_sampling_repair_status",
    "rank0_safe_diagnostic_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "stage16_boundary_status",
    "approved_stage12_13_14_chain_status",
    "no_direct_rhs_injection_status",
    "no_legacy_ibm_forcing_status",
    "no_bending_solve_status",
    "no_tension_solve_status",
    "no_wall_contact_status",
    "no_multifibre_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "final_status",
]

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

FORBIDDEN_STAGE16_PATTERNS = {
    "no_direct_rhs_injection_status": re.compile(r"stage16.*(rhs|right.hand.side).*(inject|add|apply)|direct.*rhs.*inject", re.I),
    "no_legacy_ibm_forcing_status": re.compile(r"legacy.*ibm|production_ibm|apply_ibm|ibm_forcing", re.I),
    "no_bending_solve_status": re.compile(r"bending.*(solve|operator|implicit|production)|implicit_bending", re.I),
    "no_tension_solve_status": re.compile(r"tension.*(solve|production)|tension_solver", re.I),
    "no_wall_contact_status": re.compile(r"wall.*contact|contact.*wall|collision", re.I),
    "no_multifibre_status": re.compile(r"multi[-_ ]?fibre|multifibre|multiple.*fib", re.I),
    "no_pressure_projection_modification_status": re.compile(r"pressure|projection", re.I),
    "no_poisson_modification_status": re.compile(r"poisson", re.I),
    "no_rk3_channel_forcing_modification_status": re.compile(r"rk3|channel.*forc|forcing.*channel", re.I),
}


def truthy(value: str) -> bool:
    return value == "1"


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def existing(paths: list[Path]) -> list[Path]:
    return [path for path in paths if path.exists()]


def all_files_under(root: Path, patterns: tuple[str, ...] = ("*",)) -> list[Path]:
    files: list[Path] = []
    if not root.exists():
        return files
    for pattern in patterns:
        files.extend(path for path in root.rglob(pattern) if path.is_file())
    return sorted(set(files))


def matching_text(files: list[Path]) -> str:
    chunks = []
    for path in files:
        chunks.append(f"\n### {path}\n")
        chunks.append(read_text(path))
    return "".join(chunks)


def set_status(status: dict[str, int], key: str, ok: bool, reasons: list[str], reason: str) -> None:
    status[key] = 1 if ok else 0
    if not ok:
        reasons.append(reason)


def check_rg_fallback(script: Path) -> bool:
    text = read_text(script)
    if not re.search(r"\brg\b", text):
        return True
    has_command_guard = "command -v rg" in text or "which rg" in text
    has_grep_fallback = re.search(r"\bgrep\b", text) is not None
    return has_command_guard and has_grep_fallback


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-stage14-closed", default="1")
    parser.add_argument("--require-stage15-closed", default="1")
    parser.add_argument("--accept-closed-stage15-evidence", default="1")
    parser.add_argument("--enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--wrapper-reasons-file", default="")
    args = parser.parse_args()

    repo = Path.cwd()
    out_dir = repo / "stage16_outputs"
    out_dir.mkdir(exist_ok=True)
    summary = out_dir / "fibre_stage16_0_preflight_closure_integrity.dat"
    reasons_file = out_dir / "stage16_0_preflight_closure_integrity_reasons.tmp"

    status = {key: 0 for key in REQUIRED_STATUS_KEYS}
    reasons: list[str] = []

    if args.wrapper_reasons_file:
        wrapper_reasons = [line.strip() for line in read_text(Path(args.wrapper_reasons_file)).splitlines() if line.strip()]
        reasons.extend(wrapper_reasons)

    requested_ok = truthy(args.enable) and truthy(args.diagnostic_only)
    set_status(
        status,
        "stage16_0_requested_status",
        requested_ok,
        reasons,
        "stage16_0_must_be_enabled_and_diagnostic_only",
    )

    stage14_closed = repo / "stage14_checks" / "STAGE14_CLOSED.md"
    stage15_closed = repo / "stage15_checks" / "STAGE15_CLOSED.md"
    stage14_text = read_text(stage14_closed)
    stage15_text = read_text(stage15_closed)

    stage14_closed_ok = (not truthy(args.require_stage14_closed)) or (stage14_closed.exists() and stage14_closed.stat().st_size > 0)
    set_status(status, "stage14_closed_file_status", stage14_closed_ok, reasons, "missing_or_empty_stage14_closed_file")

    stage15_closed_ok = (not truthy(args.require_stage15_closed)) or (stage15_closed.exists() and stage15_closed.stat().st_size > 0)
    set_status(status, "stage15_closed_file_status", stage15_closed_ok, reasons, "missing_or_empty_stage15_closed_file")

    stage15_content_ok = False
    if stage15_closed_ok and stage15_text:
        closed_ok = re.search(r"STAGE15_CLOSED\s*=\s*YES|Stage\s+15\s+(?:is\s+)?closed", stage15_text, re.I) is not None
        inactive_ok = all(re.search(word, stage15_text, re.I) for word in ["bending", "tension", "wall", "contact"]) and re.search(r"multi[- ]?fibre|multifibre", stage15_text, re.I) and re.search(r"inactive|not\s+active|still\s+inactive", stage15_text, re.I)
        chain_ok = re.search(r"Stage\s*15\s*(?:->|→).*Stage\s*12\s*(?:->|→).*Stage\s*13\s*(?:->|→).*Stage\s*14", stage15_text, re.I | re.S) is not None
        stage15_content_ok = bool(closed_ok and inactive_ok and chain_ok)
        if not closed_ok:
            reasons.append("stage15_closed_file_does_not_state_stage15_closed")
        if not inactive_ok:
            reasons.append("stage15_closed_file_does_not_preserve_inactive_full_physics_statement")
        if not chain_ok:
            reasons.append("stage15_closed_file_does_not_preserve_stage15_to_stage12_to_stage13_to_stage14_route")
    else:
        reasons.append("stage15_closed_content_unavailable")
    status["stage15_closed_content_status"] = 1 if stage15_content_ok else 0

    src_files = all_files_under(repo / "src", ("*.f90",))
    stage11_files = all_files_under(repo / "stage11_checks", ("*.sh", "*.md"))
    stage13_files = all_files_under(repo / "stage13_checks", ("*.sh", "*.md"))
    stage14_files = all_files_under(repo / "stage14_checks", ("*.sh", "*.md"))
    stage15_files = all_files_under(repo / "stage15_checks", ("*.sh", "*.md"))
    src_stage_text = matching_text(src_files)
    stage11_14_text = matching_text(stage11_files + stage13_files + stage14_files)
    stage15_checks_text = matching_text(stage15_files)

    forbidden_lambda_gate_absent = re.search(r"stage14_get_injection_gain\s*\(\s*\)\s*(?:==|\.eq\.)\s*0\.0", src_stage_text + stage11_14_text, re.I) is None
    stage11_diag = (repo / "stage11_checks" / "run_stage11_5_production_oneway_hook.sh").exists() and "stage11_5_production_oneway_hook" in src_stage_text + stage11_14_text
    stage13_diag = (repo / "stage13_checks" / "run_stage13_6_production_force_density_candidate.sh").exists() and "stage13_6_production_force_density_candidate" in src_stage_text + stage11_14_text
    stage14_diag = (repo / "stage14_checks" / "run_stage14_5_production_rhs_hook.sh").exists() and "stage14_5_production_rhs_hook" in src_stage_text + stage11_14_text
    small_lambda = (repo / "stage14_checks" / "run_stage14_7_small_lambda_response_np1.sh").exists() and (repo / "stage14_checks" / "run_stage14_8_parallel_small_lambda_response.sh").exists() and re.search(r"small[-_ ]lambda|lambda.*nonzero|lambda_zero_status\s*==\s*0", src_stage_text + stage11_14_text, re.I)
    stage14_regression_ok = all([forbidden_lambda_gate_absent, stage11_diag, stage13_diag, stage14_diag, bool(small_lambda)])
    status["stage14_regression_status"] = 1 if stage14_regression_ok else 0
    if not forbidden_lambda_gate_absent:
        reasons.append("stage14_forbidden_lambda_zero_registration_gate_found")
    if not small_lambda:
        reasons.append("stage14_small_lambda_rhs_hook_diagnostic_or_registration_evidence_missing")
    if not stage11_diag:
        reasons.append("stage11_5_production_oneway_hook_diagnostics_missing")
    if not stage13_diag:
        reasons.append("stage13_6_production_force_density_diagnostics_missing")
    if not stage14_diag:
        reasons.append("stage14_5_production_rhs_hook_diagnostics_missing")

    stage13_6_name_ok = "stage13_6_production_force_density_candidate" in src_stage_text + stage11_14_text and "fibre_stage13_6_production_force_density_candidate.dat" in src_stage_text + stage11_14_text
    old_stage13_5_prod_absent = "stage13_5_production_force_density_candidate" not in src_stage_text + stage11_14_text + stage15_checks_text
    status["stage13_6_diagnostic_name_status"] = 1 if stage13_6_name_ok and old_stage13_5_prod_absent else 0
    if not stage13_6_name_ok:
        reasons.append("current_stage13_6_production_force_density_diagnostic_names_missing")
    if not old_stage13_5_prod_absent:
        reasons.append("old_stage13_5_production_force_density_name_reappeared")

    stage13_source = read_text(repo / "src" / "fibre_stage13_production_force_density_candidate.f90")
    center_pattern = re.compile(r"[ijk]0\s*=\s*\(\s*lbound\s*\(\s*ux\s*,\s*[123]\s*\)\s*\+\s*ubound\s*\(\s*ux\s*,\s*[123]\s*\)\s*\)\s*/\s*2", re.I)
    sampling_repair_ok = center_pattern.search(stage13_source) is None and "lbound(ux, 1) + 2" in stage13_source and "local subdomain centre" in stage13_source
    set_status(status, "stage13_sampling_repair_status", sampling_repair_ok, reasons, "stage13_force_density_local_subdomain_center_sampling_detected_or_repair_evidence_missing")

    rank0_files = [
        repo / "src" / "fibre_stage11_production_oneway_hook.f90",
        repo / "src" / "fibre_stage13_production_force_density_candidate.f90",
        repo / "src" / "fibre_stage14_production_rhs_injection.f90",
        repo / "src" / "fibre_stage15_production_structure_hook.f90",
    ]
    rank0_ok = all("rank0_write_allowed" in read_text(path) or "MPI_COMM_RANK" in read_text(path) or "nrank" in read_text(path) for path in rank0_files)
    set_status(status, "rank0_safe_diagnostic_status", rank0_ok, reasons, "rank0_safe_diagnostic_writing_regressed_for_stage11_13_14_or_15")

    missing_stage15_docs = []
    for stage, stem in STAGE15_DOC_SCRIPT_STEMS.items():
        md = repo / "stage15_checks" / f"stage15_{stage}_{stem}.md"
        sh = repo / "stage15_checks" / f"run_stage15_{stage}_{stem}.sh"
        if not md.exists() or not sh.exists():
            missing_stage15_docs.append(f"stage15_{stage}_{stem}")
    if missing_stage15_docs:
        reasons.append("missing_stage15_docs_or_scripts:" + ",".join(missing_stage15_docs))

    evidence_ok = True
    if truthy(args.accept_closed_stage15_evidence) and stage15_content_ok:
        for stage in range(1, 12):
            if not re.search(fr"Stage\s*15\.{stage}\b|stage15_{stage}_", stage15_text, re.I):
                evidence_ok = False
                reasons.append(f"stage15_{stage}_closed_evidence_not_named_in_stage15_closed_file")
    else:
        for stage in range(1, 12):
            output_hits = list((repo / "stage15_outputs").glob(f"*stage15_{stage}*.dat")) if (repo / "stage15_outputs").exists() else []
            if not output_hits:
                evidence_ok = False
                reasons.append(f"stage15_{stage}_diagnostic_evidence_missing_and_closed_evidence_not_accepted")

    no_rg_only = all(check_rg_fallback(path) for path in all_files_under(repo / "stage15_checks", ("*.sh",)))
    set_status(status, "no_rg_only_dependency_status", no_rg_only, reasons, "stage15_or_stage16_script_has_rg_without_grep_fallback")

    stage15_11_text = read_text(repo / "stage15_checks" / "run_stage15_11_total_smoke_closure.sh")
    no_unknown_ok = "add_unmet_final_status_reasons" in stage15_11_text and "unknown_stage15_11_failure" in stage15_11_text and re.search(r"if \[ ! -s \"\$REASONS_FILE\" \].*add_unmet_final_status_reasons", stage15_11_text, re.S)
    set_status(status, "no_unknown_failure_status", bool(no_unknown_ok), reasons, "stage15_11_can_fail_without_specific_reasons")

    closure_generation_ok = "rm -f \"$CLOSURE_FILE\"" in stage15_11_text and "write_closure_file" in stage15_11_text and "STAGE15_CLOSED=YES" in stage15_11_text and "STAGE15_CLOSED=NO" in stage15_11_text
    stage15_regression_ok = not missing_stage15_docs and evidence_ok and no_rg_only and bool(no_unknown_ok) and closure_generation_ok
    status["stage15_regression_status"] = 1 if stage15_regression_ok else 0
    if not closure_generation_ok:
        reasons.append("stage15_11_closure_file_generation_not_guarded_by_pass_fail_evidence")

    stage16_paths = all_files_under(repo / "stage16_checks", ("*.sh", "*.py", "*.md")) + all_files_under(repo / "src", ("*stage16*.f90",))
    stage16_text_by_status = {key: False for key in FORBIDDEN_STAGE16_PATTERNS}
    for path in stage16_paths:
        text = read_text(path)
        if path.name in {"stage16_0_preflight_closure_integrity.md", "assert_stage16_0_preflight_closure_integrity.py"}:
            # Documentation and this parser necessarily name prohibited physics.
            # The executable wrapper and any Stage 16 source must not implement it.
            continue
        for key, pattern in FORBIDDEN_STAGE16_PATTERNS.items():
            if pattern.search(text):
                stage16_text_by_status[key] = True
                reasons.append(f"stage16_boundary_forbidden_pattern_in_{path}:{key}")
    for key, found in stage16_text_by_status.items():
        status[key] = 0 if found else 1

    stage16_src_files = all_files_under(repo / "src", ("*stage16*.f90",))
    xcompact_text = read_text(repo / "src" / "xcompact3d.f90")
    stage16_hook_absent = "stage16" not in xcompact_text.lower()
    closed_stage_write_pattern = re.compile(r"(rm\s+-f|mv\s+|cp\s+|sed\s+-i|tee\s+|>+)\s+[^\n;]*(stage1[0-5]_checks|stage1[0-5]_outputs|STAGE1[0-5]_CLOSED)", re.I)
    stage16_modifies_closed = any(
        closed_stage_write_pattern.search(read_text(path))
        for path in stage16_paths
        if path.name != "assert_stage16_0_preflight_closure_integrity.py"
    )
    required_stage16_files = (repo / "stage16_checks" / "run_stage16_0_preflight_closure_integrity.sh").exists() and (repo / "stage16_checks" / "stage16_0_preflight_closure_integrity.md").exists()
    stage16_boundary_ok = not stage16_src_files and stage16_hook_absent and not stage16_modifies_closed and required_stage16_files and all(status[key] == 1 for key in FORBIDDEN_STAGE16_PATTERNS)
    status["stage16_boundary_status"] = 1 if stage16_boundary_ok else 0
    if stage16_src_files:
        reasons.append("stage16_source_files_exist_before_stage16_0_physics_is_allowed")
    if not stage16_hook_absent:
        reasons.append("stage16_hook_detected_in_xcompact3d_f90")
    if stage16_modifies_closed:
        reasons.append("stage16_script_contains_closed_stage10_15_write_pattern")
    if not required_stage16_files:
        reasons.append("missing_stage16_0_required_doc_or_wrapper")

    chain_terms = [
        r"controlled\s+structure\s+update",
        r"Stage\s*12.*feedback",
        r"Stage\s*13.*force[- ]density",
        r"Stage\s*14.*RHS",
        r"Stage\s*12\s*(?:->|→).*Stage\s*13\s*(?:->|→).*Stage\s*14",
        r"no\s+direct.*RHS|direct.*RHS.*inactive|outside\s+approved\s+chain",
        r"legacy\s+IBM.*inactive|no\s+legacy\s+IBM|outside\s+approved\s+chain",
    ]
    chain_ok = bool(stage15_text) and all(re.search(pattern, stage15_text, re.I | re.S) for pattern in chain_terms)
    status["approved_stage12_13_14_chain_status"] = 1 if chain_ok else 0
    if not chain_ok:
        reasons.append("stage15_closed_file_does_not_preserve_complete_approved_stage12_13_14_chain_evidence")

    # Make status aliases fail closed when the boundary or chain audit failed.
    if not chain_ok:
        status["no_direct_rhs_injection_status"] = 0
        status["no_legacy_ibm_forcing_status"] = 0
    if not stage16_boundary_ok:
        for key in FORBIDDEN_STAGE16_PATTERNS:
            status[key] = min(status[key], 0)

    final_ok = all(status[key] == 1 for key in REQUIRED_STATUS_KEYS if key != "final_status")
    status["final_status"] = 1 if final_ok else 0

    with summary.open("w") as handle:
        handle.write("# Stage 16.0 preflight closure integrity summary\n")
        for key in REQUIRED_STATUS_KEYS:
            handle.write(f"{key} {status[key]}\n")
        if reasons:
            handle.write("failure_reasons_begin\n")
            for reason in reasons:
                handle.write(f"failure_reason {reason}\n")
            handle.write("failure_reasons_end\n")
    reasons_file.write_text("\n".join(reasons) + ("\n" if reasons else ""))

    if final_ok:
        print("STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS")
        print("STAGE 16.0 FINAL VERDICT: PASS")
        return 0

    print("STAGE 16.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: FAIL")
    print("STAGE 16.0 FINAL VERDICT: FAIL")
    print("Reasons:")
    for reason in reasons:
        print(f" - {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
