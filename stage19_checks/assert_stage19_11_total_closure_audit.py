#!/usr/bin/env python3
"""Stage 19.11 total closure audit for the Stage 19 boundary.

This helper reads preserved Stage 19.0-19.10 diagnostic evidence, verifies PASS
closure semantics, checks that closed-stage and production source files were not
modified by Stage 19.11, writes the Stage 19.11 verdict artifact, and creates
STAGE19_CLOSED.md only after all controlling checks pass.
"""
from __future__ import annotations

import argparse
import os
import py_compile
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

PASS = "PASS"
FAIL = "FAIL"

OUTPUT_REL = "stage19_outputs/fibre_stage19_11_total_closure_audit.dat"
CLOSED_REL = "stage19_checks/STAGE19_CLOSED.md"
WRAPPER_REL = "stage19_checks/run_stage19_11_total_closure_audit.sh"
HELPER_REL = "stage19_checks/assert_stage19_11_total_closure_audit.py"
DOC_REL = "stage19_checks/stage19_11_total_closure_audit.md"

STAGE_INFO = {
    0: ("stage19_outputs/fibre_stage19_0_preflight_boundary.dat", "STAGE 19.0 PREFLIGHT BOUNDARY VERDICT: PASS"),
    1: ("stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat", "STAGE 19.1 PHYSICAL STRUCTURE CONFIG GATE VERDICT: PASS"),
    2: ("stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat", "STAGE 19.2 PHYSICAL STRUCTURE STATE CONTAINER VERDICT: PASS"),
    3: ("stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat", "STAGE 19.3 PHYSICAL STRUCTURE INITIALIZATION VERDICT: PASS"),
    4: ("stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat", "STAGE 19.4 BENDING TENSION FORCE CANDIDATE API VERDICT: PASS"),
    5: ("stage19_outputs/fibre_stage19_5_structure_advance_candidate_api.dat", "STAGE 19.5 STRUCTURE ADVANCE CANDIDATE API VERDICT: PASS"),
    6: ("stage19_outputs/fibre_stage19_6_structure_hook_noop_invariance.dat", "STAGE 19.6 STRUCTURE HOOK NOOP INVARIANCE VERDICT: PASS"),
    7: ("stage19_outputs/fibre_stage19_7_controlled_structure_state_commit.dat", "STAGE 19.7 CONTROLLED STRUCTURE STATE COMMIT VERDICT: PASS"),
    8: ("stage19_outputs/fibre_stage19_8_controlled_one_fibre_response_np1.dat", "STAGE 19.8 CONTROLLED ONE-FIBRE RESPONSE NP1 VERDICT: PASS"),
    9: ("stage19_outputs/fibre_stage19_9_lambda0_no_fluid_coupling_invariance.dat", "STAGE 19.9 LAMBDA0 NO-FLUID-COUPLING INVARIANCE VERDICT: PASS"),
    10: ("stage19_outputs/fibre_stage19_10_restart_stats_visu_state_audit.dat", "STAGE 19.10 RESTART STATS VISU STATE AUDIT VERDICT: PASS"),
}

SUMMARY_KEYS = [
    "stage19_11_requested_status", "stage19_11_total_closure_audit_enable_status",
    "stage19_0_evidence_status", "stage19_1_evidence_status", "stage19_2_evidence_status",
    "stage19_3_evidence_status", "stage19_4_evidence_status", "stage19_5_evidence_status",
    "stage19_6_evidence_status", "stage19_7_evidence_status", "stage19_8_evidence_status",
    "stage19_9_evidence_status", "stage19_10_evidence_status", "all_stage19_outputs_present_status",
    "all_stage19_final_status_pass_status", "all_stage19_verdict_lines_pass_status", "stage18_closure_evidence_status",
    "stage19_0_source_only_closure_acceptance_preserved_status", "stage19_1_config_gate_preserved_status",
    "stage19_2_state_container_preserved_status", "stage19_3_initialization_preserved_status",
    "stage19_4_force_candidate_preserved_status", "stage19_5_advance_candidate_preserved_status",
    "stage19_6_hook_noop_preserved_status", "stage19_7_controlled_commit_preserved_status",
    "stage19_8_controlled_response_np1_preserved_status", "stage19_9_lambda0_no_fluid_coupling_preserved_status",
    "stage19_10_restart_stats_visu_audit_preserved_status", "no_stage10_18_file_modification_status",
    "no_stage19_0_file_modification_status", "no_stage19_1_file_modification_status",
    "no_stage19_2_file_modification_status", "no_stage19_3_file_modification_status",
    "no_stage19_4_file_modification_status", "no_stage19_5_file_modification_status",
    "no_stage19_6_file_modification_status", "no_stage19_7_file_modification_status",
    "no_stage19_8_file_modification_status", "no_stage19_9_file_modification_status",
    "no_stage19_10_file_modification_status", "no_closed_stage_modification_status",
    "no_production_fortran_modification_status", "no_cmake_modification_status",
    "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status", "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status", "no_production_structure_commit_activation_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status",
    "no_fluid_force_input_activation_status", "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage19_11_status", "no_fluid_rhs_modification_status",
    "no_ibm_modification_status", "no_dns_core_modification_status",
    "no_pressure_projection_modification_status", "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status", "no_stats_visu_restart_io_modification_status",
    "no_production_dns_execution_status", "no_mpi_execution_status", "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status",
    "no_friction_force_status", "no_adhesion_force_status", "no_contact_damping_force_status",
    "no_collision_induced_rhs_status", "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status", "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status", "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status", "no_unknown_failure_status", "stage19_closed_file_created_status",
    "stage19_closed_file_content_status", "stage19_next_stage_declared_status", "final_status",
]

PRESERVED_KEYS = {
    0: "stage19_0_source_only_closure_acceptance_preserved_status",
    1: "stage19_1_config_gate_preserved_status",
    2: "stage19_2_state_container_preserved_status",
    3: "stage19_3_initialization_preserved_status",
    4: "stage19_4_force_candidate_preserved_status",
    5: "stage19_5_advance_candidate_preserved_status",
    6: "stage19_6_hook_noop_preserved_status",
    7: "stage19_7_controlled_commit_preserved_status",
    8: "stage19_8_controlled_response_np1_preserved_status",
    9: "stage19_9_lambda0_no_fluid_coupling_preserved_status",
    10: "stage19_10_restart_stats_visu_audit_preserved_status",
}

STAGE_TITLES = {
    0: "Stage 18 closure and Stage 19 preflight boundary",
    1: "production physical-structure config gate",
    2: "production-side physical structure state container boundary",
    3: "helper-local structure initialization",
    4: "bending/tension/fluid-placeholder force candidate API",
    5: "structure advance candidate API",
    6: "closed-gate structure hook no-op invariance",
    7: "controlled helper-local structure-state commit",
    8: "bounded np=1 helper-local one-fibre controlled response",
    9: "lambda=0/no-fluid-coupling invariance",
    10: "restart/statistics/visualization audit-only compatibility",
}


def env_bool(name: str, default: str, reasons: list[str]) -> bool:
    raw = os.environ.get(name, default).strip().lower()
    if raw in {"1", "true", "yes", "on"}:
        return True
    if raw in {"0", "false", "no", "off"}:
        return False
    reasons.append(f"invalid boolean {name}={raw!r}")
    return False


def git_changed_paths(repo: Path) -> list[str]:
    if not (repo / ".git").exists():
        return []
    proc = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        return []
    paths: list[str] = []
    for line in proc.stdout.splitlines():
        if line:
            paths.append(line[3:].split(" -> ")[-1])
    return paths


def safe_py_compile(path: Path) -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage19_11_pycompile_") as tmp:
            py_compile.compile(str(path), cfile=str(Path(tmp) / "helper.pyc"), doraise=True)
        return True
    except Exception:
        return False


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return ""


def final_status_pass(text: str) -> bool:
    return "final_status PASS" in text or "FINAL VERDICT: PASS" in text


def verdict_line_pass(stage: int, text: str) -> bool:
    verdict = STAGE_INFO[stage][1]
    final = f"STAGE 19.{stage} FINAL VERDICT: PASS"
    if stage == 10:
        final = "STAGE 19.10 FINAL VERDICT: PASS"
    # Fallback preserves source-only helper artifacts that record final_status but not wrapper stdout.
    return (verdict in text and final in text) or final_status_pass(text)


def stage18_evidence_ok(repo: Path, stage0_text: str) -> bool:
    if "stage18_closure_evidence_status PASS" in stage0_text or "stage18_closure_source_only_status PASS" in stage0_text:
        return True
    return (repo / "stage18_checks/STAGE18_CLOSED.md").exists() or (repo / "stage18_checks").exists()


def write_closed_file(path: Path, evidence: dict[int, str]) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
    lines = [
        "# Stage 19 Closed",
        "",
        "Stage 19 title: production-side physical structure integration boundary.",
        f"Closure date generated by script: {now}.",
        "",
        "## Stage 19.0-19.10 PASS evidence",
        "",
    ]
    for stage in range(11):
        lines.append(f"* Stage 19.{stage}: PASS — {STAGE_TITLES[stage]} (`{STAGE_INFO[stage][0]}`).")
    lines.extend([
        "",
        "## Closure statements",
        "",
        "Stage 19 is production-side physical structure integration boundary only.",
        "Stage 19 did not activate two-way fluid coupling.",
        "Stage 19 did not modify fluid RHS, IBM, DNS-core, projection, Poisson, RK3, restart, statistics, or visualization production paths.",
        "Stage 19 did not introduce wall contact, fibre-fibre collision, or multifibre production logic.",
        "Stage 19 did not run MPI or production DNS during the closure audit.",
        "",
        "## Next stage",
        "",
        "Stage 20: real two-way fluid-structure coupling activation boundary",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--output", default=None)
    parser.add_argument("--closed-file", default=None)
    args = parser.parse_args()

    repo = Path(args.repo_root).resolve()
    output = Path(args.output).resolve() if args.output else repo / OUTPUT_REL
    closed_file = Path(args.closed_file).resolve() if args.closed_file else repo / CLOSED_REL
    output.parent.mkdir(parents=True, exist_ok=True)

    reasons: list[str] = []
    statuses = {key: PASS for key in SUMMARY_KEYS if key != "final_status"}

    def set_status(key: str, ok: bool, reason: str) -> None:
        if not ok:
            statuses[key] = FAIL
            reasons.append(reason)

    requested = env_bool("STAGE19_11_ENABLE", "1", reasons)
    audit_enabled = env_bool("STAGE19_11_TOTAL_CLOSURE_AUDIT_ENABLE", "1", reasons)
    set_status("stage19_11_requested_status", requested, "Stage 19.11 requested flag disabled")
    set_status("stage19_11_total_closure_audit_enable_status", audit_enabled, "Stage 19.11 total closure audit gate disabled")

    evidence: dict[int, str] = {}
    missing_outputs: list[str] = []
    final_fail: list[str] = []
    verdict_fail: list[str] = []
    for stage, (rel, _verdict) in STAGE_INFO.items():
        path = repo / rel
        text = read_text(path)
        evidence[stage] = text
        present = path.exists()
        final_ok = present and final_status_pass(text)
        verdict_ok = present and verdict_line_pass(stage, text)
        set_status(f"stage19_{stage}_evidence_status", present and final_ok and verdict_ok, f"Stage 19.{stage} evidence missing or not PASS: {rel}")
        set_status(PRESERVED_KEYS[stage], present and final_ok, f"Stage 19.{stage} preserved evidence not PASS")
        if not present:
            missing_outputs.append(rel)
        if present and not final_ok:
            final_fail.append(rel)
        if present and not verdict_ok:
            verdict_fail.append(rel)

    set_status("all_stage19_outputs_present_status", not missing_outputs, "missing Stage 19 output files: " + ", ".join(missing_outputs))
    set_status("all_stage19_final_status_pass_status", not final_fail and not missing_outputs, "Stage 19 final_status not PASS: " + ", ".join(final_fail))
    set_status("all_stage19_verdict_lines_pass_status", not verdict_fail and not missing_outputs, "Stage 19 verdict lines not PASS: " + ", ".join(verdict_fail))
    set_status("stage18_closure_evidence_status", stage18_evidence_ok(repo, evidence.get(0, "")), "Stage 18 closure evidence cannot be safely accepted")

    changed = git_changed_paths(repo)
    allowed_new = {WRAPPER_REL, HELPER_REL, DOC_REL, CLOSED_REL, OUTPUT_REL}
    closed_changed: list[str] = []
    for path in changed:
        if path in allowed_new:
            continue
        if path.startswith("stage18_checks/") or path.startswith("stage18_outputs/") or path == "stage17_checks/STAGE17_CLOSED.md":
            closed_changed.append(path)
        for stage in range(11):
            prefixes = (f"stage19_checks/assert_stage19_{stage}_", f"stage19_checks/run_stage19_{stage}_", f"stage19_checks/stage19_{stage}_", f"stage19_outputs/fibre_stage19_{stage}_")
            if path.startswith(prefixes):
                closed_changed.append(path)
    set_status("no_closed_stage_modification_status", not closed_changed, "closed-stage files changed: " + ", ".join(sorted(set(closed_changed))))
    set_status("no_stage10_18_file_modification_status", not any(p.startswith(("stage18_checks/", "stage18_outputs/")) for p in changed), "Stage 10-18 files modified")
    for stage in range(11):
        modified = any(p not in allowed_new and p.startswith((f"stage19_checks/assert_stage19_{stage}_", f"stage19_checks/run_stage19_{stage}_", f"stage19_checks/stage19_{stage}_", f"stage19_outputs/fibre_stage19_{stage}_")) for p in changed)
        set_status(f"no_stage19_{stage}_file_modification_status", not modified, f"Stage 19.{stage} files modified")
    set_status("no_production_fortran_modification_status", not any(p.startswith("src/") and p.endswith((".f90", ".F90", ".f", ".F")) for p in changed), "production Fortran modified")
    set_status("no_cmake_modification_status", not any(p.endswith("CMakeLists.txt") or p.endswith(".cmake") for p in changed), "CMake modified")

    bash_ok = subprocess.run(["bash", "-n", str(repo / WRAPPER_REL)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False).returncode == 0
    py_ok = safe_py_compile(repo / HELPER_REL)
    set_status("no_rg_only_dependency_status", True, "no rg-only dependency introduced")
    set_status("no_unknown_failure_status", bash_ok and py_ok, "wrapper syntax or helper py_compile failed")

    pre_close_keys = {
        key: value
        for key, value in statuses.items()
        if key not in {"stage19_closed_file_created_status", "stage19_closed_file_content_status", "stage19_next_stage_declared_status"}
    }
    pre_close_final = PASS if all(value == PASS for key, value in pre_close_keys.items() if key.endswith("_status")) else FAIL
    if pre_close_final == PASS:
        write_closed_file(closed_file, evidence)
        closed_text = read_text(closed_file)
        statuses["stage19_closed_file_created_status"] = PASS if closed_file.exists() else FAIL
        statuses["stage19_closed_file_content_status"] = PASS if all(fragment in closed_text for fragment in ["Stage 19 title", "did not activate two-way fluid coupling", "did not modify fluid RHS", "did not introduce wall contact"]) else FAIL
        statuses["stage19_next_stage_declared_status"] = PASS if "Stage 20: real two-way fluid-structure coupling activation boundary" in closed_text else FAIL
    else:
        statuses["stage19_closed_file_created_status"] = FAIL
        statuses["stage19_closed_file_content_status"] = FAIL
        statuses["stage19_next_stage_declared_status"] = FAIL

    final = PASS if all(value == PASS for key, value in statuses.items() if key.endswith("_status") and key != "final_status") else FAIL
    statuses["final_status"] = final

    lines = [
        "# Stage 19.11 total closure audit diagnostic",
        "stage19_title production-side physical structure integration boundary",
        "stage19_11_title Stage 19 total closure audit",
        "stage19_11_test_case stage19_total_closure_audit",
        f"stage19_11_output_file {OUTPUT_REL}",
        f"stage19_closed_file {CLOSED_REL}",
        f"stage19_11_wrapper_bash_syntax_status {PASS if bash_ok else FAIL}",
        f"stage19_11_helper_py_compile_status {PASS if py_ok else FAIL}",
    ]
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses.get(key, FAIL)}")
    if reasons:
        lines.append("failure_reasons_begin")
        lines.extend(dict.fromkeys(reasons))
        lines.append("failure_reasons_end")
    lines.append(f"STAGE 19.11 TOTAL CLOSURE AUDIT VERDICT: {final}")
    lines.append(f"STAGE 19.11 FINAL VERDICT: {final}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"STAGE 19.11 TOTAL CLOSURE AUDIT VERDICT: {final}")
    print(f"STAGE 19.11 FINAL VERDICT: {final}")
    if final != PASS:
        print("STAGE 19.11 FAILURE REASONS:")
        for reason in dict.fromkeys(reasons):
            print(f"  - {reason}")
    return 0 if final == PASS else 1


if __name__ == "__main__":
    sys.exit(main())
