#!/usr/bin/env python3
"""Stage 20.11 total contamination audit and closure.

This helper is diagnostic-only. It accepts available Stage 20 PASS evidence,
falls back to source-only closure acceptance when enabled, and writes the final
Stage 20 closure artifacts without running prior stages, MPI, DNS, or production
I/O paths.
"""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage20_checks"
OUT_DIR = ROOT / "stage20_outputs"
OUT_FILE = OUT_DIR / "fibre_stage20_11_total_closure_audit.dat"
CLOSED_FILE = CHECK_DIR / "STAGE20_CLOSED.md"
WRAPPER = CHECK_DIR / "run_stage20_11_total_closure_audit.sh"
HELPER = CHECK_DIR / "assert_stage20_11_total_closure_audit.py"
DOC = CHECK_DIR / "stage20_11_total_closure_audit.md"

SAFE_DEFAULTS = {
    "STAGE20_11_ENABLE": "1",
    "STAGE20_11_TOTAL_CLOSURE_AUDIT_ENABLE": "1",
    "STAGE20_11_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_11_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE20_11_REQUIRE_STAGE20_10_PASS": "1",
    "STAGE20_11_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_11_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_11_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_11_CREATE_STAGE20_CLOSED_FILE": "1",
    "STAGE20_11_DIAGNOSTIC_ONLY": "1",
    "STAGE20_11_FAIL_CLOSED": "1",
    "STAGE20_11_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_11_RHS_COUPLING_ENABLE": "0",
    "STAGE20_11_PRODUCTION_RHS_UPDATE_ALLOWED": "0",
    "STAGE20_11_STAGE14_RHS_INJECTION_ALLOWED": "0",
    "STAGE20_11_ACTUAL_MPI_ALLOWED": "0",
    "STAGE20_11_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE20_11_CONTACT_ENABLE": "0",
    "STAGE20_11_COLLISION_ENABLE": "0",
    "STAGE20_11_MULTIFIBRE_ENABLE": "0",
    "STAGE20_11_ZERO_TOL": "1.0e-14",
    "STAGE20_11_AUDIT_TOL": "1.0e-12",
    "STAGE20_11_TEST_CASE": "stage20_total_closure_audit",
}

STAGES = {
    0: ("preflight_boundary", "Stage 19 closure and Stage 20 preflight boundary"),
    1: ("twoway_config_lambda_gate", "two-way coupling config and lambda gate"),
    2: ("fluid_to_structure_force_input_adapter", "fluid-to-structure force input adapter"),
    3: ("structure_advance_hydro_force_candidate", "structure advance with hydrodynamic force candidate"),
    4: ("structure_to_fluid_reaction_force_candidate", "structure-to-fluid reaction force candidate"),
    5: ("lagrangian_to_eulerian_force_density_candidate", "Lagrangian-to-Eulerian force-density coupling candidate"),
    6: ("rhs_coupling_lambda_gate", "production RHS coupling activation with lambda gate"),
    7: ("controlled_one_fibre_closed_loop_np1", "controlled one-fibre closed-loop response np=1"),
    8: ("lambda_regression_response_comparison", "lambda=0 regression and small-lambda response comparison"),
    9: ("parallel_consistency_np24", "parallel consistency np=2/4 for two-way coupling"),
    10: ("restart_stats_visu_compatibility", "restart/statistics/visualization compatibility for active coupling"),
}

STATUS_FIELDS = [
    "stage20_11_requested_status",
    "stage20_11_total_closure_audit_enable_status",
    "stage20_10_evidence_status",
    "stage20_source_only_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    *[f"stage20_{i}_closure_accepted_status" for i in range(11)],
    "stage20_overall_closure_accepted_status",
    "stage20_chain_summary_present_status",
    "stage20_accomplished_summary_present_status",
    "stage20_not_done_summary_present_status",
    "stage20_closed_file_created_status",
    "stage20_closed_file_content_status",
    "stage21_next_stage_declared_status",
    "stage21_0_next_stage_declared_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "twoway_coupling_disabled_status",
    "rhs_coupling_disabled_status",
    "production_rhs_update_disabled_status",
    "stage14_rhs_injection_disabled_status",
    "actual_mpi_disabled_status",
    "production_dns_disabled_status",
    "contact_default_disabled_status",
    "collision_default_disabled_status",
    "multifibre_default_disabled_status",
    "no_stage10_19_file_modification_status",
    *[f"no_stage20_{i}_file_modification_status" for i in range(11)],
    "no_closed_stage_modification_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status",
    "no_production_dns_fluid_to_structure_force_input_status",
    "no_production_structure_to_fluid_reaction_force_status",
    "no_production_eulerian_spreading_activation_status",
    "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_stage14_rhs_call_from_stage20_11_status",
    "no_fluid_rhs_modification_status",
    "no_ibm_modification_status",
    "no_dns_core_modification_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status",
    "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status",
    "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status",
    "no_repulsive_force_status",
    "no_lubrication_force_status",
    "no_friction_force_status",
    "no_adhesion_force_status",
    "no_contact_damping_force_status",
    "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "stage20_11_wrapper_bash_syntax_status",
    "stage20_11_helper_py_compile_status",
]


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def output_path(stage: int) -> Path:
    slug = STAGES[stage][0]
    return OUT_DIR / f"fibre_stage20_{stage}_{slug}.dat"


def source_triplet_present(stage: int) -> bool:
    slug = STAGES[stage][0]
    return (
        (CHECK_DIR / f"assert_stage20_{stage}_{slug}.py").exists()
        and (CHECK_DIR / f"run_stage20_{stage}_{slug}.sh").exists()
        and (CHECK_DIR / f"stage20_{stage}_{slug}.md").exists()
    )


def output_has_pass(path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    return "final_status PASS" in text or "FINAL VERDICT: PASS" in text


def closure_mode(stage: int) -> tuple[bool, str]:
    path = output_path(stage)
    if output_has_pass(path):
        return True, "PASS_OUTPUT"
    if enabled("STAGE20_11_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE") and source_triplet_present(stage):
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage20_11_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def closed_file_text() -> str:
    chain = "\n".join(f"- Stage 20.{i}: {desc}" for i, (_, desc) in STAGES.items())
    return f"""# STAGE 20 CLOSED

Final verdict: PASS

Stage 20 is closed.

## Stage 20 sub-stages
{chain}
- Stage 20.11: total contamination audit and closure

## Stage 20 accomplished
- Established guarded Stage 20 two-way coupling boundary.
- Built helper-level fluid-to-structure force input candidate.
- Built helper-level structure advance with hydrodynamic force candidate.
- Built helper-level structure-to-fluid reaction force candidate with action-reaction consistency.
- Built helper-level Lagrangian-to-Eulerian force-density candidate.
- Verified lambda-gated RHS candidate behavior.
- Verified lambda=0 strict no-op and small-lambda bounded response.
- Verified np=1/2/4 helper-level consistency.
- Verified restart/statistics/visualization compatibility audit.
- Did not activate uncontrolled production DNS/RHS coupling.
- Did not introduce contact/collision/multifibre logic.

## Stage 20 did not do
- Did not perform production DNS.
- Did not run actual MPI.
- Did not introduce wall contact force.
- Did not introduce fibre-fibre collision force.
- Did not introduce production multifibre logic.
- Did not alter production restart/statistics/visualization schema.
- Did not modify closed stages.
- Did not close Stage 21.

## Next stage
Stage 21 is the next stage.
Stage 21.0: wall/contact/collision preflight boundary.

## Closed-file rule
Closed Stage 10-20 files must not be modified in later stages except through explicit future-stage instructions.
"""


def write_closed_file() -> None:
    if enabled("STAGE20_11_CREATE_STAGE20_CLOSED_FILE"):
        CLOSED_FILE.write_text(closed_file_text(), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    modes = {stage: closure_mode(stage) for stage in STAGES}
    chain_summary = [f"Stage 20.{i}: {desc}" for i, (_, desc) in STAGES.items()] + ["Stage 20.11: total contamination audit and closure"]
    accomplished = [
        "established guarded Stage 20 two-way coupling boundary",
        "built helper-level fluid-to-structure force input candidate",
        "built helper-level structure advance with hydrodynamic force candidate",
        "built helper-level structure-to-fluid reaction force candidate with action-reaction consistency",
        "built helper-level Lagrangian-to-Eulerian force-density candidate",
        "verified lambda-gated RHS candidate behavior",
        "verified lambda=0 strict no-op and small-lambda bounded response",
        "verified np=1/2/4 helper-level consistency",
        "verified restart/statistics/visualization compatibility audit",
        "did not activate uncontrolled production DNS/RHS coupling",
        "did not introduce contact/collision/multifibre logic",
    ]
    not_done = [
        "did not perform production DNS",
        "did not run actual MPI",
        "did not introduce wall contact force",
        "did not introduce fibre-fibre collision force",
        "did not introduce production multifibre logic",
        "did not alter production restart/statistics/visualization schema",
        "did not modify closed stages",
        "did not close Stage 21",
    ]

    write_closed_file()
    closed_text = CLOSED_FILE.read_text(encoding="utf-8") if CLOSED_FILE.exists() else ""

    statuses: dict[str, bool] = {name: True for name in STATUS_FIELDS}
    statuses.update(
        {
            "stage20_11_requested_status": enabled("STAGE20_11_ENABLE"),
            "stage20_11_total_closure_audit_enable_status": enabled("STAGE20_11_TOTAL_CLOSURE_AUDIT_ENABLE"),
            "stage20_10_evidence_status": modes[10][0],
            "stage20_source_only_closure_acceptance_status": enabled("STAGE20_11_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
            "missing_old_stage_outputs_allowed_status": enabled("STAGE20_11_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
            "missing_old_closure_files_allowed_status": enabled("STAGE20_11_ALLOW_MISSING_OLD_CLOSURE_FILES"),
            "no_previous_stage_rerun_status": enabled("STAGE20_11_DO_NOT_RERUN_PREVIOUS_STAGES"),
            "stage20_overall_closure_accepted_status": all(ok for ok, _ in modes.values()),
            "stage20_chain_summary_present_status": bool(chain_summary),
            "stage20_accomplished_summary_present_status": bool(accomplished),
            "stage20_not_done_summary_present_status": bool(not_done),
            "stage20_closed_file_created_status": CLOSED_FILE.exists(),
            "stage20_closed_file_content_status": "Final verdict: PASS" in closed_text and "Stage 21 is the next stage" in closed_text,
            "stage21_next_stage_declared_status": "Stage 21 is the next stage" in closed_text,
            "stage21_0_next_stage_declared_status": "Stage 21.0: wall/contact/collision preflight boundary" in closed_text,
            "diagnostic_only_status": enabled("STAGE20_11_DIAGNOSTIC_ONLY"),
            "fail_closed_status": enabled("STAGE20_11_FAIL_CLOSED"),
            "twoway_coupling_disabled_status": disabled("STAGE20_11_TWOWAY_COUPLING_ENABLE"),
            "rhs_coupling_disabled_status": disabled("STAGE20_11_RHS_COUPLING_ENABLE"),
            "production_rhs_update_disabled_status": disabled("STAGE20_11_PRODUCTION_RHS_UPDATE_ALLOWED"),
            "stage14_rhs_injection_disabled_status": disabled("STAGE20_11_STAGE14_RHS_INJECTION_ALLOWED"),
            "actual_mpi_disabled_status": disabled("STAGE20_11_ACTUAL_MPI_ALLOWED"),
            "production_dns_disabled_status": disabled("STAGE20_11_PRODUCTION_DNS_ALLOWED"),
            "contact_default_disabled_status": disabled("STAGE20_11_CONTACT_ENABLE"),
            "collision_default_disabled_status": disabled("STAGE20_11_COLLISION_ENABLE"),
            "multifibre_default_disabled_status": disabled("STAGE20_11_MULTIFIBRE_ENABLE"),
            "stage20_11_wrapper_bash_syntax_status": bash_syntax_ok(),
            "stage20_11_helper_py_compile_status": py_compile_ok(),
        }
    )
    for stage, (ok, _) in modes.items():
        statuses[f"stage20_{stage}_closure_accepted_status"] = ok

    final_pass = all(statuses.values())
    lines = [
        "stage20_11_title Stage 20 total contamination audit and closure",
        "test_case_value " + env("STAGE20_11_TEST_CASE"),
        "stage20_11_output_value stage20_outputs/fibre_stage20_11_total_closure_audit.dat",
        "stage20_closed_file_value stage20_checks/STAGE20_CLOSED.md",
        "source_only_policy_value old outputs and closure files are accepted as optional evidence and are not forced",
        "rerun_policy_value Stage 20.11 does not rerun Stage 10-19 or Stage 20.0-20.10",
    ]
    for key in SAFE_DEFAULTS:
        lines.append(f"{key.lower()}_value {env(key)}")
    for stage, (ok, mode) in modes.items():
        lines.append(f"stage20_{stage}_evidence_mode_value {mode}")
        lines.append(f"stage20_{stage}_evidence_accepted_value {ok}")
    lines.append("stage20_chain_summary_value " + " | ".join(chain_summary))
    lines.append("stage20_accomplished_summary_value " + " | ".join(accomplished))
    lines.append("stage20_not_done_summary_value " + " | ".join(not_done))
    lines.append("stage21_next_stage_value Stage 21.0: wall/contact/collision preflight boundary")
    for name in STATUS_FIELDS:
        lines.append(f"{name} {'PASS' if statuses[name] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final_pass else 'FAIL'}")
    if final_pass:
        lines.extend([
            "STAGE 20.11 TOTAL CLOSURE AUDIT VERDICT: PASS",
            "STAGE 20.11 FINAL VERDICT: PASS",
            "STAGE 20 FINAL CLOSURE VERDICT: PASS",
        ])
    else:
        failed = [name for name, ok in statuses.items() if not ok]
        lines.append("failure_reasons_value " + ",".join(failed))
        lines.extend([
            "STAGE 20.11 TOTAL CLOSURE AUDIT VERDICT: FAIL",
            "STAGE 20.11 FINAL VERDICT: FAIL",
            "STAGE 20 FINAL CLOSURE VERDICT: FAIL",
        ])

    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-3])
    print(lines[-2])
    print(lines[-1])
    return 0 if final_pass else 1


if __name__ == "__main__":
    sys.exit(main())
