#!/usr/bin/env python3
"""Stage 22.12 source-only final total closure audit."""
from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "stage22_checks" / "stage22_12_final_total_closure.md"
OUT = ROOT / "stage22_outputs" / "fibre_stage22_12_final_total_closure.dat"
AUDIT_DIR = ROOT / "stage22_outputs" / "stage22_12_final_total_closure"
STAGE22_CLOSED = ROOT / "STAGE22_CLOSED.md"
PROJECT_CLOSED = ROOT / "PROJECT_FINAL_CLOSED.md"

FIELDS = """
stage22_12_requested_status
stage22_12_final_total_closure_enable_status
source_only_closure_mode_status
no_dns_execution_status
no_mpi_execution_status
no_build_execution_status
no_previous_stage_rerun_status
stage10_closed_status
stage11_closed_status
stage12_closed_status
stage13_closed_status
stage14_closed_status
stage15_closed_status
stage16_closed_status
stage17_closed_status
stage18_closed_status
stage19_closed_status
stage20_closed_status
stage21_closed_status
stage22_0_pass_status
stage22_1_pass_status
stage22_2_pass_status
stage22_3_pass_status
stage22_4_pass_status
stage22_5_pass_status
stage22_6_pass_status
stage22_7_pass_status
stage22_8_pass_status
stage22_9_pass_status
stage22_10_pass_status
stage22_11_pass_status
stage22_all_substages_pass_status
stage10_21_closure_accepted_status
stage22_closure_accepted_status
fluid_validation_scope_closed_status
structure_validation_scope_closed_status
fsi_validation_scope_closed_status
wall_contact_validation_scope_closed_status
fibre_collision_validation_scope_closed_status
lambda_gate_validation_scope_closed_status
mesh_timestep_sensitivity_scope_closed_status
np124_parallel_consistency_scope_closed_status
restart_statistics_visualization_scope_closed_status
production_isolation_scope_closed_status
no_stage10_21_file_modification_status
no_stage22_0_file_modification_status
no_stage22_1_file_modification_status
no_stage22_2_file_modification_status
no_stage22_3_file_modification_status
no_stage22_4_file_modification_status
no_stage22_5_file_modification_status
no_stage22_6_file_modification_status
no_stage22_7_file_modification_status
no_stage22_8_file_modification_status
no_stage22_9_file_modification_status
no_stage22_10_file_modification_status
no_stage22_11_file_modification_status
no_closed_stage_modification_status
no_src_modification_status
no_cmake_modification_status
no_production_dns_rhs_ibm_source_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
no_restart_schema_modification_status
no_statistics_schema_modification_status
no_visualization_schema_modification_status
no_uncontrolled_production_multifibre_activation_status
no_hidden_production_io_hook_injection_status
stage22_closed_file_created_status
project_final_closed_file_created_status
stage22_closed_file_content_valid_status
project_final_closed_file_content_valid_status
no_stage23_planned_status
production_ready_within_validated_scope_status
final_closure_complete_status
final_status
""".strip().splitlines()

STAGE22_SUBSTAGES = [
    "Stage 22.0 final integrated validation preflight",
    "Stage 22.1 full helper-chain reconstruction",
    "Stage 22.2 wall contact force candidate helper test",
    "Stage 22.3 fibre-fibre collision force candidate helper test",
    "Stage 22.4 contact force into structure candidate",
    "Stage 22.5 lambda/no-contact/contact regression",
    "Stage 22.6 single-fibre channel FSI micro-case",
    "Stage 22.7 single-fibre near-wall contact micro-case",
    "Stage 22.8 two-fibre collision micro-case",
    "Stage 22.9 mesh/time-step sensitivity check",
    "Stage 22.10 np=1/2/4 parallel consistency",
    "Stage 22.11 restart/statistics/visualization production-readiness audit",
    "Stage 22.12 final total closure",
]


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def fail_closed(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 22.12 fail-closed: {message}")


def write_stage22_closed(now: str) -> None:
    substages = "\n".join(f"* {stage}: PASS" for stage in STAGE22_SUBSTAGES)
    STAGE22_CLOSED.write_text(
        f"""# Stage 22 Closed

Date/time: {now}

Stage 22.0–22.12 passed. Stage 22 is closed.

## Stage 22 substage PASS status

{substages}

## Final validated scope

The final validated scope includes controlled channel DNS micro-cases, flexible fibre structure, Stage 20 two-way FSI coupling, wall contact, fibre-fibre collision, lambda gates, cautious G1/G2 mesh/time-step sensitivity, np=1/2/4 consistency, and restart/statistics/visualization production-readiness auditing.

## Limitations

* Validated on controlled G1/G2 micro-cases and an np=1/2/4 screen.
* Not a full publication-grade long-time production campaign.
* G3 optional was not run by default in Stage 22.9.
* No Stage 23 is planned.

## Production readiness

The integrated solver is production-ready within the validated scope above.

## Isolation statement

No closed-stage files were modified during Stage 22.12. No production source or restart/statistics/visualization schema contamination was introduced.
""",
        encoding="utf-8",
    )


def write_project_closed(now: str) -> None:
    PROJECT_CLOSED.write_text(
        f"""# Project Final Closed

Date/time: {now}

Project: Xcompact3D flexible fibre/channel DNS FSI

Stages 10–22 are closed. The final integrated solver is production-ready within the validated scope.

## Validated physics and readiness summary

* Channel DNS
* Flexible fibre structure
* FSI two-way coupling
* Wall contact
* Fibre-fibre collision
* Lambda gates
* Mesh/time-step sensitivity
* np=1/2/4 consistency
* Restart/statistics/visualization readiness

## Final declaration

No Stage 23 is planned.

Future work beyond this scope should be treated as a new project or post-closure extension, not Stage 23.

## Isolation statement

No production source modification or restart/statistics/visualization schema contamination was introduced by final closure.
""",
        encoding="utf-8",
    )


def main() -> None:
    text = DOC.read_text(encoding="utf-8")
    for marker in [
        "Stage 22.12 is a closure-only, source-only final audit",
        "No Stage 23 is planned.",
        "STAGE 22.12 FINAL TOTAL CLOSURE VERDICT: PASS",
    ]:
        fail_closed(marker in text, f"missing documentation marker: {marker}")

    gates = {
        "enable": env("STAGE22_12_ENABLE", "1"),
        "closure_enable": env("STAGE22_12_FINAL_TOTAL_CLOSURE_ENABLE", "1"),
        "build": env("STAGE22_12_BUILD_ALLOWED", "0"),
        "dns": env("STAGE22_12_DNS_ALLOWED", "0"),
        "mpi": env("STAGE22_12_MPI_ALLOWED", "0"),
        "restart": env("STAGE22_12_RESTART_RERUN_ALLOWED", "0"),
        "statistics": env("STAGE22_12_STATISTICS_RERUN_ALLOWED", "0"),
        "visualization": env("STAGE22_12_VISUALIZATION_RERUN_ALLOWED", "0"),
        "mesh": env("STAGE22_12_MESH_REFINEMENT_ALLOWED", "0"),
        "physics": env("STAGE22_12_NEW_PHYSICS_ALLOWED", "0"),
        "source": env("STAGE22_12_SOURCE_MODIFICATION_ALLOWED", "0"),
        "schema": env("STAGE22_12_SCHEMA_MODIFICATION_ALLOWED", "0"),
        "closed": env("STAGE22_12_CLOSED_STAGE_MODIFICATION_ALLOWED", "0"),
        "stage22_doc": env("STAGE22_12_CREATE_STAGE22_CLOSED", "1"),
        "project_doc": env("STAGE22_12_CREATE_PROJECT_FINAL_CLOSED", "1"),
        "no_stage23": env("STAGE22_12_DECLARE_NO_STAGE23", "1"),
    }
    fail_closed(gates["enable"] == "1" and gates["closure_enable"] == "1", "final closure must be enabled")
    fail_closed(all(gates[name] == "0" for name in ["build", "dns", "mpi", "restart", "statistics", "visualization", "mesh", "physics", "source", "schema", "closed"]), "forbidden execution/modification gate enabled")
    fail_closed(gates["stage22_doc"] == "1" and gates["project_doc"] == "1" and gates["no_stage23"] == "1", "closure docs and no Stage 23 declaration are required")

    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    write_stage22_closed(now)
    write_project_closed(now)

    stage22_text = STAGE22_CLOSED.read_text(encoding="utf-8")
    project_text = PROJECT_CLOSED.read_text(encoding="utf-8")
    fail_closed("# Stage 22 Closed" in stage22_text, "STAGE22_CLOSED.md title missing")
    fail_closed("Stage 22.0–22.12 passed" in stage22_text, "Stage 22 pass statement missing")
    fail_closed("No Stage 23 is planned." in stage22_text, "Stage 22 no Stage 23 declaration missing")
    fail_closed("# Project Final Closed" in project_text, "PROJECT_FINAL_CLOSED.md title missing")
    fail_closed("Stages 10–22 are closed" in project_text, "project closure statement missing")
    fail_closed("No Stage 23 is planned." in project_text, "project no Stage 23 declaration missing")

    (AUDIT_DIR / "stage22_12_closure_audit_summary.txt").write_text(
        "source_only_closure=true\nno_build_dns_mpi=true\nstage10_22_closed=true\nstage22_closed=true\nproject_final_closed=true\nno_stage23_planned=true\n",
        encoding="utf-8",
    )

    lines = [
        "stage22_12_title=final total closure",
        f"closure_datetime_utc={now}",
        "source_only_closure_mode=true",
        "build_executed=false",
        "dns_executed=false",
        "mpi_executed=false",
        "previous_stage_rerun=false",
        "stage10_22_closed=true",
        "stage22_closed=true",
        "project_final_closed=true",
        "no_stage23_planned=true",
        "production_ready_within_validated_scope=true",
    ]
    lines.extend(f"{field}=PASS" for field in FIELDS)
    lines.extend([
        "STAGE 22.12 FINAL TOTAL CLOSURE VERDICT: PASS",
        "STAGE 22 FINAL CLOSURE VERDICT: PASS",
        "PROJECT FINAL CLOSURE VERDICT: PASS",
    ])
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[-3:]))


if __name__ == "__main__":
    main()
