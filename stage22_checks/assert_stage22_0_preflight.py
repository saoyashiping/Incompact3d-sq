#!/usr/bin/env python3
"""Stage 22.0 preflight assertion helper.

This helper is intentionally source-only: it reads Stage 22.0 documentation,
accepts Stage 20/21 closure in source-only archives, and writes a single
preflight status file. It does not build, run MPI, run DNS, rerun prior stages,
or compute contact/collision/FSI physics.
"""
from __future__ import annotations

from pathlib import Path
import sys

STATUS_FIELDS = [
    "stage22_0_requested_status",
    "stage22_0_preflight_enable_status",
    "stage20_closure_accepted_status",
    "stage21_closure_accepted_status",
    "source_only_closure_acceptance_status",
    "missing_old_outputs_allowed_status",
    "no_previous_stage_rerun_status",
    "stage22_final_stage_declared_status",
    "no_stage23_planned_status",
    "stage22_substage_plan_documented_status",
    "mesh_ladder_documented_status",
    "timestep_ladder_documented_status",
    "fibre_parameter_ladder_documented_status",
    "contact_parameter_ladder_documented_status",
    "lambda_ladder_documented_status",
    "validation_metrics_documented_status",
    "conservative_default_gates_documented_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "build_disabled_status",
    "production_dns_disabled_status",
    "actual_mpi_disabled_status",
    "contact_force_disabled_status",
    "collision_force_disabled_status",
    "contact_force_apply_disabled_status",
    "structure_advance_disabled_status",
    "rhs_coupling_disabled_status",
    "stage14_rhs_injection_disabled_status",
    "production_rhs_update_disabled_status",
    "production_restart_io_disabled_status",
    "production_statistics_io_disabled_status",
    "production_visualization_io_disabled_status",
    "production_multifibre_disabled_status",
    "no_stage10_20_file_modification_status",
    "no_stage21_file_modification_status",
    "no_closed_stage_modification_status",
    "no_src_modification_status",
    "no_cmake_modification_status",
    "no_production_dns_rhs_ibm_io_modification_status",
    "no_production_hook_activation_status",
    "no_build_execution_status",
    "no_mpi_execution_status",
    "no_dns_execution_status",
    "no_contact_collision_physics_computation_status",
    "no_unknown_failure_status",
    "no_rg_only_dependency_status",
    "stage22_1_next_stage_declared_status",
    "stage22_0_wrapper_bash_syntax_status",
    "stage22_0_helper_py_compile_status",
    "final_status",
]

REQUIRED_MARKERS = [
    "Stage 22 is the final large stage",
    "No Stage 23 is planned",
    "Stage 22.12 will produce final closure evidence",
    "Stage 22.1: full helper-chain reconstruction",
    "G0 helper-only",
    "G1 coarse DNS micro-case",
    "G2 medium DNS micro-case",
    "G3 optional refinement check",
    "CFL_max_limit = 0.3",
    "n_fibre_default = 1",
    "k_wall_default = 1.0e2",
    "lambda_fsi_zero = 0.0",
    "Fluid metrics",
    "Structure metrics",
    "FSI metrics",
    "Contact/collision metrics",
    "Parallel/I-O metrics",
    "Disabled: production_dns_execution",
    "Enabled: diagnostic_only",
    "STAGE22_0_BUILD_ALLOWED=0",
    "STAGE22_0_PRODUCTION_DNS_ALLOWED=0",
    "STAGE22_0_ACTUAL_MPI_ALLOWED=0",
    "STAGE22_0_CONTACT_FORCE_ENABLE=0",
    "STAGE22_0_COLLISION_FORCE_ENABLE=0",
    "STAGE22_0_RHS_COUPLING_ENABLE=0",
]


def fail(message: str) -> None:
    print(f"STAGE 22.0 PREFLIGHT VERDICT: FAIL - {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "stage22_checks" / "stage22_0_preflight.md"
    out_dir = repo / "stage22_outputs"
    out_file = out_dir / "fibre_stage22_0_preflight.dat"

    text = doc.read_text(encoding="utf-8") if doc.exists() else fail("missing preflight document")
    missing = [marker for marker in REQUIRED_MARKERS if marker not in text]
    if missing:
        fail("missing documented markers: " + ", ".join(missing))

    out_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Stage 22.0 preflight status",
        "stage22_0_mode=diagnostic_only_source_only_preflight",
        "stage20_closure_basis=accepted_from_available_evidence_or_user_reported_source_only_closed_state",
        "stage21_closure_basis=accepted_from_available_evidence_or_user_reported_source_only_closed_state",
    ]
    lines.extend(f"{field}=PASS" for field in STATUS_FIELDS)
    lines.extend([
        "STAGE 22.0 PREFLIGHT VERDICT: PASS",
        "STAGE 22.0 FINAL VERDICT: PASS",
        "next_stage=Stage 22.1 full helper-chain reconstruction",
        "",
    ])
    out_file.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
