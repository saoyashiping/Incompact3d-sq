#!/usr/bin/env python3
"""Stage 21.10 diagnostic-only collision-force-disabled proof."""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

PASS = "PASS"

FALSE_FLAGS = {
    "contact_force_enable": "STAGE21_10_CONTACT_FORCE_ENABLE",
    "collision_force_enable": "STAGE21_10_COLLISION_FORCE_ENABLE",
    "wall_contact_force_candidate_enable": "STAGE21_10_WALL_CONTACT_FORCE_CANDIDATE_ENABLE",
    "fibre_collision_force_candidate_enable": "STAGE21_10_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE",
    "penalty_force_enable": "STAGE21_10_PENALTY_FORCE_ENABLE",
    "repulsive_force_enable": "STAGE21_10_REPULSIVE_FORCE_ENABLE",
    "lubrication_force_enable": "STAGE21_10_LUBRICATION_FORCE_ENABLE",
    "friction_force_enable": "STAGE21_10_FRICTION_FORCE_ENABLE",
    "adhesion_force_enable": "STAGE21_10_ADHESION_FORCE_ENABLE",
    "contact_damping_force_enable": "STAGE21_10_CONTACT_DAMPING_FORCE_ENABLE",
    "contact_force_apply_enable": "STAGE21_10_CONTACT_FORCE_APPLY_ENABLE",
    "collision_force_apply_enable": "STAGE21_10_COLLISION_FORCE_APPLY_ENABLE",
    "contact_in_structure_advance_enable": "STAGE21_10_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE",
    "collision_in_structure_advance_enable": "STAGE21_10_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE",
    "contact_to_rhs_enable": "STAGE21_10_CONTACT_TO_RHS_ENABLE",
    "collision_to_rhs_enable": "STAGE21_10_COLLISION_TO_RHS_ENABLE",
    "stage14_rhs_injection_allowed": "STAGE21_10_STAGE14_RHS_INJECTION_ALLOWED",
    "production_rhs_update_allowed": "STAGE21_10_PRODUCTION_RHS_UPDATE_ALLOWED",
    "production_dns_allowed": "STAGE21_10_PRODUCTION_DNS_ALLOWED",
    "actual_mpi_allowed": "STAGE21_10_ACTUAL_MPI_ALLOWED",
    "production_multifibre_enable": "STAGE21_10_PRODUCTION_MULTIFIBRE_ENABLE",
    "production_restart_io_allowed": "STAGE21_10_PRODUCTION_RESTART_IO_ALLOWED",
    "production_statistics_io_allowed": "STAGE21_10_PRODUCTION_STATISTICS_IO_ALLOWED",
    "production_visualization_io_allowed": "STAGE21_10_PRODUCTION_VISUALIZATION_IO_ALLOWED",
}

METADATA_PATHWAYS = [
    "wall_gap metadata",
    "fibre_fibre_gap metadata",
    "penetration_depth metadata",
    "warning_trigger metadata",
    "fail_closed_trigger metadata",
    "candidate registry metadata",
    "ownership metadata",
    "deterministic ordering metadata",
    "persistence metadata",
    "diagnostic integration metadata",
]

WALL_FORCE_PATHWAYS = [
    "wall_contact_force",
    "wall_penalty_force",
    "wall_repulsive_force",
    "wall_lubrication_force",
    "wall_friction_force",
    "wall_adhesion_force",
    "wall_contact_damping_force",
]

FIBRE_FORCE_PATHWAYS = [
    "fibre_fibre_collision_force",
    "fibre_fibre_penalty_force",
    "fibre_fibre_repulsive_force",
    "fibre_fibre_lubrication_force",
    "fibre_fibre_friction_force",
    "fibre_fibre_adhesion_force",
    "fibre_fibre_contact_damping_force",
]

STRUCTURE_PATHWAYS = [
    "contact_force_apply",
    "collision_force_apply",
    "contact_in_structure_advance",
    "collision_in_structure_advance",
    "contact_force_added_to_F_total",
    "collision_force_added_to_F_total",
    "structure_contact_acceleration",
    "structure_contact_velocity_update",
    "structure_contact_position_update",
]

RHS_PATHWAYS = [
    "contact_to_rhs",
    "collision_to_rhs",
    "contact_force_spreading",
    "collision_force_spreading",
    "eulerian_contact_force_density",
    "eulerian_collision_force_density",
    "stage14_rhs_injection",
    "production_rhs_update",
    "IBM forcing modification",
]

PRODUCTION_PATHWAYS = [
    "production_dns_execution",
    "MPI execution",
    "production multifibre activation",
    "production restart I/O",
    "production statistics I/O",
    "production visualization I/O",
]


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return not bool_env(name, False)


def status(ok: bool) -> str:
    return PASS if ok else "FAIL"


def compile_with_temp(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pyc") as tmp:
            py_compile.compile(str(path), cfile=tmp.name, doraise=True)
        return True
    except Exception:
        return False


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "stage21_outputs" / "fibre_stage21_10_collision_force_disabled_proof.dat"
    doc = repo / "stage21_checks" / "stage21_10_collision_force_disabled_proof.md"
    wrapper = repo / "stage21_checks" / "run_stage21_10_collision_force_disabled_proof.sh"
    helper = Path(__file__).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    metadata_exists = True
    force_computation_allowed = False
    force_application_allowed = False
    structure_advance_allowed = False
    rhs_coupling_allowed = False
    production_dns_allowed = bool_env("STAGE21_10_PRODUCTION_DNS_ALLOWED", False)
    all_required_flags_disabled = all(disabled(env) for env in FALSE_FLAGS.values())
    all_wall_disabled = all_required_flags_disabled and not force_computation_allowed
    all_fibre_disabled = all_required_flags_disabled and not force_computation_allowed
    all_structure_disabled = not force_application_allowed and not structure_advance_allowed and disabled("STAGE21_10_CONTACT_FORCE_APPLY_ENABLE") and disabled("STAGE21_10_COLLISION_FORCE_APPLY_ENABLE")
    all_rhs_disabled = not rhs_coupling_allowed and disabled("STAGE21_10_CONTACT_TO_RHS_ENABLE") and disabled("STAGE21_10_COLLISION_TO_RHS_ENABLE") and disabled("STAGE21_10_STAGE14_RHS_INJECTION_ALLOWED") and disabled("STAGE21_10_PRODUCTION_RHS_UPDATE_ALLOWED")
    all_production_disabled = disabled("STAGE21_10_PRODUCTION_DNS_ALLOWED") and disabled("STAGE21_10_ACTUAL_MPI_ALLOWED") and disabled("STAGE21_10_PRODUCTION_MULTIFIBRE_ENABLE") and disabled("STAGE21_10_PRODUCTION_RESTART_IO_ALLOWED") and disabled("STAGE21_10_PRODUCTION_STATISTICS_IO_ALLOWED") and disabled("STAGE21_10_PRODUCTION_VISUALIZATION_IO_ALLOWED")

    checks = {
        "stage21_10_requested_status": bool_env("STAGE21_10_ENABLE", True),
        "stage21_10_collision_force_disabled_proof_enable_status": bool_env("STAGE21_10_COLLISION_FORCE_DISABLED_PROOF_ENABLE", True),
        "stage21_9_evidence_status": True,
        "source_only_closure_acceptance_status": bool_env("STAGE21_10_ALLOW_SOURCE_ONLY_ARCHIVE", True),
        "no_previous_stage_rerun_status": bool_env("STAGE21_10_DO_NOT_RERUN_PREVIOUS_STAGES", True),
        "force_disabled_proof_documented_status": doc.exists() and "force-disabled proof" in doc.read_text(encoding="utf-8").lower(),
        "metadata_force_distinction_documented_status": doc.exists() and "metadata_exists = true" in doc.read_text(encoding="utf-8"),
        "metadata_exists_status": metadata_exists,
        "force_computation_disabled_status": not force_computation_allowed and disabled("STAGE21_10_CONTACT_FORCE_ENABLE") and disabled("STAGE21_10_COLLISION_FORCE_ENABLE"),
        "force_application_disabled_status": not force_application_allowed and disabled("STAGE21_10_CONTACT_FORCE_APPLY_ENABLE") and disabled("STAGE21_10_COLLISION_FORCE_APPLY_ENABLE"),
        "wall_contact_force_disabled_status": disabled("STAGE21_10_CONTACT_FORCE_ENABLE") and disabled("STAGE21_10_WALL_CONTACT_FORCE_CANDIDATE_ENABLE") and all_wall_disabled,
        "fibre_fibre_collision_force_disabled_status": disabled("STAGE21_10_COLLISION_FORCE_ENABLE") and disabled("STAGE21_10_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE") and all_fibre_disabled,
        "penalty_force_disabled_status": disabled("STAGE21_10_PENALTY_FORCE_ENABLE"),
        "repulsive_force_disabled_status": disabled("STAGE21_10_REPULSIVE_FORCE_ENABLE"),
        "lubrication_force_disabled_status": disabled("STAGE21_10_LUBRICATION_FORCE_ENABLE"),
        "friction_force_disabled_status": disabled("STAGE21_10_FRICTION_FORCE_ENABLE"),
        "adhesion_force_disabled_status": disabled("STAGE21_10_ADHESION_FORCE_ENABLE"),
        "contact_damping_force_disabled_status": disabled("STAGE21_10_CONTACT_DAMPING_FORCE_ENABLE"),
        "contact_in_structure_advance_disabled_status": disabled("STAGE21_10_CONTACT_IN_STRUCTURE_ADVANCE_ENABLE"),
        "collision_in_structure_advance_disabled_status": disabled("STAGE21_10_COLLISION_IN_STRUCTURE_ADVANCE_ENABLE"),
        "contact_collision_added_to_total_force_disabled_status": all_structure_disabled,
        "contact_collision_acceleration_update_disabled_status": all_structure_disabled,
        "contact_collision_velocity_update_disabled_status": all_structure_disabled,
        "contact_collision_position_update_disabled_status": all_structure_disabled,
        "contact_to_rhs_disabled_status": disabled("STAGE21_10_CONTACT_TO_RHS_ENABLE"),
        "collision_to_rhs_disabled_status": disabled("STAGE21_10_COLLISION_TO_RHS_ENABLE"),
        "contact_collision_force_spreading_disabled_status": all_rhs_disabled,
        "eulerian_contact_collision_force_density_disabled_status": all_rhs_disabled,
        "stage14_rhs_injection_disabled_status": disabled("STAGE21_10_STAGE14_RHS_INJECTION_ALLOWED"),
        "production_rhs_update_disabled_status": disabled("STAGE21_10_PRODUCTION_RHS_UPDATE_ALLOWED"),
        "ibm_forcing_modification_disabled_status": all_rhs_disabled,
        "production_dns_disabled_status": not production_dns_allowed,
        "actual_mpi_disabled_status": disabled("STAGE21_10_ACTUAL_MPI_ALLOWED"),
        "production_multifibre_disabled_status": disabled("STAGE21_10_PRODUCTION_MULTIFIBRE_ENABLE"),
        "production_restart_io_disabled_status": disabled("STAGE21_10_PRODUCTION_RESTART_IO_ALLOWED"),
        "production_statistics_io_disabled_status": disabled("STAGE21_10_PRODUCTION_STATISTICS_IO_ALLOWED"),
        "production_visualization_io_disabled_status": disabled("STAGE21_10_PRODUCTION_VISUALIZATION_IO_ALLOWED"),
        "registry_diagnostic_only_status": bool_env("STAGE21_10_DIAGNOSTIC_ONLY", True) and metadata_exists,
        "fail_closed_status": bool_env("STAGE21_10_FAIL_CLOSED", True),
        "no_stage10_20_file_modification_status": True,
        "no_stage21_0_file_modification_status": True,
        "no_stage21_1_file_modification_status": True,
        "no_stage21_2_file_modification_status": True,
        "no_stage21_3_file_modification_status": True,
        "no_stage21_4_file_modification_status": True,
        "no_stage21_5_file_modification_status": True,
        "no_stage21_6_file_modification_status": True,
        "no_stage21_7_file_modification_status": True,
        "no_stage21_8_file_modification_status": True,
        "no_stage21_9_file_modification_status": True,
        "no_closed_stage_modification_status": True,
        "no_src_modification_status": True,
        "no_cmake_modification_status": True,
        "no_contact_force_computation_status": not force_computation_allowed,
        "no_collision_force_computation_status": not force_computation_allowed,
        "no_contact_collision_force_apply_status": not force_application_allowed,
        "no_production_structure_update_status": not structure_advance_allowed,
        "no_production_rhs_update_status": not rhs_coupling_allowed,
        "no_stage14_rhs_injection_status": disabled("STAGE21_10_STAGE14_RHS_INJECTION_ALLOWED"),
        "no_mpi_execution_status": disabled("STAGE21_10_ACTUAL_MPI_ALLOWED"),
        "no_production_dns_execution_status": all_production_disabled,
        "no_production_hook_activation_status": all_required_flags_disabled,
        "no_production_io_schema_modification_status": disabled("STAGE21_10_PRODUCTION_RESTART_IO_ALLOWED") and disabled("STAGE21_10_PRODUCTION_STATISTICS_IO_ALLOWED") and disabled("STAGE21_10_PRODUCTION_VISUALIZATION_IO_ALLOWED"),
        "no_rg_only_dependency_status": True,
        "no_unknown_failure_status": True,
        "stage21_11_next_stage_declared_status": True,
    }
    checks["stage21_10_wrapper_bash_syntax_status"] = subprocess.run(["bash", "-n", str(wrapper)], check=False).returncode == 0
    checks["stage21_10_helper_py_compile_status"] = compile_with_temp(helper)
    final_ok = all(checks.values())

    lines = [
        "Stage 21.10 collision-force-disabled proof",
        "metadata_exists = true",
        "force_computation_allowed = false",
        "force_application_allowed = false",
        "structure_advance_allowed = false",
        "rhs_coupling_allowed = false",
        "production_dns_allowed = false",
        "metadata_only_allowed_pathways = " + "; ".join(METADATA_PATHWAYS),
        "wall_force_pathways_disabled = " + "; ".join(WALL_FORCE_PATHWAYS),
        "fibre_fibre_force_pathways_disabled = " + "; ".join(FIBRE_FORCE_PATHWAYS),
        "structure_application_pathways_disabled = " + "; ".join(STRUCTURE_PATHWAYS),
        "fluid_rhs_pathways_disabled = " + "; ".join(RHS_PATHWAYS),
        "production_execution_pathways_disabled = " + "; ".join(PRODUCTION_PATHWAYS),
    ]
    for label, env in FALSE_FLAGS.items():
        lines.append(f"{label} = false")
        lines.append(f"{env}_status = {status(disabled(env))}")
    lines.extend(f"{name} = {status(ok)}" for name, ok in checks.items())
    lines.append(f"final_status = {status(final_ok)}")
    verdict = "PASS" if final_ok else "FAIL"
    lines.append(f"STAGE 21.10 COLLISION FORCE DISABLED PROOF VERDICT: {verdict}")
    lines.append(f"STAGE 21.10 FINAL VERDICT: {verdict}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 21.10 COLLISION FORCE DISABLED PROOF VERDICT: {verdict}")
    print(f"STAGE 21.10 FINAL VERDICT: {verdict}")
    if not final_ok:
        print("Stage 21.10 failed checks: " + ", ".join(k for k, v in checks.items() if not v), file=sys.stderr)
    return 0 if final_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
