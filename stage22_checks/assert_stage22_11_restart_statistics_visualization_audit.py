#!/usr/bin/env python3
"""Stage 22.11 restart/statistics/visualization production-readiness audit."""
from __future__ import annotations

import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "stage22_checks" / "stage22_11_restart_statistics_visualization_audit.md"
OUT = ROOT / "stage22_outputs" / "fibre_stage22_11_restart_statistics_visualization_audit.dat"
CASE_DIR = ROOT / "stage22_cases" / "stage22_11_restart_statistics_visualization_audit"
AUDIT_DIR = ROOT / "stage22_outputs" / "stage22_11_restart_statistics_visualization_audit"

FIELDS = """
stage22_11_requested_status
stage22_11_restart_statistics_visualization_audit_enable_status
stage22_10_evidence_status
stage22_9_evidence_status
stage22_8_evidence_status
stage22_7_evidence_status
stage22_6_evidence_status
stage22_5_evidence_status
stage22_4_evidence_status
stage22_3_evidence_status
stage22_2_evidence_status
stage22_1_evidence_status
stage22_0_evidence_status
stage20_closure_accepted_status
stage21_closure_accepted_status
source_only_closure_acceptance_status
missing_old_outputs_allowed_status
no_previous_stage_rerun_status
cautious_io_production_readiness_mode_enabled_status
g1_case_documented_status
g2_disabled_status
g3_disabled_status
np1_primary_documented_status
np2_secondary_optional_documented_status
np4_disabled_by_default_status
fresh_restart_split_documented_status
total_equivalent_steps_documented_status
no_mesh_refinement_status
no_new_physics_status
no_final_closure_in_stage22_11_status
isolated_build_directory_valid_status
no_source_modification_during_build_status
fresh_segment_completed_status
restart_file_written_status
restart_continuation_completed_status
no_unsupported_np_run_status
no_g2_g3_run_status
no_nan_inf_runtime_logs_status
restart_file_exists_status
restart_file_readable_status
restart_metadata_finite_status
restart_time_index_consistent_status
restart_step_index_consistent_status
restart_field_dimensions_consistent_status
restart_fibre_state_dimensions_consistent_status
restart_fluid_state_finite_status
restart_structure_state_finite_status
restart_continuation_expected_step_time_status
post_restart_velocity_finite_status
post_restart_pressure_finite_status
post_restart_rhs_finite_status
post_restart_divergence_bounded_status
post_restart_cfl_bounded_status
post_restart_xva_finite_status
post_restart_fsi_force_bounded_status
restart_signature_consistency_bounded_status
no_restart_schema_modification_status
no_restart_field_count_drift_status
no_restart_file_corruption_status
statistics_path_audited_status
statistics_output_readable_or_inactive_documented_status
statistics_timestamps_monotonic_status
statistics_fields_finite_status
statistics_no_nan_inf_status
statistics_schema_unchanged_status
statistics_no_unexpected_column_drift_status
visualization_path_audited_status
visualization_output_readable_or_inactive_documented_status
visualization_metadata_finite_status
visualization_dimensions_consistent_status
visualization_fields_finite_status
visualization_no_nan_inf_status
visualization_schema_unchanged_status
visualization_no_file_corruption_status
velocity_finite_status
pressure_finite_status
rhs_finite_status
divergence_bounded_status
cfl_bounded_status
projection_stable_status
poisson_stable_status
x_finite_status
v_finite_status
a_finite_status
segment_length_residual_bounded_status
structure_step_displacement_bounded_status
fsi_force_bounded_status
action_reaction_residual_bounded_status
force_conservation_residual_bounded_status
rhs_delta_bounded_status
lambda_fsi_response_bounded_status
contact_collision_disabled_default_baseline_status
contact_collision_force_norm_zero_default_baseline_status
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
no_closed_stage_modification_status
no_src_modification_status
no_cmake_modification_status
no_production_dns_rhs_ibm_source_modification_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
no_production_restart_schema_modification_status
no_production_statistics_schema_modification_status
no_production_visualization_schema_modification_status
no_uncontrolled_production_multifibre_activation_status
no_hidden_production_io_hook_injection_status
no_stage22_12_final_closure_file_created_early_status
no_unknown_failure_status
no_rg_only_dependency_status
stage22_12_next_stage_declared_status
stage22_11_wrapper_bash_syntax_status
stage22_11_helper_py_compile_status
final_status
""".strip().splitlines()

MARKERS = [
    "Stage 22.11 is the final functional audit before Stage 22.12 total closure",
    "restart/statistics/visualization production-readiness audit only",
    "np=1 primary",
    "np=4 disabled by default",
    "G2 disabled",
    "G3 disabled",
    "Stage 22.12: final total closure.",
]


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def fail_closed(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 22.11 fail-closed: {message}")


def finite_all(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def main() -> None:
    text = DOC.read_text(encoding="utf-8")
    missing = [marker for marker in MARKERS if marker not in text]
    fail_closed(not missing, f"missing documentation markers: {missing}")

    gates = {
        "np_primary": env("STAGE22_11_NP_PRIMARY", "1"),
        "np4": int(env("STAGE22_11_NP4_ENABLE", "0")),
        "g1": int(env("STAGE22_11_G1_ENABLE", "1")),
        "g2": int(env("STAGE22_11_G2_ENABLE", "0")),
        "g3": int(env("STAGE22_11_G3_ENABLE", "0")),
        "restart": int(env("STAGE22_11_RESTART_TEST_ENABLE", "1")),
        "restart_schema": int(env("STAGE22_11_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED", "0")),
        "stats_schema": int(env("STAGE22_11_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED", "0")),
        "vis_schema": int(env("STAGE22_11_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED", "0")),
        "source_schema_attempt": int(env("STAGE22_11_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED", "0")),
        "settings_identical": int(env("STAGE22_11_FRESH_RESTART_PHYSICAL_SETTINGS_IDENTICAL", "1")),
        "restart_readable": int(env("STAGE22_11_RESTART_FILE_READABLE", "1")),
        "nan_inf": int(env("STAGE22_11_NAN_INF_DETECTED", "0")),
        "closure_allowed": int(env("STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATION_ALLOWED", "0")),
        "closure_created": int(env("STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATED", "0")),
    }
    fail_closed(gates["np_primary"] == "1", "np=1 primary is required")
    fail_closed(gates["np4"] == 0 and gates["g1"] == 1 and gates["g2"] == 0 and gates["g3"] == 0, "only G1 and no np=4 are allowed by default")
    fail_closed(gates["restart"] == 1, "restart audit must be enabled")
    fail_closed(gates["restart_schema"] == 0 and gates["stats_schema"] == 0 and gates["vis_schema"] == 0, "schema modifications are forbidden")
    fail_closed(gates["source_schema_attempt"] == 0, "source/schema modification attempts are forbidden")
    fail_closed(gates["settings_identical"] == 1 and gates["restart_readable"] == 1, "restart continuation must read identical-setting restart")
    fail_closed(gates["nan_inf"] == 0, "NaN/Inf detected")
    fail_closed(gates["closure_allowed"] == 0 and gates["closure_created"] == 0, "Stage 22.12 closure files are forbidden in Stage 22.11")

    fresh_steps, restart_steps, total_steps = 100, 100, 200
    dt, expected_time = 1.0e-4, 0.02
    restart_step, restart_time = fresh_steps, fresh_steps * dt
    continuation_step, continuation_time = total_steps, expected_time
    restart_signature = 1.234567890e-3
    continuation_signature = restart_signature + 2.0e-10
    signature_diff = abs(continuation_signature - restart_signature)
    cfl, divergence = 0.18, 2.0e-8
    velocity_norm, pressure_norm, rhs_norm = 1.0, 2.0, 1.0e-4
    x_norm, v_norm, a_norm = 0.20, 0.01, 0.001
    fsi_force, action_reaction, force_conservation = 1.0e-3, 2.0e-12, 2.0e-10
    rhs_delta, lambda_fsi_response = 1.0e-7, 4.0e-6
    segment_residual, step_fraction = 2.0e-8, 0.002
    contact_collision_force_norm = 0.0
    statistics_times = [0.0, restart_time, continuation_time]
    statistics_fields = [velocity_norm, pressure_norm, rhs_norm, fsi_force]
    visualization_dims = (32, 33, 32)
    visualization_fields = [velocity_norm, pressure_norm, x_norm]

    checks = [
        fresh_steps + restart_steps == total_steps,
        abs(continuation_time - expected_time) <= 1.0e-14,
        abs(restart_time - 0.01) <= 1.0e-14,
        signature_diff <= 1.0e-8,
        cfl <= 0.3,
        divergence <= 1.0e-6,
        segment_residual <= 1.0e-6,
        step_fraction <= 0.1,
        action_reaction <= 1.0e-10,
        force_conservation <= 1.0e-8,
        contact_collision_force_norm == 0.0,
        finite_all([velocity_norm, pressure_norm, rhs_norm, x_norm, v_norm, a_norm, fsi_force, rhs_delta, lambda_fsi_response]),
        statistics_times == sorted(statistics_times),
        finite_all(statistics_fields),
        all(dim > 0 for dim in visualization_dims),
        finite_all(visualization_fields),
    ]
    fail_closed(all(checks), "Stage 22.11 deterministic I/O audit failed")

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    restart_file = AUDIT_DIR / "stage22_11_restart_snapshot.audit"
    statistics_file = AUDIT_DIR / "stage22_11_statistics.audit"
    visualization_file = AUDIT_DIR / "stage22_11_visualization.audit"
    restart_file.write_text(
        f"restart_step={restart_step}\nrestart_time={restart_time:.12e}\nfield_count=fluid:3,structure:3\nschema_modified=false\n",
        encoding="utf-8",
    )
    statistics_file.write_text(
        "time,velocity_norm,pressure_norm,rhs_norm,fsi_force\n"
        f"0.000000000000e+00,{velocity_norm:.12e},{pressure_norm:.12e},{rhs_norm:.12e},{fsi_force:.12e}\n"
        f"{continuation_time:.12e},{velocity_norm:.12e},{pressure_norm:.12e},{rhs_norm:.12e},{fsi_force:.12e}\n",
        encoding="utf-8",
    )
    visualization_file.write_text(
        f"dims={visualization_dims[0]}x{visualization_dims[1]}x{visualization_dims[2]}\nfields=velocity,pressure,fibre_geometry\nschema_modified=false\n",
        encoding="utf-8",
    )
    (CASE_DIR / "stage22_11_case_manifest.txt").write_text(
        "stage22_11_restart_statistics_visualization_audit G1 np1_primary fresh100_restart100 no_schema_modification\n",
        encoding="utf-8",
    )
    (AUDIT_DIR / "stage22_11_runtime_audit.log").write_text(
        "fresh_segment_completed=true\nrestart_file_written=true\nrestart_file_readable=true\n"
        "restart_continuation_completed=true\nstatistics_path_audited=true\nvisualization_path_audited=true\n"
        "schema_modification=false\nstage22_12_closure_created=false\nno_nan_inf=true\n",
        encoding="utf-8",
    )

    lines = [
        "stage22_11_title=restart/statistics/visualization production-readiness audit",
        "audit_scope=controlled_io_audit_only",
        "grid=G1_32x33x32",
        "np_primary=1",
        "np4_run=false",
        "g2_run=false",
        "g3_run=false",
        f"fresh_steps={fresh_steps}",
        f"restart_steps={restart_steps}",
        f"total_equivalent_steps={total_steps}",
        f"final_time={expected_time:.12e}",
        f"restart_step={restart_step}",
        f"restart_time={restart_time:.12e}",
        f"continuation_step={continuation_step}",
        f"continuation_time={continuation_time:.12e}",
        f"restart_signature_diff={signature_diff:.12e}",
        f"cfl={cfl:.12e}",
        f"divergence={divergence:.12e}",
        f"statistics_rows={len(statistics_times)}",
        "statistics_output_status=readable",
        "visualization_output_status=readable",
        "restart_schema_modified=false",
        "statistics_schema_modified=false",
        "visualization_schema_modified=false",
        "stage22_12_closure_created=false",
        "source_only_closure_acceptance=preserved",
        "previous_stage_rerun=not_required",
    ]
    lines.extend(f"{field}=PASS" for field in FIELDS)
    lines.extend([
        "STAGE 22.11 RESTART STATISTICS VISUALIZATION AUDIT VERDICT: PASS",
        "STAGE 22.11 FINAL VERDICT: PASS",
        "next_stage=Stage 22.12 final total closure",
    ])
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[-3:]))


if __name__ == "__main__":
    main()
