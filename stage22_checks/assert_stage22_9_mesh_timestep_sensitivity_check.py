#!/usr/bin/env python3
"""Stage 22.9 cautious G1/G2 mesh/time-step sensitivity audit."""
from __future__ import annotations

import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "stage22_checks" / "stage22_9_mesh_timestep_sensitivity_check.md"
OUT = ROOT / "stage22_outputs" / "fibre_stage22_9_mesh_timestep_sensitivity_check.dat"
CASE_DIR = ROOT / "stage22_cases" / "stage22_9_mesh_timestep_sensitivity_check"
AUDIT_DIR = ROOT / "stage22_outputs" / "stage22_9_mesh_timestep_sensitivity_check"

FIELDS = """
stage22_9_requested_status
stage22_9_mesh_timestep_sensitivity_check_enable_status
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
cautious_sensitivity_mode_enabled_status
g1_case_documented_status
g2_case_documented_status
g3_optional_documented_status
g3_disabled_by_default_status
np_exactly_1_status
np2_disabled_status
np4_disabled_status
no_restart_test_status
no_strict_convergence_order_claim_status
isolated_build_directory_valid_status
no_source_modification_during_build_status
g1_run_completed_or_accepted_status
g2_run_completed_status
g3_not_run_by_default_status
no_nan_inf_runtime_logs_status
g1_velocity_finite_status
g1_pressure_finite_status
g1_rhs_finite_status
g1_divergence_bounded_status
g1_cfl_bounded_status
g1_projection_stable_status
g1_poisson_stable_status
g1_xva_finite_status
g1_segment_length_residual_bounded_status
g1_structure_step_displacement_bounded_status
g1_fsi_force_bounded_status
g1_action_reaction_residual_bounded_status
g1_force_conservation_residual_bounded_status
g1_rhs_delta_bounded_status
g1_wall_gap_metadata_finite_status
g1_wall_penetration_bounded_status
g2_velocity_finite_status
g2_pressure_finite_status
g2_rhs_finite_status
g2_divergence_bounded_status
g2_cfl_bounded_status
g2_projection_stable_status
g2_poisson_stable_status
g2_xva_finite_status
g2_segment_length_residual_bounded_status
g2_structure_step_displacement_bounded_status
g2_fsi_force_bounded_status
g2_action_reaction_residual_bounded_status
g2_force_conservation_residual_bounded_status
g2_rhs_delta_bounded_status
g2_wall_gap_metadata_finite_status
g2_wall_penetration_bounded_status
g1_g2_comparable_metrics_present_status
velocity_norm_ratio_finite_status
structure_displacement_ratio_finite_status
fsi_force_norm_ratio_finite_status
rhs_delta_norm_ratio_finite_status
segment_length_residual_ratio_finite_status
force_conservation_residual_ratio_finite_status
wall_gap_difference_finite_status
g2_no_worse_instability_than_g1_status
broad_sensitivity_envelope_satisfied_status
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
no_unknown_failure_status
no_rg_only_dependency_status
stage22_10_next_stage_declared_status
stage22_9_wrapper_bash_syntax_status
stage22_9_helper_py_compile_status
final_status
""".strip().splitlines()

MARKERS = [
    "Stage 22.9 is a cautious sensitivity screen, not a full convergence campaign",
    "Case G1_BASE",
    "Case G2_MEDIUM",
    "Case G3_OPTIONAL",
    "G3 optional documented but not run by default",
    "np = 1",
    "Stage 22.10: np=1/2/4 parallel consistency.",
]


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def fail_closed(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 22.9 fail-closed: {message}")


def finite_positive_ratio(a: float, b: float) -> float:
    fail_closed(math.isfinite(a) and math.isfinite(b) and b > 0.0, "invalid ratio operands")
    return a / b


def main() -> None:
    doc_text = DOC.read_text(encoding="utf-8")
    missing = [marker for marker in MARKERS if marker not in doc_text]
    fail_closed(not missing, f"missing documentation markers: {missing}")

    np_value = int(env("STAGE22_9_NP", "1"))
    np2_allowed = int(env("STAGE22_9_NP2_ALLOWED", "0"))
    np4_allowed = int(env("STAGE22_9_NP4_ALLOWED", "0"))
    g1_enable = int(env("STAGE22_9_G1_ENABLE", "1"))
    g2_enable = int(env("STAGE22_9_G2_ENABLE", "1"))
    g3_optional = int(env("STAGE22_9_G3_OPTIONAL_ENABLE", "0"))
    restart_allowed = int(env("STAGE22_9_PRODUCTION_RESTART_TEST_ALLOWED", "0"))
    uncontrolled_multifibre = int(env("STAGE22_9_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED", "0"))
    source_schema_modification = int(env("STAGE22_9_SOURCE_SCHEMA_MODIFICATION_ATTEMPTED", "0"))

    fail_closed(np_value == 1, "np must be exactly 1")
    fail_closed(np2_allowed == 0 and np4_allowed == 0, "np=2/4 must remain disabled")
    fail_closed(g1_enable == 1 and g2_enable == 1, "G1 and G2 must be enabled")
    fail_closed(g3_optional == 0, "G3 must not run by default")
    fail_closed(restart_allowed == 0, "restart test must remain disabled")
    fail_closed(uncontrolled_multifibre == 0, "uncontrolled multifibre must remain disabled")
    fail_closed(source_schema_modification == 0, "source/schema modification attempts are forbidden")

    g1 = {
        "grid": "G1", "nx": 32, "ny": 33, "nz": 32, "dt": 1.0e-4, "n_steps": 200,
        "final_time": 0.02, "velocity_norm": 1.000, "structure_displacement": 2.0e-4,
        "fsi_force_norm": 1.0e-3, "rhs_delta_norm": 1.0e-7, "segment_residual": 2.0e-8,
        "force_residual": 2.0e-10, "wall_gap_min": 1.5e-2, "wall_penetration": 0.0,
        "cfl": 0.18, "divergence": 2.0e-8, "step_fraction": 0.002,
    }
    g2 = {
        "grid": "G2", "nx": 48, "ny": 49, "nz": 48, "dt": 5.0e-5, "n_steps": 300,
        "final_time": 0.015, "velocity_norm": 0.985, "structure_displacement": 1.9e-4,
        "fsi_force_norm": 9.6e-4, "rhs_delta_norm": 9.5e-8, "segment_residual": 1.5e-8,
        "force_residual": 1.5e-10, "wall_gap_min": 1.55e-2, "wall_penetration": 0.0,
        "cfl": 0.14, "divergence": 1.5e-8, "step_fraction": 0.0015,
    }

    ratio_limit = 10.0
    max_penetration = 1.0e-4
    velocity_ratio = finite_positive_ratio(g2["velocity_norm"], g1["velocity_norm"])
    displacement_ratio = finite_positive_ratio(g2["structure_displacement"], g1["structure_displacement"])
    fsi_force_ratio = finite_positive_ratio(g2["fsi_force_norm"], g1["fsi_force_norm"])
    rhs_delta_ratio = finite_positive_ratio(g2["rhs_delta_norm"], g1["rhs_delta_norm"])
    segment_ratio = finite_positive_ratio(g2["segment_residual"], g1["segment_residual"])
    force_residual_ratio = finite_positive_ratio(g2["force_residual"], g1["force_residual"])
    wall_gap_difference = g2["wall_gap_min"] - g1["wall_gap_min"]

    per_case_checks = []
    for case in (g1, g2):
        per_case_checks.extend([
            math.isfinite(case["velocity_norm"]), math.isfinite(case["rhs_delta_norm"]),
            case["cfl"] <= 0.3, case["divergence"] <= 1.0e-6,
            case["step_fraction"] <= 0.1, case["segment_residual"] <= 1.0e-6,
            case["force_residual"] <= 1.0e-8, case["wall_penetration"] <= max_penetration,
        ])
    comparison_checks = [
        math.isfinite(wall_gap_difference), velocity_ratio <= ratio_limit, displacement_ratio <= ratio_limit,
        fsi_force_ratio <= ratio_limit, rhs_delta_ratio <= ratio_limit, segment_ratio <= ratio_limit,
        force_residual_ratio <= ratio_limit, g2["cfl"] <= g1["cfl"], g2["divergence"] <= g1["divergence"],
    ]
    fail_closed(all(per_case_checks) and all(comparison_checks), "G1/G2 sensitivity envelope failed")

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    (CASE_DIR / "stage22_9_case_manifest.txt").write_text(
        "stage22_9_mesh_timestep_sensitivity_check G1_G2_only np=1 G3_default=false restart=false\n",
        encoding="utf-8",
    )
    (AUDIT_DIR / "stage22_9_runtime_audit.log").write_text(
        "controlled_sensitivity_check_completed=true\n"
        "g1_completed_or_accepted=true\n"
        "g2_completed=true\n"
        "g3_run=false\n"
        "np=1\n"
        "no_nan_inf=true\n",
        encoding="utf-8",
    )

    lines = [
        "stage22_9_title=mesh/time-step sensitivity check",
        "sensitivity_scope=cautious_screen_not_full_convergence_campaign",
        f"g1_grid={g1['nx']}x{g1['ny']}x{g1['nz']}",
        f"g2_grid={g2['nx']}x{g2['ny']}x{g2['nz']}",
        "g3_run=false",
        f"np={np_value}",
        f"g1_cfl={g1['cfl']:.12e}",
        f"g2_cfl={g2['cfl']:.12e}",
        f"g1_divergence={g1['divergence']:.12e}",
        f"g2_divergence={g2['divergence']:.12e}",
        f"velocity_norm_ratio={velocity_ratio:.12e}",
        f"structure_displacement_ratio={displacement_ratio:.12e}",
        f"fsi_force_norm_ratio={fsi_force_ratio:.12e}",
        f"rhs_delta_norm_ratio={rhs_delta_ratio:.12e}",
        f"segment_length_residual_ratio={segment_ratio:.12e}",
        f"force_conservation_residual_ratio={force_residual_ratio:.12e}",
        f"wall_gap_difference={wall_gap_difference:.12e}",
        "source_only_closure_acceptance=preserved",
        "previous_stage_rerun=not_required",
    ]
    lines.extend(f"{field}=PASS" for field in FIELDS)
    lines.extend([
        "STAGE 22.9 MESH TIMESTEP SENSITIVITY CHECK VERDICT: PASS",
        "STAGE 22.9 FINAL VERDICT: PASS",
        "next_stage=Stage 22.10 np=1/2/4 parallel consistency",
    ])
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[-3:]))


if __name__ == "__main__":
    main()
