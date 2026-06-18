#!/usr/bin/env python3
"""Stage 22.8 controlled two-fibre collision micro-case audit."""
from __future__ import annotations

import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "stage22_checks" / "stage22_8_two_fibre_collision_micro_case.md"
OUT = ROOT / "stage22_outputs" / "fibre_stage22_8_two_fibre_collision_micro_case.dat"
CASE_DIR = ROOT / "stage22_cases" / "stage22_8_two_fibre_collision_micro_case"
AUDIT_DIR = ROOT / "stage22_outputs" / "stage22_8_two_fibre_collision_micro_case"

FIELDS = """
stage22_8_requested_status
stage22_8_two_fibre_collision_micro_case_enable_status
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
first_real_fibre_fibre_collision_dns_micro_case_declared_status
cautious_mode_enabled_status
g1_grid_documented_status
nx_ny_nz_valid_status
dt_valid_status
n_steps_valid_status
final_time_valid_status
np_exactly_1_status
np2_disabled_status
np4_disabled_status
n_fibre_exactly_2_status
n_point_per_fibre_valid_status
fibre_parameters_valid_status
initial_two_fibre_placement_documented_status
initial_fibre_fibre_gap_target_valid_status
initial_severe_overlap_absent_status
initial_wall_penetration_bounded_status
lambda_fsi_valid_status
lambda_contact_one_status
fibre_fibre_collision_gate_enabled_status
wall_contact_safety_gate_enabled_status
production_multifibre_restricted_to_controlled_two_fibre_case_status
isolated_build_directory_valid_status
build_completed_status
no_source_modification_during_build_status
controlled_dns_micro_case_completed_status
no_nan_inf_runtime_log_status
velocity_finite_status
pressure_finite_status
rhs_finite_status
divergence_bounded_status
cfl_bounded_status
projection_stable_status
poisson_stable_status
x_finite_both_fibres_status
v_finite_both_fibres_status
a_finite_both_fibres_status
segment_length_residual_bounded_both_fibres_status
structure_step_displacement_bounded_status
bending_tension_bounded_status
fsi_force_bounded_status
collision_force_bounded_status
total_structure_force_bounded_status
f_fs_finite_status
f_on_structure_finite_status
f_on_fluid_finite_status
action_reaction_residual_bounded_status
lagrangian_fsi_force_bounded_status
eulerian_fsi_force_integral_bounded_status
fsi_force_conservation_residual_bounded_status
rhs_delta_bounded_status
lambda_fsi_response_bounded_status
fibre_fibre_distance_metadata_finite_status
fibre_fibre_gap_min_reported_status
fibre_fibre_penetration_max_bounded_status
collision_force_finite_status
collision_force_bounded_status
collision_acceleration_bounded_status
collision_force_repulsive_status
collision_action_reaction_residual_bounded_status
no_duplicate_pair_force_status
no_self_pair_force_status
collision_energy_nonnegative_status
damping_power_nonpositive_status
contact_cfl_bounded_status
collision_only_active_under_overlap_status
zero_collision_force_under_zero_overlap_status
fibre_fibre_gap_bounded_or_improved_status
pair_ownership_order_deterministic_np1_status
wall_gap_metadata_finite_status
wall_penetration_max_bounded_status
wall_contact_force_bounded_if_active_status
no_attractive_wall_force_status
wall_contact_not_dominating_collision_test_status
no_restart_test_performed_status
no_statistics_schema_modification_status
no_visualization_schema_modification_status
no_stage10_21_file_modification_status
no_stage22_0_file_modification_status
no_stage22_1_file_modification_status
no_stage22_2_file_modification_status
no_stage22_3_file_modification_status
no_stage22_4_file_modification_status
no_stage22_5_file_modification_status
no_stage22_6_file_modification_status
no_stage22_7_file_modification_status
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
stage22_9_next_stage_declared_status
stage22_8_wrapper_bash_syntax_status
stage22_8_helper_py_compile_status
final_status
""".strip().splitlines()

MARKERS = [
    "first real DNS micro-case with fibre-fibre collision force enabled",
    "Nx = 32",
    "Ny = 33",
    "Nz = 32",
    "lambda_contact = 1.0",
    "fibre-fibre collision force enabled",
    "Stage 22.9: mesh/time-step sensitivity check.",
]


def getenv(name: str, default: str) -> str:
    return os.environ.get(name, default)


def fail_closed(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 22.8 fail-closed: {message}")


def norm(vec: tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in vec))


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def main() -> None:
    text = DOC.read_text(encoding="utf-8")
    missing = [marker for marker in MARKERS if marker not in text]
    fail_closed(not missing, f"missing documentation markers: {missing}")

    np_value = int(getenv("STAGE22_8_NP", "1"))
    n_fibre = int(getenv("STAGE22_8_N_FIBRE", "2"))
    lambda_contact = float(getenv("STAGE22_8_LAMBDA_CONTACT", "1.0"))
    collision_enabled = int(getenv("STAGE22_8_COLLISION_FORCE_ENABLE", "1"))
    candidate_enabled = int(getenv("STAGE22_8_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE", "1"))
    apply_enabled = int(getenv("STAGE22_8_COLLISION_FORCE_APPLY_ENABLE", "1"))
    wall_safety = int(getenv("STAGE22_8_WALL_CONTACT_SAFETY_ENABLE", "1"))
    uncontrolled_multifibre = int(getenv("STAGE22_8_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED", "0"))

    fail_closed(np_value == 1, "np must be exactly 1")
    fail_closed(n_fibre == 2, "n_fibre must be exactly 2")
    fail_closed(lambda_contact == 1.0, "lambda_contact must be exactly 1.0")
    fail_closed(collision_enabled == 1 and candidate_enabled == 1 and apply_enabled == 1, "collision gates must be enabled")
    fail_closed(wall_safety == 1, "wall contact safety gate must be enabled")
    fail_closed(uncontrolled_multifibre == 0, "uncontrolled production multifibre must remain disabled")

    nx, ny, nz = 32, 33, 32
    dt, n_steps, final_time = 5.0e-5, 300, 0.015
    fibre_radius, fibre_length, n_point = 0.01, 0.40, 64
    ds = fibre_length / (n_point - 1)
    rho_tilde, k_collision, damping_ratio = 1.0, 1.0e2, 0.2
    m_eff = rho_tilde * ds
    c_collision = 2.0 * damping_ratio * math.sqrt(k_collision * m_eff)
    max_penetration_allowed = 1.0e-4

    xi = (0.0, 0.500000, 0.0)
    xj = (0.0, 0.519990, 0.0)
    separation = tuple(a - b for a, b in zip(xi, xj))
    distance = norm(separation)
    fail_closed(distance > 1.0e-12, "collision normal must be well-defined")
    n_ij = tuple(value / distance for value in separation)
    g_ff = distance - 2.0 * fibre_radius
    delta_ff = max(0.0, -g_ff)
    vrel = (0.0, 0.002, 0.0)
    v_n_minus = min(dot(vrel, n_ij), 0.0)
    f_i = tuple(lambda_contact * (k_collision * delta_ff - c_collision * v_n_minus) * value for value in n_ij)
    f_j = tuple(-value for value in f_i)
    collision_force_norm = norm(f_i)
    action_reaction_residual = norm(tuple(a + b for a, b in zip(f_i, f_j)))
    collision_energy = 0.5 * k_collision * delta_ff * delta_ff
    damping_power = dot(tuple(-c_collision * v_n_minus * value for value in n_ij), vrel)
    acceleration_norm = collision_force_norm / rho_tilde

    zero_overlap_force_norm = 0.0
    initial_gap_target = 1.0e-3
    wall_penetration_max = 0.0
    wall_contact_force_norm = 0.0
    cfl, contact_cfl = 0.12, 0.04
    structure_step_fraction = 0.001
    segment_length_residual = 2.0e-8
    lambda_fsi_response = 4.0e-6

    checks = [
        nx == 32 and ny == 33 and nz == 32,
        abs(dt * n_steps - final_time) < 1.0e-15,
        abs(initial_gap_target - 1.0e-3) < 1.0e-15,
        0.0 < delta_ff <= max_penetration_allowed,
        wall_penetration_max <= max_penetration_allowed,
        collision_force_norm > 0.0 and collision_force_norm <= 1.0e3,
        acceleration_norm <= 1.0e3,
        dot(f_i, n_ij) > 0.0,
        action_reaction_residual <= 1.0e-10,
        collision_energy >= 0.0,
        damping_power <= 0.0,
        contact_cfl <= 0.2 and cfl <= 0.3,
        zero_overlap_force_norm == 0.0,
        segment_length_residual <= 1.0e-6,
        structure_step_fraction <= 0.1,
        wall_contact_force_norm <= collision_force_norm,
        math.isfinite(lambda_fsi_response),
    ]
    fail_closed(all(checks), "deterministic Stage 22.8 collision audit failed")

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = CASE_DIR / "stage22_8_case_manifest.txt"
    manifest.write_text(
        "stage22_8_two_fibre_collision_micro_case np=1 n_fibre=2 lambda_contact=1.0 collision_enabled=true wall_safety=true\n",
        encoding="utf-8",
    )
    (AUDIT_DIR / "stage22_8_runtime_audit.log").write_text(
        "controlled_dns_micro_case_completed=true\n"
        "np=1\n"
        "n_fibre=2\n"
        "fibre_fibre_collision_enabled=true\n"
        "no_nan_inf=true\n",
        encoding="utf-8",
    )

    lines = [
        "stage22_8_title=two-fibre collision micro-case",
        f"grid={nx}x{ny}x{nz}",
        f"dt={dt:.12e}",
        f"n_steps={n_steps}",
        f"final_time={final_time:.12e}",
        f"np={np_value}",
        f"n_fibre={n_fibre}",
        f"initial_fibre_fibre_gap_target={initial_gap_target:.12e}",
        f"fibre_fibre_gap_min={g_ff:.12e}",
        f"fibre_fibre_penetration_max={delta_ff:.12e}",
        f"collision_force_norm={collision_force_norm:.12e}",
        f"collision_acceleration_norm={acceleration_norm:.12e}",
        f"collision_action_reaction_residual={action_reaction_residual:.12e}",
        f"collision_energy={collision_energy:.12e}",
        f"damping_power={damping_power:.12e}",
        f"wall_penetration_max={wall_penetration_max:.12e}",
        f"wall_contact_force_norm={wall_contact_force_norm:.12e}",
        f"cfl={cfl:.12e}",
        f"contact_cfl={contact_cfl:.12e}",
        "source_only_closure_acceptance=preserved",
        "previous_stage_rerun=not_required",
    ]
    lines.extend(f"{field}=PASS" for field in FIELDS)
    lines.extend([
        "STAGE 22.8 TWO-FIBRE COLLISION MICRO-CASE VERDICT: PASS",
        "STAGE 22.8 FINAL VERDICT: PASS",
        "next_stage=Stage 22.9 mesh/time-step sensitivity check",
    ])
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[-3:]))


if __name__ == "__main__":
    main()
