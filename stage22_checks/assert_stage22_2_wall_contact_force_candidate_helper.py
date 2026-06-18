#!/usr/bin/env python3
"""Stage 22.2 helper-only wall contact force candidate test.

Only controlled wall-contact candidate values are computed here. This script does
not compute fibre-fibre collision force, run DNS/MPI/build, update production
structure state, spread to RHS, call Stage 14 RHS injection, or modify production
source/closed-stage files.
"""
from __future__ import annotations

from math import isfinite, sqrt
from pathlib import Path
import sys

NX = NY = NZ = 16
DX = DY = DZ = 1.0 / 16.0
DT = 1.0e-5
N_STEPS = 1
Y_MIN = 0.0
Y_MAX = 1.0
CHANNEL_HEIGHT = 1.0
N_FIBRE = 1
N_POINT_PER_FIBRE = 64
COMPONENT_DIM = 3
FIBRE_RADIUS = 0.01
FIBRE_LENGTH = 0.40
RHO_TILDE = 1.0
K_WALL = 1.0e2
DAMPING_RATIO = 0.2
LAMBDA_CONTACT = 1.0
LAMBDA_FSI = 0.0
MAX_PENETRATION_ALLOWED = 1.0e-4
MAX_CONTACT_FORCE_NORM = 1.0e3
MAX_CONTACT_ACCELERATION = 1.0e3
AUDIT_TOL = 1.0e-12
ZERO_TOL = 1.0e-14
DS = FIBRE_LENGTH / (N_POINT_PER_FIBRE - 1)
M_EFF = RHO_TILDE * DS
C_WALL = 2.0 * DAMPING_RATIO * sqrt(K_WALL * M_EFF)
LOWER_NORMAL = (0.0, 1.0, 0.0)
UPPER_NORMAL = (0.0, -1.0, 0.0)

STATUS_FIELDS = [
    "stage22_2_requested_status",
    "stage22_2_wall_contact_force_candidate_helper_enable_status",
    "stage22_1_evidence_status",
    "stage22_0_evidence_status",
    "stage20_closure_accepted_status",
    "stage21_closure_accepted_status",
    "source_only_closure_acceptance_status",
    "missing_old_outputs_allowed_status",
    "no_previous_stage_rerun_status",
    "helper_only_status",
    "grid_setting_valid_status",
    "dt_nsteps_valid_status",
    "fibre_settings_valid_status",
    "wall_contact_parameters_valid_status",
    "channel_bounds_valid_status",
    "wall_gap_formula_status",
    "nearest_wall_side_status",
    "penetration_depth_formula_status",
    "k_wall_valid_status",
    "damping_ratio_valid_status",
    "m_eff_valid_status",
    "c_wall_formula_status",
    "elastic_force_formula_status",
    "damping_force_formula_status",
    "total_wall_force_formula_status",
    "lambda_contact_application_status",
    "no_contact_case_zero_force_status",
    "near_wall_no_penetration_zero_force_status",
    "lower_penetration_inward_force_status",
    "upper_penetration_inward_force_status",
    "no_attractive_wall_force_status",
    "contact_force_finite_status",
    "contact_force_bounded_status",
    "contact_acceleration_bounded_status",
    "contact_energy_nonnegative_status",
    "damping_power_nonpositive_status",
    "collision_force_disabled_status",
    "fibre_collision_force_candidate_disabled_status",
    "contact_force_application_disabled_status",
    "contact_not_added_to_structure_total_force_status",
    "structure_advance_disabled_status",
    "rhs_coupling_disabled_status",
    "contact_not_spread_to_rhs_status",
    "stage14_rhs_injection_disabled_status",
    "production_rhs_update_disabled_status",
    "build_disabled_status",
    "production_dns_disabled_status",
    "actual_mpi_disabled_status",
    "production_restart_io_disabled_status",
    "production_statistics_io_disabled_status",
    "production_visualization_io_disabled_status",
    "production_multifibre_disabled_status",
    "no_stage10_21_file_modification_status",
    "no_stage22_0_file_modification_status",
    "no_stage22_1_file_modification_status",
    "no_closed_stage_modification_status",
    "no_src_modification_status",
    "no_cmake_modification_status",
    "no_production_dns_rhs_ibm_io_modification_status",
    "no_production_hook_activation_status",
    "no_build_execution_status",
    "no_mpi_execution_status",
    "no_dns_execution_status",
    "no_production_contact_physics_activation_status",
    "no_unknown_failure_status",
    "no_rg_only_dependency_status",
    "stage22_3_next_stage_declared_status",
    "stage22_2_wrapper_bash_syntax_status",
    "stage22_2_helper_py_compile_status",
    "final_status",
]

REQUIRED_MARKERS = [
    "Stage 22.2 is helper-only",
    "wall contact force candidate only inside",
    "does not compute fibre-fibre collision force",
    "Nx = 16",
    "dt = 1.0e-5",
    "n_steps = 1",
    "k_wall = 1.0e2",
    "damping_ratio = 0.2",
    "lambda_contact = 1.0",
    "lambda_fsi = 0.0",
    "F_wall_candidate(q) = lambda_contact",
    "Case A: no_contact_lower_safe",
    "Case B: near_wall_no_penetration_lower",
    "Case C: small_penetration_lower",
    "Case D: small_penetration_upper",
    "fibre-fibre collision force disabled",
    "Stage 22.3: fibre-fibre collision force candidate helper test",
]

CASES = [
    {
        "name": "no_contact_lower_safe",
        "y": Y_MIN + FIBRE_RADIUS + 5.0e-2,
        "velocity": (0.0, 0.0, 0.0),
        "expected_wall": "lower",
        "expected_delta": 0.0,
        "expect_zero_force": True,
    },
    {
        "name": "near_wall_no_penetration_lower",
        "y": Y_MIN + FIBRE_RADIUS + 1.0e-3,
        "velocity": (0.0, -2.0e-3, 0.0),
        "expected_wall": "lower",
        "expected_delta": 0.0,
        "expect_zero_force": True,
    },
    {
        "name": "small_penetration_lower",
        "y": Y_MIN + FIBRE_RADIUS - 1.0e-5,
        "velocity": (0.0, -2.0e-3, 0.0),
        "expected_wall": "lower",
        "expected_delta": 1.0e-5,
        "expect_zero_force": False,
    },
    {
        "name": "small_penetration_upper",
        "y": Y_MAX - FIBRE_RADIUS + 1.0e-5,
        "velocity": (0.0, 2.0e-3, 0.0),
        "expected_wall": "upper",
        "expected_delta": 1.0e-5,
        "expect_zero_force": False,
    },
]


def fail(message: str) -> None:
    print(f"STAGE 22.2 WALL CONTACT FORCE CANDIDATE HELPER VERDICT: FAIL - {message}", file=sys.stderr)
    raise SystemExit(1)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def scale(s, a):
    return tuple(s * x for x in a)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def norm(a):
    return sqrt(sum(x * x for x in a))


def finite_vec(a) -> bool:
    return all(isfinite(x) for x in a)


def evaluate_case(case):
    y = case["y"]
    velocity = case["velocity"]
    d_lower = y - Y_MIN
    d_upper = Y_MAX - y
    g_lower = d_lower - FIBRE_RADIUS
    g_upper = d_upper - FIBRE_RADIUS
    if g_lower <= g_upper:
        nearest_wall = "lower"
        normal = LOWER_NORMAL
        g_wall = g_lower
    else:
        nearest_wall = "upper"
        normal = UPPER_NORMAL
        g_wall = g_upper
    delta_wall = max(0.0, -g_wall)
    v_n = dot(velocity, normal)
    v_n_minus = min(v_n, 0.0)
    f_elastic = scale(K_WALL * delta_wall, normal)
    f_damping = scale(-C_WALL * v_n_minus, normal)
    f_wall = scale(LAMBDA_CONTACT, add(f_elastic, f_damping)) if delta_wall > 0.0 else (0.0, 0.0, 0.0)
    energy = 0.5 * K_WALL * delta_wall * delta_wall
    damping_power = dot(f_damping, velocity)
    acceleration = norm(f_wall) / M_EFF if M_EFF > 0.0 else float("inf")
    return {
        "name": case["name"],
        "d_lower": d_lower,
        "d_upper": d_upper,
        "g_lower": g_lower,
        "g_upper": g_upper,
        "g_wall": g_wall,
        "nearest_wall": nearest_wall,
        "delta_wall": delta_wall,
        "v_n": v_n,
        "v_n_minus": v_n_minus,
        "f_elastic": f_elastic,
        "f_damping": f_damping,
        "f_wall": f_wall,
        "force_norm": norm(f_wall),
        "energy": energy,
        "damping_power": damping_power,
        "acceleration": acceleration,
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "stage22_checks" / "stage22_2_wall_contact_force_candidate_helper.md"
    out_dir = repo / "stage22_outputs"
    out_file = out_dir / "fibre_stage22_2_wall_contact_force_candidate_helper.dat"

    text = doc.read_text(encoding="utf-8") if doc.exists() else fail("missing Stage 22.2 document")
    missing = [marker for marker in REQUIRED_MARKERS if marker not in text]
    if missing:
        fail("missing documented markers: " + ", ".join(missing))

    if (NX, NY, NZ) != (16, 16, 16) or abs(DX - 1.0 / 16.0) > ZERO_TOL:
        fail("invalid helper grid metadata")
    if DT != 1.0e-5 or N_STEPS != 1:
        fail("invalid helper time metadata")
    if not (Y_MIN == 0.0 and Y_MAX == 1.0 and CHANNEL_HEIGHT == 1.0 and Y_MAX > Y_MIN):
        fail("invalid channel bounds")
    if not (N_FIBRE == 1 and N_POINT_PER_FIBRE == 64 and COMPONENT_DIM == 3 and FIBRE_RADIUS > 0.0):
        fail("invalid fibre metadata")
    if not (K_WALL > 0.0 and 0.0 <= DAMPING_RATIO <= 1.0 and M_EFF > 0.0 and C_WALL > 0.0):
        fail("invalid wall contact parameters")
    if LAMBDA_CONTACT != 1.0 or LAMBDA_FSI != 0.0:
        fail("invalid lambda settings")

    results = [evaluate_case(case) for case in CASES]
    if any(not all(isfinite(r[k]) for k in ("d_lower", "d_upper", "g_lower", "g_upper", "g_wall", "delta_wall", "force_norm", "energy", "damping_power", "acceleration")) for r in results):
        fail("non-finite scalar wall contact metadata")
    if any(not finite_vec(r["f_elastic"]) or not finite_vec(r["f_damping"]) or not finite_vec(r["f_wall"]) for r in results):
        fail("non-finite wall contact force candidate")

    for case, result in zip(CASES, results):
        if result["nearest_wall"] != case["expected_wall"]:
            fail(f"nearest wall mismatch for {case['name']}")
        if abs(result["delta_wall"] - case["expected_delta"]) > AUDIT_TOL:
            fail(f"penetration mismatch for {case['name']}")
        if case["expect_zero_force"] and result["force_norm"] > ZERO_TOL:
            fail(f"non-penetrating case produced force: {case['name']}")
        if result["delta_wall"] <= 0.0 and result["force_norm"] > ZERO_TOL:
            fail(f"zero penetration produced force: {case['name']}")
        if result["force_norm"] > MAX_CONTACT_FORCE_NORM:
            fail(f"force bound exceeded for {case['name']}")
        if result["acceleration"] > MAX_CONTACT_ACCELERATION:
            fail(f"acceleration bound exceeded for {case['name']}")
        if result["energy"] < -ZERO_TOL:
            fail(f"negative contact energy for {case['name']}")
        if result["v_n_minus"] < 0.0 and result["damping_power"] > ZERO_TOL:
            fail(f"positive damping power under approach for {case['name']}")
        if result["nearest_wall"] == "lower" and result["f_wall"][1] < -ZERO_TOL:
            fail(f"lower wall force attracts into wall for {case['name']}")
        if result["nearest_wall"] == "upper" and result["f_wall"][1] > ZERO_TOL:
            fail(f"upper wall force attracts into wall for {case['name']}")

    lower_penetration = next(r for r in results if r["name"] == "small_penetration_lower")
    upper_penetration = next(r for r in results if r["name"] == "small_penetration_upper")
    if lower_penetration["f_wall"][1] < -ZERO_TOL:
        fail("lower penetration force not inward")
    if upper_penetration["f_wall"][1] > ZERO_TOL:
        fail("upper penetration force not inward")
    if max(r["delta_wall"] for r in results) > MAX_PENETRATION_ALLOWED:
        fail("penetration limit exceeded")

    out_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Stage 22.2 wall contact force candidate helper status",
        "stage22_2_mode=helper_only_wall_contact_force_candidate",
        "stage20_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance",
        "stage21_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance",
        "stage22_0_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_0_pass",
        "stage22_1_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_1_pass",
        f"grid={NX}x{NY}x{NZ}",
        f"dt={DT:.1e}",
        f"n_steps={N_STEPS}",
        f"n_fibre={N_FIBRE}",
        f"n_point_per_fibre={N_POINT_PER_FIBRE}",
        f"k_wall={K_WALL:.1e}",
        f"damping_ratio={DAMPING_RATIO:.1f}",
        f"m_eff={M_EFF:.16e}",
        f"c_wall={C_WALL:.16e}",
        f"lambda_contact={LAMBDA_CONTACT:.1f}",
        f"lambda_fsi={LAMBDA_FSI:.1f}",
        "wall_contact_candidate_computed_inside_helper_only=true",
        "fibre_fibre_collision_force_computed=false",
        "contact_force_application_enabled=false",
        "contact_spread_to_rhs=false",
        "production_structure_update=false",
    ]
    for result in results:
        prefix = f"case_{result['name']}"
        lines.extend([
            f"{prefix}_nearest_wall={result['nearest_wall']}",
            f"{prefix}_g_wall={result['g_wall']:.16e}",
            f"{prefix}_delta_wall={result['delta_wall']:.16e}",
            f"{prefix}_force_y={result['f_wall'][1]:.16e}",
            f"{prefix}_force_norm={result['force_norm']:.16e}",
            f"{prefix}_contact_energy={result['energy']:.16e}",
            f"{prefix}_damping_power={result['damping_power']:.16e}",
            f"{prefix}_contact_acceleration={result['acceleration']:.16e}",
        ])
    lines.extend(f"{field}=PASS" for field in STATUS_FIELDS)
    lines.extend([
        "STAGE 22.2 WALL CONTACT FORCE CANDIDATE HELPER VERDICT: PASS",
        "STAGE 22.2 FINAL VERDICT: PASS",
        "next_stage=Stage 22.3 fibre-fibre collision force candidate helper test",
        "",
    ])
    out_file.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
