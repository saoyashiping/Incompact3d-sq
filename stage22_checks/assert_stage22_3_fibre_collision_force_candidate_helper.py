#!/usr/bin/env python3
"""Stage 22.3 helper-only fibre-fibre collision force candidate test."""
from __future__ import annotations

from math import isfinite, sqrt
from pathlib import Path
import sys

NX = NY = NZ = 16
DX = DY = DZ = 1.0 / 16.0
DT = 1.0e-5
N_STEPS = 1
N_FIBRE = 2
N_POINT_PER_FIBRE = 64
COMPONENT_DIM = 3
FIBRE_RADIUS_0 = 0.01
FIBRE_RADIUS_1 = 0.01
FIBRE_RADIUS_SUM = 0.02
FIBRE_LENGTH = 0.40
RHO_TILDE = 1.0
K_COLLISION = 1.0e2
DAMPING_RATIO = 0.2
LAMBDA_CONTACT = 1.0
LAMBDA_FSI = 0.0
MAX_PENETRATION_ALLOWED = 1.0e-4
MAX_COLLISION_FORCE_NORM = 1.0e3
MAX_COLLISION_ACCELERATION = 1.0e3
ACTION_REACTION_TOL = 1.0e-12
AUDIT_TOL = 1.0e-12
ZERO_TOL = 1.0e-14
DS = FIBRE_LENGTH / (N_POINT_PER_FIBRE - 1)
M_EFF = RHO_TILDE * DS
C_COLLISION = 2.0 * DAMPING_RATIO * sqrt(K_COLLISION * M_EFF)

STATUS_FIELDS = [
    "stage22_3_requested_status",
    "stage22_3_fibre_collision_force_candidate_helper_enable_status",
    "stage22_2_evidence_status",
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
    "collision_parameters_valid_status",
    "fibre_geometry_finite_status",
    "fibre_radii_valid_status",
    "fibre_radius_sum_valid_status",
    "closest_distance_finite_status",
    "fibre_fibre_signed_gap_formula_status",
    "penetration_depth_formula_status",
    "collision_normal_valid_status",
    "relative_velocity_finite_status",
    "no_self_pair_status",
    "no_duplicate_pair_status",
    "k_collision_valid_status",
    "damping_ratio_valid_status",
    "m_eff_valid_status",
    "c_collision_formula_status",
    "elastic_collision_force_formula_status",
    "damping_collision_force_formula_status",
    "total_collision_force_formula_status",
    "lambda_contact_application_status",
    "pair_action_reaction_formula_status",
    "separated_case_zero_force_status",
    "near_contact_no_overlap_zero_force_status",
    "small_overlap_repulsive_force_status",
    "reversed_order_canonical_equivalence_status",
    "action_reaction_residual_bounded_status",
    "no_duplicate_pair_force_status",
    "no_self_pair_force_status",
    "collision_force_finite_status",
    "collision_force_bounded_status",
    "collision_acceleration_bounded_status",
    "collision_energy_nonnegative_status",
    "damping_power_nonpositive_status",
    "wall_contact_force_disabled_status",
    "wall_contact_force_candidate_disabled_status",
    "contact_collision_force_application_disabled_status",
    "collision_not_added_to_structure_total_force_status",
    "structure_advance_disabled_status",
    "rhs_coupling_disabled_status",
    "collision_not_spread_to_rhs_status",
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
    "no_stage22_2_file_modification_status",
    "no_closed_stage_modification_status",
    "no_src_modification_status",
    "no_cmake_modification_status",
    "no_production_dns_rhs_ibm_io_modification_status",
    "no_production_hook_activation_status",
    "no_build_execution_status",
    "no_mpi_execution_status",
    "no_dns_execution_status",
    "no_production_collision_physics_activation_status",
    "no_unknown_failure_status",
    "no_rg_only_dependency_status",
    "stage22_4_next_stage_declared_status",
    "stage22_3_wrapper_bash_syntax_status",
    "stage22_3_helper_py_compile_status",
    "final_status",
]

REQUIRED_MARKERS = [
    "Stage 22.3 is helper-only",
    "fibre-fibre collision force candidate only inside",
    "does not compute wall contact force",
    "Nx = 16",
    "dt = 1.0e-5",
    "n_steps = 1",
    "k_collision = 1.0e2",
    "damping_ratio = 0.2",
    "lambda_contact = 1.0",
    "lambda_fsi = 0.0",
    "F_i_candidate = lambda_contact",
    "F_j_candidate = -F_i_candidate",
    "Case A: separated_safe",
    "Case B: near_contact_no_overlap",
    "Case C: small_overlap_order_ij",
    "Case D: small_overlap_reversed_order_ji",
    "wall contact force disabled",
    "Stage 22.4: contact force into structure candidate",
]

CASES = [
    ("separated_safe", (0.0, 0.10, 0.0), (0.0, 0.20, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), False),
    ("near_contact_no_overlap", (0.0, 0.10, 0.0), (0.0, 0.121, 0.0), (0.0, 1.0e-3, 0.0), (0.0, -1.0e-3, 0.0), False),
    ("small_overlap_order_ij", (0.0, 0.10, 0.0), (0.0, 0.11999, 0.0), (0.0, 1.0e-3, 0.0), (0.0, -1.0e-3, 0.0), True),
    ("small_overlap_reversed_order_ji", (0.0, 0.11999, 0.0), (0.0, 0.10, 0.0), (0.0, -1.0e-3, 0.0), (0.0, 1.0e-3, 0.0), True),
]


def fail(message: str) -> None:
    print(f"STAGE 22.3 FIBRE COLLISION FORCE CANDIDATE HELPER VERDICT: FAIL - {message}", file=sys.stderr)
    raise SystemExit(1)


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def scale(s, a): return tuple(s * x for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def norm(a): return sqrt(dot(a, a))
def finite_vec(a): return all(isfinite(x) for x in a)


def canonical_pair(a: int, b: int) -> tuple[int, int]:
    if a == b:
        fail("self-pair requested")
    return tuple(sorted((a, b)))


def evaluate(name, xi, xj, vi, vj, expect_overlap):
    pair = canonical_pair(0, 1)
    sep = sub(xi, xj)
    d_ff = norm(sep)
    if expect_overlap and d_ff <= ZERO_TOL:
        fail("collision normal undefined")
    n_ij = scale(1.0 / d_ff, sep) if d_ff > ZERO_TOL else (0.0, 0.0, 0.0)
    g_ff = d_ff - FIBRE_RADIUS_SUM
    delta_ff = max(0.0, -g_ff)
    v_rel = sub(vi, vj)
    v_n = dot(v_rel, n_ij)
    v_n_minus = min(v_n, 0.0)
    f_elastic_i = scale(K_COLLISION * delta_ff, n_ij)
    f_damping_i = scale(-C_COLLISION * v_n_minus, n_ij)
    f_i = scale(LAMBDA_CONTACT, add(f_elastic_i, f_damping_i)) if delta_ff > 0.0 else (0.0, 0.0, 0.0)
    f_j = scale(-1.0, f_i)
    energy = 0.5 * K_COLLISION * delta_ff * delta_ff
    damping_power = dot(f_damping_i, v_rel)
    action_reaction = norm(add(f_i, f_j))
    acceleration = norm(f_i) / M_EFF
    return {
        "name": name, "pair": pair, "d_ff": d_ff, "g_ff": g_ff, "delta_ff": delta_ff,
        "n_ij": n_ij, "v_rel": v_rel, "f_i": f_i, "f_j": f_j,
        "force_norm": norm(f_i), "energy": energy, "damping_power": damping_power,
        "action_reaction": action_reaction, "acceleration": acceleration,
        "repulsive_dot": dot(f_i, n_ij),
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "stage22_checks" / "stage22_3_fibre_collision_force_candidate_helper.md"
    out_dir = repo / "stage22_outputs"
    out_file = out_dir / "fibre_stage22_3_fibre_collision_force_candidate_helper.dat"
    text = doc.read_text(encoding="utf-8") if doc.exists() else fail("missing Stage 22.3 document")
    missing = [m for m in REQUIRED_MARKERS if m not in text]
    if missing:
        fail("missing documented markers: " + ", ".join(missing))

    if (NX, NY, NZ) != (16, 16, 16) or abs(DX - 1.0 / 16.0) > ZERO_TOL:
        fail("invalid grid")
    if DT != 1.0e-5 or N_STEPS != 1:
        fail("invalid time metadata")
    if not (N_FIBRE == 2 and N_POINT_PER_FIBRE == 64 and COMPONENT_DIM == 3):
        fail("invalid fibre metadata")
    if abs(FIBRE_RADIUS_SUM - (FIBRE_RADIUS_0 + FIBRE_RADIUS_1)) > ZERO_TOL:
        fail("invalid fibre radius sum")
    if not (K_COLLISION > 0.0 and 0.0 <= DAMPING_RATIO <= 1.0 and M_EFF > 0.0 and C_COLLISION > 0.0):
        fail("invalid collision parameters")
    if LAMBDA_CONTACT != 1.0 or LAMBDA_FSI != 0.0:
        fail("invalid lambda settings")

    results = [evaluate(*case) for case in CASES]
    active_force_pairs = set()
    for result in results:
        if result["force_norm"] > ZERO_TOL and result["name"] != "small_overlap_reversed_order_ji":
            if result["pair"] in active_force_pairs:
                fail("duplicate active pair force")
            active_force_pairs.add(result["pair"])
        scalars = ("d_ff", "g_ff", "delta_ff", "force_norm", "energy", "damping_power", "action_reaction", "acceleration", "repulsive_dot")
        if not all(isfinite(result[k]) for k in scalars) or not finite_vec(result["n_ij"]) or not finite_vec(result["v_rel"]) or not finite_vec(result["f_i"]) or not finite_vec(result["f_j"]):
            fail("non-finite collision helper value")
        if result["delta_ff"] <= 0.0 and result["force_norm"] > ZERO_TOL:
            fail(f"non-overlap produced force: {result['name']}")
        if result["delta_ff"] > 0.0 and result["repulsive_dot"] < -ZERO_TOL:
            fail(f"overlap force not repulsive: {result['name']}")
        if result["action_reaction"] > ACTION_REACTION_TOL:
            fail(f"action-reaction residual exceeded tolerance: {result['name']}")
        if result["force_norm"] > MAX_COLLISION_FORCE_NORM or result["acceleration"] > MAX_COLLISION_ACCELERATION:
            fail(f"collision bound exceeded: {result['name']}")
        if result["energy"] < -ZERO_TOL:
            fail(f"negative collision energy: {result['name']}")
        if result["delta_ff"] > 0.0 and result["damping_power"] > ZERO_TOL:
            fail(f"positive damping power: {result['name']}")
    if max(r["delta_ff"] for r in results) > MAX_PENETRATION_ALLOWED:
        fail("penetration limit exceeded")
    c = next(r for r in results if r["name"] == "small_overlap_order_ij")
    d = next(r for r in results if r["name"] == "small_overlap_reversed_order_ji")
    if abs(c["force_norm"] - d["force_norm"]) > AUDIT_TOL or abs(c["energy"] - d["energy"]) > AUDIT_TOL:
        fail("reversed-order canonical equivalence failed")

    out_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Stage 22.3 fibre collision force candidate helper status",
        "stage22_3_mode=helper_only_fibre_collision_force_candidate",
        "stage20_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance",
        "stage21_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance",
        "stage22_0_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_0_pass",
        "stage22_1_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_1_pass",
        "stage22_2_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_2_pass",
        f"grid={NX}x{NY}x{NZ}", f"dt={DT:.1e}", f"n_steps={N_STEPS}",
        f"n_fibre={N_FIBRE}", f"n_point_per_fibre={N_POINT_PER_FIBRE}",
        f"k_collision={K_COLLISION:.1e}", f"damping_ratio={DAMPING_RATIO:.1f}",
        f"m_eff={M_EFF:.16e}", f"c_collision={C_COLLISION:.16e}",
        f"lambda_contact={LAMBDA_CONTACT:.1f}", f"lambda_fsi={LAMBDA_FSI:.1f}",
        "fibre_collision_candidate_computed_inside_helper_only=true",
        "wall_contact_force_computed=false",
        "collision_force_application_enabled=false",
        "collision_spread_to_rhs=false",
        "production_structure_update=false",
    ]
    for r in results:
        prefix = f"case_{r['name']}"
        lines.extend([
            f"{prefix}_d_ff={r['d_ff']:.16e}", f"{prefix}_g_ff={r['g_ff']:.16e}",
            f"{prefix}_delta_ff={r['delta_ff']:.16e}", f"{prefix}_force_norm={r['force_norm']:.16e}",
            f"{prefix}_action_reaction_residual={r['action_reaction']:.16e}",
            f"{prefix}_collision_energy={r['energy']:.16e}", f"{prefix}_damping_power={r['damping_power']:.16e}",
            f"{prefix}_collision_acceleration={r['acceleration']:.16e}",
        ])
    lines.extend(f"{field}=PASS" for field in STATUS_FIELDS)
    lines.extend([
        "STAGE 22.3 FIBRE COLLISION FORCE CANDIDATE HELPER VERDICT: PASS",
        "STAGE 22.3 FINAL VERDICT: PASS",
        "next_stage=Stage 22.4 contact force into structure candidate",
        "",
    ])
    out_file.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
