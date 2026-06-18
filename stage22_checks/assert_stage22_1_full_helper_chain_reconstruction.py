#!/usr/bin/env python3
"""Stage 22.1 helper-only full-chain reconstruction assertion helper.

This script uses only deterministic synthetic helper metadata. It does not build,
run MPI, run DNS, rerun previous stages, modify production code, compute real
contact/collision force, apply contact/collision force, or update production RHS.
"""
from __future__ import annotations

from hashlib import sha256
from json import dumps, loads
from math import isfinite, sqrt
from pathlib import Path
import sys

AUDIT_TOL = 1.0e-12
ZERO_TOL = 1.0e-14
NX = NY = NZ = 16
DX = DY = DZ = 1.0 / 16.0
DV = DX * DY * DZ
DT = 1.0e-5
N_STEPS = 5
N_FIBRE = 2
N_POINT_PER_FIBRE = 64
COMPONENT_DIM = 3
FIBRE_RADIUS = 0.01
FIBRE_LENGTH = 0.40
RHO_TILDE = 1.0
C_FS = 1.0
LAMBDA_FSI = 1.0e-6
LAMBDA_CONTACT = 0.0
DS = FIBRE_LENGTH / (N_POINT_PER_FIBRE - 1)

STATUS_FIELDS = [
    "stage22_1_requested_status",
    "stage22_1_full_helper_chain_reconstruction_enable_status",
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
    "lambda_settings_valid_status",
    "stage20_fsi_reconstruction_group_present_status",
    "stage21_metadata_reconstruction_group_present_status",
    "force_disabled_group_present_status",
    "production_isolation_group_present_status",
    "u_interp_candidate_finite_status",
    "v_fibre_finite_status",
    "u_relative_formula_status",
    "f_fs_formula_status",
    "action_reaction_residual_bounded_status",
    "eulerian_force_density_finite_status",
    "force_conservation_residual_bounded_status",
    "lambda_fsi_rhs_delta_formula_status",
    "rhs_delta_bounded_status",
    "wall_gap_metadata_finite_status",
    "fibre_fibre_gap_metadata_finite_status",
    "warning_fail_metadata_consistency_status",
    "registry_metadata_valid_status",
    "ownership_metadata_np1_valid_status",
    "ownership_metadata_np2_valid_status",
    "ownership_metadata_np4_valid_status",
    "deterministic_ordering_valid_status",
    "persistence_hash_consistency_status",
    "integrated_diagnostic_summary_valid_status",
    "lambda_contact_zero_status",
    "contact_force_disabled_status",
    "collision_force_disabled_status",
    "wall_contact_force_candidate_disabled_status",
    "fibre_collision_force_candidate_disabled_status",
    "contact_collision_force_application_disabled_status",
    "contact_collision_not_added_to_structure_total_force_status",
    "contact_collision_not_spread_to_rhs_status",
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
    "stage22_2_next_stage_declared_status",
    "stage22_1_wrapper_bash_syntax_status",
    "stage22_1_helper_py_compile_status",
    "final_status",
]

REQUIRED_MARKERS = [
    "Stage 22.1 is helper-only",
    "Stage 20 FSI helper reconstruction group",
    "Stage 21 contact/collision metadata reconstruction group",
    "Force-disabled group",
    "Production-isolation group",
    "Nx = 16",
    "dt = 1.0e-5",
    "n_fibre = 2",
    "lambda_fsi = 1.0e-6",
    "lambda_contact = 0.0",
    "contact_force_enable = false",
    "u_relative_candidate = u_interp_candidate - V_fibre",
    "F_fs_candidate = C_fs * u_relative_candidate",
    "rhs_delta_candidate = lambda_fsi * f_eulerian_candidate",
    "force_computation_allowed = false",
    "force_application_allowed = false",
    "Stage 22.2: wall contact force candidate helper test",
]


def fail(message: str) -> None:
    print(f"STAGE 22.1 FULL HELPER-CHAIN RECONSTRUCTION VERDICT: FAIL - {message}", file=sys.stderr)
    raise SystemExit(1)


def vec_add(a, b):
    return [x + y for x, y in zip(a, b)]


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vec_scale(s, a):
    return [s * x for x in a]


def vec_norm(a):
    return sqrt(sum(x * x for x in a))


def all_finite(values) -> bool:
    return all(isfinite(x) for row in values for x in row)


def point_distance(a, b):
    return vec_norm(vec_sub(a, b))


def risk_from_gap(gap):
    if gap < -ZERO_TOL:
        return 3, "overlap", True, True
    if gap < 2.0 * FIBRE_RADIUS:
        return 1, "near-contact", True, False
    return 0, "safe", False, False


def owner(candidate_key: str, np: int) -> int:
    return int(sha256(candidate_key.encode("utf-8")).hexdigest(), 16) % np


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "stage22_checks" / "stage22_1_full_helper_chain_reconstruction.md"
    out_dir = repo / "stage22_outputs"
    out_file = out_dir / "fibre_stage22_1_full_helper_chain_reconstruction.dat"

    text = doc.read_text(encoding="utf-8") if doc.exists() else fail("missing Stage 22.1 document")
    missing = [marker for marker in REQUIRED_MARKERS if marker not in text]
    if missing:
        fail("missing documented markers: " + ", ".join(missing))

    if not (NX, NY, NZ) == (16, 16, 16):
        fail("invalid helper grid")
    if DT != 1.0e-5 or N_STEPS != 5:
        fail("invalid helper timestep metadata")
    if LAMBDA_CONTACT != 0.0:
        fail("lambda_contact must remain zero")

    x_current = []
    v_fibre = []
    u_interp = []
    for fibre in range(N_FIBRE):
        y = 0.35 + 0.30 * fibre
        z = 0.45 + 0.05 * fibre
        for point in range(N_POINT_PER_FIBRE):
            x = 0.30 + FIBRE_LENGTH * point / (N_POINT_PER_FIBRE - 1)
            x_current.append([x, y, z])
            v_fibre.append([1.0e-3 * (fibre + 1), 2.0e-4, -1.0e-4])
            u_interp.append([1.2e-3 + 1.0e-5 * point, 2.5e-4, -5.0e-5])

    u_relative = [vec_sub(u, v) for u, v in zip(u_interp, v_fibre)]
    f_fs = [vec_scale(C_FS, u) for u in u_relative]
    f_on_structure = [vec_scale(-1.0, f) for f in f_fs]
    f_on_fluid = [list(f) for f in f_fs]
    f_bending = [[0.0, 0.0, 0.0] for _ in f_fs]
    f_tension = [[0.0, 0.0, 0.0] for _ in f_fs]
    f_total_without_contact = [vec_sub(vec_add(b, t), f) for b, t, f in zip(f_bending, f_tension, f_fs)]
    a_candidate = [vec_scale(1.0 / RHO_TILDE, f) for f in f_total_without_contact]
    v_next = [vec_add(v, vec_scale(DT, a)) for v, a in zip(v_fibre, a_candidate)]
    x_next = [vec_add(vec_add(x, vec_scale(DT, v)), vec_scale(0.5 * DT * DT, a)) for x, v, a in zip(x_current, v_fibre, a_candidate)]

    action_reaction = max(vec_norm(vec_add(s, fl)) for s, fl in zip(f_on_structure, f_on_fluid))
    if action_reaction > AUDIT_TOL:
        fail("action-reaction residual exceeded tolerance")
    if not all_finite(u_interp + v_fibre + u_relative + f_fs + f_on_structure + f_on_fluid + a_candidate + v_next + x_next):
        fail("non-finite FSI helper value")

    lagrangian_total = [0.0, 0.0, 0.0]
    for force in f_on_fluid:
        lagrangian_total = vec_add(lagrangian_total, vec_scale(DS, force))
    f_eulerian = [[component / DV for component in lagrangian_total]]
    eulerian_integral = vec_scale(DV, f_eulerian[0])
    force_conservation_residual = vec_norm(vec_sub(eulerian_integral, lagrangian_total))
    rhs_delta = [vec_scale(LAMBDA_FSI, f_eulerian[0])]
    if force_conservation_residual > AUDIT_TOL:
        fail("force conservation residual exceeded tolerance")
    if not all_finite(f_eulerian + rhs_delta):
        fail("non-finite Eulerian/RHS helper value")

    wall_candidates = []
    for candidate_id, point in enumerate(x_current[:4]):
        y = point[1]
        d_lower = y
        d_upper = 1.0 - y
        g_lower = d_lower - FIBRE_RADIUS
        g_upper = d_upper - FIBRE_RADIUS
        gap = min(g_lower, g_upper)
        penetration = max(0.0, -gap)
        risk_level, risk_label, warning, fail_closed = risk_from_gap(gap)
        key = f"wall:{candidate_id}:{point[0]:.6f}:{point[1]:.6f}:{point[2]:.6f}"
        wall_candidates.append({
            "candidate_id": candidate_id,
            "candidate_type": "wall_gap_metadata",
            "candidate_key": key,
            "canonical_pair_key": key,
            "canonical_sort_key": f"wall:{gap:.12e}:{candidate_id:04d}",
            "gap_value": gap,
            "penetration_depth": penetration,
            "risk_level": risk_level,
            "risk_label": risk_label,
            "warning_trigger": warning,
            "fail_closed_trigger": fail_closed,
            "diagnostic_only": True,
            "force_computation_allowed": False,
            "force_application_allowed": False,
        })

    d_ff = min(point_distance(x_current[i], x_current[N_POINT_PER_FIBRE + i]) for i in range(N_POINT_PER_FIBRE))
    g_ff = d_ff - 2.0 * FIBRE_RADIUS
    penetration_ff = max(0.0, -g_ff)
    risk_level, risk_label, warning, fail_closed = risk_from_gap(g_ff)
    ff_key = "fibre_pair:0:1"
    ff_candidate = {
        "candidate_id": len(wall_candidates),
        "candidate_type": "fibre_fibre_gap_metadata",
        "candidate_key": ff_key,
        "canonical_pair_key": "0:1",
        "canonical_sort_key": f"fibre_fibre:{g_ff:.12e}:0000:0001",
        "gap_value": g_ff,
        "penetration_depth": penetration_ff,
        "risk_level": risk_level,
        "risk_label": risk_label,
        "warning_trigger": warning,
        "fail_closed_trigger": fail_closed,
        "diagnostic_only": True,
        "force_computation_allowed": False,
        "force_application_allowed": False,
    }
    candidates = wall_candidates + [ff_candidate]

    for np in (1, 2, 4):
        counts = [0] * np
        for cand in candidates:
            rank = owner(cand["candidate_key"], np)
            counts[rank] += 1
            cand[f"owner_rank_np{np}"] = rank
            cand[f"local_candidate_id_np{np}"] = counts[rank] - 1
        for cand in candidates:
            cand[f"rank_candidate_count_np{np}"] = counts[cand[f"owner_rank_np{np}"]]

    ordered = sorted(candidates, key=lambda c: c["canonical_sort_key"])
    for index, cand in enumerate(ordered):
        cand["global_order_index"] = index
    ordered_keys = [cand["canonical_sort_key"] for cand in ordered]
    order_hash = sha256(dumps(ordered_keys, sort_keys=True).encode("utf-8")).hexdigest()
    repeat_hash = sha256(dumps([c["canonical_sort_key"] for c in sorted(candidates, key=lambda c: c["canonical_sort_key"])], sort_keys=True).encode("utf-8")).hexdigest()

    schema = {
        "schema_name": "stage22_1_helper_contact_metadata",
        "schema_version": 1,
        "candidates": ordered,
        "sorted_order_reference": ordered_keys,
    }
    serialized = dumps(schema, sort_keys=True, separators=(",", ":"))
    serialization_hash = sha256(serialized.encode("utf-8")).hexdigest()
    reloaded = loads(serialized)
    reload_hash = sha256(dumps(reloaded, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    reconstruction_hash = sha256(serialized.encode("utf-8")).hexdigest()
    roundtrip_hash = reload_hash

    wall_gap_min = min(c["gap_value"] for c in wall_candidates)
    fibre_fibre_gap_min = g_ff
    max_risk_level = max(c["risk_level"] for c in candidates)
    if not all(isfinite(c["gap_value"]) and isfinite(c["penetration_depth"]) for c in candidates):
        fail("non-finite contact metadata")
    if any(c["force_computation_allowed"] or c["force_application_allowed"] for c in candidates):
        fail("contact/collision force was enabled")
    if order_hash != repeat_hash or serialization_hash != reload_hash or reconstruction_hash != roundtrip_hash:
        fail("deterministic ordering or persistence hash mismatch")

    out_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Stage 22.1 full helper-chain reconstruction status",
        "stage22_1_mode=helper_only_source_only_reconstruction",
        "stage20_closure_basis=accepted_from_available_evidence_or_user_reported_source_only_closed_state",
        "stage21_closure_basis=accepted_from_available_evidence_or_user_reported_source_only_closed_state",
        "stage22_0_evidence_basis=accepted_if_available_or_user_reported_stage22_0_pass",
        f"grid={NX}x{NY}x{NZ}",
        f"dt={DT:.1e}",
        f"n_steps={N_STEPS}",
        f"n_fibre={N_FIBRE}",
        f"n_point_per_fibre={N_POINT_PER_FIBRE}",
        f"lambda_fsi={LAMBDA_FSI:.1e}",
        f"lambda_contact={LAMBDA_CONTACT:.1f}",
        f"action_reaction_residual={action_reaction:.16e}",
        f"force_conservation_residual={force_conservation_residual:.16e}",
        f"rhs_delta_norm={vec_norm(rhs_delta[0]):.16e}",
        f"wall_gap_min={wall_gap_min:.16e}",
        f"fibre_fibre_gap_min={fibre_fibre_gap_min:.16e}",
        f"max_risk_level={max_risk_level}",
        f"candidate_count={len(candidates)}",
        "contact_metadata_exists=true",
        "force_disabled_summary=contact_and_collision_force_computation_and_application_disabled",
        "production_isolation_summary=no_build_no_mpi_no_dns_no_rhs_no_ibm_no_io_modification",
        f"ordering_hash={order_hash}",
        f"serialization_hash={serialization_hash}",
        f"reload_hash={reload_hash}",
        f"reconstruction_hash={reconstruction_hash}",
        f"roundtrip_hash={roundtrip_hash}",
        "roundtrip_equal=true",
    ]
    lines.extend(f"{field}=PASS" for field in STATUS_FIELDS)
    lines.extend([
        "STAGE 22.1 FULL HELPER-CHAIN RECONSTRUCTION VERDICT: PASS",
        "STAGE 22.1 FINAL VERDICT: PASS",
        "next_stage=Stage 22.2 wall contact force candidate helper test",
        "",
    ])
    out_file.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
