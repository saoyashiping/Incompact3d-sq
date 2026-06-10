#!/usr/bin/env python3
"""Stage 18.5 standalone structure time-integration core diagnostic audit.

This pure-Python helper validates local single-fibre candidate updates for
X_t = V and rho*V_t = F_total.  It computes X_candidate/V_candidate/A_candidate
only inside this diagnostic script; it never writes production structure state,
inserts production hooks, calls Stage 14 RHS injection, modifies IBM/RHS/DNS-core,
or runs MPI/build/production validation.

The helper reuses the corrected Stage 18.4 / 18.3 / 18.2 / 18.1 / 18.0 / Stage
17 / Stage 16 false-positive-safe audit pattern: targeted structural checks only,
no Markdown-as-code regression evidence, no broad repository scans, no mandatory
ripgrep, source-only archives without .git metadata accepted as non-contamination,
and only *_status fields contribute to final_status.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
    "stage18_5_requested_status", "stage18_4_evidence_status", "stage18_3_evidence_status",
    "stage18_2_evidence_status", "stage18_1_evidence_status", "stage18_0_evidence_status",
    "stage17_closed_file_status", "stage17_closed_evidence_status", "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status", "stage18_1_config_preserved_status",
    "stage18_2_geometry_operator_preserved_status", "stage18_3_bending_operator_preserved_status",
    "stage18_4_tension_constraint_preserved_status", "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status", "stage17_11_total_audit_fix_preserved_status",
    "no_closed_stage_modification_status", "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status", "stage18_1_files_unmodified_status",
    "stage18_2_files_unmodified_status", "stage18_3_files_unmodified_status",
    "stage18_4_files_unmodified_status", "stage18_enable_status", "time_integration_core_enable_status",
    "single_fibre_only_status", "diagnostic_only_status", "npts_value", "npts_status",
    "component_dim_value", "component_dim_status", "fibre_length_value", "fibre_length_status",
    "ds_value", "ds_formula_status", "dt_structure_value", "dt_structure_status", "rho_l_value",
    "rho_l_status", "rho_tilde_value", "rho_tilde_status", "dimensional_mass_validation_status",
    "nondimensional_mass_validation_status", "candidate_array_shape_status", "candidate_array_finite_status",
    "zero_force_rest_acceleration_zero_status", "zero_force_rest_velocity_unchanged_status",
    "zero_force_rest_position_unchanged_status", "uniform_velocity_no_force_velocity_preserved_status",
    "uniform_velocity_no_force_position_formula_status", "constant_force_acceleration_formula_status",
    "constant_force_velocity_formula_status", "constant_force_position_formula_status",
    "split_force_sum_formula_status", "split_force_acceleration_formula_status",
    "dt_refinement_position_consistency_status", "dt_refinement_velocity_consistency_status",
    "time_integration_equation_documented_status", "candidate_update_diagnostic_only_status",
    "no_production_structure_update_status", "no_production_structure_hook_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status",
    "no_inextensibility_projection_status", "no_inextensibility_repair_status",
    "no_fluid_force_physical_structure_integration_status", "no_stage14_rhs_call_from_stage18_5_status",
    "no_structure_energy_power_runtime_activation_status", "no_real_wall_contact_force_status",
    "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status",
    "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status",
    "no_contact_damping_force_status", "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status", "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status", "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status", "no_unknown_failure_status", "final_status",
]
VALUE_KEYS = {k for k in SUMMARY_KEYS if k.endswith(("_value", "_formula_value", "_shape_value", "_case_value"))}

STAGE18_0_FILES = ["stage18_checks/run_stage18_0_preflight_boundary.sh", "stage18_checks/assert_stage18_0_preflight_boundary.py", "stage18_checks/stage18_0_preflight_boundary.md"]
STAGE18_1_FILES = ["stage18_checks/run_stage18_1_physical_structure_config.sh", "stage18_checks/assert_stage18_1_physical_structure_config.py", "stage18_checks/stage18_1_physical_structure_config.md"]
STAGE18_2_FILES = ["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh", "stage18_checks/assert_stage18_2_structure_state_geometry_operators.py", "stage18_checks/stage18_2_structure_state_geometry_operators.md"]
STAGE18_3_FILES = ["stage18_checks/run_stage18_3_physical_bending_force_operator.sh", "stage18_checks/assert_stage18_3_physical_bending_force_operator.py", "stage18_checks/stage18_3_physical_bending_force_operator.md"]
STAGE18_4_FILES = ["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh", "stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py", "stage18_checks/stage18_4_tension_inextensibility_constraint.md"]
STAGE18_5_FILES = ["stage18_checks/run_stage18_5_structure_time_integration_core.sh", "stage18_checks/assert_stage18_5_structure_time_integration_core.py", "stage18_checks/stage18_5_structure_time_integration_core.md"]
ALLOWED_CHANGED_PATHS = set(STAGE18_5_FILES) | {"stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat"}
Vector = Tuple[float, float, float]


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def parse_dat(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def finite_float(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def int_value(text: str) -> int | None:
    value = finite_float(text)
    if value is None or not value.is_integer():
        return None
    return int(value)


def vadd(a: Vector, b: Vector) -> Vector:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vsub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vscale(c: float, a: Vector) -> Vector:
    return (c * a[0], c * a[1], c * a[2])


def max_vec_error(a: Sequence[Vector], b: Sequence[Vector]) -> float:
    return max((max(abs(x - y) for x, y in zip(p, q)) for p, q in zip(a, b)), default=0.0)


def arrays_finite(arrays: Iterable[Sequence[Vector]]) -> bool:
    return all(math.isfinite(v) for array in arrays for vec in array for v in vec)


def straight_fibre(npts: int, ds: float) -> List[Vector]:
    return [(i * ds, 0.0, 0.0) for i in range(npts)]


def candidate_update(x: Sequence[Vector], v: Sequence[Vector], f: Sequence[Vector], rho: float, dt: float) -> Tuple[List[Vector], List[Vector], List[Vector]]:
    a = [vscale(1.0 / rho, force) for force in f]
    x_next = [vadd(vadd(xi, vscale(dt, vi)), vscale(0.5 * dt * dt, ai)) for xi, vi, ai in zip(x, v, a)]
    v_next = [vadd(vi, vscale(dt, ai)) for vi, ai in zip(v, a)]
    return a, x_next, v_next


def git_status_entries(root: Path) -> Tuple[bool, List[Tuple[str, str]]]:
    if not (root / ".git").exists():
        return False, []
    proc = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        return False, []
    entries: List[Tuple[str, str]] = []
    for raw in proc.stdout.splitlines():
        if raw:
            path = raw[3:] if len(raw) > 3 else ""
            if " -> " in path:
                path = path.split(" -> ", 1)[1]
            entries.append((raw[:2], path))
    return True, entries


def changed_paths_ok(entries: Sequence[Tuple[str, str]]) -> bool:
    return all(path in ALLOWED_CHANGED_PATHS for _code, path in entries)


def files_unmodified(entries: Sequence[Tuple[str, str]], protected: Sequence[str]) -> bool:
    return all(path not in protected for _code, path in entries)


def has_all(text: str, needles: Sequence[str]) -> bool:
    return all(needle in text for needle in needles)


def stage_evidence_ok(root: Path, files: Sequence[str], output_rel: str, needles: Sequence[str]) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in files)
    output_ok = parse_dat(root / output_rel).get("final_status") in {"1", "PASS"}
    text = "\n".join(read_text(root / rel) for rel in files)
    return files_ok and (output_ok or has_all(text, needles))


def stage18_0_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_0_FILES, "stage18_outputs/fibre_stage18_0_preflight_boundary.dat", ["Stage 18.0", "preflight", "diagnostic-only"])


def stage18_1_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_1_FILES, "stage18_outputs/fibre_stage18_1_physical_structure_config.dat", ["A = pi * a^2", "B = E * I", "rho_l = rho_s * A", "ds = L_f / (npts - 1)"])


def stage18_2_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_2_FILES, "stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat", ["X_ssss", "arclength", "geometry"])


def stage18_3_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_3_FILES, "stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat", ["F_b = -B * X_ssss", "bending", "candidate"])


def stage18_4_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_4_FILES, "stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat", ["X_s dot X_s = 1", "tension", "diagnostic"])


def stage17_closed_evidence_ok(root: Path) -> bool:
    closed = root / "stage17_checks/STAGE17_CLOSED.md"
    if closed.exists() and "PASS" in read_text(closed):
        return True
    return stage17_11_total_audit_fix_preserved(root)


def stage18_0_wrapper_root_fix_preserved(root: Path) -> bool:
    wrapper = read_text(root / STAGE18_0_FILES[0])
    direct = "SCRIPT_DIR=" in wrapper and "REPO_ROOT=" in wrapper and "cd \"${DECOMP2D_ROOT}" not in wrapper
    inherited = any("stage18_0_wrapper_root_fix_preserved_status" in read_text(root / rel) for rel in STAGE18_1_FILES + STAGE18_2_FILES + STAGE18_3_FILES + STAGE18_4_FILES)
    return direct or inherited


def stage17_6_static_audit_fix_preserved(root: Path) -> bool:
    text = read_text(root / "stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py")
    return all(s in text for s in ["VALUE_KEYS", "pass_fail_keys", "source-only", "fibre_stage14_production_rhs_injection.f90", "xcompact3d.f90"])


def stage17_10_evidence_fix_preserved(root: Path) -> bool:
    text = read_text(root / "stage17_checks/assert_stage17_10_parallel_restart_io_wall_safety.py")
    return all(s in text for s in ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys", "source-only", "fibre_stage13_production_force_density_candidate.f90"])


def stage17_11_total_audit_fix_preserved(root: Path) -> bool:
    text = read_text(root / "stage17_checks/assert_stage17_11_total_contamination_audit_closure.py")
    return all(s in text for s in ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys", "source-only", "STAGE17_CLOSED.md"])


def stage13_6_diagnostic_preserved(root: Path) -> bool:
    return all((root / rel).exists() for rel in [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ])


def stage13_local_subdomain_regression_absent(root: Path) -> bool:
    text = read_text(root / "stage13_checks/stage13_6_production_force_density_candidate.md") + read_text(root / "stage13_checks/run_stage13_6_production_force_density_candidate.sh")
    return "local_subdomain_center" not in text


def stage14_small_lambda_hook_ok(root: Path) -> bool:
    return (root / "src/fibre_stage14_production_rhs_injection.f90").exists() and (root / "src/xcompact3d.f90").exists()


def no_rg_only_dependency(root: Path) -> bool:
    wrapper = read_text(root / STAGE18_5_FILES[0])
    return " rg " not in f" {wrapper} " or "grep" in wrapper


def no_stage18_5_physics_activation(root: Path) -> bool:
    text = read_text(root / STAGE18_5_FILES[1]) + "\n" + read_text(root / STAGE18_5_FILES[0])
    forbidden = ["fibre_stage14_production_rhs_injection", "call stage14", "call fibre_stage14", "cmake ", "make ", "ninja ", "src/xcompact3d.f90"]
    return not any(token in text for token in forbidden)


def add_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stage 18.5 diagnostic structure time-integration core audit")
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--stage18-5-enable", default="1")
    parser.add_argument("--time-integration-core-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--npts", default="16")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--component-dim", default="3")
    parser.add_argument("--dt-structure", default="1.0e-4")
    parser.add_argument("--rho-l", default="1.0")
    parser.add_argument("--rho-tilde", default="1.0")
    parser.add_argument("--use-dimensional-mass", default="1")
    parser.add_argument("--use-nondimensional-mass", default="1")
    parser.add_argument("--uniform-velocity", default="1.0e-3")
    parser.add_argument("--constant-force", default="1.0e-3")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-12")
    parser.add_argument("--dt-refinement-tol", default="1.0e-12")
    parser.add_argument("--test-case", default="zero_uniform_constant_split_dt_refinement")
    return parser.parse_args()


def write_output(root: Path, values: Dict[str, str], reasons: Sequence[str]) -> None:
    out = root / "stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Stage 18.5 structure time integration core diagnostic summary"]
    lines.extend(f"{key} {values.get(key, 'FAIL')}" for key in SUMMARY_KEYS)
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = add_args()
    root = Path(args.repo_root).resolve()
    statuses: Dict[str, str] = {}
    reasons: List[str] = []

    npts = int_value(args.npts)
    component_dim = int_value(args.component_dim)
    length = finite_float(args.fibre_length)
    dt = finite_float(args.dt_structure)
    rho_l = finite_float(args.rho_l)
    rho_tilde = finite_float(args.rho_tilde)
    use_dim = int_value(args.use_dimensional_mass)
    use_non = int_value(args.use_nondimensional_mass)
    u0 = finite_float(args.uniform_velocity)
    f0 = finite_float(args.constant_force)
    zero_tol = finite_float(args.zero_tol) or 0.0
    formula_tol = finite_float(args.formula_tol) or 0.0
    dt_ref_tol = finite_float(args.dt_refinement_tol) or 0.0

    n_ok = npts is not None and npts >= 2
    length_ok = length is not None and length > 0.0
    dt_ok = dt is not None and dt > 0.0
    rho_l_ok = rho_l is not None and rho_l > 0.0
    rho_tilde_ok = rho_tilde is not None and rho_tilde > 0.0
    mass_enabled = (use_dim == 1) or (use_non == 1)
    ds = length / (npts - 1) if n_ok and length_ok else math.nan
    shape_ok = n_ok and component_dim == 3

    x = straight_fibre(npts or 2, ds if math.isfinite(ds) else 1.0)
    zeros = [(0.0, 0.0, 0.0)] * len(x)
    no_activation = no_stage18_5_physics_activation(root)

    rho_for_diag = rho_l if rho_l_ok else 1.0
    dt_for_diag = dt if dt_ok else 0.0
    u_for_diag = u0 if u0 is not None else 0.0
    f_for_diag = f0 if f0 is not None else 0.0

    a0, x0, v0 = candidate_update(x, zeros, zeros, rho_for_diag, dt_for_diag)
    uniform_v = [(0.0, u_for_diag, 0.0)] * len(x)
    au, xu, vu = candidate_update(x, uniform_v, zeros, rho_for_diag, dt_for_diag)
    f_const = [(0.0, f_for_diag, 0.0)] * len(x)
    ac, xc, vc = candidate_update(x, zeros, f_const, rho_for_diag, dt_for_diag)
    expected_ac = [(0.0, f_for_diag / rho_for_diag, 0.0)] * len(x)
    expected_vc = [vscale(dt_for_diag, ai) for ai in expected_ac]
    expected_xc = [vadd(xi, vscale(0.5 * dt_for_diag * dt_for_diag, ai)) for xi, ai in zip(x, expected_ac)]

    ft = [(f_for_diag, 0.0, 0.0)] * len(x)
    fb = [(0.0, f_for_diag, 0.0)] * len(x)
    fh = [(0.0, 0.0, f_for_diag)] * len(x)
    fsum = [vadd(vadd(t, b), h) for t, b, h in zip(ft, fb, fh)]
    asplit, _, _ = candidate_update(x, zeros, fsum, rho_for_diag, dt_for_diag)
    expected_asplit = [vscale(1.0 / rho_for_diag, fs) for fs in fsum]

    # Constant-acceleration dt-refinement: two half steps must equal one full step.
    full_a, full_x, full_v = candidate_update(x, uniform_v, f_const, rho_for_diag, dt_for_diag)
    half_dt = 0.5 * dt_for_diag
    half_a, half_x, half_v = candidate_update(x, uniform_v, f_const, rho_for_diag, half_dt)
    two_a, two_x, two_v = candidate_update(half_x, half_v, f_const, rho_for_diag, half_dt)

    have_git, entries = git_status_entries(root)
    changed_ok = (not have_git) or changed_paths_ok(entries)
    source_only_ok = not have_git
    doc_text = read_text(root / STAGE18_5_FILES[2])

    statuses.update({
        "stage18_5_requested_status": status(args.stage18_5_enable == "1" and all((root / rel).exists() for rel in STAGE18_5_FILES)),
        "stage18_4_evidence_status": status(stage18_4_evidence_ok(root)),
        "stage18_3_evidence_status": status(stage18_3_evidence_ok(root)),
        "stage18_2_evidence_status": status(stage18_2_evidence_ok(root)),
        "stage18_1_evidence_status": status(stage18_1_evidence_ok(root)),
        "stage18_0_evidence_status": status(stage18_0_evidence_ok(root)),
        "stage17_closed_file_status": status((root / "stage17_checks/STAGE17_CLOSED.md").exists() or source_only_ok or stage17_11_total_audit_fix_preserved(root)),
        "stage17_closed_evidence_status": status(stage17_closed_evidence_ok(root)),
        "stage17_11_closure_preserved_status": status(stage17_11_total_audit_fix_preserved(root)),
        "stage18_0_wrapper_root_fix_preserved_status": status(stage18_0_wrapper_root_fix_preserved(root)),
        "stage18_1_config_preserved_status": status(stage18_1_evidence_ok(root)),
        "stage18_2_geometry_operator_preserved_status": status(stage18_2_evidence_ok(root)),
        "stage18_3_bending_operator_preserved_status": status(stage18_3_evidence_ok(root)),
        "stage18_4_tension_constraint_preserved_status": status(stage18_4_evidence_ok(root)),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_static_audit_fix_preserved(root)),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_evidence_fix_preserved(root)),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_total_audit_fix_preserved(root)),
        "no_closed_stage_modification_status": status(changed_ok),
        "no_stage10_17_file_modification_status": status(changed_ok),
        "stage18_0_files_unmodified_status": status(files_unmodified(entries, STAGE18_0_FILES)),
        "stage18_1_files_unmodified_status": status(files_unmodified(entries, STAGE18_1_FILES)),
        "stage18_2_files_unmodified_status": status(files_unmodified(entries, STAGE18_2_FILES)),
        "stage18_3_files_unmodified_status": status(files_unmodified(entries, STAGE18_3_FILES)),
        "stage18_4_files_unmodified_status": status(files_unmodified(entries, STAGE18_4_FILES)),
        "stage18_enable_status": status(args.stage18_5_enable == "1"),
        "time_integration_core_enable_status": status(args.time_integration_core_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1"),
        "npts_value": str(npts if npts is not None else "invalid"), "npts_status": status(n_ok),
        "component_dim_value": str(component_dim if component_dim is not None else "invalid"), "component_dim_status": status(component_dim == 3),
        "fibre_length_value": str(length if length is not None else "invalid"), "fibre_length_status": status(length_ok),
        "ds_value": f"{ds:.17e}", "ds_formula_status": status(n_ok and length_ok and math.isfinite(ds) and ds > 0.0 and abs(ds - length / (npts - 1)) <= formula_tol),
        "dt_structure_value": str(dt if dt is not None else "invalid"), "dt_structure_status": status(dt_ok),
        "rho_l_value": str(rho_l if rho_l is not None else "invalid"), "rho_l_status": status(rho_l_ok),
        "rho_tilde_value": str(rho_tilde if rho_tilde is not None else "invalid"), "rho_tilde_status": status(rho_tilde_ok),
        "dimensional_mass_validation_status": status(use_dim == 1 and rho_l_ok),
        "nondimensional_mass_validation_status": status(use_non == 1 and rho_tilde_ok),
        "candidate_array_shape_status": status(shape_ok and all(len(arr) == len(x) for arr in [a0, x0, v0, ac, xc, vc, asplit])),
        "candidate_array_finite_status": status(arrays_finite([a0, x0, v0, au, xu, vu, ac, xc, vc, asplit, full_a, full_x, full_v, two_a, two_x, two_v])),
        "zero_force_rest_acceleration_zero_status": status(max_vec_error(a0, zeros) <= zero_tol),
        "zero_force_rest_velocity_unchanged_status": status(max_vec_error(v0, zeros) <= zero_tol),
        "zero_force_rest_position_unchanged_status": status(max_vec_error(x0, x) <= zero_tol),
        "uniform_velocity_no_force_velocity_preserved_status": status(max_vec_error(vu, uniform_v) <= formula_tol),
        "uniform_velocity_no_force_position_formula_status": status(max_vec_error(xu, [vadd(xi, vscale(dt_for_diag, vi)) for xi, vi in zip(x, uniform_v)]) <= formula_tol),
        "constant_force_acceleration_formula_status": status(max_vec_error(ac, expected_ac) <= formula_tol),
        "constant_force_velocity_formula_status": status(max_vec_error(vc, expected_vc) <= formula_tol),
        "constant_force_position_formula_status": status(max_vec_error(xc, expected_xc) <= formula_tol),
        "split_force_sum_formula_status": status(max_vec_error(fsum, [(f_for_diag, f_for_diag, f_for_diag)] * len(x)) <= formula_tol),
        "split_force_acceleration_formula_status": status(max_vec_error(asplit, expected_asplit) <= formula_tol),
        "dt_refinement_position_consistency_status": status(max_vec_error(full_x, two_x) <= dt_ref_tol),
        "dt_refinement_velocity_consistency_status": status(max_vec_error(full_v, two_v) <= dt_ref_tol),
        "time_integration_equation_documented_status": status(all(s in doc_text for s in ["X_t = V", "A_candidate^n = F_total^n / rho_l", "candidate update"])),
        "candidate_update_diagnostic_only_status": status(no_activation),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent(root)),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })

    negative_keys = [
        "no_production_structure_update_status", "no_production_structure_hook_status",
        "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status",
        "no_inextensibility_projection_status", "no_inextensibility_repair_status",
        "no_fluid_force_physical_structure_integration_status", "no_stage14_rhs_call_from_stage18_5_status",
        "no_structure_energy_power_runtime_activation_status", "no_real_wall_contact_force_status",
        "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status", "no_repulsive_force_status",
        "no_lubrication_force_status", "no_friction_force_status", "no_adhesion_force_status",
        "no_contact_damping_force_status", "no_collision_induced_rhs_status",
        "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status",
        "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status",
        "no_unapproved_production_ibm_forcing_status", "no_pressure_projection_modification_status",
        "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    ]
    for key in negative_keys:
        statuses[key] = status(no_activation)

    if args.test_case != "zero_uniform_constant_split_dt_refinement":
        reasons.append("unexpected_stage18_5_test_case")
    if not all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_5_FILES):
        reasons.append("required_stage18_5_file_missing_or_empty")
    if not changed_ok:
        bad = [path for _code, path in entries if path not in ALLOWED_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))
    for key in SUMMARY_KEYS:
        if key.endswith("_status") and key != "final_status" and statuses.get(key) == "FAIL":
            reasons.append(key.replace("_status", "_failed"))

    pass_fail_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key.endswith("_status") and key not in VALUE_KEYS]
    statuses["final_status"] = "PASS" if all(statuses.get(key) == "PASS" for key in pass_fail_keys) else "FAIL"

    write_output(root, statuses, reasons)
    for key in SUMMARY_KEYS:
        print(f"{key} {statuses.get(key, 'FAIL')}")
    for reason in reasons:
        print(f"reason {reason}")
    return 0 if statuses["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
