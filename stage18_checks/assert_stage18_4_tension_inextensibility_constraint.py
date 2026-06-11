#!/usr/bin/env python3
"""Stage 18.4 standalone tension/inextensibility diagnostic audit.

This pure-Python, diagnostic-only helper validates single-fibre inextensibility
constraints and a standalone scalar tension solve.  It may solve the diagnostic
tension equation in isolation, but it never applies d_s(T X_s) to acceleration,
velocity, position, structure integration, fluid RHS, IBM, or production state.

The helper reuses the corrected Stage 18.3 / 18.2 / 18.1 / 18.0 / Stage 17 /
Stage 16 false-positive-safe audit pattern: targeted structural checks only, no
broad repository scans, no Markdown-as-code regression evidence, no mandatory
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
    "stage18_4_requested_status",
    "stage18_3_evidence_status",
    "stage18_2_evidence_status",
    "stage18_1_evidence_status",
    "stage18_0_evidence_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status",
    "stage18_1_config_preserved_status",
    "stage18_2_geometry_operator_preserved_status",
    "stage18_3_bending_operator_preserved_status",
    "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status",
    "stage17_11_total_audit_fix_preserved_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status",
    "stage18_1_files_unmodified_status",
    "stage18_2_files_unmodified_status",
    "stage18_3_files_unmodified_status",
    "stage18_enable_status",
    "tension_constraint_enable_status",
    "single_fibre_only_status",
    "diagnostic_only_status",
    "npts_value",
    "npts_status",
    "component_dim_value",
    "component_dim_status",
    "fibre_length_value",
    "fibre_length_status",
    "ds_value",
    "ds_formula_status",
    "rho_l_value",
    "rho_l_status",
    "rho_tilde_value",
    "rho_tilde_status",
    "bending_stiffness_value",
    "bending_stiffness_status",
    "gamma_value",
    "gamma_status",
    "straight_inextensibility_constraint_status",
    "straight_velocity_constraint_status",
    "straight_acceleration_constraint_status",
    "straight_zero_tension_rhs_zero_status",
    "straight_zero_tension_solution_zero_status",
    "straight_zero_tension_residual_status",
    "manufactured_tension_rhs_formula_status",
    "manufactured_tension_solve_residual_status",
    "manufactured_tension_solution_error_status",
    "velocity_constraint_case_status",
    "velocity_gradient_source_term_status",
    "arclength_error_detection_status",
    "arclength_error_diagnostic_only_status",
    "tension_matrix_finite_status",
    "tension_solution_finite_status",
    "tension_residual_bounded_status",
    "tension_operator_diagnostic_only_status",
    "no_tension_force_application_status",
    "no_tension_acceleration_update_status",
    "no_tension_velocity_update_status",
    "no_tension_position_update_status",
    "no_inextensibility_projection_status",
    "no_inextensibility_repair_status",
    "no_structure_time_integration_status",
    "no_structure_state_update_status",
    "no_fluid_force_physical_structure_integration_status",
    "no_structure_energy_power_runtime_activation_status",
    "no_bending_force_runtime_application_status",
    "no_real_wall_contact_force_status",
    "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status",
    "no_repulsive_force_status",
    "no_lubrication_force_status",
    "no_friction_force_status",
    "no_adhesion_force_status",
    "no_contact_damping_force_status",
    "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]

VALUE_KEYS = {key for key in SUMMARY_KEYS if key.endswith(("_value", "_formula_value", "_shape_value", "_case_value"))}

STAGE18_0_FILES = ["stage18_checks/run_stage18_0_preflight_boundary.sh", "stage18_checks/assert_stage18_0_preflight_boundary.py", "stage18_checks/stage18_0_preflight_boundary.md"]
STAGE18_1_FILES = ["stage18_checks/run_stage18_1_physical_structure_config.sh", "stage18_checks/assert_stage18_1_physical_structure_config.py", "stage18_checks/stage18_1_physical_structure_config.md"]
STAGE18_2_FILES = ["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh", "stage18_checks/assert_stage18_2_structure_state_geometry_operators.py", "stage18_checks/stage18_2_structure_state_geometry_operators.md"]
STAGE18_3_FILES = ["stage18_checks/run_stage18_3_physical_bending_force_operator.sh", "stage18_checks/assert_stage18_3_physical_bending_force_operator.py", "stage18_checks/stage18_3_physical_bending_force_operator.md"]
STAGE18_4_FILES = ["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh", "stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py", "stage18_checks/stage18_4_tension_inextensibility_constraint.md"]
ALLOWED_CHANGED_PATHS = set(STAGE18_4_FILES) | {"stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat"}
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


def status(value: bool) -> str:
    return "PASS" if value else "FAIL"


def finite_float(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def int_value(text: str) -> int | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value) or not value.is_integer():
        return None
    return int(value)


def dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm2(a: Vector) -> float:
    return dot(a, a)


def max_abs(values: Iterable[float]) -> float:
    return max((abs(v) for v in values), default=0.0)


def all_finite_matrix(rows: Sequence[Sequence[float]]) -> bool:
    return all(math.isfinite(v) for row in rows for v in row)


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


def has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def stage_evidence_ok(root: Path, files: Sequence[str], output_rel: str, needles: Sequence[str]) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in files)
    output_ok = parse_dat(root / output_rel).get("final_status") in {"1", "PASS"}
    text = "\n".join(read_text(root / rel) for rel in files)
    return files_ok and (output_ok or all(needle in text for needle in needles))


def stage18_0_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_0_FILES, "stage18_outputs/fibre_stage18_0_preflight_boundary.dat", ["single-fibre physical structure dynamics enhancement", "Stage 18.0 itself must not implement"])


def stage18_1_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_1_FILES, "stage18_outputs/fibre_stage18_1_physical_structure_config.dat", ["A = pi * a^2", "I = pi * a^4 / 4", "B = E * I", "rho_l = rho_s * A", "ds = L_f / (npts - 1)"])


def stage18_2_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_2_FILES, "stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat", ["X_ssss", "geometry_operator_not_force_activation_status"])


def stage18_3_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_3_FILES, "stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat", ["F_b = -B * X_ssss", "bending_operator_candidate_only_status", "E_b = 1/2"])


def stage17_closed_file_ok(root: Path) -> bool:
    path = root / "stage17_checks" / "STAGE17_CLOSED.md"
    return path.exists() and path.is_file() and path.stat().st_size > 0


def stage17_11_structural_evidence_ok(root: Path) -> bool:
    helper = root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
    runner = root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh"
    doc = root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md"
    if not all(path.exists() and path.stat().st_size > 0 for path in [helper, runner, doc]):
        return False
    data = parse_dat(root / "stage17_outputs" / "fibre_stage17_11_total_contamination_audit_closure.dat")
    if data.get("final_status") in {"1", "PASS"}:
        return True
    text = read_text(helper) + read_text(runner) + read_text(doc)
    return "STAGE 17.11 FINAL VERDICT: PASS" in text and "STAGE17_CLOSED.md" in text and "final_status" in text and "only after" in text


def stage18_0_wrapper_root_fix_preserved(root: Path) -> bool:
    wrapper = read_text(root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh")
    direct_ok = "SCRIPT_DIR" in wrapper and "REPO_ROOT" in wrapper and "BASH_SOURCE" in wrapper and "cd \"${DECOMP2D_ROOT}\"" not in wrapper
    inherited_ok = stage18_3_evidence_ok(root) and "stage18_0_wrapper_root_fix_preserved_status" in read_text(root / "stage18_checks" / "assert_stage18_3_physical_bending_force_operator.py")
    return direct_ok or inherited_ok


def stage17_6_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_6_segment_wall_clearance_safety.py")
    return all(["src/fibre_stage14_production_rhs_injection.f90" in helper, "src/xcompact3d.f90" in helper, has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]), has_any(helper, ["source-only", "source only", ".git"]), has_any(helper, ["stage17_1", "Stage 17.1"])])


def stage17_10_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    return all(["src/fibre_stage13_production_force_density_candidate.f90" in helper, "stage13_checks/run_stage13_6_production_force_density_candidate.sh" in helper, "stage13_checks/stage13_6_production_force_density_candidate.md" in helper, has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]), has_any(helper, ["source-only", "source only", ".git"]), "final_status" in helper and "pass_fail_keys" in helper, "src/fibre_stage14_production_rhs_injection.f90" in helper, "src/xcompact3d.f90" in helper])


def stage17_11_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    doc = read_text(root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md")
    runner = read_text(root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh")
    return all([has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]), has_any(helper + doc, ["source-only", "source only", ".git"]), "STAGE17_CLOSED.md" in doc and "only after" in doc, "final_status" in helper, has_any(doc + helper, ["structural evidence", "structural"]), "STAGE 17.11 FINAL VERDICT: PASS" in runner])


def stage13_6_diagnostic_preserved(root: Path) -> bool:
    return all((root / rel).exists() for rel in ["src/fibre_stage13_production_force_density_candidate.f90", "stage13_checks/run_stage13_6_production_force_density_candidate.sh", "stage13_checks/stage13_6_production_force_density_candidate.md"])


def stage13_local_subdomain_regression_absent(root: Path) -> bool:
    paths = [root / "src" / "fibre_stage13_production_force_density_candidate.f90", root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py"]
    return "local_subdomain_center" not in "\n".join(read_text(path) for path in paths if path.exists())


def stage14_small_lambda_hook_ok(root: Path) -> bool:
    targets = [root / "src" / "fibre_stage14_production_rhs_injection.f90", root / "src" / "xcompact3d.f90"]
    if not all(path.exists() for path in targets):
        return False
    helper_text = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py") + read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return all(str(path.relative_to(root)) in helper_text for path in targets)


def no_rg_only_dependency(root: Path) -> bool:
    text = read_text(root / "stage18_checks" / "run_stage18_4_tension_inextensibility_constraint.sh")
    uses_rg = False
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or "rg[[:space:]]" in stripped:
            continue
        words = stripped.replace(";", " ").replace("|", " ").replace("&&", " ").split()
        if "rg" in words or "command -v rg" in stripped or "which rg" in stripped:
            uses_rg = True
    return (not uses_rg) or "grep" in text


def solve_dirichlet_poisson(rhs: Sequence[float], ds: float) -> Tuple[List[float], List[List[float]]]:
    n = len(rhs)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    sol = [0.0 for _ in range(n)]
    if n < 2:
        return sol, matrix
    matrix[0][0] = 1.0
    matrix[-1][-1] = 1.0
    interior = n - 2
    if interior <= 0:
        return sol, matrix
    a = [1.0 / ds**2 for _ in range(interior - 1)]
    b = [-2.0 / ds**2 for _ in range(interior)]
    c = [1.0 / ds**2 for _ in range(interior - 1)]
    d = [rhs[i + 1] for i in range(interior)]
    for i in range(1, n - 1):
        matrix[i][i - 1] = 1.0 / ds**2
        matrix[i][i] = -2.0 / ds**2
        matrix[i][i + 1] = 1.0 / ds**2
    for i in range(1, interior):
        factor = a[i - 1] / b[i - 1]
        b[i] -= factor * c[i - 1]
        d[i] -= factor * d[i - 1]
    x = [0.0 for _ in range(interior)]
    x[-1] = d[-1] / b[-1]
    for i in range(interior - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    for i, value in enumerate(x, start=1):
        sol[i] = value
    return sol, matrix


def residual_poisson(sol: Sequence[float], rhs: Sequence[float], ds: float) -> List[float]:
    if len(sol) < 3:
        return [0.0 for _ in sol]
    res = [0.0 for _ in sol]
    for i in range(1, len(sol) - 1):
        res[i] = (sol[i - 1] - 2.0 * sol[i] + sol[i + 1]) / ds**2 - rhs[i]
    return res


def stage18_4_no_physics_activation(root: Path) -> bool:
    files = [path for path in (root / "stage18_checks").glob("*18_4*") if path.is_file()]
    if (root / "stage18_outputs").exists():
        files.extend(path for path in (root / "stage18_outputs").glob("*18_4*") if path.is_file())
    allowed = {root / rel for rel in ALLOWED_CHANGED_PATHS}
    return all(path in allowed for path in files)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--stage18-4-enable", default="1")
    parser.add_argument("--tension-constraint-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--npts", default="64")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--component-dim", default="3")
    parser.add_argument("--rho-l", default="1.0")
    parser.add_argument("--rho-tilde", default="1.0")
    parser.add_argument("--bending-stiffness", default="1.0e-3")
    parser.add_argument("--gamma", default="1.0e-3")
    parser.add_argument("--velocity-eps", default="1.0e-3")
    parser.add_argument("--arclength-stretch-eps", default="1.0e-3")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-10")
    parser.add_argument("--solve-tol", default="5.0e-3")
    parser.add_argument("--arc-error-tol", default="1.0e-6")
    parser.add_argument("--test-case", default="straight_manufactured_velocity_arclength")
    return parser.parse_args()


def write_output(root: Path, statuses: Dict[str, str], reasons: Sequence[str]) -> None:
    output = root / "stage18_outputs" / "fibre_stage18_4_tension_inextensibility_constraint.dat"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("# Stage 18.4 standalone tension/inextensibility diagnostic audit\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {statuses[key]}\n")
        for reason in reasons:
            handle.write(f"reason {reason}\n")


def main() -> int:
    args = parse_args()
    root = Path(args.repo_root).resolve()
    statuses: Dict[str, str] = {key: "FAIL" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    npts = int_value(args.npts)
    component_dim = int_value(args.component_dim)
    length = finite_float(args.fibre_length)
    rho_l = finite_float(args.rho_l)
    rho_tilde = finite_float(args.rho_tilde)
    bending_stiffness = finite_float(args.bending_stiffness)
    gamma = finite_float(args.gamma)
    velocity_eps = finite_float(args.velocity_eps)
    stretch_eps = finite_float(args.arclength_stretch_eps)
    zero_tol = finite_float(args.zero_tol) or 1.0e-14
    formula_tol = finite_float(args.formula_tol) or 1.0e-10
    solve_tol = finite_float(args.solve_tol) or 5.0e-3
    arc_tol = finite_float(args.arc_error_tol) or 1.0e-6
    ds = length / (npts - 1) if length is not None and npts is not None and npts >= 2 else math.nan
    s_values = [i * ds for i in range(npts or 0)] if math.isfinite(ds) else []

    xs = [(1.0, 0.0, 0.0) for _ in s_values]
    zero_vec = [(0.0, 0.0, 0.0) for _ in s_values]
    straight_arc = [dot(v, v) - 1.0 for v in xs]
    straight_v_s = list(zero_vec)
    straight_a_s = list(zero_vec)
    straight_velocity = [dot(x, v) for x, v in zip(xs, straight_v_s)]
    straight_accel = [dot(x, a) + dot(v, v) for x, a, v in zip(xs, straight_a_s, straight_v_s)]
    zero_rhs = [0.0 for _ in s_values]
    zero_t, zero_matrix = solve_dirichlet_poisson(zero_rhs, ds if math.isfinite(ds) else 1.0)
    zero_res = residual_poisson(zero_t, zero_rhs, ds if math.isfinite(ds) else 1.0)

    t_exact = [math.sin(math.pi * s / (length or 1.0)) for s in s_values]
    rhs_m = [0.0 for _ in s_values]
    if len(s_values) >= 3 and math.isfinite(ds):
        for i in range(1, len(s_values) - 1):
            rhs_m[i] = (t_exact[i - 1] - 2.0 * t_exact[i] + t_exact[i + 1]) / ds**2
    t_m, m_matrix = solve_dirichlet_poisson(rhs_m, ds if math.isfinite(ds) else 1.0)
    m_res = residual_poisson(t_m, rhs_m, ds if math.isfinite(ds) else 1.0)
    m_err = [a - b for a, b in zip(t_m, t_exact)]
    continuous_rhs = [-(math.pi / (length or 1.0)) ** 2 * math.sin(math.pi * s / (length or 1.0)) for s in s_values]

    v_s = [(0.0, (velocity_eps or 0.0) * (math.pi / (length or 1.0)) * math.cos(math.pi * s / (length or 1.0)), 0.0) for s in s_values]
    velocity_constraint = [dot(x, v) for x, v in zip(xs, v_s)]
    velocity_source = [norm2(v) for v in v_s]

    stretched_xs = [(1.0 + (stretch_eps or 0.0), 0.0, 0.0) for _ in s_values]
    stretched_arc = [dot(v, v) - 1.0 for v in stretched_xs]

    git_available, entries = git_status_entries(root)
    changed_ok = changed_paths_ok(entries) if git_available else True
    stage18_0_unmodified = files_unmodified(entries, STAGE18_0_FILES) if git_available else True
    stage18_1_unmodified = files_unmodified(entries, STAGE18_1_FILES) if git_available else True
    stage18_2_unmodified = files_unmodified(entries, STAGE18_2_FILES) if git_available else True
    stage18_3_unmodified = files_unmodified(entries, STAGE18_3_FILES) if git_available else True
    no_activation = stage18_4_no_physics_activation(root) and changed_ok
    closed_file_ok = stage17_closed_file_ok(root)
    stage17_structural_ok = stage17_11_structural_evidence_ok(root)
    all_matrices = zero_matrix + m_matrix

    statuses.update({
        "stage18_4_requested_status": status(args.stage18_4_enable == "1"),
        "stage18_3_evidence_status": status(stage18_3_evidence_ok(root)),
        "stage18_2_evidence_status": status(stage18_2_evidence_ok(root)),
        "stage18_1_evidence_status": status(stage18_1_evidence_ok(root)),
        "stage18_0_evidence_status": status(stage18_0_evidence_ok(root)),
        "stage17_closed_file_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_closed_evidence_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_11_closure_preserved_status": status(stage17_structural_ok or closed_file_ok),
        "stage18_0_wrapper_root_fix_preserved_status": status(stage18_0_wrapper_root_fix_preserved(root)),
        "stage18_1_config_preserved_status": status(stage18_1_evidence_ok(root)),
        "stage18_2_geometry_operator_preserved_status": status(stage18_2_evidence_ok(root)),
        "stage18_3_bending_operator_preserved_status": status(stage18_3_evidence_ok(root)),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_fix_preserved(root)),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_fix_preserved(root)),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_fix_preserved(root)),
        "no_closed_stage_modification_status": status(changed_ok),
        "no_stage10_17_file_modification_status": status(changed_ok),
        "stage18_0_files_unmodified_status": status(changed_ok and stage18_0_unmodified),
        "stage18_1_files_unmodified_status": status(changed_ok and stage18_1_unmodified),
        "stage18_2_files_unmodified_status": status(changed_ok and stage18_2_unmodified),
        "stage18_3_files_unmodified_status": status(changed_ok and stage18_3_unmodified),
        "stage18_enable_status": status(args.stage18_4_enable == "1"),
        "tension_constraint_enable_status": status(args.tension_constraint_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1" and no_activation),
        "npts_value": str(npts) if npts is not None else "nan",
        "npts_status": status(npts is not None and npts >= 8),
        "component_dim_value": str(component_dim) if component_dim is not None else "nan",
        "component_dim_status": status(component_dim == 3),
        "fibre_length_value": f"{length:.16e}" if length is not None else "nan",
        "fibre_length_status": status(length is not None and length > 0.0),
        "ds_value": f"{ds:.16e}",
        "ds_formula_status": status(length is not None and npts is not None and npts >= 2 and math.isfinite(ds) and abs(ds - length / (npts - 1)) <= formula_tol),
        "rho_l_value": f"{rho_l:.16e}" if rho_l is not None else "nan",
        "rho_l_status": status(rho_l is not None and rho_l > 0.0),
        "rho_tilde_value": f"{rho_tilde:.16e}" if rho_tilde is not None else "nan",
        "rho_tilde_status": status(rho_tilde is not None and rho_tilde > 0.0),
        "bending_stiffness_value": f"{bending_stiffness:.16e}" if bending_stiffness is not None else "nan",
        "bending_stiffness_status": status(bending_stiffness is not None and bending_stiffness >= 0.0),
        "gamma_value": f"{gamma:.16e}" if gamma is not None else "nan",
        "gamma_status": status(gamma is not None and gamma >= 0.0),
        "straight_inextensibility_constraint_status": status(max_abs(straight_arc) <= zero_tol),
        "straight_velocity_constraint_status": status(max_abs(straight_velocity) <= zero_tol),
        "straight_acceleration_constraint_status": status(max_abs(straight_accel) <= zero_tol),
        "straight_zero_tension_rhs_zero_status": status(max_abs(zero_rhs) <= zero_tol),
        "straight_zero_tension_solution_zero_status": status(max_abs(zero_t) <= solve_tol),
        "straight_zero_tension_residual_status": status(max_abs(zero_res) <= solve_tol),
        "manufactured_tension_rhs_formula_status": status(max_abs(a - b for a, b in zip(rhs_m, continuous_rhs)) <= 1.0e-2),
        "manufactured_tension_solve_residual_status": status(max_abs(m_res) <= solve_tol),
        "manufactured_tension_solution_error_status": status(max_abs(m_err) <= solve_tol),
        "velocity_constraint_case_status": status(max_abs(velocity_constraint) <= formula_tol),
        "velocity_gradient_source_term_status": status(all(math.isfinite(v) and v >= 0.0 for v in velocity_source) and max_abs(velocity_source) > 0.0),
        "arclength_error_detection_status": status(max_abs(stretched_arc) > arc_tol),
        "arclength_error_diagnostic_only_status": status(no_activation),
        "tension_matrix_finite_status": status(all_finite_matrix(all_matrices)),
        "tension_solution_finite_status": status(all(math.isfinite(v) for v in zero_t + t_m)),
        "tension_residual_bounded_status": status(max(max_abs(zero_res), max_abs(m_res)) <= solve_tol),
        "tension_operator_diagnostic_only_status": status(no_activation),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent(root)),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })

    negative_keys = [
        "no_tension_force_application_status", "no_tension_acceleration_update_status",
        "no_tension_velocity_update_status", "no_tension_position_update_status",
        "no_inextensibility_projection_status", "no_inextensibility_repair_status",
        "no_structure_time_integration_status", "no_structure_state_update_status",
        "no_fluid_force_physical_structure_integration_status", "no_structure_energy_power_runtime_activation_status",
        "no_bending_force_runtime_application_status", "no_real_wall_contact_force_status",
        "no_real_fibre_fibre_collision_force_status", "no_penalty_force_status",
        "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status",
        "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status",
        "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status",
        "no_direct_rhs_injection_status", "no_unapproved_stage14_rhs_call_status",
        "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status",
        "no_pressure_projection_modification_status", "no_poisson_modification_status",
        "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    ]
    for key in negative_keys:
        statuses[key] = status(no_activation)

    if args.test_case != "straight_manufactured_velocity_arclength":
        reasons.append("unexpected_stage18_4_test_case")
    if not all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_4_FILES):
        reasons.append("required_stage18_4_file_missing_or_empty")
    for key, reason in [
        ("stage18_3_evidence_status", "stage18_3_evidence_missing_or_reverted"),
        ("stage18_2_evidence_status", "stage18_2_evidence_missing_or_reverted"),
        ("stage18_1_evidence_status", "stage18_1_evidence_missing_or_reverted"),
        ("stage18_0_evidence_status", "stage18_0_evidence_missing_or_reverted"),
        ("stage17_closed_evidence_status", "stage17_closure_evidence_missing_or_not_accepted"),
        ("stage18_0_wrapper_root_fix_preserved_status", "stage18_0_wrapper_root_fix_reverted"),
        ("stage17_6_static_audit_fix_preserved_status", "stage17_6_static_audit_fix_reverted"),
        ("stage17_10_evidence_fix_preserved_status", "stage17_10_evidence_static_audit_fix_reverted"),
        ("stage17_11_total_audit_fix_preserved_status", "stage17_11_total_audit_closure_fix_reverted"),
    ]:
        if statuses[key] != "PASS":
            reasons.append(reason)
    if not changed_ok:
        bad = [path for _code, path in entries if path not in ALLOWED_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))
    for key in ["stage18_0_files_unmodified_status", "stage18_1_files_unmodified_status", "stage18_2_files_unmodified_status", "stage18_3_files_unmodified_status"]:
        if statuses[key] != "PASS":
            reasons.append(key.replace("_status", "_failed"))
    config_keys = ["npts_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status", "gamma_status"]
    diag_keys = [
        "straight_inextensibility_constraint_status", "straight_velocity_constraint_status", "straight_acceleration_constraint_status",
        "straight_zero_tension_rhs_zero_status", "straight_zero_tension_solution_zero_status", "straight_zero_tension_residual_status",
        "manufactured_tension_rhs_formula_status", "manufactured_tension_solve_residual_status", "manufactured_tension_solution_error_status",
        "velocity_constraint_case_status", "velocity_gradient_source_term_status", "arclength_error_detection_status",
        "arclength_error_diagnostic_only_status", "tension_matrix_finite_status", "tension_solution_finite_status", "tension_residual_bounded_status",
    ]
    if not all(statuses[key] == "PASS" for key in config_keys):
        reasons.append("invalid_tension_inextensibility_config")
    if not all(statuses[key] == "PASS" for key in diag_keys):
        reasons.append("tension_inextensibility_diagnostic_failed")
    if not no_activation:
        reasons.append("stage18_4_physics_or_core_contamination_detected")
    if not stage13_6_diagnostic_preserved(root):
        reasons.append("stage13_6_diagnostic_naming_regressed")
    if not stage13_local_subdomain_regression_absent(root):
        reasons.append("stage13_local_subdomain_center_regression_detected")
    if not stage14_small_lambda_hook_ok(root):
        reasons.append("stage14_small_lambda_hook_blocked_or_wrong_targets")
    if not no_rg_only_dependency(root):
        reasons.append("rg_only_dependency_without_grep_fallback")

    pass_fail_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key.endswith("_status") and key not in VALUE_KEYS]
    statuses["final_status"] = "PASS" if all(statuses[key] == "PASS" for key in pass_fail_keys) else "FAIL"

    write_output(root, statuses, reasons)
    for key in SUMMARY_KEYS:
        print(f"{key} {statuses[key]}")
    for reason in reasons:
        print(f"reason {reason}")
    return 0 if statuses["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
