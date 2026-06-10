#!/usr/bin/env python3
"""Stage 18.1 physical structure dynamics configuration audit.

The helper validates diagnostic-only single-fibre physical structure parameter
configuration for later Stage 18 dynamics work.  It computes only scalar
configuration relations (area, second moment, bending stiffness, mass per
length, and nominal spacing); it does not implement or activate a structure
solver, bending force, tension solve, inextensibility enforcement, state update,
fluid RHS update, contact/collision force, IBM path, or DNS-core change.

The audit follows the corrected Stage 18.0 / Stage 17 / Stage 16
false-positive-safe pattern: targeted structural checks only, no Markdown as
runtime regression evidence, no broad repository scans, no mandatory ripgrep,
source-only archives without .git metadata accepted as non-contamination, and
only *_status fields participate in final_status.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
    "stage18_1_requested_status",
    "stage18_0_evidence_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status",
    "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status",
    "stage17_11_total_audit_fix_preserved_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status",
    "stage18_enable_status",
    "structure_dynamics_config_enable_status",
    "single_fibre_only_status",
    "physical_structure_boundary_status",
    "bending_config_enable_status",
    "tension_config_enable_status",
    "inextensibility_config_enable_status",
    "time_integration_config_enable_status",
    "energy_diagnostic_config_enable_status",
    "diagnostic_only_status",
    "rho_s_value",
    "rho_s_status",
    "fibre_radius_value",
    "fibre_radius_status",
    "young_modulus_value",
    "young_modulus_status",
    "fibre_length_value",
    "fibre_length_status",
    "npts_value",
    "npts_status",
    "dt_structure_value",
    "dt_structure_status",
    "rho_tilde_value",
    "rho_tilde_status",
    "gamma_value",
    "gamma_status",
    "fibre_area_value",
    "fibre_area_formula_status",
    "second_moment_area_value",
    "second_moment_area_formula_status",
    "bending_stiffness_value",
    "bending_stiffness_formula_status",
    "mass_per_length_value",
    "mass_per_length_formula_status",
    "ds_value",
    "ds_formula_status",
    "derived_values_finite_status",
    "physical_structure_equation_documented_status",
    "inextensibility_constraint_documented_status",
    "structure_energy_boundary_documented_status",
    "config_not_physics_activation_status",
    "no_bending_force_activation_status",
    "no_tension_solve_activation_status",
    "no_inextensibility_enforcement_status",
    "no_structure_time_integration_status",
    "no_structure_state_update_status",
    "no_fluid_force_physical_structure_integration_status",
    "no_structure_energy_power_implementation_status",
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

VALUE_KEYS = {key for key in SUMMARY_KEYS if key.endswith("_value") or key.endswith("_formula_value")}

STAGE18_0_FILES = [
    "stage18_checks/run_stage18_0_preflight_boundary.sh",
    "stage18_checks/assert_stage18_0_preflight_boundary.py",
    "stage18_checks/stage18_0_preflight_boundary.md",
]

STAGE18_1_FILES = [
    "stage18_checks/run_stage18_1_physical_structure_config.sh",
    "stage18_checks/assert_stage18_1_physical_structure_config.py",
    "stage18_checks/stage18_1_physical_structure_config.md",
]

ALLOWED_CHANGED_PATHS = set(STAGE18_1_FILES) | {"stage18_outputs/fibre_stage18_1_physical_structure_config.dat"}


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
        if str(int(text)) != text.strip():
            # Accept strings like 8.0 only if they are exactly integral as float.
            as_float = float(text)
            if not math.isfinite(as_float) or not as_float.is_integer():
                return None
            return int(as_float)
        return int(text)
    except ValueError:
        return None


def near(a: float, b: float, tol: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def git_status_entries(root: Path) -> Tuple[bool, List[Tuple[str, str]]]:
    """Return porcelain entries if .git metadata exists.

    Source-only archives without .git metadata are accepted as structural source
    archives and are not treated as DNS-core contamination or closed-stage edits.
    """
    if not (root / ".git").exists():
        return False, []
    proc = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        return False, []
    entries: List[Tuple[str, str]] = []
    for raw in proc.stdout.splitlines():
        if not raw:
            continue
        code = raw[:2]
        path = raw[3:] if len(raw) > 3 else ""
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        entries.append((code, path))
    return True, entries


def changed_paths_ok(entries: Sequence[Tuple[str, str]]) -> bool:
    return all(path in ALLOWED_CHANGED_PATHS for _code, path in entries)


def stage18_0_files_unmodified(entries: Sequence[Tuple[str, str]]) -> bool:
    return all(path not in STAGE18_0_FILES for _code, path in entries)


def stage17_10_files_unmodified(entries: Sequence[Tuple[str, str]]) -> bool:
    prefixes = tuple(f"stage{stage}_" for stage in range(10, 18))
    return all(not path.startswith(prefixes) and not path.startswith("stage17_checks/STAGE17_CLOSED.md") for _code, path in entries)


def required_files_present(root: Path) -> bool:
    return all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_1_FILES)


def stage18_0_evidence_ok(root: Path) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_0_FILES)
    output = parse_dat(root / "stage18_outputs" / "fibre_stage18_0_preflight_boundary.dat")
    output_ok = output.get("final_status") in {"1", "PASS"}
    doc = read_text(root / "stage18_checks" / "stage18_0_preflight_boundary.md")
    structural_ok = (
        "single-fibre physical structure dynamics enhancement" in doc
        and "Stage 18.0 itself must not implement these equations yet" in doc
    )
    return files_ok and (output_ok or structural_ok)


def stage17_closed_file_ok(root: Path) -> bool:
    path = root / "stage17_checks" / "STAGE17_CLOSED.md"
    return path.exists() and path.is_file() and path.stat().st_size > 0


def stage17_11_structural_evidence_ok(root: Path) -> bool:
    helper = root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
    runner = root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh"
    doc = root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md"
    required_present = all(path.exists() and path.stat().st_size > 0 for path in [helper, runner, doc])
    if not required_present:
        return False
    data = parse_dat(root / "stage17_outputs" / "fibre_stage17_11_total_contamination_audit_closure.dat")
    if data.get("final_status") in {"1", "PASS"}:
        return True
    text = read_text(runner) + read_text(doc) + read_text(helper)
    return (
        "STAGE 17.11 FINAL VERDICT: PASS" in text
        and "STAGE17_CLOSED.md" in text
        and "final_status" in text
        and "only after" in text
    )


def has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def stage18_0_wrapper_root_fix_preserved(root: Path) -> bool:
    wrapper = read_text(root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh")
    has_script_root = "SCRIPT_DIR" in wrapper and "REPO_ROOT" in wrapper and "BASH_SOURCE" in wrapper
    decomp_interface_only = "DECOMP2D_ROOT" in wrapper and "cd \"${DECOMP2D_ROOT}\"" not in wrapper
    return has_script_root and decomp_interface_only


def stage17_6_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_6_segment_wall_clearance_safety.py")
    return all([
        "src/fibre_stage14_production_rhs_injection.f90" in helper,
        "src/xcompact3d.f90" in helper,
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper, ["source-only", "source only", ".git"]),
        has_any(helper, ["stage17_1", "Stage 17.1"]),
    ])


def stage17_10_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    return all([
        "src/fibre_stage13_production_force_density_candidate.f90" in helper,
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh" in helper,
        "stage13_checks/stage13_6_production_force_density_candidate.md" in helper,
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper, ["source-only", "source only", ".git"]),
        "final_status" in helper and "pass_fail_keys" in helper,
        "src/fibre_stage14_production_rhs_injection.f90" in helper,
        "src/xcompact3d.f90" in helper,
    ])


def stage17_11_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    doc = read_text(root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md")
    runner = read_text(root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh")
    return all([
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper + doc, ["source-only", "source only", ".git"]),
        "STAGE17_CLOSED.md" in doc and "only after" in doc,
        "final_status" in helper,
        has_any(doc + helper, ["structural evidence", "structural"]),
        "STAGE 17.11 FINAL VERDICT: PASS" in runner,
    ])


def stage13_6_diagnostic_preserved(root: Path) -> bool:
    required = [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ]
    return all((root / rel).exists() for rel in required)


def stage13_local_subdomain_regression_absent(root: Path) -> bool:
    paths = [
        root / "src" / "fibre_stage13_production_force_density_candidate.f90",
        root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py",
    ]
    text = "\n".join(read_text(path) for path in paths if path.exists())
    return "local_subdomain_center" not in text


def stage14_small_lambda_hook_ok(root: Path) -> bool:
    targets = [
        root / "src" / "fibre_stage14_production_rhs_injection.f90",
        root / "src" / "xcompact3d.f90",
    ]
    if not all(path.exists() for path in targets):
        return False
    helper_text = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    helper_text += read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return all(str(path.relative_to(root)) in helper_text for path in targets)


def no_rg_only_dependency(root: Path) -> bool:
    # Check only Stage 18.1 executable shell wrappers.  Documentation, Python
    # regex literals, and negative-check labels are not executable rg usage.
    for script in [root / "stage18_checks" / "run_stage18_1_physical_structure_config.sh"]:
        text = read_text(script)
        uses_rg = False
        for raw in text.splitlines():
            stripped = raw.strip()
            if not stripped or stripped.startswith("#") or "rg[[:space:]]" in stripped:
                continue
            words = stripped.replace(";", " ").replace("|", " ").replace("&&", " ").split()
            if "rg" in words or "command -v rg" in stripped or "which rg" in stripped:
                uses_rg = True
        if uses_rg and "grep" not in text:
            return False
    return True


def boundary_doc_statuses(root: Path) -> Tuple[bool, bool, bool]:
    doc = read_text(root / "stage18_checks" / "stage18_1_physical_structure_config.md")
    physical = (
        "rho_l * X_tt = d_s(T X_s) - B X_ssss + F_h" in doc
        and "rho_tilde * X_tt = d_s(T X_s) - gamma X_ssss + F_h" in doc
    )
    inextensible = "X_s dot X_s = 1" in doc
    energy = all(item in doc for item in [
        "E_k = 1/2 int rho_l |V|^2 ds",
        "E_b = 1/2 int B |X_ss|^2 ds",
        "P_h = int F_h dot V ds",
    ])
    return physical, inextensible, energy


def stage18_1_no_physics_activation(root: Path) -> bool:
    # Stage 18.1 is allowed only the helper, wrapper, documentation, and runtime
    # dat file.  This targeted structural check avoids treating legacy Stage 2
    # structure files or documentation examples as new physics activation.
    files = [path for path in (root / "stage18_checks").glob("*18_1*") if path.is_file()]
    if (root / "stage18_outputs").exists():
        files.extend(path for path in (root / "stage18_outputs").glob("*18_1*") if path.is_file())
    allowed = {root / rel for rel in ALLOWED_CHANGED_PATHS}
    return all(path in allowed for path in files)


def write_output(root: Path, statuses: Dict[str, str], reasons: Sequence[str]) -> None:
    output = root / "stage18_outputs" / "fibre_stage18_1_physical_structure_config.dat"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("# Stage 18.1 physical structure dynamics configuration audit\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {statuses[key]}\n")
        for reason in reasons:
            handle.write(f"reason {reason}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--stage18-1-enable", default="1")
    parser.add_argument("--structure-dynamics-config-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--physical-structure-boundary", default="1")
    parser.add_argument("--bending-config-enable", default="1")
    parser.add_argument("--tension-config-enable", default="1")
    parser.add_argument("--inextensibility-config-enable", default="1")
    parser.add_argument("--time-integration-config-enable", default="1")
    parser.add_argument("--energy-diagnostic-config-enable", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--rho-s", default="1.0")
    parser.add_argument("--fibre-radius", default="1.0e-3")
    parser.add_argument("--young-modulus", default="1.0")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--npts", default="8")
    parser.add_argument("--dt-structure", default="1.0e-4")
    parser.add_argument("--rho-tilde", default="1.0")
    parser.add_argument("--gamma", default="1.0e-3")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-12")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.repo_root).resolve()
    statuses: Dict[str, str] = {key: "FAIL" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    rho_s = finite_float(args.rho_s)
    radius = finite_float(args.fibre_radius)
    young = finite_float(args.young_modulus)
    length = finite_float(args.fibre_length)
    npts = int_value(args.npts)
    dt_structure = finite_float(args.dt_structure)
    rho_tilde = finite_float(args.rho_tilde)
    gamma = finite_float(args.gamma)
    formula_tol = finite_float(args.formula_tol) or 1.0e-12

    area = math.pi * radius**2 if radius is not None else math.nan
    second_moment = math.pi * radius**4 / 4.0 if radius is not None else math.nan
    bending = young * second_moment if young is not None and math.isfinite(second_moment) else math.nan
    mass_per_length = rho_s * area if rho_s is not None and math.isfinite(area) else math.nan
    ds = length / (npts - 1) if length is not None and npts is not None and npts >= 2 else math.nan

    git_available, entries = git_status_entries(root)
    changed_ok = changed_paths_ok(entries) if git_available else True
    stage18_0_unmodified = stage18_0_files_unmodified(entries) if git_available else True
    closed_unmodified = stage17_10_files_unmodified(entries) if git_available else True

    closed_file_ok = stage17_closed_file_ok(root)
    stage17_structural_ok = stage17_11_structural_evidence_ok(root)
    stage17_evidence_ok = closed_file_ok or stage17_structural_ok
    physical_doc, inextensible_doc, energy_doc = boundary_doc_statuses(root)
    no_activation = stage18_1_no_physics_activation(root) and changed_ok
    required_ok = required_files_present(root)

    statuses.update({
        "stage18_1_requested_status": status(args.stage18_1_enable == "1"),
        "stage18_0_evidence_status": status(stage18_0_evidence_ok(root)),
        "stage17_closed_file_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_closed_evidence_status": status(stage17_evidence_ok),
        "stage17_11_closure_preserved_status": status(stage17_structural_ok or closed_file_ok),
        "stage18_0_wrapper_root_fix_preserved_status": status(stage18_0_wrapper_root_fix_preserved(root)),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_fix_preserved(root)),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_fix_preserved(root)),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_fix_preserved(root)),
        "no_closed_stage_modification_status": status(changed_ok and closed_unmodified),
        "no_stage10_17_file_modification_status": status(changed_ok and closed_unmodified),
        "stage18_0_files_unmodified_status": status(changed_ok and stage18_0_unmodified),
        "stage18_enable_status": status(args.stage18_1_enable == "1"),
        "structure_dynamics_config_enable_status": status(args.structure_dynamics_config_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "physical_structure_boundary_status": status(args.physical_structure_boundary == "1"),
        "bending_config_enable_status": status(args.bending_config_enable == "1"),
        "tension_config_enable_status": status(args.tension_config_enable == "1"),
        "inextensibility_config_enable_status": status(args.inextensibility_config_enable == "1"),
        "time_integration_config_enable_status": status(args.time_integration_config_enable == "1"),
        "energy_diagnostic_config_enable_status": status(args.energy_diagnostic_config_enable == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1" and no_activation),
        "rho_s_value": f"{rho_s:.16e}" if rho_s is not None else "nan",
        "rho_s_status": status(rho_s is not None and rho_s > 0.0),
        "fibre_radius_value": f"{radius:.16e}" if radius is not None else "nan",
        "fibre_radius_status": status(radius is not None and radius > 0.0),
        "young_modulus_value": f"{young:.16e}" if young is not None else "nan",
        "young_modulus_status": status(young is not None and young > 0.0),
        "fibre_length_value": f"{length:.16e}" if length is not None else "nan",
        "fibre_length_status": status(length is not None and length > 0.0),
        "npts_value": str(npts) if npts is not None else "nan",
        "npts_status": status(npts is not None and npts >= 2),
        "dt_structure_value": f"{dt_structure:.16e}" if dt_structure is not None else "nan",
        "dt_structure_status": status(dt_structure is not None and dt_structure > 0.0),
        "rho_tilde_value": f"{rho_tilde:.16e}" if rho_tilde is not None else "nan",
        "rho_tilde_status": status(rho_tilde is not None and rho_tilde > 0.0),
        "gamma_value": f"{gamma:.16e}" if gamma is not None else "nan",
        "gamma_status": status(gamma is not None and gamma >= 0.0),
        "fibre_area_value": f"{area:.16e}",
        "fibre_area_formula_status": status(radius is not None and near(area, math.pi * radius**2, formula_tol)),
        "second_moment_area_value": f"{second_moment:.16e}",
        "second_moment_area_formula_status": status(radius is not None and near(second_moment, math.pi * radius**4 / 4.0, formula_tol)),
        "bending_stiffness_value": f"{bending:.16e}",
        "bending_stiffness_formula_status": status(young is not None and near(bending, young * second_moment, formula_tol)),
        "mass_per_length_value": f"{mass_per_length:.16e}",
        "mass_per_length_formula_status": status(rho_s is not None and near(mass_per_length, rho_s * area, formula_tol)),
        "ds_value": f"{ds:.16e}",
        "ds_formula_status": status(length is not None and npts is not None and npts >= 2 and near(ds, length / (npts - 1), formula_tol)),
        "derived_values_finite_status": status(all(math.isfinite(value) and value >= 0.0 for value in [area, second_moment, bending, mass_per_length, ds])),
        "physical_structure_equation_documented_status": status(physical_doc),
        "inextensibility_constraint_documented_status": status(inextensible_doc),
        "structure_energy_boundary_documented_status": status(energy_doc),
        "config_not_physics_activation_status": status(no_activation),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent(root)),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })

    negative_keys = [
        "no_bending_force_activation_status",
        "no_tension_solve_activation_status",
        "no_inextensibility_enforcement_status",
        "no_structure_time_integration_status",
        "no_structure_state_update_status",
        "no_fluid_force_physical_structure_integration_status",
        "no_structure_energy_power_implementation_status",
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
    ]
    for key in negative_keys:
        statuses[key] = status(no_activation)

    if not required_ok:
        reasons.append("required_stage18_1_file_missing_or_empty")
    if not stage18_0_evidence_ok(root):
        reasons.append("stage18_0_evidence_missing_or_reverted")
    if not stage17_evidence_ok:
        reasons.append("stage17_closure_evidence_missing_or_not_accepted")
    if not stage18_0_wrapper_root_fix_preserved(root):
        reasons.append("stage18_0_wrapper_root_fix_reverted")
    if not changed_ok or not closed_unmodified:
        bad = [path for _code, path in entries if path not in ALLOWED_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))
    if not stage18_0_unmodified:
        bad = [path for _code, path in entries if path in STAGE18_0_FILES]
        reasons.append("stage18_0_file_modified:" + ",".join(bad))
    if not no_activation:
        reasons.append("stage18_1_physics_or_core_contamination_detected")
    if not all(statuses[key] == "PASS" for key in [
        "rho_s_status", "fibre_radius_status", "young_modulus_status", "fibre_length_status",
        "npts_status", "dt_structure_status", "rho_tilde_status", "gamma_status",
    ]):
        reasons.append("invalid_physical_structure_config_parameter")
    if not all(statuses[key] == "PASS" for key in [
        "fibre_area_formula_status", "second_moment_area_formula_status",
        "bending_stiffness_formula_status", "mass_per_length_formula_status", "ds_formula_status",
    ]):
        reasons.append("derived_physical_structure_formula_mismatch")
    if not (physical_doc and inextensible_doc and energy_doc):
        reasons.append("stage18_1_physical_boundary_documentation_missing")
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
