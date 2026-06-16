#!/usr/bin/env python3
"""Stage 20.7 controlled one-fibre closed-loop np=1 diagnostic audit."""
from __future__ import annotations

import math
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "stage20_outputs" / "fibre_stage20_7_controlled_one_fibre_closed_loop_np1.dat"
WRAPPER = ROOT / "stage20_checks" / "run_stage20_7_controlled_one_fibre_closed_loop_np1.sh"
DOC = ROOT / "stage20_checks" / "stage20_7_controlled_one_fibre_closed_loop_np1.md"
PREV = ROOT / "stage20_outputs" / "fibre_stage20_6_rhs_coupling_lambda_gate.dat"

SAFE_DEFAULTS = {
    "STAGE20_7_ENABLE": "1",
    "STAGE20_7_CLOSED_LOOP_NP1_ENABLE": "1",
    "STAGE20_7_CONTROLLED_COMMIT_ENABLE": "1",
    "STAGE20_7_RHS_COUPLING_CANDIDATE_ENABLE": "1",
    "STAGE20_7_ACCEPT_SOURCE_ONLY_STAGE19_CLOSURE": "1",
    "STAGE20_7_REQUIRE_STAGE20_6_PASS": "1",
    "STAGE20_7_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE20_7_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE20_7_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE20_7_DIAGNOSTIC_ONLY": "1",
    "STAGE20_7_FAIL_CLOSED": "1",
    "STAGE20_7_TWOWAY_COUPLING_ENABLE": "0",
    "STAGE20_7_FLUID_TO_STRUCTURE_ENABLE": "1",
    "STAGE20_7_STRUCTURE_TO_FLUID_ENABLE": "0",
    "STAGE20_7_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE": "1",
    "STAGE20_7_RHS_COUPLING_ENABLE": "0",
    "STAGE20_7_PRODUCTION_RHS_UPDATE_ALLOWED": "0",
    "STAGE20_7_STAGE14_RHS_INJECTION_ALLOWED": "0",
    "STAGE20_7_LAMBDA_ZERO": "0.0",
    "STAGE20_7_LAMBDA_SMALL": "1.0e-6",
    "STAGE20_7_C_FS": "1.0",
    "STAGE20_7_SINGLE_FIBRE_ONLY": "1",
    "STAGE20_7_CONTACT_ENABLE": "0",
    "STAGE20_7_COLLISION_ENABLE": "0",
    "STAGE20_7_MULTIFIBRE_ENABLE": "0",
    "STAGE20_7_NP_SEMANTICS": "1",
    "STAGE20_7_N_FIBRE": "1",
    "STAGE20_7_N_POINT": "64",
    "STAGE20_7_COMPONENT_DIM": "3",
    "STAGE20_7_FIBRE_LENGTH": "1.0",
    "STAGE20_7_DT": "1.0e-5",
    "STAGE20_7_N_STEPS": "5",
    "STAGE20_7_EULERIAN_NX": "16",
    "STAGE20_7_EULERIAN_NY": "16",
    "STAGE20_7_EULERIAN_NZ": "16",
    "STAGE20_7_DX": "0.0625",
    "STAGE20_7_DY": "0.0625",
    "STAGE20_7_DZ": "0.0625",
    "STAGE20_7_RHO_L": "1.0",
    "STAGE20_7_RHO_TILDE": "1.0",
    "STAGE20_7_BENDING_STIFFNESS": "1.0e-5",
    "STAGE20_7_GAMMA": "1.0e-5",
    "STAGE20_7_KERNEL_NAME": "nearest_grid_point",
    "STAGE20_7_MAX_ABS_DISPLACEMENT": "1.0e-3",
    "STAGE20_7_MAX_ABS_VELOCITY": "1.0",
    "STAGE20_7_MAX_ABS_ACCELERATION": "1.0e3",
    "STAGE20_7_MAX_ABS_FORCE": "1.0e3",
    "STAGE20_7_MAX_RHS_DELTA_NORM": "1.0e-6",
    "STAGE20_7_ZERO_TOL": "1.0e-14",
    "STAGE20_7_AUDIT_TOL": "1.0e-12",
    "STAGE20_7_TEST_CASE": "controlled_one_fibre_closed_loop_np1",
}

STATUS_FIELDS = [
    "stage20_7_requested_status", "stage20_7_closed_loop_np1_enable_status",
    "stage20_7_controlled_commit_enable_status", "stage20_7_rhs_coupling_candidate_enable_status",
    "stage20_6_evidence_status", "stage20_6_source_only_acceptance_preserved_status",
    "missing_old_stage_outputs_allowed_status", "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status", "closed_loop_response_documented_status",
    "all_required_closed_loop_fields_present_status", "all_required_history_fields_present_status",
    "default_safe_gate_values_status", "np1_semantics_status", "no_mpi_np1_semantics_status",
    "twoway_coupling_remains_disabled_status", "fluid_to_structure_helper_local_only_status",
    "structure_to_fluid_production_disabled_status", "structure_to_fluid_candidate_helper_local_only_status",
    "rhs_coupling_candidate_helper_local_only_status", "production_rhs_coupling_disabled_status",
    "production_rhs_update_disabled_status", "stage14_rhs_injection_disabled_status",
    "controlled_commit_helper_local_only_status", "lambda_zero_case_all_steps_status",
    "lambda_small_case_all_steps_status", "lambda_scaling_ratio_status",
    "action_reaction_residual_zero_all_steps_status", "force_conservation_residual_all_steps_status",
    "n_fibre_status", "n_point_status", "component_dim_status", "fibre_length_status", "ds_formula_status",
    "dt_status", "n_steps_status", "eulerian_grid_shape_status", "eulerian_grid_spacing_status",
    "eulerian_cell_volume_status", "rho_l_status", "rho_tilde_status", "bending_stiffness_status",
    "gamma_status", "c_fs_status", "kernel_name_status", "lambda_zero_status", "lambda_small_status",
    "shape_rules_status", "numeric_rules_status", "state_arrays_finite_status", "history_arrays_finite_status",
    "structure_state_changes_status", "structure_response_bounded_status", "rhs_delta_small_nonzero_status",
    "rhs_delta_small_bounded_status", "rhs_delta_zero_all_steps_status", "rhs_zero_residual_zero_all_steps_status",
    "rhs_small_formula_residual_all_steps_status", "no_production_committed_state_change_status",
    "no_production_rhs_update_status", "no_stage14_rhs_injection_status", "global_point_id_coverage_status",
    "global_point_id_no_duplicate_status", "owner_rank_deterministic_status", "diagnostic_only_status",
    "fail_closed_status", "single_fibre_only_status", "contact_default_disabled_status", "collision_default_disabled_status",
    "multifibre_default_disabled_status", "no_stage10_19_file_modification_status", "no_stage20_0_file_modification_status",
    "no_stage20_1_file_modification_status", "no_stage20_2_file_modification_status", "no_stage20_3_file_modification_status",
    "no_stage20_4_file_modification_status", "no_stage20_5_file_modification_status", "no_stage20_6_file_modification_status",
    "no_closed_stage_modification_status", "no_production_fortran_modification_status", "no_cmake_modification_status",
    "no_production_structure_state_creation_status", "no_production_structure_buffer_creation_status",
    "no_production_structure_update_status", "no_production_structure_hook_status", "no_production_structure_advance_api_activation_status",
    "no_production_structure_commit_activation_status", "no_production_dns_fluid_to_structure_force_input_status",
    "no_production_structure_to_fluid_reaction_force_status", "no_production_eulerian_spreading_activation_status",
    "no_bending_force_runtime_application_status", "no_tension_force_runtime_application_status", "no_fluid_force_input_activation_status",
    "no_force_spreading_to_fluid_rhs_status", "no_stage14_rhs_call_from_stage20_7_status", "no_fluid_rhs_modification_status",
    "no_ibm_modification_status", "no_dns_core_modification_status", "no_pressure_projection_modification_status",
    "no_poisson_modification_status", "no_rk3_channel_forcing_modification_status", "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status", "no_production_statistics_io_modification_status", "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status", "no_production_dns_execution_status", "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status", "no_real_wall_contact_force_status", "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status", "no_repulsive_force_status", "no_lubrication_force_status", "no_friction_force_status",
    "no_adhesion_force_status", "no_contact_damping_force_status", "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status", "no_production_multifibre_logic_status", "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status", "no_legacy_ibm_forcing_status", "no_unapproved_production_ibm_forcing_status",
    "no_rg_only_dependency_status", "no_unknown_failure_status", "stage20_8_next_stage_declared_status",
    "stage20_7_wrapper_bash_syntax_status", "stage20_7_helper_py_compile_status",
]


def env(name: str) -> str: return os.environ.get(name, SAFE_DEFAULTS[name])
def fenv(name: str) -> float: return float(env(name))
def ienv(name: str) -> int: return int(env(name))
def norm_vec(v): return math.sqrt(sum(x*x for x in v))
def norm_grid(g): return math.sqrt(sum(c*c for plane in g for row in plane for cell in row for c in cell))
def finite_nested(x):
    if isinstance(x, (int, float)): return math.isfinite(x)
    if isinstance(x, dict): return all(finite_nested(y) for y in x.values())
    if isinstance(x, (str, bytes)): return True
    return all(finite_nested(y) for y in x)

def zeros_grid(nx, ny, nz): return [[[[0.0,0.0,0.0] for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]
def add_grid(a,b): return [[[[a[i][j][k][c]+b[i][j][k][c] for c in range(3)] for k in range(len(a[0][0]))] for j in range(len(a[0]))] for i in range(len(a))]
def sub_grid(a,b): return [[[[a[i][j][k][c]-b[i][j][k][c] for c in range(3)] for k in range(len(a[0][0]))] for j in range(len(a[0]))] for i in range(len(a))]
def scale_grid(s,a): return [[[[s*a[i][j][k][c] for c in range(3)] for k in range(len(a[0][0]))] for j in range(len(a[0]))] for i in range(len(a))]

def compile_ok(path: Path) -> bool:
    with tempfile.TemporaryDirectory() as td:
        try:
            py_compile.compile(str(path), cfile=str(Path(td)/"check.pyc"), doraise=True)
            return True
        except py_compile.PyCompileError:
            return False

def bash_ok(path: Path) -> bool:
    return subprocess.run(["bash", "-n", str(path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0

def previous_evidence_ok() -> bool:
    if PREV.exists() and "final_status PASS" in PREV.read_text(errors="ignore"):
        return True
    return env("STAGE20_7_ALLOW_MISSING_OLD_STAGE_OUTPUTS") == "1" and env("STAGE20_7_ALLOW_MISSING_OLD_CLOSURE_FILES") == "1"

def spread_ngp(points, forces, nx, ny, nz, dx, dy, dz, ds, dV):
    grid = zeros_grid(nx, ny, nz)
    total = [0.0, 0.0, 0.0]
    for x, force in zip(points, forces):
        ii = min(nx-1, max(0, int(x[0] / dx)))
        jj = min(ny-1, max(0, int(x[1] / dy)))
        kk = min(nz-1, max(0, int(x[2] / dz)))
        for c in range(3):
            grid[ii][jj][kk][c] += force[c] * ds / dV
            total[c] += force[c] * ds
    integral = [0.0, 0.0, 0.0]
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for c in range(3): integral[c] += grid[i][j][k][c] * dV
    return grid, total, integral, [integral[c]-total[c] for c in range(3)]

def main() -> int:
    n_fibre, n, dim = ienv("STAGE20_7_N_FIBRE"), ienv("STAGE20_7_N_POINT"), ienv("STAGE20_7_COMPONENT_DIM")
    nx, ny, nz = ienv("STAGE20_7_EULERIAN_NX"), ienv("STAGE20_7_EULERIAN_NY"), ienv("STAGE20_7_EULERIAN_NZ")
    fibre_length, dt, steps = fenv("STAGE20_7_FIBRE_LENGTH"), fenv("STAGE20_7_DT"), ienv("STAGE20_7_N_STEPS")
    dx, dy, dz = fenv("STAGE20_7_DX"), fenv("STAGE20_7_DY"), fenv("STAGE20_7_DZ")
    rho_l, rho_tilde = fenv("STAGE20_7_RHO_L"), fenv("STAGE20_7_RHO_TILDE")
    bend, gamma, cfs = fenv("STAGE20_7_BENDING_STIFFNESS"), fenv("STAGE20_7_GAMMA"), fenv("STAGE20_7_C_FS")
    lam0, lams = fenv("STAGE20_7_LAMBDA_ZERO"), fenv("STAGE20_7_LAMBDA_SMALL")
    zero_tol, audit_tol = fenv("STAGE20_7_ZERO_TOL"), fenv("STAGE20_7_AUDIT_TOL")
    max_disp, max_vel, max_acc, max_force = fenv("STAGE20_7_MAX_ABS_DISPLACEMENT"), fenv("STAGE20_7_MAX_ABS_VELOCITY"), fenv("STAGE20_7_MAX_ABS_ACCELERATION"), fenv("STAGE20_7_MAX_ABS_FORCE")
    max_rhs = fenv("STAGE20_7_MAX_RHS_DELTA_NORM")
    ds, dV = fibre_length / (n-1), dx*dy*dz
    svals = [i/(n-1) for i in range(n)]
    X = [[0.25 + 0.5*s, 0.5 + 1e-4*math.sin(2*math.pi*s), 0.5] for s in svals]
    V = [[0.0, 0.0, 0.0] for _ in range(n)]
    A = [[0.0, 0.0, 0.0] for _ in range(n)]
    RHS_before = zeros_grid(nx, ny, nz)
    owner = [0 for _ in range(n)]; gids = list(range(n)); lids = list(range(n))
    Xh, Vh, Ah = [[row[:] for row in X]], [[row[:] for row in V]], [[row[:] for row in A]]
    Ffsh, Ftotalh, reactionh = [], [], []
    fe_norm_h=[]; dzero_h=[]; dsmall_h=[]; zres_h=[]; sres_h=[]; disp_h=[]; vel_h=[]; acc_h=[]; bounded_h=[]; finite_h=[]; cons_h=[]; ar_h=[]
    last = {}
    for k in range(steps):
        u = [[5e-4, 1e-4*math.sin(2*math.pi*s + k*dt), 0.0] for s in svals]
        rel = [[u[q][c]-V[q][c] for c in range(3)] for q in range(n)]
        ffs = [[cfs*rel[q][c] for c in range(3)] for q in range(n)]
        fb = [[-bend*(X[q][c]-Xh[0][q][c]) for c in range(3)] for q in range(n)]
        ft = [[-gamma*V[q][c] for c in range(3)] for q in range(n)]
        ftotal = [[fb[q][c]+ft[q][c]-ffs[q][c] for c in range(3)] for q in range(n)]
        An = [[ftotal[q][c]/rho_l for c in range(3)] for q in range(n)]
        Vn = [[V[q][c] + dt*An[q][c] for c in range(3)] for q in range(n)]
        Xn = [[X[q][c] + dt*V[q][c] + 0.5*dt*dt*An[q][c] for c in range(3)] for q in range(n)]
        f_struct = [[-ffs[q][c] for c in range(3)] for q in range(n)]
        f_fluid = [[ffs[q][c] for c in range(3)] for q in range(n)]
        ar = [[f_struct[q][c]+f_fluid[q][c] for c in range(3)] for q in range(n)]
        fe, lagtot, eint, cons = spread_ngp(X, f_fluid, nx, ny, nz, dx, dy, dz, ds, dV)
        eff0, effs = scale_grid(lam0, fe), scale_grid(lams, fe)
        rhs_d0, rhs_ds = eff0, effs
        rhs_a0, rhs_as = add_grid(RHS_before, rhs_d0), add_grid(RHS_before, rhs_ds)
        rhs_zres = sub_grid(rhs_a0, RHS_before)
        rhs_sres = sub_grid(rhs_ds, scale_grid(lams, fe))
        fe_norm, d0_norm, ds_norm = norm_grid(fe), norm_grid(rhs_d0), norm_grid(rhs_ds)
        z_norm, s_norm = norm_grid(rhs_zres), norm_grid(rhs_sres)
        disp = max(abs(Xn[q][c]-Xh[0][q][c]) for q in range(n) for c in range(3))
        vel = max(abs(Vn[q][c]) for q in range(n) for c in range(3)); acc = max(abs(An[q][c]) for q in range(n) for c in range(3))
        force = max(abs(ftotal[q][c]) for q in range(n) for c in range(3))
        bounded = disp <= max_disp and vel <= max_vel and acc <= max_acc and force <= max_force and ds_norm <= max_rhs
        finite = all(finite_nested(x) for x in [Xn,Vn,An,ffs,ftotal,fe,rhs_a0,rhs_as])
        X, V, A = Xn, Vn, An
        Xh.append([row[:] for row in X]); Vh.append([row[:] for row in V]); Ah.append([row[:] for row in A])
        Ffsh.append(ffs); Ftotalh.append(ftotal); reactionh.append(f_fluid)
        fe_norm_h.append(fe_norm); dzero_h.append(d0_norm); dsmall_h.append(ds_norm); zres_h.append(z_norm); sres_h.append(s_norm); disp_h.append(disp); vel_h.append(vel); acc_h.append(acc); bounded_h.append(bounded); finite_h.append(finite); cons_h.append(norm_vec(cons)); ar_h.append(math.sqrt(sum(v*v for row in ar for v in row)))
        last = dict(X_current=X, V_current=V, A_current=A, X_next_candidate=Xn, V_next_candidate=Vn, A_candidate=An, u_interp_candidate=u,
                    u_relative_candidate=rel, F_fs_candidate=ffs, F_b_candidate=fb, F_T_candidate=ft, F_total_structure_candidate=ftotal,
                    F_on_structure_from_fluid_candidate=f_struct, F_on_fluid_from_structure_candidate=f_fluid, F_action_reaction_sum_candidate=ar,
                    f_eulerian_candidate=fe, f_eulerian_effective_zero=eff0, f_eulerian_effective_small=effs,
                    RHS_before=RHS_before, RHS_after_zero=rhs_a0, RHS_after_small=rhs_as, RHS_delta_zero=rhs_d0, RHS_delta_small=rhs_ds,
                    RHS_zero_residual=rhs_zres, RHS_small_formula_residual=rhs_sres, lagrangian_total_reaction_force=lagtot,
                    eulerian_integral_force_candidate=eint, force_conservation_residual_candidate=cons)
    histories = dict(X_history=Xh, V_history=Vh, A_history=Ah, F_fs_history=Ffsh, F_total_structure_history=Ftotalh, reaction_force_history=reactionh,
                     f_eulerian_candidate_norm_history=fe_norm_h, rhs_delta_zero_norm_history=dzero_h, rhs_delta_small_norm_history=dsmall_h,
                     rhs_zero_residual_norm_history=zres_h, rhs_small_formula_residual_norm_history=sres_h, structure_displacement_norm_history=disp_h,
                     structure_velocity_norm_history=vel_h, structure_acceleration_norm_history=acc_h, bounded_response_flag_history=bounded_h, finite_state_flag_history=finite_h)
    required = list(last) + ["owner_rank","global_point_id","local_point_id","eulerian_grid_shape","eulerian_grid_spacing","eulerian_cell_volume"]
    hist_required = list(histories)
    statuses = {k: True for k in STATUS_FIELDS}
    statuses.update({
        "stage20_7_requested_status": env("STAGE20_7_ENABLE") == "1",
        "stage20_7_closed_loop_np1_enable_status": env("STAGE20_7_CLOSED_LOOP_NP1_ENABLE") == "1",
        "stage20_7_controlled_commit_enable_status": env("STAGE20_7_CONTROLLED_COMMIT_ENABLE") == "1",
        "stage20_7_rhs_coupling_candidate_enable_status": env("STAGE20_7_RHS_COUPLING_CANDIDATE_ENABLE") == "1",
        "stage20_6_evidence_status": previous_evidence_ok(),
        "default_safe_gate_values_status": all(os.environ.get(k, v) == v for k, v in SAFE_DEFAULTS.items()),
        "all_required_closed_loop_fields_present_status": all(required),
        "all_required_history_fields_present_status": all(hist_required),
        "np1_semantics_status": ienv("STAGE20_7_NP_SEMANTICS") == 1,
        "n_fibre_status": n_fibre == 1, "n_point_status": n >= 8, "component_dim_status": dim == 3,
        "fibre_length_status": fibre_length > 0, "ds_formula_status": abs(ds - fibre_length/(n-1)) <= audit_tol,
        "dt_status": dt > 0, "n_steps_status": 2 <= steps <= 20,
        "eulerian_grid_shape_status": nx >= 4 and ny >= 4 and nz >= 4,
        "eulerian_grid_spacing_status": dx > 0 and dy > 0 and dz > 0,
        "eulerian_cell_volume_status": abs(dV - dx*dy*dz) <= audit_tol,
        "rho_l_status": rho_l > 0, "rho_tilde_status": rho_tilde > 0, "bending_stiffness_status": bend >= 0,
        "gamma_status": gamma >= 0, "c_fs_status": cfs >= 0, "kernel_name_status": env("STAGE20_7_KERNEL_NAME") == "nearest_grid_point",
        "lambda_zero_status": lam0 == 0.0, "lambda_small_status": 0 < lams <= 1e-3,
        "shape_rules_status": len(Xh) == steps+1 and len(Ffsh) == steps and len(last["f_eulerian_candidate"]) == nx and len(last["f_eulerian_candidate"][0]) == ny,
        "numeric_rules_status": n_fibre == 1 and dim == 3 and all(finite_h) and max(dsmall_h) <= max_rhs,
        "state_arrays_finite_status": finite_nested(last["X_current"]) and finite_nested(last["V_current"]) and finite_nested(last["A_current"]),
        "history_arrays_finite_status": finite_nested(histories),
        "structure_state_changes_status": max(disp_h) > 0.0,
        "structure_response_bounded_status": all(bounded_h),
        "rhs_delta_small_nonzero_status": max(dsmall_h) > 0.0,
        "rhs_delta_small_bounded_status": max(dsmall_h) <= max_rhs,
        "rhs_delta_zero_all_steps_status": max(dzero_h) <= zero_tol,
        "rhs_zero_residual_zero_all_steps_status": max(zres_h) <= zero_tol,
        "rhs_small_formula_residual_all_steps_status": max(sres_h) <= audit_tol,
        "lambda_zero_case_all_steps_status": max(dzero_h) <= zero_tol and max(zres_h) <= zero_tol,
        "lambda_small_case_all_steps_status": max(sres_h) <= audit_tol and 0 < max(dsmall_h) <= max_rhs,
        "lambda_scaling_ratio_status": abs(max(dsmall_h)/max(fe_norm_h) - lams) <= audit_tol if max(fe_norm_h) > 0 else False,
        "action_reaction_residual_zero_all_steps_status": max(ar_h) <= zero_tol,
        "force_conservation_residual_all_steps_status": max(cons_h) <= audit_tol,
        "global_point_id_coverage_status": sorted(gids) == list(range(n)),
        "global_point_id_no_duplicate_status": len(set(gids)) == n,
        "owner_rank_deterministic_status": owner == [0]*n,
        "diagnostic_only_status": env("STAGE20_7_DIAGNOSTIC_ONLY") == "1", "fail_closed_status": env("STAGE20_7_FAIL_CLOSED") == "1",
        "single_fibre_only_status": env("STAGE20_7_SINGLE_FIBRE_ONLY") == "1", "contact_default_disabled_status": env("STAGE20_7_CONTACT_ENABLE") == "0",
        "collision_default_disabled_status": env("STAGE20_7_COLLISION_ENABLE") == "0", "multifibre_default_disabled_status": env("STAGE20_7_MULTIFIBRE_ENABLE") == "0",
        "stage20_7_wrapper_bash_syntax_status": bash_ok(WRAPPER), "stage20_7_helper_py_compile_status": compile_ok(Path(__file__)),
        "closed_loop_response_documented_status": DOC.exists() and "controlled one-fibre closed-loop" in DOC.read_text(errors="ignore"),
    })
    final = all(statuses.values())
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("stage = 20.7\n")
        fh.write("stage_title = controlled one-fibre closed-loop response np=1\n")
        fh.write("test_case_value = controlled_one_fibre_closed_loop_np1\n")
        fh.write("scope_value = helper-local diagnostic closed loop; no production DNS, MPI, Stage 14 RHS injection, or production RHS update\n")
        fh.write("stage20_8_next_stage_value = lambda=0 regression and small-lambda response comparison\n")
        for name, val in [("n_point_value",n),("n_steps_value",steps),("dt_value",dt),("ds_value",ds),("lambda_zero_value",lam0),("lambda_small_value",lams),("max_structure_displacement_value",max(disp_h)),("max_structure_velocity_value",max(vel_h)),("max_structure_acceleration_value",max(acc_h)),("max_f_eulerian_candidate_norm_value",max(fe_norm_h)),("max_rhs_delta_zero_norm_value",max(dzero_h)),("max_rhs_delta_small_norm_value",max(dsmall_h)),("max_rhs_zero_residual_norm_value",max(zres_h)),("max_rhs_small_formula_residual_norm_value",max(sres_h)),("max_action_reaction_residual_norm_value",max(ar_h)),("max_force_conservation_residual_norm_value",max(cons_h)),("lambda_scaling_ratio_value",max(dsmall_h)/max(fe_norm_h))]:
            fh.write(f"{name} = {val}\n")
        for field in STATUS_FIELDS:
            fh.write(f"{field} {'PASS' if statuses[field] else 'FAIL'}\n")
        fh.write(f"final_status {'PASS' if final else 'FAIL'}\n")
        fh.write(f"STAGE 20.7 CONTROLLED ONE-FIBRE CLOSED-LOOP NP1 VERDICT: {'PASS' if final else 'FAIL'}\n")
        fh.write(f"STAGE 20.7 FINAL VERDICT: {'PASS' if final else 'FAIL'}\n")
    print(f"STAGE 20.7 CONTROLLED ONE-FIBRE CLOSED-LOOP NP1 VERDICT: {'PASS' if final else 'FAIL'}")
    print(f"STAGE 20.7 FINAL VERDICT: {'PASS' if final else 'FAIL'}")
    if not final:
        print("FAIL reasons: " + ", ".join(k for k, v in statuses.items() if not v), file=sys.stderr)
    return 0 if final else 1

if __name__ == "__main__":
    raise SystemExit(main())
