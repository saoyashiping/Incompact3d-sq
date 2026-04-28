#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys


def parse_kv(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = value.strip()
    return data


def get_float(data: dict[str, str], key: str) -> float:
    if key not in data:
        raise KeyError(key)
    return float(data[key].split()[0])


def get_int(data: dict[str, str], key: str) -> int:
    if key not in data:
        raise KeyError(key)
    return int(data[key].split()[0])


def require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"STAGE 3.6 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_structure_coupling_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "zero_slip_f_ext_norm") <= 1e-12,
            f"zero_slip_f_ext_norm={get_float(data, 'zero_slip_f_ext_norm')}")
    require(get_float(data, "zero_slip_force_structure_norm") <= 1e-12,
            f"zero_slip_force_structure_norm={get_float(data, 'zero_slip_force_structure_norm')}")
    require(get_float(data, "zero_slip_force_fluid_norm") <= 1e-12,
            f"zero_slip_force_fluid_norm={get_float(data, 'zero_slip_force_fluid_norm')}")

    require(get_float(data, "uniform_drag_final_center_velocity_x") > 0.0,
            f"uniform_drag_final_center_velocity_x={get_float(data, 'uniform_drag_final_center_velocity_x')}")
    require(get_int(data, "uniform_drag_direction_check") == 1,
            f"uniform_drag_direction_check={get_int(data, 'uniform_drag_direction_check')}")
    require(get_float(data, "uniform_drag_center_velocity_error") <= 1e-4,
            f"uniform_drag_center_velocity_error={get_float(data, 'uniform_drag_center_velocity_error')}")
    require(get_float(data, "uniform_drag_shape_error_max") <= 1e-8,
            f"uniform_drag_shape_error_max={get_float(data, 'uniform_drag_shape_error_max')}")
    require(get_float(data, "uniform_drag_length_error") <= 1e-8,
            f"uniform_drag_length_error={get_float(data, 'uniform_drag_length_error')}")
    require(get_int(data, "uniform_drag_solver_failure_count") == 0,
            f"uniform_drag_solver_failure_count={get_int(data, 'uniform_drag_solver_failure_count')}")

    require(get_float(data, "reverse_drag_final_center_velocity_x") < 0.0,
            f"reverse_drag_final_center_velocity_x={get_float(data, 'reverse_drag_final_center_velocity_x')}")
    require(get_int(data, "reverse_drag_direction_check") == 1,
            f"reverse_drag_direction_check={get_int(data, 'reverse_drag_direction_check')}")

    require(get_float(data, "force_set_error") <= 1e-14,
            f"force_set_error={get_float(data, 'force_set_error')}")
    require(get_float(data, "force_add_error") <= 1e-14,
            f"force_add_error={get_float(data, 'force_add_error')}")
    require(get_float(data, "force_clear_error") <= 1e-14,
            f"force_clear_error={get_float(data, 'force_clear_error')}")

    require(get_float(data, "nonuniform_f_ext_norm") > 0.0,
            f"nonuniform_f_ext_norm={get_float(data, 'nonuniform_f_ext_norm')}")
    require(get_float(data, "nonuniform_force_variation_norm") > 0.0,
            f"nonuniform_force_variation_norm={get_float(data, 'nonuniform_force_variation_norm')}")
    require(get_float(data, "nonuniform_force_matches_f_ext_error") <= 1e-14,
            f"nonuniform_force_matches_f_ext_error={get_float(data, 'nonuniform_force_matches_f_ext_error')}")

    require(get_int(data, "ibm_structure_coupling_status") == 1,
            f"ibm_structure_coupling_status={get_int(data, 'ibm_structure_coupling_status')}")

    print("STAGE 3.6 STRUCTURE COUPLING CHECK PASSED")


if __name__ == "__main__":
    main()
