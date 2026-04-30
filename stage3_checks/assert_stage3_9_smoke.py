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
        print(f"STAGE 3.9 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_stage3_smoke_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_int(data, "smoke_boundary_safe_count") == 33,
            f"smoke_boundary_safe_count={get_int(data, 'smoke_boundary_safe_count')}")
    require(get_int(data, "smoke_boundary_periodic_wrap_count") == 0,
            f"smoke_boundary_periodic_wrap_count={get_int(data, 'smoke_boundary_periodic_wrap_count')}")
    require(get_int(data, "smoke_boundary_unsafe_count") == 0,
            f"smoke_boundary_unsafe_count={get_int(data, 'smoke_boundary_unsafe_count')}")
    require(get_int(data, "smoke_boundary_outside_count") == 0,
            f"smoke_boundary_outside_count={get_int(data, 'smoke_boundary_outside_count')}")
    require(get_int(data, "smoke_boundary_total_unsafe_count") == 0,
            f"smoke_boundary_total_unsafe_count={get_int(data, 'smoke_boundary_total_unsafe_count')}")

    require(get_float(data, "uniform_smoke_final_center_velocity_x") > 0.0,
            f"uniform_smoke_final_center_velocity_x={get_float(data, 'uniform_smoke_final_center_velocity_x')}")
    require(get_int(data, "uniform_smoke_direction_check") == 1,
            f"uniform_smoke_direction_check={get_int(data, 'uniform_smoke_direction_check')}")
    require(get_float(data, "uniform_smoke_length_error") <= 1e-8,
            f"uniform_smoke_length_error={get_float(data, 'uniform_smoke_length_error')}")
    require(get_float(data, "uniform_smoke_shape_error_max") <= 1e-8,
            f"uniform_smoke_shape_error_max={get_float(data, 'uniform_smoke_shape_error_max')}")
    require(get_int(data, "uniform_smoke_solver_failure_count") == 0,
            f"uniform_smoke_solver_failure_count={get_int(data, 'uniform_smoke_solver_failure_count')}")
    require(get_int(data, "uniform_smoke_nan_detected") == 0,
            f"uniform_smoke_nan_detected={get_int(data, 'uniform_smoke_nan_detected')}")

    require(get_float(data, "smoke_force_conservation_error") <= 1e-10,
            f"smoke_force_conservation_error={get_float(data, 'smoke_force_conservation_error')}")
    require(get_float(data, "smoke_force_conservation_relative_error") <= 1e-10,
            f"smoke_force_conservation_relative_error={get_float(data, 'smoke_force_conservation_relative_error')}")
    require(get_float(data, "smoke_force_buffer_max_abs") > 0.0,
            f"smoke_force_buffer_max_abs={get_float(data, 'smoke_force_buffer_max_abs')}")

    require(get_float(data, "smoke_power_abs_error") <= 1e-10,
            f"smoke_power_abs_error={get_float(data, 'smoke_power_abs_error')}")
    require(get_float(data, "smoke_power_relative_error") <= 1e-10,
            f"smoke_power_relative_error={get_float(data, 'smoke_power_relative_error')}")

    require(abs(get_float(data, "smoke_power_abs_error") -
                abs(get_float(data, "smoke_power_eulerian") - get_float(data, "smoke_power_lagrangian"))) <= 1e-12,
            "smoke power abs-error consistency mismatch")
    require(abs(get_float(data, "smoke_power_recomputed_abs_error") -
                abs(get_float(data, "smoke_power_eulerian") - get_float(data, "smoke_power_lagrangian"))) <= 1e-12,
            "smoke_power_recomputed_abs_error mismatch")
    require(get_float(data, "smoke_power_error_consistency_check") <= 1e-12,
            f"smoke_power_error_consistency_check={get_float(data, 'smoke_power_error_consistency_check')}")

    require(get_float(data, "zero_slip_smoke_f_ext_norm") <= 1e-10,
            f"zero_slip_smoke_f_ext_norm={get_float(data, 'zero_slip_smoke_f_ext_norm')}")
    require(get_float(data, "zero_slip_smoke_force_buffer_norm") <= 1e-10,
            f"zero_slip_smoke_force_buffer_norm={get_float(data, 'zero_slip_smoke_force_buffer_norm')}")

    require(get_float(data, "nonuniform_smoke_f_ext_norm") > 0.0,
            f"nonuniform_smoke_f_ext_norm={get_float(data, 'nonuniform_smoke_f_ext_norm')}")
    require(get_float(data, "nonuniform_smoke_force_buffer_norm") > 0.0,
            f"nonuniform_smoke_force_buffer_norm={get_float(data, 'nonuniform_smoke_force_buffer_norm')}")
    require(get_float(data, "nonuniform_smoke_force_conservation_relative_error") <= 1e-10,
            f"nonuniform_smoke_force_conservation_relative_error={get_float(data, 'nonuniform_smoke_force_conservation_relative_error')}")
    require(get_float(data, "nonuniform_smoke_power_relative_error") <= 1e-10,
            f"nonuniform_smoke_power_relative_error={get_float(data, 'nonuniform_smoke_power_relative_error')}")

    require(get_float(data, "nonuniform_smoke_power_error_consistency_check") <= 1e-12,
            f"nonuniform_smoke_power_error_consistency_check={get_float(data, 'nonuniform_smoke_power_error_consistency_check')}")
    require(get_float(data, "nonuniform_smoke_length_error") <= 1e-8,
            f"nonuniform_smoke_length_error={get_float(data, 'nonuniform_smoke_length_error')}")
    require(get_int(data, "nonuniform_smoke_nan_detected") == 0,
            f"nonuniform_smoke_nan_detected={get_int(data, 'nonuniform_smoke_nan_detected')}")
    require(get_int(data, "nonuniform_smoke_unsafe_count") == 0,
            f"nonuniform_smoke_unsafe_count={get_int(data, 'nonuniform_smoke_unsafe_count')}")

    require(get_int(data, "stage3_smoke_status") == 1,
            f"stage3_smoke_status={get_int(data, 'stage3_smoke_status')}")

    print("STAGE 3.9 SMOKE CHECK PASSED")


if __name__ == "__main__":
    main()
