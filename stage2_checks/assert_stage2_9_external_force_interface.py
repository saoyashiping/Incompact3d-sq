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
        print(f"STAGE 2.9 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_external_force_interface_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "zero_force_preservation_error_max") <= 1e-10,
            f"zero_force_preservation_error_max={get_float(data, 'zero_force_preservation_error_max')}")
    require(get_float(data, "zero_force_final_bending_energy") <= 1e-12,
            f"zero_force_final_bending_energy={get_float(data, 'zero_force_final_bending_energy')}")
    require(get_float(data, "zero_force_final_kinetic_energy") <= 1e-12,
            f"zero_force_final_kinetic_energy={get_float(data, 'zero_force_final_kinetic_energy')}")
    require(get_float(data, "zero_force_final_length_error") <= 1e-12,
            f"zero_force_final_length_error={get_float(data, 'zero_force_final_length_error')}")
    require(get_int(data, "zero_force_solver_failure_count") == 0,
            f"zero_force_solver_failure_count={get_int(data, 'zero_force_solver_failure_count')}")

    require(get_float(data, "uniform_force_center_velocity_error_norm") <= 1e-6,
            f"uniform_force_center_velocity_error_norm={get_float(data, 'uniform_force_center_velocity_error_norm')}")
    require(get_float(data, "uniform_force_center_accel_error_norm") <= 1e-6,
            f"uniform_force_center_accel_error_norm={get_float(data, 'uniform_force_center_accel_error_norm')}")
    require(get_float(data, "uniform_force_length_error") <= 1e-10,
            f"uniform_force_length_error={get_float(data, 'uniform_force_length_error')}")
    require(get_float(data, "uniform_force_shape_error_max") <= 1e-10,
            f"uniform_force_shape_error_max={get_float(data, 'uniform_force_shape_error_max')}")
    require(get_int(data, "uniform_force_solver_failure_count") == 0,
            f"uniform_force_solver_failure_count={get_int(data, 'uniform_force_solver_failure_count')}")

    require(get_float(data, "sinusoidal_force_net_force_norm") <= 1e-10,
            f"sinusoidal_force_net_force_norm={get_float(data, 'sinusoidal_force_net_force_norm')}")
    require(get_float(data, "sinusoidal_force_final_bending_energy") > 0.0,
            f"sinusoidal_force_final_bending_energy={get_float(data, 'sinusoidal_force_final_bending_energy')}")
    require(get_float(data, "sinusoidal_force_center_drift_norm") <= 1e-6,
            f"sinusoidal_force_center_drift_norm={get_float(data, 'sinusoidal_force_center_drift_norm')}")
    require(get_float(data, "sinusoidal_force_final_length_error") <= 1e-8,
            f"sinusoidal_force_final_length_error={get_float(data, 'sinusoidal_force_final_length_error')}")
    require(get_float(data, "sinusoidal_force_max_displacement") > 0.0,
            f"sinusoidal_force_max_displacement={get_float(data, 'sinusoidal_force_max_displacement')}")
    require(get_int(data, "sinusoidal_force_solver_failure_count") == 0,
            f"sinusoidal_force_solver_failure_count={get_int(data, 'sinusoidal_force_solver_failure_count')}")

    require(get_float(data, "clear_force_norm_after_clear") <= 1e-14,
            f"clear_force_norm_after_clear={get_float(data, 'clear_force_norm_after_clear')}")

    print("STAGE 2.9 EXTERNAL FORCE INTERFACE CHECK PASSED")


if __name__ == "__main__":
    main()
