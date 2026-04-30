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
        print(f"STAGE 3.1 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_delta_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "phi_0") > 0.0, f"phi_0={get_float(data, 'phi_0')}")
    require(get_float(data, "phi_05") > 0.0, f"phi_05={get_float(data, 'phi_05')}")
    require(get_float(data, "phi_10") > 0.0, f"phi_10={get_float(data, 'phi_10')}")
    require(get_float(data, "phi_15") > 0.0, f"phi_15={get_float(data, 'phi_15')}")
    require(abs(get_float(data, "phi_20")) <= 1e-15, f"phi_20={get_float(data, 'phi_20')}")

    require(abs(get_float(data, "weight_sum_1d") - 1.0) <= 1e-12,
            f"weight_sum_1d={get_float(data, 'weight_sum_1d')}")
    require(abs(get_float(data, "first_moment_1d")) <= 1e-12,
            f"first_moment_1d={get_float(data, 'first_moment_1d')}")

    require(abs(get_float(data, "weight_sum_3d") - 1.0) <= 1e-12,
            f"weight_sum_3d={get_float(data, 'weight_sum_3d')}")
    require(abs(get_float(data, "first_moment_x_3d")) <= 1e-12,
            f"first_moment_x_3d={get_float(data, 'first_moment_x_3d')}")
    require(abs(get_float(data, "first_moment_y_3d")) <= 1e-12,
            f"first_moment_y_3d={get_float(data, 'first_moment_y_3d')}")
    require(abs(get_float(data, "first_moment_z_3d")) <= 1e-12,
            f"first_moment_z_3d={get_float(data, 'first_moment_z_3d')}")

    require(get_float(data, "constant_field_error") <= 1e-12,
            f"constant_field_error={get_float(data, 'constant_field_error')}")

    require(abs(get_float(data, "periodic_wrap_distance_x") - 0.10) <= 1e-12,
            f"periodic_wrap_distance_x={get_float(data, 'periodic_wrap_distance_x')}")
    require(abs(get_float(data, "nonperiodic_distance_x") + 1.90) <= 1e-12,
            f"nonperiodic_distance_x={get_float(data, 'nonperiodic_distance_x')}")

    require(abs(get_float(data, "phi_outside_positive")) <= 1e-15,
            f"phi_outside_positive={get_float(data, 'phi_outside_positive')}")
    require(abs(get_float(data, "phi_outside_negative")) <= 1e-15,
            f"phi_outside_negative={get_float(data, 'phi_outside_negative')}")

    require(get_int(data, "delta_kernel_status") == 1,
            f"delta_kernel_status={get_int(data, 'delta_kernel_status')}")

    print("STAGE 3.1 DELTA KERNEL CHECK PASSED")


if __name__ == "__main__":
    main()
