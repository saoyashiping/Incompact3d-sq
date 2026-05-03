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
        print(f"STAGE 3.3 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_spreading_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "zero_force_grid_norm") <= 1e-14,
            f"zero_force_grid_norm={get_float(data, 'zero_force_grid_norm')}")
    require(get_float(data, "zero_force_total_force_norm") <= 1e-14,
            f"zero_force_total_force_norm={get_float(data, 'zero_force_total_force_norm')}")

    require(get_float(data, "single_lag_total_force_error_norm") <= 1e-12,
            f"single_lag_total_force_error_norm={get_float(data, 'single_lag_total_force_error_norm')}")
    require(get_float(data, "single_lag_relative_force_error") <= 1e-12,
            f"single_lag_relative_force_error={get_float(data, 'single_lag_relative_force_error')}")
    require(get_int(data, "single_lag_nonzero_grid_count") > 0,
            f"single_lag_nonzero_grid_count={get_int(data, 'single_lag_nonzero_grid_count')}")
    require(get_int(data, "single_lag_nonzero_grid_count") <= 64,
            f"single_lag_nonzero_grid_count={get_int(data, 'single_lag_nonzero_grid_count')}")
    require(get_float(data, "single_lag_max_abs_force_density") > 0.0,
            f"single_lag_max_abs_force_density={get_float(data, 'single_lag_max_abs_force_density')}")

    require(get_float(data, "multi_lag_total_force_error_norm") <= 1e-12,
            f"multi_lag_total_force_error_norm={get_float(data, 'multi_lag_total_force_error_norm')}")
    require(get_float(data, "multi_lag_relative_force_error") <= 1e-12,
            f"multi_lag_relative_force_error={get_float(data, 'multi_lag_relative_force_error')}")

    require(get_float(data, "superposition_max_error") <= 1e-12,
            f"superposition_max_error={get_float(data, 'superposition_max_error')}")
    require(get_float(data, "superposition_l2_error") <= 1e-12,
            f"superposition_l2_error={get_float(data, 'superposition_l2_error')}")

    require(get_float(data, "periodic_lag_total_force_error_norm") <= 1e-12,
            f"periodic_lag_total_force_error_norm={get_float(data, 'periodic_lag_total_force_error_norm')}")
    require(get_float(data, "periodic_lag_relative_force_error") <= 1e-12,
            f"periodic_lag_relative_force_error={get_float(data, 'periodic_lag_relative_force_error')}")

    require(get_float(data, "clear_grid_force_norm_after_clear") <= 1e-14,
            f"clear_grid_force_norm_after_clear={get_float(data, 'clear_grid_force_norm_after_clear')}")

    require(get_int(data, "spreading_status") == 1,
            f"spreading_status={get_int(data, 'spreading_status')}")

    print("STAGE 3.3 SPREADING CHECK PASSED")


if __name__ == "__main__":
    main()
