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
        print(f"STAGE 3.2 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_interpolation_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "scalar_constant_max_error") <= 1e-12,
            f"scalar_constant_max_error={get_float(data, 'scalar_constant_max_error')}")
    require(get_float(data, "scalar_constant_weight_sum_max_error") <= 1e-12,
            f"scalar_constant_weight_sum_max_error={get_float(data, 'scalar_constant_weight_sum_max_error')}")
    require(get_float(data, "vector_constant_max_error") <= 1e-12,
            f"vector_constant_max_error={get_float(data, 'vector_constant_max_error')}")
    require(get_float(data, "scalar_linear_inner_max_error") <= 1e-11,
            f"scalar_linear_inner_max_error={get_float(data, 'scalar_linear_inner_max_error')}")
    require(get_float(data, "vector_linear_inner_max_error") <= 1e-11,
            f"vector_linear_inner_max_error={get_float(data, 'vector_linear_inner_max_error')}")
    require(get_float(data, "periodic_wrap_constant_error") <= 1e-12,
            f"periodic_wrap_constant_error={get_float(data, 'periodic_wrap_constant_error')}")
    # This test uses a coarse synthetic grid and a smoothed Peskin kernel.
    # The periodic-field check is intended to catch wrap discontinuities,
    # not to enforce high-order interpolation accuracy for sinusoidal fields.
    require(get_float(data, "periodic_wrap_periodic_field_error_p1") <= 0.15,
            f"periodic_wrap_periodic_field_error_p1={get_float(data, 'periodic_wrap_periodic_field_error_p1')}")
    require(get_float(data, "periodic_wrap_periodic_field_error_p2") <= 0.15,
            f"periodic_wrap_periodic_field_error_p2={get_float(data, 'periodic_wrap_periodic_field_error_p2')}")
    require(abs(get_float(data, "interpolation_weight_sum_min") - 1.0) <= 1e-12,
            f"interpolation_weight_sum_min={get_float(data, 'interpolation_weight_sum_min')}")
    require(abs(get_float(data, "interpolation_weight_sum_max") - 1.0) <= 1e-12,
            f"interpolation_weight_sum_max={get_float(data, 'interpolation_weight_sum_max')}")
    require(get_int(data, "interpolation_status") == 1,
            f"interpolation_status={get_int(data, 'interpolation_status')}")

    print("STAGE 3.2 INTERPOLATION CHECK PASSED")


if __name__ == "__main__":
    main()
