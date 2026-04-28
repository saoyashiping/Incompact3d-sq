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
        print(f"STAGE 3.4 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_power_consistency_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(abs(get_float(data, "zero_force_power_eulerian")) <= 1e-14,
            f"zero_force_power_eulerian={get_float(data, 'zero_force_power_eulerian')}")
    require(abs(get_float(data, "zero_force_power_lagrangian")) <= 1e-14,
            f"zero_force_power_lagrangian={get_float(data, 'zero_force_power_lagrangian')}")
    require(get_float(data, "zero_force_power_abs_error") <= 1e-14,
            f"zero_force_power_abs_error={get_float(data, 'zero_force_power_abs_error')}")

    require(get_float(data, "constant_velocity_power_abs_error") <= 1e-12,
            f"constant_velocity_power_abs_error={get_float(data, 'constant_velocity_power_abs_error')}")
    require(get_float(data, "constant_velocity_power_rel_error") <= 1e-12,
            f"constant_velocity_power_rel_error={get_float(data, 'constant_velocity_power_rel_error')}")

    require(get_float(data, "nonuniform_power_abs_error") <= 1e-12,
            f"nonuniform_power_abs_error={get_float(data, 'nonuniform_power_abs_error')}")
    require(get_float(data, "nonuniform_power_rel_error") <= 1e-12,
            f"nonuniform_power_rel_error={get_float(data, 'nonuniform_power_rel_error')}")

    require(get_float(data, "periodic_power_abs_error") <= 1e-12,
            f"periodic_power_abs_error={get_float(data, 'periodic_power_abs_error')}")
    require(get_float(data, "periodic_power_rel_error") <= 1e-12,
            f"periodic_power_rel_error={get_float(data, 'periodic_power_rel_error')}")

    require(get_float(data, "deterministic_power_abs_error") <= 1e-12,
            f"deterministic_power_abs_error={get_float(data, 'deterministic_power_abs_error')}")
    require(get_float(data, "deterministic_power_rel_error") <= 1e-12,
            f"deterministic_power_rel_error={get_float(data, 'deterministic_power_rel_error')}")

    require(get_int(data, "power_consistency_status") == 1,
            f"power_consistency_status={get_int(data, 'power_consistency_status')}")

    print("STAGE 3.4 POWER CONSISTENCY CHECK PASSED")


if __name__ == "__main__":
    main()
