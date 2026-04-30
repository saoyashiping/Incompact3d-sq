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
        print(f"STAGE 3.8 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_boundary_safety_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_int(data, "interior_safe_count") == 3,
            f"interior_safe_count={get_int(data, 'interior_safe_count')}")
    require(get_int(data, "interior_periodic_wrap_count") == 0,
            f"interior_periodic_wrap_count={get_int(data, 'interior_periodic_wrap_count')}")
    require(get_int(data, "interior_unsafe_count") == 0,
            f"interior_unsafe_count={get_int(data, 'interior_unsafe_count')}")
    require(get_int(data, "interior_outside_count") == 0,
            f"interior_outside_count={get_int(data, 'interior_outside_count')}")

    require(get_int(data, "periodic_safe_count") == 0,
            f"periodic_safe_count={get_int(data, 'periodic_safe_count')}")
    require(get_int(data, "periodic_wrap_count") == 4,
            f"periodic_wrap_count={get_int(data, 'periodic_wrap_count')}")
    require(get_int(data, "periodic_unsafe_count") == 0,
            f"periodic_unsafe_count={get_int(data, 'periodic_unsafe_count')}")
    require(get_int(data, "periodic_outside_count") == 0,
            f"periodic_outside_count={get_int(data, 'periodic_outside_count')}")

    require(get_int(data, "wall_safe_count") == 0,
            f"wall_safe_count={get_int(data, 'wall_safe_count')}")
    require(get_int(data, "wall_periodic_wrap_count") == 0,
            f"wall_periodic_wrap_count={get_int(data, 'wall_periodic_wrap_count')}")
    require(get_int(data, "wall_unsafe_count") == 4,
            f"wall_unsafe_count={get_int(data, 'wall_unsafe_count')}")
    require(get_int(data, "wall_outside_count") == 0,
            f"wall_outside_count={get_int(data, 'wall_outside_count')}")
    require(get_float(data, "wall_min_wall_distance_min") >= 0.0,
            f"wall_min_wall_distance_min={get_float(data, 'wall_min_wall_distance_min')}")
    require(get_float(data, "wall_min_wall_distance_max") < 2.0 * (2.0 / 12.0) + 1e-12,
            f"wall_min_wall_distance_max={get_float(data, 'wall_min_wall_distance_max')}")

    require(get_int(data, "outside_safe_count") == 0,
            f"outside_safe_count={get_int(data, 'outside_safe_count')}")
    require(get_int(data, "outside_periodic_wrap_count") == 0,
            f"outside_periodic_wrap_count={get_int(data, 'outside_periodic_wrap_count')}")
    require(get_int(data, "outside_unsafe_count") == 0,
            f"outside_unsafe_count={get_int(data, 'outside_unsafe_count')}")
    require(get_int(data, "outside_outside_count") == 2,
            f"outside_outside_count={get_int(data, 'outside_outside_count')}")

    require(get_int(data, "mixed_safe_count") == 1,
            f"mixed_safe_count={get_int(data, 'mixed_safe_count')}")
    require(get_int(data, "mixed_periodic_wrap_count") == 1,
            f"mixed_periodic_wrap_count={get_int(data, 'mixed_periodic_wrap_count')}")
    require(get_int(data, "mixed_unsafe_count") == 1,
            f"mixed_unsafe_count={get_int(data, 'mixed_unsafe_count')}")
    require(get_int(data, "mixed_outside_count") == 1,
            f"mixed_outside_count={get_int(data, 'mixed_outside_count')}")
    require(get_int(data, "mixed_total_unsafe_count") == 2,
            f"mixed_total_unsafe_count={get_int(data, 'mixed_total_unsafe_count')}")

    require(get_int(data, "boundary_safety_policy_status") == 1,
            f"boundary_safety_policy_status={get_int(data, 'boundary_safety_policy_status')}")
    require(get_int(data, "boundary_safety_status") == 1,
            f"boundary_safety_status={get_int(data, 'boundary_safety_status')}")

    print("STAGE 3.8 BOUNDARY SAFETY CHECK PASSED")


if __name__ == "__main__":
    main()
