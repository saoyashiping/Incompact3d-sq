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
        print(f"STAGE 3.7 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_force_buffer_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_int(data, "buffer_allocated_flag") == 1,
            f"buffer_allocated_flag={get_int(data, 'buffer_allocated_flag')}")
    require(get_int(data, "buffer_nx") == 16,
            f"buffer_nx={get_int(data, 'buffer_nx')}")
    require(get_int(data, "buffer_ny") == 12,
            f"buffer_ny={get_int(data, 'buffer_ny')}")
    require(get_int(data, "buffer_nz") == 10,
            f"buffer_nz={get_int(data, 'buffer_nz')}")
    require(get_float(data, "buffer_initial_max_abs") <= 1e-14,
            f"buffer_initial_max_abs={get_float(data, 'buffer_initial_max_abs')}")
    require(get_float(data, "buffer_initial_l2_norm") <= 1e-14,
            f"buffer_initial_l2_norm={get_float(data, 'buffer_initial_l2_norm')}")

    require(get_float(data, "single_accumulate_total_force_error") <= 1e-12,
            f"single_accumulate_total_force_error={get_float(data, 'single_accumulate_total_force_error')}")
    require(get_float(data, "single_accumulate_relative_force_error") <= 1e-12,
            f"single_accumulate_relative_force_error={get_float(data, 'single_accumulate_relative_force_error')}")
    require(get_float(data, "single_accumulate_buffer_max_abs") > 0.0,
            f"single_accumulate_buffer_max_abs={get_float(data, 'single_accumulate_buffer_max_abs')}")
    require(get_float(data, "single_accumulate_buffer_l2_norm") > 0.0,
            f"single_accumulate_buffer_l2_norm={get_float(data, 'single_accumulate_buffer_l2_norm')}")

    require(get_float(data, "clear_after_accumulate_max_abs") <= 1e-14,
            f"clear_after_accumulate_max_abs={get_float(data, 'clear_after_accumulate_max_abs')}")
    require(get_float(data, "clear_after_accumulate_l2_norm") <= 1e-14,
            f"clear_after_accumulate_l2_norm={get_float(data, 'clear_after_accumulate_l2_norm')}")

    require(get_float(data, "double_accumulate_total_force_error") <= 1e-12,
            f"double_accumulate_total_force_error={get_float(data, 'double_accumulate_total_force_error')}")
    require(get_float(data, "double_accumulate_relative_force_error") <= 1e-12,
            f"double_accumulate_relative_force_error={get_float(data, 'double_accumulate_relative_force_error')}")

    require(get_float(data, "buffer_superposition_max_error") <= 1e-12,
            f"buffer_superposition_max_error={get_float(data, 'buffer_superposition_max_error')}")
    require(get_float(data, "buffer_superposition_l2_error") <= 1e-12,
            f"buffer_superposition_l2_error={get_float(data, 'buffer_superposition_l2_error')}")

    require(get_float(data, "feedback_buffer_total_force_error") <= 1e-12,
            f"feedback_buffer_total_force_error={get_float(data, 'feedback_buffer_total_force_error')}")
    require(get_float(data, "feedback_buffer_relative_force_error") <= 1e-12,
            f"feedback_buffer_relative_force_error={get_float(data, 'feedback_buffer_relative_force_error')}")
    require(get_float(data, "feedback_buffer_max_abs") > 0.0,
            f"feedback_buffer_max_abs={get_float(data, 'feedback_buffer_max_abs')}")

    require(get_int(data, "force_buffer_status") == 1,
            f"force_buffer_status={get_int(data, 'force_buffer_status')}")

    print("STAGE 3.7 FORCE BUFFER CHECK PASSED")


if __name__ == "__main__":
    main()
