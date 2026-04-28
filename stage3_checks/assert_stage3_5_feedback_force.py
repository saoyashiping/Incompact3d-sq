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
        print(f"STAGE 3.5 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_feedback_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_float(data, "zero_slip_force_structure_norm") <= 1e-14,
            f"zero_slip_force_structure_norm={get_float(data, 'zero_slip_force_structure_norm')}")
    require(get_float(data, "zero_slip_force_fluid_norm") <= 1e-14,
            f"zero_slip_force_fluid_norm={get_float(data, 'zero_slip_force_fluid_norm')}")
    require(get_float(data, "zero_slip_max_slip_norm") <= 1e-14,
            f"zero_slip_max_slip_norm={get_float(data, 'zero_slip_max_slip_norm')}")
    require(abs(get_float(data, "zero_slip_power_total")) <= 1e-14,
            f"zero_slip_power_total={get_float(data, 'zero_slip_power_total')}")

    require(get_float(data, "constant_slip_structure_force_error") <= 1e-14,
            f"constant_slip_structure_force_error={get_float(data, 'constant_slip_structure_force_error')}")
    require(get_float(data, "constant_slip_fluid_force_error") <= 1e-14,
            f"constant_slip_fluid_force_error={get_float(data, 'constant_slip_fluid_force_error')}")
    require(get_float(data, "constant_slip_action_reaction_error") <= 1e-14,
            f"constant_slip_action_reaction_error={get_float(data, 'constant_slip_action_reaction_error')}")

    require(get_float(data, "action_reaction_total_force_error") <= 1e-14,
            f"action_reaction_total_force_error={get_float(data, 'action_reaction_total_force_error')}")

    require(get_float(data, "feedback_power_total") <= 1e-14,
            f"feedback_power_total={get_float(data, 'feedback_power_total')}")
    require(get_float(data, "feedback_expected_dissipation") > 0.0,
            f"feedback_expected_dissipation={get_float(data, 'feedback_expected_dissipation')}")
    require(get_float(data, "feedback_dissipation_identity_error") <= 1e-12,
            f"feedback_dissipation_identity_error={get_float(data, 'feedback_dissipation_identity_error')}")

    require(get_float(data, "feedback_spread_total_force_error") <= 1e-12,
            f"feedback_spread_total_force_error={get_float(data, 'feedback_spread_total_force_error')}")
    require(get_float(data, "feedback_spread_relative_force_error") <= 1e-12,
            f"feedback_spread_relative_force_error={get_float(data, 'feedback_spread_relative_force_error')}")

    require(get_float(data, "interpolated_zero_slip_force_norm") <= 1e-14,
            f"interpolated_zero_slip_force_norm={get_float(data, 'interpolated_zero_slip_force_norm')}")

    require(get_int(data, "feedback_status") == 1,
            f"feedback_status={get_int(data, 'feedback_status')}")

    print("STAGE 3.5 FEEDBACK FORCE CHECK PASSED")


if __name__ == "__main__":
    main()
