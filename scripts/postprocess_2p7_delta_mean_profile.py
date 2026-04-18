#!/usr/bin/env python3
"""Compute delta mean profile between fiber and single-phase profiles.

Usage:
  python scripts/postprocess_2p7_delta_mean_profile.py --single umean_single.dat --fiber umean_fiber.dat --output delta_umean.dat
"""

from __future__ import annotations
import argparse
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--single", required=True, help="Single-phase mean profile (y value)")
    ap.add_argument("--fiber", required=True, help="Fiber mean profile (y value)")
    ap.add_argument("--output", required=True, help="Output delta profile")
    args = ap.parse_args()

    s = np.loadtxt(args.single)
    f = np.loadtxt(args.fiber)
    if s.ndim == 1:
        s = s.reshape(1, -1)
    if f.ndim == 1:
        f = f.reshape(1, -1)

    ys, us = s[:, 0], s[:, 1]
    yf, uf = f[:, 0], f[:, 1]

    uf_interp = np.interp(ys, yf, uf)
    out = np.column_stack((ys, uf_interp - us))
    np.savetxt(args.output, out, header="y delta_mean_value", comments="")


if __name__ == "__main__":
    main()
