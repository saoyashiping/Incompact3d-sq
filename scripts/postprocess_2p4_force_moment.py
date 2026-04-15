#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path


def periodic_delta(delta: float, period: float) -> float:
    if period <= 0.0:
        return delta
    if delta > 0.5 * period:
        delta -= period
    if delta < -0.5 * period:
        delta += period
    return delta


def read_table(path: Path):
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split()
        rows = []
        for line in f:
            if not line.strip():
                continue
            vals = [float(x) for x in line.split()]
            rows.append(dict(zip(header, vals)))
    return header, rows


def parse_input_center_and_domain(input_file: Path):
    xc = [4.0, 1.0, 2.0]
    xlx, zlz = 8.0, 4.0
    with input_file.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("fiber_center"):
                vals = s.split("=")[1].split("!")[0].split(",")
                xc = [float(v) for v in vals[:3]]
            elif s.startswith("xlx"):
                xlx = float(s.split("=")[1].split("!")[0])
            elif s.startswith("zlz"):
                zlz = float(s.split("=")[1].split("!")[0])
    return xc, xlx, zlz


def compute_weights(first_pts):
    n = len(first_pts)
    if n < 2:
        return [1.0]
    ds = []
    for i in range(n - 1):
        dx = first_pts[i + 1]["x"] - first_pts[i]["x"]
        dy = first_pts[i + 1]["y"] - first_pts[i]["y"]
        dz = first_pts[i + 1]["z"] - first_pts[i]["z"]
        ds.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    ds_mean = sum(ds) / len(ds)
    w = [ds_mean] * n
    w[0] = 0.5 * ds_mean
    w[-1] = 0.5 * ds_mean
    return w


def process_case(report_dir: Path, force_out: Path, slip_out: Path):
    summary_file = report_dir / "fiber_rigid_coupling_summary.dat"
    points_file = report_dir / "fiber_rigid_coupling_points.dat"
    input_file = report_dir / "input.i3d"
    if not summary_file.exists() or not points_file.exists() or not input_file.exists():
        raise FileNotFoundError(f"missing required files in {report_dir}")

    _, sum_rows = read_table(summary_file)
    _, pts_rows = read_table(points_file)
    xc, xlx, zlz = parse_input_center_and_domain(input_file)

    by_time = {}
    for row in pts_rows:
        by_time.setdefault(int(row["itime"]), []).append(row)
    for it in by_time:
        by_time[it].sort(key=lambda r: int(r["index"]))

    first_it = sorted(by_time.keys())[0]
    weights = compute_weights(by_time[first_it])

    with force_out.open("w", encoding="utf-8") as f_force, slip_out.open("w", encoding="utf-8") as f_slip:
        f_force.write("itime time total_Fx total_Fy total_Fz total_Mx total_My total_Mz\n")
        f_slip.write("itime time slip_max slip_rms\n")

        for srow in sum_rows:
            it = int(srow["itime"])
            t = srow["time"]
            pts = by_time[it]
            mx = my = mz = 0.0
            for p, w in zip(pts, weights):
                rx = periodic_delta(p["x"] - xc[0], xlx)
                ry = p["y"] - xc[1]
                rz = periodic_delta(p["z"] - xc[2], zlz)
                fx = p["Ffs_x"]
                fy = p["Ffs_y"]
                fz = p["Ffs_z"]
                mx += (ry * fz - rz * fy) * w
                my += (rz * fx - rx * fz) * w
                mz += (rx * fy - ry * fx) * w
            f_force.write(
                f"{it:d} {t:.16e} {srow['lag_total_Fx']:.16e} {srow['lag_total_Fy']:.16e} {srow['lag_total_Fz']:.16e} "
                f"{mx:.16e} {my:.16e} {mz:.16e}\n"
            )
            f_slip.write(f"{it:d} {t:.16e} {srow['slip_max']:.16e} {srow['slip_rms']:.16e}\n")


def main():
    run_root = Path("/home/sq/runs24_2p4_round4_codefix")
    cases = {
        "parallel_report": ("force_moment_series_parallel.dat", "slip_error_parallel.dat"),
        "wallnormal_report": ("force_moment_series_wallnormal.dat", "slip_error_wallnormal.dat"),
        "deg45_report": ("force_moment_series_45deg.dat", "slip_error_45deg.dat"),
    }
    for case_dir, outputs in cases.items():
        report_dir = run_root / case_dir
        force_out = run_root / outputs[0]
        slip_out = run_root / outputs[1]
        process_case(report_dir, force_out, slip_out)


if __name__ == "__main__":
    main()
