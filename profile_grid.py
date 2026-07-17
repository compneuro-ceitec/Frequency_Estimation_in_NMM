#!/usr/bin/env python3
"""
Line-profile the grid search on a small high-Md subset rich in spurious points.

Requires line_profiler:
    pip install line_profiler

Run with:
    kernprof -l -v profile_grid.py

The script wraps the main grid functions with LineProfiler, runs a small
Md ∈ [3.8, 4.2] × Vt ∈ [-38, -36] grid, and prints per-line timing.
"""

from line_profiler import LineProfiler
import gast_hopf_grid
import sweep_Md


def main():
    profile = LineProfiler()

    # Wrap the functions we want to profile
    gast_hopf_grid.run_grid = profile(gast_hopf_grid.run_grid)
    gast_hopf_grid.iext_roots_for_Md_Vt = profile(gast_hopf_grid.iext_roots_for_Md_Vt)
    sweep_Md.solve_r_for_Md_Vt_Iext = profile(sweep_Md.solve_r_for_Md_Vt_Iext)
    sweep_Md.equilibrium_state = profile(sweep_Md.equilibrium_state)
    sweep_Md.dVdt_residual = profile(sweep_Md.dVdt_residual)
    sweep_Md.eigenfrequencies = profile(sweep_Md.eigenfrequencies)

    # Small high-Md subset where many spurious roots occur
    records = gast_hopf_grid.run_grid(
        md_min=3.8, md_max=4.2, n_md=11,
        vt_min=-38.0, vt_max=-36.0, n_vt=6,
    )
    print(f"\nFound {len(records)} records")

    profile.print_stats()


if __name__ == '__main__':
    main()
