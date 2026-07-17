#!/usr/bin/env python3
"""
Master runscript for generating the Gast-type Hopf grid-search data.

Runs the 1-D Md sweep and/or the focused 2-D (Md, Vt) grids, saves CSV/NPZ
outputs, and produces the default overview plots.

Usage examples
--------------
# Run everything
python3 run_grid_search.py

# Only the low-Md 2-D grid
python3 run_grid_search.py --only-low-md

# Only the high-Md (Md≈4) 2-D grid
python3 run_grid_search.py --only-high-md

# Only the 1-D sweep at Vt=-43
python3 run_grid_search.py --only-sweep

# Custom resolutions
python3 run_grid_search.py --n-md-low 401 --n-vt-low 101 --n-md-high 401 --n-vt-high 101
"""

import argparse
import time
from pathlib import Path

import gast_hopf_grid as grid
import sweep_Md as sweep


# -----------------------------------------------------------------------------
# Default resolutions
# -----------------------------------------------------------------------------
DEFAULTS = {
    'sweep': {'vt': -43.0, 'md_min': 0.0, 'md_max': 20.0, 'n_points': 401},
    'low': {'md_min': 0.0, 'md_max': 1.0, 'n_md': 301, 'vt_min': -55.0, 'vt_max': -35.0, 'n_vt': 101},
    'high': {'md_min': 2.5, 'md_max': 5.0, 'n_md': 301, 'vt_min': -55.0, 'vt_max': -35.0, 'n_vt': 101},
}


def run_1d_sweep(args):
    """Run the 1-D Md sweep at constant Vt."""
    print("\n=== 1-D Md sweep ===")
    print(f"Using n_jobs={args.n_jobs}")
    t0 = time.time()
    records = sweep.run_sweep(
        vt=args.vt_sweep,
        md_min=args.md_min_sweep,
        md_max=args.md_max_sweep,
        n_points=args.n_points_sweep,
        n_jobs=args.n_jobs,
    )
    print(f"Found {len(records)} records in {time.time()-t0:.1f}s")
    csv_path, npz_path = sweep.save_results(records, args.vt_sweep, args.outdir)
    fig_path = sweep.plot_results(records, args.vt_sweep, args.outdir)
    print(f"  CSV:  {csv_path}")
    print(f"  NPZ:  {npz_path}")
    print(f"  Plot: {fig_path}")


def run_2d_grid(segment, args):
    """Run one 2-D (Md, Vt) grid segment and save/plot results."""
    outdir_map = {'low': 'grid_Md0_1', 'high': 'fine_grid_Md2.5_5'}
    outdir = Path(args.outdir) / outdir_map[segment]
    outdir.mkdir(parents=True, exist_ok=True)

    if segment == 'low':
        md_min, md_max = args.md_min_low, args.md_max_low
        n_md = args.n_md_low
        vt_min, vt_max = args.vt_min_low, args.vt_max_low
        n_vt = args.n_vt_low
    else:  # high
        md_min, md_max = args.md_min_high, args.md_max_high
        n_md = args.n_md_high
        vt_min, vt_max = args.vt_min_high, args.vt_max_high
        n_vt = args.n_vt_high

    print(f"\n=== 2-D grid: Md [{md_min}, {md_max}] x {n_md}, Vt [{vt_min}, {vt_max}] x {n_vt} ===")
    print(f"Using n_jobs={args.n_jobs}")
    t0 = time.time()
    records = grid.run_grid(
        md_min=md_min, md_max=md_max, n_md=n_md,
        vt_min=vt_min, vt_max=vt_max, n_vt=n_vt,
        n_jobs=args.n_jobs,
    )
    print(f"Found {len(records)} records in {time.time()-t0:.1f}s")
    csv_path, npz_path = grid.save_results(records, outdir)
    grid_png, cuts_png = grid.plot_grid(records, outdir)
    print(f"  CSV:  {csv_path}")
    print(f"  NPZ:  {npz_path}")
    print(f"  Plot: {grid_png}")
    print(f"  Cuts: {cuts_png}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=str, default='.',
                        help='Base output directory (default: current directory)')
    parser.add_argument('--n-jobs', type=int, default=1,
                        help='Number of parallel workers for 2-D grids (default: 1; use -1 for all cores)')

    # 1-D sweep options
    parser.add_argument('--only-sweep', action='store_true',
                        help='Run only the 1-D sweep')
    parser.add_argument('--vt-sweep', type=float, default=DEFAULTS['sweep']['vt'])
    parser.add_argument('--md-min-sweep', type=float, default=DEFAULTS['sweep']['md_min'])
    parser.add_argument('--md-max-sweep', type=float, default=DEFAULTS['sweep']['md_max'])
    parser.add_argument('--n-points-sweep', type=int, default=DEFAULTS['sweep']['n_points'])

    # Low-Md grid options (Md [0,1])
    parser.add_argument('--only-low-md', action='store_true',
                        help='Run only the Md [0,1] grid')
    parser.add_argument('--md-min-low', type=float, default=DEFAULTS['low']['md_min'])
    parser.add_argument('--md-max-low', type=float, default=DEFAULTS['low']['md_max'])
    parser.add_argument('--n-md-low', type=int, default=DEFAULTS['low']['n_md'])
    parser.add_argument('--vt-min-low', type=float, default=DEFAULTS['low']['vt_min'])
    parser.add_argument('--vt-max-low', type=float, default=DEFAULTS['low']['vt_max'])
    parser.add_argument('--n-vt-low', type=int, default=DEFAULTS['low']['n_vt'])

    # High-Md grid options (Md [2.5,5])
    parser.add_argument('--only-high-md', action='store_true',
                        help='Run only the Md [2.5,5] grid')
    parser.add_argument('--md-min-high', type=float, default=DEFAULTS['high']['md_min'])
    parser.add_argument('--md-max-high', type=float, default=DEFAULTS['high']['md_max'])
    parser.add_argument('--n-md-high', type=int, default=DEFAULTS['high']['n_md'])
    parser.add_argument('--vt-min-high', type=float, default=DEFAULTS['high']['vt_min'])
    parser.add_argument('--vt-max-high', type=float, default=DEFAULTS['high']['vt_max'])
    parser.add_argument('--n-vt-high', type=int, default=DEFAULTS['high']['n_vt'])

    args = parser.parse_args()

    run_sweep = (not args.only_low_md and not args.only_high_md) or args.only_sweep
    run_low = (not args.only_sweep and not args.only_high_md) or args.only_low_md
    run_high = (not args.only_sweep and not args.only_low_md) or args.only_high_md

    t_start = time.time()

    if run_sweep:
        run_1d_sweep(args)
    if run_low:
        run_2d_grid('low', args)
    if run_high:
        run_2d_grid('high', args)

    print(f"\nTotal elapsed: {time.time()-t_start:.1f}s")


if __name__ == '__main__':
    main()
