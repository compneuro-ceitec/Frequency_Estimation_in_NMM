#!/usr/bin/env python3
"""
Plot 1-D cuts of Iext versus Md for fixed Vt values.

Reads a grid result .npz file, optionally filters true Hopf points, and produces
a multi-panel figure with one Iext-vs-Md curve per fixed Vt slice.

Usage examples
--------------
# Default Vt cuts for the low-Md grid
python3 plot_Iext_vs_Md_cuts.py grid_Md0_1/gast_hopf_grid_results.npz

# True-Hopf only, custom Vt cuts
python3 plot_Iext_vs_Md_cuts.py fine_grid_Md2.5_5/gast_hopf_grid_results.npz \
    --true-hopf --vt-cuts -42 -41 -40 -39 -38 -37 -36 -o Iext_cuts_highMd.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def min_abs_re(eigvals):
    """Return min |Re| among complex-conjugate eigenvalue pairs."""
    pairs = [abs(e.real) for e in eigvals if abs(e.imag) > 1e-10]
    return min(pairs) if pairs else np.nan


def plot_cuts(path, vt_cuts=None, true_hopf=False, hopf_tol=1e-3,
              iext_max=300.0, ncols=4, output='Iext_vs_Md_cuts.png'):
    """Create Iext-vs-Md cuts and save the figure."""
    data = np.load(path, allow_pickle=True)
    Md = data['Md']
    Vt = data['Vt']
    Iext = data['Iext']
    eigvals = data['eigvals']

    if true_hopf:
        min_re = np.array([min_abs_re(ev) for ev in eigvals])
        mask = min_re < hopf_tol
        Md = Md[mask]
        Vt = Vt[mask]
        Iext = Iext[mask]

    if vt_cuts is None:
        # Default: evenly spaced cuts across the Vt range
        vt_min = float(np.ceil(Vt.min()))
        vt_max = float(np.floor(Vt.max()))
        vt_cuts = np.linspace(vt_min, vt_max, 12)

    vt_cuts = np.asarray(vt_cuts)
    n_panels = len(vt_cuts)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    title_suffix = ' (true Hopf)' if true_hopf else ''
    fig.suptitle(f'Iext vs Md for fixed Vt{title_suffix}\n{Path(path).parent}')

    for i, vt0 in enumerate(vt_cuts):
        ax = axes[i]
        mask = np.abs(Vt - vt0) < 0.2
        if mask.any():
            ax.plot(Md[mask], Iext[mask], '.', ms=3)
            ax.set_title(f'Vt ≈ {vt0:.1f}')
        else:
            ax.set_title(f'Vt ≈ {vt0:.1f} (no data)')
        ax.set_ylabel('Iext')
        ax.set_ylim(0, iext_max)
        ax.grid(True, ls=':', alpha=0.5)

    for ax in axes[-ncols:]:
        ax.set_xlabel('Md')

    # Hide unused panels
    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f"Saved {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Grid result .npz file')
    parser.add_argument('-o', '--output', type=str, default='Iext_vs_Md_cuts.png',
                        help='Output PNG filename')
    parser.add_argument('--vt-cuts', type=float, nargs='+', default=None,
                        help='List of fixed Vt values (default: 12 evenly spaced)')
    parser.add_argument('--true-hopf', action='store_true',
                        help='Keep only points where |Re(lambda)| < tol')
    parser.add_argument('--hopf-tol', type=float, default=1e-3,
                        help='Tolerance for true-Hopf filtering (default: 1e-3)')
    parser.add_argument('--iext-max', type=float, default=300.0,
                        help='Upper y-axis limit for Iext (default: 300)')
    parser.add_argument('--ncols', type=int, default=4,
                        help='Number of columns in the figure (default: 4)')
    args = parser.parse_args()

    if not Path(args.input).exists():
        raise FileNotFoundError(args.input)

    plot_cuts(
        path=args.input,
        vt_cuts=args.vt_cuts,
        true_hopf=args.true_hopf,
        hopf_tol=args.hopf_tol,
        iext_max=args.iext_max,
        ncols=args.ncols,
        output=args.output,
    )


if __name__ == '__main__':
    main()
