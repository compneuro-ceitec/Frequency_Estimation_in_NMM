#!/usr/bin/env python3
"""
Render publication-quality static 3-D matplotlib scatter plots of Hopf
frequency in (Md, Vt, Iext) space.

This script mirrors render_3d_frequency.py but produces high-resolution
static raster (PNG) and/or vector (PDF) images suitable for a paper.

For each record the marker position is (Md, Vt, Iext) and the color is the
dominant frequency f_dom.  By default only true Hopf points
(|Re(lambda)| < tol) are shown, and each input grid is saved as a separate
file with its own colorbar.

Usage examples
--------------
# Default: separate files, one per grid
python3 render_3d_frequency_static.py grid_Md0_1/gast_hopf_grid_results.npz \
    fine_grid_Md2.5_5/gast_hopf_grid_results.npz -o figure_3d_frequency_true_hopf

# Single multi-panel figure with shared colorbar
python3 render_3d_frequency_static.py grid_Md0_1/gast_hopf_grid_results.npz \
    fine_grid_Md2.5_5/gast_hopf_grid_results.npz \
    --no-separate -o figure_3d_frequency_true_hopf

# Single combined panel with all grids
python3 render_3d_frequency_static.py grid_Md0_1/gast_hopf_grid_results.npz \
    fine_grid_Md2.5_5/gast_hopf_grid_results.npz \
    grid_Md0_20/gast_hopf_grid_results.npz \
    --combined -o figure_3d_frequency_combined

# Colour by secondary frequency (only points with valid f_sec)
python3 render_3d_frequency_static.py grid_Md0_1/gast_hopf_grid_results.npz \
    -o figure_3d_frequency_true_hopf --secondary

# Custom Vt ticks for the low-Md panel
python3 render_3d_frequency_static.py grid_Md0_1/gast_hopf_grid_results.npz \
    -o figure_3d_frequency_true_hopf --vt-ticks -35 -40 -45 -50

# Custom elevation/azimuth for the 3-D view
python3 render_3d_frequency_static.py grid_Md0_20/gast_hopf_grid_results.npz \
    --elev 20 --azim 135 -o figure_3d_frequency
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# Paper-quality defaults
_DEFAULT_DPI = 600
_DEFAULT_WIDTH_IN = 10.0
_DEFAULT_HEIGHT_IN = 7.0
_DEFAULT_MARKER_SIZE = 2.2
_DEFAULT_CMAP = 'plasma'


def min_abs_re(eigvals):
    """Return min |Re| among complex-conjugate eigenvalue pairs."""
    pairs = [abs(e.real) for e in eigvals if abs(e.imag) > 1e-10]
    return min(pairs) if pairs else np.nan


def load_segment(path, true_hopf=True, hopf_tol=1e-3):
    """Load a grid NPZ and return arrays, optionally filtering true Hopf points."""
    data = np.load(path, allow_pickle=True)
    Md = data['Md']
    Vt = data['Vt']
    Iext = data['Iext']
    r = data['r']
    f_dom = data['f_dom']
    f_sec = data['f_sec']
    eigvals = data['eigvals']

    if true_hopf:
        min_re = np.array([min_abs_re(ev) for ev in eigvals])
        mask = min_re < hopf_tol
        Md = Md[mask]
        Vt = Vt[mask]
        Iext = Iext[mask]
        r = r[mask]
        f_dom = f_dom[mask]
        f_sec = f_sec[mask]

    # Drop records with invalid dominant frequency
    valid = ~np.isnan(f_dom)
    Md = Md[valid]
    Vt = Vt[valid]
    Iext = Iext[valid]
    r = r[valid]
    f_dom = f_dom[valid]
    f_sec = f_sec[valid]

    return Md, Vt, Iext, r, f_dom, f_sec


def _make_scatter(ax, Md, Vt, Iext, f_dom, marker_size, cmap, vmin, vmax,
                  elev, azim, vt_ticks=None, xlabel=True, ylabel=True, zlabel=True):
    """Add a 3-D scatter plot to an existing Axes3D."""
    sc = ax.scatter(
        Md, Vt, Iext,
        c=f_dom,
        s=marker_size,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=1.0,
        edgecolors='none',
        rasterized=True,
    )

    # Publication font sizes (~15% smaller than the previous doubled size)
    ax.set_xlabel(r'$M_d$', labelpad=16, fontsize=17)
    ax.set_ylabel(r'$V_t$ (mV)', labelpad=18, fontsize=17)
    # Horizontal z-label for better readability and to avoid clipping
    ax.set_zlabel(r'$I_{\mathrm{ext}}$', labelpad=16, fontsize=17)
    ax.zaxis.set_rotate_label(False)
    ax.zaxis.label.set_rotation(0)

    # Tick labels
    ax.tick_params(axis='x', labelsize=15, pad=8)
    ax.tick_params(axis='y', labelsize=15, pad=8)
    ax.tick_params(axis='z', labelsize=15, pad=8)

    if vt_ticks is not None:
        ax.set_yticks(vt_ticks)

    ax.view_init(elev=elev, azim=azim)

    # Tight panes
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('0.85')
    ax.yaxis.pane.set_edgecolor('0.85')
    ax.zaxis.pane.set_edgecolor('0.85')
    ax.grid(True, ls=':', lw=0.5, color='0.75')

    return sc


def _build_single_panel_figure(seg, title, secondary=False,
                                marker_size=_DEFAULT_MARKER_SIZE,
                                cmap=_DEFAULT_CMAP,
                                width_in=_DEFAULT_WIDTH_IN,
                                height_in=_DEFAULT_HEIGHT_IN,
                                elev=20, azim=135, vt_ticks=None):
    """Build a single 3-D panel figure with its own colorbar."""
    Md, Vt, Iext, r, f_dom, f_sec = seg

    # Choose the frequency to colour by
    f_plot = f_sec if secondary else f_dom
    freq_label = 'f_{\\mathrm{sec}}' if secondary else 'f_{\\mathrm{dom}}'

    # Drop records with invalid chosen frequency
    valid = ~np.isnan(f_plot)
    Md = Md[valid]
    Vt = Vt[valid]
    Iext = Iext[valid]
    r = r[valid]
    f_plot = f_plot[valid]

    fig = plt.figure(figsize=(width_in, height_in))
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, projection='3d')

    # Manual layout: room for a single colorbar and axis labels
    fig.subplots_adjust(
        left=0.12,
        right=0.78,
        bottom=0.10,
        top=0.86,
    )

    if len(Md) == 0:
        label = "No secondary frequency records" if secondary else "No true Hopf points"
        ax.text2D(0.5, 0.5, label,
                  ha='center', va='center', transform=ax.transAxes)
        return fig

    vmin = float(np.nanmin(f_plot))
    vmax = float(np.nanmax(f_plot))

    sc = _make_scatter(
        ax, Md, Vt, Iext, f_plot,
        marker_size=marker_size,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        elev=elev,
        azim=azim,
        vt_ticks=vt_ticks,
    )

    cbar = fig.colorbar(
        sc,
        ax=ax,
        shrink=0.70,
        aspect=25,
        pad=0.16,
        label=rf'Secondary frequency ${freq_label}$ (Hz)' if secondary
              else rf'Dominant frequency ${freq_label}$ (Hz)',
    )
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(rf'Secondary frequency ${freq_label}$ (Hz)' if secondary
                   else rf'Dominant frequency ${freq_label}$ (Hz)', fontsize=15)

    return fig


def build_figure(paths, titles, combined=False, separate=False,
                 secondary=False,
                 true_hopf=True, hopf_tol=1e-3,
                 marker_size=_DEFAULT_MARKER_SIZE,
                 cmap=_DEFAULT_CMAP,
                 width_in=_DEFAULT_WIDTH_IN,
                 height_in=_DEFAULT_HEIGHT_IN,
                 elev=20, azim=135, vt_ticks=None):
    """Build a matplotlib figure with one or more 3-D subplots.

    If ``separate`` is True, return a list of single-panel figures instead of
    one multi-panel figure.  Each panel has its own colorbar.

    If ``secondary`` is True, colour by the secondary frequency ``f_sec`` and
    keep only points where it is defined.
    """
    segments = [
        load_segment(p, true_hopf=true_hopf, hopf_tol=hopf_tol)
        for p in paths
    ]

    if combined:
        Md = np.concatenate([s[0] for s in segments])
        Vt = np.concatenate([s[1] for s in segments])
        Iext = np.concatenate([s[2] for s in segments])
        r = np.concatenate([s[3] for s in segments])
        f_dom = np.concatenate([s[4] for s in segments])
        f_sec = np.concatenate([s[5] for s in segments])
        segments = [(Md, Vt, Iext, r, f_dom, f_sec)]
        titles = ['All grids']

    if not any(len(s[0]) for s in segments):
        raise RuntimeError("No records to plot after filtering.")

    if separate:
        return [
            _build_single_panel_figure(
                seg, title, secondary=secondary,
                marker_size=marker_size,
                cmap=cmap,
                width_in=width_in,
                height_in=height_in,
                elev=elev,
                azim=azim,
                vt_ticks=vt_ticks,
            )
            for seg, title in zip(segments, titles)
        ]

    # Multi-panel figure: choose frequency and filter each segment
    freq_idx = 5 if secondary else 4
    freq_label = 'f_{\\mathrm{sec}}' if secondary else 'f_{\\mathrm{dom}}'
    filtered_segments = []
    for seg in segments:
        Md, Vt, Iext, r, f_dom, f_sec = seg
        f_plot = f_sec if secondary else f_dom
        valid = ~np.isnan(f_plot)
        filtered_segments.append(
            (Md[valid], Vt[valid], Iext[valid], r[valid],
             f_dom[valid], f_sec[valid], f_plot[valid])
        )

    f_all = np.concatenate([s[-1] for s in filtered_segments])
    if len(f_all) == 0:
        raise RuntimeError("No records with the requested frequency.")
    vmin = float(np.nanmin(f_all))
    vmax = float(np.nanmax(f_all))

    n_cols = len(filtered_segments)
    fig, axes = plt.subplots(
        1, n_cols,
        figsize=(width_in * n_cols + 1.2, height_in + 0.6),
        subplot_kw={'projection': '3d'},
    )
    fig.patch.set_facecolor('white')
    if n_cols == 1:
        axes = [axes]

    # Manual layout: generous room for colorbar and axis labels
    fig.subplots_adjust(
        left=0.10,
        right=0.72,
        bottom=0.10,
        top=0.82,
        wspace=0.75,
    )

    for ax, (seg, title) in zip(axes, zip(filtered_segments, titles)):
        Md, Vt, Iext, r, f_dom, f_sec, f_plot = seg

        if len(Md) == 0:
            label = "No secondary frequency records" if secondary else "No true Hopf points"
            ax.text2D(0.5, 0.5, label,
                      ha='center', va='center', transform=ax.transAxes)
            continue

        sc = _make_scatter(
            ax, Md, Vt, Iext, f_plot,
            marker_size=marker_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            elev=elev,
            azim=azim,
            vt_ticks=vt_ticks,
        )

    # Shared colorbar placed in the right margin
    cbar = fig.colorbar(
        sc,
        ax=axes,
        shrink=0.60,
        aspect=22,
        pad=0.28,
        label=rf'Secondary frequency ${freq_label}$ (Hz)' if secondary
              else rf'Dominant frequency ${freq_label}$ (Hz)',
    )
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(rf'Secondary frequency ${freq_label}$ (Hz)' if secondary
                   else rf'Dominant frequency ${freq_label}$ (Hz)', fontsize=15)

    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'inputs',
        nargs='+',
        help='One or more grid result .npz files',
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='outputs/figure_3d_frequency_true_hopf',
        help='Output filename stem (default: outputs/figure_3d_frequency_true_hopf)',
    )
    parser.add_argument(
        '--combined',
        action='store_true',
        help='Plot all grids in a single 3-D panel instead of subplots',
    )
    parser.add_argument(
        '--separate',
        action='store_true',
        default=True,
        help='Save each grid as a separate file with its own colorbar (default: True)',
    )
    parser.add_argument(
        '--no-separate',
        dest='separate',
        action='store_false',
        help='Plot all grids in a single multi-panel figure',
    )
    parser.add_argument(
        '--true-hopf',
        action='store_true',
        default=True,
        help='Keep only points where the Hopf pair has |Re(lambda)| < tol (default: True)',
    )
    parser.add_argument(
        '--no-true-hopf',
        dest='true_hopf',
        action='store_false',
        help='Show all points, not only true Hopf points',
    )
    parser.add_argument(
        '--hopf-tol',
        type=float,
        default=1e-3,
        help='Tolerance for true-Hopf filtering (default: 1e-3)',
    )
    parser.add_argument(
        '--marker-size',
        type=float,
        default=_DEFAULT_MARKER_SIZE,
        help=f'Marker size (default: {_DEFAULT_MARKER_SIZE})',
    )
    parser.add_argument(
        '--width',
        type=float,
        default=_DEFAULT_WIDTH_IN,
        help=f'Figure width per panel in inches (default: {_DEFAULT_WIDTH_IN})',
    )
    parser.add_argument(
        '--height',
        type=float,
        default=_DEFAULT_HEIGHT_IN,
        help=f'Figure height in inches (default: {_DEFAULT_HEIGHT_IN})',
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=_DEFAULT_DPI,
        help=f'PNG resolution in dpi (default: {_DEFAULT_DPI})',
    )
    parser.add_argument(
        '--cmap',
        type=str,
        default=_DEFAULT_CMAP,
        help=f'Colormap (default: {_DEFAULT_CMAP})',
    )
    parser.add_argument(
        '--elev',
        type=float,
        default=20,
        help='3-D view elevation (default: 20)',
    )
    parser.add_argument(
        '--azim',
        type=float,
        default=120,
        help='3-D view azimuth (default: 120)',
    )
    parser.add_argument(
        '--formats',
        nargs='+',
        choices=['png', 'pdf', 'svg', 'eps'],
        default=['png', 'pdf'],
        help='Output formats (default: png pdf)',
    )
    parser.add_argument(
        '--secondary',
        action='store_true',
        help='Colour by secondary frequency f_sec and keep only points where it is defined',
    )
    parser.add_argument(
        '--vt-ticks',
        type=float,
        nargs='+',
        default=None,
        help='Custom Vt axis tick locations (default: auto)',
    )
    args = parser.parse_args()

    for p in args.inputs:
        if not Path(p).exists():
            raise FileNotFoundError(p)

    titles = [Path(p).parent.name for p in args.inputs]

    figs = build_figure(
        paths=args.inputs,
        titles=titles,
        combined=args.combined,
        separate=args.separate,
        secondary=args.secondary,
        true_hopf=args.true_hopf,
        hopf_tol=args.hopf_tol,
        marker_size=args.marker_size,
        cmap=args.cmap,
        width_in=args.width,
        height_in=args.height,
        elev=args.elev,
        azim=args.azim,
        vt_ticks=args.vt_ticks,
    )

    if not isinstance(figs, list):
        figs = [figs]
        titles = [args.output]

    out_stem = Path(args.output)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    stem = f"{out_stem.stem}_secondary" if args.secondary else out_stem.stem
    for fig, title in zip(figs, titles):
        for fmt in args.formats:
            if args.separate:
                out_path = out_stem.parent / f"{stem}_{title}.{fmt}"
            else:
                out_path = (out_stem.parent / stem).with_suffix(f'.{fmt}')
            if fmt == 'png':
                fig.savefig(out_path, dpi=args.dpi, bbox_inches='tight',
                            pad_inches=0.15, facecolor='white', edgecolor='none')
            else:
                fig.savefig(out_path, bbox_inches='tight', pad_inches=0.15,
                            facecolor='white', edgecolor='none')
            print(f"Saved {out_path}")
        plt.close(fig)


if __name__ == '__main__':
    main()
