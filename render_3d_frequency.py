#!/usr/bin/env python3
"""
Render interactive 3-D Plotly scatter plots of Hopf frequency in
(Md, Vt, Iext) space.

For each record the marker position is (Md, Vt, Iext) and the color is the
dominant frequency f_dom.  Hover text shows Md, Vt, Iext, r, f_dom and f_sec.

Usage examples
--------------
# All records from one or more grids
python3 render_3d_frequency.py grid_Md0_1/gast_hopf_grid_results.npz fine_grid_Md2.5_5/gast_hopf_grid_results.npz -o interactive_3d_frequency.html

# Only true Hopf points (|Re(lambda)| < tol)
python3 render_3d_frequency.py grid_Md0_1/gast_hopf_grid_results.npz fine_grid_Md2.5_5/gast_hopf_grid_results.npz --true-hopf -o interactive_3d_frequency_true_hopf.html
"""

import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def min_abs_re(eigvals):
    """Return min |Re| among complex-conjugate eigenvalue pairs."""
    pairs = [abs(e.real) for e in eigvals if abs(e.imag) > 1e-10]
    return min(pairs) if pairs else np.nan


def load_segment(path, true_hopf=False, hopf_tol=1e-3):
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
        eigvals = eigvals[mask]
        min_re = min_re[mask]
    else:
        min_re = np.full(len(Md), np.nan)

    return Md, Vt, Iext, r, f_dom, f_sec, min_re


def build_figure(paths, titles, true_hopf=False, hopf_tol=1e-3,
                 marker_size=2, width=1400, height=700):
    """Build a Plotly figure with one 3-D subplot per input grid."""
    fig = make_subplots(
        rows=1,
        cols=len(paths),
        subplot_titles=titles,
        specs=[[{'type': 'scatter3d'}] * len(paths)],
    )

    for col, (path, title) in enumerate(zip(paths, titles), 1):
        Md, Vt, Iext, r, f_dom, f_sec, min_re = load_segment(
            path, true_hopf=true_hopf, hopf_tol=hopf_tol
        )

        if len(Md) == 0:
            print(f"Warning: no records in {path} after filtering.")
            continue

        hover = [
            f"Md={m:.4f}<br>Vt={v:.2f}<br>Iext={ie:.4f}<br>r={rr:.6f}<br>"
            f"f_dom={fd:.3f} Hz<br>f_sec={fs:.3f} Hz"
            + (f"<br>min|Re|={mr:.3e}" if not np.isnan(mr) else "")
            for m, v, ie, rr, fd, fs, mr in zip(
                Md, Vt, Iext, r, f_dom, f_sec, min_re
            )
        ]

        scatter = go.Scatter3d(
            x=Md,
            y=Vt,
            z=Iext,
            mode='markers',
            marker=dict(
                size=marker_size,
                color=f_dom,
                colorscale='Plasma',
                colorbar=dict(
                    title='f_dom (Hz)',
                    x=0.45 if col == 1 and len(paths) > 1 else 1.0,
                ),
                cmin=np.nanmin(f_dom),
                cmax=np.nanmax(f_dom),
                showscale=True,
            ),
            text=hover,
            hoverinfo='text',
            name=title,
        )
        fig.add_trace(scatter, row=1, col=col)

        # Per-scene axis titles
        scene_key = 'scene' if col == 1 else f'scene{col}'
        fig.layout[scene_key].update(
            xaxis_title='Md',
            yaxis_title='Vt',
            zaxis_title='Iext',
        )

    title_text = (
        'Interactive 3-D: true Hopf frequency over Md, Vt, Iext'
        if true_hopf else 'Interactive 3-D: frequency over Md, Vt, Iext'
    )
    fig.update_layout(
        title=title_text,
        width=width,
        height=height,
    )
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
        default='interactive_3d_frequency.html',
        help='Output HTML filename (default: interactive_3d_frequency.html)',
    )
    parser.add_argument(
        '--true-hopf',
        action='store_true',
        help='Keep only points where the Hopf pair has |Re(lambda)| < tol',
    )
    parser.add_argument(
        '--hopf-tol',
        type=float,
        default=1e-3,
        help='Tolerance for true-Hopf filtering (default: 1e-3)',
    )
    parser.add_argument(
        '--marker-size',
        type=int,
        default=2,
        help='Marker size (default: 2)',
    )
    parser.add_argument(
        '--width',
        type=int,
        default=1400,
        help='Figure width in pixels (default: 1400)',
    )
    parser.add_argument(
        '--height',
        type=int,
        default=700,
        help='Figure height in pixels (default: 700)',
    )
    args = parser.parse_args()

    for p in args.inputs:
        if not Path(p).exists():
            raise FileNotFoundError(p)

    titles = [Path(p).parent.name for p in args.inputs]
    fig = build_figure(
        paths=args.inputs,
        titles=titles,
        true_hopf=args.true_hopf,
        hopf_tol=args.hopf_tol,
        marker_size=args.marker_size,
        width=args.width,
        height=args.height,
    )

    fig.write_html(args.output)
    print(f"Saved {args.output}")


if __name__ == '__main__':
    main()
