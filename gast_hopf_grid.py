#!/usr/bin/env python3
"""
2-D grid sweep over Md and Vt for the Gast-type neuromodulation model.

For each (Md, Vt) pair the script:
  1. extracts the quartic Hopf-manifold coefficients in Iext,
  2. finds real Iext roots in (0, 300) using NumPy (fast),
  3. solves the equilibrium quadratic for positive real r via SymPy,
  4. reconstructs the 4-D equilibrium and evaluates Jacobian eigenvalues,
  5. extracts dominant / secondary frequencies.

Outputs:
  - CSV table: gast_hopf_grid_results.csv
  - NPZ archive: gast_hopf_grid_results.npz
  - Figures: gast_hopf_grid.png (2-D overview) and
             gast_hopf_Iext_vs_Md.png (1-D cuts)
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from joblib import Parallel, delayed

import gast_coeffs as gc
import gast_hopf_frequency as ghf
from sweep_Md import (
    equilibrium_state,
    dVdt_residual,
    solve_r_for_Md_Vt_Iext,
    eigenfrequencies,
    get_jacobian,
)


# -----------------------------------------------------------------------------
# Fast extraction of Hopf-manifold quartic coefficients in Iext
# -----------------------------------------------------------------------------
def _build_hopf_coeffs_func():
    """Return a numpy callable that returns the quartic coefficients in Iext."""
    hopf_expr = ghf.get_hopf_manifold_fixed()
    poly = sp.Poly(hopf_expr, gc.Iext)
    coeffs = poly.all_coeffs()  # highest to lowest power
    return sp.lambdify((gc.Md, gc.Vt), coeffs, modules='numpy')


_HOPF_COEFFS_FUNC = _build_hopf_coeffs_func()


def iext_roots_for_Md_Vt(Md, Vt, iext_min=0.0, iext_max=300.0):
    """Return real Iext roots of Hopf_manifold(Md, Vt, Iext)=0 inside the bounds."""
    coeffs = _HOPF_COEFFS_FUNC(Md, Vt)
    roots = np.roots(coeffs)
    real_roots = []
    for root in roots:
        if abs(root.imag) < 1e-8:
            x = root.real
            if iext_min < x < iext_max:
                real_roots.append(x)
    # Unique up to tolerance
    seen = set()
    unique = []
    for x in sorted(real_roots):
        key = round(x, 12)
        if key not in seen:
            seen.add(key)
            unique.append(x)
    return unique


# -----------------------------------------------------------------------------
# Grid computation
# -----------------------------------------------------------------------------
def _process_vt_slice(vt_info, md_grid, iext_min, iext_max):
    """Process one Vt slice: return all records for (Md in md_grid, Vt fixed)."""
    vi, Vt = vt_info
    a0_expr = ghf.get_equilibrium_coeffs_fixed()['a0']
    a1_expr = ghf.get_equilibrium_coeffs_fixed()['a1']
    a2_expr = ghf.get_equilibrium_coeffs_fixed()['a2']
    jac_func = get_jacobian()

    coeffs = _HOPF_COEFFS_FUNC(md_grid, Vt)  # shape (5, n_md)
    records = []
    for mi, Md in enumerate(md_grid):
        iext_roots = iext_roots_for_Md_Vt(Md, Vt, iext_min, iext_max)
        for Iext in iext_roots:
            r_roots = solve_r_for_Md_Vt_Iext(
                Md, Vt, Iext, a0_expr, a1_expr, a2_expr
            )
            for r in r_roots:
                V, u, Sa = equilibrium_state(Md, Vt, Iext, r)
                res_dV = dVdt_residual(Md, Vt, Iext, r, V, u, Sa)
                J = jac_func(r, V, u, Sa, Md, Vt, Iext)
                f_dom, f_sec, eigvals = eigenfrequencies(J)
                records.append({
                    'Md': Md,
                    'Vt': Vt,
                    'Iext': Iext,
                    'r': r,
                    'V': V,
                    'u': u,
                    'Sa': Sa,
                    'dVdt_res': res_dV,
                    'f_dom': f_dom,
                    'f_sec': f_sec,
                    'eigvals': eigvals,
                })
    return records


def run_grid(md_min=0.0, md_max=20.0, n_md=401,
             vt_min=-55.0, vt_max=-35.0, n_vt=101,
             iext_min=0.0, iext_max=300.0,
             n_jobs=1):
    """Run the 2-D (Md, Vt) grid and return a list of records.

    Parameters
    ----------
    n_jobs : int
        Number of parallel workers.  Use -1 to use all available cores.
        Parallelism is over Vt slices.
    """
    md_grid = np.linspace(md_min, md_max, n_md)
    vt_grid = np.linspace(vt_min, vt_max, n_vt)
    vt_infos = list(enumerate(vt_grid))

    t_start = time.time()

    if n_jobs == 1:
        results = [
            _process_vt_slice(vt_info, md_grid, iext_min, iext_max)
            for vt_info in vt_infos
        ]
    else:
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_process_vt_slice)(vt_info, md_grid, iext_min, iext_max)
            for vt_info in vt_infos
        )

    records = [rec for sublist in results for rec in sublist]

    elapsed = time.time() - t_start
    n_total = n_md * n_vt
    print(f"  {n_total}/{n_total} points processed, {len(records)} records, "
          f"elapsed {elapsed:.1f}s", file=sys.stderr)

    return records


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
def save_results(records, outdir='.'):
    """Save CSV (without eigenvalue arrays) and a numpy archive (with eigenvalues)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "gast_hopf_grid_results.csv"
    npz_path = outdir / "gast_hopf_grid_results.npz"

    if records:
        keys = ['Md', 'Vt', 'Iext', 'r', 'V', 'u', 'Sa', 'dVdt_res', 'f_dom', 'f_sec']
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for rec in records:
                writer.writerow({k: rec[k] for k in keys})

        np.savez(
            npz_path,
            Md=np.array([r['Md'] for r in records]),
            Vt=np.array([r['Vt'] for r in records]),
            Iext=np.array([r['Iext'] for r in records]),
            r=np.array([r['r'] for r in records]),
            V=np.array([r['V'] for r in records]),
            u=np.array([r['u'] for r in records]),
            Sa=np.array([r['Sa'] for r in records]),
            dVdt_res=np.array([r['dVdt_res'] for r in records]),
            f_dom=np.array([r['f_dom'] for r in records]),
            f_sec=np.array([r['f_sec'] for r in records]),
            eigvals=np.array([r['eigvals'] for r in records]),
        )
    else:
        csv_path.write_text("Md,Vt,Iext,r,V,u,Sa,dVdt_res,f_dom,f_sec\n")
        np.savez(npz_path)

    return csv_path, npz_path


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_grid(records, outdir='.', vt_cuts=None):
    """Create 2-D overview plots and 1-D cuts."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    outdir = Path(outdir)

    if not records:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No Hopf solutions found in this grid",
                ha='center', va='center', transform=ax.transAxes)
        fig.savefig(outdir / "gast_hopf_grid.png", dpi=150)
        plt.close(fig)
        return outdir / "gast_hopf_grid.png", None

    Md = np.array([r['Md'] for r in records])
    Vt = np.array([r['Vt'] for r in records])
    Iext = np.array([r['Iext'] for r in records])
    f_dom = np.array([r['f_dom'] for r in records])
    f_sec = np.array([r['f_sec'] for r in records])
    r = np.array([r['r'] for r in records])

    # ---- 2-D overview figure ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"2-D Hopf grid: Md vs Vt ({len(records)} records)")

    ax = axes[0, 0]
    sc = ax.scatter(Md, Vt, c=Iext, s=8, cmap='viridis', vmin=0, vmax=300)
    ax.set_xlabel("Md")
    ax.set_ylabel("Vt")
    ax.set_title("Iext roots")
    plt.colorbar(sc, ax=ax, label="Iext")
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[0, 1]
    sc = ax.scatter(Md, Vt, c=f_dom, s=8, cmap='plasma')
    ax.set_xlabel("Md")
    ax.set_ylabel("Vt")
    ax.set_title("Dominant frequency f_dom (Hz)")
    plt.colorbar(sc, ax=ax, label="Hz")
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[1, 0]
    sc = ax.scatter(Md, Vt, c=f_sec, s=8, cmap='plasma')
    ax.set_xlabel("Md")
    ax.set_ylabel("Vt")
    ax.set_title("Secondary frequency f_sec (Hz)")
    plt.colorbar(sc, ax=ax, label="Hz")
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[1, 1]
    sc = ax.scatter(Md, Vt, c=r, s=8, cmap='cividis')
    ax.set_xlabel("Md")
    ax.set_ylabel("Vt")
    ax.set_title("Equilibrium firing rate r")
    plt.colorbar(sc, ax=ax, label="r")
    ax.grid(True, ls=':', alpha=0.5)

    fig.tight_layout()
    grid_png = outdir / "gast_hopf_grid.png"
    fig.savefig(grid_png, dpi=150)
    plt.close(fig)

    # ---- 1-D cuts: Iext vs Md for fixed Vt ----
    if vt_cuts is None:
        # Default cuts across the Vt range.
        vt_cuts = np.linspace(Vt.min(), Vt.max(), 5)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    axes = axes.ravel()
    fig.suptitle("1-D cuts: Iext vs Md for fixed Vt")

    for ax, vt0 in zip(axes, vt_cuts):
        mask = np.abs(Vt - vt0) < 1e-8
        if mask.any():
            ax.plot(Md[mask], Iext[mask], '.', ms=4)
            ax.set_title(f"Vt = {vt0:.1f}")
        else:
            ax.set_title(f"Vt = {vt0:.1f} (no data)")
        ax.set_ylabel("Iext")
        ax.set_ylim(0, 300)
        ax.grid(True, ls=':', alpha=0.5)

    for ax in axes[-2:]:
        ax.set_xlabel("Md")
    fig.tight_layout()
    cuts_png = outdir / "gast_hopf_Iext_vs_Md.png"
    fig.savefig(cuts_png, dpi=150)
    plt.close(fig)

    return grid_png, cuts_png


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--md-min', type=float, default=0.0)
    parser.add_argument('--md-max', type=float, default=20.0)
    parser.add_argument('--n-md', type=int, default=401)
    parser.add_argument('--vt-min', type=float, default=-55.0)
    parser.add_argument('--vt-max', type=float, default=-35.0)
    parser.add_argument('--n-vt', type=int, default=101)
    parser.add_argument('--iext-min', type=float, default=0.0)
    parser.add_argument('--iext-max', type=float, default=300.0)
    parser.add_argument('--outdir', type=str, default='.')
    args = parser.parse_args()

    print(f"Running 2-D grid: Md in [{args.md_min}, {args.md_max}] x {args.n_md}, "
          f"Vt in [{args.vt_min}, {args.vt_max}] x {args.n_vt}")
    records = run_grid(
        md_min=args.md_min, md_max=args.md_max, n_md=args.n_md,
        vt_min=args.vt_min, vt_max=args.vt_max, n_vt=args.n_vt,
        iext_min=args.iext_min, iext_max=args.iext_max,
    )
    print(f"Total records: {len(records)}")

    csv_path, npz_path = save_results(records, args.outdir)
    grid_png, cuts_png = plot_grid(records, args.outdir)

    print(f"Saved CSV:       {csv_path}")
    print(f"Saved NPZ:       {npz_path}")
    print(f"Saved grid plot: {grid_png}")
    if cuts_png:
        print(f"Saved cuts plot: {cuts_png}")


if __name__ == '__main__':
    main()
