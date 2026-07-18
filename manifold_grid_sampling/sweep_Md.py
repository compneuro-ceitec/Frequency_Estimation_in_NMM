#!/usr/bin/env python3
"""
1-D sweep over Md at constant Vt for the Gast-type neuromodulation model.

For each Md value the script:
  1. solves the Hopf manifold for Iext in (0, 300),
  2. solves the equilibrium quadratic for positive real r,
  3. reconstructs the full 4-D equilibrium (r, V, u, Sa),
  4. evaluates the Jacobian and its eigenvalues,
  5. extracts dominant / secondary frequencies.

Outputs:
  - CSV table: sweep_Md_Vt<vt>.csv
  - Multi-panel figure: sweep_Md_Vt<vt>.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from joblib import Parallel, delayed

import gast_coeffs as gc
import gast_hopf_frequency as ghf


# -----------------------------------------------------------------------------
# Fixed model parameters (same as gast_hopf_frequency.py)
# -----------------------------------------------------------------------------
_DEFAULT_PARAMS = {
    'B': sp.Rational(1, 5),
    'Ea': 0,
    'b': 1,
    'Tu': 50,
    'gaa': 1,
    'Eu': -90,
    'Sja': 15,
    'Ta': 6,
    'gammaa': sp.Rational(11, 10),
    'ujump': 45,
    'C': 50,
    'Vr': -60,
}
_PARAM_SUBS = {getattr(gc, k): v for k, v in _DEFAULT_PARAMS.items()}

# Short names for the fixed numeric values (used in analytic formulas)
_B = float(_DEFAULT_PARAMS['B'])
_Ea = float(_DEFAULT_PARAMS['Ea'])
_b = float(_DEFAULT_PARAMS['b'])
_Tu = float(_DEFAULT_PARAMS['Tu'])
_gaa = float(_DEFAULT_PARAMS['gaa'])
_Eu = float(_DEFAULT_PARAMS['Eu'])
_Sja = float(_DEFAULT_PARAMS['Sja'])
_Ta = float(_DEFAULT_PARAMS['Ta'])
_gammaa = float(_DEFAULT_PARAMS['gammaa'])
_ujump = float(_DEFAULT_PARAMS['ujump'])
_C = float(_DEFAULT_PARAMS['C'])
_Vr = float(_DEFAULT_PARAMS['Vr'])


# -----------------------------------------------------------------------------
# Build a lambdified Jacobian
# -----------------------------------------------------------------------------
def _build_jacobian():
    """Return a numpy callable J(r, V, u, Sa, Md, Vt, Iext)."""
    r_sym, V_sym, u_sym, Sa_sym = gc.r, gc.V, gc.u, gc.Sa
    f_r = (
        gc.gammaa * r_sym * (2 * V_sym - gc.Vr - gc.Vt)
        - r_sym * (gc.Md + gc.B) * gc.gaa * gc.Sa
    ) / gc.C
    f_V = (
        gc.gammaa * (V_sym - gc.Vr) * (V_sym - gc.Vt)
        - u_sym
        + gc.Iext
        + (gc.Md + gc.B) * gc.gaa * gc.Sa * (gc.Ea - V_sym)
        - gc.C * gc.pi * gc.r * (gc.C * gc.pi * gc.r / gc.gammaa)
    ) / gc.C
    f_u = (gc.b * (V_sym - gc.Vr) - (u_sym - gc.Eu)) / gc.Tu + gc.ujump * gc.r
    f_Sa = -Sa_sym / gc.Ta + gc.Sja * gc.r  # cexc = 0

    J_sym = sp.Matrix([f_r, f_V, f_u, f_Sa]).jacobian([r_sym, V_sym, u_sym, Sa_sym])
    J_param = J_sym.subs(_PARAM_SUBS)
    return sp.lambdify(
        (r_sym, V_sym, u_sym, Sa_sym, gc.Md, gc.Vt, gc.Iext),
        J_param,
        modules='numpy',
    )


_JACOBIAN = _build_jacobian()


def get_jacobian():
    """Return the lambdified Jacobian J(r, V, u, Sa, Md, Vt, Iext)."""
    return _JACOBIAN


def _build_r_solver():
    """Pre-compute and lambdify the symbolic roots of the equilibrium quadratic."""
    a0_expr = ghf.get_equilibrium_coeffs_fixed()['a0']
    a1_expr = ghf.get_equilibrium_coeffs_fixed()['a1']
    a2_expr = ghf.get_equilibrium_coeffs_fixed()['a2']
    poly = a2_expr * gc.r**2 + a1_expr * gc.r + a0_expr
    roots = sp.solve(poly, gc.r)
    return sp.lambdify((gc.Md, gc.Vt, gc.Iext), roots, modules='numpy')


_R_SOLVER = _build_r_solver()


# -----------------------------------------------------------------------------
# Core helpers
# -----------------------------------------------------------------------------
# The functions below are also used by gast_hopf_grid.py.
def equilibrium_state(Md, Vt, Iext, r):
    """Reconstruct (V, u, Sa) from r and the dr/dt=0, du/dt=0 equations."""
    Sa = _Ta * _Sja * r
    V = (_Vr + Vt) / 2.0 + (Md + _B) * _gaa * Sa / (2.0 * _gammaa)
    u = _b * (V - _Vr) + _Eu + _Tu * _ujump * r
    return V, u, Sa


def dVdt_residual(Md, Vt, Iext, r, V, u, Sa):
    """Check that the reconstructed equilibrium also satisfies dV/dt = 0."""
    dV = (
        _gammaa * (V - _Vr) * (V - Vt)
        - u
        + Iext
        + (Md + _B) * _gaa * Sa * (_Ea - V)
        - _C * np.pi * r * (_C * np.pi * r / _gammaa)
    ) / _C
    return dV


def solve_r_for_Md_Vt_Iext(Md, Vt, Iext, a0_expr=None, a1_expr=None, a2_expr=None):
    """Solve a2*r^2 + a1*r + a0 = 0 for positive real r.

    Uses a pre-lambdified quadratic formula. The a0/a1/a2 arguments are kept for
    backward compatibility but are no longer used.
    """
    roots = _R_SOLVER(Md, Vt, Iext)

    tol = 1e-10
    positive_real = []
    for root in roots:
        root_val = complex(root)
        if abs(root_val.imag) < tol and root_val.real > 0:
            positive_real.append(root_val.real)
    return positive_real


def eigenfrequencies(J):
    """Return (f_dom, f_sec, eigvals). Frequencies in Hz.

    Eigenvalues are in ms^-1 because all model time constants are in ms.
    Therefore |Im(lambda)|/(2*pi) is in kHz and is multiplied by 1000 to
    obtain Hz.

    The complex-conjugate pair whose real part is closest to zero is taken as
    the Hopf pair (dominant frequency); any additional complex pair gives the
    secondary frequency. This follows the convention in INSTRUCTIONS_next_steps.md.
    """
    eigvals = np.linalg.eigvals(J)
    tol = 1e-10
    # Extract one representative per complex-conjugate pair, keep |Re| and freq.
    pairs = []
    used = np.zeros(len(eigvals), dtype=bool)
    for i, ev in enumerate(eigvals):
        if used[i]:
            continue
        if abs(ev.imag) > tol:
            freq_hz = abs(ev.imag) / (2.0 * np.pi) * 1000.0
            pairs.append((abs(ev.real), freq_hz))
            # mark the conjugate
            for j in range(i + 1, len(eigvals)):
                if abs(eigvals[j] - ev.conjugate()) < tol:
                    used[j] = True
                    break
    # Sort by |Re| ascending: the Hopf pair has real part closest to zero.
    pairs.sort(key=lambda x: x[0])
    f_dom = pairs[0][1] if len(pairs) > 0 else np.nan
    f_sec = pairs[1][1] if len(pairs) > 1 else np.nan
    return f_dom, f_sec, eigvals


# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------
def _process_md_point(Md, vt, iext_min, iext_max):
    """Process one Md point in the 1-D sweep and return its records."""
    a0_expr = ghf.get_equilibrium_coeffs_fixed()['a0']
    a1_expr = ghf.get_equilibrium_coeffs_fixed()['a1']
    a2_expr = ghf.get_equilibrium_coeffs_fixed()['a2']
    hopf_func = ghf.solve_Iext_for_Md_Vt
    jac_func = get_jacobian()

    records = []
    iext_roots = hopf_func(Md, vt, Iext_min=iext_min, Iext_max=iext_max)
    for Iext in iext_roots:
        r_roots = solve_r_for_Md_Vt_Iext(Md, vt, Iext, a0_expr, a1_expr, a2_expr)
        for r in r_roots:
            V, u, Sa = equilibrium_state(Md, vt, Iext, r)
            # Verify dV/dt residual
            res_dV = dVdt_residual(Md, vt, Iext, r, V, u, Sa)
            J = jac_func(r, V, u, Sa, Md, vt, Iext)
            f_dom, f_sec, eigvals = eigenfrequencies(J)
            records.append({
                'Md': Md,
                'Vt': vt,
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


def run_sweep(vt, md_min=0.0, md_max=20.0, n_points=401,
              iext_min=0.0, iext_max=300.0, n_jobs=1):
    """Run the 1-D Md sweep and return a list of records.

    Parameters
    ----------
    n_jobs : int
        Number of parallel workers.  Use -1 to use all available cores.
    """
    md_grid = np.linspace(md_min, md_max, n_points)

    if n_jobs == 1:
        results = [
            _process_md_point(Md, vt, iext_min, iext_max)
            for Md in md_grid
        ]
    else:
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_process_md_point)(Md, vt, iext_min, iext_max)
            for Md in md_grid
        )

    return [rec for sublist in results for rec in sublist]


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
def save_results(records, vt, outdir='.'):
    """Save CSV (without eigenvalue arrays) and a numpy archive (with eigenvalues)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    suffix = f"_Vt{vt:.2f}"
    csv_path = outdir / f"sweep_Md{suffix}.csv"
    npz_path = outdir / f"sweep_Md{suffix}.npz"

    if records:
        # CSV without the eigenvalue arrays
        import csv
        keys = ['Md', 'Vt', 'Iext', 'r', 'V', 'u', 'Sa', 'dVdt_res', 'f_dom', 'f_sec']
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for rec in records:
                writer.writerow({k: rec[k] for k in keys})

        # NPZ with arrays
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


def plot_results(records, vt, outdir='outputs'):
    """Create a multi-panel figure of the sweep results."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig_path = outdir / f"sweep_Md_Vt{vt:.2f}.png"

    if not records:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No Hopf solutions found in this range",
                ha='center', va='center', transform=ax.transAxes)
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        return fig_path

    Md = np.array([r['Md'] for r in records])
    Iext = np.array([r['Iext'] for r in records])
    r = np.array([r['r'] for r in records])
    f_dom = np.array([r['f_dom'] for r in records])
    f_sec = np.array([r['f_sec'] for r in records])
    res = np.array([r['dVdt_res'] for r in records])

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    fig.suptitle(f"1-D Md sweep at Vt = {vt:.2f} mV")

    ax = axes[0, 0]
    ax.plot(Md, Iext, '.', ms=3)
    ax.set_ylabel("Iext")
    ax.set_title("Hopf locus Iext vs Md")
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[0, 1]
    ax.plot(Md, r, '.', ms=3)
    ax.set_ylabel("r")
    ax.set_title("Equilibrium firing rate vs Md")
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[1, 0]
    ax.plot(Md, f_dom, '.', ms=3, label='dominant')
    ax.plot(Md, f_sec, '.', ms=3, label='secondary')
    ax.set_xlabel("Md")
    ax.set_ylabel("frequency (Hz)")
    ax.set_title("Hopf frequency vs Md")
    ax.legend()
    ax.grid(True, ls=':', alpha=0.5)

    ax = axes[1, 1]
    ax.semilogy(Md, np.abs(res) + 1e-20, '.', ms=3)
    ax.set_xlabel("Md")
    ax.set_ylabel("|dV/dt residual|")
    ax.set_title("Equilibrium consistency check")
    ax.grid(True, ls=':', alpha=0.5)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    return fig_path


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vt', type=float, default=-43.0,
                        help='Constant Vt value (mV). Default: -43.0')
    parser.add_argument('--md-min', type=float, default=0.0,
                        help='Minimum Md. Default: 0.0')
    parser.add_argument('--md-max', type=float, default=20.0,
                        help='Maximum Md. Default: 20.0')
    parser.add_argument('--n-points', type=int, default=401,
                        help='Number of Md sample points. Default: 401')
    parser.add_argument('--iext-min', type=float, default=0.0,
                        help='Lower bound for Iext. Default: 0.0')
    parser.add_argument('--iext-max', type=float, default=300.0,
                        help='Upper bound for Iext. Default: 300.0')
    parser.add_argument('--outdir', type=str, default='outputs',
                        help='Output directory. Default: outputs')
    args = parser.parse_args()

    print(f"Running Md sweep: Md in [{args.md_min}, {args.md_max}], Vt = {args.vt}")
    records = run_sweep(
        vt=args.vt,
        md_min=args.md_min,
        md_max=args.md_max,
        n_points=args.n_points,
        iext_min=args.iext_min,
        iext_max=args.iext_max,
    )
    print(f"  found {len(records)} valid (Md, Vt, Iext, r) records")

    csv_path, npz_path = save_results(records, args.vt, args.outdir)
    fig_path = plot_results(records, args.vt, args.outdir)

    print(f"Saved CSV:  {csv_path}")
    print(f"Saved NPZ:  {npz_path}")
    print(f"Saved plot: {fig_path}")


if __name__ == '__main__':
    main()
