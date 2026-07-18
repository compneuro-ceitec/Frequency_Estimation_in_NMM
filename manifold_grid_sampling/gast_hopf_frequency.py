"""
Frequency computation for the Gast-type neuromodulation model.

This module takes the symbolic Hopf manifold from ``gast_coeffs``, substitutes
the fixed model parameters (copied from ``params.md`` plus the defaults from
``alg_notes.md``), and prepares the expressions needed to solve for the Hopf
locus and the corresponding Jacobian eigenfrequencies.
"""

import sympy as sp
import gast_coeffs as gc

# -----------------------------------------------------------------------------
# Fixed model parameters
# -----------------------------------------------------------------------------
# Values copied from params.md (Maple parameter setting).
# C and Vr are not listed in params.md, so their defaults from alg_notes.md
# are used instead.
_DEFAULT_PARAMS = {
    'B': sp.Rational(1, 5),          # 1/5
    'Ea': 0,                          # mV
    'b': 1,                           # ms^-1
    'Tu': 50,                         # ms
    'gaa': 1,                         # nS? (dimensionless in this setting)
    'Eu': -90,                        # mV
    'Sja': 15,                        # dimensionless gain
    'Ta': 6,                          # ms
    'gammaa': sp.Rational(11, 10),   # 11/10
    'ujump': 45,                      # mV ms^-1
    'C': 50,                          # from alg_notes.md default
    'Vr': -60,                        # from alg_notes.md default
}

# Verify that the parameter names match declared symbols.
_unknown = [k for k in _DEFAULT_PARAMS if not hasattr(gc, k)]
if _unknown:
    raise NameError(f"Unknown parameter keys: {_unknown}")

_PARAM_SUBS = {getattr(gc, k): v for k, v in _DEFAULT_PARAMS.items()}

# -----------------------------------------------------------------------------
# Substitute fixed parameters into the Hopf manifold and the a/b coefficients
# -----------------------------------------------------------------------------
_hopf_manifold_full = gc.get_hopf_manifold()
_hopf_manifold = _hopf_manifold_full.subs(_PARAM_SUBS)

_a0 = gc._a0.subs(_PARAM_SUBS)
_a1 = gc._a1.subs(_PARAM_SUBS)
_a2 = gc._a2.subs(_PARAM_SUBS)

_b = {f'b{i}': getattr(gc, f'_b{i}').subs(_PARAM_SUBS) for i in range(5)}

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------
def get_hopf_manifold_fixed():
    """Return Hopf manifold expression with fixed parameters substituted."""
    return _hopf_manifold


def get_equilibrium_coeffs_fixed():
    """Return a2, a1, a0 with fixed parameters substituted."""
    return {'a2': _a2, 'a1': _a1, 'a0': _a0}


def get_hopf_coeffs_fixed():
    """Return b4..b0 with fixed parameters substituted."""
    return _b


def solve_Iext_for_Md_Vt(Md_val, Vt_val, Iext_min=0.0, Iext_max=300.0):
    """
    For given numerical Md and Vt, solve Hopf_manifold = 0 for Iext.

    Returns all real roots that lie inside (Iext_min, Iext_max).  The Hopf
    manifold is a quartic in Iext, so up to four roots may be returned.
    """
    expr = _hopf_manifold.subs({gc.Md: Md_val, gc.Vt: Vt_val})
    # Numerical root finding for the univariate polynomial in Iext.
    roots = sp.nroots(sp.Poly(expr, gc.Iext))
    real_roots = []
    for r in roots:
        if abs(sp.im(r)) < 1e-8:
            x = float(sp.re(r))
            if Iext_min < x < Iext_max:
                real_roots.append(x)
    # Return sorted unique roots.
    seen = set()
    unique = []
    for x in sorted(real_roots):
        key = round(x, 12)
        if key not in seen:
            seen.add(key)
            unique.append(x)
    return unique


if __name__ == '__main__':
    print("Free symbols after substituting fixed parameters:")
    print(f"  Hopf manifold: {sorted(str(s) for s in _hopf_manifold.free_symbols)}")
    print(f"  a2: {sorted(str(s) for s in _a2.free_symbols)}")
    print(f"  a1: {sorted(str(s) for s in _a1.free_symbols)}")
    print(f"  a0: {sorted(str(s) for s in _a0.free_symbols)}")
    print(f"\nHopf manifold degree in each remaining symbol:")
    for sym in sorted(_hopf_manifold.free_symbols, key=str):
        print(f"  degree in {sym}: {sp.degree(_hopf_manifold, gen=sym)}")

    # Example: Md = 5.0, Vt = -43.0
    Md_test = 5.0
    Vt_test = -43.0
    Iext_solutions = solve_Iext_for_Md_Vt(Md_test, Vt_test)
    print(f"\nReal Iext solutions for Md={Md_test}, Vt={Vt_test} in (0, 300):")
    for sol in Iext_solutions:
        print(f"  Iext = {sol:.6f}")
