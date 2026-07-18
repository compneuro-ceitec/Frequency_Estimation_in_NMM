"""
Symbolic coefficients for the Gast-type Hopf manifold computation.

The module pre-computes (once, at import time) the equilibrium coefficients
a2, a1, a0 and the Hopf-condition coefficients b4..b0, and exposes them
through simple getter functions.

- The a-coefficients are entered directly as SymPy expressions.
- The b-coefficients are explicit SymPy expressions copied from the Maple
  ``*.math`` files (b0.math..b4.math).  Each expression already contains the
  multiplicative ``C`` factor, so the returned b's match the numerators used in
  Gast_Hopf_biproduct_Delta0.md.
"""

import sympy as sp
# -----------------------------------------------------------------------------
# Symbol declarations
# -----------------------------------------------------------------------------
r, V, u, Sa = sp.symbols('r V u Sa', real=True)

# Free grid parameters
Md, Vt, Iext = sp.symbols('Md Vt Iext', real=True)

# Model parameters (fixed during the grid search, but kept symbolic here)
# Assumptions must match those used in compute_hopf_r_bialternate.py.
C, gammaa, Vr, B, gaa, Ea, b, Eu, Tu, ujump, Sja, Ta = sp.symbols(
    'C gammaa Vr B gaa Ea b Eu Tu ujump Sja Ta', real=True
)
pi = sp.pi

# -----------------------------------------------------------------------------
# Equilibrium coefficients a2, a1, a0 (numerator of f2_EQ)
# These are taken directly from Gast_Hopf_biproduct_Delta0.md, eqs. (3.12)-(3.14).
# The external input cexc is set to zero (as in the PDF parameter setting).
# -----------------------------------------------------------------------------
_a2 = (
    -B**2 * Sja**2 * Ta**2 * gaa**2
    - 2 * B * Md * Sja**2 * Ta**2 * gaa**2
    - Md**2 * Sja**2 * Ta**2 * gaa**2
    - 4 * C**2 * pi**2
)

_a1 = (
    4 * B * Ea * Sja * Ta * gaa * gammaa
    - 2 * B * Sja * Ta * Vr * gaa * gammaa
    - 2 * B * Sja * Ta * Vt * gaa * gammaa
    + 4 * Ea * Md * Sja * Ta * gaa * gammaa
    - 2 * Md * Sja * Ta * Vr * gaa * gammaa
    - 2 * Md * Sja * Ta * Vt * gaa * gammaa
    - 2 * B * Sja * Ta * b * gaa
    - 2 * Md * Sja * Ta * b * gaa
    - 4 * Tu * gammaa * ujump
)

_a0 = (
    -Vr**2 * gammaa**2
    + 2 * Vr * Vt * gammaa**2
    - Vt**2 * gammaa**2
    + 2 * Vr * b * gammaa
    - 2 * Vt * b * gammaa
    - 4 * Eu * gammaa
    + 4 * Iext * gammaa
)

# -----------------------------------------------------------------------------
# Hopf-condition coefficients b4..b0
# Explicit SymPy expressions copied from the Maple .math files (b0.math..b4.math).
# They are written as plain sums (matching the style of the a-coefficients above).
# -----------------------------------------------------------------------------
_b0 = (
    C**4 * Ta * b
    + C**4 * Tu * b
    + C**3 * Ta**2 * b**2
)

_b1 = (
    B * C**4 * Sja * Ta**2 * gaa
    + B * C**4 * Sja * Ta * Tu * gaa
    + 2 * B * C**3 * Ea * Sja * Ta**3 * gaa * gammaa
    + 2 * B * C**3 * Ea * Sja * Ta**2 * Tu * gaa * gammaa
    - B * C**3 * Sja * Ta**3 * Vr * gaa * gammaa
    - B * C**3 * Sja * Ta**3 * Vt * gaa * gammaa
    - B * C**3 * Sja * Ta**2 * Tu * Vr * gaa * gammaa
    - B * C**3 * Sja * Ta**2 * Tu * Vt * gaa * gammaa
    - 2 * B * C**3 * Sja * Ta**2 * Tu * b * gaa
    - 2 * B * C**2 * Ea * Sja * Ta**3 * Tu * b * gaa * gammaa
    + 2 * B * C**2 * Ea * Sja * Ta**2 * Tu**2 * b * gaa * gammaa
    + B * C**2 * Sja * Ta**3 * Tu * Vr * b * gaa * gammaa
    + B * C**2 * Sja * Ta**3 * Tu * Vt * b * gaa * gammaa
    - B * C**2 * Sja * Ta**2 * Tu**2 * Vr * b * gaa * gammaa
    - B * C**2 * Sja * Ta**2 * Tu**2 * Vt * b * gaa * gammaa
    + C**4 * Md * Sja * Ta**2 * gaa
    + C**4 * Md * Sja * Ta * Tu * gaa
    + 2 * C**3 * Ea * Md * Sja * Ta**3 * gaa * gammaa
    + 2 * C**3 * Ea * Md * Sja * Ta**2 * Tu * gaa * gammaa
    - C**3 * Md * Sja * Ta**3 * Vr * gaa * gammaa
    - C**3 * Md * Sja * Ta**3 * Vt * gaa * gammaa
    - C**3 * Md * Sja * Ta**2 * Tu * Vr * gaa * gammaa
    - C**3 * Md * Sja * Ta**2 * Tu * Vt * gaa * gammaa
    - 2 * C**3 * Md * Sja * Ta**2 * Tu * b * gaa
    - 2 * C**3 * Ta * Tu**2 * gammaa * ujump
    - 2 * C**3 * Tu**3 * gammaa * ujump
    - 2 * C**2 * Ea * Md * Sja * Ta**3 * Tu * b * gaa * gammaa
    + 2 * C**2 * Ea * Md * Sja * Ta**2 * Tu**2 * b * gaa * gammaa
    + C**2 * Md * Sja * Ta**3 * Tu * Vr * b * gaa * gammaa
    + C**2 * Md * Sja * Ta**3 * Tu * Vt * b * gaa * gammaa
    - C**2 * Md * Sja * Ta**2 * Tu**2 * Vr * b * gaa * gammaa
    - C**2 * Md * Sja * Ta**2 * Tu**2 * Vt * b * gaa * gammaa
    + 2 * C**2 * Ta**3 * Tu * b * gammaa * ujump
    - 2 * C**2 * Ta**2 * Tu**2 * b * gammaa * ujump
)

_b2 = (
    - B**2 * C**3 * Sja**2 * Ta**4 * gaa**2
    - B**2 * C**3 * Sja**2 * Ta**3 * Tu * gaa**2
    + B**2 * C**3 * Sja**2 * Ta**2 * Tu**2 * gaa**2
    + 2 * B**2 * C**2 * Ea * Sja**2 * Ta**3 * Tu**2 * gaa**2 * gammaa
    - 2 * B**2 * C**2 * Ea * Sja**2 * Ta**2 * Tu**3 * gaa**2 * gammaa
    + B**2 * C**2 * Sja**2 * Ta**4 * Tu * b * gaa**2
    - B**2 * C**2 * Sja**2 * Ta**3 * Tu**2 * Vr * gaa**2 * gammaa
    - B**2 * C**2 * Sja**2 * Ta**3 * Tu**2 * Vt * gaa**2 * gammaa
    - B**2 * C**2 * Sja**2 * Ta**3 * Tu**2 * b * gaa**2
    + B**2 * C**2 * Sja**2 * Ta**2 * Tu**3 * Vr * gaa**2 * gammaa
    + B**2 * C**2 * Sja**2 * Ta**2 * Tu**3 * Vt * gaa**2 * gammaa
    - 4 * B**2 * C * Ea**2 * Sja**2 * Ta**3 * Tu**3 * gaa**2 * gammaa**2
    + 4 * B**2 * C * Ea * Sja**2 * Ta**3 * Tu**3 * Vr * gaa**2 * gammaa**2
    + 4 * B**2 * C * Ea * Sja**2 * Ta**3 * Tu**3 * Vt * gaa**2 * gammaa**2
    - B**2 * C * Sja**2 * Ta**3 * Tu**3 * Vr**2 * gaa**2 * gammaa**2
    - 2 * B**2 * C * Sja**2 * Ta**3 * Tu**3 * Vr * Vt * gaa**2 * gammaa**2
    - B**2 * C * Sja**2 * Ta**3 * Tu**3 * Vt**2 * gaa**2 * gammaa**2
    - 2 * B * C**3 * Md * Sja**2 * Ta**4 * gaa**2
    - 2 * B * C**3 * Md * Sja**2 * Ta**3 * Tu * gaa**2
    + 2 * B * C**3 * Md * Sja**2 * Ta**2 * Tu**2 * gaa**2
    + 4 * B * C**2 * Ea * Md * Sja**2 * Ta**3 * Tu**2 * gaa**2 * gammaa
    - 4 * B * C**2 * Ea * Md * Sja**2 * Ta**2 * Tu**3 * gaa**2 * gammaa
    + 2 * B * C**2 * Md * Sja**2 * Ta**4 * Tu * b * gaa**2
    - 2 * B * C**2 * Md * Sja**2 * Ta**3 * Tu**2 * Vr * gaa**2 * gammaa
    - 2 * B * C**2 * Md * Sja**2 * Ta**3 * Tu**2 * Vt * gaa**2 * gammaa
    - 2 * B * C**2 * Md * Sja**2 * Ta**3 * Tu**2 * b * gaa**2
    + 2 * B * C**2 * Md * Sja**2 * Ta**2 * Tu**3 * Vr * gaa**2 * gammaa
    + 2 * B * C**2 * Md * Sja**2 * Ta**2 * Tu**3 * Vt * gaa**2 * gammaa
    - 2 * B * C**2 * Sja * Ta**3 * Tu**2 * gaa * gammaa * ujump
    + 2 * B * C**2 * Sja * Ta**2 * Tu**3 * gaa * gammaa * ujump
    - 8 * B * C * Ea**2 * Md * Sja**2 * Ta**3 * Tu**3 * gaa**2 * gammaa**2
    + 8 * B * C * Ea * Md * Sja**2 * Ta**3 * Tu**3 * Vr * gaa**2 * gammaa**2
    + 8 * B * C * Ea * Md * Sja**2 * Ta**3 * Tu**3 * Vt * gaa**2 * gammaa**2
    + 8 * B * C * Ea * Sja * Ta**3 * Tu**3 * gaa * gammaa**2 * ujump
    - 2 * B * C * Md * Sja**2 * Ta**3 * Tu**3 * Vr**2 * gaa**2 * gammaa**2
    - 4 * B * C * Md * Sja**2 * Ta**3 * Tu**3 * Vr * Vt * gaa**2 * gammaa**2
    - 2 * B * C * Md * Sja**2 * Ta**3 * Tu**3 * Vt**2 * gaa**2 * gammaa**2
    - 4 * B * C * Sja * Ta**3 * Tu**3 * Vr * gaa * gammaa**2 * ujump
    - 4 * B * C * Sja * Ta**3 * Tu**3 * Vt * gaa * gammaa**2 * ujump
    + 4 * pi**2 * C**4 * Ta**3 * b
    + 4 * pi**2 * C**4 * Ta**2 * Tu * b
    - C**3 * Md**2 * Sja**2 * Ta**4 * gaa**2
    - C**3 * Md**2 * Sja**2 * Ta**3 * Tu * gaa**2
    + C**3 * Md**2 * Sja**2 * Ta**2 * Tu**2 * gaa**2
    + 2 * C**2 * Ea * Md**2 * Sja**2 * Ta**3 * Tu**2 * gaa**2 * gammaa
    - 2 * C**2 * Ea * Md**2 * Sja**2 * Ta**2 * Tu**3 * gaa**2 * gammaa
    + C**2 * Md**2 * Sja**2 * Ta**4 * Tu * b * gaa**2
    - C**2 * Md**2 * Sja**2 * Ta**3 * Tu**2 * Vr * gaa**2 * gammaa
    - C**2 * Md**2 * Sja**2 * Ta**3 * Tu**2 * Vt * gaa**2 * gammaa
    - C**2 * Md**2 * Sja**2 * Ta**3 * Tu**2 * b * gaa**2
    + C**2 * Md**2 * Sja**2 * Ta**2 * Tu**3 * Vr * gaa**2 * gammaa
    + C**2 * Md**2 * Sja**2 * Ta**2 * Tu**3 * Vt * gaa**2 * gammaa
    - 2 * C**2 * Md * Sja * Ta**3 * Tu**2 * gaa * gammaa * ujump
    + 2 * C**2 * Md * Sja * Ta**2 * Tu**3 * gaa * gammaa * ujump
    - 4 * C * Ea**2 * Md**2 * Sja**2 * Ta**3 * Tu**3 * gaa**2 * gammaa**2
    + 4 * C * Ea * Md**2 * Sja**2 * Ta**3 * Tu**3 * Vr * gaa**2 * gammaa**2
    + 4 * C * Ea * Md**2 * Sja**2 * Ta**3 * Tu**3 * Vt * gaa**2 * gammaa**2
    + 8 * C * Ea * Md * Sja * Ta**3 * Tu**3 * gaa * gammaa**2 * ujump
    - C * Md**2 * Sja**2 * Ta**3 * Tu**3 * Vr**2 * gaa**2 * gammaa**2
    - 2 * C * Md**2 * Sja**2 * Ta**3 * Tu**3 * Vr * Vt * gaa**2 * gammaa**2
    - C * Md**2 * Sja**2 * Ta**3 * Tu**3 * Vt**2 * gaa**2 * gammaa**2
    - 4 * C * Md * Sja * Ta**3 * Tu**3 * Vr * gaa * gammaa**2 * ujump
    - 4 * C * Md * Sja * Ta**3 * Tu**3 * Vt * gaa * gammaa**2 * ujump
    - 4 * C * Ta**3 * Tu**3 * gammaa**2 * ujump**2
)

_b3 = (
    - B**3 * C**2 * Sja**3 * Ta**4 * Tu**2 * gaa**3
    + B**3 * C**2 * Sja**3 * Ta**3 * Tu**3 * gaa**3
    + 4 * B**3 * C * Ea * Sja**3 * Ta**4 * Tu**3 * gaa**3 * gammaa
    - 2 * B**3 * C * Sja**3 * Ta**4 * Tu**3 * Vr * gaa**3 * gammaa
    - 2 * B**3 * C * Sja**3 * Ta**4 * Tu**3 * Vt * gaa**3 * gammaa
    - 3 * B**2 * C**2 * Md * Sja**3 * Ta**4 * Tu**2 * gaa**3
    + 3 * B**2 * C**2 * Md * Sja**3 * Ta**3 * Tu**3 * gaa**3
    + 12 * B**2 * C * Ea * Md * Sja**3 * Ta**4 * Tu**3 * gaa**3 * gammaa
    - 6 * B**2 * C * Md * Sja**3 * Ta**4 * Tu**3 * Vr * gaa**3 * gammaa
    - 6 * B**2 * C * Md * Sja**3 * Ta**4 * Tu**3 * Vt * gaa**3 * gammaa
    - 4 * B**2 * C * Sja**2 * Ta**4 * Tu**3 * gaa**2 * gammaa * ujump
    + 4 * pi**2 * B * C**4 * Sja * Ta**2 * Tu**2 * gaa
    + 4 * pi**2 * B * C**4 * Sja * Ta * Tu**3 * gaa
    + 8 * pi**2 * B * C**3 * Ea * Sja * Ta**3 * Tu**2 * gaa * gammaa
    + 8 * pi**2 * B * C**3 * Ea * Sja * Ta**2 * Tu**3 * gaa * gammaa
    - 4 * pi**2 * B * C**3 * Sja * Ta**3 * Tu**2 * Vr * gaa * gammaa
    - 4 * pi**2 * B * C**3 * Sja * Ta**3 * Tu**2 * Vt * gaa * gammaa
    - 4 * pi**2 * B * C**3 * Sja * Ta**2 * Tu**3 * Vr * gaa * gammaa
    - 4 * pi**2 * B * C**3 * Sja * Ta**2 * Tu**3 * Vt * gaa * gammaa
    - 3 * B * C**2 * Md**2 * Sja**3 * Ta**4 * Tu**2 * gaa**3
    + 3 * B * C**2 * Md**2 * Sja**3 * Ta**3 * Tu**3 * gaa**3
    + 12 * B * C * Ea * Md**2 * Sja**3 * Ta**4 * Tu**3 * gaa**3 * gammaa
    - 6 * B * C * Md**2 * Sja**3 * Ta**4 * Tu**3 * Vr * gaa**3 * gammaa
    - 6 * B * C * Md**2 * Sja**3 * Ta**4 * Tu**3 * Vt * gaa**3 * gammaa
    - 8 * B * C * Md * Sja**2 * Ta**4 * Tu**3 * gaa**2 * gammaa * ujump
    + 4 * pi**2 * C**4 * Md * Sja * Ta**2 * Tu**2 * gaa
    + 4 * pi**2 * C**4 * Md * Sja * Ta * Tu**3 * gaa
    + 8 * pi**2 * C**3 * Ea * Md * Sja * Ta**3 * Tu**2 * gaa * gammaa
    + 8 * pi**2 * C**3 * Ea * Md * Sja * Ta**2 * Tu**3 * gaa * gammaa
    - 4 * pi**2 * C**3 * Md * Sja * Ta**3 * Tu**2 * Vr * gaa * gammaa
    - 4 * pi**2 * C**3 * Md * Sja * Ta**3 * Tu**2 * Vt * gaa * gammaa
    - 4 * pi**2 * C**3 * Md * Sja * Ta**2 * Tu**3 * Vr * gaa * gammaa
    - 4 * pi**2 * C**3 * Md * Sja * Ta**2 * Tu**3 * Vt * gaa * gammaa
    - 8 * pi**2 * C**3 * Ta**3 * Tu**2 * gammaa * ujump
    - 8 * pi**2 * C**3 * Ta**2 * Tu**3 * gammaa * ujump
    - C**2 * Md**3 * Sja**3 * Ta**4 * Tu**2 * gaa**3
    + C**2 * Md**3 * Sja**3 * Ta**3 * Tu**3 * gaa**3
    + 4 * C * Ea * Md**3 * Sja**3 * Ta**4 * Tu**3 * gaa**3 * gammaa
    - 2 * C * Md**3 * Sja**3 * Ta**4 * Tu**3 * Vr * gaa**3 * gammaa
    - 2 * C * Md**3 * Sja**3 * Ta**4 * Tu**3 * Vt * gaa**3 * gammaa
    - 4 * C * Md**2 * Sja**2 * Ta**4 * Tu**3 * gaa**2 * gammaa * ujump
)

_b4 = (
    - B**4 * C * Sja**4 * Ta**5 * Tu**3 * gaa**4
    - 4 * B**3 * C * Md * Sja**4 * Ta**5 * Tu**3 * gaa**4
    - 4 * pi**2 * B**2 * C**3 * Sja**2 * Ta**4 * Tu**2 * gaa**2
    - 4 * pi**2 * B**2 * C**3 * Sja**2 * Ta**3 * Tu**3 * gaa**2
    - 6 * B**2 * C * Md**2 * Sja**4 * Ta**5 * Tu**3 * gaa**4
    - 8 * pi**2 * B * C**3 * Md * Sja**2 * Ta**4 * Tu**2 * gaa**2
    - 8 * pi**2 * B * C**3 * Md * Sja**2 * Ta**3 * Tu**3 * gaa**2
    - 4 * B * C * Md**3 * Sja**4 * Ta**5 * Tu**3 * gaa**4
    - 4 * pi**2 * C**3 * Md**2 * Sja**2 * Ta**4 * Tu**2 * gaa**2
    - 4 * pi**2 * C**3 * Md**2 * Sja**2 * Ta**3 * Tu**3 * gaa**2
    - C * Md**4 * Sja**4 * Ta**5 * Tu**3 * gaa**4
)

# -----------------------------------------------------------------------------
# Hopf manifold resultant in terms of the a- and b-coefficients
# This is the discriminant-like expression obtained by eliminating the
# equilibrium variable r.  See hopf_manifold.math for the maple form.
# -----------------------------------------------------------------------------
_hopf_manifold = (
    _a0**4 * _b4**2
    - _a0**3 * _a1 * _b3 * _b4
    + _a0**3 * _a2 * (-2 * _b2 * _b4 + _b3**2)
    + _a0**2 * _a1**2 * _b2 * _b4
    + _a0**2 * _a1 * _a2 * (3 * _b1 * _b4 - _b3 * _b2)
    + _a0**2 * _a2**2 * (2 * _b0 * _b4 - 2 * _b3 * _b1 + _b2**2)
    - _a0 * _a1**3 * _b1 * _b4
    + _a0 * _a1**2 * _a2 * (-4 * _b0 * _b4 + _b3 * _b1)
    + _a0 * _a1 * _a2**2 * (3 * _b0 * _b3 - _b1 * _b2)
    + _a0 * _a2**3 * (-2 * _b0 * _b2 + _b1**2)
    + _a1**4 * _b0 * _b4
    - _a1**3 * _a2 * _b0 * _b3
    + _a1**2 * _a2**2 * _b0 * _b2
    - _a1 * _a2**3 * _b0 * _b1
    + _a2**4 * _b0**2
)

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def get_equilibrium_coeffs():
    """Return dict with symbolic a2, a1, a0."""
    return {'a2': _a2, 'a1': _a1, 'a0': _a0}


def get_hopf_coeffs():
    """Return dict with symbolic b4, b3, b2, b1, b0 (1/C factor removed)."""
    return {'b4': _b4, 'b3': _b3, 'b2': _b2, 'b1': _b1, 'b0': _b0}


def get_all_coeffs():
    """Return all eight coefficients in one dict."""
    return {**get_equilibrium_coeffs(), **get_hopf_coeffs()}


def get_hopf_manifold():
    """Return the symbolic Hopf manifold expression (resultant in a/b coefficients)."""
    return _hopf_manifold


if __name__ == '__main__':
    print("Equilibrium coefficients:")
    for name, expr in get_equilibrium_coeffs().items():
        print(f"  {name} = {expr}")
    print("\nHopf coefficients (1/C removed):")
    for name, expr in get_hopf_coeffs().items():
        print(f"  {name} = {expr}")

    print(f"\nHopf manifold expression: {len(_hopf_manifold.args)}-term resultant in a0, a1, a2, b0..b4")
    print(f"  Free symbols in full parameter form: {sorted(str(s) for s in _hopf_manifold.free_symbols)}")
