"""Core objects: Ramanujan's coefficient lattice and the Thales beta family.

Symbolic objects use SymPy with exact rational parameters.  Numerical objects
use mpmath at the precision declared in :mod:`thales_ramanujan.precision`.

Notation follows the manuscript:

* ``A_n = (4n)! / (n!)^4`` is the coefficient of Ramanujan's 1/pi series.
* ``J(s) = int_0^1 [a(1-a)]^s da = B(s+1, s+1) = Gamma(s+1)^2 / Gamma(2s+2)``
  is the Thales partition beta family attached to ``a + b = 1``.
* ``I_s(alpha) = B(s+alpha+1, s-alpha+1)`` is the mirror-tilted family.
* ``A_s(alpha) = int_0^1 (2a-1) a^(s+alpha) (1-a)^(s-alpha) da``.
"""

from __future__ import annotations

import sympy as sp
from mpmath import mp, mpf

from .precision import QUAD_MAXDEGREE

# ---------------------------------------------------------------------------
# symbols
# ---------------------------------------------------------------------------

a_sym, s_sym, alpha_sym, z_sym = sp.symbols("a s alpha z")
n_sym = sp.Symbol("n", integer=True, nonnegative=True)

#: the quarter-integer residue lattice shared by both sides of the comparison
QUARTER_LATTICE = (sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4))

#: upper/lower parameters of the Ramanujan 3F2 of manuscript eq. (3)
RAMANUJAN_3F2_UPPER = QUARTER_LATTICE
RAMANUJAN_3F2_LOWER = (sp.Integer(1), sp.Integer(1))

#: upper parameters of the Clausen square root of manuscript eq. (4)
CLAUSEN_2F1_PARAMS = (sp.Rational(1, 8), sp.Rational(3, 8), sp.Integer(1))


# ---------------------------------------------------------------------------
# Ramanujan side
# ---------------------------------------------------------------------------

def A_factorial_sym(n):
    """A_n = (4n)! / (n!)^4 (symbolic)."""
    return sp.factorial(4 * n) / sp.factorial(n) ** 4


def A_pochhammer_sym(n):
    """256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3 (symbolic)."""
    return (
        256 ** n
        * sp.rf(sp.Rational(1, 4), n)
        * sp.rf(sp.Rational(1, 2), n)
        * sp.rf(sp.Rational(3, 4), n)
        / sp.rf(sp.Integer(1), n) ** 3
    )


def A_binomial_sym(n):
    """C(4n,2n) C(2n,n)^2 (symbolic) -- the central-binomial decomposition."""
    return sp.binomial(4 * n, 2 * n) * sp.binomial(2 * n, n) ** 2


def A_exact(n: int) -> int:
    """A_n as an exact integer."""
    return int(sp.factorial(4 * n) // sp.factorial(n) ** 4)


def F_series(z, terms: int = 60):
    """Partial sum of F(z) = sum A_n z^n, evaluated numerically."""
    z = mpf(z)
    total = mp.mpf(0)
    for k in range(terms):
        total += mpf(A_exact(k)) * z ** k
    return total


def F_hypergeometric(z):
    """F(z) = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z), evaluated numerically."""
    z = mpf(z)
    return mp.hyper(
        [mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4], [1, 1], 256 * z
    )


def clausen_square_root(x):
    """2F1(1/8, 3/8; 1; 256 z) with argument x = 256 z, evaluated numerically."""
    return mp.hyp2f1(mp.mpf(1) / 8, mp.mpf(3) / 8, 1, mpf(x))


def ramanujan_partial_sum(terms: int = 40):
    """Partial sum of Ramanujan's series for 1/pi (manuscript eq. (1))."""
    total = mp.mpf(0)
    for k in range(terms):
        total += mpf(A_exact(k)) * (1103 + 26390 * k) / mpf(396) ** (4 * k)
    return 2 * mp.sqrt(2) / 9801 * total


# ---------------------------------------------------------------------------
# Thales side
# ---------------------------------------------------------------------------

def J_sym(s):
    """J(s) = Gamma(s+1)^2 / Gamma(2s+2) (symbolic)."""
    return sp.gamma(s + 1) ** 2 / sp.gamma(2 * s + 2)


def J_beta_sym(s):
    """J(s) = B(s+1, s+1) (symbolic)."""
    return sp.beta(s + 1, s + 1)


def J_integrand_sym(a, s):
    """[a(1-a)]^s (symbolic)."""
    return (a * (1 - a)) ** s


def J_num(s):
    """J(s) via the beta function (closed form, full working precision)."""
    s = mpf(s)
    return mp.beta(s + 1, s + 1)


def J_quad(s):
    """J(s) by direct tanh-sinh quadrature of the defining integral."""
    s = mpf(s)
    return mp.quad(
        lambda a: (a * (1 - a)) ** s,
        [0, 1],
        method="tanh-sinh",
        maxdegree=QUAD_MAXDEGREE,
    )


def I_sym(s, alpha):
    """I_s(alpha) = B(s+alpha+1, s-alpha+1) (symbolic)."""
    return sp.beta(s + alpha + 1, s - alpha + 1)


def A_tilt_sym(s, alpha):
    """A_s(alpha) as a difference of beta functions (symbolic).

    ``2a - 1 = a - (1-a)`` gives
    ``A_s(alpha) = B(s+alpha+2, s-alpha+1) - B(s+alpha+1, s-alpha+2)``.
    """
    return sp.beta(s + alpha + 2, s - alpha + 1) - sp.beta(s + alpha + 1, s - alpha + 2)


def I_num(s, alpha):
    """I_s(alpha) via the beta function."""
    s, alpha = mpf(s), mpf(alpha)
    return mp.beta(s + alpha + 1, s - alpha + 1)


def I_quad(s, alpha):
    """I_s(alpha) by direct quadrature."""
    s, alpha = mpf(s), mpf(alpha)
    return mp.quad(
        lambda a: a ** (s + alpha) * (1 - a) ** (s - alpha),
        [0, 1],
        method="tanh-sinh",
        maxdegree=QUAD_MAXDEGREE,
    )


def A_tilt_quad(s, alpha):
    """A_s(alpha) by direct quadrature."""
    s, alpha = mpf(s), mpf(alpha)
    return mp.quad(
        lambda a: (2 * a - 1) * a ** (s + alpha) * (1 - a) ** (s - alpha),
        [0, 1],
        method="tanh-sinh",
        maxdegree=QUAD_MAXDEGREE,
    )


def log_I_alpha_derivative(s, alpha):
    """d/d(alpha) log I_s(alpha), by numerical differentiation."""
    s = mpf(s)
    return mp.diff(lambda t: mp.log(I_num(s, t)), mpf(alpha))


def digamma_response(s, alpha):
    """psi(s+alpha+1) - psi(s-alpha+1)."""
    s, alpha = mpf(s), mpf(alpha)
    return mp.digamma(s + alpha + 1) - mp.digamma(s - alpha + 1)


# ---------------------------------------------------------------------------
# arcsine / lemniscatic
# ---------------------------------------------------------------------------

def arcsine_moment_quad(n: int):
    """(1/pi) int_0^1 x^n / sqrt(x(1-x)) dx by direct quadrature."""
    return (
        mp.quad(
            lambda x: x ** n / mp.sqrt(x * (1 - x)),
            [0, 1],
            method="tanh-sinh",
            maxdegree=QUAD_MAXDEGREE,
        )
        / mp.pi
    )


def arcsine_moment_closed(n: int):
    """4^{-n} C(2n, n)."""
    return mp.binomial(2 * n, n) / mpf(4) ** n


def K_lemniscatic():
    """K(1/sqrt(2)).

    mpmath's ``ellipk`` takes the *parameter* m = k^2, so the modulus
    k = 1/sqrt(2) of the manuscript corresponds to m = 1/2.
    """
    return mp.ellipk(mp.mpf(1) / 2)


def K_lemniscatic_gamma():
    """Gamma(1/4)^2 / (4 sqrt(pi)) -- the classical lemniscatic evaluation."""
    return mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))


# ---------------------------------------------------------------------------
# generating functions of the Clausen section
# ---------------------------------------------------------------------------

def G_beta_series(z, terms: int = 40):
    """sum_{n>=0} J(n) z^n, evaluated as a partial sum with exact rational J(n)."""
    z = mpf(z)
    total = mp.mpf(0)
    for k in range(terms):
        coeff = sp.Rational(sp.factorial(k) ** 2, sp.factorial(2 * k + 1))
        total += mpf(sp.Float(coeff, mp.dps + 10)) * z ** k
    return total


def G_beta_closed(z):
    """2F1(1, 1; 3/2; z/4)."""
    z = mpf(z)
    return mp.hyp2f1(1, 1, mp.mpf(3) / 2, z / 4)


def G_arcsine_series(z, terms: int = 60):
    """sum_{n>=0} 4^{-n} C(2n,n) z^n as a partial sum."""
    z = mpf(z)
    total = mp.mpf(0)
    for k in range(terms):
        total += arcsine_moment_closed(k) * z ** k
    return total


def G_arcsine_closed(z):
    """(1 - z)^{-1/2}."""
    z = mpf(z)
    return (1 - z) ** mpf("-0.5")


def clausen_condition_holds(a, b, c) -> bool:
    """True iff c = a + b + 1/2 exactly (rational arithmetic)."""
    return sp.simplify(sp.nsimplify(c) - (sp.nsimplify(a) + sp.nsimplify(b) + sp.Rational(1, 2))) == 0
