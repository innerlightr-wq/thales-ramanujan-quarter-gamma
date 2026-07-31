"""Declared precision and tolerance policy.

Two regimes are used, and they are kept separate on purpose.

``DPS`` / ``TOL``
    Working precision for closed-form numerical evaluation (gamma, beta,
    digamma, elliptic integrals, hypergeometric sums).  These evaluations are
    accurate to essentially full working precision, so the tolerance is set
    close to it.

``QUAD_DPS`` / ``QUAD_TOL``
    Precision for *direct numerical quadrature* of the defining integrals.
    Several Thales integrands carry algebraic endpoint singularities (the
    arcsine density ``[a(1-a)]**(-1/2)`` being the extreme case).  Tanh-sinh
    quadrature handles them well but saturates well short of full working
    precision, so a separately declared, looser tolerance is used.  Reporting
    one global tolerance would either be dishonest about the quadrature or
    would needlessly weaken the closed-form audit.

No check anywhere in this package uses floating-point equality.
"""

from __future__ import annotations

from contextlib import contextmanager

from mpmath import mp, mpf

#: working decimal digits for closed-form numerical evaluation
DPS = 60

#: absolute-error tolerance for closed-form numerical checks
TOL = "1e-50"

#: working decimal digits for direct numerical quadrature
QUAD_DPS = 40

#: absolute-error tolerance for direct numerical quadrature
QUAD_TOL = "1e-20"

#: tanh-sinh maximum degree used for all quadrature in this package
QUAD_MAXDEGREE = 10


def set_precision(dps: int = DPS) -> None:
    """Set the global mpmath working precision."""
    mp.dps = dps


@contextmanager
def use_dps(dps: int):
    """Temporarily set the mpmath working precision."""
    old = mp.dps
    mp.dps = dps
    try:
        yield
    finally:
        mp.dps = old


def tol(value: str = TOL) -> mpf:
    """Return a tolerance as an exact-ish mpf at the current precision."""
    return mpf(value)


# Set a sane default at import time so that any consumer of this package gets
# high precision without having to remember to configure mpmath.
set_precision(DPS)
