"""R1: Ramanujan's quarter-parameter coefficient lattice."""

import sympy as sp
from mpmath import mp

from thales_ramanujan import core
from thales_ramanujan.checks import check_ramanujan_coefficients
from thales_ramanujan.precision import DPS, use_dps


def test_group_passes():
    for record in check_ramanujan_coefficients():
        assert record.status == "PASS", (record.check_id, record.symbolic_residual)


def test_pochhammer_form_exact():
    """A_n = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3, in exact arithmetic."""
    for k in range(0, 26):
        assert sp.simplify(
            core.A_factorial_sym(k) - core.A_pochhammer_sym(k)
        ) == 0


def test_term_ratio_identity_symbolic():
    """The induction step, for symbolic n."""
    n = core.n_sym
    lhs = sp.simplify(core.A_factorial_sym(n + 1) / core.A_factorial_sym(n))
    rhs = (
        256
        * (sp.Rational(1, 4) + n)
        * (sp.Rational(1, 2) + n)
        * (sp.Rational(3, 4) + n)
        / (1 + n) ** 3
    )
    assert sp.simplify(sp.expand(lhs - rhs)) == 0


def test_generating_function_matches_3f2():
    with use_dps(DPS):
        z = mp.mpf(1) / 1000
        assert mp.fabs(
            core.F_series(z, terms=140) - core.F_hypergeometric(z)
        ) < mp.mpf("1e-50")


def test_first_coefficients():
    assert [core.A_exact(k) for k in range(4)] == [1, 24, 2520, 369600]
