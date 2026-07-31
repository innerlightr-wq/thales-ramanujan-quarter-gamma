"""R7: the Clausen obstruction for the tested generating functions.

The scope of these tests is exactly the scope of the manuscript: the natural
elementary Thales generating functions examined here fail the Clausen
parameter condition.  Nothing here rules out every conceivable construction.
"""

import sympy as sp
from mpmath import mp

from thales_ramanujan import core
from thales_ramanujan.checks import check_clausen_obstruction, clausen_audit_rows
from thales_ramanujan.precision import DPS, TOL, use_dps


def test_group_passes():
    for record in check_clausen_obstruction():
        assert record.status == "PASS", (record.check_id, record.symbolic_residual)


def test_beta_generating_function_parameters():
    """G_beta(z) = 2F1(1, 1; 3/2; z/4): the Clausen condition fails."""
    assert not core.clausen_condition_holds(1, 1, sp.Rational(3, 2))
    assert sp.Rational(3, 2) != 1 + 1 + sp.Rational(1, 2)


def test_beta_generating_function_value():
    with use_dps(DPS):
        z = mp.mpf(1)
        assert mp.fabs(
            core.G_beta_series(z, terms=140) - core.G_beta_closed(z)
        ) < mp.mpf(TOL)


def test_beta_moment_coefficient_identity():
    """J(n) = n! / (4^n (3/2)_n)."""
    n = core.n_sym
    lhs = sp.factorial(n) ** 2 / sp.factorial(2 * n + 1)
    rhs = sp.factorial(n) / (4 ** n * sp.rf(sp.Rational(3, 2), n))
    assert sp.simplify(sp.gammasimp(sp.simplify(lhs / rhs)) - 1) == 0


def test_arcsine_generating_function():
    with use_dps(DPS):
        z = mp.mpf(1) / 4
        assert mp.fabs(
            core.G_arcsine_series(z, terms=180) - core.G_arcsine_closed(z)
        ) < mp.mpf(TOL)
        assert mp.fabs(core.G_arcsine_closed(z) ** 2 - 1 / (1 - z)) < mp.mpf(TOL)


def test_ramanujan_side_does_satisfy_clausen():
    """Positive control: the obstruction is specific, not vacuous."""
    assert core.clausen_condition_holds(
        sp.Rational(1, 8), sp.Rational(3, 8), sp.Integer(1)
    )
    with use_dps(DPS):
        x = mp.mpf(3) / 10
        lhs = core.clausen_square_root(x) ** 2
        rhs = mp.hyper([mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4], [1, 1], x)
        assert mp.fabs(lhs - rhs) < mp.mpf(TOL)


def test_audit_table_is_deterministic():
    assert clausen_audit_rows() == clausen_audit_rows()
    holds = {row["generating_function"]: row["clausen_condition_holds"]
             for row in clausen_audit_rows()}
    assert holds["G_beta(z) = sum J(n) z^n"] == "no"
    assert holds["Ramanujan square root (positive control)"] == "yes"
