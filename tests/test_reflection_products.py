"""R3: the conjugate-stratum product identity."""

import sympy as sp
from mpmath import mp

from thales_ramanujan import core
from thales_ramanujan.checks import (
    CONJUGATE_SPECIAL_CASES,
    REMOVABLE_POINTS,
    check_reflection_products,
    product_rhs_num,
    product_rhs_sym,
)
from thales_ramanujan.precision import DPS, TOL, use_dps


def test_group_passes():
    for record in check_reflection_products():
        assert record.status == "PASS", (record.check_id, record.symbolic_residual)


def test_general_identity_numerically():
    with use_dps(DPS):
        for value in ("0.11", "0.2", "0.3", "0.62", "0.9", "1.23"):
            v = mp.mpf(value)
            lhs = core.J_num(v) * core.J_num(1 - v)
            assert mp.fabs(lhs - product_rhs_num(v)) < mp.mpf(TOL)


def test_special_cases_exact():
    for s_value, exact, _ in CONJUGATE_SPECIAL_CASES:
        assert sp.simplify(
            core.J_sym(s_value) * core.J_sym(1 - s_value) - exact
        ) == 0


def test_quarter_case_is_pi_over_20():
    assert sp.simplify(
        core.J_sym(sp.Rational(1, 4)) * core.J_sym(sp.Rational(3, 4)) - sp.pi / 20
    ) == 0


def test_apparent_poles_are_removable():
    s = core.s_sym
    for point in REMOVABLE_POINTS:
        limit = sp.limit(product_rhs_sym(s), s, point)
        product = core.J_sym(point) * core.J_sym(1 - point)
        assert sp.simplify(product - limit) == 0
