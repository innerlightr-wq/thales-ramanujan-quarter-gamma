"""R2/R4/R5: the Thales beta family, its gamma strata, the lemniscatic
period, and the arcsine moments."""

import sympy as sp
from mpmath import mp

from thales_ramanujan import core
from thales_ramanujan.checks import (
    check_arcsine_moments,
    check_beta_gamma_values,
    check_lemniscatic_relation,
)
from thales_ramanujan.precision import DPS, QUAD_DPS, QUAD_TOL, TOL, use_dps


def test_groups_pass():
    for group in (
        check_beta_gamma_values,
        check_lemniscatic_relation,
        check_arcsine_moments,
    ):
        for record in group():
            assert record.status == "PASS", (record.check_id, record.symbolic_residual)


def test_beta_equals_gamma_form():
    s = core.s_sym
    assert sp.simplify(core.J_beta_sym(s) - core.J_sym(s)) == 0


def test_quarter_strata_exact():
    assert sp.simplify(
        core.J_sym(sp.Rational(1, 4))
        - sp.gamma(sp.Rational(1, 4)) ** 2 / (12 * sp.sqrt(sp.pi))
    ) == 0
    assert sp.simplify(
        core.J_sym(sp.Rational(3, 4))
        - 3 * sp.gamma(sp.Rational(3, 4)) ** 2 / (10 * sp.sqrt(sp.pi))
    ) == 0
    assert sp.simplify(core.J_sym(sp.Rational(1, 2)) - sp.pi / 8) == 0


def test_integral_matches_beta_by_quadrature():
    with use_dps(QUAD_DPS):
        for value in ("-0.5", "0.25", "0.75", "2"):
            err = mp.fabs(core.J_quad(value) - core.J_num(value))
            assert err < mp.mpf(QUAD_TOL)


def test_lemniscatic_relation():
    with use_dps(DPS):
        assert mp.fabs(core.K_lemniscatic() - core.K_lemniscatic_gamma()) < mp.mpf(TOL)
        assert mp.fabs(
            core.J_num(mp.mpf(1) / 4) - core.K_lemniscatic() / 3
        ) < mp.mpf(TOL)


def test_arcsine_moments_and_decomposition():
    n = core.n_sym
    moment = sp.beta(n + sp.Rational(1, 2), sp.Rational(1, 2)) / sp.pi
    assert sp.simplify(
        sp.combsimp(moment - sp.binomial(2 * n, n) / 4 ** n)
    ) == 0
    for k in range(0, 21):
        assert core.A_exact(k) == int(
            sp.binomial(4 * k, 2 * k) * sp.binomial(2 * k, k) ** 2
        )
