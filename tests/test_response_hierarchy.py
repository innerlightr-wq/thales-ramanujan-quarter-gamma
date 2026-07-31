"""R6/R8: the mirror-tilt diagnostic and the parameter-response hierarchy."""

import sympy as sp
from mpmath import mp

from thales_ramanujan import core
from thales_ramanujan.checks import check_mirror_tilt, check_response_hierarchy
from thales_ramanujan.precision import DPS, QUAD_DPS, QUAD_TOL, TOL, use_dps


def test_groups_pass():
    for group in (check_mirror_tilt, check_response_hierarchy):
        for record in group():
            assert record.status == "PASS", (record.check_id, record.symbolic_residual)


def test_mirror_tilt_ratio_symbolic():
    s, alpha = core.s_sym, core.alpha_sym
    ratio = sp.gammasimp(core.A_tilt_sym(s, alpha) / core.I_sym(s, alpha))
    assert sp.simplify(ratio - alpha / (s + 1)) == 0


def test_mirror_tilt_ratio_by_quadrature():
    with use_dps(QUAD_DPS):
        for s_value, a_value in (("0.75", "0.25"), ("1.5", "-0.5")):
            ratio = core.A_tilt_quad(s_value, a_value) / core.I_quad(s_value, a_value)
            expected = mp.mpf(a_value) / (mp.mpf(s_value) + 1)
            assert mp.fabs(ratio - expected) < mp.mpf(QUAD_TOL)


def test_digamma_response():
    s, alpha = core.s_sym, core.alpha_sym
    derivative = sp.diff(sp.log(core.I_sym(s, alpha)), alpha)
    expected = sp.digamma(s + alpha + 1) - sp.digamma(s - alpha + 1)
    assert sp.simplify(derivative - expected) == 0


def test_digamma_reflection():
    assert sp.simplify(
        sp.expand_func(sp.digamma(sp.Rational(3, 4)) - sp.digamma(sp.Rational(1, 4)))
        - sp.pi
    ) == 0


def test_diagonal_value():
    assert sp.simplify(
        sp.expand_func(sp.digamma(sp.Rational(3, 2)) - sp.digamma(1))
        - (2 - 2 * sp.log(2))
    ) == 0


def test_trigamma_catalan():
    with use_dps(DPS):
        assert mp.fabs(
            mp.polygamma(1, mp.mpf(1) / 4) - (mp.pi ** 2 + 8 * mp.catalan)
        ) < mp.mpf(TOL)
        assert mp.fabs(
            mp.polygamma(1, mp.mpf(1) / 4) + mp.polygamma(1, mp.mpf(3) / 4)
            - 2 * mp.pi ** 2
        ) < mp.mpf(TOL)
