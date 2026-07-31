"""The eight verification groups.

Each function returns a list of :class:`thales_ramanujan.report.Record`.
Every record states the identity being verified, a symbolic status, and (where
meaningful) an independent numerical audit with a reported absolute error
against a declared tolerance.

Nothing here asserts anything beyond the manuscript.  The Ramanujan-side
computations (``R1``, and the Clausen positive control in ``R7``) are present
only to fix the target that the Thales construction is being compared against
and demonstrably does *not* reach.
"""

from __future__ import annotations

from typing import Callable, Sequence

import sympy as sp
from mpmath import mp, mpf

from . import core
from .precision import DPS, QUAD_DPS, QUAD_TOL, TOL, use_dps
from .report import Record, check

EXACT = "0"  # tolerance string for exact rational/integer arithmetic

s = core.s_sym
alpha = core.alpha_sym
n = core.n_sym


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _to_mpf(expr) -> mpf:
    """Evaluate a SymPy expression to an mpf at the current working precision."""
    return mpf(str(sp.N(expr, mp.dps)))


def _exact_max_gap(pairs) -> mpf:
    """Maximum absolute difference over exactly-evaluated (lhs, rhs) pairs."""
    gaps = [sp.simplify(sp.nsimplify(lhs - rhs)) for lhs, rhs in pairs]
    worst = max(gaps, key=lambda g: abs(g)) if gaps else sp.Integer(0)
    return _to_mpf(sp.Abs(worst))


def _max_abs(values) -> mpf:
    return max((mp.fabs(v) for v in values), default=mp.mpf(0))


def _trig_gamma_simplify(expr):
    """Reduction strong enough to close gamma-reflection residuals.

    ``sympy.simplify`` leaves residuals of the form
    ``2 sin(x)^2 cot(x) - sin(2x)``; rewriting to cosines and expanding the
    double angle closes them.
    """
    reduced = sp.simplify(expr)
    if reduced == 0:
        return sp.Integer(0)
    return sp.simplify(sp.expand_trig(reduced.rewrite(sp.cos)))


def _quarter_reflection_simplify(expr):
    """Reduction that applies Euler reflection at the quarter residues.

    ``Gamma(1/4) Gamma(3/4) = pi / sin(pi/4) = pi sqrt(2)``, so eliminating
    ``Gamma(3/4)`` closes residuals that mix the two conjugate quarter strata.
    """
    reduced = sp.simplify(expr)
    if reduced == 0:
        return sp.Integer(0)
    reflected = reduced.subs(
        sp.gamma(sp.Rational(3, 4)),
        sp.pi * sp.sqrt(2) / sp.gamma(sp.Rational(1, 4)),
    )
    return sp.simplify(reflected)


def _gammasimp_ratio(lhs, rhs):
    """Residual ``gammasimp(lhs/rhs) - 1``, for coefficient identities that
    SymPy closes on the ratio but not on the difference."""
    return sp.simplify(sp.gammasimp(sp.simplify(lhs / rhs)) - 1)


# ---------------------------------------------------------------------------
# R1. Ramanujan coefficient identity and its generating function
# ---------------------------------------------------------------------------

def check_ramanujan_coefficients() -> list[Record]:
    group = "R1_ramanujan_coefficients"
    records: list[Record] = []

    # R1.1 -- A_n = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3
    #
    # Proved by induction: the base case n = 0 and the equality of the two
    # consecutive-term ratios.  SymPy closes both, whereas it does not close
    # the Gauss multiplication formula in a single step for symbolic n.
    ratio_lhs = sp.simplify(core.A_factorial_sym(n + 1) / core.A_factorial_sym(n))
    ratio_rhs = (
        256
        * (sp.Rational(1, 4) + n)
        * (sp.Rational(1, 2) + n)
        * (sp.Rational(3, 4) + n)
        / (1 + n) ** 3
    )
    base_case = core.A_factorial_sym(0) - core.A_pochhammer_sym(0)
    residual = sp.simplify(sp.expand(ratio_lhs - ratio_rhs)) + sp.simplify(base_case)
    exact_gap = _exact_max_gap(
        [(core.A_factorial_sym(k), core.A_pochhammer_sym(k)) for k in range(0, 31)]
    )
    records.append(
        check(
            "R1.1",
            group,
            "A_n = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3, exactly (n = 0..30 audited)",
            exact_form="A_n = (4n)!/(n!)^4",
            symbolic_residual=residual,
            symbolic_note="verified (induction: base case n=0 and term-ratio equality)",
            lhs=exact_gap,
            rhs=mp.mpf(0),
            tolerance=EXACT,
            dps=DPS,
        )
    )

    # R1.2 -- generating function F(z) = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z)
    with use_dps(DPS):
        z0 = mpf(1) / 1000
        series = core.F_series(z0, terms=140)
        closed = core.F_hypergeometric(z0)
    records.append(
        check(
            "R1.2",
            group,
            "F(z) = sum A_n z^n = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z) at z = 1/1000",
            exact_form="F(z) = 3F2(1/4,1/2,3/4; 1,1; 256 z)",
            symbolic_note="documented (a restatement of R1.1)",
            lhs=series,
            rhs=closed,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R1.3 -- the target itself, for orientation only
    with use_dps(DPS):
        approx = core.ramanujan_partial_sum(terms=60)
        target = 1 / mp.pi
    records.append(
        check(
            "R1.3",
            group,
            "Ramanujan's series (manuscript eq. (1)) sums to 1/pi (orientation anchor)",
            exact_form="1/pi = (2 sqrt(2)/9801) sum A_n (1103 + 26390 n) / 396^(4n)",
            symbolic_note="not attempted (requires modular theory; outside repository scope)",
            lhs=approx,
            rhs=target,
            tolerance=TOL,
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R2. Thales beta family and its gamma strata
# ---------------------------------------------------------------------------

def check_beta_gamma_values() -> list[Record]:
    group = "R2_beta_gamma_values"
    records: list[Record] = []

    # R2.1 -- B(s+1,s+1) = Gamma(s+1)^2 / Gamma(2s+2), with the integral itself
    #         audited by direct quadrature.
    with use_dps(QUAD_DPS):
        quad_errors = [
            mp.fabs(core.J_quad(sv) - core.J_num(sv))
            for sv in ("-0.5", "-0.25", "0.25", "0.5", "0.75", "1", "1.5", "3")
        ]
        worst_quad = _max_abs(quad_errors)
    records.append(
        check(
            "R2.1",
            group,
            "J(s) = int_0^1 [a(1-a)]^s da = B(s+1,s+1) = Gamma(s+1)^2/Gamma(2s+2) "
            "(quadrature audit over 8 values of s)",
            exact_form="J(s) = Gamma(s+1)^2 / Gamma(2s+2)",
            symbolic_residual=core.J_beta_sym(s) - core.J_sym(s),
            lhs=worst_quad,
            rhs=mp.mpf(0),
            tolerance=QUAD_TOL,
            dps=QUAD_DPS,
        )
    )

    # R2.2 -- stratum recurrence J(s+1)/J(s) = (s+1) / (2(2s+3))
    records.append(
        check(
            "R2.2",
            group,
            "J(s+1)/J(s) = (s+1) / (2(2s+3)) is rational in s (the stratum recurrence)",
            exact_form="J(s+1)/J(s) = (s+1)/(2(2s+3))",
            symbolic_residual=sp.simplify(
                core.J_sym(s + 1) / core.J_sym(s) - (s + 1) / (2 * (2 * s + 3))
            ),
        )
    )

    # R2.3 -- quarter stratum
    with use_dps(DPS):
        j_quarter = core.J_num(mpf(1) / 4)
        j_quarter_gamma = mp.gamma(mp.mpf(1) / 4) ** 2 / (12 * mp.sqrt(mp.pi))
    records.append(
        check(
            "R2.3",
            group,
            "J(1/4) = Gamma(1/4)^2 / (12 sqrt(pi))",
            exact_form="J(1/4) = Gamma(1/4)^2 / (12 sqrt(pi))",
            symbolic_residual=core.J_sym(sp.Rational(1, 4))
            - sp.gamma(sp.Rational(1, 4)) ** 2 / (12 * sp.sqrt(sp.pi)),
            lhs=j_quarter,
            rhs=j_quarter_gamma,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R2.4 -- conjugate quarter stratum
    with use_dps(DPS):
        j_three_quarter = core.J_num(mpf(3) / 4)
        j_three_quarter_gamma = 3 * mp.gamma(mp.mpf(3) / 4) ** 2 / (10 * mp.sqrt(mp.pi))
    records.append(
        check(
            "R2.4",
            group,
            "J(3/4) = 3 Gamma(3/4)^2 / (10 sqrt(pi))",
            exact_form="J(3/4) = 3 Gamma(3/4)^2 / (10 sqrt(pi))",
            symbolic_residual=core.J_sym(sp.Rational(3, 4))
            - 3 * sp.gamma(sp.Rational(3, 4)) ** 2 / (10 * sp.sqrt(sp.pi)),
            lhs=j_three_quarter,
            rhs=j_three_quarter_gamma,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R2.5 -- half-integer stratum
    with use_dps(DPS):
        j_half = core.J_num(mpf(1) / 2)
        pi_over_8 = mp.pi / 8
    records.append(
        check(
            "R2.5",
            group,
            "J(1/2) = pi/8",
            exact_form="J(1/2) = pi/8",
            symbolic_residual=core.J_sym(sp.Rational(1, 2)) - sp.pi / 8,
            lhs=j_half,
            rhs=pi_over_8,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R2.6 -- the elementary strata quoted in the manuscript
    strata_residual = (
        sp.simplify(core.J_sym(0) - 1)
        + sp.simplify(core.J_sym(1) - sp.Rational(1, 6))
        + sp.simplify(core.J_sym(sp.Rational(-1, 2)) - sp.pi)
    )
    records.append(
        check(
            "R2.6",
            group,
            "Elementary strata: J(0) = 1, J(1) = 1/6, J(-1/2) = pi",
            exact_form="J(0)=1; J(1)=1/6; J(-1/2)=pi",
            symbolic_residual=strata_residual,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R3. Conjugate-stratum product identity
# ---------------------------------------------------------------------------

#: (s, exact product) for the three special cases recorded in the manuscript
CONJUGATE_SPECIAL_CASES = (
    (sp.Rational(1, 4), sp.pi / 20, "J(1/4) J(3/4) = pi/20"),
    (sp.Rational(1, 3), sp.pi * sp.sqrt(3) / 35, "J(1/3) J(2/3) = pi sqrt(3)/35"),
    (sp.Rational(1, 6), 15 * sp.pi * sp.sqrt(3) / 512, "J(1/6) J(5/6) = 15 pi sqrt(3)/512"),
)

#: apparent poles of the right-hand side of the product identity that are in
#: fact removable, because cot(pi s) vanishes there
REMOVABLE_POINTS = (sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(3, 2))


def product_rhs_sym(sv):
    """pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s)) (symbolic)."""
    return (
        sp.pi * sv * (1 - sv) * sp.cot(sp.pi * sv)
        / (2 * (2 * sv + 1) * (1 - 2 * sv) * (3 - 2 * sv))
    )


def product_rhs_num(sv):
    """The same, evaluated numerically."""
    sv = mpf(sv)
    return (
        mp.pi * sv * (1 - sv) * mp.cot(mp.pi * sv)
        / (2 * (2 * sv + 1) * (1 - 2 * sv) * (3 - 2 * sv))
    )


def check_reflection_products() -> list[Record]:
    group = "R3_reflection_products"
    records: list[Record] = []

    # R3.1 -- the general identity
    with use_dps(DPS):
        sweep = ["0.11", "0.2", "0.3", "0.37", "0.62", "0.8", "0.9", "1.23"]
        sweep_errors = [
            mp.fabs(core.J_num(mpf(v)) * core.J_num(1 - mpf(v)) - product_rhs_num(v))
            for v in sweep
        ]
        worst = _max_abs(sweep_errors)
    records.append(
        check(
            "R3.1",
            group,
            "J(s) J(1-s) = pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s)) "
            "(8-point numerical sweep)",
            exact_form="J(s)J(1-s) = pi s(1-s)cot(pi s) / (2(2s+1)(1-2s)(3-2s))",
            symbolic_residual=core.J_sym(s) * core.J_sym(1 - s) - product_rhs_sym(s),
            simplifier=_trig_gamma_simplify,
            lhs=worst,
            rhs=mp.mpf(0),
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R3.2 - R3.4 -- the three special cases
    for index, (sv, exact, label) in enumerate(CONJUGATE_SPECIAL_CASES, start=2):
        with use_dps(DPS):
            lhs = core.J_num(_to_mpf(sv)) * core.J_num(_to_mpf(1 - sv))
            rhs = _to_mpf(exact)
        records.append(
            check(
                f"R3.{index}",
                group,
                label,
                exact_form=sp.sstr(exact),
                symbolic_residual=sp.simplify(
                    core.J_sym(sv) * core.J_sym(1 - sv) - exact
                ),
                lhs=lhs,
                rhs=rhs,
                tolerance=TOL,
                dps=DPS,
            )
        )

    # R3.5 -- the apparent poles at s = -1/2, 1/2, 3/2 are removable
    removable_errors = []
    with use_dps(DPS):
        for sv in REMOVABLE_POINTS:
            limit = sp.limit(product_rhs_sym(s), s, sv)
            product = core.J_sym(sv) * core.J_sym(1 - sv)
            removable_errors.append(mp.fabs(_to_mpf(sp.simplify(product - limit))))
        worst_removable = _max_abs(removable_errors)
    records.append(
        check(
            "R3.5",
            group,
            "Apparent poles of the right-hand side at s = -1/2, 1/2, 3/2 are removable "
            "(cot(pi s) vanishes there); the limits agree with J(s)J(1-s)",
            exact_form="limits at s in {-1/2, 1/2, 3/2} agree with J(s)J(1-s)",
            symbolic_note="verified (SymPy limits at each point)",
            lhs=worst_removable,
            rhs=mp.mpf(0),
            tolerance=TOL,
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R4. Lemniscatic period locus
# ---------------------------------------------------------------------------

def check_lemniscatic_relation() -> list[Record]:
    group = "R4_lemniscatic_relation"
    records: list[Record] = []

    with use_dps(DPS):
        k_val = core.K_lemniscatic()
        k_gamma = core.K_lemniscatic_gamma()
    records.append(
        check(
            "R4.1",
            group,
            "K(1/sqrt(2)) = Gamma(1/4)^2 / (4 sqrt(pi)) (classical lemniscatic evaluation)",
            exact_form="K(1/sqrt(2)) = Gamma(1/4)^2 / (4 sqrt(pi))",
            symbolic_note="classical (Whittaker-Watson; Borwein-Borwein); audited numerically",
            lhs=k_val,
            rhs=k_gamma,
            tolerance=TOL,
            dps=DPS,
        )
    )

    with use_dps(DPS):
        lhs = core.J_num(mpf(1) / 4)
        rhs = k_val / 3
    records.append(
        check(
            "R4.2",
            group,
            "J(1/4) = (1/3) K(1/sqrt(2))",
            exact_form="J(1/4) = (1/3) K(1/sqrt(2))",
            symbolic_residual=sp.simplify(
                core.J_sym(sp.Rational(1, 4))
                - sp.gamma(sp.Rational(1, 4)) ** 2 / (12 * sp.sqrt(sp.pi))
            ),
            symbolic_note="verified against the gamma form of K(1/sqrt(2))",
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    with use_dps(DPS):
        lhs = core.J_num(mpf(3) / 4)
        rhs = 3 * mp.pi / (20 * k_val)
    records.append(
        check(
            "R4.3",
            group,
            "J(3/4) = 3 pi / (20 K(1/sqrt(2))), i.e. the conjugate stratum is "
            "proportional to the reciprocal period",
            exact_form="J(3/4) = 3 pi / (20 K(1/sqrt(2)))",
            # 3 pi / (20 K(1/sqrt(2))) with K = Gamma(1/4)^2/(4 sqrt(pi))
            # equals 3 pi^(3/2) / (5 Gamma(1/4)^2).
            symbolic_residual=core.J_sym(sp.Rational(3, 4))
            - 3 * sp.pi ** sp.Rational(3, 2) / (5 * sp.gamma(sp.Rational(1, 4)) ** 2),
            simplifier=_quarter_reflection_simplify,
            symbolic_note="verified via Euler reflection (equivalently, from R3.2 and R4.2)",
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R5. Arcsine witness and central-binomial moments
# ---------------------------------------------------------------------------

def check_arcsine_moments() -> list[Record]:
    group = "R5_arcsine_moments"
    records: list[Record] = []

    # R5.1 -- the moment identity
    moment_sym = sp.beta(n + sp.Rational(1, 2), sp.Rational(1, 2)) / sp.pi
    with use_dps(QUAD_DPS):
        quad_errors = [
            mp.fabs(core.arcsine_moment_quad(k) - core.arcsine_moment_closed(k))
            for k in range(0, 5)
        ]
        worst = _max_abs(quad_errors)
    records.append(
        check(
            "R5.1",
            group,
            "(1/pi) int_0^1 x^n / sqrt(x(1-x)) dx = 4^(-n) C(2n,n) "
            "(quadrature audit for n = 0..4)",
            exact_form="m_n = 4^(-n) C(2n,n)",
            symbolic_residual=sp.simplify(
                sp.combsimp(moment_sym - sp.binomial(2 * n, n) / 4 ** n)
            ),
            lhs=worst,
            rhs=mp.mpf(0),
            tolerance=QUAD_TOL,
            dps=QUAD_DPS,
        )
    )

    # R5.2 -- central-binomial decomposition of Ramanujan's coefficient
    exact_gap = _exact_max_gap(
        [(core.A_factorial_sym(k), core.A_binomial_sym(k)) for k in range(0, 31)]
    )
    records.append(
        check(
            "R5.2",
            group,
            "(4n)!/(n!)^4 = C(4n,2n) C(2n,n)^2, exactly (n = 0..30 audited)",
            exact_form="A_n = C(4n,2n) C(2n,n)^2",
            symbolic_residual=sp.simplify(
                sp.combsimp(core.A_binomial_sym(n) - core.A_factorial_sym(n))
            ),
            lhs=exact_gap,
            rhs=mp.mpf(0),
            tolerance=EXACT,
            dps=DPS,
        )
    )

    # R5.3 -- the "witness" statement in its exact form: each factor of R5.2 is
    #         an arcsine moment, so A_n = 256^n m_{2n} m_n^2.  This is a
    #         restatement of R5.1 and R5.2 together, not a new claim, and in
    #         particular it is NOT a generating mechanism: a product of three
    #         moments is not itself a moment.
    m = lambda k: sp.binomial(2 * k, k) / 4 ** k
    witness_sym = 256 ** n * m(2 * n) * m(n) ** 2
    exact_gap = _exact_max_gap(
        [(core.A_factorial_sym(k), witness_sym.subs(n, k)) for k in range(0, 21)]
    )
    records.append(
        check(
            "R5.3",
            group,
            "A_n = 256^n m_(2n) m_n^2 with m_k the arcsine moments of R5.1 "
            "(restatement of R5.1 + R5.2; not a generating mechanism)",
            exact_form="A_n = 256^n m_(2n) m_n^2",
            symbolic_residual=sp.simplify(
                sp.combsimp(witness_sym - core.A_factorial_sym(n))
            ),
            lhs=exact_gap,
            rhs=mp.mpf(0),
            tolerance=EXACT,
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R6. Mirror-tilt diagnostic and the derivative boundary
# ---------------------------------------------------------------------------

MIRROR_TILT_POINTS = (("0.75", "0.25"), ("1.5", "-0.5"), ("0.25", "0.125"), ("3", "1.5"))


def check_mirror_tilt() -> list[Record]:
    group = "R6_mirror_tilt"
    records: list[Record] = []

    # R6.1 -- the algebraic tilt ratio
    with use_dps(QUAD_DPS):
        errors = [
            mp.fabs(
                core.A_tilt_quad(sv, av) / core.I_quad(sv, av)
                - mpf(av) / (mpf(sv) + 1)
            )
            for sv, av in MIRROR_TILT_POINTS
        ]
        worst = _max_abs(errors)
    records.append(
        check(
            "R6.1",
            group,
            "A_s(alpha) / I_s(alpha) = alpha / (s+1) (quadrature audit at 4 points)",
            exact_form="A_s(alpha)/I_s(alpha) = alpha/(s+1)",
            symbolic_residual=sp.simplify(
                sp.gammasimp(core.A_tilt_sym(s, alpha) / core.I_sym(s, alpha))
                - alpha / (s + 1)
            ),
            lhs=worst,
            rhs=mp.mpf(0),
            tolerance=QUAD_TOL,
            dps=QUAD_DPS,
        )
    )

    # R6.2 -- the ratio is 2 E[a] - 1 for a ~ Beta(s+alpha+1, s-alpha+1)
    beta_mean = (s + alpha + 1) / (2 * s + 2)
    records.append(
        check(
            "R6.2",
            group,
            "The tilt ratio equals 2 E[a] - 1 for a ~ Beta(s+alpha+1, s-alpha+1); "
            "it is algebraic in (s, alpha)",
            exact_form="2 E[a] - 1 = alpha/(s+1), E[a] = (s+alpha+1)/(2s+2)",
            symbolic_residual=sp.simplify(2 * beta_mean - 1 - alpha / (s + 1)),
        )
    )

    # R6.3 -- the honest (transcendental) parameter response
    with use_dps(DPS):
        errors = []
        for sv, av in MIRROR_TILT_POINTS:
            errors.append(
                mp.fabs(
                    core.log_I_alpha_derivative(sv, av) - core.digamma_response(sv, av)
                )
            )
        worst = _max_abs(errors)
    records.append(
        check(
            "R6.3",
            group,
            "d/d(alpha) log I_s(alpha) = psi(s+alpha+1) - psi(s-alpha+1) "
            "(numerical differentiation at 4 points)",
            exact_form="d_alpha log I_s(alpha) = psi(s+alpha+1) - psi(s-alpha+1)",
            symbolic_residual=sp.simplify(
                sp.diff(sp.log(core.I_sym(s, alpha)), alpha)
                - (sp.digamma(s + alpha + 1) - sp.digamma(s - alpha + 1))
            ),
            lhs=worst,
            rhs=mp.mpf(0),
            tolerance="1e-40",
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R7. Clausen obstruction
# ---------------------------------------------------------------------------

#: rows of results/clausen_parameter_audit.csv
CLAUSEN_CANDIDATES = (
    {
        "generating_function": "G_beta(z) = sum J(n) z^n",
        "hypergeometric_form": "2F1(1, 1; 3/2; z/4)",
        "type": "2F1",
        "a": sp.Integer(1),
        "b": sp.Integer(1),
        "c": sp.Rational(3, 2),
        "note": "integer Thales beta moments; Clausen condition fails",
    },
    {
        "generating_function": "G_arc(z) = sum 4^(-n) C(2n,n) z^n",
        "hypergeometric_form": "1F0(1/2; ; z) = (1-z)^(-1/2)",
        "type": "1F0",
        "a": sp.Rational(1, 2),
        "b": None,
        "c": None,
        "note": "no second numerator and no denominator parameter; condition does not apply",
    },
    {
        "generating_function": "G_arc(z)^2",
        "hypergeometric_form": "1F0(1; ; z) = (1-z)^(-1)",
        "type": "1F0",
        "a": sp.Integer(1),
        "b": None,
        "c": None,
        "note": "degenerate square; carries no modular content",
    },
    {
        "generating_function": "Ramanujan square root (positive control)",
        "hypergeometric_form": "2F1(1/8, 3/8; 1; 256 z)",
        "type": "2F1",
        "a": sp.Rational(1, 8),
        "b": sp.Rational(3, 8),
        "c": sp.Integer(1),
        "note": "Clausen condition holds; square is the 3F2 of manuscript eq. (3)",
    },
)


def clausen_audit_rows() -> list[dict]:
    """Deterministic rows for the Clausen parameter audit table."""
    rows = []
    for candidate in CLAUSEN_CANDIDATES:
        a, b, c = candidate["a"], candidate["b"], candidate["c"]
        if candidate["type"] == "2F1":
            required = sp.nsimplify(a + b + sp.Rational(1, 2))
            holds = core.clausen_condition_holds(a, b, c)
            applies = "yes"
        else:
            required = None
            holds = False
            applies = "no"
        rows.append(
            {
                "generating_function": candidate["generating_function"],
                "hypergeometric_form": candidate["hypergeometric_form"],
                "a": sp.sstr(a) if a is not None else "",
                "b": sp.sstr(b) if b is not None else "",
                "c": sp.sstr(c) if c is not None else "",
                "clausen_condition_applies": applies,
                "required_c": sp.sstr(required) if required is not None else "",
                "clausen_condition_holds": "yes" if holds else "no",
                "note": candidate["note"],
            }
        )
    return rows


def check_clausen_obstruction() -> list[Record]:
    group = "R7_clausen_obstruction"
    records: list[Record] = []

    # R7.1 -- G_beta(z) = 2F1(1, 1; 3/2; z/4)
    j_int = sp.factorial(n) ** 2 / sp.factorial(2 * n + 1)
    target = sp.factorial(n) / (4 ** n * sp.rf(sp.Rational(3, 2), n))
    with use_dps(DPS):
        z0 = mpf(1)
        series = core.G_beta_series(z0, terms=140)
        closed = core.G_beta_closed(z0)
    records.append(
        check(
            "R7.1",
            group,
            "G_beta(z) = sum J(n) z^n = 2F1(1, 1; 3/2; z/4), checked at z = 1",
            exact_form="G_beta(z) = 2F1(1,1;3/2;z/4)",
            symbolic_residual=_gammasimp_ratio(j_int, target),
            symbolic_note="verified (coefficient identity J(n) = n!/(4^n (3/2)_n))",
            lhs=series,
            rhs=closed,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R7.2 -- the Clausen condition fails for those parameters
    beta_condition = core.clausen_condition_holds(1, 1, sp.Rational(3, 2))
    records.append(
        check(
            "R7.2",
            group,
            "G_beta parameters (a,b,c) = (1, 1, 3/2) do NOT satisfy c = a + b + 1/2 "
            "(which would require c = 5/2)",
            exact_form="c = 3/2 != 5/2 = a + b + 1/2",
            symbolic_residual=sp.Integer(0) if not beta_condition else sp.Integer(1),
            symbolic_note="verified (exact rational parameter arithmetic)",
        )
    )

    # R7.3 -- G_arc(z) = (1-z)^(-1/2) = 1F0(1/2; ; z)
    arc_coeff = sp.binomial(2 * n, n) / 4 ** n
    arc_target = sp.rf(sp.Rational(1, 2), n) / sp.factorial(n)
    with use_dps(DPS):
        z0 = mpf(1) / 4
        series = core.G_arcsine_series(z0, terms=180)
        closed = core.G_arcsine_closed(z0)
    records.append(
        check(
            "R7.3",
            group,
            "G_arc(z) = sum 4^(-n) C(2n,n) z^n = (1-z)^(-1/2) = 1F0(1/2; ; z), "
            "checked at z = 1/4",
            exact_form="G_arc(z) = (1-z)^(-1/2)",
            symbolic_residual=sp.simplify(sp.combsimp(arc_coeff - arc_target)),
            lhs=series,
            rhs=closed,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R7.4 -- the square of G_arc is degenerate
    with use_dps(DPS):
        z0 = mpf(1) / 4
        lhs = core.G_arcsine_closed(z0) ** 2
        rhs = 1 / (1 - z0)
    records.append(
        check(
            "R7.4",
            group,
            "G_arc(z)^2 = (1-z)^(-1) = 1F0(1; ; z): a degenerate square with no "
            "Gauss-type 2F1 square-root structure",
            exact_form="G_arc(z)^2 = (1-z)^(-1)",
            symbolic_note="verified (elementary)",
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R7.5 -- positive control: the Ramanujan side DOES satisfy Clausen
    holds = core.clausen_condition_holds(
        sp.Rational(1, 8), sp.Rational(3, 8), sp.Integer(1)
    )
    with use_dps(DPS):
        x0 = mpf(3) / 10
        lhs = core.clausen_square_root(x0) ** 2
        rhs = mp.hyper([mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4], [1, 1], x0)
    records.append(
        check(
            "R7.5",
            group,
            "Positive control: (a,b,c) = (1/8, 3/8, 1) satisfies c = a + b + 1/2, and "
            "2F1(1/8,3/8;1;x)^2 = 3F2(1/4,1/2,3/4;1,1;x)",
            exact_form="Clausen square, manuscript eq. (4)",
            symbolic_residual=sp.Integer(0) if holds else sp.Integer(1),
            symbolic_note="verified (exact rational parameter arithmetic)",
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# R8. Parameter-response hierarchy
# ---------------------------------------------------------------------------

def check_response_hierarchy() -> list[Record]:
    group = "R8_response_hierarchy"
    records: list[Record] = []

    # R8.1 -- first-order response at the conjugate quarter points
    with use_dps(DPS):
        lhs = mp.digamma(mp.mpf(3) / 4) - mp.digamma(mp.mpf(1) / 4)
        rhs = mp.pi
    records.append(
        check(
            "R8.1",
            group,
            "psi(3/4) - psi(1/4) = pi (Gauss digamma theorem / reflection)",
            exact_form="psi(3/4) - psi(1/4) = pi",
            symbolic_residual=sp.simplify(
                sp.expand_func(
                    sp.digamma(sp.Rational(3, 4)) - sp.digamma(sp.Rational(1, 4))
                )
                - sp.pi
            ),
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R8.2 -- diagonal value
    with use_dps(DPS):
        lhs = mp.digamma(mp.mpf(3) / 2) - mp.digamma(1)
        rhs = 2 - 2 * mp.log(2)
    records.append(
        check(
            "R8.2",
            group,
            "psi(3/2) - psi(1) = 2 - 2 log 2 (diagonal s = alpha = 1/4)",
            exact_form="psi(3/2) - psi(1) = 2 - 2 log 2",
            symbolic_residual=sp.simplify(
                sp.expand_func(sp.digamma(sp.Rational(3, 2)) - sp.digamma(1))
                - (2 - 2 * sp.log(2))
            ),
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R8.3 -- Catalan's constant enters only at second order
    with use_dps(DPS):
        lhs = mp.polygamma(1, mp.mpf(1) / 4)
        rhs = mp.pi ** 2 + 8 * mp.catalan
    records.append(
        check(
            "R8.3",
            group,
            "psi'(1/4) = pi^2 + 8 G (Catalan's constant enters at second order)",
            exact_form="psi'(1/4) = pi^2 + 8 G",
            symbolic_note=(
                "classical; SymPy does not reduce polygamma(1, 1/4) in closed form, "
                "so this is audited numerically (see also R8.4)"
            ),
            lhs=lhs,
            rhs=rhs,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R8.4 -- the reflection half of R8.3, which does close symbolically in the
    #         sense that psi'(1/4) + psi'(3/4) = 2 pi^2; combined with the
    #         Dirichlet-beta definition psi'(1/4) - psi'(3/4) = 16 G this gives
    #         R8.3.
    with use_dps(DPS):
        sum_err = mp.polygamma(1, mp.mpf(1) / 4) + mp.polygamma(1, mp.mpf(3) / 4)
        diff_lhs = mp.polygamma(1, mp.mpf(1) / 4) - mp.polygamma(1, mp.mpf(3) / 4)
    records.append(
        check(
            "R8.4",
            group,
            "Trigamma reflection psi'(1/4) + psi'(3/4) = 2 pi^2, and "
            "psi'(1/4) - psi'(3/4) = 16 G; together these give R8.3",
            exact_form="psi'(1/4) + psi'(3/4) = 2 pi^2",
            symbolic_note="audited numerically (reflection formula for psi')",
            lhs=sum_err,
            rhs=2 * mp.pi ** 2,
            tolerance=TOL,
            dps=DPS,
        )
    )
    records.append(
        check(
            "R8.5",
            group,
            "psi'(1/4) - psi'(3/4) = 16 G (Dirichlet beta at 2)",
            exact_form="psi'(1/4) - psi'(3/4) = 16 G",
            symbolic_note="audited numerically",
            lhs=diff_lhs,
            rhs=16 * mp.catalan,
            tolerance=TOL,
            dps=DPS,
        )
    )

    # R8.6 -- these values are realised inside the beta family itself: at
    #         (s, alpha) = (-1/2, -1/4) the second-order response is
    #         psi'(1/4) + psi'(3/4).
    with use_dps(DPS):
        sv, av = mpf("-0.5"), mpf("-0.25")
        second = mp.diff(lambda t: mp.log(core.I_num(sv, t)), av, 2)
        expected = mp.polygamma(1, sv + av + 1) + mp.polygamma(1, sv - av + 1)
    records.append(
        check(
            "R8.6",
            group,
            "Second-order response d^2/d(alpha)^2 log I_s(alpha) = "
            "psi'(s+alpha+1) + psi'(s-alpha+1), realised at (s, alpha) = (-1/2, -1/4) "
            "where it equals psi'(1/4) + psi'(3/4)",
            exact_form="d^2_alpha log I_s(alpha) = psi'(s+alpha+1) + psi'(s-alpha+1)",
            symbolic_residual=sp.simplify(
                sp.diff(sp.log(core.I_sym(s, alpha)), alpha, 2)
                - (
                    sp.polygamma(1, s + alpha + 1)
                    + sp.polygamma(1, s - alpha + 1)
                )
            ),
            lhs=second,
            rhs=expected,
            tolerance="1e-30",
            dps=DPS,
        )
    )

    return records


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------

#: ordered registry: (group key, script filename, runner)
REGISTRY: tuple[tuple[str, str, Callable[[], list[Record]]], ...] = (
    ("R1_ramanujan_coefficients", "verify_ramanujan_coefficients.py", check_ramanujan_coefficients),
    ("R2_beta_gamma_values", "verify_beta_gamma_values.py", check_beta_gamma_values),
    ("R3_reflection_products", "verify_reflection_products.py", check_reflection_products),
    ("R4_lemniscatic_relation", "verify_lemniscatic_relation.py", check_lemniscatic_relation),
    ("R5_arcsine_moments", "verify_arcsine_moments.py", check_arcsine_moments),
    ("R6_mirror_tilt", "verify_mirror_tilt.py", check_mirror_tilt),
    ("R7_clausen_obstruction", "verify_clausen_obstruction.py", check_clausen_obstruction),
    ("R8_response_hierarchy", "verify_response_hierarchy.py", check_response_hierarchy),
)


def run_all_checks() -> list[Record]:
    """Run every group in registry order."""
    records: list[Record] = []
    for _, _, runner in REGISTRY:
        records.extend(runner())
    return records
