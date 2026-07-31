# Computational methods

## Verification protocol

Every identity produces exactly one record, and every record carries the same
five things:

1. **What is being verified** — a plain statement of the identity, stored in
   the record's `description` and `exact_form` fields.
2. **A symbolic reduction where practical** — a SymPy residual that must reduce
   to `0` in exact arithmetic.
3. **An independent numerical check** — an arbitrary-precision mpmath
   evaluation of the two sides, computed by a different route from the symbolic
   one wherever possible (quadrature of the defining integral against a closed
   beta form; a truncated power series against a hypergeometric evaluation;
   numerical differentiation against a digamma expression).
4. **A reported absolute error.**
5. **A declared tolerance, and a clear failure if the error exceeds it.**

Floating-point equality is never used anywhere in this repository. Comparisons
are either exact (rational or integer arithmetic in SymPy) or bounded by a
declared absolute tolerance.

## Symbolic reductions

Most residuals close under plain `sympy.simplify`. Three do not, and each is
handled by a named, documented reduction in `thales_ramanujan/checks.py` rather
than by weakening the check:

| Situation | Reduction | Used by |
|---|---|---|
| Gamma reflection leaves a residual of the form `2 sin(x)^2 cot(x) - sin(2x)` | `_trig_gamma_simplify`: rewrite to cosines, then `expand_trig` | `R3.1` |
| Residual mixes `Gamma(1/4)` and `Gamma(3/4)` | `_quarter_reflection_simplify`: substitute `Gamma(3/4) = pi sqrt(2) / Gamma(1/4)`, which is Euler reflection at `1/4` | `R4.3` |
| Factorial/Pochhammer identity that SymPy closes on the ratio but not the difference | `_gammasimp_ratio`: check `gammasimp(lhs/rhs) - 1 == 0` | `R7.1` |

Two identities are not reduced symbolically by SymPy at all and say so in their
record rather than pretending otherwise:

* `R1.1`, the Gauss multiplication formula at `m = 4`, is instead verified by
  induction — SymPy closes both the base case and the equality of consecutive
  term ratios — and independently by exact integer arithmetic for `n = 0..30`.
* `R8.3`, `psi'(1/4) = pi^2 + 8 G`, has no closed-form reduction in SymPy for
  `polygamma(1, 1/4)`. Its record states this, it is audited numerically to 60
  digits, and the supporting decomposition into the reflection `R8.4` and the
  Dirichlet-beta value `R8.5` is recorded separately.

The `symbolic_status` column of `results/exact_identity_summary.csv` reports
exactly which of these applies to each identity. Reading that column is the
fastest way to see what is proved symbolically and what is audited numerically.

## Precision policy

Two regimes, declared in `thales_ramanujan/precision.py` and reported per check
in `results/high_precision_checks.csv`:

| Regime | Digits | Tolerance | Used for |
|---|---|---|---|
| Closed form | `DPS = 60` | `1e-50` | gamma, beta, digamma, polygamma, elliptic integrals, hypergeometric evaluations, truncated series |
| Quadrature | `QUAD_DPS = 40` | `1e-20` | direct tanh-sinh integration of the defining integrals |

The regimes are separated deliberately. Several Thales integrands carry
algebraic endpoint singularities — the arcsine density `[a(1-a)]^(-1/2)` being
the extreme case — and although tanh-sinh quadrature handles them well, it
saturates around `1e-22` at 40 digits rather than reaching full working
precision. Reporting a single global tolerance would either overstate the
quadrature or needlessly weaken the closed-form audit. Observed quadrature
errors are typically two or more orders of magnitude inside the declared bound;
smooth integrands do far better than the bound requires.

Numerical differentiation (`R6.3`, `R8.6`) uses `mpmath.diff` and is given a
correspondingly looser tolerance (`1e-40` and `1e-30`), also declared per
check.

Truncated series are given enough terms that the tail is well below tolerance:
`F(z)` at `z = 1/1000` has ratio `0.256` per term and is summed to 140 terms;
`G_beta` and `G_arc` decay like `4^(-n)` at the tested arguments and are summed
to 140 and 180 terms respectively.

## Determinism

The result tables are byte-for-byte reproducible:

* No timestamps, hostnames, version strings, or random seeds appear in any
  output file.
* Row order follows the fixed `REGISTRY` tuple in `thales_ramanujan/checks.py`
  and the fixed `CLAUSEN_CANDIDATES` tuple; nothing is sorted at runtime and
  nothing depends on dictionary iteration order.
* Column order is fixed by explicit field lists in `scripts/run_all.py`.
* Errors are formatted through a single helper (`report.fmt_error`) to three
  significant digits, with exact zero rendered as `0e+0`.
* Files are written with `newline="\n"` regardless of platform.

Running `python scripts/run_all.py` twice and diffing `results/` should produce
no changes. This was verified during construction.

## Dependencies and why each is present

* **SymPy** — exact rational arithmetic and all symbolic reductions.
* **mpmath** — arbitrary-precision numerics. Ships with SymPy, so it is not a
  separate install burden.
* **pytest** — the test suite.

pandas is deliberately *not* a dependency. The result tables are a few dozen
rows of plain strings; the standard library `csv` module writes them
deterministically with less machinery and one fewer thing to pin.

## Layout note

The requested layout placed the mathematics directly in `scripts/`. This
repository instead keeps the mathematics in an importable `thales_ramanujan/`
package and makes the eight `scripts/verify_*.py` files thin wrappers around
it, with the requested filenames unchanged. This is what lets `pytest` exercise
exactly the same code paths that `run_all.py` runs, rather than a
reimplementation of them — a duplicated identity in a test file is a test of
the test, not of the result.

Each script inserts the repository root on `sys.path`, so the scripts run from
a bare checkout with no installation step. A root `conftest.py` does the same
for `pytest`.
