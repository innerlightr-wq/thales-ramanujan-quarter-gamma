# Mathematical scope

Definitions and statements of every identity reproduced by this repository,
with a pointer to the check that verifies it. Proofs are in the manuscript
(see [`zenodo_record.md`](zenodo_record.md)); this file fixes notation and
records what is being computed, so that the code can be read without the paper
open.

## Notation

Ramanujan side:

    A_n = (4n)! / (n!)^4
    F(z) = sum_{n>=0} A_n z^n

Thales side, attached to the binary partition `a + b = 1`:

    C(a) = a(1-a)                        (the partition product)
    h(a) = sqrt(a(1-a))                  (the altitude)
    J(s) = int_0^1 [a(1-a)]^s da         (the beta family)

Mirror-tilted family, for `|alpha| < s + 1`:

    I_s(alpha) = int_0^1 a^(s+alpha) (1-a)^(s-alpha) da
    A_s(alpha) = int_0^1 (2a-1) a^(s+alpha) (1-a)^(s-alpha) da

Arcsine moments:

    m_n = (1/pi) int_0^1 x^n / sqrt(x(1-x)) dx

Throughout, `(x)_n = Gamma(x+n)/Gamma(x)` is the rising factorial, `psi` is the
digamma function, `G` is Catalan's constant, and `K` is the complete elliptic
integral of the first kind in the modulus convention
`K(k) = int_0^(pi/2) dtheta / sqrt(1 - k^2 sin^2 theta)`.

> **Convention warning.** `mpmath.ellipk` takes the *parameter* `m = k^2`, not
> the modulus `k`. The modulus `k = 1/sqrt(2)` used in the manuscript therefore
> corresponds to `mp.ellipk(1/2)`. This is handled in
> `thales_ramanujan.core.K_lemniscatic` and is the single easiest place in this
> subject to introduce a silent factor error.

---

## 1. Ramanujan coefficient identity — `R1`

    A_n = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3                      (R1.1)

Follows from the Gauss multiplication formula at `m = 4`. Verified here by
induction — base case `n = 0` plus equality of the consecutive-term ratios,
both of which SymPy closes for symbolic `n` — and independently by exact
integer arithmetic for `n = 0..30`.

Consequently the generating function is the classical `3F2`

    F(z) = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z)                             (R1.2)

which is where the quarter-integer residue lattice `{1/4, 1/2, 3/4} (mod 1)`
enters. By Clausen's identity `F(z) = [2F1(1/8, 3/8; 1; 256 z)]^2`; this
square factorization is the entry point to the modular theory and is verified
here only as a positive control (`R7.5`).

Check `R1.3` evaluates Ramanujan's series itself to confirm convergence to
`1/pi`. It is an orientation anchor. Nothing in this repository derives it.

## 2. Thales beta family and gamma strata — `R2`

    J(s) = B(s+1, s+1) = Gamma(s+1)^2 / Gamma(2s+2),   Re(s) > -1      (R2.1)

Immediate from `[a(1-a)]^s = a^s (1-a)^s`. The recurrence
`Gamma(x+1) = x Gamma(x)` gives

    J(s+1) / J(s) = (s+1) / (2(2s+3)),                                  (R2.2)

rational in `s`, so the arithmetic nature of `J(s)` depends only on the residue
of `s` mod 1. This yields four strata:

    s in Z_{>=0}      =>  J(s) in Q
    s in Z + 1/2      =>  J(s) in Q pi
    s in Z + 1/4      =>  J(s) in Q Gamma(1/4)^2 / sqrt(pi)
    s in Z + 3/4      =>  J(s) in Q Gamma(3/4)^2 / sqrt(pi)

with the special values

    J(1/4) = Gamma(1/4)^2 / (12 sqrt(pi))                              (R2.3)
    J(3/4) = 3 Gamma(3/4)^2 / (10 sqrt(pi))                            (R2.4)
    J(1/2) = pi/8                                                      (R2.5)
    J(0) = 1,   J(1) = 1/6,   J(-1/2) = pi                             (R2.6)

The quarter strata are the locus of contact: the upper parameters of the `3F2`
in `R1.2` are exactly the residues that select them. This is a statement about
residue classes and the field `Q(pi, Gamma(1/4))`, not about a shared analytic
operation.

## 3. Conjugate-stratum product identity — `R3`

    J(s) J(1-s) = pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s))         (R3.1)

proved by Euler reflection applied twice, once to
`Gamma(s+1) Gamma(2-s)` and once to `Gamma(2s) Gamma(1-2s)`. Special cases:

    J(1/4) J(3/4) = pi/20                                              (R3.2)
    J(1/3) J(2/3) = pi sqrt(3)/35                                      (R3.3)
    J(1/6) J(5/6) = 15 pi sqrt(3)/512                                  (R3.4)

The quarter case collapses to a rational multiple of `pi` because `1/4` is the
self-complementary residue, `cot(pi/4) = 1`, and `Gamma(1/4) Gamma(3/4) = pi
sqrt(2)`. The third and sixth residues contribute `sqrt(3)` through
`cot(pi/3)` and `cot(pi/6)`, landing in the equianharmonic field
`Q(sqrt(-3))` rather than the lemniscatic `Q(i)`.

The apparent poles of the right-hand side at `s = -1/2, 1/2, 3/2` are removable
(`R3.5`), since `cot(pi s)` vanishes at exactly those points.

## 4. Lemniscatic-period relation — `R4`

With the classical evaluation

    K(1/sqrt(2)) = Gamma(1/4)^2 / (4 sqrt(pi))                         (R4.1)

comparison with `R2.3` gives

    J(1/4) = (1/3) K(1/sqrt(2))                                        (R4.2)

and, via `R3.2`,

    J(3/4) = 3 pi / (20 K(1/sqrt(2))).                                 (R4.3)

The quarter stratum is proportional to the lemniscatic period and the conjugate
stratum to its reciprocal; the product identity eliminates the period and
returns `pi/20`. This supplies a period only — no multiplier, no quasi-period —
and it selects the CM field `Q(i)`, not the singular modulus fixing the
constants of Ramanujan's series.

## 5. Arcsine moments — `R5`

The reciprocal altitude `1/h(x)` is the (unnormalized) arcsine density, with
total mass `B(1/2, 1/2) = pi = J(-1/2)`, and

    m_n = (1/pi) int_0^1 x^n / sqrt(x(1-x)) dx = 4^(-n) C(2n,n)        (R5.1)

Ramanujan's coefficient decomposes into central binomials,

    A_n = (4n)!/(n!)^4 = C(4n,2n) C(2n,n)^2                            (R5.2)

so each factor is an arcsine moment, equivalently

    A_n = 256^n m_(2n) m_n^2.                                          (R5.3)

`R5.3` is a restatement of `R5.1` and `R5.2` together, recorded because it
makes the "witness" language exact. It is emphatically **not** a generating
mechanism: a product of three moments is not itself a moment, reconstructing
`A_n` requires a threefold Mellin convolution against three distinct arcsine
measures, and the assembly into the `1/pi` series lies entirely outside the
elementary construction.

## 6. Mirror-tilt diagnostic — `R6`

    A_s(alpha) / I_s(alpha) = alpha / (s+1)                            (R6.1)

which is exactly `2 E[a] - 1` for `a ~ Beta(s+alpha+1, s-alpha+1)`, whose mean
is `(s+alpha+1)/(2s+2)` (`R6.2`). The ratio is algebraic in `(s, alpha)` even
when the normalizing mass `I_s(alpha)` lies in a transcendental stratum.

It is **not** an analogue of the Euler operator `theta = z d/dz`: the factor
`2a - 1` acts on the coordinate `a`, whereas `theta` acts on the index `n`.
The honest parameter response is

    d/d(alpha) log I_s(alpha) = psi(s+alpha+1) - psi(s-alpha+1)        (R6.3)

which is transcendental, and which belongs to beta/digamma analysis.

## 7. Clausen obstruction — `R7`

    G_beta(z) = sum_{n>=0} J(n) z^n = 2F1(1, 1; 3/2; z/4)              (R7.1)
    G_arc(z)  = sum_{n>=0} m_n z^n   = (1-z)^(-1/2) = 1F0(1/2; ; z)    (R7.3)

The first follows from `J(n) = (n!)^2/(2n+1)! = n! / (4^n (3/2)_n)`, using
`(2n+1)! = 4^n n! (3/2)_n`. Neither satisfies the Clausen parameter condition
`c = a + b + 1/2`: for `G_beta` the condition would require `c = 5/2` against
the actual `c = 3/2` (`R7.2`), and for `G_arc` there is no second numerator
parameter for the condition to apply to, its square `(1-z)^(-1)` being
degenerate (`R7.4`). The positive control `(1/8, 3/8, 1)` does satisfy it
(`R7.5`).

See [`scope_and_nonclaims.md`](scope_and_nonclaims.md) for the exact scope of
this negative result.

## 8. Parameter-response hierarchy — `R8`

    psi(3/4) - psi(1/4) = pi                                           (R8.1)
    psi(3/2) - psi(1) = 2 - 2 log 2                                    (R8.2)
    psi'(1/4) = pi^2 + 8 G                                             (R8.3)

with the supporting pair

    psi'(1/4) + psi'(3/4) = 2 pi^2                                     (R8.4)
    psi'(1/4) - psi'(3/4) = 16 G                                       (R8.5)

and the second-order response of the family itself,

    d^2/d(alpha)^2 log I_s(alpha) = psi'(s+alpha+1) + psi'(s-alpha+1)  (R8.6)

realised at `(s, alpha) = (-1/2, -1/4)`, where it equals
`psi'(1/4) + psi'(3/4)`.

The hierarchy is therefore: `pi` at first order via reflection of `psi`, the
logarithmic constant `2 - 2 log 2` on the diagonal, and Catalan's constant only
at second order. All of it is beta-, digamma-, and trigamma-level data, on the
elementary side of the boundary throughout. None of it is the modular
derivative responsible for the factor `1103 + 26390 n`.
