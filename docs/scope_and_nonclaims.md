# Scope and non-claims

This document draws the boundary precisely. It is the single most important
file in the repository: everything computational here exists to support the
division below, not to soften it.

The manuscript's own summary is that **the shared structure is a residue
lattice, not a derivation mechanism.** The code reproduces both halves of that
statement — the exact contact and the demonstrated obstruction.

---

## Established contact

Each item is exact, elementary, and reproduced by this repository.

| Item | Statement | Checks |
|---|---|---|
| Quarter-parameter residue lattice | The upper parameters of `3F2(1/4, 1/2, 3/4; 1, 1; 256 z)`, generating `A_n = (4n)!/(n!)^4`, are exactly the residues mod 1 that select the transcendental strata of `J(s)`. | `R1.1`, `R1.2`, `R2.2` |
| Gamma strata | `J(s)` is rational for `s` a non-negative integer, a rational multiple of `pi` for `s` in `Z + 1/2`, and a rational multiple of `Gamma(1/4)^2/sqrt(pi)` resp. `Gamma(3/4)^2/sqrt(pi)` for `s` in `Z + 1/4` resp. `Z + 3/4`. | `R2.1`–`R2.6` |
| Reflection products | `J(s) J(1-s) = pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s))`, giving `pi/20`, `pi sqrt(3)/35`, and `15 pi sqrt(3)/512` at the quarter, third, and sixth residues. | `R3.1`–`R3.5` |
| Lemniscatic period | `J(1/4) = (1/3) K(1/sqrt(2))`, with `K(1/sqrt(2)) = Gamma(1/4)^2/(4 sqrt(pi))`; and `J(3/4) = 3 pi / (20 K(1/sqrt(2)))`. | `R4.1`–`R4.3` |
| Arcsine and central-binomial moments | `(1/pi) int_0^1 x^n/sqrt(x(1-x)) dx = 4^(-n) C(2n,n)`, and `A_n = C(4n,2n) C(2n,n)^2 = 256^n m_(2n) m_n^2`. | `R5.1`–`R5.3` |

Two remarks on how far the contact reaches.

**Two CM fields, not one.** The reflection product lands in `Q pi` at the
self-complementary residue `1/4` (the lemniscatic field `Q(i)`, curve
`y^2 = x^3 - x`), and in `Q pi sqrt(3)` at the residues `1/3` and `1/6` (the
equianharmonic field `Q(sqrt(-3))`, curve `y^2 = x^3 - 1`). The construction is
therefore not confined to a single quadratic field. This *strengthens* the
statement of the obstruction rather than weakening it: reaching a second field
costs nothing and still supplies no quasi-period content.

**Removable poles.** The right-hand side of the product identity appears to
have poles at `s = -1/2, 1/2, 3/2`. All three are removable, because
`cot(pi s)` vanishes at exactly those points; check `R3.5` verifies the limits
agree with `J(s) J(1-s)`. The manuscript's qualifier "for all `s` avoiding the
poles of the factors" is therefore conservative, and the identity extends by
continuity.

---

## Demonstrated boundary

Each item is a computationally demonstrated negative, not an assertion.

### Failure of the Clausen condition

Clausen's identity requires the parameter condition `c = a + b + 1/2` for a
`2F1(a, b; c; x)` whose square is a Gauss-type `3F2`. Ramanujan's coefficient
satisfies it: `(a, b, c) = (1/8, 3/8, 1)` gives `a + b + 1/2 = 1 = c`, and the
square is the `3F2` of the main text. The tested Thales generating functions do
not:

* `G_beta(z) = sum_n J(n) z^n = 2F1(1, 1; 3/2; z/4)`. Here `a = b = 1` and
  `c = 3/2`, while the condition demands `c = 5/2`. It fails. (`R7.1`, `R7.2`)
* `G_arc(z) = sum_n 4^(-n) C(2n,n) z^n = (1-z)^(-1/2) = 1F0(1/2; ; z)`. This has
  no second numerator parameter and no denominator parameter at all, so the
  condition does not even apply; its square `(1-z)^(-1) = 1F0(1; ; z)` is
  degenerate and carries no modular content. (`R7.3`, `R7.4`)

The parameter table is regenerated as `results/clausen_parameter_audit.csv`,
including the positive control, so that the obstruction is visibly specific
rather than vacuous.

**Scope of this negative, stated exactly.** This is a demonstrated obstruction
*for the natural elementary Thales generating functions tested here*. It is not
a proof that no conceivable construction built from binary partition data could
ever produce a Clausen square. The repository makes the narrower claim and only
the narrower claim.

### Absence of modular multipliers and quasi-period data

The lemniscatic relation `J(1/4) = (1/3) K(1/sqrt(2))` supplies a **period**.
It does not supply a multiplier, a singular modulus, or a quasi-period
response. The distinction is load-bearing:

* The `pi` produced by the reflection product (`R3`) comes from Euler
  reflection of the gamma function. The `pi` produced on the Ramanujan side
  comes from the Legendre relation, pairing a period against a quasi-period at
  a complementary pair. The two mechanisms are structurally parallel and are
  not the same.
* The normalized tilt `A_s(alpha)/I_s(alpha) = alpha/(s+1)` (`R6.1`) is
  algebraic and acts on the *coordinate* `a`. The Euler operator
  `theta = z d/dz` acts on the *summation index* `n`. Conflating them would
  misidentify a coordinate first moment as a quasi-period response.
* The honest parameter response is
  `d_alpha log I_s(alpha) = psi(s+alpha+1) - psi(s-alpha+1)` (`R6.3`). It is
  transcendental, and it is digamma data. `pi` enters at first order via
  `psi(3/4) - psi(1/4) = pi`; Catalan's constant enters only at second order
  via `psi'(1/4) = pi^2 + 8 G` (`R8`). None of these is the modular derivative
  producing `1103 + 26390 n`.

### No derivation of Ramanujan's constants

Nothing in the elementary construction selects or produces `1103`, `26390`,
`396`, or `9801`. These are fixed by a singular modulus reached only through
the modular column. Check `R1.3` evaluates Ramanujan's series to confirm that
it does converge to `1/pi` — this is present purely to fix the target, and no
script in this repository derives it.

---

## Not claimed

The repository does **not** claim any of the following, and no output should be
read as supporting them:

1. **A new proof of Ramanujan's formula for `1/pi`.** The proof requires
   modular equations, singular moduli, multipliers, the Legendre relation, and
   derivative/quasi-period data. None of that is here.
2. **A derivation of modular equations from Thales geometry.** No modular
   equation appears anywhere in this repository.
3. **A universal relationship between binary partitions and modular forms.**
   The contact demonstrated is a shared residue lattice and a shared period,
   at specific residues.
4. **New transcendence or algebraic-independence results.** Every transcendence
   statement used here (`Gamma(1/4)`, `pi`, Catalan's constant) is classical
   and is used only as vocabulary.
5. **That the Clausen obstruction is universal.** See the scope statement
   above. It is demonstrated for the generating functions actually tested.

---

## One-line summary

The strongest defensible statement extractable from the comparison is the pair

    J(1/4) = (1/3) K(1/sqrt(2))          (an exact lemniscatic period relation)
    J(s) J(1-s) = pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s))

— the quarter stratum is one third of the lemniscatic period, and the collapse
to `pi/20` is a reflection phenomenon structurally parallel to, but not an
instance of, the Legendre relation that yields `pi` on the Ramanujan side.
