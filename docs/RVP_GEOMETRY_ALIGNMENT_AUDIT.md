# Geometry Beyond the Scalar Strain–Rotation Partition: An RVP Alignment and Ellipticity Audit

**Status.** Documentation of a completed, falsification-first computational audit.
This note is a research record, not a new exploratory round: every mathematical
claim below was checked in `scratch/rvp_ellipticity_audit.py` and reported in
`scratch/rvp_ellipticity_audit_report.txt`, both in this repository. This note
does not alter, and is not part of, the Ramanujan/Thales manuscript this
repository otherwise documents (see the main [README](../README.md)), and it
does not modify the RVP manuscript itself. RVP — "Strain–Vorticity Interaction
and Rotational Coherence" — is a separate technical note by the same author;
its PDF is included in this repository for the audit's own use, not as a
manuscript this repository publishes.

**Central question.** Does ellipticity, strain-eigenvalue anisotropy, or the
Thales–Legendre complementary-modulus geometry developed elsewhere in this
repository provide independent, non-redundant information about
strain–vorticity interaction beyond RVP's scalar magnitude coordinate ζ?

**Answer, in one line.** Strain-eigenvalue shape and vorticity alignment do
carry independent information beyond ζ (proved exactly below); the
Thales–Legendre elliptic-integral apparatus specifically does not — the one
natural, non-circular mapping from RVP's own data onto a Thales `(k,k')` pair
is an exact but information-free reparametrization of ζ.

---

## 1. Purpose and scope

RVP.pdf decomposes the velocity-gradient tensor as `∇u = S + Ω` and defines

```
ζ = (‖Ω‖_F² − ‖S‖_F²) / (‖Ω‖_F² + ‖S‖_F²) = 2Q/D
```

an exact affine rescaling of the second invariant `Q` (Hunt–Wray–Moin
convention), and, since its September 2026 interaction-response revision,
replaces any strain-versus-rotation magnitude criterion with the exact
interaction term appearing in the enstrophy equation,

```
P = ω · S ω = |ω|² σ_eff,     σ_eff = n̂ · S n̂,     n̂ = ω/|ω|.
```

RVP.pdf itself states plainly that it no longer proposes a universal
strain-rotation corridor or threshold, and that ζ is retained only "as a point
of comparison," not as a criterion. This audit asks a narrower, falsifiable
question left open by that revision: does the elliptic/modular geometry
developed in this repository's Thales–Legendre material (`k`, `k'`,
`k²+k'²=1`, complete elliptic integrals `K(k)`, `K(k')`, modular self-duality
at `k=k'=1/√2`) supply anything ζ is missing?

---

## 2. Central result: ζ does not determine P

**Claim (exact mathematics).** Fixing `‖S‖_F`, `‖Ω‖_F` (hence ζ) does not fix
`P`; `P` can be driven from a strictly positive value to a strictly negative
value by changing only the orientation of ω, with `S` and `|ω|` held
completely unchanged.

**Exact counterexample** (`scratch/rvp_ellipticity_audit.py`, "Counterexample
B"). Fix

```
S = [[2, 0, 0],
     [0,-1, 0],
     [0, 0,-1]]        (traceless, ‖S‖_F² = 6)
```

and fix `|ω|² = 6` (so `‖Ω‖_F² = |ω|²/2 = 3`, using the standard identity
`Ω_ij = −½ε_ijk ω_k ⇒ ‖Ω‖_F² = |ω|²/2`, re-verified symbolically in the audit
script). This fixes `D = ‖S‖_F² + ‖Ω‖_F² = 9` and

```
ζ = (‖Ω‖_F² − ‖S‖_F²)/D = −1/3
```

**identically**, for every choice of ω's direction. Varying only that
direction:

| ω direction | P = ω·Sω |
|---|---|
| aligned with the distinct (extensional) eigenvector, e₁ | **+12** |
| aligned with the degenerate-plane eigenvector, e₂ | **−6** |
| 45° mix of e₁, e₂ | +3 |
| (1,1,1)/√3 | 0 |

`S` is never altered — its eigenvalue shape is fixed throughout — and `|ω|`
never changes. Only alignment changes, and `P` swings from +12 to −6 while
`ζ = −1/3` stays fixed. **This one example already refutes any claim that ζ,
or any function of `‖S‖_F, ‖Ω‖_F` alone, determines the sign or magnitude of
`P`.**

A second, independent counterexample (`scratch/rvp_ellipticity_audit.py`,
"Counterexample A") shows the same conclusion via a *different* mechanism —
holding `ω` fixed and changing `S`'s eigenvalue shape at fixed `‖S‖_F`:
`S₁ = diag(2,−1,−1)` and `S₂ = diag(√3,0,−√3)` (rescaled to the same
`‖S‖_F = √6`) give, for `ω=(1,0,0)`, `P₁ = 2` and `P₂ = √3`, exactly.

**Therefore: magnitude balance ≠ interaction geometry.** ζ is a function of
`‖S‖_F, ‖Ω‖_F` only; `P` additionally depends on `S`'s eigenvalue shape and on
the orientation of ω relative to `S`'s eigenframe, neither of which ζ
records. Classification: **exact mathematics**, not an approximation, not a
numerical tendency.

---

## 3. Eigenframe decomposition

This is standard material, restated here (following RVP.pdf §4.1, itself
elementary) to fix notation for what follows. Writing `S eᵢ = λᵢ eᵢ` for `S`'s
real orthonormal eigenbasis and `ω = Σᵢ ωᵢ eᵢ`,

```
P = ω · S ω = Σᵢ λᵢ ωᵢ².
```

For incompressible flow, `tr(S) = ∇·u = 0`, i.e. `λ₁ + λ₂ + λ₃ = 0`. `P`
therefore depends on four logically separate ingredients:

1. **strain magnitude** — `‖S‖_F = √(Σᵢλᵢ²)`;
2. **strain eigenvalue shape** — the one residual degree of freedom once
   `‖S‖_F` is fixed and `λ₁+λ₂+λ₃=0` is imposed (a single angle on a circle;
   the classical *Lode angle* of continuum mechanics, equivalently the
   Lumley/Lumley-triangle invariants of turbulence anisotropy);
3. **vorticity magnitude** — `|ω|`;
4. **vorticity orientation** — the direction cosines of ω in `S`'s eigenframe
   (two independent angles once `|ω|` is fixed).

ζ is built entirely from (1) and (3) (via `‖S‖_F, ‖Ω‖_F = |ω|/√2`); it
contains no information from (2) or (4). That vorticity's alignment with the
strain eigenframe — not merely strain magnitude — controls vortex stretching
is long-established in the turbulence literature, independent of this audit:
Betchov (1956) derived an exact isotropic-turbulence relation between mean
strain self-amplification and mean enstrophy production that already presumes
strain and vorticity are dynamically coupled through alignment, not through a
magnitude comparison; and Tsinober, Shtilman & Vaisburd (1997), combining
direct numerical simulation and laboratory measurement, report that enstrophy
production is governed by the alignment of vorticity with the strain-rate
eigenvectors — in particular a preferential alignment with the intermediate
eigenvector `λ₂` — and explicitly that regions of weak vorticity are "not
structureless," i.e. magnitude alone does not organize the production field.
(Both sources are checked against primary bibliographic metadata in this
audit; see References.)

---

## 4. The Thales mapping — exact but redundant

The most natural, non-circular way to apply this repository's Thales–Legendre
`(k,k')` construction to RVP's own data is to normalize the two Frobenius
energies that already sum to `D`:

```
a := ‖Ω‖_F² / D,     b := ‖S‖_F² / D,     a + b = 1   (exact, by definition of D)
k := 2√(ab),          k' := |a − b|.
```

`a, b` are not arbitrarily imposed: they are already meaningful fractions of
the local velocity-gradient energy, and `a+b=1` follows automatically from the
definition of `D`, not from a chosen rescaling. This is therefore a fair test
of the construction, not a straw one.

**Exact result** (verified symbolically, SymPy, in
`scratch/rvp_ellipticity_audit.py`):

```
k'  =  |ζ|
k   =  √(1 − ζ²)
```

**identically.** `k² + k'² − 1`, `k' − |ζ|`, and `k² − (1−ζ²)` all simplify to
exactly `0`. The Thales complementary-modulus pair `(k,k')`, applied this way,
is a smooth, bijective reparametrization of the single existing scalar ζ. It
therefore supplies **no additional state information beyond ζ** — and, per
§2, ζ itself already does not determine `P`.

**Classification: EXACT REPARAMETRIZATION / REDUNDANT FOR RVP DYNAMICS.**

A second, shape-based Thales-style mapping was also tested — `a :=
λ_max/(λ_max−λ_min)`, `b := −λ_min/(λ_max−λ_min)`, using only `S`'s two
extreme eigenvalues (also satisfies `a+b=1` automatically, for any nonzero
traceless `S`). This one is *not* redundant with ζ — it carries genuine shape
information — but it is a coarsened restatement of the same Lode angle already
named in §3, and it involves no complete elliptic integral, complementary
period, or modular structure of any kind; it is elementary eigenvalue algebra.
Neither Thales-style mapping resolves the alignment gap of §2.

The self-dual point of this construction, `k=k'=1/√2`, corresponds under the
magnitude-based mapping to `ζ = ±1/√2 ≈ ±0.70711`. No governing equation, no
canonical flow (solid rotation, pure strain, simple shear, Lamb–Oseen,
Burgers), and no sharp bound tested in this audit selects this value
dynamically; any one of these flows can obviously be *tuned*, via its own free
parameter, to pass through `ζ = ±1/√2` at some particular parameter value —
exactly as it passes through every other value of ζ in its continuous range.
Passing through a value under parameter tuning is not the same as that value
being dynamically selected. **No physical self-dual point was found.**

---

## 5. Two centers: magnitude balance and interaction neutrality are different centers

Two natural "zero" conditions exist in RVP's state space, and they are not the
same condition:

- **`ζ = 0`** — equal Frobenius-norm contribution from strain and rotation
  (`‖S‖_F = ‖Ω‖_F`); a statement purely about magnitudes.
- **`P = ω·Sω = 0`** — zero instantaneous strain–vorticity
  production/compression at that point; a statement about alignment and shape
  as much as magnitude (`P=0` whenever ω lies in the null cone of the
  quadratic form `S`, regardless of `‖S‖_F, ‖Ω‖_F`).

Neither condition generally implies the other:

- **Simple shear**: `ζ = 0` *and* `P = 0` — but this is a coincidence of this
  one flow (RVP.pdf §5.3 states this explicitly in its own text), not a
  general identity: `S`'s single off-diagonal shear component happens to be
  exactly orthogonal, in the quadratic-form sense, to the shear flow's own
  vorticity direction.
- **Burgers vortex**: `P = a ω_z² > 0` for *every* `a > 0` (RVP.pdf Eq. 8,
  exact, independent of the swirl profile's radial shape) — `P` never reaches
  zero for a sustained vortex — while `ζ` is driven monotonically toward `−1`
  as the sustaining strain rate `a` grows. The two loci are approached from
  opposite ends of the parameter range and never coincide for `a>0`.
- **General state-space argument**: `{ζ=0}` is a codimension-1 hypersurface
  depending only on the two Frobenius norms; `{P=0}` is a codimension-1
  hypersurface depending on shape and alignment as well. Generic surfaces
  built from only partially overlapping data intersect transversally; they do
  not coincide as sets.

This is an **RVP state-space distinction**, derivable directly from RVP's own
definitions, independent of any external framework. It is analogous *in
spirit* — not in mathematical structure — to the separate repository note
distinguishing an equipartition center from a modular self-dual center for an
abstract Thales–Legendre partition (`a=b` vs. `k=k'`): both cases warn against
conflating two differently-defined "symmetric" points. The RVP two centers are
codimension-1 surfaces in an ordinary tensor state space; the Thales two
centers are two points on a single one-dimensional modulus line inside genuine
modular/elliptic-period geometry. **The two structures are not being
identified here — only their cautionary shape is being noted as similar.**

---

## 6. Sharp-bound geometry: axisymmetry, not ellipticity

RVP.pdf (Theorem 4.1, a stated specialization of the Wolkowicz–Styan (1980)
eigenvalue-trace bound) proves, for real traceless symmetric `S`,

```
|ω · S ω|  ≤  √(2/3) ‖S‖_F |ω|²,
```

with equality **if and only if `S` has a repeated eigenvalue** — eigenvalues
proportional to `(2,−1,−1)` up to sign and permutation — and ω aligned with
the non-repeated eigenvector. This audit re-derives the same equality
condition independently (Lagrange multipliers on the constraint set
`λ₁+λ₂+λ₃=0`, `Σλᵢ²=‖S‖_F²`; exact symbolic ratio at the extremum = `√6/3 =
√(2/3)`).

The equality geometry is a discrete **eigenvalue-degeneracy** condition (two
eigenvalues coincide), not a continuously variable eccentricity. Its standard
fluid-mechanics name is **tensorial axisymmetry** — the strain-rate tensor has
an axis of rotational symmetry — exactly the two straight edges of the Lumley
turbulence-anisotropy triangle. A geometric ellipsoid of revolution (a
spheroid) can legitimately be pictured at this state, but nothing about
complete elliptic integrals, a complementary modulus, or modular geometry
enters the derivation anywhere; it is elementary constrained optimization on
three real numbers.

**This bound is not reinterpreted through the Thales modulus.** `√(2/3) ≈
0.8165` and the self-dual value `1/√2 ≈ 0.7071` (§4) are unrelated numbers
answering unrelated questions; no identity connects them, and none is claimed.

---

## 7. Elliptical vortices: what the classical literature actually establishes

Literature verified against primary bibliographic sources during this
documentation pass (DOI metadata confirmed via CrossRef; abstract content
confirmed where noted; see References for exact status of each item).

**What is established.** Kida (1981a), "Motion of an Elliptic Vortex in a
Uniform Shear Flow" (*J. Phys. Soc. Japan* **50**(10), 3517–3520), gives an
exact solution for a uniform-vorticity elliptical vortex patch immersed in a
uniform external strain/shear. Its own abstract (verified via the journal's
J-STAGE listing) states: "The elliptic shape is preserved and the area of the
vortex is conserved but the axis ratio of the ellipse changes in general,"
with the vortex exhibiting rotation and nutation depending on flow
parameters, becoming elongated under strong strain, and admitting stationary
elliptic solutions under weak strain. **This is a legitimate, primary-source
example in which vortex-core aspect ratio and orientation relative to an
external strain field are dynamically meaningful** — real content, not
analogy, and it directly supports the general RVP lesson that shape and
orientation carry information a scalar magnitude ratio cannot (§2–3).

**What is not established by this same paper**: the fetched abstract makes no
mention of complete elliptic integrals, a complementary modulus, or modular
structure in its solution. Absence from an abstract summary is not proof of
absence in the paper's internal equations (a nutation period could still be
expressed via `K(k)` in the body even where the abstract does not say so);
this was **not resolved** by full-text access in this session and is recorded
here as an open literature question, not asserted either way.

**A genuinely confirmed elliptic-integral occurrence — in a different
problem.** Kida (1981b), "A vortex filament moving without change of form"
(*J. Fluid Mech.* **112**, 397–409), solves for thin vortex-filament shapes
under the localized-induction approximation that translate/rotate without
changing form; per verified bibliographic summary, these solutions —
including helical filaments and forms related to travelling-wave solutions of
the nonlinear Schrödinger equation — are "expressed through elliptic
integrals." **This is a real, primary-source-adjacent confirmation that
complete elliptic integrals arise naturally in *some* classical vortex
problem.** It is not the same problem as Kida (1981a)'s elliptical patch, and
it is not RVP's pointwise tensor quantity `P = ω·Sω`. No mapping from this
filament problem's modulus to RVP's state, or to this repository's
Thales–Legendre modulus, was constructed or is claimed.

**Three concepts kept explicit and separate, per the audit's own
terminology discipline:**

1. **elliptical vortex geometry** — Kida (1981a): an actual elliptical patch
   boundary, dynamically meaningful, confirmed.
2. **elliptic-integral mathematical solution structure** — Kida (1981b): `K`,
   `E`, `Π`-type integrals arising in a filament-shape solution, confirmed in
   a *different* problem.
3. **Thales complementary-modulus geometry** — this repository's `k, k',
   K(k'), τ`, modular self-duality: confirmed nowhere in the fluid literature
   checked here, and not established as governing either (1) or (2).

**Elliptic instability** (Kerswell, *Annu. Rev. Fluid Mech.* **34**, 83–113,
2002, metadata verified; content partially verified — the review's own
framing describes "the linear instability mechanism that tends to break up
regions of elliptical streamlines in a rotating flow" across three
independently-discovered contexts) is a short-wavelength parametric
instability of elliptical *streamlines*, tracing to Moore & Saffman (1971) and
studied by Bayly (1986), Pierrehumbert (1986), and Waleffe (1990). Its
classical weak-strain growth-rate result is, by standard characterization,
algebraic in the strain rate/eccentricity rather than built from `K(k)` — this
specific mathematical-form claim was **not independently re-verified against
primary text in this session** (the relevant pages were paywalled in every
route attempted); it is reported here as the standard textbook
characterization, explicitly flagged as not re-confirmed this session, per
the instruction to say so plainly when a claim is not verified.

Dritschel (1990), "The stability of elliptical vortices in an external
straining flow" (*J. Fluid Mech.* **210**, 223–261, metadata verified), was
identified as a further primary reference for elliptical-patch stability in
strain; its content was not independently fetched in this session and no
claim is attributed to it beyond its existence and subject matter.

Kirchhoff's classical steady-rotation result for a uniform-vorticity ellipse
(`Ω = ω₀ ab/(a+b)²`, algebraic, no elliptic integral needed) is widely
reported in standard references (e.g. Lamb, *Hydrodynamics*) but was **not
independently re-verified against a primary or tertiary source in this
session**; it is mentioned here for context only, explicitly flagged as
unverified this session, and is not used to support any claim in this note.

**Conclusion of §7, stated as instructed:** the existence of elliptic
integrals in *some* classical vortex problem (Kida 1981b) does **not**
establish that the Thales–Legendre complementary modulus of this repository
controls RVP, or controls any of the elliptical-patch/elliptic-instability
literature checked here. These remain separate claims.

---

## 8. Information hierarchy

```
scalar magnitude balance                magnitude + strain shape          magnitude + strain shape
        ζ                        <         (Lode angle)             <     + vorticity alignment      →   P = ω · S ω
                                                                       (local interaction)
```

This is a **descriptive ordering of what information §2–3 show is needed**,
not a universal minimal-coordinate theorem. In the general (non-axisymmetric)
case, the full completion beyond ζ requires the shape angle (1 scalar) *and*
two independent alignment angles (the orientation of ω in `S`'s eigenframe) —
three additional scalars, not one. Only in the special axisymmetric case
(§6) does the shape freedom collapse, leaving a single alignment angle. For
the narrower question of whether `P` can change sign at fixed ζ, alignment
alone (two angles), with shape held fixed, already suffices (§2,
Counterexample B) — ellipticity/shape is not logically necessary to produce
that degeneracy, only to set its amplitude.

---

## 9. Relation to the OpenAI Navier–Stokes construction

OpenAI's *Finite Time Blowup for Navier–Stokes* preprint (cited identically to
RVP.pdf's own reference [15]; `navier-stokes.pdf` in this repository)
constructs, for every positive viscosity, a forced solution developing
unbounded velocity in finite time from rest with uniformly bounded kinetic
energy. Its own text (§2.1, §3.1) gives radial and axial core length scales

```
ℓ_r ≍ τ^{1/2},    ℓ_z ≍ τ^{1/2−h}    (0 < h < 1/100 fixed, τ = 1−t),
```

so the aspect ratio `q(t) := ℓ_r(t)/ℓ_z(t) ≍ τ^h → 0` as the singular time is
approached — the paper's own words describe the core becoming "an
increasingly slender column."

**Supported lesson only:** this is a genuine, well-documented example of
anisotropic geometric collapse in a rigorous Navier–Stokes construction,
reinforcing — structurally, not quantitatively — the general point of this
note that directional/shape structure is a separate axis of information from
scalar strain/rotation magnitude. RVP.pdf's own text (§8.1) already notes the
construction's background mechanism is structurally the same
stretching mechanism examined for the Burgers vortex, and its own §8.2/§9.5
leave the sign of `σ_eff` for this construction **explicitly open** — this
audit did not attempt to resolve that, and does not change that status.

**Explicitly not claimed, per this audit's own findings and instructions:**

- the OpenAI construction does **not** validate RVP;
- its anisotropic exponent `h` (or `q(t)`) is **not** a Thales modulus;
- its aspect ratio does **not** obey, or relate to, the RVP `(a,b)` magnitude
  partition of §4;
- it does **not** select, or relate to, the Thales self-dual point
  `k=k'=1/√2`.

`log q(t) = h · log τ` is a similarity-scaling exponent evaluated as a system
approaches a singularity in *time* — a different kind of quantity from either
a strain-eigenvalue anisotropy ratio or the Thales–Legendre modulus `k`,
evaluated at fixed time. No exponent-matching argument is advanced here.

---

## 10. Claim table

| Claim | Status | Evidence |
|---|---|---|
| ζ does not determine P | **EXACT** | §2, Counterexamples A & B (symbolic, `scratch/rvp_ellipticity_audit.py`) |
| Alignment changes P at fixed ζ (S and \|ω\| unchanged) | **EXACT** | §2, Counterexample B: P = +12 → −6 at fixed ζ = −1/3 |
| Strain eigengeometry contains information absent from ζ | **EXACT / STANDARD FLUID MECHANICS** | §3; Betchov (1956); Tsinober, Shtilman & Vaisburd (1997) |
| Natural (magnitude-based) Thales mapping adds new information | **FALSE / REDUNDANT** | §4: k' = \|ζ\| and k = √(1−ζ²) exactly (symbolic) |
| Thales self-dual point is physically selected in RVP | **NO EVIDENCE** | §4: no governing equation, bound, or canonical flow selects ζ = 1/√2 |
| RVP sharp-bound equality is axisymmetric, not elliptic-integral geometry | **EXACT** | §6: equality ⟺ repeated eigenvalue of S; re-derived symbolically |
| Elliptical vortex aspect ratio can be dynamically meaningful | **KNOWN FLUID MECHANICS (verified)** | §7: Kida (1981a), abstract confirmed via J-STAGE |
| Elliptic integrals occur in some vortex problems | **KNOWN MATHEMATICS/FLUID MECHANICS (verified, different problem)** | §7: Kida (1981b), vortex-filament (LIA) shape solutions |
| Thales complementary periods govern vortex stretching | **NOT ESTABLISHED** | §4, §7: no derivation found or attempted; explicitly not claimed |
| Elliptic instability growth rate is algebraic, not K(k)-based | **STANDARD CHARACTERIZATION, NOT RE-VERIFIED THIS SESSION** | §7: Kerswell (2002) metadata/partial content verified; growth-rate form not re-confirmed (paywalled) |

---

## 11. What this audit does not establish

- **No new universal vortex threshold.** RVP.pdf already withdrew its earlier
  strain-rotation corridor claim; nothing here reinstates any threshold, and
  no new one is proposed.
- **No physical significance for `k = 1/√2`.** §4 finds no fluid-mechanical
  quantity, bound, or flow that dynamically selects the Thales self-dual
  modulus.
- **No evidence that `K(k)`/`K(k')` controls enstrophy production.** No
  derivation connecting complete elliptic integrals to `P = ω·Sω` was found or
  constructed anywhere in this audit.
- **No identification of Thales modular geometry with vortex dynamics.** The
  one confirmed elliptic-integral occurrence in classical vortex theory (Kida
  1981b, vortex filaments) is a different problem from both RVP's tensor
  interaction and this repository's Thales–Legendre partition; no bridge
  between them is asserted.
- **No proof that one scalar ellipticity coordinate completes ζ.** §8: the
  general completion needs shape *and* two alignment angles; only the special
  axisymmetric case reduces to one further scalar.
- **No claim that the OpenAI Navier–Stokes construction validates RVP.** §9:
  used only as a documented anisotropic-collapse example; its sign of
  `σ_eff` remains open, exactly as RVP.pdf itself already states.

---

## 12. References

Full computational detail, exact symbolic derivations, and the complete
falsification-first audit trail (including the literature-access limitations
of the earlier pass) are recorded in `scratch/rvp_ellipticity_audit.py` and
`scratch/rvp_ellipticity_audit_report.txt`. Below, sources already cited by
RVP.pdf are marked *(RVP)*; sources checked and added specifically for this
geometry/alignment audit are marked *(new)*, with their verification status
in this session stated explicitly.

**Sources already cited in RVP.pdf, reused here *(RVP)*:**

1. J. C. R. Hunt, A. A. Wray and P. Moin, "Eddies, streams, and convergence
   zones in turbulent flows," in *Studying Turbulence Using Numerical
   Simulation Databases, 2* — Proc. 1988 Summer Program, Center for
   Turbulence Research, 193–208 (1988).
2. C. Liu, Y. Wang, Y. Yang and Z. Duan, "New omega vortex identification
   method," *Sci. China Phys. Mech. Astron.* **59**, 684711 (2016).
   DOI:10.1007/s11433-016-0022-6.
3. A. Okubo, "Horizontal dispersion of floatable particles in the vicinity of
   velocity singularities such as convergences," *Deep-Sea Res. Oceanogr.
   Abstr.* **17**(3), 445–454 (1970). DOI:10.1016/0011-7471(70)90059-8.
4. J. Weiss, "The dynamics of enstrophy transfer in two-dimensional
   hydrodynamics," *Physica D* **48**(2–3), 273–294 (1991).
   DOI:10.1016/0167-2789(91)90088-Q.
5. R. Betchov, "An inequality concerning the production of vorticity in
   isotropic turbulence," *J. Fluid Mech.* **1**(5), 497–504 (1956).
   DOI:10.1017/S0022112056000317. — metadata re-verified via CrossRef this
   session; used here (§3) only for the claim it actually supports (an exact
   isotropic-turbulence relation coupling mean strain self-amplification to
   mean enstrophy production), not re-derived or extended.
6. W. T. Ashurst, A. R. Kerstein, R. M. Kerr and C. H. Gibson, "Alignment of
   vorticity and scalar gradient with strain rate in simulated Navier–Stokes
   turbulence," *Phys. Fluids* **30**, 2343–2353 (1987).
7. H. Wolkowicz and G. P. H. Styan, "Bounds for eigenvalues using traces,"
   *Linear Algebra Appl.* **29**, 471–506 (1980).
8. J. T. Beale, T. Kato and A. Majda, "Remarks on the breakdown of smooth
   solutions for the 3-D Euler equations," *Commun. Math. Phys.* **94**(1),
   61–66 (1984).
9. J. M. Burgers, "A mathematical model illustrating the theory of
   turbulence," *Adv. Appl. Mech.* **1**, 171–199 (1948).
10. N. Rott, "On the viscous core of a line vortex," *Z. Angew. Math. Phys.*
    **9**(5), 543–553 (1958).
15. OpenAI, *Finite Time Blowup for Navier–Stokes*, preprint (September
    2026). https://openai.com/research.

**Sources checked and added specifically for this geometry/alignment audit
*(new)*, with verification status:**

- S. Kida, "Motion of an Elliptic Vortex in a Uniform Shear Flow," *J. Phys.
  Soc. Japan* **50**(10), 3517–3520 (1981). DOI:10.1143/JPSJ.50.3517. —
  **Metadata and abstract content verified this session** (title/journal/
  volume/pages/year confirmed via CrossRef; abstract content confirmed via
  J-STAGE listing: elliptical shape preserved, area conserved, axis ratio
  varies, rotation/nutation, strain-dependent stability). Used in §7 for
  exactly this claim.
- S. Kida, "A vortex filament moving without change of form," *J. Fluid
  Mech.* **112**, 397–409 (1981). DOI:10.1017/S0022112081000475. —
  **Metadata verified via CrossRef; content (elliptic-integral solution
  structure for filament shapes) verified via a secondary bibliographic
  summary, not the full primary text** (Cambridge Core returned HTTP 403 in
  every route attempted this session). Used in §7 only for the claim that
  elliptic integrals arise in this (distinct) vortex-filament problem; the
  "first, second, and third kind" breakdown specifically was not
  independently re-confirmed against primary text this session.
- A. Tsinober, L. Shtilman and H. Vaisburd, "A study of properties of vortex
  stretching and enstrophy generation in numerical and laboratory
  turbulence," *Fluid Dynamics Research* **21**(6), 477–494 (1997).
  DOI:10.1016/S0169-5983(97)00022-1. — **Metadata and content verified this
  session** (title/journal/volume/page/year confirmed; content confirmed via
  IOPscience listing: enstrophy production depends on alignment of vorticity
  with strain-rate eigenvectors, notably λ₂; weak-vorticity regions are not
  structureless). Used in §3 for exactly this claim.
- R. R. Kerswell, "Elliptical Instability," *Annu. Rev. Fluid Mech.* **34**,
  83–113 (2002). DOI:10.1146/annurev.fluid.34.081701.171829. — **Metadata and
  partial content verified this session** (title/journal/volume/pages/year
  confirmed via CrossRef; the review's characterization of the instability
  mechanism as breaking up elliptical-streamline regions, discovered
  independently in three contexts, confirmed). **Not verified this session:**
  the specific mathematical form (algebraic vs. elliptic-integral) of the
  classical growth-rate formula — the relevant pages were paywalled in every
  route attempted. Reported in §7 as the standard textbook characterization
  only, explicitly flagged as unconfirmed here.
- D. W. Moore and P. G. Saffman, "Structure of a Line Vortex in an Imposed
  Strain," in *Aircraft Wake Turbulence and Its Detection*, Springer US,
  339–354 (1971). DOI:10.1007/978-1-4684-8346-8_20. — **Metadata verified via
  CrossRef only**; content not independently fetched this session. Cited in
  §7 only for its existence and subject matter (the origin of the
  strain-resonance mechanism later named elliptic instability), not for any
  specific formula.
- D. G. Dritschel, "The stability of elliptical vortices in an external
  straining flow," *J. Fluid Mech.* **210**, 223–261 (1990).
  DOI:10.1017/S0022112090001276. — **Metadata verified via CrossRef only**;
  content not independently fetched this session. Cited in §7 only for its
  existence and subject matter.
- Kirchhoff's classical steady-rotation result for a uniform-vorticity
  elliptical patch (`Ω = ω₀ ab/(a+b)²`) — **not independently verified
  against any primary or tertiary source in this session**; mentioned in §7
  for context only, and no claim in this note depends on it.

**Explicitly not cited:** no search-result snippet, and no AI system (this
model or any other), is cited as a source anywhere in this note, per the
audit's own instructions.
