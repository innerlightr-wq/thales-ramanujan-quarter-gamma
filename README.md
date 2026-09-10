# thales-ramanujan-quarter-gamma

Reproducibility and verification companion for the technical note *Ramanujan's
Formula for 1/π and the Thales Integral Witness*.

The note identifies a **limited but exact** overlap between the
quarter-parameter hypergeometric lattice underlying Ramanujan's fastest
classical series for 1/π and the gamma strata of the Thales partition beta
family

```
J(s) = ∫₀¹ [a(1-a)]ˢ da = B(s+1, s+1) = Γ(s+1)² / Γ(2s+2).
```

Both objects are organized by the same quarter-integer residue lattice
`{1/4, 1/2, 3/4} (mod 1)`. This repository verifies every identity that makes
that contact exact, **and** verifies the demonstrated obstruction at which the
analogy stops. The contribution is organizational and diagnostic: the shared
lattice is not a mechanism.

---

## Manuscript

> De Jesús, Elias. (2026). *Ramanujan's Formula for 1/π and the Thales Integral
> Witness*. Zenodo. https://doi.org/10.5281/zenodo.20665588

**Read the paper here → [https://doi.org/10.5281/zenodo.20665588](https://doi.org/10.5281/zenodo.20665588)**
(concept DOI; always resolves to the latest version).

> **Version reproduced.** The identities verified in this repository correspond
> to the deposited version with DOI
> [10.5281/zenodo.20672985](https://doi.org/10.5281/zenodo.20672985). Cite the
> concept DOI above for the work; cite the version DOI when you need to state
> exactly which deposit a computation corresponds to.

**This repository does not contain the manuscript PDF or LaTeX source. Zenodo
is the authoritative archival location for the paper.**

---

## Scope of this repository

GitHub serves only as the reproducibility, verification, documentation, and
code companion. Concretely, it:

* reproduces the exact mathematical identities of the note, symbolically where
  practical and numerically at high precision as an independent audit;
* records the demonstrated limits of the analogy, including a positive control
  so that the negative result is visibly specific rather than vacuous;
* regenerates three deterministic result tables under `results/`.

It does not restate the manuscript's prose, and it introduces no mathematical
claim beyond it.

---

## Exact results reproduced

34 checks across eight groups. Every one currently passes.

| Group | Content | Script |
|---|---|---|
| `R1` | `Aₙ = (4n)!/(n!)⁴ = 256ⁿ (1/4)ₙ(1/2)ₙ(3/4)ₙ / (1)ₙ³`; generating function `₃F₂(1/4, 1/2, 3/4; 1, 1; 256z)` | `verify_ramanujan_coefficients.py` |
| `R2` | `J(s) = B(s+1,s+1) = Γ(s+1)²/Γ(2s+2)`; the four gamma strata; `J(1/4) = Γ(1/4)²/(12√π)`, `J(3/4) = 3Γ(3/4)²/(10√π)`, `J(1/2) = π/8` | `verify_beta_gamma_values.py` |
| `R3` | `J(s)J(1-s) = π s(1-s) cot(πs) / (2(2s+1)(1-2s)(3-2s))`; the evaluations `π/20`, `π√3/35`, `15π√3/512` | `verify_reflection_products.py` |
| `R4` | `K(1/√2) = Γ(1/4)²/(4√π)`; `J(1/4) = ⅓ K(1/√2)`; `J(3/4) = 3π/(20 K(1/√2))` | `verify_lemniscatic_relation.py` |
| `R5` | arcsine moments `4⁻ⁿ C(2n,n)`; the decomposition `Aₙ = C(4n,2n) C(2n,n)²` | `verify_arcsine_moments.py` |
| `R6` | mirror tilt `A_s(α)/I_s(α) = α/(s+1)`; the true response `∂_α log I_s(α) = ψ(s+α+1) − ψ(s−α+1)` | `verify_mirror_tilt.py` |
| `R7` | `G_beta(z) = ₂F₁(1,1;3/2;z/4)` and `G_arc(z) = (1−z)^(−1/2)`, neither satisfying `c = a + b + ½` | `verify_clausen_obstruction.py` |
| `R8` | `ψ(3/4) − ψ(1/4) = π`; `ψ(3/2) − ψ(1) = 2 − 2 log 2`; `ψ'(1/4) = π² + 8G` | `verify_response_hierarchy.py` |

Full statements, notation, and a manuscript-to-check mapping:
[`docs/mathematical_scope.md`](docs/mathematical_scope.md) and
[`docs/zenodo_record.md`](docs/zenodo_record.md).

---

## Demonstrated limits of the analogy

The exact contact is real, and it is narrow. Stated plainly:

* **The shared quarter-gamma lattice is exact.** The upper parameters of
  Ramanujan's `₃F₂` are precisely the residues that select the transcendental
  strata of `J(s)`.
* **The lemniscatic-period relation is exact.** `J(1/4) = ⅓ K(1/√2)`.
* **The reflection-product identities are exact.** Including `J(1/4)J(3/4) =
  π/20`, which is Euler reflection, not the Legendre relation.
* **The Thales construction does not derive Ramanujan's 1/π series.**
* **It does not produce the constants 1103, 26390, 396, or 9801.**
* **It does not supply singular moduli, multipliers, the Legendre relation, or
  quasi-period data.**
* **The Clausen obstruction concerns the natural beta- and arcsine-moment
  generating functions tested here.** It is a demonstrated negative for those
  functions. It is *not* a proof that every conceivable construction is ruled
  out, and this repository does not claim otherwise.

Two further points worth stating, both reproduced in code:

* The construction reaches **two** CM fields, `Q(i)` at the quarter residue and
  `Q(√−3)` at the third and sixth residues, by the same elementary device. The
  failure of the analogy is therefore not that it selects the "wrong" quadratic
  field; the obstruction lies one level deeper, in the absence of quasi-period
  content in every stratum.
* The apparent poles of the product identity at `s = −1/2, 1/2, 3/2` are
  **removable**, since `cot(πs)` vanishes there (check `R3.5`).

The full division into *established contact*, *demonstrated boundary*, and *not
claimed* is in [`docs/scope_and_nonclaims.md`](docs/scope_and_nonclaims.md).
Read that file before drawing any conclusion from the outputs.

---

## Related material

This repository also hosts a self-contained technical note testing a separate
question: whether this repository's own elliptic/modular-period geometry (or
ordinary strain-eigenvalue anisotropy) adds information beyond the scalar
strain–rotation partition used in the author's separate fluid-mechanics note,
RVP. It is a documentation record of a completed falsification-first audit,
not part of the Ramanujan/Thales manuscript above, and it does not modify any
claim in it. See
[`docs/RVP_GEOMETRY_ALIGNMENT_AUDIT.md`](docs/RVP_GEOMETRY_ALIGNMENT_AUDIT.md) —
tests whether ellipticity, strain-eigenvalue shape, or Thales
complementary-modulus geometry adds information beyond the scalar
strain–rotation partition.

---

## Installation

Python 3.11 or later.

```bash
git clone https://github.com/innerlightr-wq/thales-ramanujan-quarter-gamma.git   # ← replace innerlightr-wq
cd thales-ramanujan-quarter-gamma
pip install -r requirements.txt
```

With conda:

```bash
conda env create -f environment.yml
conda activate thales-ramanujan
```

Optionally, install the package itself (not required — the scripts and tests
both bootstrap `sys.path` and run from a bare checkout):

```bash
pip install -e .
```

---

## Reproduction

Run every check and regenerate all result tables:

```bash
python scripts/run_all.py
```

Exit status is `0` if every check passes and `1` otherwise, so the runner can
be used directly in CI.

Run the test suite:

```bash
pytest
```

Run an individual group:

```bash
python scripts/verify_ramanujan_coefficients.py
python scripts/verify_beta_gamma_values.py
python scripts/verify_reflection_products.py
python scripts/verify_lemniscatic_relation.py
python scripts/verify_arcsine_moments.py
python scripts/verify_mirror_tilt.py
python scripts/verify_clausen_obstruction.py
python scripts/verify_response_hierarchy.py
```

Each script prints, for every identity: the statement being verified, the
symbolic status, the absolute error, and the declared tolerance. Numerical
checks run at 60 decimal digits (closed-form evaluations) or 40 digits (direct
quadrature of singular integrands); both regimes and their tolerances are
declared in `thales_ramanujan/precision.py` and reported per check. Exact
rational and integer arithmetic is used wherever possible, and floating-point
equality is never used. See
[`docs/computational_methods.md`](docs/computational_methods.md).

Outputs are deterministic — no timestamps, fixed row and column order — so
running the runner twice and diffing `results/` yields no changes.

---

## Repository structure

```
thales-ramanujan-quarter-gamma/
├── README.md
├── CITATION.cff
├── LICENSE
├── conftest.py                  # makes the package importable to pytest
├── environment.yml
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── thales_ramanujan/            # the mathematics, importable
│   ├── __init__.py
│   ├── precision.py             # declared precision and tolerance policy
│   ├── core.py                  # symbolic and numeric objects
│   ├── checks.py                # the eight verification groups
│   └── report.py                # records, formatting, deterministic CSV output
├── scripts/
│   ├── verify_ramanujan_coefficients.py
│   ├── verify_beta_gamma_values.py
│   ├── verify_reflection_products.py
│   ├── verify_lemniscatic_relation.py
│   ├── verify_arcsine_moments.py
│   ├── verify_mirror_tilt.py
│   ├── verify_clausen_obstruction.py
│   ├── verify_response_hierarchy.py
│   └── run_all.py
├── tests/
│   ├── test_coefficients.py
│   ├── test_beta_gamma.py
│   ├── test_reflection_products.py
│   ├── test_clausen_obstruction.py
│   └── test_response_hierarchy.py
├── results/
│   ├── README.md
│   ├── exact_identity_summary.csv
│   ├── high_precision_checks.csv
│   └── clausen_parameter_audit.csv
└── docs/
    ├── mathematical_scope.md
    ├── computational_methods.md
    ├── scope_and_nonclaims.md
    ├── zenodo_record.md
    └── RVP_GEOMETRY_ALIGNMENT_AUDIT.md   # separate RVP geometry/alignment audit record
```

The mathematics lives in the `thales_ramanujan/` package and the eight
`scripts/verify_*.py` files are thin wrappers around it. This is a small
departure from a flat `scripts/`-only layout, and the reason is that it lets
`pytest` exercise exactly the code paths `run_all.py` runs rather than a
parallel reimplementation of them.

---

## Citation

To cite **this repository**, see [`CITATION.cff`](CITATION.cff) (GitHub renders
a "Cite this repository" button from it).

To cite the **work**, cite the manuscript:

```bibtex
@misc{dejesus2026thales_ramanujan,
  author       = {De Jes{\'u}s, Elias},
  title        = {Ramanujan's Formula for 1/$\pi$ and the Thales Integral Witness},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20665588},
  url          = {https://doi.org/10.5281/zenodo.20665588},
  note         = {Concept DOI, resolving to the latest version.
                  Version reproduced by this repository: 10.5281/zenodo.20672985}
}
```

---

## License

Code and documentation in this repository are released under the MIT License;
see [`LICENSE`](LICENSE). MIT is applied to the entire repository, including
the generated tables under `results/`, to keep licensing unambiguous.

The manuscript is not included here and remains governed by its own Zenodo
license (CC BY 4.0).

---

## AI-assistance acknowledgment

Consistent with the manuscript: AI assistance was used for symbolic checking,
skeptical review, literature organization, and manuscript drafting, and it was
used in constructing this repository. All conceptual development,
interpretation, and final responsibility for the content belong to the author.
