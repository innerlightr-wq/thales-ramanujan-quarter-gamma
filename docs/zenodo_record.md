# Zenodo record

The manuscript is archived on Zenodo. **Zenodo is the authoritative archival
location for the paper.** This repository contains no manuscript PDF and no
LaTeX source, by design; it is the reproducibility, verification,
documentation, and code companion only.

## Citation

> De Jesús, Elias. (2026). *Ramanujan's Formula for 1/π and the Thales Integral
> Witness*. Zenodo. https://doi.org/10.5281/zenodo.20665588

ORCID: [0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)

## Identifiers

| Kind | DOI | Resolves to |
|---|---|---|
| **Concept DOI** (preferred) | [10.5281/zenodo.20665588](https://doi.org/10.5281/zenodo.20665588) | all versions; always the latest |
| Version DOI (reproduced here) | [10.5281/zenodo.20672985](https://doi.org/10.5281/zenodo.20672985) | the specific deposited version whose identities this repository verifies |

Use the **concept DOI** when citing the work. Use the **version DOI** when you
need to state exactly which deposit a computation corresponds to — which is the
case for this repository, since the check identifiers below are tied to the
numbered results of that specific version.

## Mapping from manuscript results to checks

| Manuscript | Statement | Checks |
|---|---|---|
| Lemma 1, eq. (2) | `A_n = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3` | `R1.1` |
| eq. (3) | `F(z) = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z)` | `R1.2` |
| eq. (1) | Ramanujan's series for `1/π` (orientation anchor) | `R1.3` |
| Lemma 2, eq. (6) | `J(s) = B(s+1,s+1) = Γ(s+1)²/Γ(2s+2)` | `R2.1` |
| §3, strata (7)–(10) | stratum recurrence and the four residue strata | `R2.2`, `R2.6` |
| eq. (11) | `J(1/4)`, `J(3/4)` in closed form | `R2.3`, `R2.4` |
| §3 | `J(1/2) = π/8` | `R2.5` |
| Proposition 1, eq. (12) | conjugate-stratum product identity | `R3.1`, `R3.5` |
| eq. (13), Proposition 2, eq. (20) | `π/20`, `π√3/35`, `15π√3/512` | `R3.2`–`R3.4` |
| eq. (14) | `K(1/√2) = Γ(1/4)²/(4√π)` | `R4.1` |
| eq. (15) | `J(1/4) = (1/3)K(1/√2)`; `J(3/4) = 3π/(20K)` | `R4.2`, `R4.3` |
| Lemma 3, eq. (16) | arcsine moments | `R5.1` |
| eq. (17) | `A_n = C(4n,2n) C(2n,n)²` | `R5.2`, `R5.3` |
| Lemma 4, eq. (18) | mirror-tilt ratio `α/(s+1)` | `R6.1`, `R6.2` |
| eq. (19) | digamma parameter response | `R6.3` |
| Proposition 3, eq. (22) | `G_beta`, `G_arc`, and the Clausen failure | `R7.1`–`R7.4` |
| eq. (4), eq. (21) | Clausen square (positive control) | `R7.5` |
| eq. (23) | `ψ(3/4) − ψ(1/4) = π` | `R8.1` |
| §A.3 | `ψ(3/2) − ψ(1) = 2 − 2 log 2` | `R8.2` |
| eq. (24) | `ψ'(1/4) = π² + 8G` | `R8.3`–`R8.5` |
| §7, §A.3 | second-order response of the family | `R8.6` |

## Licensing relationship

Code and documentation in this repository are MIT licensed (see `LICENSE`).
The manuscript itself remains governed by its Zenodo license (CC BY 4.0) and is
not redistributed here.
