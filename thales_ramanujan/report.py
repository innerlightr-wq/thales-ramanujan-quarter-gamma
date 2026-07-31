"""Uniform result records, formatting, and deterministic CSV output.

Every verified identity produces exactly one :class:`Record`, carrying both a
symbolic status and (where meaningful) an independent numerical audit with a
reported absolute error and a declared tolerance.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field, asdict
from typing import Iterable, Optional, Sequence

import sympy as sp
from mpmath import mp, mpf, nstr

from .precision import DPS, TOL

PASS = "PASS"
FAIL = "FAIL"

SYMBOLIC_VERIFIED = "verified"
SYMBOLIC_FAILED = "failed"
SYMBOLIC_NA = "not applicable"


@dataclass
class Record:
    """One verified (or failed) identity."""

    check_id: str
    group: str
    description: str
    exact_form: str = ""
    symbolic_status: str = SYMBOLIC_NA
    symbolic_residual: str = ""
    numeric_dps: int = 0
    tolerance: str = ""
    abs_error: str = ""
    status: str = PASS

    def as_row(self) -> dict:
        return asdict(self)


def fmt_error(err) -> str:
    """Deterministic short scientific formatting of an absolute error."""
    if err is None:
        return ""
    err = mp.fabs(mpf(err))
    if err == 0:
        return "0e+0"
    return nstr(err, 3)


def _symbolic_status(residual, simplifier=None) -> tuple[str, str]:
    """Reduce a sympy residual to (status, printed residual).

    ``simplifier`` may supply a stronger reduction than :func:`sympy.simplify`
    (for example one that applies ``expand_trig`` after rewriting cotangents),
    for residuals that SymPy does not close on its own.
    """
    if residual is None:
        return SYMBOLIC_NA, ""
    simplified = (simplifier or sp.simplify)(residual)
    if simplified == 0:
        return SYMBOLIC_VERIFIED, "0"
    return SYMBOLIC_FAILED, sp.sstr(simplified)


def check(
    check_id: str,
    group: str,
    description: str,
    exact_form: str = "",
    symbolic_residual=None,
    symbolic_note: Optional[str] = None,
    simplifier=None,
    lhs=None,
    rhs=None,
    tolerance: str = TOL,
    dps: int = DPS,
) -> Record:
    """Build a :class:`Record` from an optional symbolic residual and an
    optional independent numerical pair ``(lhs, rhs)``.

    The record fails if the symbolic residual does not reduce to zero, or if
    the absolute numerical error exceeds the declared tolerance.  Numerical
    comparison is by absolute error against a declared bound: floating-point
    equality is never used.
    """
    status = PASS

    if symbolic_note is not None:
        sym_status, sym_residual = symbolic_note, ""
        if symbolic_residual is not None:
            sym_status_auto, sym_residual = _symbolic_status(symbolic_residual, simplifier)
            if sym_status_auto == SYMBOLIC_FAILED:
                status = FAIL
                sym_status = f"{symbolic_note} [residual nonzero]"
    else:
        sym_status, sym_residual = _symbolic_status(symbolic_residual, simplifier)
        if sym_status == SYMBOLIC_FAILED:
            status = FAIL

    err_str = ""
    tol_str = ""
    used_dps = 0
    if lhs is not None and rhs is not None:
        used_dps = dps
        tol_str = tolerance
        err = mp.fabs(mpf(lhs) - mpf(rhs))
        err_str = fmt_error(err)
        if not (err <= mpf(tolerance)):
            status = FAIL

    return Record(
        check_id=check_id,
        group=group,
        description=description,
        exact_form=exact_form,
        symbolic_status=sym_status,
        symbolic_residual=sym_residual,
        numeric_dps=used_dps,
        tolerance=tol_str,
        abs_error=err_str,
        status=status,
    )


def format_records(records: Sequence[Record]) -> str:
    """Human-readable console block for a list of records."""
    lines = []
    width = max((len(r.check_id) for r in records), default=8)
    for r in records:
        bits = [f"[{r.status}] {r.check_id.ljust(width)}  {r.description}"]
        detail = [f"symbolic: {r.symbolic_status}"]
        if r.abs_error:
            detail.append(f"abs err: {r.abs_error} (tol {r.tolerance}, {r.numeric_dps} dps)")
        bits.append("           " + " | ".join(detail))
        if r.symbolic_residual not in ("", "0"):
            bits.append(f"           residual: {r.symbolic_residual}")
        lines.append("\n".join(bits))
    return "\n".join(lines)


def summarise(records: Sequence[Record]) -> str:
    failed = [r for r in records if r.status != PASS]
    return (
        f"{len(records)} checks, {len(records) - len(failed)} passed, "
        f"{len(failed)} failed"
    )


def any_failed(records: Iterable[Record]) -> bool:
    return any(r.status != PASS for r in records)


CSV_FIELDS = [
    "check_id",
    "group",
    "description",
    "exact_form",
    "symbolic_status",
    "numeric_dps",
    "tolerance",
    "abs_error",
    "status",
]


def write_csv(path: str, records: Sequence[Record], fields: Sequence[str] = CSV_FIELDS) -> None:
    """Write records deterministically: fixed field order, fixed row order,
    LF line endings, no timestamps."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        for record in records:
            row = record.as_row()
            writer.writerow({key: row.get(key, "") for key in fields})


def write_rows(path: str, fields: Sequence[str], rows: Sequence[dict]) -> None:
    """Write arbitrary dict rows deterministically (used for the Clausen audit)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def repo_root() -> str:
    """Absolute path of the repository root (parent of this package)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def results_dir() -> str:
    return os.path.join(repo_root(), "results")
