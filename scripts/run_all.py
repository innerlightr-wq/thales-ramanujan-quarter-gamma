#!/usr/bin/env python3
"""Run every verification group and regenerate the result tables.

Run with:

    python scripts/run_all.py

Outputs (written deterministically, no timestamps, fixed row order):

    results/exact_identity_summary.csv   one row per identity: exact form and
                                         symbolic status
    results/high_precision_checks.csv    one row per identity: declared
                                         precision, declared tolerance, and
                                         reported absolute error
    results/clausen_parameter_audit.csv  hypergeometric parameters of each
                                         tested generating function against
                                         the Clausen condition c = a + b + 1/2

Exit status is 0 if every check passes and 1 otherwise.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import REGISTRY, clausen_audit_rows  # noqa: E402
from thales_ramanujan.report import (  # noqa: E402
    any_failed,
    format_records,
    results_dir,
    summarise,
    write_csv,
    write_rows,
)

EXACT_FIELDS = [
    "check_id",
    "group",
    "description",
    "exact_form",
    "symbolic_status",
    "status",
]

NUMERIC_FIELDS = [
    "check_id",
    "group",
    "description",
    "numeric_dps",
    "tolerance",
    "abs_error",
    "status",
]

CLAUSEN_FIELDS = [
    "generating_function",
    "hypergeometric_form",
    "a",
    "b",
    "c",
    "clausen_condition_applies",
    "required_c",
    "clausen_condition_holds",
    "note",
]


def main() -> int:
    records = []
    for group, script, runner in REGISTRY:
        print("=" * 78)
        print(f"{group}   [{script}]")
        print("=" * 78)
        group_records = runner()
        print(format_records(group_records))
        print(summarise(group_records))
        print()
        records.extend(group_records)

    out = results_dir()
    exact_path = os.path.join(out, "exact_identity_summary.csv")
    numeric_path = os.path.join(out, "high_precision_checks.csv")
    clausen_path = os.path.join(out, "clausen_parameter_audit.csv")

    write_csv(exact_path, records, EXACT_FIELDS)
    write_csv(
        numeric_path,
        [r for r in records if r.abs_error != ""],
        NUMERIC_FIELDS,
    )
    write_rows(clausen_path, CLAUSEN_FIELDS, clausen_audit_rows())

    print("=" * 78)
    print("wrote:")
    for path in (exact_path, numeric_path, clausen_path):
        print("  " + os.path.relpath(path, os.path.dirname(out)))
    print()
    print("TOTAL: " + summarise(records))

    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
