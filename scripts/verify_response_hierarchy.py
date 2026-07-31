#!/usr/bin/env python3
"""Verifies the parameter-response hierarchy of the beta family.

Identities verified:

  R8.1  psi(3/4) - psi(1/4) = pi
  R8.2  psi(3/2) - psi(1) = 2 - 2 log 2
  R8.3  psi'(1/4) = pi^2 + 8 G  (G = Catalan's constant)
  R8.4  psi'(1/4) + psi'(3/4) = 2 pi^2
  R8.5  psi'(1/4) - psi'(3/4) = 16 G
  R8.6  d^2/d(alpha)^2 log I_s(alpha) = psi'(s+alpha+1) + psi'(s-alpha+1)

These are beta-, digamma-, and trigamma-level identities.  None of them is
the modular derivative responsible for the factor 1103 + 26390 n.

Run with:

    python scripts/verify_response_hierarchy.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_response_hierarchy  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_response_hierarchy()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
