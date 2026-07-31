#!/usr/bin/env python3
"""Verifies the arcsine moments and the central-binomial decomposition.

Identities verified:

  R5.1  (1/pi) int_0^1 x^n / sqrt(x(1-x)) dx = 4^(-n) C(2n,n)
  R5.2  (4n)!/(n!)^4 = C(4n,2n) C(2n,n)^2
  R5.3  A_n = 256^n m_(2n) m_n^2, with m_k the arcsine moments of R5.1

These witness the building blocks of Ramanujan's coefficient.  They do not
derive the 1/pi series: a product of three moments is not itself a moment,
and the assembly requires the modular theory.

Run with:

    python scripts/verify_arcsine_moments.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_arcsine_moments  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_arcsine_moments()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
