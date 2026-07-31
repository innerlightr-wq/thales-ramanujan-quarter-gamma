#!/usr/bin/env python3
"""Verifies the Clausen obstruction for the tested Thales generating functions.

Identities verified:

  R7.1  G_beta(z) = sum J(n) z^n = 2F1(1, 1; 3/2; z/4)
  R7.2  (a,b,c) = (1, 1, 3/2) does not satisfy c = a + b + 1/2
  R7.3  G_arc(z) = sum 4^(-n) C(2n,n) z^n = (1-z)^(-1/2) = 1F0(1/2; ; z)
  R7.4  G_arc(z)^2 = (1-z)^(-1) is a degenerate square
  R7.5  positive control: (1/8, 3/8, 1) does satisfy the Clausen condition

Scope: this is a demonstrated obstruction for the natural elementary Thales
generating functions tested here.  It does not rule out every conceivable
construction.

Run with:

    python scripts/verify_clausen_obstruction.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_clausen_obstruction  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_clausen_obstruction()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
