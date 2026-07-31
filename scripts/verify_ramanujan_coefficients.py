#!/usr/bin/env python3
"""Verifies Ramanujan's quarter-parameter coefficient lattice.

Identities verified:

  R1.1  A_n = (4n)!/(n!)^4 = 256^n (1/4)_n (1/2)_n (3/4)_n / (1)_n^3
  R1.2  F(z) = sum_n A_n z^n = 3F2(1/4, 1/2, 3/4; 1, 1; 256 z)
  R1.3  Ramanujan's series (manuscript eq. (1)) sums to 1/pi

R1.3 is an orientation anchor only: it fixes the target that the elementary
Thales construction is compared against and does not reach.

Run with:

    python scripts/verify_ramanujan_coefficients.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_ramanujan_coefficients  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_ramanujan_coefficients()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
