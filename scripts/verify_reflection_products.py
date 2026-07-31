#!/usr/bin/env python3
"""Verifies the conjugate-stratum product identity.

Identities verified:

  R3.1  J(s) J(1-s) = pi s(1-s) cot(pi s) / (2 (2s+1)(1-2s)(3-2s))
  R3.2  J(1/4) J(3/4) = pi/20
  R3.3  J(1/3) J(2/3) = pi sqrt(3)/35
  R3.4  J(1/6) J(5/6) = 15 pi sqrt(3)/512
  R3.5  the apparent poles at s = -1/2, 1/2, 3/2 are removable

The collapse to a rational multiple of pi is Euler reflection.  It is not an
instance of the Legendre relation and carries no quasi-period content.

Run with:

    python scripts/verify_reflection_products.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_reflection_products  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_reflection_products()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
