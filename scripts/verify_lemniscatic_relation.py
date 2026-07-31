#!/usr/bin/env python3
"""Verifies the lemniscatic-period relation for the quarter stratum.

Identities verified:

  R4.1  K(1/sqrt(2)) = Gamma(1/4)^2 / (4 sqrt(pi))
  R4.2  J(1/4) = (1/3) K(1/sqrt(2))
  R4.3  J(3/4) = 3 pi / (20 K(1/sqrt(2)))

This supplies a period only.  It supplies no multiplier and no quasi-period
response, and it selects the CM field Q(i), not the singular modulus that
fixes the constants of Ramanujan's series.

Run with:

    python scripts/verify_lemniscatic_relation.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_lemniscatic_relation  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_lemniscatic_relation()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
