#!/usr/bin/env python3
"""Verifies the mirror-tilt diagnostic and the derivative boundary.

Identities verified:

  R6.1  A_s(alpha) / I_s(alpha) = alpha / (s+1)
  R6.2  the same ratio is 2 E[a] - 1 for a ~ Beta(s+alpha+1, s-alpha+1)
  R6.3  d/d(alpha) log I_s(alpha) = psi(s+alpha+1) - psi(s-alpha+1)

The algebraic tilt (R6.1) acts on the coordinate a; the Euler operator
theta = z d/dz of the Ramanujan side acts on the summation index n.  They
must not be conflated.  R6.3 is the honest parameter response, and it is
beta/digamma data throughout.

Run with:

    python scripts/verify_mirror_tilt.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_mirror_tilt  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_mirror_tilt()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
