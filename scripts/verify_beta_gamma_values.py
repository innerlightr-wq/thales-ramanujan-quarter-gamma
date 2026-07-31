#!/usr/bin/env python3
"""Verifies the Thales beta family and its gamma strata.

Identities verified:

  R2.1  J(s) = int_0^1 [a(1-a)]^s da = B(s+1,s+1) = Gamma(s+1)^2/Gamma(2s+2)
  R2.2  J(s+1)/J(s) = (s+1)/(2(2s+3))  (the stratum recurrence)
  R2.3  J(1/4) = Gamma(1/4)^2 / (12 sqrt(pi))
  R2.4  J(3/4) = 3 Gamma(3/4)^2 / (10 sqrt(pi))
  R2.5  J(1/2) = pi/8
  R2.6  J(0) = 1, J(1) = 1/6, J(-1/2) = pi

Run with:

    python scripts/verify_beta_gamma_values.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thales_ramanujan.checks import check_beta_gamma_values  # noqa: E402
from thales_ramanujan.report import any_failed, format_records, summarise  # noqa: E402


def main() -> int:
    print(__doc__.strip().split("Run with:")[0].strip())
    print("-" * 78)
    records = check_beta_gamma_values()
    print(format_records(records))
    print("-" * 78)
    print(summarise(records))
    return 1 if any_failed(records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
