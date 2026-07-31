"""Reproducibility companion for

    De Jesus, Elias (2026). *Ramanujan's Formula for 1/pi and the Thales
    Integral Witness*. Zenodo. https://doi.org/10.5281/zenodo.20665588

This package verifies the exact identities recorded in the technical note and
the demonstrated boundary at which the Thales/Ramanujan analogy stops.  It
makes no mathematical claim beyond the manuscript.
"""

from .precision import DPS, QUAD_DPS, QUAD_TOL, TOL, set_precision
from .report import Record, any_failed, format_records, summarise

__all__ = [
    "DPS",
    "QUAD_DPS",
    "QUAD_TOL",
    "TOL",
    "Record",
    "any_failed",
    "format_records",
    "set_precision",
    "summarise",
]

__version__ = "1.0.0"
