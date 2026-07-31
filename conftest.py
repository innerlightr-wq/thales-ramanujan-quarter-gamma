"""Make the package importable when running ``pytest`` from a checkout that
has not been installed with ``pip install -e .``."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
