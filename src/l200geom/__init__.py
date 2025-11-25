from __future__ import annotations

import sys
import warnings

import pygeoml200  # noqa: F401

sys.modules[__name__] = sys.modules["pygeoml200"]

warnings.warn("Please use `pygeoml200` instead of `l200geom`.", FutureWarning, stacklevel=2)
