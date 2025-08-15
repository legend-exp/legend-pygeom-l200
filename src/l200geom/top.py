from __future__ import annotations

import logging

from pyg4ometry import geant4

from . import core
from .utils import _read_model

log = logging.getLogger(__name__)

TOP_PLATE_THICKNESS = 3


def place_top_plate(b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 copper top plate."""
    plate = _read_model("TopPlate.stl", "top_plate", b.materials.metal_copper, b)
    if plate is None:
        return
    plate.pygeom_color_rgba = (0.72, 0.45, 0.2, 0.2)

    geant4.PhysicalVolume([0, 0, 0], [0, 0, b.top_plate_z_pos], plate, plate.name, b.mother_lv, b.registry)
