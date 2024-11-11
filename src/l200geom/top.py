from __future__ import annotations

import logging
from importlib import resources

import pyg4ometry
from pyg4ometry import geant4

from . import core

log = logging.getLogger(__name__)

TOP_PLATE_THICKNESS = 3


def place_top_plate(b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 HPGe strings."""
    plate_file = resources.files("l200geom") / "models" / "TopPlate.stl"
    plate = pyg4ometry.stl.Reader(
        plate_file, solidname="top_plate", centre=False, registry=b.registry
    ).getSolid()
    plate = geant4.LogicalVolume(plate, b.materials.metal_copper, "top_plate", b.registry)
    plate.pygeom_color_rgba = (0.72, 0.45, 0.2, 0.2)

    geant4.PhysicalVolume([0, 0, 0], [0, 0, b.top_plate_z_pos], plate, plate.name, b.mother_lv, b.registry)
