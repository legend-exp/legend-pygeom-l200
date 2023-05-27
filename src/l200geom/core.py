from __future__ import annotations

import pathlib

from legendmeta import LegendMetadata
from pyg4ometry import geant4

from . import cryo, gearray, materials

lmeta = LegendMetadata()
config = pathlib.Path(__file__).parent.resolve() / "configs"


def construct() -> geant4.Registry:
    """Construct the LEGEND-200 geometry and return the pyg4ometry Registry containing the world volume."""
    reg = geant4.Registry()
    mats = materials.OpticalMaterialRegistry(reg)

    # Create the world volume
    world_material = geant4.MaterialPredefined("G4_Galactic")
    world = geant4.solid.Box("world", 20, 20, 20, reg, "m")
    world_lv = geant4.LogicalVolume(world, world_material, "world", reg)
    reg.setWorld(world_lv)

    # TODO: Shift the global coordinate system that z=0 is a reasonable value for defining hit positions.
    coordinate_z_displacement = 0

    # Create basic structure with argon and cryostat.
    cryostat_lv = cryo.construct_cryostat(mats.metal_steel, reg)
    cryo.place_cryostat(cryostat_lv, world_lv, coordinate_z_displacement, reg)

    lar_lv = cryo.construct_argon(mats.liquidargon, reg)
    cryo.place_argon(lar_lv, cryostat_lv, coordinate_z_displacement, reg)

    # Place the germanium detector array inside the liquid argon
    metadata = lmeta.hardware.detectors.germanium.diodes
    temp_ch_map = lmeta.hardware.configuration.channelmaps[r"l200-p03-r%-T%-all-config"]
    temp_array_config = config / "dummy_array_config.json"

    gearray.place_gearray(metadata, temp_array_config, temp_ch_map, 1950, lar_lv, reg)

    return reg
