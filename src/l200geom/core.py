from __future__ import annotations

from importlib import resources

from legendmeta import JsonDB, LegendMetadata
from pyg4ometry import geant4

from . import cryo, hpge_strings, materials, wlsr

lmeta = LegendMetadata()
configs = JsonDB(resources.files("l200geom") / "configs")


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

    # Place the WLSR into the cryostat.
    wlsr_lvs = wlsr.construct_wlsr(
        mats.metal_copper,
        mats.tetratex,
        mats.tpb_on_tetratex,
        reg,
    )
    wlsr.place_wlsr(*wlsr_lvs, lar_lv, 3 * 180, reg)
    # TODO: implement and enable optical surfaces.
    # wlsr.add_surfaces_wlsr(*wlsr_lvs[1:], lar_lv, mats, reg)

    # Place the germanium detector array inside the liquid argon
    hpge_channelmap = lmeta.channelmap("20230311T235840Z")
    hpge_string_config = configs.on("20230311T235840Z")

    hpge_strings.place_hpge_strings(
        hpge_channelmap, hpge_string_config, 1950, lar_lv, reg
    )

    return reg
