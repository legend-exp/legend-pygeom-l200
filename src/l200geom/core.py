from __future__ import annotations

from importlib import resources

from legendmeta import LegendMetadata, TextDB
from pyg4ometry import geant4

from . import calibration, cryo, fibers, hpge_strings, materials, top, wlsr

lmeta = LegendMetadata()
configs = TextDB(resources.files("l200geom") / "configs")

DEFINED_ASSEMBLIES = ["wlsr", "strings", "calibration", "fibers", "top"]


def construct(
    assemblies: list[str] = DEFINED_ASSEMBLIES,
    use_detailed_fiber_model: bool = False,
) -> geant4.Registry:
    """Construct the LEGEND-200 geometry and return the pyg4ometry Registry containing the world volume."""
    if set(assemblies) - set(DEFINED_ASSEMBLIES) != set():
        msg = "invalid geometrical assembly specified"
        raise ValueError(msg)

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
    lar_pv = cryo.place_argon(lar_lv, cryostat_lv, coordinate_z_displacement, reg)

    if "wlsr" in assemblies:
        # Place the WLSR into the cryostat.
        wlsr_lvs = wlsr.construct_wlsr(
            mats.metal_copper,
            mats.tetratex,
            mats.tpb_on_tetratex,
            reg,
        )
        wlsr_pvs = wlsr.place_wlsr(*wlsr_lvs, lar_lv, 3 * 180, reg)
        wlsr.add_surfaces_wlsr(*wlsr_pvs[1:], lar_lv, mats, reg)

    channelmap = lmeta.channelmap("20230311T235840Z")

    # Place the germanium detector array inside the liquid argon
    hpge_string_config = configs.on("20230311T235840Z")

    # top of the top plate, this is still a dummy value!
    top_plate_z_pos = 1700

    if "strings" in assemblies:
        hpge_strings.place_hpge_strings(channelmap, hpge_string_config, top_plate_z_pos, lar_lv, mats, reg)
    if "calibration" in assemblies:
        calibration.place_calibration_system(top_plate_z_pos, lar_lv, mats, reg)
    if "top" in assemblies:
        top.place_top_plate(top_plate_z_pos, lar_lv, mats, reg)

    # build fiber modules
    if "fibers" in assemblies:
        fiber_modules = lmeta.hardware.detectors.lar.fibers
        fibers.place_fiber_modules(
            fiber_modules, channelmap, top_plate_z_pos, lar_lv, lar_pv, mats, reg, use_detailed_fiber_model
        )

    return reg
