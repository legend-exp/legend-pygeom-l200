from __future__ import annotations

from importlib import resources
from typing import Callable, NamedTuple

from legendmeta import AttrsDict, LegendMetadata, TextDB
from pyg4ometry import geant4

from . import calibration, cryo, det_utils, fibers, hpge_strings, materials, top, vis_utils, wlsr

lmeta = LegendMetadata()
configs = TextDB(resources.files("l200geom") / "configs")

DEFINED_ASSEMBLIES = ["wlsr", "strings", "calibration", "fibers", "top"]


class InstrumentationData(NamedTuple):
    mother_lv: geant4.LogicalVolume
    """LogicalVolume instance in which all components are to be placed."""
    mother_pv: geant4.PhysicalVolume
    """PhysicalVolume instance in which all components are to be placed."""
    materials: materials.OpticalMaterialRegistry
    """Material properties for common materials"""
    registry: geant4.Registry
    """pyg4ometry registry instance."""

    channelmap: AttrsDict
    """LEGEND-200 channel map containing germanium/spms detectors configuration in the string
    and their geometry."""
    special_metadata: AttrsDict
    """LEGEND-200 special geometry metadata file. Used to reconstruct the spatial position of each
    string, detector and calibration tube."""
    runtime_config: AttrsDict
    """Volatile runtime config, settings that are not tied to a specific detector configuration."""

    top_plate_z_pos: float
    """The z coordinate of the top face of the array top plate."""


def construct(
    assemblies: list[str] = DEFINED_ASSEMBLIES,
    use_detailed_fiber_model: bool = False,
    config: dict | None = None,
) -> geant4.Registry:
    """Construct the LEGEND-200 geometry and return the pyg4ometry Registry containing the world volume."""
    if set(assemblies) - set(DEFINED_ASSEMBLIES) != set():
        msg = "invalid geometrical assembly specified"
        raise ValueError(msg)

    config = config if config is not None else {}

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

    # top of the top plate, this is still a dummy value!
    top_plate_z_pos = 1700

    timestamp = config.get("metadata_timestamp", "20230311T235840Z")
    channelmap = _load_map_from_config(config, "channelmap", lambda: lmeta.channelmap(timestamp))
    special_metadata = _load_map_from_config(config, "special_metadata", lambda: configs.on(timestamp))
    instr = InstrumentationData(
        lar_lv, lar_pv, mats, reg, channelmap, special_metadata, AttrsDict(config), top_plate_z_pos
    )

    if "wlsr" in assemblies:
        # Place the WLSR into the cryostat.
        wlsr_lvs = wlsr.construct_wlsr(
            mats.metal_copper,
            mats.tetratex,
            mats.tpb_on_tetratex,
            reg,
        )
        # TODO: the z offset here is still a dummy value?
        wlsr_pvs = wlsr.place_wlsr(*wlsr_lvs, lar_lv, 3 * 180, reg)
        wlsr.add_surfaces_wlsr(*wlsr_pvs[1:], instr)

    # Place all other instrumentation into the liquid argon
    if "strings" in assemblies:
        hpge_strings.place_hpge_strings(instr)
    if "calibration" in assemblies:
        calibration.place_calibration_system(instr)
    if "top" in assemblies:
        top.place_top_plate(instr)
    if "fibers" in assemblies:
        fibers.place_fiber_modules(lmeta.hardware.detectors.lar.fibers, instr, use_detailed_fiber_model)

    _assign_common_copper_surface(instr)

    det_utils.write_detector_auxvals(reg)
    vis_utils.write_color_auxvals(reg)

    return reg


def _load_map_from_config(config: dict, key: str, default: Callable[[], AttrsDict]) -> AttrsDict:
    m = config.get(key)
    if isinstance(m, str):
        import json
        from pathlib import Path

        with Path(m).open() as jfile:
            return AttrsDict(json.load(jfile))
    elif isinstance(m, dict):
        return AttrsDict(m)
    return default()


def _assign_common_copper_surface(b: InstrumentationData) -> None:
    if hasattr(b.materials, "_metal_copper") is None:
        return
    surf = None
    cu_mat = b.materials.metal_copper

    for _, pv in b.registry.physicalVolumeDict.items():
        if pv.motherVolume != b.mother_lv or pv.logicalVolume.material != cu_mat:
            continue
        if surf is None:
            surf = b.materials.surfaces.to_copper

        geant4.BorderSurface("bsurface_lar_cu_" + pv.name, b.mother_pv, pv, surf, b.registry)
        geant4.BorderSurface("bsurface_cu_lar_" + pv.name, pv, b.mother_pv, surf, b.registry)
