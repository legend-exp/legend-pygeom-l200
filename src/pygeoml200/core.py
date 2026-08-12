from __future__ import annotations

import logging
from dataclasses import dataclass
from importlib import resources

from dbetto import AttrsDict, TextDB
from pyg4ometry import geant4
from pygeomtools.utils import load_dict_from_config

from . import calibration, cryo, fibers, hpge_strings, materials, top, watertank, wlsr
from .metadata import resolve_metadata

log = logging.getLogger(__name__)

configs = TextDB(resources.files("pygeoml200") / "configs" / "extra_meta")

DEFAULT_ASSEMBLIES = {"wlsr", "strings", "calibration", "fibers", "top"}
DEFINED_ASSEMBLIES = DEFAULT_ASSEMBLIES | {"watertank"}


@dataclass(frozen=True, slots=True, kw_only=True)
class InstrumentationData:
    mother_lv: geant4.LogicalVolume
    """Argon LogicalVolume instance in which all components are to be placed."""
    mother_pv: geant4.PhysicalVolume
    """Argon PhysicalVolume instance in which all components are to be placed."""
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
    assemblies: list[str] | set[str] = DEFAULT_ASSEMBLIES,
    use_detailed_fiber_model: bool = False,
    config: dict | None = None,
    public_geometry: bool = False,
) -> geant4.Registry:
    """Construct the LEGEND-200 geometry and return the pyg4ometry Registry containing the world volume."""
    if set(assemblies) - set(DEFINED_ASSEMBLIES) != set():
        msg = "invalid geometrical assembly specified"
        raise ValueError(msg)

    config = config if config is not None else {}

    meta = resolve_metadata(
        config,
        public_geometry,
        "cannot construct geometry from public testdata only, if not explicitly instructed",
    )

    reg = geant4.Registry()
    mats = materials.OpticalMaterialRegistry(reg)

    # Create the world volume
    world_material = geant4.MaterialPredefined("G4_Galactic")
    world = geant4.solid.Box("world", 20, 20, 20, reg, "m")
    world_lv = geant4.LogicalVolume(world, world_material, "world", reg)
    reg.setWorld(world_lv)

    # TODO: Shift the global coordinate system that z=0 is a reasonable value for defining hit positions.
    cryo_z_displacement: float = 0

    cryo_parent = world_lv
    if "watertank" in assemblies:
        # TODO: Shift the global coordinate system that z=0 is a reasonable value for defining hit positions.
        tank_z_displacement = 0.0
        cryo_z_displacement = (
            watertank.water_height / 2
            - cryo.cryo_access_height
            - (cryo.cryo_tub_height / 2 + cryo.cryo_top_height)
            - cryo.access_overlap / 2
            - 1e-9  # safety
        )  # -153
        tank_z_displacement = -cryo_z_displacement

        water_lv, water_pv, _ = watertank.insert_muon_veto(
            reg,
            world_lv,
            tank_z_displacement,
            cryo_z_displacement,
            mats,
        )
        cryo_parent = water_lv

    # Create basic structure with argon and cryostat.
    cryostat_lv = cryo.construct_cryostat(mats.metal_steel, reg)
    cryostat_pv = cryo.place_cryostat(cryostat_lv, cryo_parent, cryo_z_displacement, reg)

    if "watertank" in assemblies:
        geant4.BorderSurface(
            "water_cryo_surface", water_pv, cryostat_pv, mats.surfaces.vm2000_reflective_border, reg
        )

    argon_z_displacement = 0  # center argon in cryostat
    lar_lv, lar_neck_height = cryo.construct_argon(mats.liquidargon, reg)
    lar_pv = cryo.place_argon(
        lar_lv, cryostat_lv, cryostat_pv, argon_z_displacement, mats.surfaces.to_cryostat_steel, reg
    )
    gar_lv = cryo.construct_ullage_argon(mats.gaseousargon, reg)
    cryo.place_ullage_argon(gar_lv, cryostat_lv, argon_z_displacement, reg)

    array_total_height = 1488  # 1484 to 1490 mm array height (OB bottom to copper plate top).
    top_plate_z_pos_relative_to_neck = (
        7300  # end position meterdrive reading.
        - 1641  # meterdrive reading when OB touches shutter.
        - (1222 + 440 + 195.15 + 354 + 510)  # distance to shutter bottom flange.
        - (74 + (180 - 74) / 2)  # distance to the actual shutter surface.
        - array_total_height
    )
    top_plate_z_pos = lar_neck_height - top_plate_z_pos_relative_to_neck

    log.info(
        "displacement from cryostat center (positive to top): %f mm", top_plate_z_pos - array_total_height / 2
    )

    special_metadata = load_dict_from_config(config, "special_metadata", lambda: configs.on(meta.timestamp))
    special_metadata = meta.update_special_metadata(special_metadata)

    instr = InstrumentationData(
        mother_lv=lar_lv,
        mother_pv=lar_pv,
        materials=mats,
        registry=reg,
        channelmap=meta.channelmap,
        special_metadata=special_metadata,
        runtime_config=AttrsDict(config),
        top_plate_z_pos=top_plate_z_pos,
    )

    # Place all instrumentation into the liquid argon
    if "wlsr" in assemblies:
        # height below the lower end of the neck (even though this intended dimension is quite certainly
        # not really met in reality, P. Krause estimates ~cm uncertainty).
        wlsr.place_wlsr(instr, lar_neck_height - 1247.41)

    if "strings" in assemblies:
        hpge_strings.place_hpge_strings(meta.diodes, instr)
    if "calibration" in assemblies:
        calibration.place_calibration_system(instr)
    if "top" in assemblies:
        top.place_top_plate(instr)
    if "fibers" in assemblies:
        fibers.place_fiber_modules(meta.fibers, instr, use_detailed_fiber_model)

    _assign_common_copper_surface(instr)

    return reg


def _assign_common_copper_surface(b: InstrumentationData) -> None:
    """Assign a common copper surface to all copper parts in the LAr volume."""
    # check that the copper material has been instantiated (i.e. we have a copper part).
    if not hasattr(b.materials, "_metal_copper"):
        return

    surf = None
    cu_mat = b.materials.metal_copper

    # collect all existing border surfaces once, so that the check below stays O(1) per part.
    existing_borders = {
        (s.physref1, s.physref2)
        for s in b.registry.surfaceDict.values()
        if isinstance(s, geant4.BorderSurface)
    }

    for pv in b.registry.physicalVolumeDict.values():
        if (
            pv.motherVolume != b.mother_lv
            or not hasattr(pv.logicalVolume, "material")
            or pv.logicalVolume.material != cu_mat
        ):
            continue

        # only lazy-load the copper surface when we encounter a copper part.
        if surf is None:
            surf = b.materials.surfaces.to_copper

        # check that we do not have another surface already at this boundary.
        if (b.mother_pv, pv) in existing_borders:
            continue

        geant4.BorderSurface("bsurface_lar_cu_" + pv.name, b.mother_pv, pv, surf, b.registry)
