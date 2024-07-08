from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import numpy as np
import pyg4ometry
from legendhpges import make_hpge
from legendmeta import AttrsDict
from pyg4ometry import geant4

from . import materials
from .det_utils import RemageDetectorInfo

log = logging.getLogger(__name__)


def place_hpge_strings(
    channelmap: str | dict | AttrsDict,
    string_config: str | dict | AttrsDict,
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> None:
    """Construct LEGEND-200 HPGe strings.

    Parameters
    ----------
    channelmap
        LEGEND-200 HPGe channel map containing germanium detectors
        configuration in the string and their geometry.
    string_config
        LEGEND-200 germanium detector string configuration file.
        Used to reconstruct the spatial position of each string.
    z0
        The z coordinate of the top face of the array top plate.
    mothervolume
        pyg4ometry Geant4 LogicalVolume instance in which the strings
        are to be placed.
    registry
        pyg4ometry Geant4 registry instance.
    """
    if registry is None:
        msg = "registry cannot be None"
        raise ValueError(msg)

    if channelmap is None:
        msg = "configuration metadata file cannot be None"
        raise ValueError(msg)

    if string_config is None:
        msg = "string configuration cannot be None"
        raise ValueError(msg)

    if not isinstance(channelmap, (dict, AttrsDict)):
        with Path(channelmap).open() as jfile:
            ch_map = AttrsDict(json.load(jfile))
    else:
        ch_map = AttrsDict(channelmap)

    if not isinstance(string_config, (dict, AttrsDict)):
        with Path(string_config).open() as jfile:
            hpge_string_config = AttrsDict(json.load(jfile))
    else:
        hpge_string_config = AttrsDict(string_config)

    # derive the strings from the channelmap.
    ch_map = ch_map.map("system", unique=False).geds.values()

    strings_to_build = {}

    for hpge_meta in ch_map:
        # Temporary fix for gedet with null enrichment value
        if hpge_meta.production.enrichment is None:
            log.warning("%s has no enrichment in metadata - setting to dummy value 0.86!", hpge_meta.name)
            hpge_meta.production.enrichment = 0.86

        hpge_string_id = str(hpge_meta.location.string)
        hpge_unit_id_in_string = hpge_meta.location.position

        if hpge_string_id not in strings_to_build:
            strings_to_build[hpge_string_id] = {}

        hpge_extra_meta = hpge_string_config.hpges[hpge_meta.name]
        strings_to_build[hpge_string_id][hpge_unit_id_in_string] = HPGeDetUnit(
            hpge_meta.name,
            hpge_meta.production.manufacturer,
            hpge_meta.daq.rawid,
            make_hpge(hpge_meta, registry),
            hpge_meta.geometry.height_in_mm,
            hpge_extra_meta["baseplate"],
            hpge_extra_meta["rodlength_in_mm"],
        )

    # now, build all strings.
    for string_id, string in strings_to_build.items():
        _place_hpge_string(string_id, string, hpge_string_config, z0, mothervolume, materials, registry)


@dataclass
class HPGeDetUnit:
    name: str
    manufacturer: str
    rawid: int
    lv: geant4.LogicalVolume
    height: float
    baseplate: str
    rodlength: float


def _place_hpge_string(
    string_id: str,
    string_slots: list,
    hpge_string_config: AttrsDict,
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
):
    """
    Place a single HPGe detector string.

    This includes all PEN plates and the nylon shroud around the string."""
    string_meta = hpge_string_config.hpge_string[string_id]

    angle_in_rad = math.pi * string_meta.angle_in_deg / 180
    x_pos = string_meta.radius_in_mm * math.cos(angle_in_rad)
    y_pos = -string_meta.radius_in_mm * math.sin(angle_in_rad)
    # rotation angle for anything in the string.
    string_rot = -np.pi + angle_in_rad
    string_rot_m = np.array(
        [[np.sin(string_rot), np.cos(string_rot)], [np.cos(string_rot), -np.sin(string_rot)]]
    )

    # offset the height of the string by the length of the string support rod.
    support = _get_support_structure(materials, registry)
    # z0_string is the upper z coordinate of the topmost detector unit.
    z0_string = z0 - 422.11  # from CAD model.

    # deliberately use max and range here. The code does not support sparse strings (i.e. with
    # unpopulated slots, that are _not_ at the end. In those cases it should produce a KeyError.
    max_unit_id = max(string_slots.keys())
    total_rod_length = 0
    for hpge_unit_id_in_string in range(1, max_unit_id + 1):
        det_unit = string_slots[hpge_unit_id_in_string]

        # convert the "warm" length of the rod to the (shorter) length in the cooled down state.
        total_rod_length += det_unit.rodlength * 0.997

        z_unit_bottom = z0_string - total_rod_length
        # - notes for those comparing this to MaGe (those offsets are not from there, but from the
        #   CAD model): the detector unit (DU)-local z coordinates are inverted in comparison to
        #   the coordinates here, as well as to the string coordinates in MaGe.
        # - In MaGe, the end of the three support rods is at +11.1 mm, the PEN plate at +4 mm, the
        #   diode at -diodeHeight/2-0.025 mm, so that the crystal contact is at DU-z 0 mm.
        z_unit_pen = z_unit_bottom + 3.7 + 1.5 / 2  # from CAD model; 1.5 mm is the PEN thickness.
        # - note from CAD model: the distance between PEN plate and detector bottom face varies a
        #   lot between different diodes (i.e. BEGe's mostly (all?) have 2.1 mm face-to-face; for
        #   PPCs this varies between 2.5 and 4 mm.)
        z_pos_det = z_unit_pen + 4

        det_pv = geant4.PhysicalVolume(
            [0, 0, 0],
            [x_pos, y_pos, z_pos_det],
            det_unit.lv,
            det_unit.name,
            mothervolume,
            registry,
        )
        det_pv.pygeom_active_dector = RemageDetectorInfo("germanium", det_unit.rawid)
        det_unit.lv.pygeom_color_rgba = (0, 1, 1, 1)

        # a lot of Ortec detectors have modified medium plates.
        if (
            det_unit.name.startswith("V")
            and det_unit.baseplate == "medium"
            and det_unit.manufacturer == "Ortec"
        ):
            # TODO: what is with "V01389A"?
            det_unit.baseplate = "medium_ortec"
        pen_plate = _get_pen_plate(det_unit.baseplate, materials, registry)
        geant4.PhysicalVolume(
            [0, 0, string_rot],
            [x_pos, y_pos, z_unit_pen],
            pen_plate,
            det_unit.name + "_pen",
            mothervolume,
            registry,
        )

        # (Majorana) PPC detectors have a top PEN ring.
        if det_unit.name.startswith("P"):
            assert det_unit.baseplate == "small"
            pen_plate = _get_pen_plate("ppc_small", materials, registry)
            geant4.PhysicalVolume(
                [0, 0, string_rot],
                [x_pos, y_pos, z_pos_det + det_unit.height + 1.5 / 2],
                pen_plate,
                det_unit.name + "_pen_top",
                mothervolume,
                registry,
            )

    # TODO: offset 6 is from MaGe. This is quite certainly incorrect, the mini shrouds extend above the string!
    shroud_length = total_rod_length + 6
    ms = _get_nylon_mini_shroud(
        string_meta.minishroud_radius_in_mm, shroud_length, False, materials, registry
    )
    geant4.PhysicalVolume(
        [0, 0, 0],
        [x_pos, y_pos, z0_string - shroud_length / 2 + 0.1],  # add the shroud thickness to avoid overlaps.
        ms,
        ms.name + "_string_" + string_id,
        mothervolume,
        registry,
    )

    geant4.PhysicalVolume(
        [0, 0, np.deg2rad(30) + string_rot],
        [x_pos, y_pos, z0_string + 12],  # this offset of 12 is measured from the CAD file.
        support,
        support.name + "_string_" + string_id,
        mothervolume,
        registry,
    )

    copper_rod_r = string_meta.rod_radius_in_mm
    assert copper_rod_r < string_meta.minishroud_radius_in_mm - 0.75
    copper_rod_name = f"string_{string_id}_cu_rod"
    copper_rod_length = total_rod_length + 3.5  # the copper rod is slightly longer after the last detector.
    # the rod has a radius of 1.5 mm, but this would overlap with the coarse model of the PPC top PEN ring.
    copper_rod = geant4.solid.Tubs(copper_rod_name, 0, 1.43, copper_rod_length, 0, 2 * math.pi, registry)
    copper_rod = geant4.LogicalVolume(copper_rod, materials.metal_copper, copper_rod_name, registry)
    copper_rod.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    for i in range(3):
        copper_rod_th = np.deg2rad(-30 - i * 120)
        delta = copper_rod_r * string_rot_m @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        geant4.PhysicalVolume(
            [0, 0, 0],
            [x_pos + delta[0], y_pos + delta[1], z0_string - copper_rod_length / 2],
            copper_rod,
            f"{copper_rod_name}_{i}",
            mothervolume,
            registry,
        )


_pen_plate_cache = {}
_minishroud_cache = {}


def _get_pen_plate(
    size: str,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    if size not in ["small", "medium", "medium_ortec", "large", "xlarge", "ppc_small"]:
        msg = f"Invalid PEN-plate size {size}"
        raise ValueError(msg)

    # just for vis purposes...
    colors = {
        "small": (1, 0, 0, 1),
        "medium": (0, 1, 0, 1),
        "medium_ortec": (1, 0, 1, 1),
        "large": (0, 0, 1, 1),
        "xlarge": (1, 1, 0, 1),
        "ppc_small": (1, 0, 0, 1),
    }

    if size not in _pen_plate_cache:
        if size != "ppc_small":
            pen_file = resources.files("l200geom") / "models" / f"BasePlate_{size}.stl"
        else:
            pen_file = resources.files("l200geom") / "models" / "TopPlate_ppc.stl"

        pen_solid = pyg4ometry.stl.Reader(
            pen_file, solidname=f"pen_{size}", centre=False, registry=registry
        ).getSolid()
        _pen_plate_cache[size] = geant4.LogicalVolume(pen_solid, materials.pen, f"pen_{size}", registry)
        _pen_plate_cache[size].pygeom_color_rgba = colors[size]

    return _pen_plate_cache[size]


def _get_support_structure(
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    """This model's origin is a the top face of the tripod structure."""
    if "string_support_structure" not in registry.solidDict:
        support_file = resources.files("l200geom") / "models" / "StringSupportStructure.stl"
        support_solid = pyg4ometry.stl.Reader(
            support_file, solidname="string_support_structure", centre=False, registry=registry
        ).getSolid()
        support_lv = geant4.LogicalVolume(
            support_solid, materials.metal_copper, "string_support_structure", registry
        )
        support_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        support_solid = registry.solidDict["string_support_structure"]
        support_lv = registry.logicalVolumeDict["string_support_structure"]

    return support_lv


def _get_nylon_mini_shroud(
    radius: int,
    length: int,
    top_open: bool,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    """Create a nylon/TPB funnel of the given outer dimensions, which will be closed at the top/bottom.

    .. note:: this can also be used for calibration tubes.
    """
    shroud_name = f"minishroud_{radius}x{length}"
    if shroud_name not in _minishroud_cache:
        shroudThickness = 0.125  # mm
        endThickness = 2 * shroudThickness
        outer = geant4.solid.Tubs(f"{shroud_name}_outer", 0, radius, length, 0, 2 * math.pi, registry)
        inner = geant4.solid.Tubs(
            f"{shroud_name}_inner",
            0,
            radius - shroudThickness,
            # at the top/bottom, the NMS has essentially two layers.
            length - (0 if top_open else 2 * endThickness),
            0,
            2 * math.pi,
            registry,
        )
        # subtract the slightly smaller solid from the larger one, to get a hollow and closed volume.
        inner_z = (1 if top_open else 0) * endThickness
        shroud = geant4.solid.Subtraction(shroud_name, outer, inner, [[0, 0, 0], [0, 0, inner_z]], registry)
        _minishroud_cache[shroud_name] = geant4.LogicalVolume(
            shroud,
            materials.tpb_on_nylon,
            shroud_name,
            registry,
        )
        _minishroud_cache[shroud_name].pygeom_color_rgba = (1, 0.86, 0.86, 0.2)

    # TODO: implement optical surfaces
    return _minishroud_cache[shroud_name]
