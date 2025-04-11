from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from importlib import resources

import numpy as np
import pyg4ometry
from dbetto import AttrsDict, TextDB
from legendhpges import make_hpge
from pyg4ometry import geant4
from pygeomtools import RemageDetectorInfo
from scipy.spatial.transform import Rotation

from . import core, materials

log = logging.getLogger(__name__)


def place_hpge_strings(hpge_metadata: TextDB, b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 HPGe strings."""
    # derive the strings from the channelmap.
    ch_map = b.channelmap.map("system", unique=False).get("geds", {}).values()
    strings_to_build = {}

    for ch_meta in ch_map:
        # ch_meta might be a full channelmap entry, i.e. containing the merged hardware meta from
        # lmeta.channelmap(), or a shallow dict with only channel data. So combine them with the hardware
        # data again.
        hpge_meta = hpge_metadata[ch_meta.name]
        assert hpge_meta.name == ch_meta.name
        full_meta = ch_meta | hpge_meta

        # Temporary fix for gedet with null enrichment value
        if hpge_meta.production.enrichment is None:
            log.warning("%s has no enrichment in metadata - setting to dummy value 0.86!", hpge_meta.name)
            hpge_meta.production.enrichment = 0.86

        hpge_string_id = str(ch_meta.location.string)
        hpge_unit_id_in_string = ch_meta.location.position

        if hpge_string_id not in strings_to_build:
            strings_to_build[hpge_string_id] = {}

        hpge_extra_meta = b.special_metadata.hpges[hpge_meta.name]
        strings_to_build[hpge_string_id][hpge_unit_id_in_string] = HPGeDetUnit(
            hpge_meta.name,
            hpge_meta.production.manufacturer,
            ch_meta.daq.rawid,
            make_hpge(full_meta, b.registry),
            hpge_meta.geometry.height_in_mm,
            hpge_meta.geometry.radius_in_mm,
            hpge_extra_meta["baseplate"],
            hpge_extra_meta["rodlength_in_mm"],
            full_meta,
        )

    for string_id, string_meta in b.special_metadata.hpge_string.items():
        if string_meta.get("empty_string_content") is None:
            continue
        if string_id in strings_to_build:
            msg = f"string {string_id} has empty_string_content and detectors"
            raise RuntimeError(msg)
        _place_empty_string(string_id, b)

    # now, build all strings.
    for string_id, string in strings_to_build.items():
        _place_hpge_string(string_id, string, b)


@dataclass
class HPGeDetUnit:
    name: str
    manufacturer: str
    rawid: int
    lv: geant4.LogicalVolume
    height: float
    radius: float
    baseplate: str
    rodlength: float
    meta: AttrsDict


def _place_front_end_and_insulators(
    det_unit: HPGeDetUnit,
    unit_length: float,
    string_info: dict,
    b: core.InstrumentationData,
    z_pos: dict,
    thickness: dict,
    parts_origin: dict,
):
    # add cable and clamp
    signal_cable, signal_clamp, signal_asic = _get_signal_cable_and_asic(
        det_unit.name,
        thickness["cable"],
        thickness["clamp"],
        unit_length,
        b.materials,
        b.mother_lv,
        b.registry,
    )
    signal_cable.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    signal_clamp.pygeom_color_rgba = (0.3, 0.3, 0.3, 1)
    signal_asic.pygeom_color_rgba = (0.73, 0.33, 0.4, 1)

    angle_signal = math.pi * 1 / 2.0 - string_info["string_rot"]
    x_clamp, y_clamp = np.array([string_info["x_pos"], string_info["y_pos"]]) + parts_origin["signal"][
        "clamp"
    ] * np.array([np.sin(string_info["string_rot"]), np.cos(string_info["string_rot"])])
    x_cable, y_cable = np.array([string_info["x_pos"], string_info["y_pos"]]) + parts_origin["signal"][
        "cable"
    ] * np.array([np.sin(string_info["string_rot"]), np.cos(string_info["string_rot"])])
    x_asic, y_asic = np.array([string_info["x_pos"], string_info["y_pos"]]) + parts_origin["signal"][
        "asic"
    ] * np.array([np.sin(string_info["string_rot"]), np.cos(string_info["string_rot"])])

    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_cable, y_cable, z_pos["cable"]],  # this offset of 12 is measured from the CAD file.
        signal_cable,
        signal_cable.name + "_string_" + string_info["string_id"],
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_clamp, y_clamp, z_pos["clamp"]],  # this offset of 12 is measured from the CAD file.
        signal_clamp,
        signal_clamp.name + "_string_" + string_info["string_id"],
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [
            x_asic,
            y_asic,
            z_pos["cable"] - thickness["cable"] - 0.5,
        ],  # this offset of 12 is measured from the CAD file.
        signal_asic,
        signal_asic.name + "_string_" + string_info["string_id"],
        b.mother_lv,
        b.registry,
    )

    hv_cable, hv_clamp = _get_hv_cable(
        det_unit.name,
        thickness["cable"],
        thickness["clamp"],
        unit_length,
        b.materials,
        b.mother_lv,
        b.registry,
    )
    hv_cable.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    hv_clamp.pygeom_color_rgba = (0.3, 0.3, 0.3, 1)

    angle_hv = math.pi * 1 / 2.0 + string_info["string_rot"]
    x_clamp, y_clamp = np.array([string_info["x_pos"], string_info["y_pos"]]) - parts_origin["hv"][
        "clamp"
    ] * np.array([np.sin(string_info["string_rot"]), np.cos(string_info["string_rot"])])
    x_cable, y_cable = np.array([string_info["x_pos"], string_info["y_pos"]]) - parts_origin["hv"][
        "cable"
    ] * np.array([np.sin(string_info["string_rot"]), np.cos(string_info["string_rot"])])

    geant4.PhysicalVolume(
        [0, 0, angle_hv],
        [x_clamp, y_clamp, z_pos["cable"]],
        hv_cable,
        hv_cable.name + "_string_" + string_info["string_id"],
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [0, 0, angle_hv],
        [x_clamp, y_clamp, z_pos["clamp"]],
        hv_clamp,
        hv_clamp.name + "_string_" + string_info["string_id"],
        b.mother_lv,
        b.registry,
    )

    insulator_top_length = string_info["string_meta"].rod_radius_in_mm - det_unit.radius + 1.5

    click, insulator = _get_click_and_insulator(
        det_unit,
        thickness["click"],
        thickness["insulator"],
        insulator_top_length,
        b.materials,
        b.registry,
    )
    click.pygeom_color_rgba = (0.6, 0.6, 0.6, 1)
    insulator.pygeom_color_rgba = (0.6, 0.6, 0.6, 1)

    for i in range(3):
        copper_rod_th = np.deg2rad(-30 - i * 120)
        pieces_th = string_info["string_rot"] + np.deg2rad(-(i + 1) * 120)
        delta_click = (
            (string_info["string_meta"].rod_radius_in_mm - 5.6)
            * string_info["string_rot_m"]
            @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        )
        delta_insulator = (
            (string_info["string_meta"].rod_radius_in_mm - (16.5 / 2.0 - 1.5))
            * string_info["string_rot_m"]
            @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        )
        geant4.PhysicalVolume(
            [0, 0, pieces_th],
            [string_info["x_pos"] + delta_click[0], string_info["y_pos"] + delta_click[1], z_pos["click"]],
            click,
            f"{click.name}_{i}",
            b.mother_lv,
            b.registry,
        )
        geant4.PhysicalVolume(
            [0, 0, pieces_th],
            [
                string_info["x_pos"] + delta_insulator[0],
                string_info["y_pos"] + delta_insulator[1],
                z_pos["insulator"],
            ],
            insulator,
            f"{insulator.name}_{i}",
            b.mother_lv,
            b.registry,
        )


def _place_hpge_unit(
    z_unit_bottom: float,
    det_unit: HPGeDetUnit,
    unit_length: float,
    string_info: dict,
    thicknesses: dict,
    b: core.InstrumentationData,
):
    safety_margin = 0.001  # 0.001 # 1 micro meter

    pen_offset = -0.15  # mm

    z_pos = {
        "det": z_unit_bottom,
        "insulator": z_unit_bottom - thicknesses["insulator"] / 2.0 - safety_margin,
        "pen": z_unit_bottom
        - thicknesses["insulator"]
        - thicknesses["pen"] / 2.0
        - pen_offset
        - safety_margin * 2,
        "click": z_unit_bottom
        - thicknesses["insulator"]
        - thicknesses["pen"]
        - thicknesses["click"] / 2.0
        - safety_margin * 3,
        "cable": z_unit_bottom
        - thicknesses["insulator"]
        - thicknesses["pen"]
        - thicknesses["cable"] / 2.0
        - safety_margin * 3,
        "clamp": z_unit_bottom
        - thicknesses["insulator"]
        - thicknesses["pen"]
        - thicknesses["cable"]
        - thicknesses["clamp"] / 2.0
        - safety_margin * 4,
        "pen_top": z_unit_bottom
        + det_unit.height
        + thicknesses["insulator"]
        + thicknesses["pen"]
    }

    det_pv = geant4.PhysicalVolume(
        [0, 0, 0],
        [string_info["x_pos"], string_info["y_pos"], z_pos["det"]],
        det_unit.lv,
        det_unit.name,
        b.mother_lv,
        b.registry,
    )
    det_pv.pygeom_active_detector = RemageDetectorInfo("germanium", det_unit.rawid, det_unit.meta)
    det_unit.lv.pygeom_color_rgba = (0, 1, 1, 1)

    # add germanium reflective surface.
    geant4.BorderSurface(
        "bsurface_lar_ge_" + det_pv.name,
        b.mother_pv,
        det_pv,
        b.materials.surfaces.to_germanium,
        b.registry,
    )

    baseplate = det_unit.baseplate
    # a lot of Ortec detectors have modified medium plates.
    if (
        det_unit.name.startswith("V")
        and det_unit.baseplate == "medium"
        and det_unit.manufacturer == "Ortec"
    ):
        # TODO: what is with "V01389A"?
        baseplate = "medium_ortec"
    pen_plate = _get_pen_plate(baseplate, b.materials, b.registry)

    # This rotation is not physical, but gets us closer to the real model of the PEN plates.
    # In the CAD model, most plates are mirrored, compared to reality (some are also correct in the
    # first place), i.e. how the plates in PGT were produced. So the STL mesh is also mirrored, so
    # flip it over.
    # note/TODO: this rotation should be replaced by a correct mesh, so that the counterbores are
    # on the correct side. This might be necessary to fit in other parts!
    pen_rot = Rotation.from_euler("XZ", [-math.pi, string_info["string_rot"]]).as_euler("xyz")
    pen_pv = geant4.PhysicalVolume(
        list(pen_rot),
        [string_info["x_pos"], string_info["y_pos"], z_pos["pen"]],
        pen_plate,
        det_unit.name + "_pen",
        b.mother_lv,
        b.registry,
    )
    _add_pen_surfaces(pen_pv, b.mother_pv, b.materials, b.registry)

    # (Majorana) PPC detectors have a top PEN ring.
    if det_unit.name.startswith("P"):
        assert det_unit.baseplate == "small"
        pen_plate = _get_pen_plate("ppc_small", b.materials, b.registry)
        pen_pv = geant4.PhysicalVolume(
            list(pen_rot),
            [string_info["x_pos"], string_info["y_pos"], z_pos["pen_top"]],
            pen_plate,
            "pen_top_" + det_unit.name,
            b.mother_lv,
            b.registry,
        )
        _add_pen_surfaces(pen_pv, b.mother_pv, b.materials, b.registry)

    front_enc_and_insulator_parts_origin = {
        "signal": {
            "clamp": 2.5
            + 4.0
            + 1.5
            + 5 / 2,  # position from center of detector to center of volume center
            "cable": 2.5 + 4.0 + 16 / 2,
            "asic": 2.5 + 4.0 + 11 + 1 / 2.0,
        },
        "hv": {"clamp": 2.5 + 29.5 + 3.5 + 5 / 2, "cable": 2.5 + 29.5 + 2.0 + 8 / 2},
    }

    _place_front_end_and_insulators(
        det_unit, unit_length, string_info, b, z_pos, thicknesses, front_enc_and_insulator_parts_origin
    )





def _place_hpge_string(
    string_id: str,
    string_slots: list,
    b: core.InstrumentationData,
):
    """Place a single HPGe detector string (with at least one detector).

    This includes all PEN plates and the nylon shroud around the string."""
    string_meta = b.special_metadata.hpge_string[string_id]

    angle_in_rad = math.pi * string_meta.angle_in_deg / 180
    x_pos = string_meta.radius_in_mm * math.cos(angle_in_rad)
    y_pos = -string_meta.radius_in_mm * math.sin(angle_in_rad)
    # rotation angle for anything in the string.
    string_rot = -np.pi + angle_in_rad
    string_rot_m = np.array(
        [[np.sin(string_rot), np.cos(string_rot)], [np.cos(string_rot), -np.sin(string_rot)]]
    )

    # offset the height of the string by the length of the string support rod.
    # z0_string is the upper z coordinate of the topmost detector unit.
    # TODO: real measurements (slides of M. Bush on 2024-07-08) show an additional offset -0.6 mm.
    # TODO: this is also still a warm length.
    z0_string = b.top_plate_z_pos - 410.1 - 12  # from CAD model.

    # deliberately use max and range here. The code does not support sparse strings (i.e. with
    # unpopulated slots, that are _not_ at the end. In those cases it should produce a KeyError.
    max_unit_id = max(string_slots.keys())
    total_rod_length = 0
    for hpge_unit_id_in_string in range(1, max_unit_id + 1):
        det_unit = string_slots[hpge_unit_id_in_string]

        # convert the "warm" length of the rod to the (shorter) length in the cooled down state.
        total_rod_length += det_unit.rodlength * 0.997

        z_unit_bottom = z0_string - total_rod_length

        unit_length = det_unit.rodlength * 0.997
        string_info = {
            "string_id": string_id,
            "string_rot": string_rot,
            "string_rot_m": string_rot_m,
            "string_meta": string_meta,
            "x_pos": x_pos,
            "y_pos": y_pos,
        }

        thicknesses = {
            "pen": 1.5,  # mm
            "cable": 0.076,  # mm
            "clamp": 1.8,  # mm
            "click": 1.5,  # mm flap thickness
            "insulator": 2.4,  # mm flap thickness
        }

        _place_hpge_unit(z_unit_bottom, det_unit, unit_length, string_info, thicknesses, b)




    # the copper rod is slightly longer after the last detector.
    copper_rod_length_from_z0 = total_rod_length + 3.5
    copper_rod_length = copper_rod_length_from_z0 + 12

    minishroud_length = MINISHROUD_LENGTH[0] + string_meta.get("minishroud_delta_length_in_mm", 0)
    assert total_rod_length < minishroud_length
    nms = _get_nylon_mini_shroud(
        string_meta.minishroud_radius_in_mm, minishroud_length, True, b.materials, b.registry
    )
    z_nms = z0_string - copper_rod_length_from_z0 + minishroud_length / 2 - MINISHROUD_END_THICKNESS
    nms_pv = geant4.PhysicalVolume(
        [0, 0, 0],
        [x_pos, y_pos, z_nms],
        nms,
        nms.name + "_string_" + string_id,
        b.mother_lv,
        b.registry,
    )
    _add_nms_surfaces(nms_pv, b.mother_pv, b.materials, b.registry)
    nms_top = _get_nylon_mini_shroud(
        string_meta.minishroud_radius_in_mm - MINISHROUD_END_THICKNESS,
        MINISHROUD_LENGTH[1],
        True,
        b.materials,
        b.registry,
        min_radius=10,
    )
    nms_pv = geant4.PhysicalVolume(
        [0, 0, 0],
        [x_pos, y_pos, z0_string + 15 + MINISHROUD_LENGTH[1] / 2],
        nms_top,
        nms_top.name + "_string_" + string_id,
        b.mother_lv,
        b.registry,
    )
    _add_nms_surfaces(nms_pv, b.mother_pv, b.materials, b.registry)

    support, tristar = _get_support_structure(string_slots[1].baseplate, b.materials, b.registry)
    geant4.PhysicalVolume(
        [0, 0, np.deg2rad(30) + string_rot],
        [x_pos, y_pos, z0_string + 12],  # this offset of 12 is measured from the CAD file.
        support,
        support.name + "_string_" + string_id,
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [0, 0, string_rot],
        [x_pos, y_pos, z0_string + 12 - 1e-6],  # this offset of 12 is measured from the CAD file.
        tristar,
        tristar.name + "_string_" + string_id,
        b.mother_lv,
        b.registry,
    )

    copper_rod_r = string_meta.rod_radius_in_mm
    assert copper_rod_r < string_meta.minishroud_radius_in_mm - 0.75
    copper_rod_name = f"string_{string_id}_cu_rod"
    # the rod has a radius of 1.5 mm, but this would overlap with the coarse model of the PPC top PEN ring.
    copper_rod = geant4.solid.Tubs(copper_rod_name, 0, 1.43, copper_rod_length, 0, 2 * math.pi, b.registry)
    copper_rod = geant4.LogicalVolume(copper_rod, b.materials.metal_copper, copper_rod_name, b.registry)
    copper_rod.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    for i in range(3):
        copper_rod_th = np.deg2rad(-30 - i * 120)
        delta = copper_rod_r * string_rot_m @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        geant4.PhysicalVolume(
            [0, 0, 0],
            [x_pos + delta[0], y_pos + delta[1], z0_string + 12 - copper_rod_length / 2],
            copper_rod,
            f"{copper_rod_name}_{i}",
            b.mother_lv,
            b.registry,
        )


def _place_empty_string(string_id: str, b: core.InstrumentationData):
    """Place an empty string (i.e. with no HPGe detectors), optionally with a counterweight."""
    string_meta = b.special_metadata.hpge_string[string_id]

    angle_in_rad = math.pi * string_meta.angle_in_deg / 180
    x_pos = string_meta.radius_in_mm * math.cos(angle_in_rad)
    y_pos = -string_meta.radius_in_mm * math.sin(angle_in_rad)
    # rotation angle for anything in the string.
    string_rot = -np.pi + angle_in_rad

    # offset the height of the string by the length of the string support rod.
    # TODO: this is also still a warm length.
    z0_string = b.top_plate_z_pos - 410.1  # from CAD model.

    if "string_support_structure_short" not in b.registry.logicalVolumeDict:
        support_lv = _read_model(
            "StringSupportStructure-short.stl",
            "string_support_structure_short",
            b.materials.metal_copper,
            b.registry,
        )
        support_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        support_lv = b.registry.logicalVolumeDict["string_support_structure_short"]

    geant4.PhysicalVolume(
        [0, 0, np.deg2rad(30) + string_rot],
        [x_pos, y_pos, z0_string],
        support_lv,
        support_lv.name + "_string_" + string_id,
        b.mother_lv,
        b.registry,
    )

    # add the optional steel counterweight to the empty string.
    string_content = string_meta.get("empty_string_content", [])
    if len(string_content) == 0:
        return
    if len(string_content) != 1 or string_content[0] not in ("counterweight", "counterweight_ttx"):
        msg = f"invalid empty string content {string_content}"
        raise ValueError(msg)
    has_counterweight = string_content[0] in ("counterweight", "counterweight_ttx")
    wrap_tetratex = has_counterweight and string_content[0] == "counterweight_ttx"

    if has_counterweight:
        counterweight_height = 513  # mm
        counterweight_name = "counterweight" + ("_wrapped" if wrap_tetratex else "")
        if counterweight_name not in b.registry.logicalVolumeDict:
            counterweight = geant4.solid.Tubs(
                counterweight_name, 0, 77 / 2, counterweight_height, 0, 2 * math.pi, b.registry, "mm"
            )
            counterweight = geant4.LogicalVolume(
                counterweight, b.materials.metal_steel, counterweight_name, b.registry
            )
            counterweight.pygeom_color_rgba = [1, 1, 1, 1] if wrap_tetratex else [0.5, 0.5, 0.5, 1]

        # account for the shorter hanger (compared to an active string), and the distance between copper
        # hanger and weight (the latter is estimated from photos).
        counterweight_z = z0_string + 130.5 - 30 - counterweight_height / 2
        counterweight_pv = geant4.PhysicalVolume(
            [0, 0, 0],
            [x_pos, y_pos, counterweight_z],
            b.registry.logicalVolumeDict[counterweight_name],
            f"{counterweight_name}_{string_id}",
            b.mother_lv,
            b.registry,
        )

        if wrap_tetratex:
            # note: no volume that actually has tetratex material, here. The surface alone should be fine
            # (propagation of light into the volume will not occur with this surface).
            geant4.BorderSurface(
                f"bsurface_lar_ttx_{string_id}",
                b.mother_pv,
                counterweight_pv,
                b.materials.surfaces.to_tetratex,
                b.registry,
            )


# Those dimensions are from an email from A. Lubashevskiy to L. Varriano on Dec 12, 2023; on the NMS made at
# TUM in May 2022.
MINISHROUD_THICKNESS = 0.125  # mm
MINISHROUD_END_THICKNESS = 2 * MINISHROUD_THICKNESS
MINISHROUD_LENGTH = (1000, 20)


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

    pen_lv_name = f"pen_{size}"
    if pen_lv_name not in registry.logicalVolumeDict:
        pen_file = f"BasePlate_{size}.stl" if size != "ppc_small" else "TopPlate_ppc.stl"
        pen_lv = _read_model(pen_file, pen_lv_name, materials.pen, registry)
        pen_lv.pygeom_color_rgba = colors[size]

    return registry.logicalVolumeDict[pen_lv_name]


def _get_support_structure(
    size: str,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> tuple[geant4.LogicalVolume, geant4.LogicalVolume]:
    """Get the (simplified) support structure and the tristar of the requested size.

    .. note :: Both models' coordinate origins are a the top face of the tristar structure."""
    if "string_support_structure" not in registry.logicalVolumeDict:
        support_lv = _read_model(
            "StringSupportStructure.stl", "string_support_structure", materials.metal_copper, registry
        )
        support_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        support_lv = registry.logicalVolumeDict["string_support_structure"]

    tristar_lv_name = f"tristar_{size}"
    if tristar_lv_name not in registry.logicalVolumeDict:
        tristar_lv = _read_model(f"TriStar_{size}.stl", f"tristar_{size}", materials.pen, registry)
        tristar_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        tristar_lv = registry.logicalVolumeDict[tristar_lv_name]

    return support_lv, tristar_lv


def _get_nylon_mini_shroud(
    radius: int,
    length: int,
    top_open: bool,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
    min_radius: int = 0,
) -> geant4.LogicalVolume:
    """Create a nylon/TPB funnel of the given outer dimensions, which will be closed at the bottom.

    .. note:: this can also be used for calibration tubes.
    """
    assert top_open  # just for b/c of this shared interface. remove in future.
    shroud_name = f"minishroud_{radius}x{length}"
    if shroud_name not in registry.logicalVolumeDict:
        outer = geant4.solid.Tubs(f"{shroud_name}_outer", min_radius, radius, length, 0, 2 * np.pi, registry)
        inner = geant4.solid.Tubs(
            f"{shroud_name}_inner",
            0,
            radius - MINISHROUD_THICKNESS,
            # at the top/bottom, the NMS has essentially two layers.
            length - (0 if top_open else 2 * MINISHROUD_END_THICKNESS),
            0,
            2 * np.pi,
            registry,
        )
        # subtract the slightly smaller solid from the larger one, to get a hollow and closed volume.
        inner_z = (1 if top_open else 0) * MINISHROUD_END_THICKNESS
        shroud = geant4.solid.Subtraction(shroud_name, outer, inner, [[0, 0, 0], [0, 0, inner_z]], registry)
        nms_lv = geant4.LogicalVolume(shroud, materials.tpb_on_nylon, shroud_name, registry)
        nms_lv.pygeom_color_rgba = (1, 0.86, 0.86, 0.2)

    return registry.logicalVolumeDict[shroud_name]


def _add_pen_surfaces(
    pen_pv: geant4.PhysicalVolume,
    mother_pv: geant4.LogicalVolume,
    mats: materials.OpticalMaterialRegistry,
    reg: geant4.Registry,
):
    # between LAr and PEN we need a surface in both directions.
    geant4.BorderSurface("bsurface_lar_pen_" + pen_pv.name, mother_pv, pen_pv, mats.surfaces.lar_to_pen, reg)
    geant4.BorderSurface("bsurface_tpb_pen_" + pen_pv.name, pen_pv, mother_pv, mats.surfaces.lar_to_pen, reg)


def _add_nms_surfaces(
    nms_pv: geant4.PhysicalVolume,
    mother_pv: geant4.LogicalVolume,
    mats: materials.OpticalMaterialRegistry,
    reg: geant4.Registry,
):
    # between LAr and the NMS we need a surface in both directions.
    geant4.BorderSurface("bsurface_lar_nms_" + nms_pv.name, mother_pv, nms_pv, mats.surfaces.lar_to_tpb, reg)
    geant4.BorderSurface("bsurface_nms_lar_" + nms_pv.name, nms_pv, mother_pv, mats.surfaces.lar_to_tpb, reg)


def _read_model(
    file: str, name: str, material: geant4.Material, registry: geant4.Registry
) -> geant4.LogicalVolume:
    file = resources.files("l200geom") / "models" / file
    solid = pyg4ometry.stl.Reader(file, solidname=name, centre=False, registry=registry).getSolid()
    return geant4.LogicalVolume(solid, material, name, registry)


def _get_hv_cable(
    name: str,
    cable_thickness: float,
    clamp_thickness: float,
    cable_length: float,
    materials: materials.OpticalMaterialRegistry,
    mother_pv: geant4.LogicalVolume,
    reg: geant4.Registry,
):
    safety_margin = 1  # mm
    cable_length -= safety_margin

    hv_cable_under_clamp = geant4.solid.Box(
        name + "_hv_cable_under_clamp",
        8,
        13,
        cable_thickness,
        reg,
        "mm",
    )
    hv_cable_clamp_to_curve = geant4.solid.Box(
        name + "_hv_cable_clamp_to_curve",
        5.5,
        2,
        cable_thickness,
        reg,
        "mm",
    )

    hv_cable_curve = geant4.solid.Tubs(
        name + "_hv_cable_curve", 3.08, 3.08 + cable_thickness, 2.0, 0, math.pi / 2.0, reg, "mm"
    )

    hv_cable_along_unit = geant4.solid.Box(
        name + "_hv_along_unit",
        cable_thickness,
        2.0,
        cable_length,
        reg,
        "mm",
    )

    hv_cable_part1 = geant4.solid.Union(
        name + "_hv_cable_part1",
        hv_cable_under_clamp,
        hv_cable_clamp_to_curve,
        [[0, 0, 0], [8 / 2.0 + 5.5 / 2.0, 0, 0]],
        reg,
    )

    hv_cable_part2 = geant4.solid.Union(
        name + "_hv_cable_part2",
        hv_cable_part1,
        hv_cable_curve,
        [[-np.pi / 2, 0, 0], [8 / 2.0 + 5.5, 0, 3.08 + cable_thickness / 2.0]],
        reg,
    )

    hv_cable = geant4.solid.Union(
        name + "_hv_cable",
        hv_cable_part2,
        hv_cable_along_unit,
        [[0, 0, 0], [8 / 2.0 + 5.5 + 3.08 + cable_thickness / 2.0, 0, 3.08 + cable_length / 2.0]],
        reg,
    )

    hv_clamp = geant4.solid.Box(
        name + "_hv_clamp",
        5,
        13,
        clamp_thickness,
        reg,
        "mm",
    )

    hv_cable_lv = geant4.LogicalVolume(
        hv_cable,
        materials.metal_copper,
        name + "_hv_cable",
        reg,
    )

    hv_clamp_lv = geant4.LogicalVolume(
        hv_clamp,
        materials.ultem,
        name + "_hv_clamp",
        reg,
    )

    return hv_cable_lv, hv_clamp_lv


def _get_signal_cable_and_asic(
    name: str,
    cable_thickness: float,
    clamp_thickness: float,
    cable_length: float,
    materials: materials.OpticalMaterialRegistry,
    mother_pv: geant4.LogicalVolume,
    reg: geant4.Registry,
):
    safety_margin = 1  # mm
    cable_length -= safety_margin

    signal_cable_under_clamp = geant4.solid.Box(
        name + "_signal_cable_under_clamp",
        16,
        13,
        cable_thickness,
        reg,
        "mm",
    )
    signal_cable_clamp_to_curve = geant4.solid.Box(
        name + "_signal_cable_clamp_to_curve",
        23.25,
        2,
        cable_thickness,
        reg,
        "mm",
    )
    signal_cable_curve = geant4.solid.Tubs(
        name + "_signal_cable_curve", 3.08, 3.08 + cable_thickness, 2.0, 0, math.pi / 2.0, reg, "mm"
    )
    signal_cable_along_unit = geant4.solid.Box(
        name + "_signal_along_unit",
        cable_thickness,
        2.0,
        cable_length,
        reg,
        "mm",
    )
    signal_cable_part1 = geant4.solid.Union(
        name + "_signal_cable_part1",
        signal_cable_under_clamp,
        signal_cable_clamp_to_curve,
        [[0, 0, 0], [16 / 2.0 + 23.25 / 2.0, 0, 0]],
        reg,
    )
    signal_cable_part2 = geant4.solid.Union(
        name + "_signal_cable_part2",
        signal_cable_part1,
        signal_cable_curve,
        [[np.pi / 2, 0, 0], [16 / 2.0 + 23.25, 0, -3.08 - cable_thickness / 2.0]],
        reg,
    )
    signal_cable = geant4.solid.Union(
        name + "_signal_cable",
        signal_cable_part2,
        signal_cable_along_unit,
        [[0, 0, 0], [16 / 2.0 + 23.25 + 3.08 + cable_thickness / 2.0, 0, -3.08 - cable_length / 2.0]],
        reg,
    )

    signal_clamp_part1 = geant4.solid.Box(
        name + "_signal_clamp_part1",
        5,
        13,
        clamp_thickness,
        reg,
        "mm",
    )
    signal_clamp_part2 = geant4.solid.Box(
        name + "_signal_clamp_part2",
        9,
        2.5,
        clamp_thickness,
        reg,
        "mm",
    )
    signal_clamp_part3 = geant4.solid.Union(
        name + "_signal_clamp_part3",
        signal_clamp_part1,
        signal_clamp_part2,
        [[0, 0, 0], [5 / 2.0 + 9 / 2.0, 13 / 2.0 - 2.5 / 2.0, 0]],
        reg,
    )
    signal_clamp = geant4.solid.Union(
        name + "_signal_clamp",
        signal_clamp_part3,
        signal_clamp_part2,
        [[0, 0, 0], [5 / 2.0 + 9 / 2.0, -13 / 2.0 + 2.5 / 2.0, 0]],
        reg,
    )

    signal_asic = geant4.solid.Box(
        name + "_signal_asic",
        1,
        1,
        0.5,
        reg,
        "mm",
    )

    signal_cable_lv = geant4.LogicalVolume(
        signal_cable,
        materials.metal_copper,
        name + "_signal_cable",
        reg,
    )

    signal_clamp_lv = geant4.LogicalVolume(
        signal_clamp,
        materials.ultem,
        name + "_signal_clamp",
        reg,
    )

    signal_asic_lv = geant4.LogicalVolume(
        signal_asic,
        materials.silica,
        name + "_signal_asic",
        reg,
    )

    return signal_cable_lv, signal_clamp_lv, signal_asic_lv


def _get_click_and_insulator(
    det_unit: HPGeDetUnit,
    click_top_flap_thickness: float,
    insulator_du_holder_flap_thickness: float,
    insulator_top_length: float,
    materials: materials.OpticalMaterialRegistry,
    reg: geant4.Registry,
):
    safety_margin = 0.1
    click_top_flap = geant4.solid.Box(
        det_unit.name + "_click_top_flap",
        20.8,
        5,
        click_top_flap_thickness,
        reg,
        "mm",
    )

    click_top_clamp = geant4.solid.Box(
        det_unit.name + "_click_top_clamp",
        7.8,
        5,
        2.2,
        reg,
        "mm",
    )

    # Union the flap and clamp
    click_top_without_hole = geant4.solid.Union(
        det_unit.name + "_click_top_without_hole",
        click_top_flap,
        click_top_clamp,
        [[0, 0, 0], [20.8 / 2.0 - 7.8 / 2.0, 0, -2.2 / 2.0 - click_top_flap_thickness / 2.0]],
        reg,
    )

    click_top_carving_hole = geant4.solid.Tubs(
        det_unit.name + "_click_top_carving_hole",
        0,
        1.5 + safety_margin,
        2 * (click_top_flap_thickness + 2.2),
        0,
        math.pi * 2,
        reg,
        "mm",
    )

    # Perform subtraction only once
    click_top = geant4.solid.Subtraction(
        det_unit.name + "_click_top",
        click_top_without_hole,
        click_top_carving_hole,
        [[0, 0, 0], [5.60, 0, 0]],  # Adjust the position of the hole as needed
        reg,
    )

    insulator_du_holder_flap = geant4.solid.Box(
        det_unit.name + "_insulator_du_holder_flap",
        16.5,
        7,
        insulator_du_holder_flap_thickness,
        reg,
        "mm",
    )

    safety_margin_touching_detector = 0.25

    insulator_du_holder_clamp = geant4.solid.Box(
        det_unit.name + "_insulator_du_holder_clamp",
        insulator_top_length - safety_margin_touching_detector,
        7,
        5.5 - insulator_du_holder_flap_thickness,
        reg,
        "mm",
    )

    # Union the flap and clamp
    insulator_du_holder_without_hole = geant4.solid.Union(
        det_unit.name + "_insulator_du_holder_without_hole",
        insulator_du_holder_flap,
        insulator_du_holder_clamp,
        [
            [0, 0, 0],
            [
                16.5 / 2.0 - (insulator_top_length - safety_margin_touching_detector) / 2.0,
                0,
                (5.5 - insulator_du_holder_flap_thickness) / 2.0 + insulator_du_holder_flap_thickness / 2.0,
            ],
        ],
        reg,
    )

    insulator_du_holder_carving_hole = geant4.solid.Tubs(
        det_unit.name + "_insulator_du_holder_carving_hole",
        0,
        1.5 + safety_margin,
        3 * 5.5,
        0,
        math.pi * 2,
        reg,
        "mm",
    )

    # Perform subtraction only once
    insulator_du_holder = geant4.solid.Subtraction(
        det_unit.name + "_insulator_du_holder",
        insulator_du_holder_without_hole,
        insulator_du_holder_carving_hole,
        [[0, 0, 0], [16.5 / 2.0 - 1.5, 0, 0]],  # Adjust the position of the hole as needed
        reg,
    )

    click_top_lv = geant4.LogicalVolume(
        click_top,
        materials.ultem,
        det_unit.name + "_click_top",
        reg,
    )

    insulator_du_holder_lv = geant4.LogicalVolume(
        insulator_du_holder,
        materials.ultem,
        det_unit.name + "_insulator_du_holder",
        reg,
    )

    return click_top_lv, insulator_du_holder_lv
