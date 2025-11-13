from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np
from dbetto import AttrsDict, TextDB
from legendhpges import make_hpge
from pyg4ometry import geant4
from pygeomtools import RemageDetectorInfo
from scipy.spatial.transform import Rotation

from . import core, materials
from .utils import _read_model

log = logging.getLogger(__name__)


def place_hpge_strings(hpge_metadata: TextDB, b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 HPGe strings."""
    # derive the strings from the channelmap.
    mass_total = 0
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
        if hpge_meta.production.enrichment.val is None:
            log.warning("%s has no enrichment in metadata - setting to dummy value 0.9!", hpge_meta.name)
            hpge_meta.production.enrichment = 0.9

        hpge_string_id = str(ch_meta.location.string)
        hpge_unit_id_in_string = ch_meta.location.position

        if hpge_string_id not in strings_to_build:
            strings_to_build[hpge_string_id] = {}

        hpge_extra_meta = b.special_metadata.hpges[hpge_meta.name]

        log.debug("Building %s", hpge_meta.name)

        hpge = make_hpge(full_meta, b.registry)
        strings_to_build[hpge_string_id][hpge_unit_id_in_string] = HPGeDetUnit(
            hpge_meta.name,
            hpge_meta.production.manufacturer,
            ch_meta.daq.rawid,
            hpge,
            hpge_meta.geometry.height_in_mm,
            hpge_meta.geometry.radius_in_mm,
            hpge_extra_meta["baseplate"],
            hpge_extra_meta["rodlength_in_mm"],
            # convert the "warm" length of the rod to the (shorter) length in the cooled down state.
            hpge_extra_meta["rodlength_in_mm"] * 0.997,
            hpge_extra_meta.get("has_top_insulators", False),
            full_meta,
        )
        mass_total += hpge.mass.to("kg").m

    log.info("Total HPGe mass %.2f kg", mass_total)

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
    rodlength_cold: float
    has_top_insulators: bool
    meta: AttrsDict


def _place_front_end_and_insulators(
    det_unit: HPGeDetUnit,
    string_info: dict,
    b: core.InstrumentationData,
    z_pos: dict,
    thickness: dict,
    parts_origin: dict,
):
    string_rot_v = np.array([np.sin(string_info.rot), np.cos(string_info.rot)])
    string_pos_v = np.array([string_info.x, string_info.y])

    cu_pin = _get_cu_pin(thickness["cu_pin"], b)
    phbr_washer = _get_phbr_washer(thickness["washer"], b)
    phbr_spring = _get_phbr_spring(b)

    # add cable and clamp
    signal_cable = _get_signal_cable(thickness["cable"], det_unit.rodlength_cold, b)
    signal_clamp, signal_lmfe, signal_holes = _get_signal_clamp_and_lmfe(thickness["clamp"], b)
    signal_cable.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    signal_clamp.pygeom_color_rgba = (0.64, 0.54, 0.31, 1)
    signal_lmfe.pygeom_color_rgba = (0.64, 0.54, 0.31, 0.5)

    angle_signal = math.pi * 1 / 2.0 - string_info.rot
    x_clamp, y_clamp = string_pos_v + parts_origin["signal"] * string_rot_v
    x_cable, y_cable = string_pos_v + (parts_origin["signal"] + 7.5 / 2 + 0.1) * string_rot_v
    x_spring, y_spring = string_pos_v + (parts_origin["signal"] - (13 - 7.5) / 2) * string_rot_v
    lmfe_origin = parts_origin["signal"] + (7.5 + 16 + 0.1) / 2
    x_lmfe, y_lmfe = string_pos_v + lmfe_origin * string_rot_v

    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_cable, y_cable, z_pos["clamp"]],
        signal_cable,
        f"{signal_cable.name}_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_clamp, y_clamp, z_pos["clamp"]],
        signal_clamp,
        f"{signal_clamp.name}_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_lmfe, y_lmfe, z_pos["clamp"]],
        signal_lmfe,
        f"{signal_lmfe.name}_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [math.pi, 0, angle_signal],
        [x_spring, y_spring, z_pos["clamp"]],
        phbr_spring,
        f"{phbr_spring.name}_signal_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )

    def place_clamp_details(r0, name, holes, z_clamp, *, sign=1, z_sign=1, rot_v=string_rot_v):
        rot_v_perp = np.cross([*rot_v, 0], [0, 0, 1])[0:2]
        for hole_idx, hole in enumerate(holes):
            x_hole, y_hole = string_pos_v + sign * (r0 + hole[0]) * rot_v + hole[1] * rot_v_perp
            pin_outside = thickness["cu_pin"] - thickness["pen"] - thickness["clamp"]
            geant4.PhysicalVolume(
                [0, 0, 0],
                [x_hole, y_hole, z_clamp + z_sign * thickness["pen"] / 2 - z_sign * pin_outside / 2],
                cu_pin,
                f"{cu_pin.name}_{name}_{det_unit.name}_{hole_idx}",
                b.mother_lv,
                b.registry,
            )
            delta_z_washer = (thickness["clamp"] + thickness["washer"]) / 2 + thickness["safety"]
            geant4.PhysicalVolume(
                [0, 0, 0],
                [x_hole, y_hole, z_clamp - z_sign * delta_z_washer],
                phbr_washer,
                f"{phbr_washer.name}_{name}_{det_unit.name}_{hole_idx}",
                b.mother_lv,
                b.registry,
            )

    place_clamp_details(parts_origin["signal"], "signal", signal_holes, z_pos["clamp"])

    # shorter HV cable for top contact on PPCs.
    hv_cable_length = det_unit.rodlength_cold if not det_unit.name.startswith("P") else 15
    hv_cable = _get_hv_cable(thickness["cable"], hv_cable_length, b)
    hv_clamp, hv_holes = _get_hv_clamp(thickness["clamp"], b)
    hv_cable.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    hv_clamp.pygeom_color_rgba = (0.64, 0.54, 0.31, 1)

    angle_hv = math.pi / 2 + string_info.rot
    hv_rot_v = string_rot_v
    hv_washer_z_sign = 1
    if det_unit.name.startswith("P"):
        hv_rot_offset = -math.pi * 1 / 3
        angle_hv += hv_rot_offset
        hv_rot_v = np.array(
            [
                np.sin(string_info.rot + hv_rot_offset),
                np.cos(string_info.rot + hv_rot_offset),
            ]
        )
        hv_washer_z_sign = -1
    hv_z_pos = z_pos["clamp" if not det_unit.name.startswith("P") else "clamp_top"]

    x_clamp, y_clamp = string_pos_v - parts_origin["hv"] * hv_rot_v
    x_cable, y_cable = string_pos_v - (parts_origin["hv"] - 3 + 2 * 2e-9) * hv_rot_v
    x_spring, y_spring = string_pos_v - (parts_origin["hv"] - 3 + 2e-9) * hv_rot_v

    geant4.PhysicalVolume(
        [0, 0, angle_hv],
        [x_cable, y_cable, hv_z_pos],
        hv_cable,
        f"{hv_cable.name}_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [0, 0, angle_hv],
        [x_clamp, y_clamp, hv_z_pos],
        hv_clamp,
        f"{hv_clamp.name}_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )
    geant4.PhysicalVolume(
        [0, 0, angle_hv],
        [x_spring, y_spring, hv_z_pos],
        phbr_spring,
        f"{phbr_spring.name}_hv_{det_unit.name}",
        b.mother_lv,
        b.registry,
    )

    place_clamp_details(
        parts_origin["hv"], "hv", hv_holes, hv_z_pos, sign=-1, z_sign=hv_washer_z_sign, rot_v=hv_rot_v
    )

    # this is a heuristic to get the value of "dimension A"; in reality this is a step function.
    insulator_top_length = string_info.meta.rod_radius_in_mm - det_unit.radius + 1.5

    weldment = _get_weldment(thickness["weldment"], b)
    insulator = _get_insulator(thickness["insulator"], insulator_top_length, b)

    for i in range(3):
        copper_rod_th = np.deg2rad(-30 - i * 120)
        pieces_th = string_info.rot + np.deg2rad(-(i + 1) * 120)
        delta_weldment = (
            (string_info.meta.rod_radius_in_mm - 7)
            * string_info.rot_m
            @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        )
        delta_insulator = (
            (string_info.meta.rod_radius_in_mm - (16.5 / 2.0 - 1.5))
            * string_info.rot_m
            @ np.array([np.cos(copper_rod_th), np.sin(copper_rod_th)])
        )
        geant4.PhysicalVolume(
            [0, 0, pieces_th],
            [
                string_info.x + delta_weldment[0],
                string_info.y + delta_weldment[1],
                z_pos["weldment"],
            ],
            weldment,
            f"{weldment.name}_{det_unit.name}_{i}",
            b.mother_lv,
            b.registry,
        )
        geant4.PhysicalVolume(
            [0, 0, pieces_th],
            [
                string_info.x + delta_insulator[0],
                string_info.y + delta_insulator[1],
                z_pos["insulator"],
            ],
            insulator,
            f"{insulator.name}_{det_unit.name}_{i}",
            b.mother_lv,
            b.registry,
        )
        if det_unit.has_top_insulators:
            assert det_unit.name.startswith("V")
            insulator_top_rot = Rotation.from_euler("XZ", [-np.pi, pieces_th]).as_euler("xyz")
            geant4.PhysicalVolume(
                list(insulator_top_rot),
                [
                    string_info.x + delta_insulator[0],
                    string_info.y + delta_insulator[1],
                    z_pos["insulator_top"],
                ],
                insulator,
                f"{insulator.name}_{det_unit.name}_{i}_top",
                b.mother_lv,
                b.registry,
            )


def _hpge_unit_get_z(bottom: float, det_unit: HPGeDetUnit) -> tuple[dict, dict]:
    t = {
        "pen": 1.5,  # mm
        "cable": 0.076,  # mm
        "clamp": 3.7,  # mm, but no constant thickness (HV +0.5 mm)
        "weldment": 1.5,  # mm flap thickness
        "insulator": 2.4,  # mm flap thickness
        "washer": 0.3,  # PhBr washers
        "cu_pin": 7.8,
    }

    # TODO: the large PEN plate model has an unexpected rib at one of the spots for the insulator,
    # so the real spacing cannot be used yet in the model.
    # the mirrored medium_ortec plate would also require a special treatment.
    insulator_spacer = t["insulator"]  # - 0.3

    # - notes for those comparing this to MaGe (those offsets are not from there, but from the
    #   CAD model): the detector unit (DU)-local z coordinates are inverted in comparison to
    #   the coordinates here, as well as to the string coordinates in MaGe.
    # - In MaGe, the end of the three support rods is at +11.1 mm, the PEN plate at +4 mm, the
    #   diode at -diodeHeight/2-0.025 mm, so that the crystal contact is at DU-z 0 mm.
    z_det = bottom + 3.7 + 1.3 + t["pen"] + insulator_spacer

    safety = 0.001  # 1 micro meter
    t["safety"] = safety
    # - note from CAD model: the distance between PEN plate top and detector bottom face varies
    #   a lot between different diodes (i.e. BEGe's/IC's all(?) use a single standard insulator
    #   type, and have a distance of 2.1 mm; for PPCs this varies between ca. 2.5 and 4 mm.)
    z = {
        "det": z_det,
        "insulator": z_det - t["insulator"] / 2.0 - safety,
        "pen": z_det - insulator_spacer - t["pen"] / 2.0 - safety * 2,
        "weldment": z_det - insulator_spacer - t["pen"] - t["weldment"] / 2.0 - safety * 3,
        "clamp": z_det - insulator_spacer - t["pen"] - t["clamp"] / 2.0 - safety * 4,
        "pen_top": z_det + det_unit.height + t["pen"] / 2 + safety,
        "insulator_top": z_det + det_unit.height + t["insulator"] / 2.0 + safety,
        "clamp_top": z_det + det_unit.height + t["pen"] + t["clamp"] / 2.0 + safety * 3,
    }
    return t, z


def _place_hpge_unit(
    z_unit_bottom: float,
    det_unit: HPGeDetUnit,
    string_info: AttrsDict,
    b: core.InstrumentationData,
):
    thicknesses, z_pos = _hpge_unit_get_z(z_unit_bottom, det_unit)

    det_pv = geant4.PhysicalVolume(
        [0, 0, 0],
        [string_info.x, string_info.y, z_pos["det"]],
        det_unit.lv,
        det_unit.name,
        b.mother_lv,
        b.registry,
    )
    det_pv.pygeom_active_detector = RemageDetectorInfo("germanium", det_unit.rawid, det_unit.meta)
    det_unit.lv.pygeom_color_rgba = (0.5, 0.5, 0.5, 1)

    # add germanium reflective surface.
    geant4.BorderSurface(
        "bsurface_lar_ge_" + det_pv.name,
        b.mother_pv,
        det_pv,
        b.materials.surfaces.to_germanium,
        b.registry,
    )

    baseplate = det_unit.baseplate
    pen_rot = [0, 0, string_info.rot]
    # a lot of Ortec detectors have modified medium plates.
    if det_unit.name.startswith("V") and det_unit.baseplate == "medium" and det_unit.manufacturer == "Ortec":
        # TODO: what is with "V01389A"?
        baseplate = "medium_ortec"
        # This rotation is not physical, but gets us closer to the real model of the PEN plates.
        # In the CAD model, most plates are mirrored, compared to reality (some are also correct in the
        # first place), i.e. how the plates in PGT were produced. So the STL mesh is also mirrored, so
        # flip it over.
        # note/TODO: this rotation should be replaced by a correct mesh, so that the counterbores are
        # on the correct side. This might be necessary to fit in other parts!
        pen_rot = Rotation.from_euler("XZ", [-np.pi, string_info.rot]).as_euler("xyz")
    pen_plate = _get_pen_plate(baseplate, b)

    if pen_plate is not None:
        pen_pv = geant4.PhysicalVolume(
            list(pen_rot),
            [string_info.x, string_info.y, z_pos["pen"]],
            pen_plate,
            "pen_" + det_unit.name,
            b.mother_lv,
            b.registry,
        )
        _add_pen_surfaces(pen_pv, b.mother_pv, b.materials, b.registry)

    # (Majorana) PPC detectors have a top PEN ring.
    if det_unit.name.startswith("P"):
        assert det_unit.baseplate == "small"
        pen_plate = _get_pen_plate("ppc_small", b)
        if pen_plate is not None:
            # this is a physical rotation (i.e. the model should be right, it is just rotated in the file).
            pen_top_rot = Rotation.from_euler("XZ", [-np.pi, string_info.rot + np.pi * 2 / 3]).as_euler("xyz")
            pen_pv = geant4.PhysicalVolume(
                list(pen_top_rot),
                [string_info.x, string_info.y, z_pos["pen_top"]],
                pen_plate,
                "pen_top_" + det_unit.name,
                b.mother_lv,
                b.registry,
            )
            _add_pen_surfaces(pen_pv, b.mother_pv, b.materials, b.registry)

    # positions from center of detector to center of volume center
    # note for finding these values: the holes of the HV receptacle have 3 mm distance to the side,
    # 4 mm for signal receptacle.
    if det_unit.baseplate == "small":
        fe_ins_origins = {"signal": 8.8, "hv": 32.35}
    elif det_unit.baseplate in {"medium", "large"}:
        fe_ins_origins = {"signal": 14.25, "hv": 37.5}
    elif det_unit.baseplate == "xlarge":
        fe_ins_origins = {"signal": 17.3, "hv": 40.35}

    _place_front_end_and_insulators(det_unit, string_info, b, z_pos, thicknesses, fe_ins_origins)


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
    string_info = AttrsDict(
        {"rot": string_rot, "rot_m": string_rot_m, "meta": string_meta, "x": x_pos, "y": y_pos}
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

        total_rod_length += det_unit.rodlength_cold
        z_unit_bottom = z0_string - total_rod_length

        _place_hpge_unit(z_unit_bottom, det_unit, string_info, b)

    # the copper rod is slightly longer after the last detector (estimate from CAD model, probably does
    # not match reality).
    copper_rod_length_from_z0 = total_rod_length + 3.5
    copper_rod_length = copper_rod_length_from_z0 + 12
    copper_rod_r = string_meta.rod_radius_in_mm

    if string_meta.minishroud_radius_in_mm is not None:
        minishroud_length = MINISHROUD_LENGTH[0] + string_meta.get("minishroud_delta_length_in_mm", 0)
        assert total_rod_length < minishroud_length
        assert copper_rod_r < string_meta.minishroud_radius_in_mm - 0.75
        nms = _get_nylon_mini_shroud(
            string_meta.minishroud_radius_in_mm,
            minishroud_length,
            True,
            b.materials,
            b.registry,
            minishroud_name="minishroud_side",
        )
        z_nms = z0_string - copper_rod_length_from_z0 + minishroud_length / 2 - MINISHROUD_END_THICKNESS - 0.1
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
            minishroud_name="minishroud_top",
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

    support, tristar = _get_support_structure(string_slots[1].baseplate, b)
    if support is not None:
        geant4.PhysicalVolume(
            [0, 0, np.deg2rad(30) + string_rot],
            [x_pos, y_pos, z0_string + 12],  # this offset of 12 is measured from the CAD file.
            support,
            support.name + "_string_" + string_id,
            b.mother_lv,
            b.registry,
        )
    if tristar is not None:
        geant4.PhysicalVolume(
            [0, 0, string_rot],
            [x_pos, y_pos, z0_string + 12 - 1e-6],  # this offset of 12 is measured from the CAD file.
            tristar,
            tristar.name + "_string_" + string_id,
            b.mother_lv,
            b.registry,
        )

    copper_rod_name = f"hpge_support_copper_string_{string_id}_cu_rod"
    # the rod has a radius of 1.5 mm, but this would overlap with the coarse model of the PPC top PEN ring.
    copper_rod = geant4.solid.Tubs(copper_rod_name, 0, 1.40, copper_rod_length, 0, 2 * math.pi, b.registry)
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

    if "hpge_support_copper_string_support_structure_short" not in b.registry.logicalVolumeDict:
        support_lv = _read_model(
            "StringSupportStructure-short.stl",
            "hpge_support_copper_string_support_structure_short",
            b.materials.metal_copper,
            b,
        )
        if support_lv is not None:
            support_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        support_lv = b.registry.logicalVolumeDict["hpge_support_copper_string_support_structure_short"]

    if support_lv is not None:
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
    b: core.InstrumentationData,
) -> geant4.LogicalVolume:
    if size not in ["small", "medium", "medium_ortec", "large", "xlarge", "ppc_small"]:
        msg = f"Invalid PEN-plate size {size}"
        raise ValueError(msg)

    pen_lv_name = f"pen_{size}"
    if pen_lv_name not in b.registry.logicalVolumeDict:
        pen_file = f"BasePlate_{size}.stl" if size != "ppc_small" else "TopPlate_ppc.stl"
        pen_lv = _read_model(pen_file, pen_lv_name, b.materials.pen, b)
        if pen_lv is not None:
            pen_lv.pygeom_color_rgba = (1, 1, 1, 0.3)

    return b.registry.logicalVolumeDict.get(pen_lv_name)


def _get_support_structure(
    size: str,
    b: core.InstrumentationData,
) -> tuple[geant4.LogicalVolume, geant4.LogicalVolume]:
    """Get the (simplified) support structure and the tristar of the requested size.

    .. note :: Both models' coordinate origins are a the top face of the tristar structure."""
    if "hpge_support_copper_string_support_structure" not in b.registry.logicalVolumeDict:
        support_lv = _read_model(
            "StringSupportStructure.stl",
            "hpge_support_copper_string_support_structure",
            b.materials.metal_copper,
            b,
        )
        if support_lv is not None:
            support_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        support_lv = b.registry.logicalVolumeDict["hpge_support_copper_string_support_structure"]

    tristar_lv_name = f"hpge_support_copper_tristar_{size}"
    if tristar_lv_name not in b.registry.logicalVolumeDict:
        tristar_lv = _read_model(f"TriStar_{size}.stl", tristar_lv_name, b.materials.metal_copper, b)
        if tristar_lv is not None:
            tristar_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    else:
        tristar_lv = b.registry.logicalVolumeDict[tristar_lv_name]

    return support_lv, tristar_lv


def _get_nylon_mini_shroud(
    radius: int,
    length: int,
    top_open: bool,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
    min_radius: int = 0,
    minishroud_name: str = "minishroud",
) -> geant4.LogicalVolume:
    """Create a nylon/TPB funnel of the given outer dimensions, which will be closed at the bottom.

    .. note:: this can also be used for calibration tubes.
    """
    assert top_open  # just for b/c of this shared interface. remove in future.
    shroud_name = f"{minishroud_name}_{radius}x{length}"
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
        nms_lv.pygeom_color_rgba = (0.55, 0.79, 0.97, 0.1)

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


def _get_hv_cable(
    cable_thickness: float,
    cable_length: float,
    b: core.InstrumentationData,
):
    cable_name = f"cable_hv_{cable_length}"
    if cable_name in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict[cable_name]

    safety_margin = 1  # mm
    cable_length -= safety_margin

    hv_cable_radius = 3.08
    hv_cable_curve = geant4.solid.Tubs(
        f"{cable_name}_curve",
        hv_cable_radius,
        hv_cable_radius + cable_thickness,
        2.0,
        0,
        math.pi / 2.0,
        b.registry,
        "mm",
    )

    hv_cable_along_unit = geant4.solid.Box(
        f"{cable_name}_along_unit",
        cable_thickness,
        2.0,
        cable_length,
        b.registry,
        "mm",
    )

    hv_cable = geant4.solid.MultiUnion(
        cable_name,
        [hv_cable_curve, hv_cable_along_unit],
        [
            [[-np.pi / 2, 0, 0], [8 / 2.0 + 5.5, 0, hv_cable_radius + cable_thickness / 2.0]],
            [
                [0, 0, 0],
                [
                    8 / 2.0 + 5.5 + hv_cable_radius + cable_thickness / 2.0,
                    0,
                    hv_cable_radius + cable_length / 2.0,
                ],
            ],
        ],
        b.registry,
    )

    return geant4.LogicalVolume(
        hv_cable,
        b.materials.metal_copper,
        cable_name,
        b.registry,
    )


def _get_hv_clamp(clamp_thickness: float, b: core.InstrumentationData):
    holes = [[-3.45, 4, 0], [-3.45, -4, 0]]
    if "ultem_clamp_hv" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["ultem_clamp_hv"], holes

    hv_clamp_bulk = geant4.solid.Box("ultem_clamp_hv_bulk", 13, 13, clamp_thickness, b.registry, "mm")

    clamp_hole = geant4.solid.Tubs(
        "ultem_clamp_hv_hole", 0, 1.5, clamp_thickness, 0, 2 * math.pi, b.registry, "mm"
    )
    clamp_springhole = geant4.solid.Box("ultem_clamp_signal_hole", 13.1, 0.901, 0.251, b.registry, "mm")
    clamp_holes = geant4.solid.MultiUnion(
        "ultem_clamp_hv_holes",
        [clamp_hole, clamp_hole, clamp_springhole],
        [[[0, 0, 0], holes[0]], [[0, 0, 0], holes[1]], [[0, 0, 0], [0, 0, 0]]],
        b.registry,
    )
    hv_clamp = geant4.solid.Subtraction(
        "ultem_clamp_hv", hv_clamp_bulk, clamp_holes, [[0, 0, 0], [0, 0, 0]], b.registry
    )
    hv_clamp = geant4.LogicalVolume(hv_clamp, b.materials.ultem, "ultem_clamp_hv", b.registry)
    return hv_clamp, holes


def _get_signal_cable(
    cable_thickness: float,
    cable_length: float,
    b: core.InstrumentationData,
):
    cable_name = f"cable_signal_{cable_length}"
    if cable_name in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict[cable_name]

    safety_margin = 1  # mm
    cable_length -= safety_margin

    signal_cable_radius = 3.08
    signal_cable_clamp_to_curve = geant4.solid.Box(
        f"{cable_name}_clamp_to_curve",
        23.25 / 3,
        2,
        cable_thickness,
        b.registry,
        "mm",
    )
    signal_cable_curve = geant4.solid.Tubs(
        f"c{cable_name}_curve",
        signal_cable_radius,
        signal_cable_radius + cable_thickness,
        2.0,
        0,
        math.pi / 2.0,
        b.registry,
        "mm",
    )
    signal_cable_along_unit = geant4.solid.Box(
        f"{cable_name}_along_unit",
        cable_thickness,
        2.0,
        cable_length,
        b.registry,
        "mm",
    )
    signal_cable = geant4.solid.MultiUnion(
        cable_name,
        [signal_cable_clamp_to_curve, signal_cable_curve, signal_cable_along_unit],
        [
            [[0, 0, 0], [16 + 23.25 / 3 / 2.0, 0, 0]],
            [[+np.pi / 2, 0, 0], [16 + 23.25 / 3, 0, -signal_cable_radius - cable_thickness / 2.0]],
            [
                [0, 0, 0],
                [
                    16 + 23.25 / 3 + signal_cable_radius + cable_thickness / 2.0,
                    0,
                    -signal_cable_radius - cable_length / 2.0,
                ],
            ],
        ],
        b.registry,
    )

    return geant4.LogicalVolume(
        signal_cable,
        b.materials.metal_copper,
        cable_name,
        b.registry,
    )


def _get_signal_clamp_and_lmfe(
    clamp_thickness: float,
    b: core.InstrumentationData,
):
    holes = [[0.25, 5.05, 0], [0.25, -5.05, 0]]
    if "lmfe" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["ultem_clamp_signal"], b.registry.logicalVolumeDict["lmfe"], holes

    signal_clamp_mid = geant4.solid.Box(
        "ultem_clamp_signal_mid",
        7.5,
        9.5,
        clamp_thickness,
        b.registry,
        "mm",
    )
    signal_clamp_side = geant4.solid.Box(
        "ultem_clamp_signal_side",
        27,
        3,
        clamp_thickness,
        b.registry,
        "mm",
    )
    signal_clamp_bulk = geant4.solid.MultiUnion(
        "ultem_clamp_signal_bulk",
        [signal_clamp_mid, signal_clamp_side, signal_clamp_side],
        [
            [[0, 0, 0], [0, 0, 0]],
            [
                [0, 0, 0],
                [
                    (signal_clamp_side.pX - signal_clamp_mid.pX) / 2,
                    +(signal_clamp_mid.pY + signal_clamp_side.pY) / 2,
                    0,
                ],
            ],
            [
                [0, 0, 0],
                [
                    (signal_clamp_side.pX - signal_clamp_mid.pX) / 2,
                    -(signal_clamp_mid.pY + signal_clamp_side.pY) / 2,
                    0,
                ],
            ],
        ],
        b.registry,
    )

    clamp_hole = geant4.solid.Tubs(
        "signal_clamp_hole", 0, 1.5, clamp_thickness, 0, 2 * math.pi, b.registry, "mm"
    )
    clamp_springhole = geant4.solid.Box("signal_clamp_springhole", 10, 0.901, 0.251, b.registry, "mm")
    clamp_holes = geant4.solid.MultiUnion(
        "signal_clamp_holes",
        [clamp_hole, clamp_hole, clamp_springhole],
        [[[0, 0, 0], holes[0]], [[0, 0, 0], holes[1]], [[0, 0, 0], [0, 0, 0]]],
        b.registry,
    )
    signal_clamp = geant4.solid.Subtraction(
        "ultem_clamp_signal", signal_clamp_bulk, clamp_holes, [[0, 0, 0], [0, 0, 0]], b.registry
    )

    signal_lmfe = geant4.solid.Box("lmfe", 16, 8, 0.5, b.registry, "mm")

    signal_clamp_lv = geant4.LogicalVolume(signal_clamp, b.materials.ultem, "ultem_clamp_signal", b.registry)

    signal_lmfe_lv = geant4.LogicalVolume(
        signal_lmfe, b.materials.silica, "lmfe", b.registry
    )  # suprasil is quartz glass.

    return signal_clamp_lv, signal_lmfe_lv, holes


def _get_weldment(
    weldment_flap_thickness: float,
    b: core.InstrumentationData,
):
    if "hpge_support_copper_weldment" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["hpge_support_copper_weldment"]

    safety_margin = 0.1
    weldment_flap = geant4.solid.Box(
        "hpge_support_copper_weldment_flap",
        19,
        5,
        weldment_flap_thickness,
        b.registry,
        "mm",
    )

    weldment_clamp = geant4.solid.Box(
        "hpge_support_copper_weldment_clamp",
        6.7,
        5,
        1.5,  # 3 mm total height around the rod
        b.registry,
        "mm",
    )

    # Union the flap and clamp
    weldment_without_hole = geant4.solid.Union(
        "hpge_support_copper_weldment_without_hole",
        weldment_flap,
        weldment_clamp,
        [[0, 0, 0], [19 / 2.0 - 6.7 / 2.0, 0, -1.5 / 2.0 - weldment_flap_thickness / 2.0]],
        b.registry,
    )

    weldment_carving_hole = geant4.solid.Tubs(
        "hpge_support_copper_weldment_carving_hole",
        0,
        1.5 + safety_margin,
        2 * (weldment_flap_thickness + 2.2),
        0,
        math.pi * 2,
        b.registry,
        "mm",
    )

    # Perform subtraction only once
    weldment_top = geant4.solid.Subtraction(
        "hpge_support_copper_weldment",
        weldment_without_hole,
        weldment_carving_hole,
        [[0, 0, 0], [7, 0, 0]],  # Adjust the position of the hole as needed
        b.registry,
    )

    weldment_top_lv = geant4.LogicalVolume(
        weldment_top,
        b.materials.metal_copper,
        "hpge_support_copper_weldment",
        b.registry,
    )
    weldment_top_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)

    return weldment_top_lv


def _get_insulator(
    insulator_du_holder_flap_thickness: float,
    insulator_top_length: float,
    b: core.InstrumentationData,
):
    name_prefix = f"ultem_insulator_{insulator_top_length}"
    if name_prefix in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict[name_prefix]

    insulator_du_holder_flap = geant4.solid.Box(
        f"{name_prefix}_flap",
        16.5,
        7,
        insulator_du_holder_flap_thickness,
        b.registry,
        "mm",
    )

    safety_margin = 0.1
    safety_margin_touching_detector = 0.25
    top_length = insulator_top_length - safety_margin_touching_detector

    insulator_du_holder_clamp = geant4.solid.Box(
        f"{name_prefix}_clamp",
        top_length,
        7,
        5.5 - insulator_du_holder_flap_thickness,
        b.registry,
        "mm",
    )

    # Union the flap and clamp
    insulator_du_holder_without_hole = geant4.solid.Union(
        f"{name_prefix}_without_hole",
        insulator_du_holder_flap,
        insulator_du_holder_clamp,
        [
            [0, 0, 0],
            [
                16.5 / 2.0 - top_length / 2.0,
                0,
                (5.5 - insulator_du_holder_flap_thickness) / 2.0 + insulator_du_holder_flap_thickness / 2.0,
            ],
        ],
        b.registry,
    )

    insulator_du_holder_carving_hole = geant4.solid.Tubs(
        f"{name_prefix}_carving_hole",
        0,
        1.5 + safety_margin,
        3 * 5.5,
        0,
        math.pi * 2,
        b.registry,
        "mm",
    )

    # Perform subtraction only once
    insulator_du_holder = geant4.solid.Subtraction(
        name_prefix,
        insulator_du_holder_without_hole,
        insulator_du_holder_carving_hole,
        [[0, 0, 0], [16.5 / 2.0 - 1.5, 0, 0]],  # Adjust the position of the hole as needed
        b.registry,
    )

    insulator_du_holder_lv = geant4.LogicalVolume(
        insulator_du_holder,
        b.materials.ultem,
        name_prefix,
        b.registry,
    )
    insulator_du_holder_lv.pygeom_color_rgba = (0.64, 0.54, 0.31, 1)

    return insulator_du_holder_lv


def _get_cu_pin(length: float, b: core.InstrumentationData):
    if f"hpge_du_pin_{length}" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict[f"hpge_du_pin_{length}"]

    pin = geant4.solid.Tubs(f"hpge_du_pin_{length}", 0, 1.3 - 2e-10, length, 0, 2 * math.pi, b.registry)
    pin = geant4.LogicalVolume(pin, b.materials.metal_copper, f"hpge_du_pin_{length}", b.registry)
    pin.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)
    return pin


def _get_phbr_washer(thickness: float, b: core.InstrumentationData):
    if "phbr_washer" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["phbr_washer"]

    pin = geant4.solid.Tubs("phbr_washer", 1.4, 2.5, thickness, 0, 2 * math.pi, b.registry)
    pin = geant4.LogicalVolume(pin, b.materials.metal_phosphor_bronze, "phbr_washer", b.registry)
    pin.pygeom_color_rgba = (0.72, 0.5, 0.2, 1)
    return pin


def _get_phbr_spring(b: core.InstrumentationData):
    if "phbr_spring" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["phbr_spring"]

    spring = geant4.solid.Box("phbr_spring", 13, 0.9, 0.25, b.registry)
    spring = geant4.LogicalVolume(spring, b.materials.metal_phosphor_bronze, "phbr_spring", b.registry)
    spring.pygeom_color_rgba = (0.72, 0.5, 0.2, 1)
    return spring
