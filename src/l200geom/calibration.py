from __future__ import annotations

import logging
import math
from typing import Literal

import numpy as np
from dbetto import AttrsDict
from pyg4ometry import geant4

from . import core, hpge_strings

log = logging.getLogger(__name__)


def place_calibration_system(b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 calibration system."""

    # overridable config for copper cap dimensions.
    cu_absorber_config = {
        "height": 15 - 0.01,  # mm, drawing from S. Schönert
        "inner_height": 12.5,  # mm
        "inner_radius": 7.6 / 2,  # mm
        "outer_radius": (14 - 0.01) / 2,  # mm, drawing from S. Schönert
    }
    cu_absorber_config = b.runtime_config.get("cu_absorber", AttrsDict(cu_absorber_config))

    # add source outside SIS, if requested by user.
    if hasattr(b.runtime_config, "extra_source"):
        source_cfg = b.runtime_config.extra_source
        source_spec = _parse_source_spec(source_cfg.source)
        _place_source(
            b,
            source_cfg.get("name", ""),
            np.array(source_cfg.position_in_mm[0:2]),
            source_cfg.position_in_mm[2],
            source_type=source_spec["type"],
            cu_absorber=(cu_absorber_config if source_spec["has_cu"] else None),
            bare=source_spec["bare"],
        )

    # place calibration tubes.
    if "calibration" not in b.special_metadata or len(b.special_metadata.calibration) == 0:
        return

    calib_tubes = {}
    calib_tube_length = []
    calib_tube_xy = np.empty((2, 4))

    sis_cfg = b.runtime_config.get("sis", {})

    for i, tube in b.special_metadata.calibration.items():
        idx = int(i) - 1
        if tube is None:
            continue
        tube_cfg = sis_cfg.get(i, {}) or {}

        if tube.length_in_mm not in calib_tubes:
            calib_tubes[tube.length_in_mm] = hpge_strings._get_nylon_mini_shroud(
                tube.tube_radius_in_mm, tube.length_in_mm, True, b.materials, b.registry
            )
            calib_tube_length.append(tube.length_in_mm)

        # allow for an offset to place properly the sis
        phi = np.deg2rad(tube.angle_in_deg)
        phi += np.deg2rad(tube_cfg.get("phi_offset", 0))
        tube.radius_in_mm += tube_cfg.get("r_offset", 0)

        # add a very small z-offset to prevent overlaps if we moved a cal tube
        off = 2 if ("phi_offset" in tube_cfg or "r_offset" in tube_cfg) else 0

        calib_tube_xy[:, idx] = np.array([tube.radius_in_mm * np.cos(phi), -tube.radius_in_mm * np.sin(phi)])

        nms_pv = geant4.PhysicalVolume(
            [0, 0, 0],
            [*calib_tube_xy[:, idx], b.top_plate_z_pos - tube.length_in_mm / 2 - off],
            calib_tubes[tube.length_in_mm],
            f"calibration_tube_{i}",
            b.mother_lv,
            b.registry,
        )
        hpge_strings._add_nms_surfaces(nms_pv, b.mother_pv, b.materials, b.registry)

    # check if we have one shared length of calibration tubes.
    calib_tube_length = (
        calib_tube_length[0] if all(length == calib_tube_length[0] for length in calib_tube_length) else None
    )

    # place the actual calibration sources inside.
    if not hasattr(b.runtime_config, "sis"):
        return

    for i, tube in b.special_metadata.calibration.items():
        if tube is None or i not in sis_cfg or sis_cfg[i] is None:
            continue
        idx = int(i) - 1

        # SIS reading to our coordinates. This marks the top of the torlon initialization pin in our
        # (pygeom) coordinates.

        sis_z = sis_cfg[i].sis_z - sis_cfg[i].get("offset", 0)
        sis_xy = calib_tube_xy[:, idx]

        pin_top = _sis_to_pygeoml200(sis_z)

        if len(sis_cfg[i].sources) != 4:
            msg = f"Invalid number of sources in config of SIS{i}"
            raise ValueError(msg)

        # z offsets from top of pin to bottom of source.
        delta_z = (-271, -171, -71, 42 + source_inside_holder)

        # always place the Ta absorber, irrespective if it holds a source.
        _place_ta_absorber(b, f"_sis{i}", sis_xy, pin_top + delta_z[3] - source_inside_holder)

        for si in range(4):
            if sis_cfg[i].sources[si] is None:
                continue
            source_spec = _parse_source_spec(sis_cfg[i].sources[si])
            _place_source(
                b,
                f"_sis{i}_source{si}",
                sis_xy,
                pin_top + delta_z[si] - source_inside_holder,
                source_type=source_spec["type"],
                cu_absorber=(cu_absorber_config if source_spec["has_cu"] else None),
                bare=source_spec["bare"],
            )


# outer dimensions of steel container:
source_height = 17.6  # mm
source_radius_outer = 6.4 / 2  # mm
# inner dimension steel container (i.e. the actual source size):
source_th_height_inner = 5  # mm
source_th_radius_inner = 3 / 2  # mm
source_th_top_inner = 2.3  # mm

source_height_inner = 4  # mm
source_radius_inner = 4 / 2  # mm
source_top_inner = 1.2  # mm

ABSORBER_HEIGHT = 37.5  # mm
# overlap of steel container with Ta absorber/source holder:
source_outside_holder = 10.6  # mm
source_inside_holder = source_height - source_outside_holder

safety = 1e-9  # mm


def _place_source(
    b: core.InstrumentationData,
    suffix: str,
    xy,
    delta_z: float,
    source_type: Literal["Th228", "Ra"],
    cu_absorber: AttrsDict | None,
    bare: bool = False,
) -> None:
    """Place a single source container.

    delta_z
        to source container holder top from top plate top
    source_type
        controls the interior design of the source container
    cu_absorber
        include a copper absorber cap of the given dimensions. The dimensions have to be the same for all
        sources in this geometry.
    bare
        Do not encapsulate the source. only use if you know what you do; this does not correspond to any
        physical source geometry used in LEGEND.
    """
    z0 = b.top_plate_z_pos - delta_z

    source_z = z0 + source_height / 2 - source_inside_holder
    # outer = steel container.
    if not bare:
        if "source_outer" not in b.registry.solidDict:
            geant4.solid.Tubs(
                "source_outer",
                0,
                source_radius_outer - safety,
                source_height - safety,
                0,
                2 * math.pi,
                b.registry,
            )
        if f"source_outer{suffix}" not in b.registry.logicalVolumeDict:
            geant4.LogicalVolume(
                b.registry.solidDict["source_outer"],
                b.materials.metal_steel,
                f"source_outer{suffix}",
                b.registry,
            )
        source_outer = b.registry.logicalVolumeDict[f"source_outer{suffix}"]
        geant4.PhysicalVolume(
            [0, 0, 0],
            [*xy, source_z],
            source_outer,
            f"source_outer{suffix}",
            b.mother_lv,
            b.registry,
        )

    # inner = contains actual source material.
    inner_dims = (source_radius_inner, source_height_inner, source_top_inner)
    if source_type == "Th228":
        inner_dims = (source_th_radius_inner, source_th_height_inner, source_th_top_inner)

    if f"source_inner_{source_type}" not in b.registry.solidDict:
        source_inner_solid = geant4.solid.Tubs(
            f"source_inner_{source_type}", 0, *inner_dims[0:2], 0, 2 * math.pi, b.registry
        )

        if source_type == "Th228":
            source_material = b.materials.metal_caps_gold  # for Th source
        elif source_type == "Ra":
            source_material = geant4.MaterialPredefined("G4_SILICON_DIOXIDE")  # for Ra source
        else:
            msg = f"unknown source type {source_type}"
            raise ValueError(msg)

        source_inner = geant4.LogicalVolume(
            source_inner_solid, source_material, f"source_inner_{source_type}", b.registry
        )
        source_inner.pygeom_color_rgba = (1, 0.843, 0, 1)

    source_inner_z = source_height / 2 - inner_dims[1] / 2 - inner_dims[2]
    geant4.PhysicalVolume(
        [0, 0, 0],
        [0, 0, source_inner_z] if not bare else [*xy, source_z + source_inner_z],
        b.registry.logicalVolumeDict[f"source_inner_{source_type}"],
        f"source_inner{suffix}",
        source_outer if not bare else b.mother_lv,
        b.registry,
    )

    # copper absorber cap, if requested.
    if cu_absorber:
        cu_absorber_cu, cu_absorber_lar = _get_cu_cap(b, cu_absorber)
        cu_absorber_z0 = z0 + max([0, source_outside_holder - cu_absorber.inner_height])
        cu_absorber_z = cu_absorber_z0 + cu_absorber.height / 2 + safety
        cu_absorber_lar_z = cu_absorber_z0 + cu_absorber.inner_height / 2 + safety
        geant4.PhysicalVolume(
            [0, 0, 0],
            [*xy, cu_absorber_z],
            cu_absorber_cu,
            f"cu_absorber{suffix}",
            b.mother_lv,
            b.registry,
        )
        if cu_absorber_lar is not None:
            geant4.PhysicalVolume(
                [0, 0, 0],
                [*xy, cu_absorber_lar_z],
                cu_absorber_lar,
                f"cu_absorber_lar_inactive{suffix}",
                b.mother_lv,
                b.registry,
            )


def _get_cu_cap(
    b: core.InstrumentationData, dimens: AttrsDict
) -> tuple[geant4.LogicalVolume, geant4.LogicalVolume | None]:
    if "cu_absorber" in b.registry.logicalVolumeDict:
        return (
            b.registry.logicalVolumeDict["cu_absorber"],
            b.registry.logicalVolumeDict.get("cu_absorber_lar_inactive", None),
        )

    if (
        dimens.inner_radius >= dimens.outer_radius
        or dimens.inner_height >= dimens.height
        or dimens.inner_radius < source_radius_outer
    ):
        msg = "invalid copper cap configuration."
        raise ValueError(msg)

    cu_absorber_outer = geant4.solid.Tubs(
        "cu_absorber_outer", 0, dimens.outer_radius, dimens.height, 0, 2 * math.pi, b.registry
    )
    cu_absorber_inner = geant4.solid.Tubs(
        "cu_absorber_inner", 0, dimens.inner_radius, dimens.inner_height, 0, 2 * math.pi, b.registry
    )
    cu_absorber = geant4.solid.Subtraction(
        "cu_absorber",
        cu_absorber_outer,
        cu_absorber_inner,
        [[0, 0, 0], [0, 0, -(dimens.height - dimens.inner_height) / 2]],
        b.registry,
    )
    cu_absorber = geant4.LogicalVolume(cu_absorber, b.materials.metal_copper, "cu_absorber", b.registry)
    cu_absorber.pygeom_color_rgba = (0.72, 0.45, 0.2, 0.3)

    cu_absorber_lar = None
    if dimens.inner_radius != source_radius_outer:
        cu_absorber_lar_outer = geant4.solid.Tubs(
            "cu_absorber_lar_inactive_outer",
            0,
            dimens.inner_radius - safety,
            dimens.inner_height - 2 * safety,
            0,
            2 * math.pi,
            b.registry,
        )
        cu_absorber_lar_inner = geant4.solid.Tubs(
            "cu_absorber_lar_inactive_inner",
            0,
            source_radius_outer,
            source_outside_holder,
            0,
            2 * math.pi,
            b.registry,
        )
        cu_absorber_lar = geant4.solid.Subtraction(
            "cu_absorber_lar_inactive",
            cu_absorber_lar_outer,
            cu_absorber_lar_inner,
            [[0, 0, 0], [0, 0, -(dimens.inner_height - 2 * safety - source_outside_holder) / 2]],
            b.registry,
        )
        cu_absorber_lar = geant4.LogicalVolume(
            cu_absorber_lar, b.materials.liquidargon, "cu_absorber_lar_inactive", b.registry
        )
        cu_absorber_lar.pygeom_color_rgba = (1, 1, 1, 0.0001)

    return cu_absorber, cu_absorber_lar


def _place_ta_absorber(
    b: core.InstrumentationData,
    suffix: str,
    xy,
    delta_z: float,
) -> None:
    """Place tantalum absorber plus source container.

    delta_z
        to absorber top from top plate top
    """
    z0 = b.top_plate_z_pos - delta_z - ABSORBER_HEIGHT / 2

    ta_absorber_lv = _get_ta_absorber(b)
    ta_absorber_lv.pygeom_color_rgba = (0.5, 0.5, 0.5, 0.9)
    geant4.PhysicalVolume(
        [0, 0, 0], [*xy, z0], ta_absorber_lv, f"ta_absorber{suffix}", b.mother_lv, b.registry
    )

    if "peek_holder" not in b.registry.logicalVolumeDict:
        peek_outside = geant4.solid.Box("peek_outside", 33.1, 9, 25, b.registry)
        peek_inside = geant4.solid.Box("peek_inside", 14, 9, 15, b.registry)  # Drawing from S. Schönert
        peek_holder = geant4.solid.Subtraction(
            "peek_holder", peek_outside, peek_inside, [[0, 0, 0], [0, 0, -10 / 2]], b.registry
        )
        peek_holder = geant4.LogicalVolume(peek_holder, b.materials.peek, "peek_holder", b.registry)
        peek_holder.pygeom_color_rgba = (0.5, 0.32, 0, 1)

    peek_holder_z = z0 + ABSORBER_HEIGHT / 2 + 25 / 2 + safety
    geant4.PhysicalVolume(
        [0, 0, 0],
        [*xy, peek_holder_z],
        b.registry.logicalVolumeDict["peek_holder"],
        f"peek_holder{suffix}",
        b.mother_lv,
        b.registry,
    )


def _get_ta_absorber(b: core.InstrumentationData):
    if "ta_absorber" in b.registry.logicalVolumeDict:
        return b.registry.logicalVolumeDict["ta_absorber"]

    ta_absorber_outer = geant4.solid.Tubs(
        "ta_absorber_outer",
        0,
        16.2,  # estimate!
        ABSORBER_HEIGHT,
        0,
        2 * math.pi,
        b.registry,
    )
    ta_absorber_inner = geant4.solid.Tubs(
        "ta_absorber_inner",
        0,
        source_radius_outer + safety,
        source_height - source_outside_holder + 2 * safety,
        0,
        2 * math.pi,
        b.registry,
    )
    ta_absorber = geant4.solid.Subtraction(
        "ta_absorber",
        ta_absorber_outer,
        ta_absorber_inner,
        [[0, 0, 0], [0, 0, (ABSORBER_HEIGHT - source_inside_holder) / 2 - safety]],
        b.registry,
    )
    return geant4.LogicalVolume(ta_absorber, b.materials.metal_tantalum, "ta_absorber", b.registry)


def _sis_to_pygeoml200(sis_coord: float) -> float:
    sis_ta_on_lar = 5315  # mm
    h_absorber = 82.5  # mm - from Pin-Jung
    h_funnel_to_lar_surface = 7333 - 5207  # mm (meterdrive readings from Matt's slides)
    h_funnel = 90  # mm - from Matt's slides; _not_ including the Cu top plate thickness.
    sis_top_plate = sis_ta_on_lar + h_absorber + h_funnel_to_lar_surface + h_funnel
    return sis_coord - sis_top_plate


def _parse_source_spec(src: str) -> dict:
    src = src.split("+")
    if src[0] not in ("Th228", "Ra"):
        msg = f"Invalid source type {src[0]} in source spec"
        raise ValueError(msg)
    if set(src[1:]) - {"Cu", "_bare"} != set():
        msg = f"Unknown extra in source spec {src[1:]}"
        raise ValueError(msg)
    return {"type": src[0], "has_cu": "Cu" in src, "bare": "_bare" in src}
