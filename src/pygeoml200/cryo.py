"""Construct the LEGEND-200/GERDA cryostat including the liquid argon volume.

Dimensions from figure 3 and table 2 of [Knoepfle2022]_, and from P. Krause. The two vessels
are built from torispherical heads of 3976 mm and 4160 mm inner diameter on cylindrical shells
of 3900 mm and 4149.2 mm, with necks of 789 mm and 964 mm inner diameter; the cryostat is
4200 mm across.

.. note::
    The shells here approximate the heads by half-ellipsoids of the same
    bounding box. This has faster navigation than
    a polycone (depending on the number of faces).

.. [Knoepfle2022] T. Knöpfle and B. Schwingenheuer "Design and Performance of the GERDA
   Low-Background Cryostat for Operation in Water" In: Journal of Instrumentation 17 P02038
   (2022). https://doi.org/10.1088/1748-0221/17/02/P02038
"""

from __future__ import annotations

from math import pi

import pyg4ometry.geant4 as g4
from pygeomtools import RemageDetectorInfo

from .utils import COLORS

cryo_radius = 3976 / 2
cryo_wall = 12
cryo_tub_height = 3900
cryo_top_height = 826
cryo_bottom_height = 829

# The outer vessel is not a scaled copy of the inner one: its cylindrical shell is longer
# (4149.2 vs 3900 mm) and its heads have their own heights, which is how the gap between the
# two vessels ends up ~150 mm axially while being only 80 mm radially.
cryo_vacuum_gap = (4160 - 4000) / 2  # inner vessel outer wall -> outer vessel bore
cryo_outer_wall = 20
cryo_outer_tub_height = 4149.2
cryo_gap_top_height = 862
cryo_gap_bottom_height = 870

# Necks: 789 mm inner diameter for the inner vessel, 964 mm for the outer one.
cryo_access_radius = 789 / 2
cryo_access_wall = 10
cryo_access_gap = (964 - 809) / 2  # inner neck outer wall -> outer neck bore
cryo_access_outer_wall = (1000 - 964) / 2
cryo_access_height = 1720
access_overlap = 200

# Outer envelope of the cryostat, i.e. the surfaces in contact with the water. These are
# annotated explicitly because pygeoml200 has an import cycle (utils -> core -> watertank ->
# utils), and mypy cannot infer the type of a derived constant across one.
cryo_outer_radius: float = cryo_radius + cryo_wall + cryo_vacuum_gap + cryo_outer_wall
cryo_outer_top_height: float = cryo_gap_top_height + cryo_outer_wall
cryo_outer_bottom_height: float = cryo_gap_bottom_height + cryo_outer_wall
cryo_access_outer_radius: float = (
    cryo_access_radius + cryo_access_wall + cryo_access_gap + cryo_access_outer_wall
)
cryo_neck_top: float = cryo_tub_height / 2 + cryo_top_height + cryo_access_height
# lowest point of the cryostat, relative to the centre of the inner vessel.
cryo_outer_bottom_z: float = -(cryo_outer_tub_height / 2 + cryo_outer_bottom_height)

lar_ullage_height = 800
lar_ullage_safety = 0.0001  # avoid surface overlaps


def cryostat_shell(
    name: str,
    reg: g4.Registry,
    *,
    r: float,
    h_cyl: float,
    h_top: float,
    h_bot: float,
    r_neck: float,
    z_neck: float,
    neck_height: float | None = None,
) -> g4.solid.Union:
    """Construct a cylinder with two half-ellipsoid heads and a neck, as a single union.

    Used for the three cryostat shells, for its VM2000 wrapper and for the argon. The shells
    are concentric but not similar: ``h_cyl`` is that shell's own cylindrical length, so the
    outer vessel's heads start further out than the inner vessel's. ``z_neck`` is the top of
    the neck; it defaults to the full height of the cryostat neck.
    """
    if neck_height is None:
        neck_height = cryo_access_height + access_overlap
    tub = g4.solid.Tubs(f"{name}_tub", 0, r, h_cyl, 0, 2 * pi, reg, "mm")
    top = g4.solid.Ellipsoid(f"{name}_top", r, r, h_top, 0, h_top, reg, "mm")
    bottom = g4.solid.Ellipsoid(f"{name}_bottom", r, r, h_bot, 0, h_bot, reg, "mm")
    neck = g4.solid.Tubs(f"{name}_neck", 0, r_neck, neck_height, 0, 2 * pi, reg, "mm")

    shell1 = g4.solid.Union(f"{name}1", tub, top, [[0, 0, 0], [0, 0, h_cyl / 2]], reg)
    shell2 = g4.solid.Union(f"{name}2", shell1, bottom, [[0, pi, 0], [0, 0, -h_cyl / 2]], reg)
    return g4.solid.Union(name, shell2, neck, [[0, 0, 0], [0, 0, z_neck - neck_height / 2]], reg)


def construct_cryostat(
    cryostat_material: g4.Material, vacuum_material: g4.Material, reg: g4.Registry
) -> g4.LogicalVolume:
    """Construct the double-walled cryostat.

    The returned outer vessel contains the insulation vacuum, which in turn contains the inner
    vessel. The argon is placed into the inner vessel, see :func:`place_argon`.
    """
    inner = cryostat_shell(
        "cryostat_inner_wall",
        reg,
        r=cryo_radius + cryo_wall,
        h_cyl=cryo_tub_height,
        h_top=cryo_top_height + cryo_wall,
        h_bot=cryo_bottom_height + cryo_wall,
        r_neck=cryo_access_radius + cryo_access_wall,
        z_neck=cryo_neck_top - 2e-6,
    )
    gap = cryostat_shell(
        "cryostat_vacuum_gap",
        reg,
        r=cryo_radius + cryo_wall + cryo_vacuum_gap,
        h_cyl=cryo_outer_tub_height,
        h_top=cryo_gap_top_height,
        h_bot=cryo_gap_bottom_height,
        r_neck=cryo_access_radius + cryo_access_wall + cryo_access_gap,
        z_neck=cryo_neck_top - 1e-6,
    )
    outer = cryostat_shell(
        "cryostat_outer_wall",
        reg,
        r=cryo_outer_radius,
        h_cyl=cryo_outer_tub_height,
        h_top=cryo_outer_top_height,
        h_bot=cryo_outer_bottom_height,
        r_neck=cryo_access_outer_radius,
        z_neck=cryo_neck_top,
    )

    outer_lv = g4.LogicalVolume(outer, cryostat_material, "cryostat_outer_wall", reg)
    gap_lv = g4.LogicalVolume(gap, vacuum_material, "cryostat_vacuum_gap", reg)
    inner_lv = g4.LogicalVolume(inner, cryostat_material, "cryostat_inner_wall", reg)
    gap_lv.pygeom_color_rgba = False
    inner_lv.pygeom_color_rgba = COLORS["steel"]

    g4.PhysicalVolume([0, 0, 0], [0, 0, 0], gap_lv, "cryostat_vacuum_gap", outer_lv, reg)
    inner_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], inner_lv, "cryostat_inner_wall", gap_lv, reg)
    # the argon goes into the inner vessel, not into the returned outer one.
    outer_lv.pygeom_cryostat_inner = (inner_lv, inner_pv)

    return outer_lv


def place_cryostat(
    cryostat_lv: g4.LogicalVolume, wl: g4.LogicalVolume, reg: g4.Registry
) -> g4.PhysicalVolume:
    """Place the cryostat at the origin of ``wl``"""
    cryostat_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], cryostat_lv, "cryostat_outer_wall", wl, reg)
    cryostat_lv.pygeom_color_rgba = COLORS["steel"]
    return cryostat_pv


def construct_argon(lar_material: g4.Material, reg: g4.Registry) -> tuple[g4.LogicalVolume, float]:
    """Construct an approximate LEGEND-200 argon volume.

    Returns
    -------
    logical volume instance and height of the cryostat neck relative to the origin of the argon volume.

    .. note::
        the constructed volume's center (i.e. for children placed at 0,0,0) is not the barycenter of the
        volume, but the center of the central tubular section of the cryostat.
    """
    lar_access_height = cryo_access_height - lar_ullage_height
    lar = cryostat_shell(
        "liquid_argon",
        reg,
        r=cryo_radius,
        h_cyl=cryo_tub_height,
        h_top=cryo_top_height,
        h_bot=cryo_bottom_height,
        r_neck=cryo_access_radius,
        z_neck=cryo_tub_height / 2 + cryo_top_height + lar_access_height,
        neck_height=lar_access_height + access_overlap,
    )

    lar_neck_z = (
        cryo_tub_height / 2 + cryo_top_height - 20
    )  # offset is below the "virtual" top point of the round segment (see technical drawing)
    return g4.LogicalVolume(lar, lar_material, "liquid_argon", reg), lar_neck_z


def construct_ullage_argon(gar_material: g4.Material, reg: g4.Registry) -> g4.LogicalVolume:
    """Construct the gaseous argon volume above the LAr."""
    lar_ullage = g4.solid.Tubs(
        "gaseous_argon",
        0,
        cryo_access_radius,
        lar_ullage_height - 4 * lar_ullage_safety,
        0,
        2 * pi,
        reg,
        "mm",
    )
    return g4.LogicalVolume(lar_ullage, gar_material, "gaseous_argon", reg)


def place_argon(
    lar_lv: g4.LogicalVolume,
    cryostat_lv: g4.LogicalVolume,
    cryostat_displacement_z: float,
    to_cryostat_steel: g4.solid.OpticalSurface,
    reg: g4.Registry,
) -> g4.PhysicalVolume:
    """Place the liquid argon volume in the inner vessel of the cryostat.

    Also adds an optical surface in between and registers the argon as active detector."""
    inner_lv, inner_pv = cryostat_lv.pygeom_cryostat_inner
    lar_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, cryostat_displacement_z], lar_lv, "liquid_argon", inner_lv, reg
    )
    lar_lv.pygeom_color_rgba = [0, 0, 0, 0.03]

    # set lar as active with det id 0
    lar_pv.set_pygeom_active_detector(RemageDetectorInfo("scintillator", 0, {}))

    # add surface argon->steel
    g4.BorderSurface("bsurface_lar_steel", lar_pv, inner_pv, to_cryostat_steel, reg)

    return lar_pv


def place_ullage_argon(
    gar_lv: g4.LogicalVolume,
    cryostat_lv: g4.LogicalVolume,
    cryostat_displacement_z: float,
    reg: g4.Registry,
) -> g4.PhysicalVolume:
    z_pos = (
        cryo_tub_height / 2 + cryo_top_height + cryo_access_height - lar_ullage_height / 2 + lar_ullage_safety
    )
    inner_lv, _ = cryostat_lv.pygeom_cryostat_inner
    gar_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, z_pos + cryostat_displacement_z], gar_lv, "gaseous_argon", inner_lv, reg
    )
    gar_lv.pygeom_color_rgba = False

    return gar_pv
