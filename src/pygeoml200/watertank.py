"""Construct the LEGEND-200/GERDA watertank including the water volume and the reflective foil VM2000."""

from __future__ import annotations

import math
import warnings
from math import pi, sqrt

import numpy as np
import pyg4ometry.geant4 as g4
from pygeomtools import RemageDetectorInfo
from scipy.spatial.transform import Rotation

from . import cryo, materials
from .utils import COLORS

pmt_id = np.array(
    [
        "ch29_312",
        "ch18_301",
        "ch16_701",
        "ch19_302",
        "ch20_303",
        "ch21_304",
        "ch22_305",
        "ch23_306",
        "ch17_703",
        "ch24_307",
        "ch25_308",
        "ch26_309",
        "ch27_310",
        "ch28_311",
        "ch15_208",
        "ch10_201",
        "ch11_202",
        "ch12_203",
        "ch13_706",
        "ch14_206",
        "ch01_704",
        "ch03_705",
        "ch06_709",
        "ch08_710",
        "ch30_401",
        "ch31_402",
        "ch32_403",
        "ch33_404",
        "ch34_409",
        "ch35_410",
        "ch43_510",
        "ch36_501",
        "ch37_502",
        "ch38_503",
        "ch39_504",
        "ch40_507",
        "ch41_508",
        "ch42_509",
        "ch44_602",
        "ch45_603",
        "ch46_605",
        "ch47_606",
        "ch48_607",
        "ch49_608",
        "ch50_609",
        "ch51_610",
        "ch52_702",
        "ch00_101",
        "ch02_102",
        "ch04_708",
        "ch05_104",
        "ch07_105",
        "ch09_707",
    ]
)

pmt_rawids = np.array(
    [
        2003219,
        2003208,
        2003206,
        2003209,
        2003210,
        2003211,
        2003212,
        2003213,
        2003207,
        2003214,
        2003215,
        2003216,
        2003217,
        2003218,
        2003205,
        2003200,
        2003201,
        2003202,
        2003203,
        2003204,
        2001605,
        2001607,
        2001610,
        2001612,
        2004800,
        2004801,
        2004802,
        2004803,
        2004804,
        2004805,
        2004813,
        2004806,
        2004807,
        2004808,
        2004809,
        2004810,
        2004811,
        2004812,
        2004814,
        2004815,
        2004816,
        2004817,
        2004818,
        2004819,
        2004820,
        2004821,
        2004822,
        2001604,
        2001606,
        2001608,
        2001609,
        2001611,
        2001613,
    ]
)


# Water tank with water and air buffer
water_tank_thickness = 7.0
inner_tank_height = 8900.0
inner_radius = 0.0

# Reflective foil
reflective_foil_thickness = 0.04

# The tank cavity is lined with VM2000, and the water sits inside that liner.
water_tank_inner_radius = 5000.0
water_tank_inner_height = inner_tank_height - 2 * water_tank_thickness
water_radius = water_tank_inner_radius - reflective_foil_thickness
water_height = water_tank_inner_height - 2 * reflective_foil_thickness

# Pillbox: the cylindrical skirt the cryostat stands on. [Knoepfle2022]_ gives "12 mm thick
# stainless steel [holding] the cryostat 1280 mm above floor", the drawing "ø4200 x 12".
shielding_foot_or = cryo.cryo_outer_radius
shielding_foot_thickness = 12.0
shielding_foot_ir = shielding_foot_or - shielding_foot_thickness
manhole_outer_radius = 400.0
cryo_clearance = 1280.0  # water floor to the lowest point of the cryostat

# The skirt is welded to the outer vessel where its cylindrical shell ends, and the bottom head
# then curves away *inside* it. Stop it just below
# the point where the head, foil included, has narrowed clear of the skirt's inner surface.
pillbox_head_clearance = 0.5  # radial play between the skirt and the bottom head.
_head_clear_radius = (
    shielding_foot_ir - 2 * reflective_foil_thickness - pillbox_head_clearance
)  # skirt bore, its own liner, the cryostat's foil, and the play in between
cryo_bottom_height = (
    cryo_clearance
    + cryo.cryo_outer_bottom_height
    - cryo.cryo_outer_bottom_height * sqrt(1 - (_head_clear_radius / cryo.cryo_outer_radius) ** 2)
    - reflective_foil_thickness  # the skirt's wrapper reaches one foil thickness higher
)

# This is what fixes the water tank in the world
cryo_z_displacement = (
    -water_height / 2 + cryo_clearance + reflective_foil_thickness - cryo.cryo_outer_bottom_z
)
# centre of the pillbox's VM2000 wrapper, which stands on the water floor and is one foil
# thickness taller than the skirt at each end (the skirt sits inside it).
pillbox_foil_height = cryo_bottom_height + 2 * reflective_foil_thickness
pillbox_offset = -water_height / 2 + 0.5 * pillbox_foil_height + 1e-9

# Air buffer
outer_water_tank_radius = water_radius + water_tank_thickness
air_buffer_radius = water_radius - 1e-9
air_buffer_height = 486.0


# z-axis Offsets
air_buffer_offset = 0.5 * (water_height - air_buffer_height)


# PMTs

angle_segment = 69.449  # PMT cap adaption

pmt_starting_angle = 0.0  # in degrees
pmt_ending_angle = 2 * pi  # in degrees

photocathode_inner_radius = 109.84
photocathode_outer_radius = 110.00
photocathode_theta_start = 0.0  # in degrees
photocathode_theta_end = (angle_segment / 180.0) * pi  # in degrees
photocatode_height = photocathode_outer_radius * (1 - math.cos(photocathode_theta_end))
photocathode_height_difference = photocathode_outer_radius - photocatode_height

pmt_outer_radius = 103.0
borosilikat_glass_thickness = 1.0

pmt_steel_cone_thickness = 4.0
pmt_steel_cone_height = 300.0
pmt_steel_cone_upper_rmin = pmt_outer_radius
pmt_steel_cone_upper_rmax = pmt_outer_radius + pmt_steel_cone_thickness
pmt_steel_cone_lower_rmin = pmt_outer_radius / 3.0
pmt_steel_cone_lower_rmax = pmt_steel_cone_lower_rmin + pmt_steel_cone_thickness

pmt_borosilikat_glass_outer_radius = photocathode_outer_radius + borosilikat_glass_thickness

pmt_air_outer_radius = 114.08

acryl_inner_radius = photocathode_inner_radius
acryl_outer_radius = pmt_air_outer_radius + 6.0
acryl_theta_start = 0.0  # in degrees
acryl_theta_end = (angle_segment / 180.0) * pi  # in degrees


pmt_steel_bottom_height = 30.0

pmt_cathode_offset = (
    -water_height / 2
    + reflective_foil_thickness
    - photocathode_height_difference
    + 0.5 * pmt_steel_cone_height
    + pmt_steel_bottom_height
)
pmt_cone_offset = (
    -water_height / 2
    + reflective_foil_thickness
    + 0.25 * pmt_steel_cone_height
    + pmt_steel_bottom_height
    + 2e-9
)
pmt_bottom_offset = -water_height / 2 + reflective_foil_thickness + 0.5 * pmt_steel_bottom_height + 2e-9

distance_pmt_base_tank = (
    water_radius
    - reflective_foil_thickness
    - np.sqrt((water_radius - reflective_foil_thickness) ** 2 - pmt_steel_cone_lower_rmax**2)
)

distance_pmt_base_pillbox = (
    shielding_foot_ir
    - reflective_foil_thickness
    - np.sqrt((shielding_foot_ir - reflective_foil_thickness) ** 2 - pmt_steel_cone_lower_rmax**2)
)


def construct_tank(reg: g4.Registry, tank_material: g4.Material) -> g4.LogicalVolume:
    water_tank_wall = g4.solid.Tubs(
        "water_tank_wall", inner_radius, outer_water_tank_radius, inner_tank_height, 0, 2 * pi, reg
    )

    tank_lv = g4.LogicalVolume(water_tank_wall, tank_material, "water_tank_lv", reg)
    tank_lv.pygeom_color_rgba = COLORS["steel"]
    return tank_lv


def place_tank(
    reg: g4.Registry,
    water_tank_lv: g4.LogicalVolume,
    world_lv: g4.LogicalVolume,
    tank_offset: float,
) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0, 0, 0], [0, 0, tank_offset], water_tank_lv, "water_tank", world_lv, reg)


def construct_tank_foil(reg: g4.Registry, vm2000_material: g4.Material) -> g4.LogicalVolume:
    """VM2000 lining the whole tank cavity; the water is placed inside it."""
    foil = g4.solid.Tubs(
        "water_tank_reflection_foil",
        inner_radius,
        water_tank_inner_radius,
        water_tank_inner_height,
        0,
        2 * pi,
        reg,
    )
    foil_lv = g4.LogicalVolume(foil, vm2000_material, "water_tank_reflection_foil_lv", reg)
    foil_lv.pygeom_color_rgba = COLORS["vm2000"]
    return foil_lv


def place_tank_foil(
    reg: g4.Registry, foil_lv: g4.LogicalVolume, tank_lv: g4.LogicalVolume
) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0, 0, 0], [0, 0, 0], foil_lv, "water_tank_reflection_foil_pv", tank_lv, reg)


def construct_water(reg: g4.Registry, water_material: g4.Material) -> g4.LogicalVolume:
    water_solid = g4.solid.Tubs("water_solid", inner_radius, water_radius, water_height, 0, 2 * pi, reg)
    water_lv = g4.LogicalVolume(water_solid, water_material, "water_lv", reg)
    water_lv.pygeom_color_rgba = COLORS["water"]
    return water_lv


def place_water(
    reg: g4.Registry, water_lv: g4.LogicalVolume, tank_foil_lv: g4.LogicalVolume
) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0, 0, 0], [0, 0, 0], water_lv, "water_pv", tank_foil_lv, reg)


def construct_air_buffer(reg: g4.Registry, air_material: g4.Material) -> g4.LogicalVolume:
    air_buffer = g4.solid.Tubs(
        "air_buffer",
        cryo.cryo_access_outer_radius + reflective_foil_thickness + 1e-9,
        air_buffer_radius,
        air_buffer_height,
        0,
        2 * pi,
        reg,
    )
    air_buffer_lv = g4.LogicalVolume(air_buffer, air_material, "air_buffer_lv", reg)
    air_buffer_lv.pygeom_color_rgba = COLORS["air"]
    return air_buffer_lv


def place_air_buffer(
    reg: g4.Registry, air_buffer_lv: g4.LogicalVolume, water_lv: g4.LogicalVolume
) -> g4.PhysicalVolume:
    return g4.PhysicalVolume(
        [0, 0, 0], [0, 0, air_buffer_offset - 1e-9], air_buffer_lv, "air_buffer_pv", water_lv, reg
    )


def _manhole_solid(name: str, radius: float, reg: g4.Registry) -> g4.solid.Union:
    """The arch-shaped manhole cut-out: a half-cylinder on a rectangle."""
    # Make the subtractions slightly deeper than the pillbox depth for sanity
    height = 2 * (shielding_foot_or + reflective_foil_thickness + 1)
    arc = g4.solid.Tubs(name + "_arc", 0, radius, height, 0, 2 * math.pi, reg)
    box = g4.solid.Box(name + "_box", 2 * radius, radius, height, reg)
    return g4.solid.Union(name, arc, box, [[0, 0, 0], [0, -0.5 * radius, 0]], reg)


def construct_pillbox(
    reg: g4.Registry, pillbox_material: g4.Material | str
) -> tuple[g4.LogicalVolume, list, float]:
    x_rot_drehvol = 0
    y_rot_drehvol = np.pi / 2.0
    z_rot_drehvol = np.pi / 2.0

    # rotation matrix around Z, Y and X
    rot_z = Rotation.from_euler("z", z_rot_drehvol, degrees=False)
    rot_y = Rotation.from_euler("y", y_rot_drehvol, degrees=False)
    rot_x = Rotation.from_euler("x", x_rot_drehvol, degrees=False)

    # combine rotations
    combined_rotation = rot_y * rot_z * rot_x

    # rotation matrix in euler angle
    with warnings.catch_warnings(action="ignore"):
        euler_angles = combined_rotation.as_euler("xyz", degrees=False)
    x_rot_global, y_rot_global, z_rot_global = euler_angles

    pillbox_tube = g4.solid.Tubs(
        "pillbox_tube",
        shielding_foot_ir,
        shielding_foot_or,
        cryo_bottom_height,
        0,
        2 * pi,
        reg,
    )  # outer steel cylinder

    manhole_rotation = [x_rot_global, y_rot_global, z_rot_global]
    manhole_pillbox = _manhole_solid("manhole_union", manhole_outer_radius, reg)

    # Subtract the manhole (half-tube) from the pillbox
    man_hole_offset = 0.5 * cryo_bottom_height - manhole_outer_radius
    pillbox = g4.solid.Subtraction(
        "pillbox_subtraction1",
        pillbox_tube,
        manhole_pillbox,
        [manhole_rotation, [0, 0, 0 - man_hole_offset]],
        reg,
    )
    pillbox_lv = g4.LogicalVolume(pillbox, pillbox_material, "pillbox_lv", reg)
    pillbox_lv.pygeom_color_rgba = COLORS["steel"]

    return pillbox_lv, manhole_rotation, man_hole_offset


def place_pillbox(
    reg: g4.Registry, pillbox_lv: g4.LogicalVolume, pillbox_foil_lv: g4.LogicalVolume
) -> g4.PhysicalVolume:
    """The skirt sits inside its VM2000 wrapper, which is already at ``pillbox_offset``."""
    return g4.PhysicalVolume([0, 0, 0], [0, 0, 0], pillbox_lv, "pillbox_pv", pillbox_foil_lv, reg)


def insert_vm2000(
    reg: g4.Registry,
    vm2000_material: g4.Material,
    surfaces: materials.surfaces.OpticalSurfaceRegistry,
    water_lv: g4.LogicalVolume,
    water_pv: g4.PhysicalVolume,
    manhole_rotation: list,
    man_hole_offset: float,
    cryo_displacement_z: float,
) -> tuple[g4.PhysicalVolume, ...]:
    # Pillbox: one VM2000 wrapper with the skirt placed inside it (see place_pillbox). The
    # manhole cut in the foil is one foil thickness smaller, so the foil wraps around the
    # opening's edge instead of sharing a face with the steel.
    pillbox_reflection_foil_tube = g4.solid.Tubs(
        "pillbox_reflection_foil_tube",
        shielding_foot_ir - reflective_foil_thickness,
        shielding_foot_or + reflective_foil_thickness,
        pillbox_foil_height,
        0,
        2 * pi,
        reg,
    )
    pillbox_reflection_foil = g4.solid.Subtraction(
        "pillbox_reflection_foil",
        pillbox_reflection_foil_tube,
        _manhole_solid("manhole_foil", manhole_outer_radius - reflective_foil_thickness, reg),
        [manhole_rotation, [0, 0, 0 - man_hole_offset]],
        reg,
    )
    pillbox_reflection_foil_lv = g4.LogicalVolume(
        pillbox_reflection_foil, vm2000_material, "pillbox_reflection_foil_lv", reg
    )
    pillbox_reflection_foil_lv.pygeom_color_rgba = COLORS["vm2000"]
    pillbox_reflection_foil_pv = g4.PhysicalVolume(
        [0, 0, 0],
        [0, 0, pillbox_offset],
        pillbox_reflection_foil_lv,
        "pillbox_reflection_foil_pv",
        water_lv,
        reg,
    )

    # VM2000 wrapped onto the cryostat: a copy of the cryostat grown by the foil thickness.
    # The cryostat itself is placed inside it (see core.py), so the VM2000 fills exactly the
    # 0.04 mm skin and there is no water in between.
    cryo_reflection_foil = cryo.cryostat_shell(
        "cryo_reflection_foil",
        reg,
        r=cryo.cryo_outer_radius + reflective_foil_thickness,
        h_cyl=cryo.cryo_outer_tub_height,
        h_top=cryo.cryo_outer_top_height + reflective_foil_thickness,
        h_bot=cryo.cryo_outer_bottom_height + reflective_foil_thickness,
        r_neck=cryo.cryo_access_outer_radius + reflective_foil_thickness,
        # the neck is cut off at the water surface: end the foil above the cryostat so that the
        # two do not share a face and the wrapper really does enclose the cryostat.
        z_neck=cryo.cryo_neck_top + reflective_foil_thickness,
    )
    cryo_reflection_foil_lv = g4.LogicalVolume(
        cryo_reflection_foil, vm2000_material, "cryo_reflection_foil_lv", reg
    )
    cryo_reflection_foil_lv.pygeom_color_rgba = COLORS["vm2000"]
    cryo_reflection_foil_pv = g4.PhysicalVolume(
        [0, 0, 0],
        [0, 0, cryo_displacement_z],
        cryo_reflection_foil_lv,
        "cryo_reflection_foil_pv",
        water_lv,
        reg,
    )

    def _vm2000_surfaces(pv: g4.PhysicalVolume) -> None:
        name = pv.name.replace("_pv", "")
        # The idea is to have the entire interface between water and VM2000 transparent.
        # VM2000 -> water
        g4.BorderSurface(name + "_water_border_surface", pv, water_pv, surfaces.vm2000_water_border, reg)
        # water -> VM2000
        g4.BorderSurface("water_" + name + "_border_surface", water_pv, pv, surfaces.vm2000_water_border, reg)

    _vm2000_surfaces(pillbox_reflection_foil_pv)
    _vm2000_surfaces(cryo_reflection_foil_pv)

    return (
        pillbox_reflection_foil_lv,
        pillbox_reflection_foil_pv,
        cryo_reflection_foil_lv,
        cryo_reflection_foil_pv,
    )


def insert_pmts(
    reg: g4.Registry,
    pmt_steel_material: g4.Material | str,
    cathode_al: g4.Material | str,
    pmt_air_material: g4.Material,
    surfaces: materials.surfaces.OpticalSurfaceRegistry,
    water_lv: g4.LogicalVolume,
    water_pv: g4.PhysicalVolume,
    acryl_material: g4.Material,
    borosilicate_material: g4.Material,
):
    # Photocathode and PMT encapsulation
    acryl = g4.solid.Sphere(
        "acryl",
        acryl_inner_radius,
        acryl_outer_radius,
        pmt_starting_angle,
        pmt_ending_angle,
        acryl_theta_start,
        acryl_theta_end,  # Is not a half circle, therefore gap i think
        reg,
    )
    pmt_air = g4.solid.Sphere(
        "pmt_air",
        photocathode_inner_radius,
        pmt_air_outer_radius,
        pmt_starting_angle,
        pmt_ending_angle,
        photocathode_theta_start,
        photocathode_theta_end,
        reg,
    )
    pmt_borosilikat_glass = g4.solid.Sphere(
        "pmt_borosilikat_glass",
        photocathode_inner_radius,
        pmt_borosilikat_glass_outer_radius,
        pmt_starting_angle,
        pmt_ending_angle,
        photocathode_theta_start,
        photocathode_theta_end,
        reg,
    )
    photocathode = g4.solid.Sphere(
        "photocathode",
        photocathode_inner_radius,
        photocathode_outer_radius,
        pmt_starting_angle,
        pmt_ending_angle,
        photocathode_theta_start,
        photocathode_theta_end,
        reg,
    )

    optical_steel_surface = surfaces.to_pmt_steel
    optical_pmt_surface = surfaces.to_photocathode

    # PMT encapsulation steel cone for Cherenkov veto
    pmt_steel_cone = g4.solid.Cons(
        "pmt_steel_cone",
        pmt_steel_cone_lower_rmin,
        pmt_steel_cone_lower_rmax,
        pmt_steel_cone_upper_rmin,
        pmt_steel_cone_upper_rmax,
        pmt_steel_cone_height * 0.5,
        pmt_starting_angle,
        pmt_ending_angle,
        reg,
    )
    pmt_steel_cone_lv = g4.LogicalVolume(pmt_steel_cone, pmt_steel_material, "pmt_steel_cone_lv", reg)
    pmt_steel_cone_lv.pygeom_color_rgba = COLORS["steel"]
    g4.SkinSurface("pmt_cone_optical_surface", pmt_steel_cone_lv, optical_steel_surface, reg)

    # PMT encapsulation bottom for Cherenkov veto
    pmt_steel_bottom = g4.solid.Tubs(
        "pmt_steel_bottom",
        0,
        pmt_steel_cone_lower_rmax,
        pmt_steel_bottom_height,
        pmt_starting_angle,
        pmt_ending_angle,
        reg,
    )
    pmt_steel_bottom_lv = g4.LogicalVolume(pmt_steel_bottom, pmt_steel_material, "pmt_steel_bottom_lv", reg)
    pmt_steel_bottom_lv.pygeom_color_rgba = COLORS["steel"]
    g4.SkinSurface("pmt_bottom_optical_surface", pmt_steel_bottom_lv, optical_steel_surface, reg)

    def build_pmt(
        broken: bool,
        index: int,
        xpos: float,
        ypos: float,
        zpos: float | None = None,
        xpos_cone: float | None = None,
        ypos_cone: float | None = None,  # The wall PMTs need extra stuff
        xpos_bottom: float | None = None,
        ypos_bottom: float | None = None,
        x_rot: float = 0.0,
        y_rot: float = 0.0,
        z_rot: float = 0.0,
    ) -> None:  # Rotations for wall PMTs
        prefix = "broken_" if broken else ""

        suffix = f"_{index}" if broken else f"_{pmt_id[index]}"

        name_pmt_acryl_lv = f"{prefix}pmt_acryl_lv{suffix}"
        name_pmt_acryl = f"{prefix}pmt_acryl{suffix}"

        name_pmt_air_lv = f"{prefix}pmt_air_lv{suffix}"
        name_pmt_air = f"{prefix}pmt_air{suffix}"

        name_pmt_borosilikat_lv = f"{prefix}pmt_borosilikat_lv{suffix}"
        name_pmt_borosilikat = f"{prefix}pmt_borosilikat{suffix}"

        namephotocathode_lv = f"{prefix}pmt_cathode_lv{suffix}"
        namephotocathode = f"{prefix}pmt_cathode{suffix}"

        namesteelcone = f"{prefix}pmt_cone{suffix}"
        namesteelbottom = f"{prefix}pmt_bottom{suffix}"

        # PMT first has acryl, then air, then borosilikat glass, then photocathode
        acryl_lv = g4.LogicalVolume(acryl, acryl_material, name_pmt_acryl_lv, reg)
        acryl_lv.pygeom_color_rgba = COLORS["acrylic"]
        # Place at different height depending if this is a wall pmt which specifies a zpos
        g4.PhysicalVolume(
            [x_rot, y_rot, z_rot],
            [xpos, ypos, pmt_cathode_offset if zpos is None else zpos],
            acryl_lv,
            name_pmt_acryl,
            water_lv,
            reg,
        )

        pmt_air_lv = g4.LogicalVolume(pmt_air, pmt_air_material, name_pmt_air_lv, reg)
        # the air would fill the PMT no matter how transparent the glass around it is
        pmt_air_lv.pygeom_color_rgba = False
        g4.PhysicalVolume([0, 0, 0], [0, 0, 0], pmt_air_lv, name_pmt_air, acryl_lv, reg)

        borosilikat_lv = g4.LogicalVolume(
            pmt_borosilikat_glass, borosilicate_material, name_pmt_borosilikat_lv, reg
        )
        borosilikat_lv.pygeom_color_rgba = COLORS["pmt_window"]
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            borosilikat_lv,
            name_pmt_borosilikat,
            pmt_air_lv,
            reg,
        )
        photocathode_lv = g4.LogicalVolume(photocathode, cathode_al, namephotocathode_lv, reg)
        photocathode_lv.pygeom_color_rgba = COLORS["pmt_cathode"]
        if broken:
            g4.PhysicalVolume([0, 0, 0], [0, 0, 0], photocathode_lv, namephotocathode, borosilikat_lv, reg)
        else:
            photocathode_pv = g4.PhysicalVolume(
                [0, 0, 0], [0, 0, 0], photocathode_lv, namephotocathode, borosilikat_lv, reg
            )
            photocathode_pv.set_pygeom_active_detector(RemageDetectorInfo("optical", pmt_rawids[index]))

        g4.SkinSurface(f"{prefix}pmt_cathode_skin_surface{suffix}", photocathode_lv, optical_pmt_surface, reg)

        # PMT steel cone and bottom
        # The wall PMTs should specify extra positions for cone and bottom
        g4.PhysicalVolume(
            [x_rot, y_rot, z_rot],
            [
                xpos if xpos_cone is None else xpos_cone,
                ypos if ypos_cone is None else ypos_cone,
                pmt_cone_offset if zpos is None else zpos,
            ],
            pmt_steel_cone_lv,
            namesteelcone,
            water_lv,
            reg,
        )
        g4.PhysicalVolume(
            [x_rot, y_rot, z_rot],
            [
                xpos if xpos_bottom is None else xpos_bottom,
                ypos if ypos_bottom is None else ypos_bottom,
                pmt_bottom_offset if zpos is None else zpos,
            ],
            pmt_steel_bottom_lv,
            namesteelbottom,
            water_lv,
            reg,
        )

    ############################################ PMT positions #########################################
    num_pmts = 0  # basically just an index
    working_pmts = 0  #  basically just an index

    # ------------------ Bottom outer ring
    n_pmt_per_ring_bo = 24  # spaces, in total 14 PMTs are actually placed
    r_pos = 4250.0
    dphi = 2.0 * np.pi / n_pmt_per_ring_bo

    empty_spaces = np.array([3, 5, 7, 9, 11, 15, 17, 19, 21, 23])  # These are empty spots without PMTs

    for k in range(n_pmt_per_ring_bo):
        if k in empty_spaces:
            continue
        dphi_c = dphi * k
        xpos = r_pos * np.cos(dphi_c)
        ypos = r_pos * np.sin(dphi_c)

        # all working
        build_pmt(False, working_pmts, xpos, ypos)
        working_pmts += 1
        num_pmts += 1

    # ------------------ Bottom inner ring
    n_pmt_per_ring_bi = 8  # 8 PMTs in inner ring
    r_pos = 2750.0  # radius of PMT positions
    dphi = 2.0 * np.pi / n_pmt_per_ring_bi

    broken_but_inside = np.array(
        [4, 7]
    )  # These pmts are inside but already broke during GERDA, so we have no DAQ ids for them

    # place PMTs
    for k in range(n_pmt_per_ring_bi):
        dphi_c = dphi * k
        xpos = r_pos * np.cos(dphi_c)
        ypos = r_pos * np.sin(dphi_c)

        build_pmt(k in broken_but_inside, num_pmts if k in broken_but_inside else working_pmts, xpos, ypos)
        working_pmts += 0 if k in broken_but_inside else 1
        num_pmts += 1

    # --------------- Pillbox Bottom ring
    n_pmt_per_ring_bpb = 6
    r_pos = (
        shielding_foot_ir
        - 0.5 * pmt_steel_cone_height
        - photocathode_height_difference
        - pmt_steel_bottom_height
        - reflective_foil_thickness
    )
    dphi = 2.0 * np.pi / n_pmt_per_ring_bpb

    empty_spaces = np.array([0, 3])  # These positions are empty

    # PMTs
    for k in range(n_pmt_per_ring_bpb):
        if k in empty_spaces:
            continue

        dphi_c = dphi * k
        xpos = r_pos * np.cos(dphi_c)
        ypos = r_pos * np.sin(dphi_c)

        build_pmt(False, working_pmts, xpos, ypos)
        working_pmts += 1
        num_pmts += 1

    # --------------------- Wall mounted PMT rings

    distancetobottom = 200.0  # distance of pillbox wall PMTs to bottom
    r_pos = (
        water_radius
        - reflective_foil_thickness
        - distance_pmt_base_tank
        - pmt_steel_cone_height * 0.5
        - pmt_steel_bottom_height
        - 2e-9
    )

    for i in [0, 1, 2, 3, 4]:
        if i == 0:
            # PMT row (2m)
            n_pmt_per_ring = 10
            dphi = 2.0 * np.pi / n_pmt_per_ring
            zpos = -2450.0

            empty_spaces = np.array([])
            broken_but_inside = np.array([4, 5, 6, 7])

        elif i == 1:
            # PMT row (3.5 m)
            n_pmt_per_ring = 10
            dphi = 2.0 * np.pi / n_pmt_per_ring
            zpos = -950.0

            empty_spaces = np.array([])
            broken_but_inside = np.array([4, 5])
        elif i == 2:
            # PMT row (5 m)
            n_pmt_per_ring = 10
            dphi = 2.0 * np.pi / n_pmt_per_ring
            zpos = 550.0

            empty_spaces = np.array([])
            broken_but_inside = np.array([0, 3])
        elif i == 3:
            # PMT row (6.5 m)
            n_pmt_per_ring = 10
            dphi = 2.0 * np.pi / n_pmt_per_ring
            zpos = 2050.0
            empty_spaces = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9])
            broken_but_inside = np.array([])
        elif i == 4:
            # Pillbox wall mounted PMTs
            empty_spaces = np.array([])
            broken_but_inside = np.array([])
            n_pmt_per_ring = 6
            r_pos = (
                shielding_foot_ir
                - reflective_foil_thickness
                - pmt_steel_cone_height * 0.5
                - pmt_steel_bottom_height
                - distance_pmt_base_pillbox
            )
            zpos = -(water_height / 2.0) + distancetobottom * 3
            dphi = 2 * np.pi / n_pmt_per_ring

        x_rot_drehvol = 0
        y_rot_drehvol = np.pi / 2.0

        # placing of PMTs in ring
        for k in range(n_pmt_per_ring):
            if k in empty_spaces:
                continue

            dphi_c = dphi * k + 0.5 * dphi if i in (0, 2, 4) else dphi * k

            xpos = (r_pos + photocathode_height_difference) * np.cos(dphi_c)
            ypos = (r_pos + photocathode_height_difference) * np.sin(dphi_c)

            if i in (0, 2, 4):
                z_rot_drehvol = -0.5 * dphi + n_pmt_per_ring * dphi - k * dphi
            else:
                z_rot_drehvol = n_pmt_per_ring * dphi - k * dphi

            # rotations around Z, Y and X
            rot_z = Rotation.from_euler("z", z_rot_drehvol, degrees=False)
            rot_y = Rotation.from_euler("y", y_rot_drehvol, degrees=False)
            rot_x = Rotation.from_euler("x", x_rot_drehvol, degrees=False)

            combined_rotation = rot_y * rot_z * rot_x

            with warnings.catch_warnings(action="ignore"):
                euler_angles = combined_rotation.as_euler("xyz", degrees=False)
            x_rot_global, y_rot_global, z_rot_global = euler_angles

            # PMT cone position
            xpos_cone = (r_pos + pmt_steel_cone_height * 0.25) * np.cos(dphi_c)
            ypos_cone = (r_pos + pmt_steel_cone_height * 0.25) * np.sin(dphi_c)

            # PMT bottom position
            xpos_bottom = (r_pos + pmt_steel_cone_height * 0.5 + pmt_steel_bottom_height * 0.5) * np.cos(
                dphi_c
            )
            ypos_bottom = (r_pos + pmt_steel_cone_height * 0.5 + pmt_steel_bottom_height * 0.5) * np.sin(
                dphi_c
            )

            build_pmt(
                k in broken_but_inside,
                index=num_pmts if k in broken_but_inside else working_pmts,
                xpos=xpos,
                ypos=ypos,
                zpos=zpos,
                xpos_cone=xpos_cone,
                ypos_cone=ypos_cone,
                xpos_bottom=xpos_bottom,
                ypos_bottom=ypos_bottom,
                x_rot=x_rot_global,
                y_rot=y_rot_global,
                z_rot=z_rot_global,
            )

            working_pmts += 0 if k in broken_but_inside else 1
            num_pmts += 1


def insert_muon_veto(
    reg: g4.Registry,
    world_lv: g4.LogicalVolume,
    mats: materials.OpticalMaterialRegistry,
    skip_pmts: bool = False,
) -> tuple[g4.LogicalVolume, g4.PhysicalVolume]:
    """Build the water tank, its VM2000 liner, the pillbox and the PMTs.

    Returns the cryostat's VM2000 wrapper, into which the cryostat itself is placed. The tank is
    offset by ``-cryo_z_displacement`` so that the cryostat ends up at the global origin.
    """
    water_tank_lv = construct_tank(reg, "G4_STAINLESS-STEEL")
    water_tank_pv = place_tank(reg, water_tank_lv, world_lv, -cryo_z_displacement)

    # tank steel -> VM2000 liner -> water
    water_tank_foil_lv = construct_tank_foil(reg, mats.vm2000)
    water_tank_foil_pv = place_tank_foil(reg, water_tank_foil_lv, water_tank_lv)

    water_lv = construct_water(reg, mats.water)
    water_pv = place_water(reg, water_lv, water_tank_foil_lv)

    air_buffer_lv = construct_air_buffer(reg, mats.nitrogen_air)
    place_air_buffer(reg, air_buffer_lv, water_lv)

    pillbox_lv, manhole_rotation, manhole_offset = construct_pillbox(reg, "G4_STAINLESS-STEEL")

    pillbox_foil_lv, pillbox_foil_pv, cryo_reflection_foil_lv, cryo_reflection_foil_pv = insert_vm2000(
        reg,
        mats.vm2000,
        mats.surfaces,
        water_lv,
        water_pv,
        manhole_rotation,
        manhole_offset,
        cryo_z_displacement,
    )
    pillbox_pv = place_pillbox(reg, pillbox_lv, pillbox_foil_lv)

    # The VM2000 <--> water boundary is transparent, so photons cross the foil and are then
    # reflected where the foil meets steel.
    # Optical path: Water --> VM2000 --> Water_tank (reflection at this surface)
    g4.BorderSurface(
        "water_tank_surface", water_tank_foil_pv, water_tank_pv, mats.surfaces.vm2000_reflective_border, reg
    )
    g4.BorderSurface(
        "water_pillbox_surface", pillbox_foil_pv, pillbox_pv, mats.surfaces.vm2000_reflective_border, reg
    )
    g4.BorderSurface(
        "water_tank_foil_water_border_surface",
        water_tank_foil_pv,
        water_pv,
        mats.surfaces.vm2000_water_border,
        reg,
    )
    g4.BorderSurface(
        "water_water_tank_foil_border_surface",
        water_pv,
        water_tank_foil_pv,
        mats.surfaces.vm2000_water_border,
        reg,
    )
    # The VM2000 -> Cryo surface has to be in core.py, where the cryostat is placed.

    if not skip_pmts:
        insert_pmts(
            reg,
            "G4_STAINLESS-STEEL",
            "G4_Al",
            mats.pmt_air,
            mats.surfaces,
            water_lv,
            water_pv,
            mats.acryl,
            mats.borosilicate,
        )

    return cryo_reflection_foil_lv, cryo_reflection_foil_pv
