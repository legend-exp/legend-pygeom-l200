"""Construct the LEGEND-200 Wavelength-Shifting Reflector (WLSR).

Dimensions mainly from P. Krause.
"""

from __future__ import annotations

from math import pi

import pyg4ometry.geant4 as g4

from . import core


def construct_wlsr(
    structure_material: g4.Material,
    tetratex_material: g4.Material,
    tpb_material: g4.Material,
    reg: g4.Registry,
) -> tuple[g4.LogicalVolume]:
    wlsr_outer_diameter = 1400 / 2  # (Patrick)
    wlsr_tpb_diameter = 1374 / 2
    wlsr_tpb_thickness = 600 * 1e-6  # 600nm of TPB coating (arXiv:2112.06675)
    wlsr_thickness = (wlsr_outer_diameter - wlsr_tpb_diameter - wlsr_tpb_thickness) / 2  # dummy
    wlsr_height = 3000
    wlsr_outer = g4.solid.Tubs(
        "wlsr_outer",
        wlsr_tpb_diameter,
        wlsr_outer_diameter,
        wlsr_height,
        0,
        2 * pi,
        reg,
        "mm",
    )
    wlsr_ttx = g4.solid.Tubs(
        "wlsr_ttx",
        wlsr_tpb_diameter,
        wlsr_outer_diameter - wlsr_thickness,
        wlsr_height,
        0,
        2 * pi,
        reg,
        "mm",
    )
    wlsr_tpb = g4.solid.Tubs(
        "wlsr_tpb",
        wlsr_tpb_diameter,
        wlsr_outer_diameter - wlsr_thickness - wlsr_tpb_thickness,
        wlsr_height,
        0,
        2 * pi,
        reg,
        "mm",
    )

    wlsr_outer_lv = g4.LogicalVolume(wlsr_outer, structure_material, "wlsr_outer", reg)
    wlsr_ttx_lv = g4.LogicalVolume(wlsr_ttx, tetratex_material, "wlsr_ttx", reg)
    wlsr_tpb_lv = g4.LogicalVolume(wlsr_tpb, tpb_material, "wlsr_tpb", reg)

    return wlsr_outer_lv, wlsr_ttx_lv, wlsr_tpb_lv


def place_wlsr(
    wlsr_outer_lv: g4.LogicalVolume,
    wlsr_ttx_lv: g4.LogicalVolume,
    wlsr_tpb_lv: g4.LogicalVolume,
    mother_lv: g4.LogicalVolume,
    z_displacement: float,
    reg: g4.Registry,
) -> tuple[g4.PhysicalVolume]:
    wlsr_outer_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, z_displacement], wlsr_outer_lv, "wlsr_outer", mother_lv, reg
    )
    wlsr_ttx_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], wlsr_ttx_lv, "wlsr_ttx", wlsr_outer_lv, reg)
    wlsr_tpb_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], wlsr_tpb_lv, "wlsr_tpb", wlsr_ttx_lv, reg)

    wlsr_ttx_lv.pygeom_color_rgba = [1, 1, 1, 0.5]
    wlsr_tpb_lv.pygeom_color_rgba = False
    wlsr_outer_lv.pygeom_color_rgba = False

    return wlsr_outer_pv, wlsr_ttx_pv, wlsr_tpb_pv


def add_surfaces_wlsr(
    wlsr_ttx_pv: g4.PhysicalVolume,
    wlsr_tpb_pv: g4.PhysicalVolume,
    b: core.InstrumentationData,
):
    mats, reg = b.materials, b.registry
    # between TPB and TTX, only one surface should be enough.
    g4.BorderSurface("bsurface_tpb_ttx", wlsr_tpb_pv, wlsr_ttx_pv, mats.surfaces.wlsr_tpb_to_tetratex, reg)

    # between LAr and TPB we need a surface in both directions.
    g4.BorderSurface("bsurface_wlsr_tpb_lar", b.mother_pv, wlsr_tpb_pv, mats.surfaces.lar_to_tpb, reg)
    g4.BorderSurface("bsurface_wlsr_lar_tpb", wlsr_tpb_pv, b.mother_pv, mats.surfaces.lar_to_tpb, reg)
