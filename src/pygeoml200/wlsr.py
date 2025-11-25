"""Construct the LEGEND-200 Wavelength-Shifting Reflector (WLSR).

Dimensions mainly from [Krause2023]_ (general construction) and [Araujo2022]_ (TPB layer thickness).

.. [Krause2023] P. Krause "Shining Light on Backgrounds: An Advanced Liquid Argon Scintillation Light
    Detector for Boosting Background Suppression in LEGEND's Neutrinoless Double Beta Decay Search".
    PhD thesis (2023).
.. [Araujo2022] G. R. Araujo et al. “R&D of wavelength-shifting reflectors and characterization of
    the quantum eﬀiciency of tetraphenyl butadiene and polyethylene naphthalate in
    liquid argon.” In: The European Physical Journal C 82.5 (May 2022).
    https://doi.org/10.1140/epjc/s10052-022-10383-0
"""

from __future__ import annotations

from math import pi

import pyg4ometry.geant4 as g4

from . import core, materials

wlsr_tpb_radius = 1374 / 2
wlsr_ttx_thickness = 254 * 1e-3  # 254 um Tetratex foil
wlsr_cu_thickness = 50 * 1e-3  # 50 um copper foil for structure
wlsr_height = 3000

wlsr_tpb_thickness = 600 * 1e-6  # 600 nm of TPB coating


def _construct_wlsr(
    mats: materials.OpticalMaterialRegistry,
    reg: g4.Registry,
) -> tuple[g4.LogicalVolume, ...]:
    spacing = 2e-9  # spacing larger than the geometry tolerance (0.01 nm = 1e-11 m)
    wlsr_outer = g4.solid.Tubs(
        "wlsr_outer",
        wlsr_tpb_radius + wlsr_tpb_thickness + wlsr_ttx_thickness + 2 * spacing,
        wlsr_tpb_radius + wlsr_tpb_thickness + wlsr_ttx_thickness + 2 * spacing + wlsr_cu_thickness,
        wlsr_height,
        0,
        2 * pi,
        reg,
        "mm",
    )
    wlsr_ttx = g4.solid.Tubs(
        "wlsr_ttx",
        wlsr_tpb_radius + wlsr_tpb_thickness,
        wlsr_tpb_radius + wlsr_tpb_thickness + wlsr_ttx_thickness,
        wlsr_height - 2 * spacing,
        0,
        2 * pi,
        reg,
        "mm",
    )
    wlsr_tpb = g4.solid.Tubs(
        "wlsr_tpb",
        wlsr_tpb_radius,
        wlsr_tpb_radius + wlsr_tpb_thickness + wlsr_ttx_thickness + spacing,
        wlsr_height,
        0,
        2 * pi,
        reg,
        "mm",
    )

    wlsr_outer_lv = g4.LogicalVolume(wlsr_outer, mats.metal_copper, "wlsr_outer", reg)
    wlsr_ttx_lv = g4.LogicalVolume(wlsr_ttx, mats.tetratex, "wlsr_ttx", reg)
    wlsr_tpb_lv = g4.LogicalVolume(wlsr_tpb, mats.tpb_on_tetratex, "wlsr_tpb", reg)

    return wlsr_outer_lv, wlsr_ttx_lv, wlsr_tpb_lv


def place_wlsr(
    b: core.InstrumentationData,
    z_displacement: float,
    reg: g4.Registry,
) -> tuple[g4.PhysicalVolume]:
    wlsr_outer_lv, wlsr_ttx_lv, wlsr_tpb_lv = _construct_wlsr(b.materials, b.registry)

    wlsr_outer_pv = g4.PhysicalVolume(
        [0, 0, 0],
        [0, 0, z_displacement - wlsr_height / 2],
        wlsr_outer_lv,
        "wlsr_outer",
        b.mother_lv,
        reg,
    )
    wlsr_tpb_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], wlsr_tpb_lv, "wlsr_tpb", b.mother_lv, reg)
    wlsr_ttx_pv = g4.PhysicalVolume([0, 0, 0], [0, 0, 0], wlsr_ttx_lv, "wlsr_ttx", wlsr_tpb_lv, reg)

    wlsr_ttx_lv.pygeom_color_rgba = [1, 1, 1, 0.2]
    wlsr_tpb_lv.pygeom_color_rgba = False
    wlsr_outer_lv.pygeom_color_rgba = False

    _add_surfaces_wlsr(wlsr_ttx_pv, wlsr_tpb_pv, b)

    return wlsr_outer_pv, wlsr_ttx_pv, wlsr_tpb_pv


def _add_surfaces_wlsr(
    wlsr_ttx_pv: g4.PhysicalVolume,
    wlsr_tpb_pv: g4.PhysicalVolume,
    b: core.InstrumentationData,
):
    mats, reg = b.materials, b.registry
    # between TPB and TTX, only one surface should be enough.
    g4.BorderSurface("bsurface_tpb_ttx", wlsr_tpb_pv, wlsr_ttx_pv, mats.surfaces.to_tetratex, reg)

    # between LAr and TPB we need a surface in both directions.
    g4.BorderSurface("bsurface_wlsr_tpb_lar", b.mother_pv, wlsr_tpb_pv, mats.surfaces.lar_to_tpb, reg)
    g4.BorderSurface("bsurface_wlsr_lar_tpb", wlsr_tpb_pv, b.mother_pv, mats.surfaces.lar_to_tpb, reg)
