"""Subpackage to provide all implemented materials and their (optical) material properties."""

from __future__ import annotations

import numpy as np
import pyg4ometry.geant4 as g4
from pygeomtools.materials import LegendMaterialRegistry, cached_property

from .surfaces import OpticalSurfaceRegistry


class OpticalMaterialRegistry(LegendMaterialRegistry):
    def __init__(self, g4_registry: g4.Registry):
        super().__init__(g4_registry)

        self.surfaces = OpticalSurfaceRegistry(g4_registry)

    @cached_property
    def metal_caps_gold(self) -> g4.Material:
        """Gold for the Th228 calibration source described in https://doi.org/10.1088/1748-0221/18/02/P02001.

        .. note:: modified density in order to have the equivalent of the gold foils inside the source.
        """
        from ..calibration import source_th_height_inner, source_th_radius_inner

        # quoting https://doi.org/10.1088/1748-0221/18/02/P02001:
        # After the deposition, the external part of the foil with no 228Th activity was cut off, and the
        # foil rolled.
        # from private communication with Ralph, this does not mean the round section, but the larger
        # quadratic foil area shown in figure 2 in the paper. The inner source dimensions are guessed from
        # photos (i.e. figure 2 on the right).

        # 1/2” diameter (measured from figure 2), 50 um thickness
        volume_of_foil = (0.5 * 2.54) ** 2 * 50e-4  # cm^3

        # volume of the implemented source region
        volume_of_inner = np.pi * (source_th_radius_inner * 0.1) ** 2 * source_th_height_inner * 0.1  # cm^3

        # scale down density of the gold block to have the same number of gold atoms.
        density = 19.3 * volume_of_foil / volume_of_inner

        _metal_caps_gold = g4.Material(
            name="metal_caps_gold",
            density=density,
            number_of_components=1,
            registry=self.g4_registry,
        )
        _metal_caps_gold.add_element_natoms(self.get_element("Au"), natoms=1)

        return _metal_caps_gold
