"""Subpackage to provide all implemented materials and their (optical) material properties."""

from __future__ import annotations

import legendoptics.fibers
import legendoptics.lar
import legendoptics.nylon
import legendoptics.pen
import legendoptics.pmts
import legendoptics.silica
import legendoptics.tpb
import legendoptics.ultem
import legendoptics.vm2000
import legendoptics.water
import numpy as np
import pint
import pyg4ometry.geant4 as g4
from pygeomtools.materials import BaseMaterialRegistry, cached_property

from .surfaces import OpticalSurfaceRegistry


class OpticalMaterialRegistry(BaseMaterialRegistry):
    def __init__(self, g4_registry: g4.Registry):
        super().__init__(g4_registry)

        self.lar_temperature = 88.8
        self.surfaces = OpticalSurfaceRegistry(g4_registry)

    @cached_property
    def liquidargon(self) -> g4.Material:
        """LEGEND liquid argon."""
        _liquidargon = g4.Material(
            name="liquid_argon",
            density=1.390,  # g/cm3
            number_of_components=1,
            state="liquid",
            temperature=self.lar_temperature,  # K
            pressure=1.0 * 1e5,  # pascal
            registry=self.g4_registry,
        )
        _liquidargon.add_element_natoms(self.get_element("Ar"), natoms=1)

        u = pint.get_application_registry().get()
        legendoptics.lar.pyg4_lar_attach_rindex(
            _liquidargon,
            self.g4_registry,
        )
        legendoptics.lar.pyg4_lar_attach_attenuation(
            _liquidargon,
            self.g4_registry,
            self.lar_temperature * u.K,
        )
        legendoptics.lar.pyg4_lar_attach_scintillation(
            _liquidargon,
            self.g4_registry,
            triplet_lifetime_method="legend200-llama",
        )

        return _liquidargon

    @cached_property
    def metal_steel(self) -> g4.Material:
        """Stainless steel of the GERDA cryostat."""
        _metal_steel = g4.Material(
            name="metal_steel",
            density=7.9,
            number_of_components=5,
            registry=self.g4_registry,
        )
        _metal_steel.add_element_massfraction(self.get_element("Si"), massfraction=0.01)
        _metal_steel.add_element_massfraction(self.get_element("Cr"), massfraction=0.20)
        _metal_steel.add_element_massfraction(self.get_element("Mn"), massfraction=0.02)
        _metal_steel.add_element_massfraction(self.get_element("Fe"), massfraction=0.67)
        _metal_steel.add_element_massfraction(self.get_element("Ni"), massfraction=0.10)

        return _metal_steel

    @cached_property
    def metal_silicon(self) -> g4.Material:
        """Silicon."""
        _metal_silicon = g4.Material(
            name="metal_silicon",
            density=2.330,
            number_of_components=1,
            registry=self.g4_registry,
        )
        _metal_silicon.add_element_natoms(self.get_element("Si"), natoms=1)

        return _metal_silicon

    @cached_property
    def metal_tantalum(self) -> g4.Material:
        """Tantalum."""
        _metal_tantalum = g4.Material(
            name="metal_tantalum",
            density=16.69,
            number_of_components=1,
            registry=self.g4_registry,
        )
        _metal_tantalum.add_element_natoms(self.get_element("Ta"), natoms=1)

        return _metal_tantalum

    @cached_property
    def metal_copper(self) -> g4.Material:
        """Copper structures.

        .. warning:: For full optics support, a reflective surface is needed, see
            :py:func:`surfaces.OpticalSurfaceRegistry.to_copper`.
        """
        _metal_copper = g4.Material(
            name="metal_copper",
            density=8.960,
            number_of_components=1,
            registry=self.g4_registry,
        )
        _metal_copper.add_element_natoms(self.get_element("Cu"), natoms=1)

        return _metal_copper

    @cached_property
    def metal_phosphor_bronze(self) -> g4.Material:
        """Copper structures."""
        _metal_phbr = g4.Material(
            name="metal_phosphor_bronze",
            density=8.960,
            number_of_components=1,
            registry=self.g4_registry,
        )
        _metal_phbr.add_element_massfraction(self.get_element("Cu"), massfraction=0.948)
        _metal_phbr.add_element_massfraction(self.get_element("Sn"), massfraction=0.05)
        _metal_phbr.add_element_massfraction(self.get_element("P"), massfraction=0.002)

        return _metal_phbr

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

    @cached_property
    def peek(self) -> g4.Material:
        """PEEK for the SIS absorber holder."""
        _peek = g4.Material(name="peek", density=1.320, number_of_components=3, registry=self.g4_registry)
        _peek.add_element_natoms(self.get_element("C"), natoms=19)
        _peek.add_element_natoms(self.get_element("H"), natoms=12)  # TODO: MaGe uses C19H1203??
        _peek.add_element_natoms(self.get_element("O"), natoms=3)

        return _peek

    @cached_property
    def pmma(self) -> g4.Material:
        """PMMA for the inner fiber cladding layer."""
        _pmma = g4.Material(name="pmma", density=1.2, number_of_components=3, registry=self.g4_registry)
        _pmma.add_element_natoms(self.get_element("H"), natoms=8)
        _pmma.add_element_natoms(self.get_element("C"), natoms=5)
        _pmma.add_element_natoms(self.get_element("O"), natoms=2)

        legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(_pmma, self.g4_registry)

        return _pmma

    @cached_property
    def pmma_out(self) -> g4.Material:
        """PMMA for the outer fiber cladding layer."""
        _pmma_out = g4.Material(
            name="pmma_cl2",
            density=1.2,
            number_of_components=3,
            registry=self.g4_registry,
        )
        _pmma_out.add_element_natoms(self.get_element("H"), natoms=8)
        _pmma_out.add_element_natoms(self.get_element("C"), natoms=5)
        _pmma_out.add_element_natoms(self.get_element("O"), natoms=2)

        legendoptics.fibers.pyg4_fiber_cladding2_attach_rindex(_pmma_out, self.g4_registry)

        return _pmma_out

    @cached_property
    def ps_fibers(self) -> g4.Material:
        """Polystyrene for the fiber core."""
        _ps_fibers = g4.Material(
            name="ps_fibers",
            density=1.05,
            number_of_components=2,
            registry=self.g4_registry,
        )
        _ps_fibers.add_element_natoms(self.get_element("H"), natoms=8)
        _ps_fibers.add_element_natoms(self.get_element("C"), natoms=8)

        legendoptics.fibers.pyg4_fiber_core_attach_rindex(_ps_fibers, self.g4_registry)
        legendoptics.fibers.pyg4_fiber_core_attach_absorption(_ps_fibers, self.g4_registry)
        legendoptics.fibers.pyg4_fiber_core_attach_wls(_ps_fibers, self.g4_registry)
        legendoptics.fibers.pyg4_fiber_core_attach_scintillation(_ps_fibers, self.g4_registry)

        return _ps_fibers

    def _tpb(self, name: str, **wls_opts) -> g4.Material:
        t = g4.Material(
            name=name,
            density=1.08,
            number_of_components=2,
            state="solid",
            registry=self.g4_registry,
        )
        t.add_element_natoms(self.get_element("H"), natoms=22)
        t.add_element_natoms(self.get_element("C"), natoms=28)

        legendoptics.tpb.pyg4_tpb_attach_rindex(t, self.g4_registry)
        legendoptics.tpb.pyg4_tpb_attach_wls(t, self.g4_registry, **wls_opts)

        return t

    @cached_property
    def tpb_on_fibers(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (evaporated on fibers)."""
        return self._tpb("tpb_on_fibers")

    @cached_property
    def tpb_on_tetratex(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (evaporated on Tetratex)."""
        return self._tpb("tpb_on_tetratex")

    @cached_property
    def tpb_on_nylon(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (in nylon matrix)."""
        # as a base, use the normal TPB properties.
        _tpb_on_nylon = self._tpb(
            "tpb_on_nylon",
            # For 30% TPB 70% PS the WLS light yield is reduced by 30% [Alexey]
            quantum_efficiency=0.7 * legendoptics.tpb.tpb_quantum_efficiency(),
            # the emission spectrum differs significantly.
            emission_spectrum="polystyrene_matrix",
        )

        # add absorption length from nylon.
        legendoptics.nylon.pyg4_nylon_attach_absorption(_tpb_on_nylon, self.g4_registry)

        return _tpb_on_nylon

    @cached_property
    def tetratex(self) -> g4.Material:
        """Tetratex diffuse reflector.

        .. warning:: For full optics support, a reflective surface is needed, see
            :py:func:`surfaces.OpticalSurfaceRegistry.wlsr_tpb_to_tetratex`.
        """
        _tetratex = g4.Material(
            name="tetratex",
            density=0.35,
            number_of_components=2,
            registry=self.g4_registry,
        )
        _tetratex.add_element_massfraction(self.get_element("F"), massfraction=0.76)
        _tetratex.add_element_massfraction(self.get_element("C"), massfraction=0.24)

        return _tetratex

    @cached_property
    def nylon(self) -> g4.Material:
        """Nylon (from Borexino)."""
        _nylon = g4.Material(
            name="nylon",
            density=1.15,
            number_of_components=4,
            registry=self.g4_registry,
        )
        _nylon.add_element_natoms(self.get_element("H"), natoms=2)
        _nylon.add_element_natoms(self.get_element("N"), natoms=2)
        _nylon.add_element_natoms(self.get_element("O"), natoms=3)
        _nylon.add_element_natoms(self.get_element("C"), natoms=13)

        legendoptics.nylon.pyg4_nylon_attach_rindex(_nylon, self.g4_registry)
        legendoptics.nylon.pyg4_nylon_attach_absorption(_nylon, self.g4_registry)

        return _nylon

    @cached_property
    def pen(self) -> g4.Material:
        """PEN wavelength-shifter and scintillator."""
        _pen = g4.Material(
            name="pen",
            density=1.3,
            number_of_components=3,
            registry=self.g4_registry,
        )
        _pen.add_element_natoms(self.get_element("C"), natoms=14)
        _pen.add_element_natoms(self.get_element("H"), natoms=10)
        _pen.add_element_natoms(self.get_element("O"), natoms=4)

        legendoptics.pen.pyg4_pen_attach_rindex(_pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_attenuation(_pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_wls(_pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_scintillation(_pen, self.g4_registry)

        return _pen

    @cached_property
    def water(self) -> g4.Material:
        """High purity water of the watertank."""
        _water = g4.MaterialCompound(
            name="Water",  # written "Water" to use Geant4 intern way of handling Rayleigh scattering with water,
            # see Geant4 BookForApplicationDevelopers pg. 270
            density=1.0,
            number_of_components=2,
            registry=self.g4_registry,
        )

        _water.add_element_natoms(self.get_element("H"), natoms=2)
        _water.add_element_natoms(self.get_element("O"), natoms=1)

        legendoptics.water.pyg4_water_attach_rindex(_water, self.g4_registry)
        legendoptics.water.pyg4_water_attach_absorption(_water, self.g4_registry)

        return _water

    @cached_property
    def vm2000(self) -> g4.Material:
        """Material for the reflective foil VM2000 based on nylon (e.g. MaGe)."""
        _vm2000 = g4.MaterialCompound(
            name="vm2000",
            density=1.15,
            number_of_components=4,
            registry=self.g4_registry,
        )

        # Add elements with their mass fractions
        _vm2000.add_element_natoms(self.get_element("H"), natoms=2)
        _vm2000.add_element_natoms(self.get_element("N"), natoms=2)
        _vm2000.add_element_natoms(self.get_element("O"), natoms=3)
        _vm2000.add_element_natoms(self.get_element("C"), natoms=13)

        legendoptics.vm2000.pyg4_vm2000_attach_absorption_length(_vm2000, self.g4_registry)
        legendoptics.vm2000.pyg4_vm2000_attach_rindex(_vm2000, self.g4_registry)
        legendoptics.vm2000.pyg4_vm2000_attach_wls(_vm2000, self.g4_registry)
        # VM2000 seem to consist of PMMA and PEN layers https://iopscience.iop.org/article/10.1088/1748-0221/12/06/P06017/pdf
        legendoptics.pen.pyg4_pen_attach_scintillation(_vm2000, self.g4_registry)
        legendoptics.vm2000.pyg4_vm2000_attach_particle_scintillationyields(_vm2000, self.g4_registry)

        return _vm2000

    @cached_property
    def acryl(self) -> g4.Material:
        """Material for the acryl cap of the PMT encapsulation."""
        _acryl = g4.MaterialCompound(
            name="acryl",
            density=1.18,
            number_of_components=2,
            registry=self.g4_registry,
        )

        _acryl.add_element_natoms(self.get_element("H"), natoms=2)
        _acryl.add_element_natoms(self.get_element("C"), natoms=1)

        legendoptics.pmts.pyg4_pmt_attach_acryl_rindex(_acryl, self.g4_registry)
        legendoptics.pmts.pyg4_pmt_attach_acryl_absorption_length(_acryl, self.g4_registry)

        return _acryl

    @cached_property
    def pmt_air(self) -> g4.Material:
        """Material for the air in between Acryl cap and PMT."""
        _pmt_air = g4.MaterialCompound(
            name="PMT_air",
            density=0.001225,
            number_of_components=2,
            registry=self.g4_registry,
        )

        _pmt_air.add_element_natoms(self.get_element("N"), natoms=3)
        _pmt_air.add_element_natoms(self.get_element("O"), natoms=1)

        legendoptics.pmts.pyg4_pmt_attach_air_rindex(_pmt_air, self.g4_registry)
        legendoptics.pmts.pyg4_pmt_attach_air_absorption_length(_pmt_air, self.g4_registry)

        return _pmt_air

    @cached_property
    def borosilicate(self) -> g4.Material:
        """Material for the borosilicate glass of the PMT."""
        _borosilicate = g4.MaterialCompound(
            name="borosilicate",
            density=2.23,
            number_of_components=4,
            registry=self.g4_registry,
        )

        _borosilicate.add_element_massfraction(self.get_element("Si"), 0.376)
        _borosilicate.add_element_massfraction(self.get_element("O"), 0.543)
        _borosilicate.add_element_massfraction(self.get_element("B"), 0.04)
        _borosilicate.add_element_massfraction(self.get_element("Na"), 0.029)
        _borosilicate.add_element_massfraction(self.get_element("Al"), 0.012)

        legendoptics.pmts.pyg4_pmt_attach_borosilicate_rindex(_borosilicate, self.g4_registry)
        legendoptics.pmts.pyg4_pmt_attach_borosilicate_absorption_length(_borosilicate, self.g4_registry)

        return _borosilicate

    @cached_property
    def ultem(self) -> g4.Material:
        """Ultem for the receptacles and insulators."""
        _ultem = g4.Material(
            name="ultem",
            density=1.27,
            number_of_components=4,
            registry=self.g4_registry,
        )
        _ultem.add_element_natoms(self.get_element("C"), natoms=37)
        _ultem.add_element_natoms(self.get_element("H"), natoms=24)
        _ultem.add_element_natoms(self.get_element("O"), natoms=6)
        _ultem.add_element_natoms(self.get_element("N"), natoms=2)

        legendoptics.ultem.pyg4_ultem_attach_rindex(_ultem, self.g4_registry)
        legendoptics.ultem.pyg4_ultem_attach_absorption(_ultem, self.g4_registry)

        return _ultem

    @cached_property
    def silica(self) -> g4.Material:
        """Silica for the LMFE."""
        _silica = g4.Material(
            name="silica",
            density=2.2,
            number_of_components=2,
            registry=self.g4_registry,
        )
        _silica.add_element_natoms(self.get_element("Si"), natoms=1)
        _silica.add_element_natoms(self.get_element("O"), natoms=2)

        legendoptics.silica.pyg4_silica_attach_rindex(_silica, self.g4_registry)

        return _silica
