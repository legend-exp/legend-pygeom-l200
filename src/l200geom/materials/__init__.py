"""Subpackage to provide all implemented materials and their (optical) material properties."""

from __future__ import annotations

import legendoptics.fibers
import legendoptics.lar
import legendoptics.nylon
import legendoptics.pen
import legendoptics.pyg4utils
import legendoptics.tpb
import numpy as np
import pint
import pyg4ometry.geant4 as g4

from . import vm2000
from .surfaces import OpticalSurfaceRegistry

u = pint.get_application_registry()


class OpticalMaterialRegistry:
    def __init__(self, g4_registry: g4.Registry):
        self.g4_registry = g4_registry
        self.lar_temperature = 88.8

        self._elements = {}
        self._elements_cb = {}
        self._define_elements()

        self.surfaces = OpticalSurfaceRegistry(g4_registry)

    def get_element(self, symbol: str) -> g4.Element:
        if (symbol in self._elements_cb) and (symbol not in self._elements):
            self._elements[symbol] = (self._elements_cb[symbol])()
        return self._elements[symbol]

    def _add_element(self, name: str, symbol: str, z: int, a: float) -> None:
        """Lazily define an element on the current registry."""
        assert symbol not in self._elements_cb
        self._elements_cb[symbol] = lambda: g4.ElementSimple(
            name=name, symbol=symbol, Z=z, A=a, registry=self.g4_registry
        )

    def _define_elements(self) -> None:
        """Lazily define all used elements."""
        self._add_element(name="Hydrogen", symbol="H", z=1, a=1.00794)
        self._add_element(name="Boron", symbol="B", z=5, a=10.811)
        self._add_element(name="Carbon", symbol="C", z=6, a=12.011)
        self._add_element(name="Nitrogen", symbol="N", z=7, a=14.01)
        self._add_element(name="Oxygen", symbol="O", z=8, a=16.00)
        self._add_element(name="Fluorine", symbol="F", z=9, a=19.00)
        self._add_element(name="Sodium", symbol="Na", z=11, a=22.99)
        self._add_element(name="Aluminium", symbol="Al", z=13, a=26.981539)
        self._add_element(name="Silicon", symbol="Si", z=14, a=28.09)
        self._add_element(name="argon", symbol="Ar", z=18, a=39.95)
        self._add_element(name="Chromium", symbol="Cr", z=24, a=51.9961)
        self._add_element(name="Manganese", symbol="Mn", z=25, a=54.93805)
        self._add_element(name="Iron", symbol="Fe", z=26, a=55.845)
        self._add_element(name="Indium", symbol="In", z=49, a=114.82)
        self._add_element(name="Cobalt", symbol="Co", z=27, a=58.9332)
        self._add_element(name="Nickel", symbol="Ni", z=28, a=58.6934)
        self._add_element(name="Copper", symbol="Cu", z=29, a=63.55)
        self._add_element(name="Tantalum", symbol="Ta", z=73, a=180.94)
        self._add_element(name="Gold", symbol="Au", z=79, a=196.967)

    @property
    def liquidargon(self) -> g4.Material:
        """LEGEND liquid argon."""
        if hasattr(self, "_liquidargon"):
            return self._liquidargon

        self._liquidargon = g4.Material(
            name="liquid_argon",
            density=1.390,  # g/cm3
            number_of_components=1,
            state="liquid",
            temperature=self.lar_temperature,  # K
            pressure=1.0 * 1e5,  # pascal
            registry=self.g4_registry,
        )
        self._liquidargon.add_element_natoms(self.get_element("Ar"), natoms=1)

        u = pint.get_application_registry().get()
        legendoptics.lar.pyg4_lar_attach_rindex(
            self._liquidargon,
            self.g4_registry,
        )
        legendoptics.lar.pyg4_lar_attach_attenuation(
            self._liquidargon,
            self.g4_registry,
            self.lar_temperature * u.K,
        )
        legendoptics.lar.pyg4_lar_attach_scintillation(
            self._liquidargon,
            self.g4_registry,
            triplet_lifetime_method="legend200-llama",
        )

        return self._liquidargon

    @property
    def metal_steel(self) -> g4.Material:
        """Stainless steel of the GERDA cryostat."""
        if hasattr(self, "_metal_steel"):
            return self._metal_steel

        self._metal_steel = g4.Material(
            name="metal_steel",
            density=7.9,
            number_of_components=5,
            registry=self.g4_registry,
        )
        self._metal_steel.add_element_massfraction(self.get_element("Si"), massfraction=0.01)
        self._metal_steel.add_element_massfraction(self.get_element("Cr"), massfraction=0.20)
        self._metal_steel.add_element_massfraction(self.get_element("Mn"), massfraction=0.02)
        self._metal_steel.add_element_massfraction(self.get_element("Fe"), massfraction=0.67)
        self._metal_steel.add_element_massfraction(self.get_element("Ni"), massfraction=0.10)

        return self._metal_steel

    @property
    def metal_silicon(self) -> g4.Material:
        """Silicon."""
        if hasattr(self, "_metal_silicon"):
            return self._metal_silicon

        self._metal_silicon = g4.Material(
            name="metal_silicon",
            density=2.330,
            number_of_components=1,
            registry=self.g4_registry,
        )
        self._metal_silicon.add_element_natoms(self.get_element("Si"), natoms=1)

        return self._metal_silicon

    @property
    def metal_tantalum(self) -> g4.Material:
        """Tantalum."""
        if hasattr(self, "_metal_tantalum"):
            return self._metal_tantalum

        self._metal_tantalum = g4.Material(
            name="metal_tantalum",
            density=16.69,
            number_of_components=1,
            registry=self.g4_registry,
        )
        self._metal_tantalum.add_element_natoms(self.get_element("Ta"), natoms=1)

        return self._metal_tantalum

    @property
    def metal_copper(self) -> g4.Material:
        """Copper structures.

        .. warning:: For full optics support, a reflective surface is needed, see
            :py:func:`surfaces.OpticalSurfaceRegistry.to_copper`.
        """
        if hasattr(self, "_metal_copper"):
            return self._metal_copper

        self._metal_copper = g4.Material(
            name="metal_copper",
            density=8.960,
            number_of_components=1,
            registry=self.g4_registry,
        )
        self._metal_copper.add_element_natoms(self.get_element("Cu"), natoms=1)

        return self._metal_copper

    @property
    def metal_caps_gold(self) -> g4.Material:
        """Gold for calibration source.

        .. note:: modified density in order to have the equivalent of two 2x2cm gold
            foils, with 20 um thickness.
        """
        if hasattr(self, "_metal_caps_gold"):
            return self._metal_caps_gold

        # quoting https://doi.org/10.1088/1748-0221/18/02/P02001:
        # After the deposition, the external part of the foil with no 228Th
        # activity was cut off, and the foil rolled

        volume_of_foil = np.pi * (1 / 8 * 2.54) ** 2 * 50e-4  # 1/4” diameter, 50 um thickness
        volume_of_inner = np.pi * 0.2**2 * 0.4  # 2 cm radius, 4 cm height
        self._metal_caps_gold = g4.Material(
            name="metal_caps_gold",
            density=19.3 * volume_of_foil / volume_of_inner,
            number_of_components=1,
            registry=self.g4_registry,
        )
        self._metal_caps_gold.add_element_natoms(self.get_element("Au"), natoms=1)

        return self._metal_caps_gold

    @property
    def peek(self) -> g4.Material:
        """PEEK for the SIS absorber holder."""
        if hasattr(self, "_peek"):
            return self._peek

        self._peek = g4.Material(
            name="peek", density=1.320, number_of_components=3, registry=self.g4_registry
        )
        self._peek.add_element_natoms(self.get_element("C"), natoms=19)
        self._peek.add_element_natoms(self.get_element("H"), natoms=12)  # TODO: MaGe uses C19H1203??
        self._peek.add_element_natoms(self.get_element("O"), natoms=3)

        return self._peek

    @property
    def pmma(self) -> g4.Material:
        """PMMA for the inner fiber cladding layer."""
        if hasattr(self, "_pmma"):
            return self._pmma

        self._pmma = g4.Material(name="pmma", density=1.2, number_of_components=3, registry=self.g4_registry)
        self._pmma.add_element_natoms(self.get_element("H"), natoms=8)
        self._pmma.add_element_natoms(self.get_element("C"), natoms=5)
        self._pmma.add_element_natoms(self.get_element("O"), natoms=2)

        legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(
            self._pmma,
            self.g4_registry,
        )

        return self._pmma

    @property
    def pmma_out(self) -> g4.Material:
        """PMMA for the outer fiber cladding layer."""
        if hasattr(self, "_pmma_out"):
            return self._pmma_out

        self._pmma_out = g4.Material(
            name="pmma_cl2",
            density=1.2,
            number_of_components=3,
            registry=self.g4_registry,
        )
        self._pmma_out.add_element_natoms(self.get_element("H"), natoms=8)
        self._pmma_out.add_element_natoms(self.get_element("C"), natoms=5)
        self._pmma_out.add_element_natoms(self.get_element("O"), natoms=2)

        legendoptics.fibers.pyg4_fiber_cladding2_attach_rindex(
            self._pmma_out,
            self.g4_registry,
        )

        return self._pmma_out

    @property
    def ps_fibers(self) -> g4.Material:
        """Polystyrene for the fiber core."""
        if hasattr(self, "_ps_fibers"):
            return self._ps_fibers

        self._ps_fibers = g4.Material(
            name="ps_fibers",
            density=1.05,
            number_of_components=2,
            registry=self.g4_registry,
        )
        self._ps_fibers.add_element_natoms(self.get_element("H"), natoms=8)
        self._ps_fibers.add_element_natoms(self.get_element("C"), natoms=8)

        legendoptics.fibers.pyg4_fiber_core_attach_rindex(
            self._ps_fibers,
            self.g4_registry,
        )
        legendoptics.fibers.pyg4_fiber_core_attach_absorption(
            self._ps_fibers,
            self.g4_registry,
        )
        legendoptics.fibers.pyg4_fiber_core_attach_wls(
            self._ps_fibers,
            self.g4_registry,
        )

        return self._ps_fibers

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

    @property
    def tpb_on_fibers(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (evaporated on fibers)."""
        if hasattr(self, "_tpb_on_fibers"):
            return self._tpb_on_fibers

        self._tpb_on_fibers = self._tpb("tpb_on_fibers")

        return self._tpb_on_fibers

    @property
    def tpb_on_tetratex(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (evaporated on Tetratex)."""
        if hasattr(self, "_tpb_on_tetratex"):
            return self._tpb_on_tetratex

        self._tpb_on_tetratex = self._tpb("tpb_on_tetratex")

        return self._tpb_on_tetratex

    @property
    def tpb_on_nylon(self) -> g4.Material:
        """Tetraphenyl-butadiene wavelength shifter (in nylon matrix)."""
        if hasattr(self, "_tpb_on_nylon"):
            return self._tpb_on_nylon

        # as a base, use the normal TPB properties.
        self._tpb_on_nylon = self._tpb(
            "tpb_on_nylon",
            # For 30% TPB 70% PS the WLS light yield is reduced by 30% [Alexey]
            quantum_efficiency=0.7 * legendoptics.tpb.tpb_quantum_efficiency(),
            # the emission spectrum differs significantly.
            emission_spectrum="polystyrene_matrix",
        )

        # add absorption length from nylon.
        legendoptics.nylon.pyg4_nylon_attach_absorption(self._tpb_on_nylon, self.g4_registry)

        return self._tpb_on_nylon

    @property
    def tetratex(self) -> g4.Material:
        """Tetratex diffuse reflector.

        .. warning:: For full optics support, a reflective surface is needed, see
            :py:func:`surfaces.OpticalSurfaceRegistry.wlsr_tpb_to_tetratex`.
        """
        if hasattr(self, "_tetratex"):
            return self._tetratex

        self._tetratex = g4.Material(
            name="tetratex",
            density=0.35,
            number_of_components=2,
            registry=self.g4_registry,
        )
        self._tetratex.add_element_massfraction(self.get_element("F"), massfraction=0.76)
        self._tetratex.add_element_massfraction(self.get_element("C"), massfraction=0.24)

        return self._tetratex

    @property
    def nylon(self) -> g4.Material:
        """Nylon (from Borexino)."""
        if hasattr(self, "_nylon"):
            return self._nylon

        self._nylon = g4.Material(
            name="nylon",
            density=1.15,
            number_of_components=4,
            registry=self.g4_registry,
        )
        self._nylon.add_element_natoms(self.get_element("H"), natoms=2)
        self._nylon.add_element_natoms(self.get_element("N"), natoms=2)
        self._nylon.add_element_natoms(self.get_element("O"), natoms=3)
        self._nylon.add_element_natoms(self.get_element("C"), natoms=13)

        legendoptics.nylon.pyg4_nylon_attach_rindex(self._nylon, self.g4_registry)
        legendoptics.nylon.pyg4_nylon_attach_absorption(self._nylon, self.g4_registry)

        return self._nylon

    @property
    def pen(self) -> g4.Material:
        """PEN wavelength-shifter and scintillator."""
        if hasattr(self, "_pen"):
            return self._pen

        self._pen = g4.Material(
            name="pen",
            density=1.3,
            number_of_components=3,
            registry=self.g4_registry,
        )
        self._pen.add_element_natoms(self.get_element("C"), natoms=14)
        self._pen.add_element_natoms(self.get_element("H"), natoms=10)
        self._pen.add_element_natoms(self.get_element("O"), natoms=4)

        legendoptics.pen.pyg4_pen_attach_rindex(self._pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_attenuation(self._pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_wls(self._pen, self.g4_registry)
        legendoptics.pen.pyg4_pen_attach_scintillation(self._pen, self.g4_registry)

        return self._pen

    @property
    def water(self) -> g4.Material:
        """High purity water of the watertank."""
        if hasattr(self, "_water"):
            return self._water

        self._water = g4.MaterialCompound(
            name="Water",  # written "Water" to use Geant4 intern way of handling Rayleigh scattering with water,
            # see Geant4 BookForApplicationDevelopers pg. 270
            density=1.0,
            number_of_components=2,
            registry=self.g4_registry,
        )

        self._water.add_element_natoms(self.get_element("H"), natoms=2)
        self._water.add_element_natoms(self.get_element("O"), natoms=1)

        # add refraction index
        photon_energy = [1.0, 6.0]
        refractive_index = [1.33, 1.33]

        # add attenuation length
        # Photon energy absorption corresponding to the wavelengths
        photon_energy_water = [
            1.239841939 / 0.6,  # ~206.6 nm
            1.239841939 / 0.55,  # ~224.5 nm
            1.239841939 / 0.50,  # ~248.0 nm
            1.239841939 / 0.45,  # ~275.5 nm
            1.239841939 / 0.40,  # ~310 nm
            1.239841939 / 0.35,  # ~354.0 nm
            1.239841939 / 0.30,  # ~413.3 nm
            1.239841939 / 0.25,  # ~496.0 nm
            1.239841939 / 0.20,  # ~620 nm
            1.239841939 / 0.19,  # ~652.6 nm
            1.239841939 / 0.10,  # ~1240 nm
        ]

        # Corresponding attenuation lengths (in mm)
        absorption_lengths = [
            10 * 1000,  # 10 m for 206.6 nm
            20 * 1000,  # 20 m for 224.5 nm
            50 * 1000,  # 50 m for 248.0 nm
            100 * 1000,  # 100 m for 275.5 nm
            100 * 1000,  # 100 m for 310 nm
            100 * 1000,  # 100 m for 354 nm
            90 * 1000,  # 90 m for 413.3 nm
            20 * 1000,  # 20 m for 496.0 nm
            1 * 1000,  # 1 m for 620 nm
            0.001,  # 0.001 mm for 652.6 nm
            0.0001,  # 0.0001 mm for 1240 nm
        ] * u.m

        with u.context("sp"):
            self._water.addVecPropertyPint("ABSLENGTH", photon_energy_water, absorption_lengths)
            self._water.addVecPropertyPint("RINDEX", photon_energy, refractive_index)

        return self._water

    @property
    def vm2000(self) -> g4.Material:
        """Material for the reflective foil VM2000 based on nylon (e.g. MaGe)."""
        if hasattr(self, "_vm2000"):
            return self._vm2000

        self._vm2000 = g4.MaterialCompound(
            name="vm2000",
            density=1.15,
            number_of_components=4,
            registry=self.g4_registry,
        )

        # Add elements with their mass fractions
        self._vm2000.add_element_natoms(self.get_element("H"), natoms=2)
        self._vm2000.add_element_natoms(self.get_element("N"), natoms=2)
        self._vm2000.add_element_natoms(self.get_element("O"), natoms=3)
        self._vm2000.add_element_natoms(self.get_element("C"), natoms=13)

        vm2000_energy_range, _, _, wls_absorption, wls_emission = vm2000.vm2000_parameters()

        refraction = np.ones_like(vm2000_energy_range) * 1.15  # Estimated refractive index
        absorptionl = np.ones_like(vm2000_energy_range) * 50.0 * u.m

        with u.context("sp"):
            self._vm2000.addVecPropertyPint("RINDEX", vm2000_energy_range, refraction)
            self._vm2000.addVecPropertyPint("ABSLENGTH", vm2000_energy_range, absorptionl)
            self._vm2000.addVecPropertyPint("WLSABSLENGTH", vm2000_energy_range, wls_absorption)
            self._vm2000.addVecPropertyPint("WLSCOMPONENT", vm2000_energy_range, wls_emission)

        # VM2000 seem to consist of PMMA and PEN layers https://iopscience.iop.org/article/10.1088/1748-0221/12/06/P06017/pdf
        legendoptics.pen.pyg4_pen_attach_scintillation(self._vm2000, self.g4_registry)
        self._vm2000.addConstProperty("WLSTIMECONSTANT", 0.5 * 10e-3)  # ns

        return self._vm2000

    @property
    def pmt_air(self) -> g4.Material:
        """Material for the air in between Acryl cap and PMT."""
        if hasattr(self, "_pmt_air"):
            return self._pmt_air

        self._pmt_air = g4.MaterialCompound(
            name="PMT_air",
            density=0.001225,
            number_of_components=2,
            registry=self.g4_registry,
        )

        self._pmt_air.add_element_natoms(self.get_element("N"), natoms=3)
        self._pmt_air.add_element_natoms(self.get_element("O"), natoms=1)

        photon_energy_air = [1.0, 6.0]
        refractive_index_air = [1.0, 1.0]
        absorption_length_air = [100.0, 100.0] * u.m

        self._pmt_air.addVecPropertyPint("RINDEX", photon_energy_air, refractive_index_air)
        self._pmt_air.addVecPropertyPint("ABSLENGTH", photon_energy_air, absorption_length_air)

        return self._pmt_air

    @property
    def acryl(self) -> g4.Material:
        """Material for the acryl cap of the PMT encapsulation."""
        if hasattr(self, "_acryl"):
            return self._acryl

        self._acryl = g4.MaterialCompound(
            name="acryl",
            density=1.18,
            number_of_components=2,
            registry=self.g4_registry,
        )

        self._acryl.add_element_natoms(self.get_element("H"), natoms=2)
        self._acryl.add_element_natoms(self.get_element("C"), natoms=1)

        photon_energy_acryl = np.array([1.0, 6.0])
        refractive_index_acryl = [1.489, 1.489]
        absorption_length_acryl = [2.5, 3.5] * u.m  # 2,5 m up to 3,5 m

        self._acryl.addVecPropertyPint("RINDEX", photon_energy_acryl, refractive_index_acryl)
        self._acryl.addVecPropertyPint("ABSLENGTH", photon_energy_acryl, absorption_length_acryl)

        return self._acryl

    @property
    def borosilicate(self) -> g4.Material:
        """Material for the borosilicate glass of the PMT."""
        if hasattr(self, "_borosilicate"):
            return self._borosilicate

        self._borosilicate = g4.MaterialCompound(
            name="borosilicate",
            density=2.23,
            number_of_components=4,
            registry=self.g4_registry,
        )

        self._borosilicate.add_element_massfraction(self.get_element("Si"), 0.376)
        self._borosilicate.add_element_massfraction(self.get_element("O"), 0.543)
        self._borosilicate.add_element_massfraction(self.get_element("B"), 0.04)
        self._borosilicate.add_element_massfraction(self.get_element("Na"), 0.029)
        self._borosilicate.add_element_massfraction(self.get_element("Al"), 0.012)

        photon_energy_cathode = np.array([1.0, 6.0])
        refractive_index_cathode = [1.49, 1.49]
        absorption_length_cathode = [2.0, 3.0] * u.m

        self._borosilicate.addVecPropertyPint("RINDEX", photon_energy_cathode, refractive_index_cathode)
        self._borosilicate.addVecPropertyPint("ABSLENGTH", photon_energy_cathode, absorption_length_cathode)

        return self._borosilicate
