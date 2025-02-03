"""Subpackage to provide all implemented materials and their (optical) material properties."""

from __future__ import annotations

import legendoptics.fibers
import legendoptics.lar
import legendoptics.nylon
import legendoptics.pen
import legendoptics.tpb
import numpy as np
import pint
import pyg4ometry.gdml.Defines as defines
import pyg4ometry.geant4 as g4

from . import VM2000
from .surfaces import OpticalSurfaceRegistry


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
            name="Water",
            density=1.0,
            number_of_components=2,
            registry=self.g4_registry,
        )

        self._water.add_element_natoms(self.get_element("H"), natoms=2)
        self._water.add_element_natoms(self.get_element("O"), natoms=1)

        # add refraction index
        photon_energy = np.array([1.0, 6.0])
        refractive_index = [1.33, 1.33]
        refractive_matrix = defines.MatrixFromVectors(
            photon_energy, refractive_index, "RefractiveMatrix", self.g4_registry, "eV", ""
        )

        # add attenuation length
        # Photon energy values corresponding to the wavelengths (in eV)
        photon_energy_water = np.array(
            [
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
        )

        # Corresponding attenuation lengths (in mm)
        absorption_lengths = np.array(
            [
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
            ]
        )

        attenuation_length_matrix = defines.MatrixFromVectors(
            photon_energy_water, absorption_lengths, "AbsorptionLengthMatrix", self.g4_registry, "eV", "mm"
        )

        self._water.addProperty("ABSLENGTH", attenuation_length_matrix)
        self._water.addProperty("RINDEX", refractive_matrix)

        return self._water

    @property
    def nylon_VM2000(self) -> g4.Material:
        """Material for the reflective foil VM2000."""
        if hasattr(self, "_nylon"):
            return self._nylon

        self._nylon = g4.MaterialCompound(
            name="nylon",
            density=1.15,
            number_of_components=4,
            registry=self.g4_registry,
        )

        # Add elements with their mass fractions
        self._nylon.add_element_natoms(self.get_element("H"), natoms=2)
        self._nylon.add_element_natoms(self.get_element("N"), natoms=2)
        self._nylon.add_element_natoms(self.get_element("O"), natoms=3)
        self._nylon.add_element_natoms(self.get_element("C"), natoms=13)

        num1 = 251
        Refraction = np.ones(num1) * 1.15  # Estimated refractive index
        AbsorptionL = np.ones(num1) * 50.0  # 50 m
        Refraction[0] = Refraction[1]
        AbsorptionL[0] = AbsorptionL[1]

        params = VM2000.VM2000_parameters()
        VM2000_energy_range, WLS_absorption, WLS_emission = params[0], params[3], params[4]

        # Make matrices for the properties
        rindex_nylon = defines.MatrixFromVectors(
            VM2000_energy_range, Refraction, "RindexNylon", self.g4_registry, "eV", ""
        )
        abslength_nylon = defines.MatrixFromVectors(
            VM2000_energy_range, AbsorptionL, "AbslengthNylon", self.g4_registry, "eV", ""
        )
        wlsabslength_nylon = defines.MatrixFromVectors(
            VM2000_energy_range, WLS_absorption, "WLSAbslengthNylon", self.g4_registry, "eV", ""
        )
        wlscomponent_nylon = defines.MatrixFromVectors(
            VM2000_energy_range, WLS_emission, "WLSComponentNylon", self.g4_registry, "eV", ""
        )

        self._nylon.addProperty("RINDEX", rindex_nylon)
        self._nylon.addProperty("ABSLENGTH", abslength_nylon)
        self._nylon.addProperty("WLSABSLENGTH", wlsabslength_nylon)
        self._nylon.addProperty("WLSCOMPONENT", wlscomponent_nylon)
        legendoptics.pen.pyg4_pen_attach_scintillation(self._nylon, self.g4_registry)
        self._nylon.addConstProperty("WLSTIMECONSTANT", 0.5 * 10e-3)  # ns

        return self._nylon

    @property
    def PMT_air(self) -> g4.Material:
        """Material for the air in between Acryl cap and PMT."""
        if hasattr(self, "_PMT_air"):
            return self._PMT_air

        self._PMT_air = g4.MaterialCompound(
            name="PMT_air",
            density=0.001225,
            number_of_components=2,
            registry=self.g4_registry,
        )

        self._PMT_air.add_element_natoms(self.get_element("N"), natoms=3)
        self._PMT_air.add_element_natoms(self.get_element("O"), natoms=1)

        photon_energy_air = np.array([1.0, 6.0])
        refractive_index_air = [1.0, 1.0]
        refractive_matrix_air = defines.MatrixFromVectors(
            photon_energy_air, refractive_index_air, "RefractiveMatrixAir", self.g4_registry, "eV", ""
        )
        absorption_length_air = [100000.0, 100000.0]
        absorption_matrix_air = defines.MatrixFromVectors(
            photon_energy_air, absorption_length_air, "AbsorptionMatrixAir", self.g4_registry, "eV", "mm"
        )
        self._PMT_air.addProperty("ABSLENGTH", absorption_matrix_air)
        self._PMT_air.addProperty("RINDEX", refractive_matrix_air)

        return self._PMT_air

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
        refractive_matrix_acryl = defines.MatrixFromVectors(
            photon_energy_acryl, refractive_index_acryl, "RefractiveMatrixAcryl", self.g4_registry, "eV", ""
        )
        absorption_length_acryl = [2500.0, 3500.0]  # in mm, 3,5 m for all energies
        absorption_matrix_acryl = defines.MatrixFromVectors(
            photon_energy_acryl,
            absorption_length_acryl,
            "AbsorptionMatrixAcryl",
            self.g4_registry,
            "eV",
            "mm",
        )
        self._acryl.addProperty("ABSLENGTH", absorption_matrix_acryl)
        self._acryl.addProperty("RINDEX", refractive_matrix_acryl)

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

        # SiO2: Silicon and Oxygen
        self._borosilicate.add_element_massfraction(self.get_element("Si"), 0.806 * 28.09 / (28.09 + 2 * 16))
        self._borosilicate.add_element_massfraction(self.get_element("O"), 0.806 * 2 * 16 / (28.09 + 2 * 16))

        # B2O3: Boron and Oxygen
        self._borosilicate.add_element_massfraction(
            self.get_element("B"), 0.130 * 2 * 10.81 / (2 * 10.81 + 3 * 16)
        )
        self._borosilicate.add_element_massfraction(
            self.get_element("O"), 0.130 * 3 * 16 / (2 * 10.81 + 3 * 16)
        )

        # Na2O: Sodium and Oxygen
        self._borosilicate.add_element_massfraction(
            self.get_element("Na"), 0.040 * 2 * 22.99 / (2 * 22.99 + 16)
        )
        self._borosilicate.add_element_massfraction(self.get_element("O"), 0.040 * 16 / (2 * 22.99 + 16))

        # Al2O3: Aluminum and Oxygen
        self._borosilicate.add_element_massfraction(
            self.get_element("Al"), 0.023 * 2 * 26.98 / (2 * 26.98 + 3 * 16)
        )
        self._borosilicate.add_element_massfraction(
            self.get_element("O"), 0.023 * 3 * 16 / (2 * 26.98 + 3 * 16)
        )

        photon_energy_cathode = np.array([1.0, 6.0])
        refractive_index_cathode = [1.49, 1.49]
        refractive_matrix_cathode = defines.MatrixFromVectors(
            photon_energy_cathode,
            refractive_index_cathode,
            "RefractiveMatrixCathode",
            self.g4_registry,
            "eV",
            "",
        )

        absorption_length_cathode = [2000.0, 3000.0]
        absorption_matrix_cathode = defines.MatrixFromVectors(
            photon_energy_cathode,
            absorption_length_cathode,
            "AbsorptionMatrixCathode",
            self.g4_registry,
            "eV",
            "mm",
        )

        self._borosilicate.addProperty("RINDEX", refractive_matrix_cathode)
        self._borosilicate.addProperty("ABSLENGTH", absorption_matrix_cathode)

        return self._borosilicate
