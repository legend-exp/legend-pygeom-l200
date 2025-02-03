"""Subpackage to provide all implemented optical surfaces and their properties."""

from __future__ import annotations

from pathlib import Path

import legendoptics.copper
import legendoptics.germanium
import legendoptics.silicon
import legendoptics.tetratex
import numpy as np
import pint
import pyg4ometry.gdml.Defines as defines
import pyg4ometry.geant4 as g4

from . import VM2000
from .ketek_sipm import ketek_sipm_efficiency

u = pint.get_application_registry()


class OpticalSurfaceRegistry:
    """Register and define optical surfaces.

    Note on Models
    --------------

    * UNIFIED model:
        `value` is the `sigma_alpha` parameter, the stddev of the newly chosen facet normal direction.
        For details on this model and its parameters, see `UNIFIED model diagram`_.
    * GLISUR model:
        `value` as smoothness, in range [0,1] (0=rough, 1=perfectly smooth).

    UNIFIED is more comprehensive, but is not directly equivalent to GLISUR. One notable difference is
    that UNIFIED/ground surfaces w/o specular probabilities set will not perform total internal reflection
    according to alpha1=alpha2, whereas GFLISUR/ground will do! Polished surfaces should behave similar
    between UNIFIED and GLISUR.

    .. _UNIFIED model diagram: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/_images/UNIFIED_model_diagram.png
    """

    def __init__(self, reg: g4.Registry):
        self.g4_registry = reg
        # do not change the surface model w/o also changing all surface values below!
        self._model = "unified"

    @property
    def to_copper(self) -> g4.solid.OpticalSurface:
        """Reflective surface for copper structure."""
        if hasattr(self, "_to_copper"):
            return self._to_copper

        self._to_copper = g4.solid.OpticalSurface(
            "surface_to_copper",
            finish="ground",
            model=self._model,
            surf_type="dielectric_metal",
            value=0.5,  # rad. converted from 0.5, probably a GLISUR smoothness parameter, in MaGe.
            registry=self.g4_registry,
        )

        legendoptics.copper.pyg4_copper_attach_reflectivity(
            self._to_copper,
            self.g4_registry,
        )

        return self._to_copper

    @property
    def to_germanium(self) -> g4.solid.OpticalSurface:
        """Reflective surface for germanium detectors."""
        if hasattr(self, "_to_germanium"):
            return self._to_germanium

        self._to_germanium = g4.solid.OpticalSurface(
            "surface_to_germanium",
            finish="ground",
            model=self._model,
            surf_type="dielectric_metal",
            value=0.3,  # rad. converted from 0.5, probably a GLISUR smoothness parameter, in MaGe.
            registry=self.g4_registry,
        )

        legendoptics.germanium.pyg4_germanium_attach_reflectivity(
            self._to_germanium,
            self.g4_registry,
        )

        return self._to_germanium

    @property
    def to_tetratex(self) -> g4.solid.OpticalSurface:
        """Reflective surface Tetratex diffuse reflector."""
        if hasattr(self, "_to_tetratex"):
            return self._to_tetratex

        self._to_tetratex = g4.solid.OpticalSurface(
            "surface_to_tetratex",
            finish="groundfrontpainted",  # only lambertian reflection
            model=self._model,
            surf_type="dielectric_dielectric",
            value=0,  # rad. perfectly lambertian reflector.
            registry=self.g4_registry,
        )

        legendoptics.tetratex.pyg4_tetratex_attach_reflectivity(
            self._to_tetratex,
            self.g4_registry,
        )

        return self._to_tetratex

    @property
    def to_sipm_silicon(self) -> g4.solid.OpticalSurface:
        """Reflective surface for KETEK SiPM."""
        if hasattr(self, "_to_sipm_silicon"):
            return self._to_sipm_silicon

        self._to_sipm_silicon = g4.solid.OpticalSurface(
            "surface_to_sipm_silicon",
            finish="ground",
            model=self._model,
            surf_type="dielectric_metal",
            value=0.05,  # converted from 0.9, probably a GLISUR smoothness parameter, in MaGe.
            registry=self.g4_registry,
        )

        legendoptics.silicon.pyg4_silicon_attach_complex_rindex(
            self._to_sipm_silicon,
            self.g4_registry,
        )

        # add custom efficiency for the KETEK SiPMs. This is not part of legendoptics.
        λ, eff = ketek_sipm_efficiency()
        with u.context("sp"):
            self._to_sipm_silicon.addVecPropertyPint("EFFICIENCY", λ.to("eV"), eff)

        return self._to_sipm_silicon

    @property
    def lar_to_tpb(self) -> g4.solid.OpticalSurface:
        """Optical surface between LAr and TBP wavelength shifting coating."""
        if hasattr(self, "_lar_to_tpb"):
            return self._lar_to_tpb

        self._lar_to_tpb = g4.solid.OpticalSurface(
            "surface_lar_to_tpb",
            finish="ground",
            model="unified",
            surf_type="dielectric_dielectric",
            value=0.3,  # rad. converted from 0.5, probably a GLISUR smoothness parameter, in MaGe.
            registry=self.g4_registry,
        )

        return self._lar_to_tpb

    @property
    def lar_to_pen(self) -> g4.solid.OpticalSurface:
        """Optical surface between LAr and PEN scintillator/wavelength shifting coating."""
        if hasattr(self, "_lar_to_pen"):
            return self._lar_to_pen

        self._lar_to_pen = g4.solid.OpticalSurface(
            "surface_lar_to_pen",
            finish="ground",
            model="unified",
            surf_type="dielectric_dielectric",
            # sigma_alpha corresponds to the value from L. Manzanillas et al 2022 JINST 17 P09007.
            value=0.01,  # rad.
            registry=self.g4_registry,
        )

        # set specular spike/lobe probabilities. From a presentation by L. Manzanillas ("Update on PEN optical
        # parameters for Geant4 studies", 29.04.2021); slide 12 "Recommended data for PEN simulations", both are 0.5.
        # note: MaGe has specularlobe=0.4, specularspike=0.6; but commented out. Luis' last code uses the same values
        # as MaGe:
        # https://github.com/lmanzanillas/AttenuationPenSetup/blob/11ff9664e3b2da3c0ebf726b60ad96111e9b2aaa/src/DetectorConstruction.cc#L1771-L1786
        λ = np.array([650.0, 115.0]) * u.nm
        specular_lobe = np.array([0.4, 0.4])
        specular_spike = np.array([0.6, 0.6])
        with u.context("sp"):
            self._lar_to_pen.addVecPropertyPint("SPECULARSPIKECONSTANT", λ.to("eV"), specular_spike)
            self._lar_to_pen.addVecPropertyPint("SPECULARLOBECONSTANT", λ.to("eV"), specular_lobe)

        return self._lar_to_pen

    @property
    def to_VM2000(self) -> g4.solid.OpticalSurface:
        """Reflective surface for VM2000."""
        if hasattr(self, "_to_VM2000"):
            return self._to_VM2000

        # Create material properties table for VM2000 surface
        self._to_VM2000 = g4.solid.OpticalSurface(
            name="WaterTankFoilSurface",
            finish="polished",
            model="glisur",
            surf_type="dielectric_metal",
            value=0.01,
            registry=self.g4_registry,
        )

        params = VM2000.VM2000_parameters()
        VM2000_energy_range, VM2000_Reflectivity, VM2000_Efficiency = params[0], params[1], params[2]

        # Make matrices for the properties
        reflectivity_matrix_foil = defines.MatrixFromVectors(
            VM2000_energy_range, VM2000_Reflectivity, "ReflectivityMatrixFoil", self.g4_registry, "eV", ""
        )
        efficiency_matrix_foil = defines.MatrixFromVectors(
            VM2000_energy_range, VM2000_Efficiency, "EfficiencyMatrixFoil", self.g4_registry, "eV", ""
        )
        # Add the matrices to the optical properties
        self._to_VM2000.addProperty(
            "REFLECTIVITY", reflectivity_matrix_foil
        )  # der VM2000 zugewandte Seite --> reflectivity=0
        self._to_VM2000.addProperty("EFFICIENCY", efficiency_matrix_foil)

        return self._to_VM2000

    @property
    def water_to_VM2000(self) -> g4.solid.OpticalSurface:
        """Optical surface between water and VM2000."""
        if hasattr(self, "_water_to_VM2000"):
            return self._water_to_VM2000

        # Create material properties table for VM2000 border surface
        self._water_to_VM2000 = g4.solid.OpticalSurface(
            name="WaterTankFoilBorder",
            finish="polished",
            model="glisur",
            surf_type="dielectric_metal",
            value=0.01,
            registry=self.g4_registry,
        )

        params = VM2000.VM2000_parameters()
        VM2000_energy_range, VM2000_Reflectivity, VM2000_Efficiency = params[0], params[1], params[2]

        Reflectivity_front = VM2000_Reflectivity * 0
        Efficiency_Border = VM2000_Efficiency * 0
        Transmittance_Border = [1.0] * len(VM2000_energy_range)
        reflectivity_matrix_front = defines.MatrixFromVectors(
            VM2000_energy_range, Reflectivity_front, "ReflectivityMatrixFoilFront", self.g4_registry, "eV", ""
        )

        # Make matrices for the properties
        efficiency_matrix_border_foil = defines.MatrixFromVectors(
            VM2000_energy_range, Efficiency_Border, "EfficiencyMatrixBorderFoil", self.g4_registry, "eV", ""
        )
        transmittance_matrix_border_foil = defines.MatrixFromVectors(
            VM2000_energy_range,
            Transmittance_Border,
            "TransmittanceMatrixBorderFoil",
            self.g4_registry,
            "eV",
            "",
        )
        # Add the matrices to the optical properties
        self._water_to_VM2000.addProperty("REFLECTIVITY", reflectivity_matrix_front)  # --> facing to VM2000
        self._water_to_VM2000.addProperty("EFFICIENCY", efficiency_matrix_border_foil)
        self._water_to_VM2000.addProperty("TRANSMITTANCE", transmittance_matrix_border_foil)

        return self._water_to_VM2000

    @property
    def to_steel(self) -> g4.solid.OpticalSurface:
        """Optical surface of steel."""
        if hasattr(self, "_to_steel"):
            return self._to_steel

        self._to_steel = g4.solid.OpticalSurface(
            name="PMTSteelSurface",
            finish="polished",
            model="glisur",
            surf_type="dielectric_metal",
            value=0.9,
            registry=self.g4_registry,
        )

        photon_energy = np.array([1.0, 6.0])
        reflectivity_steel = [0.9, 0.9]
        efficiency_steel = np.array([1.0, 1.0])

        reflectivity_steel_matrix = defines.MatrixFromVectors(
            photon_energy, reflectivity_steel, "ReflectiveSteelMatrix", self.g4_registry, "eV", ""
        )
        efficiency_steel_matrix = defines.MatrixFromVectors(
            photon_energy, efficiency_steel, "EfficiencySteelMatrix", self.g4_registry, "eV", ""
        )

        self._to_steel.addProperty("REFLECTIVITY", reflectivity_steel_matrix)
        self._to_steel.addProperty("EFFICIENCY", efficiency_steel_matrix)

        return self._to_steel

    @property
    def to_photocathode(self) -> g4.solid.OpticalSurface:
        """Optical surface of the PMT photocathode."""
        if hasattr(self, "_to_photocathode"):
            return self._to_photocathode

        path = Path(__file__).resolve().parent

        file_path = path / "PMT_QE.csv"

        data = np.loadtxt(file_path, delimiter=",")  # load data

        # Split the data into two arrays: wavelengths and efficiencies
        wavelengths = data[:, 0]  # First column: wavelengths
        PMT_quantum_efficiencies = data[:, 1]  # Define wavelength range (nm) and corresponding efficiencies
        PMT_quantum_efficiencies = PMT_quantum_efficiencies * 0.01  # in percent

        # Convert wavelengths to energy using hc/λ (in eV)
        h_planck = 4.1357e-15  # eV·s
        c_speed = 299792458  # m/s
        photon_energy = h_planck * c_speed / (wavelengths * 1e-9)

        # Detector Surface
        self._to_photocathode = g4.solid.OpticalSurface(
            name="PMTCathodeSurface",
            finish="polished",  # Finish of surface
            model="glisur",  # Model of surface
            surf_type="dielectric_metal",  # Type of surface
            value=0.9,  # parameter ((max. absorbance?)
            registry=self.g4_registry,
        )

        collection_efficiency = 0.85
        # Make matrices for the properties
        efficiency_matrix = defines.MatrixFromVectors(
            photon_energy,
            PMT_quantum_efficiencies * collection_efficiency,
            "EfficiencyMatrix",
            self.g4_registry,
            "eV",
            "",
        )

        # Add the matrices to the optical properties
        self._to_photocathode.addProperty("EFFICIENCY", efficiency_matrix)

        return self._to_photocathode

    @property
    def acryl_to_air(self) -> g4.solid.OpticalSurface:
        """Optical surface between acryl and air."""
        if hasattr(self, "_acryl_to_air"):
            return self._acryl_to_air

        self._acryl_to_air = g4.solid.OpticalSurface(
            name="PMTAirSurface",
            finish="polished",  # Finish of surface
            model="glisur",  # Model of surface
            surf_type="dielectric_dielectric",  # Type of surface
            value=1.0,  # parameter ((max. absorbance?)
            registry=self.g4_registry,
        )

        photon_energy_air = np.array([1.0, 6.0])
        reflectivity_air = [0.0386, 0.0386]
        transmittance_air = [1.0 - 0.0386, 1.0 - 0.0386]

        reflectivity_matrix_air = defines.MatrixFromVectors(
            photon_energy_air, reflectivity_air, "ReflectivityMatrixAir", self.g4_registry, "eV", ""
        )
        transmittance_matrix_air = defines.MatrixFromVectors(
            photon_energy_air, transmittance_air, "TransmittanceMatrixAir", self.g4_registry, "eV", ""
        )

        self._acryl_to_air.addProperty("REFLECTIVITY", reflectivity_matrix_air)
        self._acryl_to_air.addProperty("TRANSMITTANCE", transmittance_matrix_air)

        return self._acryl_to_air

    @property
    def water_to_acryl(self) -> g4.solid.OpticalSurface:
        """Optical surface between water and acryl."""
        if hasattr(self, "_water_to_acryl"):
            return self._water_to_acryl

        self._water_to_acryl = g4.solid.OpticalSurface(
            name="PMTAcrylSurface",
            finish="polished",  # Finish of surface
            model="glisur",  # Model of surface
            surf_type="dielectric_dielectric",  # Type of surface
            value=0.01,  # parameter ((max. absorbance?) --> bei glisur 0.0 vollständig diffuse 1.0 vollständig spekulare reflektion
            registry=self.g4_registry,
        )

        photon_energy_acryl = np.array([1.0, 6.0])
        reflectivity_acryl = [0.00318, 0.00318]
        transmittance_acryl = [1.0 - 0.00318, 1.0 - 0.00318]

        reflectivity_matrix_acryl = defines.MatrixFromVectors(
            photon_energy_acryl, reflectivity_acryl, "ReflectivityMatrixAcryl", self.g4_registry, "eV", ""
        )
        efficiency_matrix_acryl = defines.MatrixFromVectors(
            photon_energy_acryl, transmittance_acryl, "TransmittanceMatrixAcryl", self.g4_registry, "eV", ""
        )

        self._water_to_acryl.addProperty("REFLECTIVITY", reflectivity_matrix_acryl)
        self._water_to_acryl.addProperty("TRANSMITTANCE", efficiency_matrix_acryl)

        return self._water_to_acryl

    @property
    def air_to_borosilicate(self) -> g4.solid.OpticalSurface:
        """Optical surface between water and acryl."""
        if hasattr(self, "_air_to_borosilicate"):
            return self._air_to_borosilicate

        self._air_to_borosilicate = g4.solid.OpticalSurface(
            name="PMTBorosilikatSurface",
            finish="polished",  # Finish of surface
            model="glisur",  # Model of surface
            surf_type="dielectric_dielectric",  # Type of surface
            value=0.01,  # parameter ((max. absorbance?) --> bei glisur 0.0 vollständig diffuse 1.0 vollständig spekulare reflektion
            registry=self.g4_registry,
        )

        # Optional: Reflexions- und Transmissionsparameter (nicht immer nötig für 'dielectric_dielectric')
        photon_energy_borosilicate = np.array([1.0, 6.0])  # Beispiel-Energiebereich in eV
        transmittance_borosilicate = [1.0 - 0.036, 1.0 - 0.036]  # 100% Transmission

        path = Path(__file__).resolve().parent

        file_path = path / "PMT_QE.csv"

        data = np.loadtxt(file_path, delimiter=",")

        # Split the data into two arrays: wavelengths and efficiencies
        wavelengths = data[:, 0]  # First column: wavelengths

        # Convert wavelengths to energy using hc/λ (in eV)
        h_planck = 4.1357e-15  # eV·s
        c_speed = 299792458  # m/s
        photon_energy = h_planck * c_speed / (wavelengths * 1e-9)
        reflectivity_max = ((1 - 1.49) / (1 + 1.49)) ** 2  # n=1.49 borosilicate
        reflectivity = [reflectivity_max - 0.01] * len(wavelengths)
        reflectivity_matrix_borosilicate = defines.MatrixFromVectors(
            photon_energy, reflectivity, "ReflectivityMatrixBorosilicate", self.g4_registry, "eV", ""
        )
        efficiency_matrix_borosilicate = defines.MatrixFromVectors(
            photon_energy_borosilicate,
            transmittance_borosilicate,
            "TransmittanceMatrixBorosilicate",
            self.g4_registry,
            "eV",
            "",
        )

        self._air_to_borosilicate.addProperty("REFLECTIVITY", reflectivity_matrix_borosilicate)
        self._air_to_borosilicate.addProperty("TRANSMITTANCE", efficiency_matrix_borosilicate)

        return self._air_to_borosilicate
