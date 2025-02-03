"""Implementation used from MaGe and improved with results of https://iopscience.iop.org/article/10.1088/1748-0221/12/06/P06017/pdf."""

from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import pint
from legendoptics import store
from pint import Quantity

log = logging.getLogger(__name__)
u = pint.get_application_registry()

################################# VM2000 ##########################################


def initialize_vm2000_spectrum():
    frequencyV = np.zeros(251)
    efficiencyV = np.zeros(251)
    successfulInitialization = False
    npoints = 251  # Set number of points
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir / "VM2000_em_spec.dat"

    try:
        with file_path.open() as file:  # "VM2000_em_spec.dat" TPBOnVM2000Emission
            for i in range(npoints):
                line = file.readline().split()
                if line:
                    gg, hh = map(float, line)
                    frequencyV[i] = gg * 1e-9  # Convert nanometers to meters
                    efficiencyV[i] = hh
            successfulInitialization = True

    except FileNotFoundError:
        for i in range(251):
            frequencyV[i] = (i + 350) * 1e-9  # Convert from mm to nanometers
            efficiencyV[i] = 20
        successfulInitialization = True

    return frequencyV, efficiencyV, successfulInitialization, npoints


frequencyV, efficiencyV, successfulInitialization, npoints = initialize_vm2000_spectrum()


def VM2000EmissionSpectrum(frequencyV, efficiencyV, successfulInitialization, npoints, energy):
    """
    Calculates the emission spectrum efficiency based on the input photon energy.

    :param energy: Photon energy in GeV.
    :return: Efficiency of the spectrum for the given energy.
    """
    if successfulInitialization:
        # Convert energy to wavelength using the LambdaE constant
        LambdaE = 2 * np.pi * 1.973269602e-7  # m * eV
        targetf = LambdaE / energy

        # Check if the wavelength is within the bounds of our frequencyV array
        if targetf < frequencyV[0] or targetf > frequencyV[npoints - 1]:
            return 0.0

        # Find the index j such that frequencyV[j] <= targetf <= frequencyV[j+1]
        j = np.searchsorted(frequencyV, targetf) - 1
        if j < 0 or j >= npoints - 1:
            return 0.0  # Ensure index is within bounds

        # Perform linear interpolation between points j and j+1
        eff = (
            (targetf - frequencyV[j])
            * (efficiencyV[j + 1] - efficiencyV[j])
            / (frequencyV[j + 1] - frequencyV[j])
        )
        eff += efficiencyV[j]
        return eff
    # Return default efficiency if initialization failed
    return 0.2


def CalculateWLSmfp(yield_value):
    # Set total path length (currently hardcoded to 1 mm)
    totalPath = 1.0e-3  # 1 mm

    # Handle edge cases
    if yield_value == 0:
        return 10.0  # 10 m,  Large mean free path, no absorption
    if yield_value == 1:
        return 0.01e-3  # 0.01 mm, Very small mean free path, 100% absorption

    # Calculate mean free path for valid yield values
    help_value = math.log(1.0 - yield_value)

    return -totalPath / help_value


@store.register_pluggable
def VM2000_parameters() -> tuple[Quantity, Quantity, Quantity, Quantity, Quantity]:
    """Wavelength-shifting parameters for the reflective foil VM2000."""
    # Constants
    # Energy goes from UV (115nm) to green (650 nm)
    WLSyield = 0.075  # 0.6 MaGe, 0.075 XENON paper

    LambdaE = 2 * np.pi * 1.973269602e-7  # m*eV

    PPSCHighE = LambdaE / (115e-9)  # 115 nm
    PPSCLowE = LambdaE / (650e-9)  # 650 nm

    num = 250
    num1 = 251
    dee = (PPSCHighE - PPSCLowE) / (num - 1)

    # Create arrays for energy and optical properties
    VM2000_energy_range = np.zeros(num1)
    VM2000_Reflectivity = np.zeros(num1)
    VM2000_Efficiency = np.zeros(num1)
    WLS_absorption = np.zeros(num1)
    WLS_emission = np.zeros(num1)

    # Populate VM2000_energy_range array with energy values
    for ji in range(1, num1):
        ee = PPSCLowE + ji * dee
        VM2000_energy_range[ji] = ee
    VM2000_energy_range[0] = 1.8  # eV

    # Set reflectivity, absorption, and emission
    for ji in range(1, num1):
        if VM2000_energy_range[ji] < (LambdaE / (370e-9)):  # 370 nm < (related to energy)
            VM2000_Reflectivity[ji] = 0.95  # Visible light 0.95, 0.99
        else:
            VM2000_Reflectivity[ji] = 0.12  # UV light 0.15, 0.3 (paper)

        VM2000_Efficiency[ji] = 0.0

        if VM2000_energy_range[ji] > 3.35:  # 5 eV 3.35
            WLS_absorption[ji] = CalculateWLSmfp(WLSyield)  # Absorbs UV
        else:
            WLS_absorption[ji] = 1.0  # 1 m, Imperturbed, no absorption of visible light

        WLS_emission[ji] = VM2000EmissionSpectrum(
            frequencyV, efficiencyV, successfulInitialization, npoints, VM2000_energy_range[ji]
        )

    # Copy the first element to 0th position
    VM2000_Reflectivity[0] = VM2000_Reflectivity[1]
    VM2000_Efficiency[0] = VM2000_Efficiency[1]
    WLS_absorption[0] = WLS_absorption[1]  # depending on path length in foil --> angle
    WLS_emission[0] = WLS_emission[1]

    return VM2000_energy_range, VM2000_Reflectivity, VM2000_Efficiency, WLS_absorption, WLS_emission
