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
    frequencyv = np.zeros(251)
    efficiencyv = np.zeros(251)
    successful_initialization = False
    npoints = 251  # Set number of points
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir / "vm2000_em_spec.dat"

    with file_path.open() as file:  # "VM2000_em_spec.dat"
        for i in range(npoints):
            line = file.readline().split()
            if line:
                gg, hh = map(float, line)
                frequencyv[i] = gg * 1e-9  # Convert nanometers to meters
                efficiencyv[i] = hh
        successful_initialization = True

    return frequencyv, efficiencyv, successful_initialization, npoints


frequencyv, efficiencyv, successful_initialization, npoints = initialize_vm2000_spectrum()


def vm2000_emission_spectrum(frequencyv, efficiencyv, successful_initialization, npoints, energy):
    """
    Calculates the emission spectrum efficiency based on the input photon energy.

    :param energy: Photon energy in GeV.
    :return: Efficiency of the spectrum for the given energy.
    """
    if successful_initialization:
        # Convert energy to wavelength using the LambdaE constant
        lambda_e = 2 * np.pi * 1.973269602e-7  # m * eV
        targetf = lambda_e / energy

        # Check if the wavelength is within the bounds of our frequencyV array
        if targetf < frequencyv[0] or targetf > frequencyv[npoints - 1]:
            return 0.0

        # Find the index j such that frequencyV[j] <= targetf <= frequencyV[j+1]
        j = np.searchsorted(frequencyv, targetf) - 1
        if j < 0 or j >= npoints - 1:
            return 0.0  # Ensure index is within bounds

        # Perform linear interpolation between points j and j+1
        eff = (
            (targetf - frequencyv[j])
            * (efficiencyv[j + 1] - efficiencyv[j])
            / (frequencyv[j + 1] - frequencyv[j])
        )
        eff += efficiencyv[j]
        return eff
    # Return default efficiency if initialization failed
    return 0.2


def calculate_wls_mfp(yield_value):
    # Set total path length (currently hardcoded to 1 mm)
    total_path = 1.0e-3  # 1 mm

    # Handle edge cases
    if yield_value == 0:
        return 10.0  # 10 m,  Large mean free path, no absorption
    if yield_value == 1:
        return 0.01e-3  # 0.01 mm, Very small mean free path, 100% absorption

    # Calculate mean free path for valid yield values
    help_value = math.log(1.0 - yield_value)

    return -total_path / help_value


@store.register_pluggable
def vm2000_parameters() -> tuple[Quantity, Quantity, Quantity, Quantity, Quantity]:
    """Wavelength-shifting parameters for the reflective foil VM2000."""
    # Constants
    # Energy goes from UV (115nm) to green (650 nm)
    wls_yield = 0.075  # 0.6 MaGe, 0.075 XENON paper

    lambda_e = 2 * np.pi * 1.973269602e-7  # m*eV

    ppsci_high_e = lambda_e / (115e-9)  # 115 nm
    ppsc_low_e = lambda_e / (650e-9)  # 650 nm

    num = 250
    num1 = 251
    dee = (ppsci_high_e - ppsc_low_e) / (num - 1)

    # Create arrays for energy and optical properties
    vm2000_energy_range = np.zeros(num1)
    vm2000_reflectivity = np.zeros(num1)
    vm2000_efficiency = np.zeros(num1)
    wls_absorption = np.zeros(num1)
    wls_emission = np.zeros(num1)

    # Populate VM2000_energy_range array with energy values
    for ji in range(1, num1):
        ee = ppsc_low_e + ji * dee
        vm2000_energy_range[ji] = ee
    vm2000_energy_range[0] = 1.8  # eV

    # Set reflectivity, absorption, and emission
    for ji in range(1, num1):
        if vm2000_energy_range[ji] < (lambda_e / (370e-9)):  # 370 nm < (related to energy)
            vm2000_reflectivity[ji] = 0.95  # Visible light 0.95, 0.99
        else:
            vm2000_reflectivity[ji] = 0.12  # UV light 0.15, 0.3 (paper)

        vm2000_efficiency[ji] = 0.0

        if vm2000_energy_range[ji] > 3.35:  # 5 eV 3.35
            wls_absorption[ji] = calculate_wls_mfp(wls_yield)  # Absorbs UV
        else:
            wls_absorption[ji] = 1.0  # 1 m, Imperturbed, no absorption of visible light

        wls_emission[ji] = vm2000_emission_spectrum(
            frequencyv, efficiencyv, successful_initialization, npoints, vm2000_energy_range[ji]
        )

    # Copy the first element to 0th position
    vm2000_reflectivity[0] = vm2000_reflectivity[1]
    vm2000_efficiency[0] = vm2000_efficiency[1]
    wls_absorption[0] = wls_absorption[1]  # depending on path length in foil --> angle
    wls_emission[0] = wls_emission[1]

    return vm2000_energy_range, vm2000_reflectivity, vm2000_efficiency, wls_absorption, wls_emission
