"""VM2000 reflective foil

Implementation used from MaGe and improved with results of https://iopscience.iop.org/article/10.1088/1748-0221/12/06/P06017/pdf."""

from __future__ import annotations

import logging
import math

import numpy as np
import pint
from legendoptics import store
from legendoptics.utils import InterpolatingGraph, readdatafile
from pint import Quantity

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def calculate_wls_mfp(yield_value):
    # Set total path length (currently hardcoded to 1 mm)
    total_path = 1.0e-3 * u.m

    # Handle edge cases
    if yield_value == 0:
        return 10.0  # 10 m,  Large mean free path, no absorption
    if yield_value == 1:
        return 0.01e-3  # 0.01 mm, Very small mean free path, 100% absorption

    # Calculate mean free path for valid yield values
    help_value = math.log(1.0 - yield_value)

    return -total_path / help_value


@store.register_pluggable
@u.with_context("sp")
def vm2000_parameters() -> tuple[Quantity, Quantity, Quantity, Quantity, Quantity]:
    """Wavelength-shifting parameters for the reflective foil VM2000."""
    # Constants
    wls_yield = 0.075  # 0.6 MaGe, 0.075 XENON paper

    # Populate VM2000_energy_range array with energy values
    ppsci_high_e = (115 * u.nm).to("eV")
    ppsc_low_e = (650 * u.nm).to("eV")

    num1 = 251
    dee = (ppsci_high_e - ppsc_low_e) / (num1 - 2)
    vm2000_energy_range = np.zeros(num1) * u.eV
    for ji in range(1, num1):
        vm2000_energy_range[ji] = ppsc_low_e + ji * dee
        vm2000_energy_range[0] = 1.8 * u.eV

    # Create arrays for energy and optical properties
    vm2000_reflectivity = np.zeros(num1)
    vm2000_efficiency = np.zeros(num1)
    wls_absorption = np.zeros(num1) * u.m

    # Set reflectivity, absorption, and emission
    for ji in range(num1):
        if vm2000_energy_range[ji] < (370 * u.nm).to("eV"):  # 370 nm < (related to energy)
            vm2000_reflectivity[ji] = 0.95  # Visible light 0.95, 0.99
        else:
            vm2000_reflectivity[ji] = 0.12  # UV light 0.15, 0.3 (paper)

        if vm2000_energy_range[ji] > 3.35 * u.eV:  # 5 eV 3.35
            # depending on path length in foil --> angle
            wls_absorption[ji] = calculate_wls_mfp(wls_yield)  # Absorbs UV
        else:
            wls_absorption[ji] = 1.0 * u.m  # Imperturbed, no absorption of visible light

    g = InterpolatingGraph(*readdatafile("vm2000_em_spec.dat", pkg="l200geom.materials"), zero_outside=True)
    wls_emission = g(vm2000_energy_range.to("nm")).to("dimensionless")

    # Copy the first element to 0th position
    wls_absorption[0] = wls_absorption[1]  # depending on path length in foil --> angle
    wls_emission[0] = wls_emission[1]

    return vm2000_energy_range, vm2000_reflectivity, vm2000_efficiency, wls_absorption, wls_emission
