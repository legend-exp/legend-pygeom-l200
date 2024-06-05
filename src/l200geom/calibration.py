"""Construct the LEGEND-200 calibration system (Gloria Senatore and Sandro Gälli at UZH contributed to this code).
"""

from __future__ import annotations

import math

from pyg4ometry import geant4

from . import materials, hpge_strings

def place_calibration_system(
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
)-> None:
    """Construct LEGEND-200 calibration system.

    Parameters
    ----------
    z0
        The position of the top most slot of the strings in the z-axis.
    mothervolume (see hpge_strings.py)
        pyg4ometry Geant4 LogicalVolume instance in which the strings
        are to be placed (see hpge_strings.py).
    registry
        pyg4ometry Geant4 registry instance.
    """
    if registry is None:
        msg = "registry cannot be None"
        raise ValueError(msg)
        
    # place calibration tubes.
    calib_tube_length = 1400  # note: just a rough guess from MaGe
    calib_tube = hpge_strings._get_nylon_mini_shroud(20, calib_tube_length, materials, registry)
    calib_tube_z = z0 - calib_tube_length / 2

    # all positions from MaGe, might be incorrect!
    geant4.PhysicalVolume(
        [0, 0, 0],
        [121.472, -96.277, calib_tube_z],
        calib_tube,
        "calibration_tube_1",
        mothervolume,
        registry,
    )
    geant4.PhysicalVolume(
        [0, 0, 0],
        [-120.9667, -96.9126, calib_tube_z],
        calib_tube,
        "calibration_tube_2",
        mothervolume,
        registry,
    )
    geant4.PhysicalVolume(
        [0, 0, 0],
        [-121.304, 96.48977, calib_tube_z],
        calib_tube,
        "calibration_tube_3",
        mothervolume,
        registry,
    )
    geant4.PhysicalVolume(
        [0, 0, 0],
        [121.135, 96.70, calib_tube_z],
        calib_tube,
        "calibration_tube_4",
        mothervolume,
        registry,
    )
    
    #build and place the calibration sources and absorbers
    calib_tube_z0 = z0 - calib_tube_length #starting z position of the calib tube in pygeometry reference system

    absorber_height = 37.5
    inner_radius_source = 1.9  # mm
    source_radius = 2  # mm
    source_height = 4  # mm
    absorber_radius = 16
    
    height = {}
    height[1] = calib_tube_z0 + 1/2 * absorber_height #this z positions of the absorbers are only a guess
    height[2] = calib_tube_z0 + 500
    height[3] = calib_tube_z0 + 1/2 * absorber_height
    height[4] = calib_tube_z0 + 500

    sourceLV = _get_calibration_source(inner_radius_source, source_radius, source_height, materials, registry)
    absorberLV = _get_calibration_absorber(absorber_radius, absorber_height, materials, registry) 
    
    base_position = {}
    # Define the z position for the lowest sample and the x,y position as well according to l200. 
    base_position[1] = [121.472, -96.277, height[1] + 1/2 * absorber_height + 1/2 * source_height]
    base_position[2] = [-120.9667, -96.9126, height[2] + 1/2 * absorber_height + 1/2 * source_height]
    base_position[3] = [-121.304, 96.48977, height[3] + 1/2 * absorber_height + 1/2 * source_height]
    base_position[4] = [121.135, 96.70, height[4] + 1/2 * absorber_height + 1/2 * source_height]

    second_lowest_source_pos = 98 + 1.2 + source_height #I add the 1/2 * absorber_height two times. 1.2 mm is the source holder high
    second_highest_source_pos = second_lowest_source_pos + 100 + 1.2 + source_height
    highest_source_pos = second_highest_source_pos + 100 + 1.2 + source_height
        
	# Place cylinders in the world where the middle of the absorber is set to be 0 -> now the middle of the absorber in SIS1 is at 300
    for i in range(1, 5, 1):
        geant4.PhysicalVolume([0, 0, 0], [base_position[i][0], base_position[i][1], base_position[i][2]], sourceLV, f'Lowest Source PV{i}', mothervolume, registry)
        geant4.PhysicalVolume([0, 0, 0], [base_position[i][0], base_position[i][1], base_position[i][2] + second_lowest_source_pos], sourceLV, f'Second Lowest Source PV{i}', mothervolume, registry)
        geant4.PhysicalVolume([0, 0, 0], [base_position[i][0], base_position[i][1], base_position[i][2] + second_highest_source_pos], sourceLV, f'Second Highest Source PV{i}', mothervolume, registry)
        geant4.PhysicalVolume([0, 0, 0], [base_position[i][0], base_position[i][1], base_position[i][2] + highest_source_pos], sourceLV, f'Highest Source PV{i}', mothervolume, registry)
        geant4.PhysicalVolume([0, 0, 0], [base_position[i][0], base_position[i][1], height[i]], absorberLV, f'Absorber PV{i}', mothervolume, registry)
        
        
        
def _get_calibration_source(
    inner_radius: float,
    outer_radius: float,
    length: float,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    """Geometry for calibration sources.
    """

    source = geant4.solid.Tubs("Calibration source", inner_radius, outer_radius, length, 0, 2 * math.pi, registry)
    source_lv = geant4.LogicalVolume(source, materials.radio_thorium, "Calibration source", registry)
  
    source_lv.pygeom_color_rgba = (255, 0, 255, 1)
    
    return source_lv
    
def _get_calibration_absorber(
    outer_radius: float,
    length: float,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    """Geometry for calibration absorbers.
    """

    absorber = geant4.solid.Tubs("Calibration absorber", 0, outer_radius, length, 0, 2 * math.pi, registry)
    absorber_lv = geant4.LogicalVolume(absorber, materials.metal_tantalum, "Calibration absorber", registry)
    
    absorber_lv.pygeom_color_rgba = (0.5, 0.5, 0, 1)
    
    return absorber_lv
