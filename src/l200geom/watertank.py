"""Construct the LEGEND-200/GERDA watertank including the water volume and the reflective foil VM2000.

"""


from __future__ import annotations

from math import pi
import math
import numpy as np
import pyg4ometry as pyg4
import pyg4ometry.geant4 as g4
from scipy.spatial.transform import Rotation as R
from . import materials, cryo


PMT_ID = np.array(["Ch29_312", "Ch18_301", "Ch16_701", "Ch19_302", "Ch20_303", "Ch21_304", "Ch22_305", "Ch23_306", "Ch17_703", "Ch24_307", "Ch25_308", "Ch26_309", "Ch27_310", "Ch28_311", \
                   "Ch15_208", "Ch10_201", "Ch11_202", "Ch12_203", "Ch13_706", "Ch14_206", \
                    "Ch01_704", "Ch03_705", "Ch06_709", "Ch08_710", \
                    "Ch30_401", "Ch31_402", "Ch32_403", "Ch33_404", "Ch34_409", "Ch35_410", \
                    "Ch43_510", "Ch36_501", "Ch37_502", "Ch38_503", "Ch39_504", "Ch40_507", "Ch41_508", "Ch42_509", \
                    "Ch44_602", "Ch45_603", "Ch46_605", "Ch47_606", "Ch48_607", "Ch49_608", "Ch50_609", "Ch51_610", \
                    "Ch52_702", \
                    "Ch00_101", "Ch02_102", "Ch04_708", "Ch05_104", "Ch07_105", "Ch09_707"])


# Water tank with water and air buffer
WaterTankThickness = 7.0
WaterTankHeight = 8900.0
InnerTankHeight = WaterTankHeight - 2 * WaterTankThickness
InnerRadius = 0.0
WaterRadius = 5000.0
WaterHeight = 8400.0

# Reflective foil
ReflectiveFoilThickness = 0.04

# Pillbox
ShieldingFootOR = 2000.0 
ShieldingFootThickness = 1.2
ShieldingFootIR = ShieldingFootOR - ShieldingFootThickness - ReflectiveFoilThickness
CryoBottomHeight = 1300.0 # K. Freund https://publikationen.uni-tuebingen.de/xmlui/bitstream/handle/10900/57937/KFreundPhD.pdf?sequence=1&isAllowed=y 
PillboxOffset = - WaterHeight/2 + 0.5 * CryoBottomHeight
PillboxTubeFoilOffset = PillboxOffset

# Air buffer
OuterWaterTankRadius = WaterRadius + WaterTankThickness
AirBufferRadius = WaterRadius
AirBufferHeight = WaterTankHeight - 2.0 * WaterTankThickness - WaterHeight


# z-axis Offsets
AirBufferOffset = 0.5 * (InnerTankHeight - AirBufferHeight)
TankOffset = 0.5 * AirBufferHeight
BottomFoilOffset = - 0.5 * WaterHeight + 0.5 * ReflectiveFoilThickness


# PMTs
PMTStartingAngle = 0.0  # in degrees
PMTEndingAngle = 2 * pi  # in degrees

PhotocathodeInnerRadius = 100.0
PhotocathodeOuterRadius = 100.16
PhotocathodeThetaStart = 0.0  # in degrees
PhotocathodeThetaEnd = (80.0 / 180.0) * pi  # in degrees
PhotocatodeHeight = PhotocathodeOuterRadius * (1 - math.cos(PhotocathodeThetaEnd))
PhotocathodeHeightDifference = PhotocathodeOuterRadius - PhotocatodeHeight

PMTInnerRadius = 100.16
PMTOuterRadius = 101.0
PMTThetaIn = 0.0  # in degrees
PMTThetaEnd = 0.5 * pi  # in degrees
BorosilikatGlassThickness = 1.0

PMTSteelConeThickness = 4.0 
PMTSteelConeHeight = 300.0 
PMTSteelConeUpperRmin = PMTOuterRadius
PMTSteelConeUpperRmax = PMTOuterRadius + PMTSteelConeThickness
PMTSteelConeLowerRmin = PMTOuterRadius / 3.0
PMTSteelConeLowerRmax = PMTSteelConeLowerRmin + PMTSteelConeThickness

PMTBorosilikatGlassInnerRadius = PhotocathodeInnerRadius
PMTBorosilikatGlassOuterRadius = PhotocathodeOuterRadius + BorosilikatGlassThickness
PMTBorosilikatGlassThetaStart = 0.0  # in degrees
PMTBorosilikatGlassThetaEnd = (80.0 / 180.0) * pi  # in degrees
PMTBorosilikatGlassHeight = PMTBorosilikatGlassOuterRadius * (1 - math.cos(PMTBorosilikatGlassThetaEnd))
PMTBorosilikatGlassHeightDifference = PMTBorosilikatGlassOuterRadius - PMTBorosilikatGlassHeight

PMTAirInnerRadius = PhotocathodeInnerRadius
PMTAirOuterRadius = PMTSteelConeUpperRmax
PMTAirThetaStart = 0.0  # in degrees
PMTAirThetaEnd = (80.0 / 180.0) * pi  # in degrees
PMTAirHeight = PMTAirOuterRadius * (1 - math.cos(PMTAirThetaEnd))
PMTAirHeightDifference = PMTAirOuterRadius - PMTAirHeight

AcrylInnerRadius = PhotocathodeInnerRadius
AcrylOuterRadius = PMTSteelConeUpperRmax + 5.5
AcrylThetaStart = 0.0  # in degrees
AcrylThetaEnd = (80.0 / 180.0) * pi  # in degrees
AcrylHeight = AcrylOuterRadius * (1 - math.cos(AcrylThetaEnd))
AcrylHeightDifference = AcrylOuterRadius - AcrylHeight


PMTSteelBottomHeight = 30.0 

PMTHeight = PhotocatodeHeight + 0.5*PMTSteelConeHeight + PMTSteelBottomHeight

PMTCathodeOffset = - WaterHeight/2 + ReflectiveFoilThickness - PhotocathodeHeightDifference + 0.5 * PMTSteelConeHeight + PMTSteelBottomHeight
PMTConeOffset = - WaterHeight/2 + ReflectiveFoilThickness + 0.25 * PMTSteelConeHeight + PMTSteelBottomHeight
PMTBottomOffset = - WaterHeight/2 + ReflectiveFoilThickness + 0.5 * PMTSteelBottomHeight

DistancePMTBaseTank = WaterRadius - ReflectiveFoilThickness - np.sqrt((WaterRadius - ReflectiveFoilThickness)**2 - PMTSteelConeLowerRmax**2)


def construct_tank(reg: g4.Registry, tank_material: g4.Material)-> g4.LogicalVolume:
    WaterTankWall = g4.solid.Tubs("WaterTankWall", WaterRadius, OuterWaterTankRadius, InnerTankHeight, 0, 2*pi, reg) 
    WaterTankFloor = g4.solid.Tubs("WaterTankFloor", InnerRadius, OuterWaterTankRadius, WaterTankThickness, 0, 2*pi, reg)

    tank_union1_transform = [[0, 0, 0], [0, 0, -0.5 * (InnerTankHeight + WaterTankThickness)]]
    tank_union2_transform = [[0, 0, 0], [0, 0, 0.5 * (InnerTankHeight + WaterTankThickness)]]

    WaterTankUnionSolid1 = g4.solid.Union("WaterTankUnionSolid1", WaterTankWall, WaterTankFloor, tank_union1_transform, reg)
    WaterTankUnionSolid2 = g4.solid.Union("WaterTankUnionSolid2", WaterTankUnionSolid1, WaterTankFloor, tank_union2_transform, reg)
    return g4.LogicalVolume(WaterTankUnionSolid2, tank_material, "WaterTankUnionSolid_lv", reg) 


def place_tank(reg:g4.Registry, WaterTankUnionSolid_lv: g4.LogicalVolume, world_lv: g4.LogicalVolume, TankOffset: float)-> g4.PhysicalVolume:
    return g4.PhysicalVolume([0,0,0], [0,0,TankOffset], WaterTankUnionSolid_lv, "WaterTankUnionSolid_pv", world_lv, reg) 

def construct_water(reg: g4.Registry, water_material: g4.Material)-> g4.LogicalVolume:
    water_solid = g4.solid.Tubs("water", InnerRadius, WaterRadius, WaterHeight, 0, 2*pi, reg) 
    water_lv = g4.LogicalVolume(water_solid, water_material, "water_lv", reg) 
    water_lv.pygeom_color_rgba = [0, 0, 1, 0.08]
    return water_lv

def place_water(reg: g4.Registry, water_lv: g4.LogicalVolume, world_lv: g4.LogicalVolume) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0,0,0], [0,0,-AirBufferHeight/2], water_lv, "water_pv", world_lv, reg) 

def construct_air_buffer(reg: g4.Registry, air_material: g4.Material)-> g4.LogicalVolume:
    air_buffer = g4.solid.Tubs("air_buffer", InnerRadius, AirBufferRadius, AirBufferHeight, 0, 2*pi, reg) 
    return g4.LogicalVolume(air_buffer, air_material, "air_buffer_lv", reg) 

def place_air_buffer(reg: g4.Registry, air_buffer_lv: g4.LogicalVolume, world_lv: g4.LogicalVolume) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0,0,0], [0,0,AirBufferOffset], air_buffer_lv, "air_buffer_pv", world_lv, reg) 

def construct_pillbox(reg: g4.Registry, pillbox_material: g4.Material)-> g4.LogicalVolume:
    ManholeOuterRadius = 400.0
    x_rot_drehvol = 0
    y_rot_drehvol = (np.pi/2.0)
    z_rot_drehvol = (np.pi/2.0)


    # rotation matrix around Z, Y and X
    rot_z = R.from_euler('z', z_rot_drehvol, degrees=False)
    rot_y = R.from_euler('y', y_rot_drehvol, degrees=False)
    rot_x = R.from_euler('x', x_rot_drehvol, degrees=False)

    # combine rotations
    combined_rotation = rot_y * rot_z * rot_x

    # rotation matrix in euler angle
    euler_angles = combined_rotation.as_euler('xyz', degrees=False)
    x_rot_global, y_rot_global, z_rot_global = euler_angles


    PillboxTube = g4.solid.Tubs("PillboxTube", ShieldingFootIR, ShieldingFootOR, CryoBottomHeight-ShieldingFootThickness, 0, 2*pi, reg) # outer steel cylinder

    # Define parameters for the semi-cylinder (half-tube) for the manhole
    ManholeInnerRadius = 0  # No inner radius for the manhole
    ManholeHeight = 2*(ShieldingFootOR + ReflectiveFoilThickness)
    ManholeAngle = math.pi  # Half-circle (180 degrees)

    # Create the half-tube (semi-cylinder) for the manhole
    ManholePillboxArc = g4.solid.Tubs("ManholePillbox", ManholeInnerRadius, ManholeOuterRadius, ManholeHeight, 0, ManholeAngle, reg)
    ManholePillboxBox = g4.solid.Box("ManholePillboxBox", 2*ManholeOuterRadius, ManholeOuterRadius, ManholeHeight, reg)

    # Position the first manhole (half-tube) along the x-axis
    # Rotate the manhole to align it with the x-axis (rotating by 90 degrees around the y-axis)
    ManholeRotation = [x_rot_global, y_rot_global, z_rot_global]
    ManholePosition1 = [ShieldingFootOR / 2, 0, 0]  # Adjust based on the geometry
    union_transform = [[0, 0, 0], [0, -0.5 * ManholeOuterRadius, 0]]
    ManholePillbox = g4.solid.Union("ManholeUnion", ManholePillboxArc, ManholePillboxBox, union_transform, reg)

    # Subtract the first manhole (half-tube) from the pillbox
    ManHoleOffset=(0.5*CryoBottomHeight-ManholeOuterRadius)
    pillbox = g4.solid.Subtraction("pillbox_subtraction1", PillboxTube, ManholePillbox, [ManholeRotation, [0,0,0-ManHoleOffset]], reg)
    pillbox_lv = g4.LogicalVolume(pillbox, pillbox_material, "pillbox_lv", reg)

    return  pillbox_lv, ManholePillbox, ManholeRotation, ManHoleOffset

def place_pillbox(reg: g4.Registry, pillbox_lv: g4.LogicalVolume, water_lv: g4.LogicalVolume) -> g4.PhysicalVolume:
    return g4.PhysicalVolume([0,0,0], [0,0,PillboxOffset], pillbox_lv, "pillbox_pv", water_lv, reg) 


def insert_VM2000(reg: g4.Registry, VM2000_material: g4.Material, water_lv: g4.LogicalVolume, water_pv: g4.PhysicalVolume, ManholePillbox: g4.solid.Union, ManholeRotation: list , ManHoleOffset:float, cryo_displacement_z: float):

    # VM2000 at inside of water tank tube
    WaterTankReflectionFoilTube = g4.solid.Tubs("WaterTankReflectionFoilTube", WaterRadius - ReflectiveFoilThickness, WaterRadius, WaterHeight, 0, 2*pi, reg) 
    WaterTankReflectionFoilTube_lv = g4.LogicalVolume(WaterTankReflectionFoilTube, VM2000_material, "WaterTankReflectionFoilTube_lv", reg) 
    WaterTankReflectionFoilTube_pv = g4.PhysicalVolume([0,0,0], [0,0,0], WaterTankReflectionFoilTube_lv, "WaterTankReflectionFoilTube_pv", water_lv, reg) 

    # VM2000 at bottom of water tank tube
    WaterTankReflectionFoilBottom = g4.solid.Tubs("WaterTankReflectionFoilBottom", ShieldingFootOR + ReflectiveFoilThickness, WaterRadius - ReflectiveFoilThickness, ReflectiveFoilThickness, 0, 2*pi, reg) 
    WaterTankReflectionFoilBottom_lv = g4.LogicalVolume(WaterTankReflectionFoilBottom, VM2000_material, "WaterTankReflectionFoilBottom_lv", reg)
    WaterTankReflectionFoilBottom_pv = g4.PhysicalVolume([0,0,0], [0,0,BottomFoilOffset], WaterTankReflectionFoilBottom_lv, "WaterTankReflectionFoilBottom_pv", water_lv, reg) 

    # Pillbox
    # VM2000 at inside of water tank tube
    PillboxOuterReflectionFoilTubeSubtraction1 = g4.solid.Tubs("PillboxOuterReflectionFoilTubeSubtraction1", ShieldingFootOR, ShieldingFootOR + ReflectiveFoilThickness, CryoBottomHeight, 0, 2*pi, reg) 
    PillboxOuterReflectionFoilTube = g4.solid.Subtraction("PillboxOuterReflectionFoilTube ", PillboxOuterReflectionFoilTubeSubtraction1, ManholePillbox, [ManholeRotation, [0,0,0-ManHoleOffset]], reg)
    PillboxOuterReflectionFoilTube_lv = g4.LogicalVolume(PillboxOuterReflectionFoilTube, VM2000_material, "PillboxOuterReflectionFoilTube_lv", reg) 
    PillboxOuterReflectionFoilTube_pv = g4.PhysicalVolume([0,0,0], [0,0,PillboxOffset], PillboxOuterReflectionFoilTube_lv, "PillboxOuterReflectionFoilTube_pv", water_lv, reg) 

    # VM2000 at inside of water tank tube
    PillboxInnerReflectionFoilTubeSubtraction1 = g4.solid.Tubs("PillboxInnerReflectionFoilTubeSubtraction1", ShieldingFootIR - ReflectiveFoilThickness, ShieldingFootIR, CryoBottomHeight-ShieldingFootThickness, 0, 2*pi, reg)
    PillboxInnerReflectionFoilTube = g4.solid.Subtraction("PillboInnerReflectionFoilTube ", PillboxInnerReflectionFoilTubeSubtraction1, ManholePillbox, [ManholeRotation, [0,0,0-ManHoleOffset]], reg)
    PillboxInnerReflectionFoilTube_lv = g4.LogicalVolume(PillboxInnerReflectionFoilTube, VM2000_material, "PillboxInnerReflectionFoilTube_lv", reg)
    PillboxInnerReflectionFoilTube_pv = g4.PhysicalVolume([0,0,0], [0,0,PillboxOffset], PillboxInnerReflectionFoilTube_lv, "PillboxInnerReflectionFoilTube_pv", water_lv, reg) 

    # VM2000 at bottom of water tank tube
    PillboxReflectionFoilTop = g4.solid.Tubs("PillboxReflectionFoilTop", InnerRadius, ShieldingFootIR, ReflectiveFoilThickness, 0, 2*pi, reg)
    PillboxReflectionFoilTop_lv = g4.LogicalVolume(PillboxReflectionFoilTop, VM2000_material, "PillboxReflectionFoilTop_lv", reg)
    PillboxReflectionFoilTop_pv = g4.PhysicalVolume([0,0,0], [0,0,BottomFoilOffset+CryoBottomHeight-2*ReflectiveFoilThickness], PillboxReflectionFoilTop_lv, "PillboxReflectionFoilTop_pv", water_lv, reg) 

    # VM2000 at bottom of water tank tube
    PillboxReflectionFoilBottom = g4.solid.Tubs("PillboxReflectionFoilBottom", InnerRadius, ShieldingFootIR, ReflectiveFoilThickness, 0, 2*pi, reg)
    PillboxReflectionFoilBottom_lv = g4.LogicalVolume(PillboxReflectionFoilBottom, VM2000_material, "PillboxReflectionFoilBottom_lv", reg)
    PillboxReflectionFoilBottom_pv = g4.PhysicalVolume([0,0,0], [0,0,BottomFoilOffset], PillboxReflectionFoilBottom_lv, "PillboxReflectionFoilBottom_pv", water_lv, reg)

    #TODO: adjust height and z position of CryoReflectionFoil to close the gap in between PillboxOuterReflectionFoilTube and CryoReflectionFoil
    CryoReflectionFoil = g4.solid.Tubs("CryoReflectionFoil", cryo.cryo_radius + cryo.cryo_wall, cryo.cryo_radius + cryo.cryo_wall + ReflectiveFoilThickness, cryo.cryo_tub_height+cryo.cryo_top_height+cryo.cryo_bottom_height+2*cryo.cryo_wall, 0, 2*pi, reg)
    CryoReflectionFoil_lv = g4.LogicalVolume(CryoReflectionFoil, VM2000_material, "CryoReflectionFoil_lv", reg) 
    CryoReflectionFoil_pv = g4.PhysicalVolume([0,0,0], [0,0,cryo_displacement_z+cryo.access_overlap+2*cryo.cryo_wall], CryoReflectionFoil_lv, "CryoReflectionFoil_pv", water_lv, reg) 

    
    VM2000BorderOptTable = materials.surfaces.OpticalSurfaceRegistry(reg).water_to_VM2000
    VM2000OptTable = materials.surfaces.OpticalSurfaceRegistry(reg).to_VM2000

    #Border Surfaces
    g4.BorderSurface("PillboxOuterTubeFoilBorderSurface", PillboxOuterReflectionFoilTube_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("PillboxInneTubeFoilBorderSurface", PillboxInnerReflectionFoilTube_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("PillboxBottomFoilBorderSurface", PillboxReflectionFoilBottom_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("PillboxTopFoilBorderSurface", PillboxReflectionFoilTop_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("WaterTankTubeFoilBorderSurface", WaterTankReflectionFoilTube_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("WaterTankBottomFoilBorderSurface", WaterTankReflectionFoilBottom_pv, water_pv, VM2000BorderOptTable, reg)
    g4.BorderSurface("CryoBorderSurface", CryoReflectionFoil_pv, water_pv, VM2000BorderOptTable, reg)
    

    g4.SkinSurface("PillboxOuterTubeFoilSkinSurface", PillboxOuterReflectionFoilTube_lv, VM2000OptTable, reg)
    g4.SkinSurface("PillboxInneTubeFoilSkinSurface", PillboxInnerReflectionFoilTube_lv, VM2000OptTable, reg)
    g4.SkinSurface("PillboxBottomFoilSkinSurface", PillboxReflectionFoilBottom_lv, VM2000OptTable, reg)
    g4.SkinSurface("PillboxTopFoilSkinSurface", PillboxReflectionFoilTop_lv, VM2000OptTable, reg)
    g4.SkinSurface("WaterTankTubeFoilSkinSurface", WaterTankReflectionFoilTube_lv, VM2000OptTable, reg)
    g4.SkinSurface("WaterTankBottomFoilSkinSurface", WaterTankReflectionFoilBottom_lv, VM2000OptTable, reg) 
    g4.SkinSurface("CryoSkinFoilSurface", CryoReflectionFoil_lv, VM2000OptTable, reg)

    return WaterTankReflectionFoilTube_pv, WaterTankReflectionFoilBottom_pv, PillboxOuterReflectionFoilTube_pv, PillboxInnerReflectionFoilTube_pv, \
        PillboxReflectionFoilTop_pv, PillboxReflectionFoilBottom_pv
    

def insert_PMTs(reg: g4.Registry, PMT_steel_material: g4.Material, cathodeAl: g4.Material, PMT_air_material: g4.Material, water_lv: g4.LogicalVolume,
                water_pv: g4.PhysicalVolume, acryl_material: g4.Material, borosilicate_material: g4.Material, PMT_configuration: str = 'LEGEND200'):

    # Photocathode and PMT encapsulation
    Acryl = g4.solid.Sphere("Acryl",AcrylInnerRadius,AcrylOuterRadius,PMTStartingAngle,PMTEndingAngle,AcrylThetaStart,AcrylThetaEnd, reg)
    PMTAir = g4.solid.Sphere("PMTAir",PhotocathodeInnerRadius,PMTAirOuterRadius,PMTStartingAngle,PMTEndingAngle,PhotocathodeThetaStart,PhotocathodeThetaEnd, reg)
    PMTBorosilikatGlass = g4.solid.Sphere("PMTBorosilikatGlass",PhotocathodeInnerRadius,PMTBorosilikatGlassOuterRadius,PMTStartingAngle,PMTEndingAngle,PhotocathodeThetaStart,PhotocathodeThetaEnd, reg)
    Photocathode = g4.solid.Sphere("Photocathode",PhotocathodeInnerRadius,PhotocathodeOuterRadius,PMTStartingAngle,PMTEndingAngle,PhotocathodeThetaStart,PhotocathodeThetaEnd, reg)

    optical_steel_surface = materials.surfaces.OpticalSurfaceRegistry(reg).to_steel
    optical_PMT_surface = materials.surfaces.OpticalSurfaceRegistry(reg).to_photocathode
    optical_border_air_acryl = materials.surfaces.OpticalSurfaceRegistry(reg).acryl_to_air
    optical_border_water_acryl = materials.surfaces.OpticalSurfaceRegistry(reg).water_to_acryl
    optical_border_air_borosilicate = materials.surfaces.OpticalSurfaceRegistry(reg).air_to_borosilicate


    # PMT encapsulation steel cone for Cherenkov veto
    PMTSteelCone = g4.solid.Cons("PMTSteelCone",PMTSteelConeLowerRmin,PMTSteelConeLowerRmax,PMTSteelConeUpperRmin,PMTSteelConeUpperRmax,PMTSteelConeHeight * 0.5,PMTStartingAngle,
                                    PMTEndingAngle, reg)
    PMTSteelCone_lv = g4.LogicalVolume(PMTSteelCone,PMT_steel_material,"PMTSteelCone_lv", reg)
    g4.SkinSurface("PMTConeOpticalSurface", PMTSteelCone_lv, optical_steel_surface, reg)

    # PMT encapsulation bottom for Cherenkov veto
    PMTSteelBottom = g4.solid.Tubs("PMTSteelBottom",0,PMTSteelConeLowerRmax,PMTSteelBottomHeight,PMTStartingAngle,PMTEndingAngle, reg)
    PMTSteelBottom_lv = g4.LogicalVolume(PMTSteelBottom,PMT_steel_material,"PMTSteelBottom_lv", reg)
    g4.SkinSurface("PMTBottomOpticalSurface", PMTSteelBottom_lv, optical_steel_surface, reg)


    ############################################ PMT positions #########################################
    Num_PMTs = 0
    Working_PMTs = 0

    # Bottom outer ring
    N_PMT_per_ring_bo = 24  # spaces, in total 14 PMTs
    R_pos = 4250.0
    dPhi = 2.0 * np.pi / N_PMT_per_ring_bo
    dPhi_c = dPhi

    no_PMT_in_outer_bottom_ring = np.array([3, 5, 7, 9, 11, 15, 17, 19, 21, 23])
    broken_but_inside_PMT_in_outer_bottom_ring_GERDA = np.array([])
    PMT_amount_outer_bottom_ring = N_PMT_per_ring_bo - len(no_PMT_in_outer_bottom_ring)

    # PMTs
    for k in range(0, N_PMT_per_ring_bo):
        if k in no_PMT_in_outer_bottom_ring:
            continue
        dPhi_c = dPhi *k
        xpos = R_pos * np.cos(dPhi_c)
        ypos = R_pos * np.sin(dPhi_c)

        # get names
        if k in broken_but_inside_PMT_in_outer_bottom_ring_GERDA:
            namephotocathode_lv = f"GERDA_PMTcathode_lv_{Num_PMTs}"
            namephotocathode = f"GERDA_PMTcathode_{Num_PMTs}"
            namesteelcone = f"GERDA_PMTcone_{Num_PMTs}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)
            namePMTair_lv = f"GERDA_PMTair_lv_{Num_PMTs}"
            namePMTair = f"GERDA_PMTair_{Num_PMTs}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"GERDA_PMTacryl_lv_{Num_PMTs}"
            namePMTacryl = f"GERDA_PMTacryl_{Num_PMTs}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)  
            namesteelbottom = f"GERDA_PMTbottom_{Num_PMTs}"
            namePMTborosilikat_lv = f"GERDA_PMTborosilikat_lv_{Num_PMTs}"
            namePMTborosilikat = f"GERDA_PMTborosilikat_{Num_PMTs}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg) 
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
        else:
            namephotocathode_lv = f"PMTcathode_lv_{PMT_ID[Working_PMTs]}"
            namephotocathode = f"PMTcathode_{PMT_ID[Working_PMTs]}"
            namesteelcone = f"PMTcone_{PMT_ID[Working_PMTs]}"
            namesteelbottom = f"PMTbottom_{PMT_ID[Working_PMTs]}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)
            namePMTair_lv = f"PMTair_lv_{PMT_ID[Working_PMTs]}"
            namePMTair = f"PMTair_{PMT_ID[Working_PMTs]}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"PMTacryl_lv_{PMT_ID[Working_PMTs]}"
            namePMTacryl = f"PMTacryl_{PMT_ID[Working_PMTs]}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)   

            namePMTborosilikat_lv = f"PMTborosilikat_lv_{PMT_ID[Working_PMTs]}"
            namePMTborosilikat = f"PMTborosilikat_{PMT_ID[Working_PMTs]}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg)   
            g4.BorderSurface(f"AirAcrylSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, Acryl_pv, optical_border_air_acryl, reg)
            g4.BorderSurface(f"WaterAcrylSurface{PMT_ID[Working_PMTs]}", water_pv, Acryl_pv, optical_border_water_acryl, reg)
            g4.BorderSurface(f"AirBorosilicateSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, borosilikat_pv, optical_border_air_borosilicate, reg)
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
            g4.SkinSurface(f"PMTCathodeSkinSurface{PMT_ID[Working_PMTs]}", Photocathode_lv, optical_PMT_surface, reg)
            Working_PMTs += 1
        
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTConeOffset], PMTSteelCone_lv, namesteelcone, water_lv, reg)
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTBottomOffset], PMTSteelBottom_lv, namesteelbottom, water_lv, reg)
        
        dPhi_c += dPhi
        Num_PMTs += 1

    # Bottom inner ring
    N_PMT_per_ring_bi = 8  # 8 PMTs in inner ring
    R_pos = 2750.0 # radius of PMT positions
    dPhi = 2.0 * np.pi / N_PMT_per_ring_bi  
    dPhi_c = dPhi

    broken_but_inside_PMT_in_inner_bottom_ring_GERDA = np.array([4, 7])
    no_PMT_in_inner_bottom_ring = np.array([])
    PMT_amount_inner_bottom_ring = N_PMT_per_ring_bi - len(no_PMT_in_inner_bottom_ring)

    # place PMTs
    for k in range(0, N_PMT_per_ring_bi):
        if k in no_PMT_in_inner_bottom_ring:
            continue
        dPhi_c = dPhi *k
        xpos = R_pos * np.cos(dPhi_c)
        ypos = R_pos * np.sin(dPhi_c)

        if k in broken_but_inside_PMT_in_inner_bottom_ring_GERDA:
            if PMT_configuration == 'LEGEND200':
                continue
            namephotocathode_lv = f"GERDA_PMTcathode_lv_{Num_PMTs}"
            namephotocathode = f"GERDA_PMTcathode_{Num_PMTs}"
            namesteelcone = f"GERDA_PMTcone_{Num_PMTs}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)
            namePMTair_lv = f"GERDA_PMTair_lv_{Num_PMTs}"
            namePMTair = f"GERDA_PMTair_{Num_PMTs}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"GERDA_PMTacryl_lv_{Num_PMTs}"
            namePMTacryl = f"GERDA_PMTacryl_{Num_PMTs}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)  
            namesteelbottom = f"GERDA_PMTbottom_{Num_PMTs}"
            namePMTborosilikat_lv = f"GERDA_PMTborosilikat_lv_{Num_PMTs}"
            namePMTborosilikat = f"GERDA_PMTborosilikat_{Num_PMTs}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg) 
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
        else:
            namephotocathode_lv = f"PMTcathode_lv_{PMT_ID[Working_PMTs]}"
            namephotocathode = f"PMTcathode_{PMT_ID[Working_PMTs]}"
            namesteelcone = f"PMTcone_{PMT_ID[Working_PMTs]}"
            namesteelbottom = f"PMTbottom_{PMT_ID[Working_PMTs]}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)

            namePMTair_lv = f"PMTair_lv_{PMT_ID[Working_PMTs]}"
            namePMTair = f"PMTair_{PMT_ID[Working_PMTs]}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"PMTacryl_lv_{PMT_ID[Working_PMTs]}"
            namePMTacryl = f"PMTacryl_{PMT_ID[Working_PMTs]}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)


            namePMTborosilikat_lv = f"PMTborosilikat_lv_{PMT_ID[Working_PMTs]}"
            namePMTborosilikat = f"PMTborosilikat_{PMT_ID[Working_PMTs]}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg)   
            g4.BorderSurface(f"AirAcrylSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, Acryl_pv, optical_border_air_acryl, reg)
            g4.BorderSurface(f"WaterAcrylSurface{PMT_ID[Working_PMTs]}", water_pv, Acryl_pv, optical_border_water_acryl, reg)
            g4.BorderSurface(f"AirBorosilicateSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, borosilikat_pv, optical_border_air_borosilicate, reg)
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
            g4.SkinSurface(f"PMTCathodeSkinSurface{PMT_ID[Working_PMTs]}", Photocathode_lv, optical_PMT_surface, reg)
            Working_PMTs += 1
        
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTConeOffset], PMTSteelCone_lv, namesteelcone, water_lv, reg)
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTBottomOffset], PMTSteelBottom_lv, namesteelbottom, water_lv, reg)
        
        dPhi_c += dPhi
        Num_PMTs += 1

    # Pillbox Bottom ring
    N_PMT_per_ring_bpb = 6  
    R_pos = ShieldingFootIR - 0.5*PMTSteelConeHeight - PhotocathodeHeightDifference - PMTSteelBottomHeight - ReflectiveFoilThickness
    dPhi = 2.0 * np.pi / N_PMT_per_ring_bpb 
    dPhi_c = dPhi

    broken_but_inside_PMT_in_pillbox_bottom_ring = np.array([])
    no_PMT_in_pillbox_bottom_ring = np.array([0,3])
    PMT_amount_pillbox_bottom_ring = N_PMT_per_ring_bpb - len(no_PMT_in_pillbox_bottom_ring)
    # PMTs
    for k in range(0, N_PMT_per_ring_bpb):
        if k in no_PMT_in_pillbox_bottom_ring:
            continue

        dPhi_c = dPhi *k
        xpos = R_pos * np.cos(dPhi_c)
        ypos = R_pos * np.sin(dPhi_c)


        if k in broken_but_inside_PMT_in_pillbox_bottom_ring:
            namephotocathode_lv = f"GERDA_PMTcathode_lv_{Num_PMTs}"
            namephotocathode = f"GERDA_PMTcathode_{Num_PMTs}"
            namesteelcone = f"GERDA_PMTcone_{Num_PMTs}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)
            namePMTair_lv = f"GERDA_PMTair_lv_{Num_PMTs}"
            namePMTair = f"GERDA_PMTair_{Num_PMTs}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"GERDA_PMTacryl_lv_{Num_PMTs}"
            namePMTacryl = f"GERDA_PMTacryl_{Num_PMTs}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)  
            namesteelbottom = f"GERDA_PMTbottom_{Num_PMTs}"
            namePMTborosilikat_lv = f"GERDA_PMTborosilikat_lv_{Num_PMTs}"
            namePMTborosilikat = f"GERDA_PMTborosilikat_{Num_PMTs}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg) 
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
        else:
            namephotocathode_lv = f"PMTcathode_lv_{PMT_ID[Working_PMTs]}"
            namephotocathode = f"PMTcathode_{PMT_ID[Working_PMTs]}"
            namesteelcone = f"PMTcone_{PMT_ID[Working_PMTs]}"
            namesteelbottom = f"PMTbottom_{PMT_ID[Working_PMTs]}"
            Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)

            namePMTair_lv = f"PMTair_lv_{PMT_ID[Working_PMTs]}"
            namePMTair = f"PMTair_{PMT_ID[Working_PMTs]}"
            PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
            PMTAir_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], PMTAir_lv, namePMTair, water_lv, reg) 
                    
            namePMTacryl_lv = f"PMTacryl_lv_{PMT_ID[Working_PMTs]}"
            namePMTacryl = f"PMTacryl_{PMT_ID[Working_PMTs]}"
            Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
            Acryl_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], Acryl_lv, namePMTacryl, water_lv, reg)


            namePMTborosilikat_lv = f"PMTborosilikat_lv_{PMT_ID[Working_PMTs]}"
            namePMTborosilikat = f"PMTborosilikat_{PMT_ID[Working_PMTs]}"
            borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
            borosilikat_pv = g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTCathodeOffset], borosilikat_lv, namePMTborosilikat, water_lv, reg)   
            g4.BorderSurface(f"AirAcrylSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, Acryl_pv, optical_border_air_acryl, reg)
            g4.BorderSurface(f"WaterAcrylSurface{PMT_ID[Working_PMTs]}", water_pv, Acryl_pv, optical_border_water_acryl, reg)
            g4.BorderSurface(f"AirBorosilicateSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, borosilikat_pv, optical_border_air_borosilicate, reg)
            PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
            g4.SkinSurface(f"PMTCathodeSkinSurface{PMT_ID[Working_PMTs]}", Photocathode_lv, optical_PMT_surface, reg)
            Working_PMTs += 1
            
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTConeOffset], PMTSteelCone_lv, namesteelcone, water_lv, reg)
        g4.PhysicalVolume([0,0,0], [xpos,ypos,PMTBottomOffset], PMTSteelBottom_lv, namesteelbottom, water_lv, reg)
        
        dPhi_c += dPhi
        Num_PMTs +=1


    NzPMT = 4  # number of PMT rings at wall
    N_PMT_per_ring = 10  # PMTs per ring

    distancetobottom = 200.0  # distance of pillbox wall PMTs to bottom
    R_pos = WaterRadius - ReflectiveFoilThickness - DistancePMTBaseTank - PMTSteelConeHeight*0.5 - PMTSteelBottomHeight 
    dPhi = 2.0 * np.pi / N_PMT_per_ring


    for i in [0, 1, 2, 3, 4]:
        if i == 0:
            # PMT row (2m)
            N_PMT_per_ring = 10
            dPhi = 2.0 * np.pi / N_PMT_per_ring
            no_PMT_in_wall_ring = np.array([])
            broken_but_inside_PMT_in_wall_ring_GERDA = np.array([4, 5, 6, 7])
            zpos = -2450.0
            dPhi_c = dPhi * 0.5
        elif i == 1:
            # PMT row (3.5 m)
            N_PMT_per_ring = 10
            dPhi = 2.0 * np.pi / N_PMT_per_ring
            no_PMT_in_wall_ring = np.array([])
            broken_but_inside_PMT_in_wall_ring_GERDA = np.array([4, 5])
            zpos = -950.0
            dPhi_c = 0.0
        elif i == 2:
            # PMT row (5 m)
            N_PMT_per_ring = 10
            dPhi = 2.0 * np.pi / N_PMT_per_ring
            no_PMT_in_wall_ring = np.array([])
            broken_but_inside_PMT_in_wall_ring_GERDA = np.array([0, 3])
            zpos = 550.0
            #zpos=-2450.0
            dPhi_c = dPhi * 0.5
        elif i == 3:
            # PMT row (6.5 m)
            N_PMT_per_ring = 10
            dPhi = 2.0 * np.pi / N_PMT_per_ring
            no_PMT_in_wall_ring = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9])
            brokManhole_heighten_but_inside_PMT_in_wall_ring_GERDA = np.array([])
            zpos = 2050.0
            dPhi_c = 0.0
        elif i == 4:
            # Pillbox wall mounted PMTs
            no_PMT_in_wall_ring = np.array([])
            broken_but_inside_PMT_in_wall_ring_GERDA = np.array([])
            N_PMT_per_ring = 6
            R_pos = ShieldingFootIR - ReflectiveFoilThickness - PMTSteelConeHeight*0.5 - PMTSteelBottomHeight
            zpos = -(WaterHeight / 2.0) + distancetobottom * 3
            dPhi = 2 * np.pi / N_PMT_per_ring
            dPhi_c = dPhi * 0.5

        
        x_rot_drehvol = 0
        y_rot_drehvol = (np.pi/2.0)
    

        # placing of PMTs in ring
        for k in range(0, N_PMT_per_ring):
            if k in no_PMT_in_wall_ring:
                continue

            if i == 0 or i==2 or i==4:
                dPhi_c = dPhi *k + 0.5*dPhi
            else:
                dPhi_c = dPhi *k

            xpos = (R_pos+ PhotocathodeHeightDifference) * np.cos(dPhi_c) 
            ypos = (R_pos+ PhotocathodeHeightDifference)  * np.sin(dPhi_c)

            if i == 0 or i==2 or i==4:
                z_rot_drehvol = -0.5*dPhi + N_PMT_per_ring*dPhi - k*dPhi
            else:
                z_rot_drehvol = N_PMT_per_ring*dPhi - k*dPhi


            # rotations around Z, Y and X
            rot_z = R.from_euler('z', z_rot_drehvol, degrees=False)
            rot_y = R.from_euler('y', y_rot_drehvol, degrees=False)
            rot_x = R.from_euler('x', x_rot_drehvol, degrees=False)

            combined_rotation = rot_y * rot_z * rot_x
            rot_matrix = combined_rotation.as_matrix()

            euler_angles = combined_rotation.as_euler('xyz', degrees=False)
            x_rot_global, y_rot_global, z_rot_global = euler_angles

            if k in broken_but_inside_PMT_in_wall_ring_GERDA:
                if PMT_configuration == 'LEGEND200':
                    continue
                namephotocathode_lv = f"GERDA_PMTcathode_lv_{Num_PMTs}"
                namephotocathode = f"GERDA_PMTcathode_{Num_PMTs}"
                namesteelcone = f"GERDA_PMTcone_{Num_PMTs}"
                Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)
                namePMTair_lv = f"GERDA_PMTair_lv_{Num_PMTs}"
                namePMTair = f"GERDA_PMTair_{Num_PMTs}"
                PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
                PMTAir_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], PMTAir_lv, namePMTair, water_lv, reg) 
                        
                namePMTacryl_lv = f"GERDA_PMTacryl_lv_{Num_PMTs}"
                namePMTacryl = f"GERDA_PMTacryl_{Num_PMTs}"
                Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
                Acryl_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], Acryl_lv, namePMTacryl, water_lv, reg)  
                namesteelbottom = f"GERDA_PMTbottom_{Num_PMTs}"
                namePMTborosilikat_lv = f"GERDA_PMTborosilikat_lv_{Num_PMTs}"
                namePMTborosilikat = f"GERDA_PMTborosilikat_{Num_PMTs}"
                borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
                borosilikat_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], borosilikat_lv, namePMTborosilikat, water_lv, reg) 
                PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
            else:
                namephotocathode_lv = f"PMTcathode_lv_{PMT_ID[Working_PMTs]}"
                namephotocathode = f"PMTcathode_{PMT_ID[Working_PMTs]}"
                namesteelcone = f"PMTcone_{PMT_ID[Working_PMTs]}"
                namesteelbottom = f"PMTbottom_{PMT_ID[Working_PMTs]}"
                Photocathode_lv = g4.LogicalVolume(Photocathode,cathodeAl,namephotocathode_lv, reg)

                namePMTair_lv = f"PMTair_lv_{PMT_ID[Working_PMTs]}"
                namePMTair = f"PMTair_{PMT_ID[Working_PMTs]}"
                PMTAir_lv = g4.LogicalVolume(PMTAir,PMT_air_material,namePMTair_lv, reg)
                PMTAir_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], PMTAir_lv, namePMTair, water_lv, reg) 
                        
                namePMTacryl_lv = f"PMTacryl_lv_{PMT_ID[Working_PMTs]}"
                namePMTacryl = f"PMTacryl_{PMT_ID[Working_PMTs]}"
                Acryl_lv = g4.LogicalVolume(Acryl,acryl_material,namePMTacryl_lv, reg)
                Acryl_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], Acryl_lv, namePMTacryl, water_lv, reg)

                namePMTborosilikat_lv = f"PMTborosilikat_lv_{PMT_ID[Working_PMTs]}"
                namePMTborosilikat = f"PMTborosilikat_{PMT_ID[Working_PMTs]}"
                borosilikat_lv = g4.LogicalVolume(PMTBorosilikatGlass,borosilicate_material,namePMTborosilikat_lv, reg)
                borosilikat_pv = g4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos, ypos, zpos], borosilikat_lv, namePMTborosilikat, water_lv, reg)   
                g4.BorderSurface(f"AirAcrylSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, Acryl_pv, optical_border_air_acryl, reg)
                g4.BorderSurface(f"WaterAcrylSurface{PMT_ID[Working_PMTs]}", water_pv, Acryl_pv, optical_border_water_acryl, reg)
                g4.BorderSurface(f"AirBorosilicateSurface{PMT_ID[Working_PMTs]}", PMTAir_pv, borosilikat_pv, optical_border_air_borosilicate, reg)
                PMT_cathode_pv = g4.PhysicalVolume([0,0,0], [0,0,0], Photocathode_lv, namephotocathode, borosilikat_lv, reg) 
                g4.SkinSurface(f"PMTCathodeSkinSurface{PMT_ID[Working_PMTs]}", Photocathode_lv, optical_PMT_surface, reg)
                Working_PMTs += 1
    
            # PMT cone
            xpos_cone = (R_pos + PMTSteelConeHeight * 0.25) * np.cos(dPhi_c)
            ypos_cone = (R_pos + PMTSteelConeHeight * 0.25) * np.sin(dPhi_c)
            pyg4.geant4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos_cone, ypos_cone, zpos], PMTSteelCone_lv, namesteelcone, water_lv, reg)
            
            # PMT bottom
            xpos_bottom = (R_pos + PMTSteelConeHeight *0.5 + PMTSteelBottomHeight * 0.5) * np.cos(dPhi_c)
            ypos_bottom = (R_pos + PMTSteelConeHeight*0.5+ PMTSteelBottomHeight * 0.5) * np.sin(dPhi_c)
            
            pyg4.geant4.PhysicalVolume([x_rot_global, y_rot_global, z_rot_global], [xpos_bottom, ypos_bottom, zpos], PMTSteelBottom_lv, namesteelbottom, water_lv, reg)
        
            dPhi_c += dPhi
            Num_PMTs += 1