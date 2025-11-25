from __future__ import annotations

import itertools
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from dbetto import TextDB
from pyg4ometry import geant4 as g4
from pygeomtools import RemageDetectorInfo

from . import core, top


def place_fiber_modules(
    fiber_metadata: TextDB,
    b: core.InstrumentationData,
    use_detailed_fiber_model: bool = False,
) -> None:
    """Construct LEGEND-200 liquid argon instrumentation.

    Parameters
    ----------
    use_detailed_fiber_model
        Switch between an implementation of single fibers (“detailed”) or
        slabs of fiber material (“segmented”).
    """
    # Unroll the provided metadata into a structure better suited for the next steps.
    # The geometry here is based on physical modules and not on channels.
    modules = {}
    ch_map = b.channelmap.map("system", unique=False).get("spms", {})
    for ch in ch_map.values():
        mod = modules.get(ch.location.fiber)

        if mod is None:
            # initialize a new module object if we don't have one yet.
            mod = FiberModuleData(
                name=ch.location.fiber,
                barrel=fiber_metadata[ch.location.fiber].type,
                tpb_thickness=fiber_metadata[ch.location.fiber].geometry.tpb.thickness_in_nm,
            )
            modules[ch.location.fiber] = mod

        assert getattr(mod, f"channel_{ch.location.position}_name") is None
        setattr(mod, f"channel_{ch.location.position}_name", ch.name)
        setattr(mod, f"channel_{ch.location.position}_rawid", ch.daq.rawid)

    factory = ModuleFactorySingleFibers if use_detailed_fiber_model else ModuleFactorySegment

    z_displacement_fiber_assembly = (
        # avoid the overlap of the top SiPMs with the top plate.
        b.top_plate_z_pos
        - top.TOP_PLATE_THICKNESS
        - ModuleFactoryBase.SIPM_HEIGHT
        - ModuleFactoryBase.SIPM_OUTER_EXTRA
        - 25  # the SiPMs are below the top plate. dummy value (guessed from photos).
    )

    # note: actually the radius is only 150mm and another short straight segment of 60mm is following after
    # the bend. to simplify things here, those two are combined to one bent shape, to have at least the same
    # covered solid angle.
    # we do not add the full straight part here (this would make the radius to large), but just a tiny delta
    # that makes the bottom of the OB cover the whole area between the straight OB and IB fibers.
    ob_radius_delta = 10
    ob_radius_mm = 155 + ob_radius_delta
    ob_inner_straight = 60 - ob_radius_delta
    # this OB fiber length is derived from measurements in the CAD model, und might not be totally correct.
    ob_fiber_length_mm = 1630
    ob_factory = factory(
        radius_mm=290 - 6,
        fiber_length_mm=ob_fiber_length_mm - math.pi * ob_radius_mm / 2 - ob_inner_straight,
        bend_radius_mm=ob_radius_mm,
        fiber_count_per_module=81,
        number_of_modules=20,
        zero_angle_module=_module_name_to_num("OB015016"),
        z_displacement_mm=z_displacement_fiber_assembly,
        b=b,
        barrel="outer",
    )

    ib_fiber_length_mm = 1400
    ib_delta_z = 35
    ib_factory = factory(
        radius_mm=126,
        fiber_length_mm=ib_fiber_length_mm,
        bend_radius_mm=None,
        fiber_count_per_module=81,
        number_of_modules=9,
        zero_angle_module=_module_name_to_num("IB013014"),
        z_displacement_mm=z_displacement_fiber_assembly - ib_delta_z,
        b=b,
        barrel="inner",
    )

    for mod in modules.values():
        if mod.barrel == "outer":
            ob_factory.create_module(mod)
        if mod.barrel == "inner":
            ib_factory.create_module(mod)

    if any(mod.barrel == "outer" for mod in modules.values()):
        create_fiber_support_outer(b, b.top_plate_z_pos - 730)
    if any(mod.barrel == "inner" for mod in modules.values()):
        create_fiber_support_inner(b, b.top_plate_z_pos - 730 - ib_delta_z)


def _module_name_to_num(mod_name: str) -> int:
    m0, m1 = int(mod_name[2:5]), int(mod_name[5:8]) - 1
    if m0 != m1:
        msg = f"Invalid fiber module name {mod_name}"
        raise ValueError(msg)
    return int((m0 - 1) / 2)


@dataclass
class FiberModuleData:
    barrel: str
    name: str
    tpb_thickness: float
    channel_top_name: str | None = None
    channel_bottom_name: str | None = None
    channel_top_rawid: str | None = None
    channel_bottom_rawid: str | None = None


class ModuleFactoryBase(ABC):
    FIBER_DIM = 1  # mm
    FIBER_THICKNESS_CL1 = 0.04 * FIBER_DIM  # (BCF-91A document)
    FIBER_THICKNESS_CL2 = 0.02 * FIBER_DIM  # (BCF-91A document)

    # The "SiPM" is only a dummy implementation here, the fibers are not bent to the real geometry of 3x3
    # fibers coupled to one physical SiPM!
    SIPM_HEIGHT = 1  # mm, dummy
    # There is a LAr gap between the fiber end and SiPM:
    SIPM_GAP = 0.05  # mm
    # Because of this gap, we cannot use surfaces between fibers and the "SiPM". To stop stray light from
    # entering the "SiPM" from the other sides, we add an outer envelope that blocks light. On each outer
    # side of the SiPM volume, an additional solid of size "outer extra" is added:
    SIPM_OUTER_EXTRA = 0.2  # mm
    # To also stop stray light from directions more close to the fibers, the envelope extends a bit more
    # along the fibers:
    SIPM_OVERLAP = 0.3  # mm
    SIPM_GAP_SIDE = 0.01  # mm, for fitting problems with round "SiPMs" and square fibers.

    ANGLE_SAFETY = 1e-9  # rad

    def __init__(
        self,
        radius_mm: float,
        fiber_length_mm: float,
        fiber_count_per_module: int,
        bend_radius_mm: float | None,
        number_of_modules: int,
        zero_angle_module: int,
        z_displacement_mm: float,
        barrel: str,
        b: core.InstrumentationData,
    ):
        """
        Create a fiber module factory.

        Parameters
        ----------
        radius_mm
            radius of the fiber barrel
        fiber_length_mm
            length of the straight section of this fiber module
        fiber_count_per_module
            number of single fibers per module
        bend_radius_mm
            radius of the bottom bend, or None if the fibers are not bent at the bottom end.
        number_of_modules
            number of modules that cover the full circle
        zero_angle_module
            module number of the module with a zero angle in polar coordinates (at the center of the module).
        z_displacement_mm
            displacement of the top of the fiber barrel, relative to the global zero point.
        barrel
            barrel name
        """
        self.radius = radius_mm
        self.fiber_length = fiber_length_mm
        self.fiber_count_per_module = fiber_count_per_module
        self.bend_radius_mm = bend_radius_mm
        self.number_of_modules = number_of_modules
        self.zero_angle_module = zero_angle_module
        self.z_displacement = z_displacement_mm
        self.barrel = barrel
        self.b = b

    def _cached_sipm_volumes(self) -> None:
        """Creates (dummy) SiPM volumes for use at the top/bottom of straight fiber sections."""
        v_suffix = f"_r{self.radius}_nmod{self.number_of_modules}"
        v_name = f"sipm{v_suffix}"
        if v_name in self.b.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square
        fiber_segment = 2 * math.pi / self.number_of_modules

        sipm = g4.solid.Tubs(
            v_name,
            self.radius - sipm_dim / 2,
            self.radius + sipm_dim / 2,
            self.SIPM_HEIGHT,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        self.sipm_lv = g4.LogicalVolume(sipm, self.b.materials.metal_silicon, v_name, self.b.registry)

        sipm_outer1 = g4.solid.Tubs(
            f"sipm_outer1{v_suffix}",
            self.radius - sipm_dim / 2 - self.SIPM_OUTER_EXTRA,
            self.radius + sipm_dim / 2 + self.SIPM_OUTER_EXTRA,
            self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Tubs(
            f"sipm_outer2{v_suffix}",
            self.radius - sipm_dim / 2 - 1e-9,
            self.radius + sipm_dim / 2 + 1e-9,
            self.SIPM_HEIGHT + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        sipm_outer_top = g4.solid.Subtraction(
            f"sipm_outer_top{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, -self.SIPM_OUTER_EXTRA / 2]],
            self.b.registry,
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, +self.SIPM_OUTER_EXTRA / 2]],
            self.b.registry,
        )
        self.sipm_outer_top_lv = g4.LogicalVolume(
            sipm_outer_top,
            self.b.materials.tetratex,
            f"sipm_outer_top{v_suffix}",
            self.b.registry,
        )
        self.sipm_outer_bottom_lv = g4.LogicalVolume(
            sipm_outer_bottom,
            self.b.materials.tetratex,
            f"sipm_outer_bottom{v_suffix}",
            self.b.registry,
        )
        # note: this is tetratex here, but do not add a surface here. I do not want to have the reflectivity
        # on the inside.

        # TODO: implement partial modules with end envelopes for SiPM.
        # sipm_outer_end = g4.solid.Box(
        #    f"sipm_outer_end{v_suffix}",
        #    sipm_dim + self.SIPM_OUTER_EXTRA * 2,
        #    self.SIPM_OUTER_EXTRA,
        #    self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
        #    self.b.registry,
        # )
        # g4.LogicalVolume(
        #    sipm_outer_end,
        #    self.b.materials.tetratex,
        #    f"sipm_outer_end{v_suffix}",
        #    self.b.registry,
        # )

    def start_angle(self, module_num: int) -> float:
        return 2 * math.pi / self.number_of_modules * (module_num - self.zero_angle_module - 0.5)

    @abstractmethod
    def create_module(self, mod: FiberModuleData) -> None:
        raise NotImplementedError()

    def _create_sipm(
        self,
        module_num: int,
        fibers: list[g4.PhysicalVolume],
        is_top: bool,
        sipm_name: str,
        sipm_detector_id: int,
        z_displacement_straight: float,
    ) -> None:
        """Creates a (dummy) SiPM physical volume for use at the top/bottom of straight fiber sections."""
        z = +self.fiber_length / 2 + self.SIPM_HEIGHT / 2 + self.SIPM_GAP  # add small gap
        z_outer = (
            z + self.SIPM_OUTER_EXTRA / 2 - self.SIPM_OVERLAP / 2 - self.SIPM_GAP + 1e-9
            if is_top
            else -z - self.SIPM_OUTER_EXTRA / 2 + self.SIPM_OVERLAP / 2 + self.SIPM_GAP - 1e-9
        )
        z = z if is_top else -z
        z += z_displacement_straight
        z_outer += z_displacement_straight
        sipm_pv = g4.PhysicalVolume(
            [0, 0, -self.start_angle(module_num)],
            [0, 0, z],
            self.sipm_lv,
            sipm_name,
            self.b.mother_lv,
            self.b.registry,
        )
        sipm_pv.set_pygeom_active_detector(RemageDetectorInfo("optical", sipm_detector_id))
        # Add border surface to mother volume.
        g4.BorderSurface(
            f"bsurface_lar_{sipm_name}",
            self.b.mother_pv,
            sipm_pv,
            self.b.materials.surfaces.to_sipm_silicon(
                self.b.runtime_config,
                sipm_name,
            ),
            self.b.registry,
        )

        g4.PhysicalVolume(
            [0, 0, -self.start_angle(module_num)],
            [0, 0, z_outer],
            self.sipm_outer_top_lv if is_top else self.sipm_outer_bottom_lv,
            f"{sipm_name}_wrap",
            self.b.mother_lv,
            self.b.registry,
        )

    def _add_tpb_surfaces(self, fiber_pvs: list[g4.PhysicalVolume]):
        # between LAr and PEN we need a surface in both directions.
        for tpb_pv in fiber_pvs:
            g4.BorderSurface(
                "bsurface_lar_tpb_" + tpb_pv.name,
                self.b.mother_pv,
                tpb_pv,
                self.b.materials.surfaces.lar_to_tpb,
                self.b.registry,
            )
            g4.BorderSurface(
                "bsurface_tpb_lar_" + tpb_pv.name,
                tpb_pv,
                self.b.mother_pv,
                self.b.materials.surfaces.lar_to_tpb,
                self.b.registry,
            )


class ModuleFactorySingleFibers(ModuleFactoryBase):
    # for bent detailed fibers, the fibers would overlap a lot near the bottom SiPMs. To avoid
    # this, use a staggered design of the fibers.
    ALLOWED_DELTA_LENGTHS = np.arange(-4, 5) * 1.9

    def _cached_sipm_volumes_bend(self) -> None:
        """Creates (dummy) SiPM volumes for use at the bottom of bent fiber sections."""
        v_suffix = "_bend"  # this dummy SiPM is for one single fiber, so we do not need any attributes here.
        v_name = f"sipm{v_suffix}"
        if v_name in self.b.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square

        self.sipm_bend = g4.solid.Box(
            v_name,
            self.SIPM_HEIGHT,
            sipm_dim,
            sipm_dim,
            self.b.registry,
            "mm",
        )
        self.sipm_lv_bend = g4.LogicalVolume(
            self.sipm_bend, self.b.materials.metal_silicon, v_name, self.b.registry
        )

        sipm_outer1 = g4.solid.Box(
            f"sipm_outer1{v_suffix}",
            self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            self.b.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Box(
            f"sipm_outer2{v_suffix}",
            self.SIPM_HEIGHT + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            sipm_dim + 1e-8,
            sipm_dim + 1e-8,
            self.b.registry,
            "mm",
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [+self.SIPM_OUTER_EXTRA / 2, 0, 0]],
            self.b.registry,
        )
        self.sipm_outer_bottom_lv_bend = g4.LogicalVolume(
            sipm_outer_bottom,
            self.b.materials.tetratex,
            f"sipm_outer_bottom{v_suffix}",
            self.b.registry,
        )
        # note: this is tetratex here, but do not add a surface here. I do not want to have the reflectivity
        # on the inside.

    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}_b{self.bend_radius_mm}"
        if f"fiber_cl2{v_suffix}" in self.b.registry.solidDict:
            return

        fibers_to_gen = [(v_suffix, self.fiber_length)]
        if self.bend_radius_mm is not None:
            for delta_length in self.ALLOWED_DELTA_LENGTHS:
                if delta_length == 0:
                    continue
                fibers_to_gen += [
                    (
                        f"_l{self.fiber_length + delta_length}_b{self.bend_radius_mm}",
                        self.fiber_length + delta_length,
                    )
                ]

        # create solids
        fiber_cl2 = {}
        fiber_cl1 = {}
        fiber_core = {}
        dim_cl1 = self.FIBER_DIM - 2 * self.FIBER_THICKNESS_CL2
        dim_core = self.FIBER_DIM - 2 * (self.FIBER_THICKNESS_CL1 + self.FIBER_THICKNESS_CL2)
        for [fiber_name, fiber_length] in fibers_to_gen:
            fiber_cl2[fiber_length] = g4.solid.Box(
                f"fiber_cl2{fiber_name}",
                self.FIBER_DIM,
                self.FIBER_DIM,
                fiber_length,
                self.b.registry,
                "mm",
            )
            fiber_cl1[fiber_length] = g4.solid.Box(
                f"fiber_cl1{fiber_name}",
                dim_cl1,
                dim_cl1,
                fiber_length,
                self.b.registry,
                "mm",
            )
            fiber_core[fiber_length] = g4.solid.Box(
                f"fiber_core{fiber_name}",
                dim_core,
                dim_core,
                fiber_length,
                self.b.registry,
                "mm",
            )
        if self.bend_radius_mm is not None:
            fiber_cl2_bend = g4.solid.Tubs(
                f"fiber_cl2_bend{v_suffix}",
                self.bend_radius_mm - self.FIBER_DIM / 2,
                self.bend_radius_mm + self.FIBER_DIM / 2,
                self.FIBER_DIM,
                0,
                math.pi / 2,
                self.b.registry,
                "mm",
            )
            fiber_cl1_bend = g4.solid.Tubs(
                f"fiber_cl1_bend{v_suffix}",
                self.bend_radius_mm - dim_cl1 / 2,
                self.bend_radius_mm + dim_cl1 / 2,
                dim_cl1,
                0,
                math.pi / 2,
                self.b.registry,
                "mm",
            )
            fiber_core_bend = g4.solid.Tubs(
                f"fiber_core_bend{v_suffix}",
                self.bend_radius_mm - dim_core / 2,
                self.bend_radius_mm + dim_core / 2,
                dim_core,
                0,
                math.pi / 2,
                self.b.registry,
                "mm",
            )

        self.fiber_cl2_lv = {}
        fiber_cl1_lv = {}
        fiber_core_lv = {}
        for [fiber_name, fiber_length] in fibers_to_gen:
            self.fiber_cl2_lv[fiber_length] = g4.LogicalVolume(
                fiber_cl2[fiber_length], self.b.materials.pmma_out, f"fiber_cl2{fiber_name}", self.b.registry
            )
            self.fiber_cl2_lv[fiber_length].pygeom_color_rgba = False
            fiber_cl1_lv[fiber_length] = g4.LogicalVolume(
                fiber_cl1[fiber_length], self.b.materials.pmma, f"fiber_cl1{fiber_name}", self.b.registry
            )
            fiber_cl1_lv[fiber_length].pygeom_color_rgba = False
            fiber_core_lv[fiber_length] = g4.LogicalVolume(
                fiber_core[fiber_length],
                self.b.materials.ps_fibers,
                f"fiber_core{fiber_name}",
                self.b.registry,
            )
            fiber_core_lv[fiber_length].pygeom_color_rgba = False
        if self.bend_radius_mm is not None:
            self.fiber_cl2_bend_lv = g4.LogicalVolume(
                fiber_cl2_bend, self.b.materials.pmma_out, f"fiber_cl2_bend{v_suffix}", self.b.registry
            )
            self.fiber_cl2_bend_lv.pygeom_color_rgba = False
            fiber_cl1_bend_lv = g4.LogicalVolume(
                fiber_cl1_bend, self.b.materials.pmma, f"fiber_cl1_bend{v_suffix}", self.b.registry
            )
            fiber_cl1_bend_lv.pygeom_color_rgba = False
            fiber_core_bend_lv = g4.LogicalVolume(
                fiber_core_bend, self.b.materials.ps_fibers, f"fiber_core_bend{v_suffix}", self.b.registry
            )
            fiber_core_bend_lv.pygeom_color_rgba = False

        for [fiber_name, fiber_length] in fibers_to_gen:
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_cl1_lv[fiber_length],
                f"fiber_{self.barrel}_cladding1{fiber_name}",
                self.fiber_cl2_lv[fiber_length],
                self.b.registry,
            )
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_core_lv[fiber_length],
                f"fiber_{self.barrel}_core{fiber_name}",
                fiber_cl1_lv[fiber_length],
                self.b.registry,
            )
        if self.bend_radius_mm is not None:
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_cl1_bend_lv,
                f"fiber_{self.barrel}_cladding1_bend{v_suffix}",
                self.fiber_cl2_bend_lv,
                self.b.registry,
            )
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_core_bend_lv,
                f"fiber_{self.barrel}_core_bend{v_suffix}",
                fiber_cl1_bend_lv,
                self.b.registry,
            )

    def _cached_tpb_coating_volume(
        self, tpb_thickness_nm: float, bend: bool = False, delta_length: int = 0
    ) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        if delta_length != 0 and bend:
            msg = "creating a bent volume with delta_length!=0 is not possible"
            raise ValueError(msg)
        if delta_length not in self.ALLOWED_DELTA_LENGTHS:
            msg = f"invalid delta length {delta_length}"
            raise ValueError(msg)
        fiber_length = self.fiber_length + delta_length

        v_suffix = f"{'_bend' if bend else ''}_l{fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.b.registry.solidDict:
            return self.b.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + 2 * tpb_thickness_nm / 1e6
        if not bend:
            coating = g4.solid.Box(v_name, coating_dim, coating_dim, fiber_length, self.b.registry, "mm")
            inner_lv = self.fiber_cl2_lv[fiber_length]
        else:
            coating = g4.solid.Tubs(
                v_name,
                self.bend_radius_mm - coating_dim / 2,
                self.bend_radius_mm + coating_dim / 2,
                coating_dim,
                0,
                math.pi / 2,
                self.b.registry,
                "mm",
            )
            inner_lv = self.fiber_cl2_bend_lv
        coating_lv = g4.LogicalVolume(coating, self.b.materials.tpb_on_fibers, v_name, self.b.registry)
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            inner_lv,
            f"fiber_{self.barrel}_cladding2{v_suffix}",
            coating_lv,
            self.b.registry,
        )

        coating_lv.pygeom_color_rgba = [0, 1, 0.165, 0.07]  # 520 nm

        return coating_lv

    def create_module(self, mod: FiberModuleData) -> None:
        assert mod.barrel == self.barrel
        module_num = _module_name_to_num(mod.name)
        if module_num < 0 or module_num >= self.number_of_modules:
            msg = f"invalid module number {module_num} for a maximum of {self.number_of_modules}-1 modules."
            raise ValueError(msg)

        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        if self.bend_radius_mm is not None:
            self._cached_sipm_volumes_bend()
            coating_lv_bend = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=True)

        start_angle = self.start_angle(module_num)

        z_displacement_straight = self.z_displacement - self.fiber_length / 2

        fibers = []
        sipm_transforms = []  # transforms for SiPMs in the MultiUnion structure.
        for n in range(self.fiber_count_per_module):
            delta_length = 0
            if self.bend_radius_mm is not None:
                # for bent detailed fibers, the fibers would overlap a lot near the bottom SiPMs. To avoid
                # this, use a staggered design of the fibers. This is certainly not the best/most
                # space-efficient packing for squares, but it is simple to implement. From above the solid
                # angle coverage is the same, but there are some holes between the fibers in the lower parts
                # if a photon arrives from the 'right' direction...
                delta_length = self.ALLOWED_DELTA_LENGTHS[n % len(self.ALLOWED_DELTA_LENGTHS)]

            coating_lv = self._cached_tpb_coating_volume(
                mod.tpb_thickness, bend=False, delta_length=delta_length
            )

            th = start_angle + 2 * math.pi / self.number_of_modules / self.fiber_count_per_module * (n + 0.5)
            x = self.radius * math.cos(th)
            y = self.radius * math.sin(th)
            fibers.append(
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [x, y, z_displacement_straight - delta_length / 2],
                    coating_lv,
                    f"fiber_{self.barrel}_tpb_{mod.name}_{n}",
                    self.b.mother_lv,
                    self.b.registry,
                )
            )
            if self.bend_radius_mm is not None:
                x2 = (self.radius - self.bend_radius_mm) * math.cos(th)
                y2 = (self.radius - self.bend_radius_mm) * math.sin(th)
                fibers.append(
                    g4.PhysicalVolume(
                        # this is an extrinsic rotation of pi/2 around X and -th around Z; expressed in intrinsic euler angles.
                        [math.pi / 2, th, 0],
                        [x2, y2, z_displacement_straight - self.fiber_length / 2 - delta_length],
                        coating_lv_bend,
                        f"fiber_{self.barrel}_tpb_{mod.name}_bend_{n}",
                        self.b.mother_lv,
                        self.b.registry,
                    )
                )

                # add per-fiber SiPMs (I do not know any other way...)
                sipm_placement_r = self.radius - self.bend_radius_mm - self.SIPM_GAP - self.SIPM_HEIGHT / 2
                x2 = sipm_placement_r * math.cos(th)
                y2 = sipm_placement_r * math.sin(th)
                z = z_displacement_straight - self.fiber_length / 2 - delta_length - self.bend_radius_mm
                sipm_transforms.append([[0, 0, th], [x2, y2, z]])

                sipm_placement_outer_r = (
                    sipm_placement_r
                    - self.SIPM_OUTER_EXTRA / 2
                    + self.SIPM_OVERLAP / 2
                    - self.SIPM_GAP
                    - 1e-9
                )
                x2 = sipm_placement_outer_r * math.cos(th)
                y2 = sipm_placement_outer_r * math.sin(th)
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [x2, y2, z],
                    self.sipm_outer_bottom_lv_bend,
                    f"{mod.channel_bottom_name}_{n}_wrap",
                    self.b.mother_lv,
                    self.b.registry,
                )

        self._add_tpb_surfaces(fibers)

        # create SiPMs and attach to fibers
        self._create_sipm(
            module_num,
            fibers,
            True,
            mod.channel_top_name,
            mod.channel_top_rawid,
            z_displacement_straight,
        )
        if self.bend_radius_mm is None:
            self._create_sipm(
                module_num,
                fibers,
                False,
                mod.channel_bottom_name,
                mod.channel_bottom_rawid,
                z_displacement_straight,
            )

        if self.bend_radius_mm is not None:
            sipm_vols = [self.sipm_bend] * len(sipm_transforms)
            sipm_mu = g4.solid.MultiUnion(
                mod.channel_bottom_name, sipm_vols, sipm_transforms, self.b.registry
            )
            sipm_mu_lv = g4.LogicalVolume(
                sipm_mu, self.b.materials.metal_silicon, mod.channel_bottom_name, self.b.registry
            )
            sipm_pv = g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                sipm_mu_lv,
                mod.channel_bottom_name,
                self.b.mother_lv,
                self.b.registry,
            )
            sipm_pv.set_pygeom_active_detector(RemageDetectorInfo("optical", mod.channel_bottom_rawid))
            # Add border surface to mother volume.
            g4.BorderSurface(
                f"bsurface_lar_{mod.channel_bottom_name}",
                self.b.mother_pv,
                sipm_pv,
                self.b.materials.surfaces.to_sipm_silicon(
                    self.b.runtime_config,
                    mod.channel_bottom_name,
                ),
                self.b.registry,
            )


class ModuleFactorySegment(ModuleFactoryBase):
    def _get_bend_polycone(
        self, inner_r: float, outer_r: float
    ) -> tuple[np.typing.ArrayLike, np.typing.ArrayLike]:
        """In the segmented model, there is no fundamental shape for the fiber bent available, so we
        use a polycone as a replacement.
        """
        delta_r_mm = (outer_r - inner_r) / 2
        bend_r_outer = self.bend_radius_mm + delta_r_mm
        bend_r_inner = self.bend_radius_mm - delta_r_mm

        angles = np.linspace(0, np.pi / 2, 100)
        z1 = bend_r_outer * np.sin(angles)
        r1 = bend_r_outer * np.cos(angles)
        z2 = bend_r_inner * np.sin(angles)
        r2 = bend_r_inner * np.cos(angles)

        z = self.bend_radius_mm - np.concatenate((z1, np.flip(z2)))
        # offset by the radius at the inner end of the bend.
        r = (outer_r - bend_r_outer) + np.concatenate((r1, np.flip(r2)))

        return z, r

    def _cached_sipm_volumes_bend(self) -> None:
        """Creates (dummy) SiPM volumes for use at the bottom of bent fiber sections."""
        v_suffix = f"_bend{self.bend_radius_mm}_r{self.radius}_nmod{self.number_of_modules}"
        v_name = f"sipm{v_suffix}"
        if v_name in self.b.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square
        fiber_segment = 2 * math.pi / self.number_of_modules
        # radius of the inner circle at the bottom, already including the small gap between fibers and SiPM.
        inner_radius = self.radius - self.bend_radius_mm - self.SIPM_GAP

        sipm = g4.solid.Tubs(
            v_name,
            inner_radius - self.SIPM_HEIGHT,
            inner_radius,
            sipm_dim,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        self.sipm_lv_bend = g4.LogicalVolume(sipm, self.b.materials.metal_silicon, v_name, self.b.registry)

        sipm_outer1 = g4.solid.Tubs(
            f"sipm_outer1{v_suffix}",
            inner_radius - self.SIPM_HEIGHT - self.SIPM_OUTER_EXTRA,
            inner_radius + self.SIPM_OVERLAP,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Tubs(
            f"sipm_outer2{v_suffix}",
            inner_radius - self.SIPM_HEIGHT - 1e-9,
            inner_radius + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            sipm_dim + 1e-9,
            self.ANGLE_SAFETY,
            fiber_segment - self.ANGLE_SAFETY,
            self.b.registry,
            "mm",
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, 0]],
            self.b.registry,
        )
        self.sipm_outer_bottom_lv_bend = g4.LogicalVolume(
            sipm_outer_bottom,
            self.b.materials.tetratex,
            f"sipm_outer_bottom{v_suffix}",
            self.b.registry,
        )
        # note: this is tetratex here, but do not add a surface here. I do not want to have the reflectivity
        # on the inside.

    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}"
        if f"fiber_cl2{v_suffix}" in self.b.registry.solidDict:
            return

        # create solids
        angle = 2 * np.pi / self.number_of_modules
        dim_cl2 = self.FIBER_DIM
        fiber_cl2 = g4.solid.Tubs(
            f"fiber_cl2{v_suffix}",
            self.radius - dim_cl2 / 2,
            self.radius + dim_cl2 / 2,
            self.fiber_length,
            0,
            angle,
            self.b.registry,
            "mm",
        )
        dim_cl1 = self.FIBER_DIM - 2 * self.FIBER_THICKNESS_CL1
        fiber_cl1 = g4.solid.Tubs(
            f"fiber_cl1{v_suffix}",
            self.radius - dim_cl1 / 2,
            self.radius + dim_cl1 / 2,
            self.fiber_length,
            0,
            angle,
            self.b.registry,
            "mm",
        )
        dim_core = self.FIBER_DIM - 2 * (self.FIBER_THICKNESS_CL1 + self.FIBER_THICKNESS_CL2)
        fiber_core = g4.solid.Tubs(
            f"fiber_core{v_suffix}",
            self.radius - dim_core / 2,
            self.radius + dim_core / 2,
            self.fiber_length,
            0,
            angle,
            self.b.registry,
            "mm",
        )
        if self.bend_radius_mm is not None:
            z, r = self._get_bend_polycone(self.radius - dim_cl2 / 2, self.radius + dim_cl2 / 2)
            fiber_cl2_bend = g4.solid.GenericPolycone(
                f"fiber_cl2_bend{v_suffix}", 0, angle, r, z, self.b.registry, "mm"
            )
            z, r = self._get_bend_polycone(self.radius - dim_cl1 / 2, self.radius + dim_cl1 / 2)
            fiber_cl1_bend = g4.solid.GenericPolycone(
                f"fiber_cl1_bend{v_suffix}", 0, angle, r, z, self.b.registry, "mm"
            )
            z, r = self._get_bend_polycone(self.radius - dim_core / 2, self.radius + dim_core / 2)
            fiber_core_bend = g4.solid.GenericPolycone(
                f"fiber_core_bend{v_suffix}", 0, angle, r, z, self.b.registry, "mm"
            )

        self.fiber_cl2_lv = g4.LogicalVolume(
            fiber_cl2, self.b.materials.pmma_out, f"fiber_cl2{v_suffix}", self.b.registry
        )
        self.fiber_cl2_lv.pygeom_color_rgba = False
        fiber_cl1_lv = g4.LogicalVolume(
            fiber_cl1, self.b.materials.pmma, f"fiber_cl1{v_suffix}", self.b.registry
        )
        fiber_cl1_lv.pygeom_color_rgba = False
        fiber_core_lv = g4.LogicalVolume(
            fiber_core, self.b.materials.ps_fibers, f"fiber_core{v_suffix}", self.b.registry
        )
        fiber_core_lv.pygeom_color_rgba = False
        if self.bend_radius_mm is not None:
            self.fiber_cl2_bend_lv = g4.LogicalVolume(
                fiber_cl2_bend, self.b.materials.pmma_out, f"fiber_cl2_bend{v_suffix}", self.b.registry
            )
            self.fiber_cl2_bend_lv.pygeom_color_rgba = False
            fiber_cl1_bend_lv = g4.LogicalVolume(
                fiber_cl1_bend, self.b.materials.pmma, f"fiber_cl1_bend{v_suffix}", self.b.registry
            )
            fiber_cl1_bend_lv.pygeom_color_rgba = False
            fiber_core_bend_lv = g4.LogicalVolume(
                fiber_core_bend, self.b.materials.ps_fibers, f"fiber_core_bend{v_suffix}", self.b.registry
            )
            fiber_core_bend_lv.pygeom_color_rgba = False

        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            fiber_cl1_lv,
            f"fiber_{self.barrel}_cladding1{v_suffix}",
            self.fiber_cl2_lv,
            self.b.registry,
        )
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            fiber_core_lv,
            f"fiber_{self.barrel}_core{v_suffix}",
            fiber_cl1_lv,
            self.b.registry,
        )
        if self.bend_radius_mm is not None:
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_cl1_bend_lv,
                f"fiber_{self.barrel}_cladding1_bend{v_suffix}",
                self.fiber_cl2_bend_lv,
                self.b.registry,
            )
            g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                fiber_core_bend_lv,
                f"fiber_{self.barrel}_core_bend{v_suffix}",
                fiber_cl1_bend_lv,
                self.b.registry,
            )

    def _cached_tpb_coating_volume(self, tpb_thickness_nm: float, bend: bool = False) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        v_suffix = f"{'_bend' if bend else ''}_l{self.fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.b.registry.solidDict:
            return self.b.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + 2 * tpb_thickness_nm / 1e6
        if not bend:
            coating = g4.solid.Tubs(
                v_name,
                self.radius - coating_dim / 2,
                self.radius + coating_dim / 2,
                self.fiber_length,
                0,
                2 * math.pi / self.number_of_modules,
                self.b.registry,
                "mm",
            )
            inner_lv = self.fiber_cl2_lv
        else:
            angle = 2 * np.pi / self.number_of_modules
            z, r = self._get_bend_polycone(self.radius - coating_dim / 2, self.radius + coating_dim / 2)
            coating = g4.solid.GenericPolycone(v_name, 0, angle, r, z, self.b.registry, "mm")
            inner_lv = self.fiber_cl2_bend_lv
        coating_lv = g4.LogicalVolume(coating, self.b.materials.tpb_on_fibers, v_name, self.b.registry)
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            inner_lv,
            f"fiber_{self.barrel}_cladding2{v_suffix}",
            coating_lv,
            self.b.registry,
        )

        coating_lv.pygeom_color_rgba = [0, 1, 0.165, 0.07]  # 520 nm

        return coating_lv

    def create_module(self, mod: FiberModuleData) -> None:
        assert mod.barrel == self.barrel
        module_num = _module_name_to_num(mod.name)
        if module_num < 0 or module_num >= self.number_of_modules:
            msg = f"invalid module number {module_num} for a maximum of {self.number_of_modules}-1 modules."
            raise ValueError(msg)

        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        coating_lv = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=False)
        if self.bend_radius_mm is not None:
            self._cached_sipm_volumes_bend()
            coating_lv_bend = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=True)

        start_angle = self.start_angle(module_num)
        z_displacement_straight = self.z_displacement - self.fiber_length / 2

        fibers = []

        th = start_angle
        fibers.append(
            g4.PhysicalVolume(
                [0, 0, -th],
                [0, 0, z_displacement_straight],
                coating_lv,
                f"fiber_{self.barrel}_tpb_{mod.name}_s",
                self.b.mother_lv,
                self.b.registry,
            )
        )
        if self.bend_radius_mm is not None:
            fibers.append(
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [0, 0, z_displacement_straight - self.fiber_length / 2 - self.bend_radius_mm],
                    coating_lv_bend,
                    f"fiber_{self.barrel}_tpb_{mod.name}_bend_s",
                    self.b.mother_lv,
                    self.b.registry,
                )
            )

        self._add_tpb_surfaces(fibers)

        # create SiPMs and attach to fibers
        self._create_sipm(
            module_num,
            fibers,
            True,
            mod.channel_top_name,
            mod.channel_top_rawid,
            z_displacement_straight,
        )
        if self.bend_radius_mm is None:
            self._create_sipm(
                module_num,
                fibers,
                False,
                mod.channel_bottom_name,
                mod.channel_bottom_rawid,
                z_displacement_straight,
            )
        else:
            z = z_displacement_straight - self.fiber_length / 2 - self.bend_radius_mm
            sipm_pv = g4.PhysicalVolume(
                [0, 0, -start_angle],
                [0, 0, z],
                self.sipm_lv_bend,
                mod.channel_bottom_name,
                self.b.mother_lv,
                self.b.registry,
            )
            sipm_pv.set_pygeom_active_detector(RemageDetectorInfo("optical", mod.channel_bottom_rawid))
            # Add border surface to mother volume.
            g4.BorderSurface(
                f"bsurface_lar_{mod.channel_bottom_name}",
                self.b.mother_pv,
                sipm_pv,
                self.b.materials.surfaces.to_sipm_silicon(
                    self.b.runtime_config,
                    mod.channel_bottom_name,
                ),
                self.b.registry,
            )

            g4.PhysicalVolume(
                [0, 0, -start_angle],
                [0, 0, z],
                self.sipm_outer_bottom_lv_bend,
                f"{mod.channel_bottom_name}_wrap",
                self.b.mother_lv,
                self.b.registry,
            )


def create_fiber_support_inner(b: core.InstrumentationData, z_pos: float) -> g4.LogicalVolume:
    inner_radius = 127.5  # mm, from CAD model.
    outer_radius = inner_radius + 6.5  # mm
    ring_thickness = 3  # mm
    rod_radius = 2.5  # mm

    vols = []
    tras = []

    # Create the rings
    ring = g4.solid.Tubs(
        "fiber_support_inner_ring", inner_radius, outer_radius, ring_thickness, 0, 2 * np.pi, b.registry
    )
    z_ring = (-700, -600, -300, 0, 300, 600, 700)  # mm
    for z in z_ring:
        vols.append(ring)
        tras.append([0, [0, 0, z]])

    # Create the rods
    radius_rod = (inner_radius + outer_radius) / 2
    for i in range(3):
        rl = 0
        for rings in itertools.pairwise(z_ring):
            rod_length = rings[1] - rings[0]
            rod_name = f"fiber_support_inner_rod_{rod_length}"
            if rod_name not in b.registry.solidDict:
                g4.solid.Tubs(
                    rod_name, 0, rod_radius, rod_length - ring_thickness - 2e-9, 0, 2 * np.pi, b.registry
                )
            rod = b.registry.solidDict[rod_name]

            vols.append(rod)
            phi = i * 2 * np.pi / 3
            tras.append([0, [radius_rod * np.cos(phi), radius_rod * np.sin(phi), rings[0] + rod_length / 2]])
            rl += rod_length
        assert rl == 1400

    # Combine rings and rods
    for idx, (vol, tra) in enumerate(zip(vols, tras, strict=True)):
        vol_lv = g4.LogicalVolume(
            vol,
            b.materials.metal_copper,
            f"fiber_support_inner_{idx}",
            b.registry,
        )
        vol_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)

        g4.PhysicalVolume(
            [0, 0, -tra[0]],
            np.array([0, 0, z_pos]) + np.array(tra[1]),
            vol_lv,
            f"fiber_support_inner_{idx}",
            b.mother_lv,
            b.registry,
        )


def create_fiber_support_outer(b: core.InstrumentationData, z_pos: float) -> g4.LogicalVolume:
    vols = []
    tras = []

    radius = 283 + 2  # mm. in CAD model 283 mm, enlarged to avoid fiber overlaps.
    radius_out = radius + 7
    thinring = g4.solid.Tubs("fiber_support_outer_ring", radius, radius_out, 2, 0, 2 * np.pi, b.registry)
    topring = g4.solid.Tubs("fiber_support_outer_topring", radius, radius_out, 3, 0, 2 * np.pi, b.registry)
    bottomring = g4.solid.Tubs("fiber_support_outer_bottomring", 73, 80, 2, 0, 2 * np.pi, b.registry)

    # add the 20 guiding fins.
    fin_radius = 155 + 10 + 20
    fin_x = 2
    fin_y = 8
    fin = g4.solid.Box("fiber_support_outer_fin_box", fin_x, fin_y, 1320, b.registry)
    curvedfin = g4.solid.Tubs(
        "fiber_support_outer_fin_curved",
        fin_radius - fin_y / 2,
        fin_radius + fin_y / 2,
        fin_x,
        0,
        np.pi / 2,
        b.registry,
    )

    fin = g4.solid.Union(
        "fiber_support_outer_fin",
        fin,
        curvedfin,
        [[0, np.pi / 2, 0], [0, -fin_radius, -450 - 200 - 10]],
        b.registry,
    )

    radius_fins = radius_out + fin_y / 2 + 0.01  # offset, but does not render in pyg4ometry without.
    for i in range(20):
        # Each fin needs to be rotated by 18 degrees to make the curved portion radial.
        vols.append(fin)
        tras.append(
            [
                i * 2 * np.pi / 20 - np.pi / 2,
                [radius_fins * np.cos(i * 2 * np.pi / 20), radius_fins * np.sin(i * 2 * np.pi / 20), 55 - 10],
            ]
        )

    # place the 7 rings at a spacing of 100?,300,300,300,150,150,160
    rings_z = [-600, -450, -300, 0, 300, 600, 700]
    rings_thickness = [2, 2, 2, 2, 2, 2, 3]
    for z, thickness in zip(rings_z, rings_thickness, strict=True):
        vols.append(thinring if thickness == 2 else topring)
        tras.append([0, [0, 0, z]])

    vols.append(bottomring)
    tras.append([0, [0, 0, -600 - fin_radius - 10]])

    # add the three support rods.
    radius_rod = (radius + radius_out) / 2
    for i in range(4):
        rl = 0
        for rings in itertools.pairwise(zip(rings_z, rings_thickness, strict=True)):
            rod_length = rings[1][0] - rings[0][0]
            rod_delta = max([rings[1][1], rings[0][1]])
            rod_name = f"fiber_support_outer_rod_{rod_length}"
            if rod_name not in b.registry.solidDict:
                g4.solid.Tubs(rod_name, 0, 2.5, rod_length - rod_delta - 2e-9, 0, 2 * np.pi, b.registry)
            rod = b.registry.solidDict[rod_name]

            vols.append(rod)
            phi = i * 2 * np.pi / 4
            tras.append(
                [0, [radius_rod * np.cos(phi), radius_rod * np.sin(phi), rings[0][0] + rod_length / 2]]
            )
            rl += rod_length
        assert rl == 1300

    # Combine rings and rods
    for idx, (vol, tra) in enumerate(zip(vols, tras, strict=True)):
        vol_lv = g4.LogicalVolume(
            vol,
            b.materials.metal_copper,
            f"fiber_support_outer_{idx}",
            b.registry,
        )
        vol_lv.pygeom_color_rgba = (0.72, 0.45, 0.2, 1)

        g4.PhysicalVolume(
            [0, 0, -tra[0]],
            np.array([0, 0, z_pos]) + np.array(tra[1]),
            vol_lv,
            f"fiber_support_outer_{idx}",
            b.mother_lv,
            b.registry,
        )
