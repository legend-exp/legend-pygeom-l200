from __future__ import annotations

import logging
import warnings
from collections.abc import Container, Sequence
from importlib import resources

import numpy as np
import pyg4ometry
from pyg4ometry import geant4
from scipy.spatial.transform import Rotation

from . import core

log = logging.getLogger(__name__)

COLORS = {
    "steel": (0.5, 0.5, 0.5, 0.05),
    "water": (0, 0, 1, 0.08),
    "air": (0.1, 0.1, 0.1, 0.025),
    "acrylic": (0.9, 0.9, 0.9, 0.05),
    "vm2000": (0.9, 0.9, 0.9, 0.05),
    "tetratex": (1, 1, 1, 0.1),
    # 520 nm. The l200 shroud is far denser than the l1000 curtain, so it needs a lower
    # alpha than l1000 uses to keep the array behind it visible.
    "fiber_coating": (0, 1, 0.165, 0.07),
    "pmt_window": (0.9, 0.8, 0.5, 0.05),
    "pmt_cathode": (0.545, 0.271, 0.074, 0.05),
}


def _read_model(
    file: str, name: str, material: geant4.Material, b: core.InstrumentationData
) -> geant4.LogicalVolume | None:
    """
    Construct a logical volume for an STL mesh.

    .. note::
        This function honours the ``no_meshes`` runtime configuration, which can either be ``True``
        to disable all meshes, or be a list of logical volume names to disable mesh loading.

    Returns
    -------
    A :class:`geant4.LogicalVolume` for the mesh or ``None``, if loading of this mesh is disabled.
    """
    # this is an (undocumented) option to remove meshes; either all or from a list (for performance tests).
    no_meshes = b.runtime_config.get("no_meshes", False)
    if (isinstance(no_meshes, Container) and (name in no_meshes)) or no_meshes is True:
        log.warning("skipping mesh %s", name)
        return None

    res = resources.files("pygeoml200") / "models" / file
    solid = pyg4ometry.stl.Reader(res, solidname=name, centre=False, registry=b.registry).getSolid()
    return geant4.LogicalVolume(solid, material, name, b.registry)


def place_disjoint_union(
    name: str,
    objects: Sequence[geant4.LogicalVolume],
    transformations: Sequence,
    mother_lv: geant4.LogicalVolume,
    registry: geant4.Registry,
    rotation: Sequence[float] = (0, 0, 0),
    position: Sequence[float] = (0, 0, 0),
) -> list[geant4.PhysicalVolume]:
    """Place a "disjoint union" of logical volumes as separate physical volumes.

    This is a drop-in replacement for constructing and placing a :class:`geant4.solid.MultiUnion`
    of disjoint solids.

    .. note::
        It is not checked whether the volumes are actually disjoint.

    Parameters
    ----------
    name
        name prefix of the physical volumes. The volumes will be named ``{name}_{idx}``.
    objects
        logical volumes to place.
    transformations
        ``[[rot1, tra1], [rot2, tra2], ...]``, with the same convention as for the nodes of a
        :class:`geant4.solid.MultiUnion`, i.e. relative to the frame of the union.
    mother_lv
        logical volume to place the union into.
    registry
        the registry to add the physical volumes to.
    rotation
        rotation of the whole union in the mother volume.
    position
        position of the whole union in the mother volume.

    Returns
    -------
    The list of created physical volumes, in the order of ``objects``.
    """
    if len(objects) != len(transformations):
        msg = "objects and transformations must have the same length"
        raise ValueError(msg)

    def _eval_vec(v):
        return np.array([float(x) for x in v])

    # Note: Geant4 uses inverse conventions for the rotations of physical volumes and of MultiUnion nodes:
    # GDML physvol angles (a, b, c) correspond to the active rotation R(a, b, c)^-1, while the same angles
    # for a MultiUnion node correspond to the active rotation R(a, b, c), with R(a, b, c) = Rz(c) Ry(b) Rx(a).
    top_rot = Rotation.from_euler("xyz", _eval_vec(rotation))
    top_pos = _eval_vec(position)

    pvs = []
    for idx, (lv, (node_rotvec, node_pos)) in enumerate(zip(objects, transformations, strict=True)):
        node_rot = Rotation.from_euler("xyz", _eval_vec(node_rotvec))
        # the active transform of the node in the mother is x -> top_rot^-1 (node_rot x + node_pos) + top_pos.
        with warnings.catch_warnings(action="ignore"):  # ignore gimbal lock warnings.
            pv_rot = (node_rot.inv() * top_rot).as_euler("xyz")
        pv_pos = top_pos + top_rot.inv().apply(_eval_vec(node_pos))

        pvs.append(
            geant4.PhysicalVolume(list(pv_rot), list(pv_pos), lv, f"{name}_{idx}", mother_lv, registry)
        )

    return pvs
