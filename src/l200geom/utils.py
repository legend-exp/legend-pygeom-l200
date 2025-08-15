from __future__ import annotations

import logging
from collections.abc import Container
from importlib import resources

import pyg4ometry
from pyg4ometry import geant4

from . import core

log = logging.getLogger(__name__)


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

    file = resources.files("l200geom") / "models" / file
    solid = pyg4ometry.stl.Reader(file, solidname=name, centre=False, registry=b.registry).getSolid()
    return geant4.LogicalVolume(solid, material, name, b.registry)
