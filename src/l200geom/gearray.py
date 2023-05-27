from __future__ import annotations

import json
import math

from legendhpges import PPC, BEGe, InvertedCoax, SemiCoax
from legendmeta.jsondb import AttrsDict, JsonDB
from pyg4ometry import geant4


def place_gearray(
    metadata: str | JsonDB,
    array_config: str | dict | AttrsDict,
    ch_map: str | dict | AttrsDict,
    z0: float,
    mothervolume: geant4.LogicalVolume,
    registry: geant4.Registry,
) -> None:
    """LEGEND-200 germanium detector array constructor.

    Parameters
    ----------
    metadata
        Path to the LEGEND-200 HPGe configuration metadata file describing the
        detector shape.
    array_config
        LEGEND-200 germanium detector array configuration file. Used to reconstruct the spatial position of each array.
    ch_map
        LEGEND-200 germanium detector channel map containing the position of each detector in each array.
    z0
        The position of the top array slot in the z-axis.
    motherVolume
        pyg4ometry Geant4 LogicalVolume instance (or its subclass) in which the arrays are to be placed.
    registry
        pyg4ometry Geant4 registry instance.
    """
    if registry is None:
        raise ValueError("registry cannot be None")

    if metadata is None:
        raise ValueError("metadata cannot be None")

    if array_config is None:
        raise ValueError("array configuration cannot be None")

    if ch_map is None:
        raise ValueError("channel map cannot be None")

    if not isinstance(metadata, JsonDB):
        metadata = JsonDB(metadata)

    if not isinstance(array_config, (dict, AttrsDict)):
        with open(array_config) as jfile:
            a_config = AttrsDict(json.load(jfile))
    else:
        a_config = AttrsDict(array_config)

    if not isinstance(ch_map, (dict, AttrsDict)):
        with open(ch_map) as jfile:
            cmap = AttrsDict(json.load(jfile))
    else:
        cmap = AttrsDict(ch_map)

    def _sin(a):
        return math.sin(math.pi * a / 180)

    def _cos(a):
        return math.cos(math.pi * a / 180)

    for name in metadata.keys():
        # Check if the detector is contained in the channel map.
        if name not in cmap.keys():
            continue

        array_idx = str(cmap[name].location.string)
        ver_pos = cmap[name].location.position

        x_pos = a_config.array[array_idx].radius_in_mm * _cos(
            a_config.array[array_idx].angle_in_deg
        )

        y_pos = a_config.array[array_idx].radius_in_mm * _sin(
            a_config.array[array_idx].angle_in_deg
        )

        z_pos = (
            z0
            - sum(a_config.array[array_idx].heights_in_mm[: ver_pos - 1])
            - a_config.array[array_idx].heights_in_mm[ver_pos - 1] / 2
        )

        if name[0] == "P":
            detector = PPC(metadata[name], registry=registry)

        elif name[0] == "B":
            detector = BEGe(metadata[name], registry=registry)

        elif name[0] == "V":
            detector = InvertedCoax(metadata[name], registry=registry)

        elif name[0] == "C":
            detector = SemiCoax(metadata[name], registry=registry)

        geant4.PhysicalVolume(
            [
                0,
                0,
                0,
            ],
            [x_pos, y_pos, z_pos],
            detector,
            name,
            mothervolume,
            registry,
        )
