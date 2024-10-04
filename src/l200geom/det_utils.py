from __future__ import annotations

import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pyg4ometry.geant4 as g4

log = logging.getLogger(__name__)


@dataclass
class RemageDetectorInfo:
    detector_type: Literal["optical", "germanium"]
    uid: int


def walk_detectors(
    pv: g4.PhysicalVolume | g4.LogicalVolume | g4.Registry,
) -> Generator[tuple[g4.PhysicalVolume, RemageDetectorInfo], None, None]:
    if isinstance(pv, g4.PhysicalVolume) and hasattr(pv, "pygeom_active_dector"):
        det = pv.pygeom_active_dector
        assert isinstance(det, RemageDetectorInfo)
        yield pv, det

    if isinstance(pv, g4.LogicalVolume):
        next_v = pv
    if isinstance(pv, g4.PhysicalVolume):
        next_v = pv.logicalVolume
    elif isinstance(pv, g4.Registry):
        next_v = pv.worldVolume
    else:
        msg = "invalid type encountered in walk_detectors volume tree"
        raise TypeError(msg)

    for dv in next_v.daughterVolumes:
        if dv.type == "placement":
            yield from walk_detectors(dv)


def generate_detector_macro(registry: g4.Registry, filename: str) -> None:
    """Create a Geant4 macro file containing the defined active detector volumes for use in remage."""

    macro_lines = {}

    for pv, det in walk_detectors(registry):
        if pv.name in macro_lines:
            return
        mac = f"/RMG/Geometry/RegisterDetector {det.detector_type.title()} {pv.name} {det.uid}\n"
        macro_lines[pv.name] = mac

    macro_contents = "".join(macro_lines.values())

    with Path.open(filename, "w", encoding="utf-8") as f:
        f.write(macro_contents)
