from __future__ import annotations

import os
import re

import pygeomtools
import pytest

public_geom = os.getenv("LEGEND_METADATA", "") == ""

pytestmark = [
    pytest.mark.xfail(run=True, reason="requires a remage installation"),
    pytest.mark.needs_remage,
]


@pytest.fixture
def gdml_file(tmp_path):
    from pygeoml200 import core

    registry = core.construct(config={}, public_geometry=public_geom)

    gdml_file = tmp_path / "l200-default.gdml"
    pygeomtools.write_pygeom(registry, gdml_file)

    return gdml_file


def _extract_stats(text):
    pattern = r"average event processing time.*?=\s*([\d.]+)\s*events/second"
    m = re.search(pattern, text, flags=re.S | re.I)
    assert m is not None
    event_rate = float(m.group(1))

    pattern = r"""
    run\ time\ was
    \s*
    (\d+)\ days?,\s*
    (\d+)\ hours?,\s*
    (\d+)\ minutes?\s*and\s*
    (\d+)\ seconds?
    """

    m = re.search(pattern, text, flags=re.I | re.X)
    assert m is not None

    days, hours, minutes, seconds = map(float, m.groups())
    runtime_seconds = days * 24 * 3600 + hours * 3600 + minutes * 60 + seconds

    print(f"runtime was: {runtime_seconds} s")
    print(f"event rate was: {event_rate} event/s")

    return runtime_seconds, event_rate


def _benchmark(macro, gdml_file, capfd):
    from remage import remage_run

    remage_run([m.strip() for m in macro.split("\n")], gdml_files=str(gdml_file))
    # remage sends to stderr
    stderr = capfd.readouterr().err

    return _extract_stats(stderr)


def test_performance(gdml_file, capfd):
    macro = """
    /RMG/Geometry/GDMLDisableOverlapCheck

    /run/initialize

    /RMG/Generator/Select GPS
    /gps/particle geantino
    /gps/ang/type iso

    /run/beamOn 100000
    """

    runtime, event_rate = _benchmark(macro, gdml_file, capfd)

    assert runtime > 1
    assert event_rate > 1_000

    macro = """
    /RMG/Geometry/GDMLDisableOverlapCheck

    /run/initialize

    /RMG/Generator/Select GPS
    /gps/particle geantino
    /gps/ang/type iso

    /RMG/Generator/Confine Volume
    /RMG/Generator/Confinement/Physical/AddVolume V.*

    /run/beamOn 50000
    """

    runtime, event_rate = _benchmark(macro, gdml_file, capfd)

    assert runtime > 1
    assert event_rate > 1_000


def test_overlaps(gdml_file):
    from remage import remage_run

    macro = [
        "/RMG/Geometry/RegisterDetectorsFromGDML Germanium",
        "/RMG/Geometry/RegisterDetectorsFromGDML Scintillator",
        "/RMG/Geometry/RegisterDetectorsFromGDML Optical",
        "/run/initialize",
    ]

    remage_run(macro, gdml_files=str(gdml_file), raise_on_error=True, raise_on_warning=True)
