from __future__ import annotations

import os
import re

import pygeomtools
import pytest

public_geom = os.getenv("LEGEND_METADATA", "") == ""

pytestmark = [
    pytest.mark.xfail(run=True, reason="requires a remage installation"),
    pytest.mark.needs_remage,
    pytest.mark.filterwarnings("ignore:lgdo.lh5 has moved to its own package:DeprecationWarning"),
]


@pytest.fixture
def gdml_file(tmp_path):
    from pygeoml200 import core

    registry = core.construct(config={}, public_geometry=public_geom)

    gdml_file = tmp_path / "l200-default.gdml"
    pygeomtools.write_pygeom(registry, gdml_file)

    return gdml_file


@pytest.fixture
def gdml_file_for_surface_check(tmp_path):
    from pygeoml200 import core

    registry = core.construct(
        assemblies=core.DEFINED_ASSEMBLIES - {"fibers"},
        config={
            "watertank_no_pmts": True,
        },
        public_geometry=public_geom,
    )

    gdml_file = tmp_path / "l200-surface-overlap-check.gdml"
    pygeomtools.write_pygeom(registry, gdml_file)

    return gdml_file


def _extract_stats(text):
    pattern = r"average event processing time.*?=\s*([\d.]+)\s*events/second"
    m = re.search(pattern, text, flags=re.DOTALL | re.IGNORECASE)
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

    m = re.search(pattern, text, flags=re.IGNORECASE | re.VERBOSE)
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


def test_surface_overlaps(gdml_file_for_surface_check):
    from remage import remage_run

    macro = [
        "/RMG/Output/ActivateOutputScheme GeometryCheck",
        "/run/initialize",
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/SampleOnSurface",
        "/RMG/Generator/Confinement/FirstSamplingVolume Geometrical",
        "/RMG/Generator/Confinement/Geometrical/AddSolid Box",
        "/RMG/Generator/Confinement/Geometrical/CenterPositionX 0 m",
        "/RMG/Generator/Confinement/Geometrical/CenterPositionY 0 m",
        "/RMG/Generator/Confinement/Geometrical/CenterPositionZ 0 m",
        "/RMG/Generator/Confinement/Geometrical/Box/XLength 10 m",
        "/RMG/Generator/Confinement/Geometrical/Box/YLength 10 m",
        "/RMG/Generator/Confinement/Geometrical/Box/ZLength 10 m",
        "/RMG/Generator/Select GPS",
        "/gps/particle     geantino",
        "/gps/energy       1 MeV",
        "/gps/ang/type     iso",
        "/run/beamOn       100000",
    ]

    remage_run(
        macro, gdml_files=str(gdml_file_for_surface_check), raise_on_error=True, raise_on_warning=False
    )
