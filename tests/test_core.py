from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from pyg4ometry import gdml
from pygeomtools import detectors

public_geom = os.getenv("LEGEND_METADATA", "") == ""


def test_import():
    import l200geom  # noqa: F401


@pytest.mark.xfail(
    public_geom,
    reason="public geometry has different detector count",
    raises=AssertionError,
)
def test_construct(tmp_path):
    from l200geom import core

    registry = core.construct(use_detailed_fiber_model=True, public_geometry=public_geom)
    # verify that we get the expected channel counts.
    ch_count = Counter([d.detector_type for f, d in detectors.walk_detectors(registry)])
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)
    assert ch_count["germanium"] == 101  # from channelmap @ 20230311T235840Z
    det_file_detailed = tmp_path / "det-detailed.mac"
    detectors.generate_detector_macro(registry, det_file_detailed, public_geometry=public_geom)

    registry = core.construct(use_detailed_fiber_model=False)
    # verify that we get the expected channel counts.
    ch_count = Counter([d.detector_type for f, d in detectors.walk_detectors(registry)])
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)
    assert ch_count["germanium"] == 101  # from channelmap @ 20230311T235840Z
    det_file_segmented = tmp_path / "det-segmented.mac"
    detectors.generate_detector_macro(registry, det_file_segmented)

    with (
        Path(det_file_detailed).open(encoding="utf-8") as f_det,
        Path(det_file_segmented).open(encoding="utf-8") as f_seg,
    ):
        assert f_det.readlines() == f_seg.readlines()


def test_read_back(tmp_path):
    from l200geom import core

    registry = core.construct(use_detailed_fiber_model=False, public_geometry=public_geom)
    # write a GDML file.
    gdml_file_detailed = tmp_path / "segmented.gdml"
    w = gdml.Writer()
    w.addDetector(registry)
    w.write(gdml_file_detailed)
    # try to read it back.
    gdml.Reader(gdml_file_detailed)


def test_material_store():
    # replacing material properties is _not_ a core functionality of this package, but
    # we have to make sure that replaced material properties from the optics package are
    # propagated correctly to the generated GDML files.

    from legendoptics import store
    from legendoptics.fibers import fiber_core_refractive_index

    from l200geom import core

    # test that replaced material properties are reflected in the GDML.
    fiber_core_refractive_index.replace_implementation(lambda: 1234)
    reg = core.construct(use_detailed_fiber_model=False, public_geometry=public_geom)
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1234, 1234])

    # test that after the reset, the created GDML contains the original values again.
    store.reset_all_to_original()
    reg = core.construct(use_detailed_fiber_model=False, public_geometry=public_geom)
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1.6, 1.6])
