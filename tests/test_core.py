from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from pyg4ometry import gdml
from pygeomtools import detectors, write_pygeom

public_geom = os.getenv("LEGEND_METADATA", "") == ""


def test_import():
    import l200geom  # noqa: F401


@pytest.fixture(scope="session")
def conctruct_fiber_variants():
    from l200geom import core

    assemblies = core.DEFINED_ASSEMBLIES
    reg_detailed = core.construct(assemblies, use_detailed_fiber_model=True, public_geometry=public_geom)
    reg_segmented = core.construct(assemblies, use_detailed_fiber_model=False, public_geometry=public_geom)

    return reg_detailed, reg_segmented


def test_construct(tmp_path, conctruct_fiber_variants):
    # do nothing here, we just want to test the fixture.
    _reg_detailed, _reg_segmented = conctruct_fiber_variants


@pytest.mark.xfail(
    public_geom,
    reason="public geometry has different detector count",
    raises=AssertionError,
)
def test_detector_count(tmp_path, conctruct_fiber_variants):
    reg_detailed, reg_segmented = conctruct_fiber_variants

    # verify that we get the expected channel counts.
    ch_count = Counter([d.detector_type for f, d in detectors.walk_detectors(reg_detailed)])
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)
    assert ch_count["germanium"] == 101  # from channelmap @ 20230311T235840Z
    det_file_detailed = tmp_path / "det-detailed.mac"
    detectors.generate_detector_macro(reg_detailed, det_file_detailed)

    # verify that we get the expected channel counts.
    ch_count = Counter([d.detector_type for f, d in detectors.walk_detectors(reg_segmented)])
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)
    assert ch_count["germanium"] == 101  # from channelmap @ 20230311T235840Z
    det_file_segmented = tmp_path / "det-segmented.mac"
    detectors.generate_detector_macro(reg_segmented, det_file_segmented)

    with (
        Path(det_file_detailed).open(encoding="utf-8") as f_det,
        Path(det_file_segmented).open(encoding="utf-8") as f_seg,
    ):
        assert f_det.readlines() == f_seg.readlines()


def test_read_back(tmp_path, conctruct_fiber_variants):
    _, reg_segmented = conctruct_fiber_variants

    # write a GDML file.
    gdml_file_detailed = tmp_path / "segmented.gdml"
    write_pygeom(reg_segmented, gdml_file_detailed)

    # try to read it back and check detector info.
    reader = gdml.Reader(gdml_file_detailed)
    reg = reader.getRegistry()
    ch_count = Counter([d.detector_type for d in detectors.get_all_sensvols(reg).values()])
    assert (
        ch_count["germanium"] > 90
    )  # the number of germanium detectors is not constant (public/private geometry).
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)


def test_material_store():
    # replacing material properties is _not_ a core functionality of this package, but
    # we have to make sure that replaced material properties from the optics package are
    # propagated correctly to the generated GDML files.

    from legendoptics import store
    from legendoptics.fibers import fiber_core_refractive_index

    from l200geom import core

    # test that replaced material properties are reflected in the GDML.
    fiber_core_refractive_index.replace_implementation(lambda: 1234)
    reg = core.construct(
        use_detailed_fiber_model=False,
        assemblies=["fibers"],
        public_geometry=public_geom,
    )
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1234, 1234])

    # test that after the reset, the created GDML contains the original values again.
    store.reset_all_to_original()
    reg = core.construct(
        use_detailed_fiber_model=False,
        assemblies=["fibers"],
        public_geometry=public_geom,
    )
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1.6, 1.6])


def test_material_store_cli(change_dir, tmp_path):
    # replacing material properties is _not_ a core functionality of this package, but
    # we have to make sure that replaced material properties from the optics package are
    # propagated correctly to the generated GDML files.

    from legendoptics import store

    from l200geom import cli

    output_file = tmp_path / "matstore.gdml"

    assert not output_file.exists()
    cli.dump_gdml_cli(
        [
            "--pygeom-optics-plugin",
            "test_cfg/matprop_change.py",
            "--assemblies=wlsr",
            str(output_file),
            "--public-geom",
        ]
    )
    assert output_file.exists()

    # test that replaced material properties are reflected in the GDML.
    reader = gdml.Reader(output_file)
    reg = reader.getRegistry()
    rindex = reg.defineDict["tpb_on_tetratex_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1234, 1234])

    # this is not a real CLI invocation, so we need to reset the store for the next tests.
    store.reset_all_to_original()


@pytest.fixture
def change_dir(request):
    os.chdir(request.fspath.dirname)
    yield
    os.chdir(request.config.invocation_params.dir)


def test_special(change_dir, tmp_path):
    from l200geom import cli

    output_file = tmp_path / "special.gdml"

    assert not output_file.exists()
    cli.dump_gdml_cli(["--config", "test_cfg/cfg_central.yaml", str(output_file), "--public-geom"])
    assert output_file.exists()

    # try to read it back and check detector info.
    reader = gdml.Reader(output_file)
    reg = reader.getRegistry()
    ch_count = Counter([d.detector_type for d in detectors.get_all_sensvols(reg).values()])
    assert ch_count["germanium"] == 0  # no germanium in channelmap.
    assert ch_count["optical"] == 2 * (9 + 20)  # 2*(IB+OB)

    # we should have only short hangers and wrapped counterweights.
    assert "counterweight_wrapped" in reg.solidDict
    assert "hpge_support_copper_string_support_structure_short" in reg.logicalVolumeDict
    assert "hpge_support_copper_string_support_structure" not in reg.logicalVolumeDict
