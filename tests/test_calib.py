from __future__ import annotations

import os

public_geom = os.getenv("LEGEND_METADATA", "") == ""


def test_construct_with_sis(tmp_path):
    from l200geom import core

    # test single source.
    cfg = {
        "sis": {
            "1": {
                "sis_z": 8250,
                "sources": [None, None, None, "Th228"],
            },
            "2": None,
            "3": None,
            "4": None,
        },
    }

    registry = core.construct(assemblies=["calibration"], config=cfg, public_geometry=public_geom)
    assert "source_inner_sis1_source3" in registry.physicalVolumeDict

    # test single source with Cu cap.
    cfg["sis"]["1"]["sources"][3] = "Ra+Cu"
    registry = core.construct(assemblies=["calibration"], config=cfg, public_geometry=public_geom)
    assert "source_inner_sis1_source3" in registry.physicalVolumeDict

    # test multiple sources.
    cfg = {
        "sis": {
            "1": {
                "sis_z": 8250,
                "sources": ["Th228", "Th228", "Th228", "Th228"],
            },
            "2": {
                "sis_z": 8250,
                "sources": ["Ra+Cu", "Ra+Cu", "Ra+Cu", "Ra+Cu"],
            },
            "3": None,
            "4": None,
        },
    }

    registry = core.construct(assemblies=["calibration"], config=cfg, public_geometry=public_geom)
    for i in range(3):
        assert f"source_inner_sis1_source{i}" in registry.physicalVolumeDict
        assert f"source_inner_sis2_source{i}" in registry.physicalVolumeDict

    # test offsets
    cfg = {
        "sis": {
            "1": {"sis_z": 8250, "sources": ["Th228", "Th228", "Th228", "Th228"], "r_offset": 1},
            "2": {"sis_z": 8250, "sources": ["Ra+Cu", "Ra+Cu", "Ra+Cu", "Ra+Cu"], "phi_offset": 3},
            "3": None,
            "4": None,
        },
    }

    # check it works
    registry = core.construct(assemblies=["calibration"], config=cfg, public_geometry=public_geom)

    for i in range(3):
        assert f"source_inner_sis1_source{i}" in registry.physicalVolumeDict
        assert f"source_inner_sis2_source{i}" in registry.physicalVolumeDict
