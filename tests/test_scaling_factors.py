from __future__ import annotations

import os

import numpy as np
from dbetto import AttrsDict

from l200geom import core

public_geom = os.getenv("LEGEND_METADATA", "") == ""


def test_construct_with_efficiencies():
    cfg = AttrsDict(
        {
            "sipm_efficiencies": AttrsDict({"S002": 0.5}),
            "sipm_use_pde_curve": False,
        }
    )

    registry = core.construct(assemblies=["fibers"], config=cfg, public_geometry=public_geom)
    assert "surface_to_sipm_silicon_S002_EFFICIENCY" in registry.defineDict

    # Get the optical surface object
    surf_obj = registry.defineDict["surface_to_sipm_silicon_S002_EFFICIENCY"]

    # Read the efficiency property
    eff_values = surf_obj.eval()[0, ::][1::2]

    # Assert that all values are 0.5
    assert np.allclose(eff_values, 0.5)
