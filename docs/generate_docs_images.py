#!/bin/env python3
from __future__ import annotations

import logging

from pyg4ometry import config as meshconfig
from pygeomtools import viewer, write_pygeom

from pygeoml200 import core

logging.basicConfig()
meshconfig.setGlobalMeshSliceAndStack(100)

images = {
    "fibers": {"assemblies": ["fibers"]},
    "top_plate": {"assemblies": ["top"]},
    "holders": {
        "assemblies": ["strings"],
        "overrides": {
            "minishroud_.*": False,
            "[BVPC].*": False,
            "pen_.*": [0, 0, 1, 1],
        },
    },
    "hpge_strings": {
        "assemblies": ["strings"],
        "overrides": {
            "pen_.*": False,
            "minishroud_.*": False,
            "hpge_support_copper_.*": False,
            "cable.*": False,
            "ultem_.*": False,
        },
    },
    "nylon": {
        "assemblies": ["strings", "calibration"],
        "overrides": {
            "pen_.*": False,
            "hpge_support_copper_.*": False,
            "cable.*": False,
            "ultem_.*": False,
            "[BVPC].*": False,
        },
    },
    "wlsr": {
        "assemblies": ["strings", "calibration", "fibers", "wlsr", "top"],
        "overrides": {"pen_.*": [0, 0, 1, 1]},
        "default": {
            # "focus": [0, 0, 0],
            # "up": [0.45, 0, 0.89],
            # "camera": [-6885.44, 64.16, 3470.46],
            "focus": [131, 0, 259],
            "up": [0.45, 0, 0.89],
            "camera": [-4572.28, 0, 2629.65],
        },
        "window_size": [571, 1000],
    },
}


def export_image(fn: str, extra: dict) -> None:
    vis_default = extra.get(
        "default",
        {
            "focus": [292.20, 0, 574.37],
            "up": [0.45, 0, 0.89],
            "camera": [-2910.76, -65.44, 2208.86],
        },
    )
    vis_scene = {
        "window_size": extra.get("window_size", [400, 700]),
        "default": vis_default,
        "color_overrides": {"lar": False, **extra.get("overrides", {})},
        "export_scale": 1,
        "export_and_exit": f"source/images/{fn}.png",
    }

    registry = core.construct(
        assemblies=extra["assemblies"],
        use_detailed_fiber_model=True,
        public_geometry=True,
    )
    write_pygeom(registry, None)
    viewer.visualize(registry, vis_scene)


for fn, extra in images.items():
    export_image(fn, extra)
