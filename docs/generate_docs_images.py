#!/usr/bin/env python3
from __future__ import annotations

import logging
import sys

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
        },
    },
    "hpge_strings": {
        "assemblies": ["strings"],
        "overrides": {
            "pen_.*": False,
            "minishroud_.*": False,
            "hpge_(string_)?support_(.*_)?copper.*": False,
            "hpge_string_.*_board_copper.*": False,
            "hpge_du_pin_.*": False,
            "cable.*": False,
            "(^|.*_)ultem_.*": False,
            "(^|.*_)phbr_.*": False,
            "lmfe": False,
        },
    },
    "nylon": {
        "assemblies": ["strings", "calibration"],
        "overrides": {
            "pen_.*": False,
            "hpge_(string_)?support_(.*_)?copper.*": False,
            "hpge_string_.*_board_copper.*": False,
            "hpge_du_pin_.*": False,
            "cable.*": False,
            "(^|.*_)ultem_.*": False,
            "(^|.*_)phbr_.*": False,
            "lmfe": False,
            "[BVPC].*": False,
        },
    },
    "wlsr": {
        "assemblies": ["strings", "calibration", "fibers", "wlsr", "top"],
        "overrides": {".*_pen_.*": [0, 0, 1, 1]},
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
        # none of these renderings show the cryostat, which now carries a (faint) color of its own
        "color_overrides": {
            "liquid_argon": False,
            "cryostat_outer_wall": False,
            "cryostat_inner_wall": False,
            **extra.get("overrides", {}),
        },
        "export_scale": 2,
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
    if len(sys.argv) > 1 and fn not in sys.argv[1:]:
        continue
    export_image(fn, extra)
