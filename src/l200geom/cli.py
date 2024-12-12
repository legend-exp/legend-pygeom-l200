from __future__ import annotations

import argparse
import logging

from pyg4ometry import config as meshconfig
from pyg4ometry import gdml
from pygeomtools import detectors, utils, visualization

from . import _version, core

log = logging.getLogger(__name__)


def dump_gdml_cli() -> None:
    parser = argparse.ArgumentParser(
        prog="legend-pygeom-l200",
        description="%(prog)s command line interface",
    )

    # global options
    parser.add_argument(
        "--version",
        action="version",
        help="""Print %(prog)s version and exit""",
        version=_version.__version__,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="""Increase the program verbosity""",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        help="""Increase the program verbosity to maximum""",
    )
    parser.add_argument(
        "--visualize",
        "-V",
        nargs="?",
        const=True,
        help="""Open a VTK visualization of the generated geometry (with optional scene file)""",
    )
    parser.add_argument(
        "--vis-macro-file",
        action="store",
        help="""Filename to write a Geant4 macro file containing visualization attributes""",
    )
    parser.add_argument(
        "--det-macro-file",
        action="store",
        help="""Filename to write a Geant4 macro file containing active detectors (to be used with remage)""",
    )
    parser.add_argument(
        "--check-overlaps",
        action="store_true",
        help="""Check for overlaps with pyg4ometry (note: this might not be accurate)""",
    )

    # options for geometry generation.
    geom_opts = parser.add_argument_group("geometry options")
    geom_opts.add_argument(
        "--assemblies",
        action="store",
        default=",".join(core.DEFINED_ASSEMBLIES),
        help="""Select the assemblies to generate in the output. (default: %(default)s)""",
    )
    geom_opts.add_argument(
        "--fiber-modules",
        action="store",
        choices=("segmented", "detailed"),
        default="segmented",
        help="""Select the fiber shroud model, either coarse segments or single fibers. (default: %(default)s)""",
    )
    geom_opts.add_argument(
        "--config",
        action="store",
        help="""Select a config file to read geometry config from.""",
    )
    geom_opts.add_argument(
        "--public-geom",
        action="store_true",
        help="""Create a geometry from public testdata only.""",
    )

    parser.add_argument(
        "filename",
        default="",
        nargs="?",
        help="""File name for the output GDML geometry.""",
    )

    args = parser.parse_args()

    if not args.visualize and args.filename == "":
        parser.error("no output file and no visualization specified")
    if (args.vis_macro_file or args.det_macro_file) and args.filename == "":
        parser.error("writing macro file(s) without gdml file is not possible")

    if args.verbose:
        logging.getLogger("l200geom").setLevel(logging.DEBUG)
    if args.debug:
        logging.root.setLevel(logging.DEBUG)

    config = {}
    if args.config:
        config = utils.load_dict(args.config)

    vis_scene = {}
    if isinstance(args.visualize, str):
        vis_scene = utils.load_dict(args.visualize)
        if vis_scene.get("fine_mesh", False):
            meshconfig.setGlobalMeshSliceAndStack(100)

    registry = core.construct(
        assemblies=[a for a in args.assemblies.split(",") if a != ""],
        use_detailed_fiber_model=args.fiber_modules == "detailed",
        config=config,
        public_geometry=args.public_geom,
    )

    if args.check_overlaps:
        msg = "checking for overlaps"
        log.info(msg)
        registry.worldVolume.checkOverlaps(recursive=True)

    if args.filename != "":
        msg = f"exporting GDML geometry to {args.filename}"
        log.info(msg)
        w = gdml.Writer()
        w.addDetector(registry)
        w.write(args.filename)

    if args.det_macro_file:
        detectors.generate_detector_macro(registry, args.det_macro_file)

    if args.vis_macro_file:
        visualization.generate_color_macro(registry, args.vis_macro_file)

    if args.visualize:
        log.info("visualizing...")
        from pygeomtools import viewer

        viewer.visualize(registry, vis_scene)
