from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from pyg4ometry import gdml

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
        action="store_true",
        help="""Open a VTK visualization of the generated geometry""",
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
        with Path.open(args.config) as config_file:
            config = json.load(config_file)

    registry = core.construct(
        assemblies=args.assemblies.split(","),
        use_detailed_fiber_model=args.fiber_modules == "detailed",
        config=config,
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
        from . import det_utils

        det_utils.generate_detector_macro(registry, args.det_macro_file)

    if args.vis_macro_file:
        from . import vis_utils

        vis_utils.generate_color_macro(registry, args.vis_macro_file)

    if args.visualize:
        log.info("visualizing...")
        from . import vis_utils

        vis_utils.visualize(registry)
