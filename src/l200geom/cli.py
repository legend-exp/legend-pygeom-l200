from __future__ import annotations

import argparse
import logging

from pyg4ometry import gdml

from . import _version
from .core import construct

log = logging.getLogger(__name__)


def dump_gdml_cli() -> None:
    parser = argparse.ArgumentParser(
        prog="legend-pygeom-l200",
        description="legend-pygeom-l200 command line interface",
    )

    # global options
    parser.add_argument(
        "--version",
        action="version",
        help="""Print %(prog)s version and exit""",
        version=_version.__version__,  # noqa: T201
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
        "filename",
        default="l200.gdml",
        help="""File name for the output GDML geometry.""",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger("l200geom").setLevel(logging.DEBUG)
    if args.debug:
        logging.root.setLevel(logging.DEBUG)

    log.info(f"exporting GDML geometry to {args.filename}")
    w = gdml.Writer()
    registry = construct()
    w.addDetector(registry)
    w.write(args.filename)

    if args.visualize:
        log.info("visualizing...")
        from pyg4ometry import visualisation

        v = visualisation.VtkViewer()
        v.addLogicalVolume(registry.worldVolume)
        v.addAxes(length=5000)
        v.view()
