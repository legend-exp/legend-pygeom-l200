from __future__ import annotations

import argparse
import logging
from collections.abc import Iterable

from dbetto import utils
from legendoptics.store import load_user_material_code
from pyg4ometry import config as meshconfig
from pygeomtools import detectors, geometry, visualization, write_pygeom

from . import _version, core

log = logging.getLogger(__name__)


def dump_gdml_cli(argv: list[str] | None = None) -> None:
    args, config = _parse_cli_args(argv)

    logging.basicConfig()
    if args.verbose:
        logging.getLogger("l200geom").setLevel(logging.DEBUG)
    if args.debug:
        logging.root.setLevel(logging.DEBUG)

    vis_scene = {}
    if isinstance(args.visualize, str):
        vis_scene = utils.load_dict(args.visualize)

    if vis_scene.get("fine_mesh", False) or args.check_overlaps:
        meshconfig.setGlobalMeshSliceAndStack(100)

    # load custom module to change material properties.
    if args.pygeom_optics_plugin:
        load_user_material_code(args.pygeom_optics_plugin)

    registry = core.construct(
        assemblies=args.assemblies,
        pmt_configuration_mv=args.pmt_config,
        use_detailed_fiber_model=args.fiber_modules == "detailed",
        config=config,
        public_geometry=args.public_geom,
    )

    if args.print_volumes:
        geometry.print_volumes(registry, args.print_volumes)

    if args.check_overlaps:
        msg = "checking for overlaps"
        log.info(msg)
        registry.worldVolume.checkOverlaps(recursive=True)

    # commit auxvals, and write to GDML file if requested.
    if args.filename is not None:
        log.info("exporting GDML geometry to %s", args.filename)
    write_pygeom(registry, args.filename)

    if args.det_macro_file:
        detectors.generate_detector_macro(registry, args.det_macro_file)

    if args.vis_macro_file:
        visualization.generate_color_macro(registry, args.vis_macro_file)

    if args.visualize:
        log.info("visualizing...")
        from pygeomtools import viewer

        viewer.visualize(registry, vis_scene)


def _parse_cli_args(argv: list[str] | None = None) -> tuple[argparse.Namespace, dict]:
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
    parser.add_argument(
        "--print-volumes",
        action="store",
        choices=("logical", "physical", "detector"),
        help="""Print a list of logical or physical volume names (from the pyg4ometry registry)""",
    )

    parser.add_argument(
        "--pygeom-optics-plugin",
        action="store",
        help="""Execute the python module given by this path before constructing the geometry""",
    )

    # options for geometry generation.
    #
    # geometry options can also be specified in the config file, so the "default" argument of the argparse
    # options cannot be used - we need to distinguish between an unspecified option and an explicitly set
    # default option.
    geom_opts = parser.add_argument_group("geometry options")
    extra_assemblies = core.DEFINED_ASSEMBLIES - core.DEFAULT_ASSEMBLIES
    geom_opts.add_argument(
        "--assemblies",
        action="store",
        help=(
            f"""Select the assemblies to generate in the output.
            (default: {",".join(core.DEFAULT_ASSEMBLIES)};
            additionally available: {",".join(extra_assemblies)})"""
        ),
    )
    fiber_modules_default = "detailed"
    geom_opts.add_argument(
        "--fiber-modules",
        action="store",
        choices=("segmented", "detailed"),
        help=f"""Select the fiber shroud model, either coarse segments or single fibers. (default: {fiber_modules_default})""",
    )
    pmt_config_default = "LEGEND200"
    geom_opts.add_argument(
        "--pmt-config",
        action="store",
        choices=core.PMT_CONFIGURATIONS,
        help=f"""Select the PMT configuration of the muon veto. (default: {pmt_config_default})""",
    )
    geom_opts.add_argument(
        "--public-geom",
        action="store_true",
        default=None,
        help="""Create a geometry from public testdata only.""",
    )
    geom_opts.add_argument(
        "--config",
        action="store",
        help="""Select a config file to read geometry config from.""",
    )

    parser.add_argument(
        "filename",
        default=None,
        nargs="?",
        help="""File name for the output GDML geometry.""",
    )

    args = parser.parse_args(argv)

    config = {}
    if args.config is not None:
        config = utils.load_dict(args.config)

    # also load geometry options from config file.
    _config_or_cli_arg(args, config, "assemblies", None)
    _config_or_cli_arg(args, config, "fiber_modules", fiber_modules_default)
    _config_or_cli_arg(args, config, "pmt_config", pmt_config_default)
    _config_or_cli_arg(args, config, "public_geom", False)

    # process assembly list after loading all parameter sources.
    args.assemblies = _parse_assemblies(args.assemblies)

    if args.fiber_modules not in ("segmented", "detailed"):
        msg = f"invalid fiber module type {args.fiber_modules}"
        raise ValueError(msg)

    if not args.visualize and args.filename == "":
        parser.error("no output file and no visualization specified")
    if (args.vis_macro_file or args.det_macro_file) and args.filename == "":
        parser.error("writing macro file(s) without gdml file is not possible")

    return args, config


def _parse_assemblies(arg: str | Iterable[str] | None) -> set[str]:
    """Parse an argument string into a set of assemblies to build the geometry for.

    Parameters
    ----------
    arg
        if it is a string, it will be split on commas. Parts can be prefixed by "set operators" '+' or '-',
        to add or remove items from the list of default assemblies. But either all parts have to be prefixed
        with such an operator, or none.
    """
    if arg is None or len(arg) == 0:
        return core.DEFAULT_ASSEMBLIES

    parts = [a.strip() for a in arg.split(",") if a.strip() != ""] if isinstance(arg, str) else arg

    with_no_op = [a[0] not in ("+", "-") for a in parts]
    if any(with_no_op) and not all(with_no_op):
        msg = "either all or no assemblies can be prefixed by the operators '+' or '-'"
        raise ValueError(msg)

    if not any(with_no_op):  # all have operators
        assemblies = set(core.DEFAULT_ASSEMBLIES)  # make a copy
        for p in parts:
            if p[0] == "-":
                assemblies -= {p[1:]}
            elif p[0] == "+":
                assemblies |= {p[1:]}
    else:
        assemblies = set(parts)

    return assemblies


def _config_or_cli_arg(args: argparse.Namespace, config: dict, name: str, default) -> None:
    """Fallback of cli args, to config file, and to default value (in this order)."""
    val_cfg = config.get(name)
    val_attrs = getattr(args, name, None)
    val = val_cfg if val_attrs is None else val_attrs
    val = default if val is None else val
    setattr(args, name, val)
