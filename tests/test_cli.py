from __future__ import annotations

from pathlib import Path


def test_cli():
    from l200geom import core
    from l200geom.cli import _parse_cli_args

    p = Path(__file__).parent.resolve() / "test_cfg"

    args, config = _parse_cli_args([])
    assert args.fiber_modules == "segmented"
    assert args.assemblies == core.DEFAULT_ASSEMBLIES

    args, config = _parse_cli_args(["--config", str(p / "cfg_geom.yaml")])
    assert args.fiber_modules == "detailed"
    assert args.assemblies == ["wlsr", "fibers"]
    args, config = _parse_cli_args(["--config", str(p / "cfg_geom.yaml"), "--fiber-modules", "segmented"])
    assert args.fiber_modules == "segmented"
    assert args.assemblies == ["wlsr", "fibers"]
    args, config = _parse_cli_args(
        ["--config", str(p / "cfg_geom.yaml"), "--assemblies", "strings,calibration"]
    )
    assert args.fiber_modules == "detailed"
    assert args.assemblies == ["strings", "calibration"]
