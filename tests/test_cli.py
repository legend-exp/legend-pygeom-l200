from __future__ import annotations

from pathlib import Path

import pytest


def test_cli():
    from l200geom.cli import _parse_cli_args
    from l200geom.core import DEFAULT_ASSEMBLIES

    p = Path(__file__).parent.resolve() / "test_cfg"

    args, _config = _parse_cli_args([])
    assert args.fiber_modules == "detailed"
    assert args.assemblies == DEFAULT_ASSEMBLIES

    args, _config = _parse_cli_args(["--config", str(p / "cfg_geom.yaml")])
    assert args.fiber_modules == "detailed"
    assert args.assemblies == {"wlsr", "fibers"}
    args, _config = _parse_cli_args(["--config", str(p / "cfg_geom.yaml"), "--fiber-modules", "segmented"])
    assert args.fiber_modules == "segmented"
    assert args.assemblies == {"wlsr", "fibers"}
    args, _config = _parse_cli_args(
        ["--config", str(p / "cfg_geom.yaml"), "--assemblies", "strings,calibration"]
    )
    assert args.fiber_modules == "detailed"
    assert args.assemblies == {"strings", "calibration"}


def test_assemblies():
    from l200geom.cli import _parse_assemblies
    from l200geom.core import DEFAULT_ASSEMBLIES

    assert _parse_assemblies(None) == DEFAULT_ASSEMBLIES
    assert _parse_assemblies("watertank,calibration") == {"watertank", "calibration"}
    assert _parse_assemblies(["watertank", "calibration"]) == {"watertank", "calibration"}
    assert _parse_assemblies("+watertank") == DEFAULT_ASSEMBLIES | {"watertank"}
    assert _parse_assemblies("+watertank,-fibers") == (DEFAULT_ASSEMBLIES | {"watertank"}) - {"fibers"}
    assert _parse_assemblies(["+watertank", "-fibers"]) == (DEFAULT_ASSEMBLIES | {"watertank"}) - {"fibers"}

    with pytest.raises(ValueError, match="all or no assemblies can be prefixed"):
        _parse_assemblies("+watertank,fibers")
