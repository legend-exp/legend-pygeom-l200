# Welcome to l200geom's documentation!

```{warning}
This is an early version of the LEGEND-200 geometry implemented with the
python-based simulation stack. It is not a drop-in replacement for MaGe, and
still under heavy development!
```

## Installation and usage

This package requires a working setup of
[`legend-metadata`](https://github.com/legend-exp/legend-metadata) before usage.

Following a git checkout, the package and its other python dependencies can be
installed with:

```
pip install -e .
```

If you do not intend to edit the python code in this geometry package, you can
omit the `-e` option.

After installation, the CLI utility `legend-pygeom-l200` is provided on your
PATH. This CLI utility is the primary way to interact with this package. For
now, you can find usage docs by running `legend-pygeom-l200 -h`.

Some geometry options can both be set on the CLI utility and on the config file.
Those are described in [the docs on geometry configuration](cfg-geometry), but
the descriptions similarly apply to the CLI options.

```{toctree}
:maxdepth: 2

runtime-cfg
vis
coordinate_systems
```

```{toctree}
:maxdepth: 1
:caption: Development

geom-dev
naming
Package API reference <api/modules>
```
