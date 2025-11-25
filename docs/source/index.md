# Welcome to pygeoml200's documentation!

<img src="_images/wlsr.png" alt="" align="right" style="height: 350px; padding: 2em">

Python package containing the Monte Carlo geometry implementation of the
[LEGEND-200 experiment](https://legend-exp.org).

This geometry can be used as an input to the
[remage](https://remage.readthedocs.io/en/stable/) simulation software.

This package is based on {doc}`pyg4ometry <pyg4ometry:index>`,
{doc}`legend-pygeom-hpges <pygeomhpges:index>` (implementation of HPGe
detectors), {doc}`legend-pygeom-optics <pygeomoptics:index>` (optical properties
of materials) and {doc}`legend-pygeom-tools <pygeomtools:index>`.

:::{warning}

This is an early version of the LEGEND-200 geometry implemented with the
python-based simulation stack. It is not a drop-in replacement for MaGe, and
still under heavy development!

:::

## Installation

:::{important}

For using all its features, this package requires a working setup of
[`legend-metadata`](https://github.com/legend-exp/legend-metadata) (_private
repository_) before usage. A limited public geometry is also implemented.

:::

The latest tagged version and all its dependencies can be installed from PyPI:
`pip install legend-pygeom-l200`.

Alternatively, the packages's development version can be installed from a git
checkout: `pip install -e .` (in the directory of the git checkout).

## Usage as CLI tool

After installation, the CLI utility `legend-pygeom-l200` is provided on your
`$PATH`. This CLI utility is the primary way to interact with this package. For
now, you can find usage docs by running `legend-pygeom-l200 -h`.

In the simplest case, you can create a usable geometry file with:

```
$ legend-pygeom-l200 --fiber-modules=detailed l200.gdml
```

The generated geometry can be customized with a large number of options. Some
geometry options can both be set on the CLI utility and on the config file.
Those are described in {doc}`cfg-geometry`, but the descriptions similarly
applies to the CLI options.

:::{note}

In the new simulation flow architecture introduced with _remage_, the geometry
GDML file contains a single static geometry. For each combination of geometry
options (common example: different
{doc}`calibration source positions <cfg-calibration>`), a separate GDML file has
to be created.

:::

```{toctree}
:maxdepth: 3

runtime-cfg
vis
description
coordinate_systems
```

```{toctree}
:maxdepth: 1
:caption: Development

geom-dev
naming
Package API reference <api/modules>
```
