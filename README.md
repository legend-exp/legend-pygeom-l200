# legend-pygeom-l200

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> [!WARNING]
>
> This is an early version of the LEGEND-200 geometry implemented with the
> python-based simulation stack. It is not a drop-in replacement for MaGe, and
> still under heavy development!

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
