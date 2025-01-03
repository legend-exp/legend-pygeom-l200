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

## Runtime configuration

For often-changing details of the geometry are configured using a runtime
configuration JSON file. This file is specified using the `--config $FILE`
parameter.

Detailed information about the configurable subsystems is available:

- [Calibration system (SIS)](docs/source/calibration.md)

## Visualization of the geometry

### Visualization with `legend-pygeom-l200`

Simply use `legend-pygeom-l200 -V [...]` to visualize the full geometry.

If you want to exclude components from the 3D rendering, append
`--assemblies=...`. Possible values are:

- `strings` (the whole HPGe array)
- `fibers`. It is highly recommended to also append the argument
  `--fiber-modules=segmented` to avoid rendering all single fibers, if you only
  need to see the overall shape.
- `calibration` (calibration tubes and sources, if any)
- `top` (copper top plate)
- `wlsr`

Multiple values can be combined with commas. Example:
`--assemblies=strings,calibration`.

The cryostat and LAr volumes are always part of the output.

### Visualizing with Geant4/[`remage`](https://github.com/legend-exp/remage) (_advanced_)

The visualization can be exported to Geant4 by using `--vis-macro-file=`:
`legend-pygeom-l200 --vis-macro-file=l200-vis.mac l200.gdml [...]`.

This generated macro does not start any visualization on its own, it just sets
the colors. To use it, create a file `vis.mac` in the same directory:

```
/run/initialize

/vis/open OGL
/vis/drawVolume lar

/vis/viewer/set/defaultColour black
/vis/viewer/set/background white
/vis/viewer/set/viewpointVector -3 -2 1
/vis/viewer/set/upVector 0 0 1
/vis/viewer/set/rotationStyle freeRotation
/vis/viewer/set/lineSegmentsPerCircle 100

/vis/scene/add/trajectories smooth
/vis/scene/endOfEventAction accumulate

# import the auto-generated visualization attributes from legend-pygeom-l200.
/control/execute l200-vis.mac
```

and use it with remage `remage vis.mac -i -g l200.gdml`. It will validate that
the given GDML file can be read by Geant4 and show a visualization from it.

It is also possible to use `--assemblies=` as described above. This will remove
any non-specified assembly from the output GDML file. Make sure that you do not
overwrite any "production" geometry with this command. Using a file with
stripped-down assemblies for a simulation will probably give wrong results.

### Adjusting the visualization from python

See the
[legend-pygeom-tools docs](https://legend-pygeom-tools.readthedocs.io/en/latest/).

## Further features (for developers)

### Registering detectors for use with [`remage`](https://github.com/legend-exp/remage)

See the
[legend-pygeom-tools docs](https://legend-pygeom-tools.readthedocs.io/en/latest/).

This information can be exported by using `--det-macro-file=l200-dets.mac` as an
additional CLI option. This macro then should be `/control/execute`d in your
main macro.

### Checking for overlaps

Using `--check-overlaps` might yield wrong results (it uses the coarsely
tessellated volumes also used for visualization); also it is very slow. Using
Geant4 to load the generated GDML file will give you correct results.

Create a file called `check-overlaps.mac` with the following contents:

```
/RMG/Manager/Logging/LogLevel error
/run/initialize
```

and use it with remage `remage check-overlaps.mac -g $PATH_TO_YOUR_GDML_FILE`.
It will validate that the given GDML file can be read by Geant4 and that it has
no overlaps.
