# legend-pygeom-l200

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> [!WARNING]
>
> This is an early version of the LEGEND-200 geometry implemented with the
> python-based simulation stack. It is not a drop-in replacement for MaGe, and
> still under heavy development!

## Installation and usage

This package rerquires a working setup of
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

## Visualization of the geometry

### Adjusting the visualization from `legend-pygeom-l200`

On a logical volume instance, you can set `pygeom_color_rgba`, e.g.

```python3
lv = g4.LogicalVolume(...)

# hide this volume in the visualization
lv.pygeom_color_rgba = False

# set the vis coloring to the given RGBA value. All 4 numbers should be given in the range 0–1.
lv.pygeom_color_rgba = (r, g, b, a)
```

### Visualizing with Geant4/[`remage`](https://github.com/legend-exp/remage)

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

and use it with remage `remage vis.mac -g l200.gdml`. It will validate that the
given GDML file can be read by Geant4 and show a visualization from it.

## Further features

### Registering detectors for use with [`remage`](https://github.com/legend-exp/remage)

On a physical volume instance, you can set `pygeom_active_dector`, e.g.

```python3
from .det_utils import RemageDetectorInfo

pv = g4.PhysicalVolume(...)

# attach an active detector to this physical volume.
lv.pygeom_active_dector = RemageDetectorInfo("optical", 1)
```

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
