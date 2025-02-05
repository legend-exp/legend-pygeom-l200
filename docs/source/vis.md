# Visualization of the geometry

## Visualization with `legend-pygeom-l200`

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

## Visualizing with Geant4/[`remage`](https://github.com/legend-exp/remage) (_advanced_)

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

## Adjusting the visualization from python

See the
[legend-pygeom-tools docs](https://legend-pygeom-tools.readthedocs.io/en/latest/).
