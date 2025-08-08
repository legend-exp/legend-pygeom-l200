# Visualization of the geometry

## Visualization with `legend-pygeom-l200`

Simply use `legend-pygeom-l200 -V [...]` to visualize the full geometry. See
{ref}`pygeomtools:gdml-viewer` for details on how to use the interactive viewer.

If you want to exclude/include components from the 3D rendering, append
`--assemblies=...`. See the section on {doc}`cfg-geometry`.

The cryostat and LAr volumes are always part of the output, but can be hidden in
a scene file, if necessary.

The visualization with the VTK-based viewer can be customized with a scene file,
including options for

- changing volume colors, transparency & hiding volumes
- better looking renderings (lighting, camera positions, file export)
- adding a clipper of the geometry (e.g., hiding one half of all volumes)
- and more tools, also for basic event visualization.

See the {ref}`pygeomtools:scene-file` for a reference of the scene file. The
scene file can be passed as an argument after `-V`:
`legend-pygeom-l200 -V scene.yaml [...]`.

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
