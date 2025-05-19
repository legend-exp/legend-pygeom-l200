# Features for geometry development

## Registering detectors for use with [`remage`](https://github.com/legend-exp/remage)

See the
[legend-pygeom-tools docs](https://legend-pygeom-tools.readthedocs.io/en/latest/).

This information can be exported by using `--det-macro-file=l200-dets.mac` as an
additional CLI option. This macro then should be `/control/execute`d in your
main macro.

## Checking for overlaps

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
