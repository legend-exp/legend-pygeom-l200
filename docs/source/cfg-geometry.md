# Geometry options

geometry options are both available in the CLI interface and the config file
interface. Provided CLI options override the values from the config file.

```{important}
Always use `fiber_modules: detailed` (in the config file) or
`--fiber-modules detailed` (on CLI) for optical simulations.

**The simplified geometry has a different (and certainly wrong) light-guiding
efficiency** - it is only suitable for visualisation or simulations not requiring
full optical physics.
```

example config (not the default values!):

```yaml
# select parts of the geometry to be built.
assemblies:
  - fibers
  - strings

# Select the fiber shroud model, either coarse segments or single fibers.
fiber_modules: detailed # or "segmented"

# Select the PMT configuration of the muon veto.
pmt_config: LEGEND200 # or "GERDA"

# If true, build from public testdata only.
public_geom: false
```

The list of assemblies provided in the `assmblies` config option or the
`--assemblies` CLI option can be either:

- a list of assemblies to build: `strings,calibration` will build the argon
  cryostat (this is always bullt), the HPGe strings and the calibration system
  only.
- a list of operations modifying the default assembly list: `+watertank` will
  build all default assemblies (i.e. the argon cryostat and anything in it), and
  dd also the muon veto watertank.

Specifying assemblies in the config file works similarly, but requires the list
to be passed as an actual JSON/YAML list instead of a comma-delimited string.
