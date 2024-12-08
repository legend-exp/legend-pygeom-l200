# Naming conventions

This document describes our naming convention for geometry parts (solids,
volumes, materials, surfaces, ...).

- **volume names**: snake_case, e.g. `lar` or `
  - exception: **detector names** are used as-is, e.g. `V01234A` or `S035`
- the names of **corresponding solids, logical and physical volumes** in GDML
  should usually be the same
  - exception: multiple placements of one LV need unique multiple physical
    volume names
  - this does not always translate to _python variable names_. It might be
    necessary to use unique names for logical and physical volume instances in
    code (i.e. to attach surfaces)
  - python variables names are not generally expected to follow the geometry
    names, but should if possible.
- **surfaces**
  - `surface_{from}_to_{to}` for OpticalSurfaces (property definition)
  - `bsurface_{from}_{to}` for border surfaces `* ssurface_{to}` for skin
    surfaces
- **materials**
  - snake_case similar to volumes, e.g. `metal_copper`.
  - Elements use capitalized names, e.g. `Hydrogen`
