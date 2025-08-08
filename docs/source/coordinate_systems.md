# Coordinate systems

## Internal coordinate system inside LAr volume

In the `hpge_strings`, `top` (and partly also `calibration`) modules, a special
internal coordinate system is used for some calculations. Unlike in many other
conventions, this is a **left handed coordinate system**! It is a cylindrical
coordinate system with an angle {math}`\phi`, a radius {math}`r` and a {math}`z`
coordinate.

- {math}`z = 0` denotes the top of the top plate in the ’birds nest‘ The
  positive {math}`z` direction extends upwards (in the same direction as the
  GDML Z-axis)
- {math}`\phi` is the **clockwise** angle. {math}`\phi = 0` denotes the
  direction of string 11, i.e. ’warm north‘.
- {math}`r` is measured to the central axis, as one would expect.

To convert this coordinate system, one can use:

```{math}
Z = z + z0

X = r \cdot \cos\phi

Y = - r \cdot \sin\phi
```

X/Y/Z is then a typical right handed cartesian coordinate system that is
centered at the barycenter of the argon volume. This is an intrinsic property of
the GDML file format and cannot be changed.

Important: The `fibers` module internally uses a normal, "right handed"
coordinate system. But this module does not consume external data. In technical
documents concerning the fiber module positioning, the above left handed
coordinate system is used.

## Global coordinate system

TBD
