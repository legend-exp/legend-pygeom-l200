# Description of the geometry

This section describes briefly the geometry, with a particular focus on the
names of the physical volumes of the various sources.

:::{note} The renderings on this page use the public geometry (see
{doc}`cfg-geometry`). :::

The geometry is divided in various assemblies, as described in {doc}`vis`.

## HPGe strings

The `strings` assembly consists of the HPGe detectors and their support
structure.

This is shown below:

```{subfigure} ABC
:subcaptions: above

:::{image} ./images/hpge_strings.png
:height: 400px
:alt: HPGe detector strings.
:::

:::{image} ./images/holders.png
:height: 400px
:alt: PEN and Copper HPGe support structure.
:::

:::{image} ./images/nylon.png
:height: 400px
:alt: Nylon minishrouds and calibration tubes.
:::

```

&nbsp;

The left figure above shows the HPGe detectors, as explained in {doc}`naming`
these are given physical volume names of the **detector name** (i.e. `V99000A`).

The center shows the detector holder structure. This is divided into several
parts:

- the pen baseplates for each detector, with physical volume names `pen_{NAME}`
  where `{NAME}` is the detector name :::{note} These are shown in yellow,
  green, red and blue in the rendering depending on the size :::
- for some detectors there are also top pen rings named `pen_top_{NAME}`,
- the HPGe detectors are supported by a copper support structure, all parts of
  this structure are made of electroformed copper and have physical volumes
  prefixed with `hpge_support_structure`.
- nylon minishrouds surrounding each string, these have a name starting in
  `minishroud`, and calibration tubes which have names `calibration_tube_{IDX}`
  where `{IDX}` is the index of the SIS (see {doc}{cfg-`calibration`}).

:::{warning} Some components of the HPGe readout chain, (electronics cables
etc.) are not yet implemented! :::

The hpge copper support structure consists of three components:

- a copper rod supporting each string (shown at the top of the rendering), these
  have names:

```
hpge_support_copper_string_support_structure_string_{STRING}
```

where `{STRING}` is the string number,

- a triangular copper support for each string (or "tristar") named

```
hpge_support_copper_tristar_{SIZE}_string_{STRING}
```

where `{SIZE}` is the size of the string (`small`,`medium`, `large`, or
`xlarge`).

- copper rods for each HPGe detector which have names:

```
hpge_support_copper_string_{STRING}_cu_rod_{IDX}
```

where `{IDX}` is an index of the rod.

:::{tip}

To select all copper parts in _remage_ you can use a wildcard
`hpge_support_copper.*` and similarly to select all copper rods, tristar or
string support structures. :::

## Top plate

The `top` assembly consists of the top plate holding the CC4 electronics,
currently this is a single physical volume called `top_plate`.

This is shown on the left figure of the rendering below.

```{subfigure} AB
:subcaptions: above

:::{image} ./images/top_plate.png
:height: 400px
:alt: Copper top plate
:::

:::{image} ./images/fibers.png
:height: 400px
:alt: Fiber shrouds.
:::

```

&nbsp;

## Fiber's and SiPMs

The fiber shroud's for the LAr readout are shown in the right figure above.

As mentioned in {doc}`runtime-cfg` there are two modes for the optical fibers,
either individual fibers or a segmented option.

In both cases the fiber volumes are divided into 4 parts:

- an outer TPB coating,
- two layers of cladding,
- the fiber core.

The optical fiber system consists of a large number of physical volumes, to
enable concise _remage_ macros the names are first prefixed with:

- `fiber_inner`: for the inner barrel fibers (inner cylinder),
- `fiber_outer`: for the outer barrel fibers.

Next the name contains an identifier of the part of the fiber either, `tpb`,
`cladding1`, `cladding2` or `core`.

The rest of the names give further information on the fiber, and uses the length
of the fiber, whether it is part of the lower bend for the outer barrel and the
fiber index to obtain a unique physical volume name.

:::{tip} For more users it is expected to use wildcards to select groups of
optical fibers. :::

The SiPMs are named after the **detector name** (ie `S001`), and there are also
physical volumes of a wrapping around the SiPMs with names `{NAME}_wrap` (where
`{NAME}` is the SiPM name).

## Wavelength shifting reflector

The array is emersed in liquid argon (physical volume named `lar`) and
surrounded by a tetratex and TPB lined copper foil for reflecting scintillation
light (WLSR).

This is shown as the outer cylinder below, which shows all the main components
of the experiment. The WLSR is the shown as the outer grey component. :::{image}
./images/wlsr.png :height: 1000px :alt: Rendering of the full experiment
including outer WLSR. :::

This consists of three physical volumes:

- `wlsr_outer` the copper foil,
- `wlsr_tetratex` the tetratex coating,
- `wlsr_tpb` the TPB coating.
