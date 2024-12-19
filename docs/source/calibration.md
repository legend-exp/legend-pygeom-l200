# Calibration sources

**legend-pygeom-l200** can easily be used to include the calibration sources
deployed with the SIS.

As the position and types of sources frequently change, they are not hard-coded
in the python code, but can be configured in a runtime configuration file.

**TL;DR**: a working example:

```yaml
sis:
  "1":
    sis_z: 8250
    sources: ["Th228", null, null, "Th228"]
  "2": null
  "3": null
  "4": null
```

## calibration source configuration

The `sis` config object contains objects describing the deployed sources in each
of the four SIS tube (or `null`, if no sources are deployed in that tube). The
numbering equals the "official" SIS numbering scheme.

The `sis_z` coordinate is the SIS reading of the deployed SIS strings. The
coordinate is transformed into simulation coordinates as described on [the
Confluence page][confluence-coord].

Each SIS string has four slots for sources, that can be filled differently with
this tool. The `sources` array contains the four slots from the top to the
bottom. The bottom source is seated on top of the tantalum absorber.

Different types of sources can be included in the slots:

- `Th228` — a normal LEGEND calibration source (as described in [L. Baudis _et
  al_ 2023 _JINST_ 18 P02001][citation-source]).
- `Ra` — a special calibration source.
- ...`+Cu` — add a copper absorber cap to any other source. The dimensions of
  the cap can be seen in the code.
- `null` — no source in this slot. The tantalum absorber is placed irrespective
  of whether a source is placed inside in the slot.

> [!NOTE]
>
> The generated geometry does **not contain the requested source material**. The
> decaying isotope has to be configured in the user's Geant4/remage macro file.
>
> The volumes named `source_inner_sis{SIS number}_source{slot}` can be used as
> the confinement volumes in remage (i.e. with a regex `^source_inner_.*`, if
> all sources share an isotope).

## extra source outside SIS

It is also possible to add an extra source, that is not part of any SIS string,
via a config file. all notes on SIS sources above also apply to this extra
source.

```yaml
extra_source:
  position_in_mm: [0, 0, 0]
  name: "_central" # to identify in geometry, this produces a volume `source_inner{name}`
  source: "Th228+Cu"
```

[confluence-coord]:
  https://legend-exp.atlassian.net/wiki/spaces/LEGEND/pages/1111785478/Calibration+simulations#Source-geometry-%2F-position
[citation-source]: https://doi.org/10.1088/1748-0221/18/02/P02001
