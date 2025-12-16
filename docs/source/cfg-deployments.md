# LEGEND-200 deployments

_legend-pygeom-l200_ supports the geometry of multiple LEGEND-200 deployments.

If you manage your own _legend-metadata_ checkout, it is strictly necessary to
have a _legend-metadata_ checkout that already contains the channelmap for this
deployment. For still untagged or pre-relase versions of _legend-metadata_, you
have to manually perform the checkout.

|             | timestamp after    | legend-metadata |
| ----------- | ------------------ | --------------- |
| **p03–p09** | `20230311T235840Z` | >= v1.0.0       |
| **p10**     | `20240226T104834Z` | >= v1.0.0       |
| **p11**     | `20240411T000000Z` | >= v1.0.0       |
| **p13**     | `20241210T225016Z` | >= v1.0.1       |
| **p14**     | `20250425T180115Z` | >= v1.0.2       |
| **p15–p16** | `20250716T161517Z` | >= v1.0.2       |
| **p18–**    | `20251108T002705Z` | >= v1.1.0       |

the timestamp can be selected in the config file:

```yaml
metadata_timestamp: 20250425T180115Z # example: p14
```

:::{tip}

adding the geometries for new deployments to _legend-pygeom-l200_ require only a
few number of changes:

1. a channelmap in _legend-metadata_ (not managed by the simulation team)
2. a new `extra_meta` file in _legend-pygeom-l200_ that contains
   - string positions & parameters such as for minishrouds
   - copper rod lengths for each detector unit
   - pen plate sizes
   - calibration tube parameters
3. an implementation of any geometry change that cannot be configured in the
   extra metadata.

:::
