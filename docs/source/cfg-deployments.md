# LEGEND-200 deployments

_legend-pygeom-l200_ supports the geometry of multiple LEGEND-200 deployments.

If you manage your own _legend-metadata_ checkout, it is strictly necessary to
have a _legend-metadata_ checkout that already contains the channelmap for this
deployment. For still untagged or pre-release versions of _legend-metadata_, you
have to manually perform the checkout.

The configuration is adapted by selecting a timestamp which controls the
metadata applied to construct the geometry:

```yaml
metadata_timestamp: 20250425T180115Z # example: l200cfg08 / p14
```

(supported-deployments)=

## supported deployments `l200cfgNN`

The naming convention is `l200cfgNN`, where `cfg` stands for "configuration" and
`NN` is a zero-padded integer. The table summarizes the defined configurations,
together with the DAQ periods are just included for reference and should not be
used to refer to simulation configurations.

| Sim HW config | DAQ periods   | timestamp after    | legend-metadata |
| ------------- | ------------- | ------------------ | --------------- |
| **l200cfg01** | **p03–p09**   | `20230311T235840Z` | ≥ v1.0.0        |
| **l200cfg02** | **p10**       | `20240226T104834Z` | ≥ v1.0.0        |
| **l200cfg03** | **p11**       | `20240411T000000Z` | ≥ v1.0.0        |
| **l200cfg04** | **p13-r000**  | `20241202T150000Z` | ≥ v1.0.1        |
| **l200cfg05** | **p13-r001**  | `20241210T225016Z` | ≥ v1.0.1        |
| **l200cfg06** | **p13-r002+** | `20241217T095149Z` | ≥ v1.0.1        |
| **l200cfg07** | **—**         | not implemented    |                 |
| **l200cfg08** | **p14**       | `20250425T180115Z` | ≥ v1.0.2        |
| **l200cfg09** | **p15–p16**   | `20250716T161517Z` | ≥ v1.0.2        |
| **—**         | **p18–**      | `20251108T002705Z` | ≥ v1.1.0        |

- `l200cfg01`: first LEGEND-200 configuration with 142 kg of enriched germanium.
  Data taking periods p03 (start: March 2023) to p09
  ([PRL 2025](https://doi.org/10.1103/25tk-nctn)).

- `l200cfg02`: same as `cfg01`, but with outer fiber barrel removed. Data taking
  period p10 (start: Feb 2024).

- `l200cfg03`: same as `cfg02`, but with nylon mini-shrouds removed in addition.
  Data taking period p11 (start: Apr 2024).

- `l200cfg04`: configuration with _only the inner barrel of the liquid argon
  instrumentation deployed_, no HPGe array. Data taking period p13, `r000` (Dec
  2024).

- `l200cfg05`: configuration with _both inner and outer barrels of the liquid
  argon instrumentation deployed_, no HPGe array. Data taking period p13, `r001`
  (Dec 2024).

- `l200cfg06`: same as `cfg05`, but with some new SiPM modules. Since new
  detectors are installed, the configuration is technically different and we
  bump the experiment label. Data taking period p13, from `r002` (start:
  Dec 2024) on.

- `l200cfg07`: 114 kg HPGe array (start: March 2025). No good data was taken
  with this configuration.

- `l200cfg08`: 99 kg HPGe array. Data taking period p14 (start: Apr 2025).

- `l200cfg09`: 138 kg HPGe array. Data taking period from p15 (start: July 2025)
  onwards. Starting with p18, some NMSs were exchanged with a slightly different
  length (no new configuration number was issued).

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
