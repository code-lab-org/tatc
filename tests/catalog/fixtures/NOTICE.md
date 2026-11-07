# Test Fixture Data

These fixtures are small samples derived from public data, used to detect
regressions without network access. The repository's license does not cover the
underlying data, whose providers' terms apply.

| Fixture | Derived from | Provider | Written by |
| --- | --- | --- | --- |
| `abi.npz` | GOES-19 ABI full-disk products (ACM, SST, LST) of 3 October 2026 | [NOAA Open Data Dissemination](https://registry.opendata.aws/noaa-goes/) | `docs/validation/validate_catalog_abi.py` |
| `atms.npz` | NOAA-20 ATMS SDR geolocation, 12 August 2026 | [NOAA Open Data Dissemination](https://registry.opendata.aws/noaa-jpss/) | `docs/validation/validate_catalog_atms.py` |
| `viirs.npz` | NOAA-20 VIIRS moderate-resolution geolocation, 12 August 2026 | [NOAA Open Data Dissemination](https://registry.opendata.aws/noaa-jpss/) | `docs/validation/validate_catalog_viirs.py` |
| `modis.npz` | Terra MODIS geolocation (`MOD03`, Collection 7), 5 November 2026 | [NASA LAADS DAAC](https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD03) | `docs/validation/validate_catalog_modis.py` |
| `gmi.npz` | `GPM_1CGPMGMI` version 08, 20 July 2026; footprint sizes published in the GMI L1B ATBD and TC PRIMED documentation | [NASA GES DISC](https://disc.gsfc.nasa.gov/datasets/GPM_1CGPMGMI_08/summary) | `docs/validation/validate_catalog_gmi.py` |

Each fixture also holds two-line element sets, used to propagate the satellite's
orbit.
