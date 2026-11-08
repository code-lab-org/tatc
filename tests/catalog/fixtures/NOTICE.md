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
| `oli.npz` | Imaging times of points in Landsat 9 scenes of 3 October 2026, from the angle coefficient files of Landsat Collection 2 Level-2 products | [USGS](https://www.usgs.gov/landsat-missions), through [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/dataset/landsat-c2-l2) | `docs/validation/validate_catalog_oli.py` |
| `tropomi.npz` | Sentinel-5P TROPOMI aerosol index product (near real time), 5 November 2026. Contains modified Copernicus Sentinel data 2026. | [Sentinel-5P Level 2 on the Registry of Open Data on AWS](https://registry.opendata.aws/sentinel5p/) | `docs/validation/validate_catalog_tropomi.py` |
| `iasi.npz` | Pixel positions of the Metop-B IASI sea surface temperature product, 4 November 2026. Contains modified EUMETSAT OSI SAF IASI sea surface temperature data 2026, under the CC BY 4.0 licence. | [NASA PO.DAAC](https://podaac.jpl.nasa.gov/dataset/IASI_SST_METOP_B-OSISAF-L2P-v1.0) | `docs/validation/validate_catalog_iasi.py` |
| `msi.npz` | Imaging times of points of the viewing angle grids of Sentinel-2B Level-2A tiles, 12 August 2026. Contains modified Copernicus Sentinel data 2026. | [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/dataset/sentinel-2-l2a) | `docs/validation/validate_catalog_msi.py` |
| `csar.npz` | Zero-Doppler times of points of the geolocation grids of Sentinel-1C Level-1 GRD products, 12 August 2026. Contains modified Copernicus Sentinel data 2026. | [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/dataset/sentinel-1-grd) | `docs/validation/validate_catalog_csar.py` |
| `gmi.npz` | `GPM_1CGPMGMI` version 08, 20 July 2026; footprint sizes published in the GMI L1B ATBD and TC PRIMED documentation | [NASA GES DISC](https://disc.gsfc.nasa.gov/datasets/GPM_1CGPMGMI_08/summary) | `docs/validation/validate_catalog_gmi.py` |

Each fixture also holds two-line element sets, used to propagate the satellite's
orbit.
