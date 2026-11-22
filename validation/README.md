# Catalog Validation

Scripts that validate the entries of the instrument catalog (`src/tatc/resources/instruments/*.yml`) against reference data from operational missions, following the [validation standard](../docs/validation/standard.rst), and the systems of the catalog (`src/tatc/resources/systems/*.yml`: satellite programs, radar networks, and surface networks) against element sets and station metadata, following its system standard. Their results are published on the catalog and systems pages of the documentation, which `docs/validation/catalog_page.py` and `systems_page.py` generate from the entries and their records.

The validation notebooks, which test orbits, analyses, and the modeling choices behind several entries, are in [`docs/validation`](../docs/validation).

## Layout

- `catalog/validate_catalog_<entry>.py`: one script per entry (some cover several entries' platforms or modes).
- `catalog/validation_io/`: helpers shared by every script, by topic:
  - `cli.py`: the shared command line (`parse_run`), with defaults from the catalog;
  - `archives.py`: HTTP with retries, range requests, the JPSS archives, the NASA Earthdata login, and element sets (`get_element_set`);
  - `geometry.py`: coordinates, and the lattices of samples that scanners image (`ScanLattice`);
  - `scoring.py`: comparisons through public TAT-C methods (passes, epochs, sample sizes, swath areas);
  - `records.py`: the validation records and reduced test fixtures (`write_results`).
- `catalog/scanner.py`: cross-track scanners whose products give each sample's position (reading NASA DAAC products, fitting the scan, and the sample size and swath tests, `validate_beams` and `validate_swath`).
- `catalog/conical.py`: conical scanners (AMSR2, GMI, SMAP, SSMIS, TMI).
- `catalog/sar.py`: side-looking radars (C-SAR and L-SAR).
- `catalog/earthcare_io.py`: reading EarthCARE products from ESA MAAP, and the profiler test.
- `catalog/compare_swath_time_steps.py`: compares swaths at different time steps on the unit-test fixtures (no reference data).
- `catalog/elements/`: element sets: `element_sets.yml`, those used on days for which no archive or file is used, and files of the history of some platforms' element sets.

The scripts import each other by module name, so run them as scripts (not as a package), and keep them in one folder.

## Running

Install the `validation` extra (`pip install -e ".[validation]"`), then run a script with a data directory:

    python validation/catalog/validate_catalog_atms.py validation/data

Downloaded data are cached in the data directory; `validation/data` is ignored by git. Every script takes the same options:

- `DAY ...`: the day of the data (UTC, such as `2026-10-02`, or a time with a zone for geostationary scans, such as `2026-10-03T12:00Z`). By default, the day of the entry's record, so that a run without a day repeats the recorded test; with `--fit`, the day of the data to which the values were fitted.
- `--platform P`: the platform whose data are compared (one of the entry's). On a platform without configurations of its own, the run is a transfer: its results are added to the record of the configurations it tests, and do not count toward the tier. `--transfer` makes a run a transfer on a platform with configurations of its own (such as NOAA-20's configurations on NOAA-21).
- `--configuration KEY`: tests only that configuration (repeatable).
- `--fit`: measures the fitted values of the configurations, printing them without writing.
- `--write`: writes the entry's validation record (`src/tatc/resources/validation/<entry>.json`) and the reduced fixture from which unit tests recompute its results (`tests/catalog/fixtures/<entry>*.npz`). It warns if the day is fewer than 30 days from the fitting data, which the entry's validation reports as a gap.
- `--elements FILE`: takes the platform's element set from a file of element sets.

Without `--write`, a script prints its results. Some scripts have options of their own (such as `--fit-view` or `--download-only`), which `--help` lists with the script's description: what it compares, the days on which its entry was fitted and tested, and what it needs.

Element sets come from, in order: `--elements`; the set registered for the platform and day in `catalog/elements/element_sets.yml`; NOAA's operational element set of a JPSS satellite (whose archive may no longer hold days of earlier records); or the file in `catalog/elements/` that holds the platform's sets (the set nearest the middle of the day). To test a day for which none of these has a set, register one, or give a file.

## Systems

Systems have scripts of their own, with a command line like the instruments' (`validation_io/cli.py`, `parse_system_run`): `SCRIPT SYSTEM DIRECTORY [DAY ...] [--configuration KEY ...] [--fit] [--write]`, which caches downloads in `DIRECTORY` (such as `validation/data/elements`), writes records to `src/tatc/resources/validation/systems/<system>.json` with `--write` (warning if fitted values are fewer than 30 days from the data), and with `--fit` fits values or builds a table instead of testing:

- `catalog/validate_system_orbits.py SYSTEM DIRECTORY [DAY ...] [--step N] [--end DAY]`: the orbits of a satellite system against element sets, each standing for the day of its own epoch, on each spacecraft's days at least 30 days after its fit to `--end` (by default yesterday, or a retired spacecraft's last day): the registered days, the days of the files of element sets, and, for the JPSS satellites, every `--step` days of NOAA's archives, which keep about the last 90 days, weekdays only). Its metrics are in `validation_io/orbits.py`: the mean elements of each set, the mean longitude of a geostationary satellite, and the ascending nodes within 6 hours of each set's epoch compared with those of an orbit maintained on a repeat ground track; with `--fit`, the node of such an orbit is the middle of the nodes of the sets of its fitting days (`fit_ascending_node`).
- `catalog/validate_system_nexrad.py nexrad DIRECTORY [DAY ...]`: on the week from 30 days after the table's snapshot by default, the NEXRAD table (built with `--fit` and a day from the radars' own Level II metadata, for the WSR-88D stations of NCEI's table) against NCEI's table and the stations recording Level II volumes, and each configuration against sampled whole volumes of its volume coverage patterns (about 1 GB). Helpers are in `validation_io/stations.py`.
- `catalog/validate_system_surface.py SYSTEM DIRECTORY` (`uscrn`, `asos`, or `awos`, as of the day of the run): the table (built with `--fit` from NCEI's HOMR) against NCEI's GHCN-Daily station inventory and the stations the Iowa Environmental Mesonet reports online on the day of the test.

To test a satellite on a day without archived or registered element sets, register one in `catalog/elements/element_sets.yml` (the current sets of CelesTrak, for example).

## Adding an entry

Write the entry's YAML file, then a script `validate_catalog_<entry>.py` that reads the reference data and compares them with the entry through public TAT-C methods: `run = parse_run("<entry>", __doc__)` gives the run, `run.element_set()` its element set, and `write_results(run, key, reference=..., metrics=..., fixture=...)` writes each configuration's results when `run.write`. Scanners whose products give each sample's position can use `scanner.Setup` with `validate_beams` and `validate_swath`, as the scripts of ATMS and DPR do. Then add the fixture's tests to `tests/catalog`.

## Accounts

- NOAA Open Data Dissemination (JPSS, GOES, Himawari), the NEXRAD Level II archive on AWS, NCEI, the Iowa Environmental Mesonet, and Microsoft Planetary Computer: anonymous.
- NASA DAACs: a [NASA Earthdata Login](https://urs.earthdata.nasa.gov/) account, as a token in `EARTHDATA_TOKEN` or, failing that, the credentials of `urs.earthdata.nasa.gov` in `~/.netrc`.
- EarthCARE (ATLID, BBR, CPR, MSI): an ESA MAAP offline token in `ESA_MAAP_KEY`.

Each script's docstring states what it needs.

## Geoid

`geoid/build_geoid_grid.py` builds the 0.25-degree grid of EGM2008 geoid heights bundled with TAT-C (`src/tatc/resources/egm2008_15.npz`, read by `tatc.utils.compute_geoid_height`) from the 2.5-minute EGM2008 grid of the National Geospatial-Intelligence Agency (public domain), as distributed by PROJ (https://cdn.proj.org/us_nga_egm08_25.tif). It reports the interpolation errors of grids of other spacings and compares the bundled grid with PROJ's transformation of heights:

    python validation/geoid/build_geoid_grid.py path/to/us_nga_egm08_25.tif
