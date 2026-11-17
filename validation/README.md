# Catalog Validation

Scripts that validate the entries of the instrument catalog (`src/tatc/resources/instruments/*.yml`) against reference data from operational missions, following the [validation standard](../docs/validation/standard.rst). Their results are published on the catalog page of the documentation, which `docs/validation/catalog_page.py` generates from the entries and their records.

The validation notebooks, which test orbits, analyses, and the modeling choices behind several entries, are in [`docs/validation`](../docs/validation).

## Layout

- `catalog/validate_catalog_<entry>.py`: one script per entry (some cover several entries' platforms or modes).
- `catalog/validation_io.py`: shared archive access, element sets, geometry, scoring, and the record and fixture writers.
- `catalog/gesdisc_scanner.py`: reading and fitting scanner swaths through `earthaccess` (NASA DAACs).
- `catalog/earthcare_io.py`: reading EarthCARE products from ESA MAAP.
- `catalog/compare_swath_time_steps.py`: compares swaths at different time steps on the unit-test fixtures (no reference data).
- `catalog/elements/`: element sets the scripts read when the archives do not supply them.

The scripts import each other by module name, so run them as scripts (not as a package), and keep them in one folder.

## Running

Install the `validation` extra (`pip install -e ".[validation]"`), then run a script with a data directory and, for most, a day:

    python validation/catalog/validate_catalog_atms.py validation/data 2026-10-02

The usage of each script, including the days on which its entry was fitted and tested, is in its docstring. Downloaded data are cached in the data directory; `validation/data` is ignored by git.

Common options:

- `--fit`: measures the fitted values of an entry's configurations on the day of the fit, printing them without writing.
- `--write`: writes the entry's validation record (`src/tatc/resources/validation/<entry>.json`) and the reduced fixture from which unit tests recompute its results (`tests/catalog/fixtures/<entry>.npz`).

Without `--write`, a script prints its results. The exception is `validate_catalog_abi.py`, which takes only a data directory and always writes.

## Accounts

- NOAA Open Data Dissemination (JPSS, GOES, Himawari) and Microsoft Planetary Computer: anonymous.
- NASA DAACs: a [NASA Earthdata Login](https://urs.earthdata.nasa.gov/) account, as a token in `EARTHDATA_TOKEN` or, failing that, the credentials of `urs.earthdata.nasa.gov` in `~/.netrc`.
- EarthCARE (ATLID, BBR, CPR, MSI): an ESA MAAP offline token in `ESA_MAAP_KEY`.

Each script's docstring states what it needs.
