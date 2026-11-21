"""
Validates the catalog's surface networks (systems `uscrn`, `asos`, and `awos`)
against NCEI's station metadata and the stations reporting observations, as of
the day of the run (UTC).

* Siting: the table's positions and ground elevations compared with NCEI's
  inventory of the Global Historical Climatology Network - Daily (GHCN-D), a
  product separate from the table's source, though compiled from the same
  station histories (NCEI's Historical Observing Metadata Repository, HOMR).
* Membership: the stations that the Iowa Environmental Mesonet (IEM) reports
  online (observations received recently: USCRN, or the METAR stations of the
  United States' networks) among the network's stations as listed on the day of
  the run (as the table is built), compared with the table's stations. As the
  table is built from the same listing, recall measures the stations added to the
  network since the table's snapshot, and precision the table's stations that no
  longer report, so a test is held out from the snapshot by
  `MIN_HOLDOUT_GAP_DAYS`.

With `--fit`, it builds the table instead (`resources/systems/tables/<system>.csv`)
from HOMR: for USCRN, the operational stations of NCEI's USCRN list (with a WBAN
identifier), at their HOMR positions and ground elevations; for ASOS and AWOS,
the stations of the United States and its territories in HOMR's tables, with
known elevations. Elevations are meters above mean sea level. Downloads are
cached in a folder of the day.

Usage:
    python validate_system_surface.py uscrn DIRECTORY --fit
    python validate_system_surface.py uscrn DIRECTORY [--write]

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import io
import json
import pathlib
import re
from datetime import date, datetime, timezone

import numpy as np
import pandas as pd

from tatc.catalog import SystemClass, SystemEntry
from validation_io.cli import SystemRun, parse_system_run
from validation_io.records import (
    REPOSITORY,
    p95_metric,
    system_result,
    value_metric,
    write_system_record,
)
from validation_io.stations import (
    FEET,
    cached_download,
    horizontal_errors,
    membership,
    read_homr_table,
)

TABLES = REPOSITORY / "src" / "tatc" / "resources" / "systems" / "tables"
USCRN_LIST = "https://www.ncei.noaa.gov/pub/data/uscrn/products/stations.tsv"
HOMR_SERVICE = "https://www.ncei.noaa.gov/access/homr/services/station/search"
GHCND = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
# GHCN-D's elevation of a station without one
GHCND_MISSING = -999.9
IEM_NETWORKS = "https://mesonet.agron.iastate.edu/api/1/networks.json"
IEM_NETWORK = "https://mesonet.agron.iastate.edu/geojson/network/{}.geojson"
# the countries of HOMR's tables in the United States and its territories
UNITED_STATES = {
    "UNITED STATES",
    "PUERTO RICO",
    "GUAM",
    "VIRGIN ISLANDS (U.S.)",
    "AMERICAN SAMOA",
    "NORTHERN MARIANA ISL",
}
REFERENCE = (
    "NCEI's GHCN-Daily station inventory (positions and ground elevations) and "
    "the stations the Iowa Environmental Mesonet reports online"
)


def today() -> date:
    """The day of the run (UTC)."""
    return datetime.now(timezone.utc).date()


def default_days(entry: SystemEntry, fitting: bool) -> list[date]:
    """Gets the day of a run: today, as the stations reporting are those of now."""
    del entry, fitting
    return [today()]


def read_ghcnd(directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads the GHCN-D station inventory, indexed by its identifier, with missing
    elevations as `NaN`.
    """
    text = cached_download(GHCND, directory / "ghcnd-stations.txt").decode("utf-8")
    table = pd.read_fwf(
        io.StringIO(text),
        colspecs=[(0, 11), (12, 20), (21, 30), (31, 37)],
        names=["ghcnd", "latitude", "longitude", "ground"],
    )
    table["ground"] = table["ground"].mask(np.isclose(table["ground"], GHCND_MISSING))
    return table.set_index("ghcnd")


def read_uscrn_listing(directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads NCEI's list of USCRN stations: the operational stations with a WBAN
    identifier.
    """
    listing = pd.read_csv(
        io.BytesIO(cached_download(USCRN_LIST, directory / "uscrn-stations.tsv")),
        sep="\t",
        dtype=str,
    )
    listing = listing[
        (listing.NETWORK == "USCRN") & (listing.OPERATION == "Operational")
    ]
    without = ~listing.WBAN.astype(str).str.isdigit()
    for station in listing[without].itertuples():
        print(f"  {station.LOCATION}: no WBAN identifier, not included")
    return listing[~without]


def read_homr_stations(system: str, directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads HOMR's table of ASOS or AWOS stations: those of the United States and
    its territories with known elevations, identified by their call signs (or
    WBAN or NCDC identifiers, if they have none).
    """
    table = read_homr_table(system, directory)
    table = table[table.COUNTRY.isin(UNITED_STATES)].copy()
    # elevations of -99999 ft are unknown
    unknown = table.ELEV.astype(float) < -1000
    print(
        f"  {unknown.sum()} of {len(table)} stations without an elevation, not included"
    )
    table = table[~unknown]
    table["id"] = table.CALL.fillna("WBAN" + table.WBAN.fillna(table.NCDCID))
    return table.drop_duplicates("id")


def list_station_ids(system: str, directory: pathlib.Path) -> set[str]:
    """Lists the identifiers of a network's stations, as its listing gives them."""
    if system == "uscrn":
        return set(read_uscrn_listing(directory).WBAN)
    return set(read_homr_stations(system, directory).id)


def homr_position(wban: str, directory: pathlib.Path) -> dict:
    """Gets a station's position and ground elevation from HOMR's service."""
    record = json.loads(
        cached_download(
            HOMR_SERVICE,
            directory / "homr" / f"{wban}.json",
            qid=f"WBAN:{wban}",
            qidMod="is",
        )
    )["stationCollection"]["stations"][0]
    location = record["location"]
    ground = [
        e for e in location.get("elevations", []) if e.get("elevationType") == "GROUND"
    ]
    return {
        "latitude": float(record["header"]["latitude_dec"]),
        "longitude": float(record["header"]["longitude_dec"]),
        "ground": float(ground[0]["elevationMeters"]) if ground else np.nan,
    }


def build_table(system: str, directory: pathlib.Path) -> pd.DataFrame:
    """Builds the table of a network's stations (see `--fit`)."""
    if system == "uscrn":
        rows = [
            {
                "id": station.WBAN,
                "name": f"{station.LOCATION} {station.VECTOR}, {station.STATE}",
                **homr_position(station.WBAN, directory),
            }
            for station in read_uscrn_listing(directory).itertuples()
        ]
        return pd.DataFrame(rows).sort_values("id")
    table = read_homr_stations(system, directory)
    return pd.DataFrame(
        {
            "id": table.id,
            "name": table.NAME.str.title(),
            "latitude": table.LAT.astype(float),
            "longitude": table.LON.astype(float),
            "ground": (table.ELEV.astype(float) * FEET).round(1),
            "wban": table.WBAN,
        }
    ).sort_values("id")


def get_reporting(system: str, directory: pathlib.Path) -> set[str]:
    """Gets the stations IEM reports online, as of the run."""
    if system == "uscrn":
        networks = ["USCRN"]
    else:
        listing = cached_download(IEM_NETWORKS, directory / "iem-networks.json")
        ids = [n["id"] for n in json.loads(listing)["data"]]
        networks = [n for n in ids if re.fullmatch(r"[A-Z]{2}_ASOS", n)]
    online = set()
    for network in networks:
        features = json.loads(
            cached_download(
                IEM_NETWORK.format(network), directory / "iem" / f"{network}.geojson"
            )
        )["features"]
        online |= {
            f["properties"]["sid"] for f in features if f["properties"]["online"]
        }
    if system != "uscrn":
        # IEM names the stations of Alaska, Hawaii, and the territories by their
        # ICAO identifiers (such as PANC), which HOMR names without the prefix
        online |= {sid[1:] for sid in online if len(sid) == 4 and sid[0] in "KPT"}
    return online


def get_ghcnd_ids(system: str, table: pd.DataFrame) -> pd.Series:
    """
    Gets the GHCN-D identifiers of a table's stations, from their WBAN
    identifiers (empty for a station without one).
    """
    wban = table.id if system == "uscrn" else table.wban
    wban = wban.fillna("").astype(str)
    return ("USW000" + wban.str.zfill(5)).where(wban.str.isdigit(), "")


def test(system: str, table: pd.DataFrame, directory: pathlib.Path) -> list[dict]:
    """Tests the siting and membership of a network's table."""
    ghcnd = read_ghcnd(directory)
    indexed = table.assign(ghcnd=get_ghcnd_ids(system, table).to_numpy()).set_index(
        "ghcnd"
    )
    indexed = indexed[indexed.index.isin(ghcnd.index)]
    horizontal = horizontal_errors(indexed, ghcnd)
    vertical = (indexed.ground - ghcnd.loc[indexed.index, "ground"]).dropna()
    print(f"  {len(indexed)} of {len(table)} stations in GHCN-D")
    # the stations of the network reporting, as its list (or HOMR's table)
    # gives them on the day of the run
    reporting = get_reporting(system, directory) & list_station_ids(system, directory)
    recall, precision = membership(set(table.id), reporting)
    print(f"  missing from the table: {sorted(reporting - set(table.id))[:20]}")
    print(f"  not reporting: {len(set(table.id) - reporting)} stations")
    return [
        p95_metric("horizontal_error_p95", horizontal),
        p95_metric("ground_elevation_error_p95", vertical),
        value_metric("membership_recall", recall),
        value_metric("membership_precision", precision),
    ]


def fit(run: SystemRun, directory: pathlib.Path) -> None:
    """Builds the table of a network's stations."""
    table = build_table(run.entry.id, directory)
    path = TABLES / run.entry.stations.file
    table.to_csv(path, index=False, lineterminator="\n")
    print(
        f"wrote {len(table)} stations to {path}; set the table's as_of and the "
        f"period of its provenance to {today()} in {run.entry.id}.yml"
    )


def main():
    run = parse_system_run(
        __doc__,
        system_classes=[SystemClass.SURFACE_NETWORK],
        default_days=default_days,
    )
    if run.days != [today()]:
        raise SystemExit("the stations reporting are those of today; give no days")
    directory = run.directory / f"{today()}"
    if run.fit:
        fit(run, directory)
        return
    metrics = test(run.entry.id, run.entry.stations.load(), directory)
    for metric in metrics:
        print(f"  {metric['name']}: {metric['value']}")
    results = {}
    for key in run.configurations:
        if run.write:
            run.warn_holdout(key)
        results[key] = system_result(
            run.entry,
            key,
            source=run.source,
            reference=REFERENCE,
            period=run.period,
            platforms=[run.entry.id],
            metrics=metrics,
        )
    if run.write:
        write_system_record(run.entry, results)


if __name__ == "__main__":
    main()
