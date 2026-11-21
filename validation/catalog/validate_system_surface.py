"""
Validates the catalog's surface networks (systems `uscrn`, `asos`, and `awos`)
against NCEI's station metadata and the stations reporting observations.

* Siting: the table's positions and ground elevations compared with NCEI's
  inventory of the Global Historical Climatology Network - Daily (GHCN-D), a
  product separate from the table's source, though compiled from the same
  station histories (NCEI's Historical Observing Metadata Repository, HOMR).
* Membership: the stations that the Iowa Environmental Mesonet (IEM) reports
  online (observations received recently: USCRN, or the METAR stations of the
  United States' networks) among the network's stations as listed on the test
  (as the table is built), compared with the table's stations.

With `--fit`, it builds the table instead (`resources/systems/tables/<system>.csv`)
from HOMR: for USCRN, the operational stations of NCEI's USCRN list (with a WBAN
identifier), at their HOMR positions and ground elevations; for ASOS and AWOS,
the stations of the United States and its territories in HOMR's tables, with
known elevations. Elevations are meters above mean sea level.

Usage:
    python validate_system_surface.py uscrn --fit
    python validate_system_surface.py uscrn [--write]

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import io
import json
import pathlib
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import tatc
from tatc.catalog import SYSTEM_STANDARD_VERSION, get_system
from validation_io import get_with_retries
from validation_io.records import REPOSITORY, source_of
from validation_io.stations import FEET, horizontal_errors, membership
from validate_system_orbits import write_record

TABLES = REPOSITORY / "src" / "tatc" / "resources" / "systems" / "tables"
USCRN_LIST = "https://www.ncei.noaa.gov/pub/data/uscrn/products/stations.tsv"
HOMR_FILE = "https://www.ncei.noaa.gov/access/homr/file/{}-stations.txt"
HOMR_SERVICE = "https://www.ncei.noaa.gov/access/homr/services/station/search"
GHCND = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
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


def _cached(url: str, directory: pathlib.Path, name: str, **params) -> bytes:
    """Downloads a file once, into a directory."""
    path = directory / name
    if not path.exists():
        directory.mkdir(parents=True, exist_ok=True)
        path.write_bytes(get_with_retries(url, params=params or None).content)
    return path.read_bytes()


def read_homr_table(kind: str, directory: pathlib.Path) -> pd.DataFrame:
    """Reads HOMR's fixed-width table of ASOS or AWOS stations."""
    text = _cached(HOMR_FILE.format(kind), directory, f"{kind}-stations.txt")
    lines = [line for line in text.decode("utf-8").splitlines() if line.strip()]
    spans = [(m.start(), m.end()) for m in re.finditer(r"-+", lines[1])]
    return pd.read_fwf(
        io.StringIO("\n".join(lines)), colspecs=spans, skiprows=[1], dtype=str
    )


def read_ghcnd(directory: pathlib.Path) -> pd.DataFrame:
    """Reads the GHCN-D station inventory, indexed by its identifier."""
    text = _cached(GHCND, directory, "ghcnd-stations.txt").decode("utf-8")
    table = pd.read_fwf(
        io.StringIO(text),
        colspecs=[(0, 11), (12, 20), (21, 30), (31, 37)],
        names=["ghcnd", "latitude", "longitude", "ground"],
    )
    return table.set_index("ghcnd")


def homr_position(wban: str, directory: pathlib.Path) -> dict:
    """Gets a station's position and ground elevation from HOMR's service."""
    record = json.loads(
        _cached(
            HOMR_SERVICE,
            directory / "homr",
            f"{wban}.json",
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
        listing = pd.read_csv(
            io.BytesIO(_cached(USCRN_LIST, directory, "uscrn-stations.tsv")),
            sep="\t",
            dtype=str,
        )
        listing = listing[
            (listing.NETWORK == "USCRN") & (listing.OPERATION == "Operational")
        ]
        rows = []
        for station in listing.itertuples():
            if not str(station.WBAN).isdigit():
                print(f"  {station.LOCATION}: no WBAN identifier, not included")
                continue
            position = homr_position(station.WBAN, directory)
            rows.append(
                {
                    "id": station.WBAN,
                    "name": f"{station.LOCATION} {station.VECTOR}, {station.STATE}",
                    **position,
                }
            )
        return pd.DataFrame(rows).sort_values("id")
    table = read_homr_table(system, directory)
    table = table[table.COUNTRY.isin(UNITED_STATES)].copy()
    # elevations of -99999 ft are unknown
    unknown = table.ELEV.astype(float) < -1000
    print(
        f"  {unknown.sum()} of {len(table)} stations without an elevation, not included"
    )
    table = table[~unknown]
    table["id"] = table.CALL.fillna("WBAN" + table.WBAN.fillna(table.NCDCID))
    table = table.drop_duplicates("id")
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
    """Gets the stations IEM reports online, as of the test."""
    if system == "uscrn":
        networks = ["USCRN"]
    else:
        ids = [
            n["id"]
            for n in json.loads(_cached(IEM_NETWORKS, directory, "iem-networks.json"))[
                "data"
            ]
        ]
        networks = [n for n in ids if re.fullmatch(r"[A-Z]{2}_ASOS", n)]
    online = set()
    for network in networks:
        features = json.loads(
            _cached(
                IEM_NETWORK.format(network), directory / "iem", f"{network}.geojson"
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


def test(system: str, table: pd.DataFrame, directory: pathlib.Path) -> list[dict]:
    """Tests the siting and membership of a network's table."""
    ghcnd = read_ghcnd(directory)
    if system == "uscrn":
        keys = "USW000" + table.id.str.zfill(5)
    else:
        keys = "USW000" + table.wban.fillna("").astype(str).str.zfill(5)
    indexed = table.assign(ghcnd=keys.to_numpy()).set_index("ghcnd")
    indexed = indexed[indexed.index.isin(ghcnd.index)]
    horizontal = horizontal_errors(indexed, ghcnd)
    vertical = (indexed.ground - ghcnd.loc[indexed.index, "ground"]).dropna()
    print(f"  {len(indexed)} of {len(table)} stations in GHCN-D")
    # the stations of the network, as its list (or HOMR's table) gives them on
    # the test
    reporting = get_reporting(system, directory) & set(
        build_table(system, directory / "test").id
    )
    recall, precision = membership(set(table.id), reporting)
    print(f"  missing from the table: {sorted(reporting - set(table.id))[:20]}")
    print(f"  not reporting: {len(set(table.id) - reporting)} stations")
    return [
        _metric("horizontal_error_p95", np.percentile(np.abs(horizontal), 95)),
        _metric("ground_elevation_error_p95", np.percentile(np.abs(vertical), 95)),
        _metric("membership_recall", recall),
        _metric("membership_precision", precision),
    ]


def _metric(name: str, value: float) -> dict:
    return {"name": name, "value": round(float(value), 4)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("system", choices=["uscrn", "asos", "awos"])
    parser.add_argument("--fit", action="store_true", help="Build the table.")
    parser.add_argument("--write", action="store_true", help="Write the record.")
    parser.add_argument(
        "--directory",
        default=str(REPOSITORY / "validation" / "data" / "surface"),
        help="Where to cache downloads (delete it to download them again).",
    )
    args = parser.parse_args()
    today = datetime.now(timezone.utc).date()
    directory = pathlib.Path(args.directory) / f"{today}"
    if args.fit:
        table = build_table(args.system, directory)
        path = TABLES / f"{args.system}.csv"
        table.to_csv(path, index=False, lineterminator="\n")
        print(f"wrote {len(table)} stations to {path}")
        return
    entry = get_system(args.system)
    metrics = test(args.system, entry.stations.load(), directory)
    for metric in metrics:
        print(f"  {metric['name']}: {metric['value']}")
    results = {
        key: {
            "digest": configuration.digest,
            "standard_version": SYSTEM_STANDARD_VERSION,
            "tatc_version": tatc.__version__,
            "source": source_of(__file__),
            "reference": REFERENCE,
            "period": f"{today}",
            "platforms": [entry.id],
            "metrics": metrics,
        }
        for key, configuration in entry.configurations.items()
    }
    if args.write:
        print(f"wrote {write_record(entry, results)}")


if __name__ == "__main__":
    main()
