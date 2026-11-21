"""
Validates the catalog's NEXRAD network (system `nexrad`) against the radars' own
Level II data on AWS (no account needed) and NCEI's table of NEXRAD stations.

* Siting: the table's positions and antenna altitudes (ground plus feedhorn)
  compared with NCEI's table, as downloaded on the test.
* Membership: the WSR-88D radars of NCEI's table that recorded Level II volumes
  on the test's days, compared with the table's stations.
* Coverage: for each configuration, the elevation angles scanned by sampled
  volumes of its volume coverage patterns (VCPs), and the range of their
  reflectivity, compared with its lowest (or a station's own) and highest
  elevation angles and maximum range.

With `--fit DAY`, it builds the table instead (`resources/systems/tables/
nexrad.csv`): every station of NCEI's table that recorded a volume on the day,
with the position, ground, and feedhorn height of its first volume of the day,
and its own lowest elevation angle if its pattern scans below 0.45 degrees.

Usage:
    python validate_system_nexrad.py --fit 2026-10-15
    python validate_system_nexrad.py [--days 2026-11-14 ... ] [--write]

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
import random
from datetime import date, timedelta

import numpy as np
import pandas as pd

import tatc
from tatc.catalog import SYSTEM_STANDARD_VERSION, get_system
from tatc.catalog.system_schemas import RadarConfiguration
from validation_io import stations
from validation_io.records import REPOSITORY, source_of
from validate_system_orbits import write_record

TABLE = REPOSITORY / "src" / "tatc" / "resources" / "systems" / "tables" / "nexrad.csv"
# the volume coverage patterns of each configuration
PATTERNS = {"precipitation": (12, 112, 212, 215), "clear_air": (31, 32, 35)}
# the elevation angle (degrees) below which a station's pattern is its own
OWN_LOWEST = 0.45
REFERENCE = (
    "NEXRAD Level II volumes (unidata-nexrad-level2 on AWS) and NCEI's table of "
    "NEXRAD stations"
)


def build_table(day: date, directory: pathlib.Path) -> pd.DataFrame:
    """Builds the table of stations from the volumes of a day (see `--fit`)."""
    ncei = stations.read_ncei_nexrad(directory)
    reporting = set(stations.list_level_2_stations(day))
    rows = []
    for icao in sorted(set(ncei.index) & reporting):
        volume = None
        for key in stations.list_level_2_volumes(day, icao)[:3]:
            try:
                volume = stations.describe_level_2(
                    stations.fetch_level_2(key, directory, records=6)
                )
                break
            except (ValueError, IndexError, KeyError) as error:
                print(f"  {key}: {error}")
        if volume is None:
            continue
        lowest = min(volume["commanded"])
        rows.append(
            {
                "id": icao,
                "name": str(ncei.loc[icao, "NAME"]).title(),
                "latitude": round(volume["latitude"], 6),
                "longitude": round(volume["longitude"], 6),
                "ground": volume["ground"],
                "height": volume["height"],
                "min_elevation_angle": lowest if lowest < OWN_LOWEST else np.nan,
            }
        )
        print(f"  {icao}: {rows[-1]}")
    return pd.DataFrame(rows)


def test_siting(table: pd.DataFrame, days: list[date], directory) -> list[dict]:
    """Tests the positions, antenna altitudes, and membership of the table."""
    ncei = stations.read_ncei_nexrad(directory)
    indexed = table.set_index("id")
    reference = ncei.rename(columns={"LAT": "latitude", "LON": "longitude"})
    horizontal = stations.horizontal_errors(indexed, reference)
    common = indexed.index.intersection(ncei.index)
    vertical = (
        indexed.loc[common, "ground"] + indexed.loc[common, "height"]
    ) - ncei.loc[common, "antenna"]
    reporting = set()
    for day in days:
        reporting |= set(stations.list_level_2_stations(day))
    reporting &= set(ncei.index)
    recall, precision = stations.membership(set(indexed.index), reporting)
    print(f"  missing from the table: {sorted(reporting - set(indexed.index))}")
    print(f"  not reporting: {sorted(set(indexed.index) - reporting)}")
    return [
        _metric("horizontal_error_p95", np.percentile(np.abs(horizontal), 95)),
        _metric("antenna_altitude_error_p95", np.percentile(np.abs(vertical), 95)),
        _metric("membership_recall", recall),
        _metric("membership_precision", precision),
    ]


def test_coverage(
    key: str,
    configuration: RadarConfiguration,
    table: pd.DataFrame,
    days: list[date],
    samples: int,
    directory,
) -> list[dict]:
    """
    Tests a configuration against whole volumes of its patterns, sampled from
    random stations (with a fixed seed) at random times of the days.
    """
    generator = random.Random(f"{key}-{days[0]}")
    indexed = table.set_index("id")
    lowest, highest, ranges, tried = [], [], [], 0
    while len(lowest) < samples and tried < 20 * samples:
        tried += 1
        icao = generator.choice(sorted(indexed.index))
        day = generator.choice(days)
        keys = stations.list_level_2_volumes(day, icao)
        if not keys:
            continue
        chosen = generator.choice(keys)
        try:
            head = stations.describe_level_2(
                stations.fetch_level_2(chosen, directory, records=6)
            )
        except (ValueError, IndexError, KeyError):
            continue
        if head["vcp"] not in PATTERNS[key]:
            continue
        volume = stations.describe_level_2(
            stations.fetch_level_2(chosen, directory, records=None)
        )
        own = indexed.loc[icao, "min_elevation_angle"]
        expected = configuration.min_elevation_angle if pd.isna(own) else own
        lowest.append(volume["lowest"] - expected)
        highest.append(volume["highest"] - configuration.max_elevation_angle)
        ranges.append(
            100 * (volume["range"] - configuration.max_range) / configuration.max_range
        )
        print(
            f"  {chosen}: VCP {volume['vcp']}, lowest {volume['lowest']:.2f}, "
            f"highest {volume['highest']:.2f}, range {volume['range'] / 1e3:.0f} km"
        )
    return [
        _metric("min_elevation_error_p95", np.percentile(np.abs(lowest), 95)),
        _metric("max_elevation_error_p95", np.percentile(np.abs(highest), 95)),
        _metric("range_error_p95", np.percentile(np.abs(ranges), 95)),
    ]


def _metric(name: str, value: float) -> dict:
    return {"name": name, "value": round(float(value), 4)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--fit", help="Build the table from a day's volumes.")
    parser.add_argument(
        "--days",
        nargs="*",
        default=[
            (date(2026, 11, 14) + timedelta(days=k)).isoformat() for k in range(7)
        ],
        help="The days (ISO 8601) of the test.",
    )
    parser.add_argument("--samples", type=int, default=24, help="Volumes per pattern.")
    parser.add_argument("--configuration", nargs="*", help="The configurations.")
    parser.add_argument("--write", action="store_true", help="Write the record.")
    parser.add_argument(
        "--directory",
        default=str(REPOSITORY / "validation" / "data" / "nexrad"),
        help="Where to cache downloads.",
    )
    args = parser.parse_args()
    directory = pathlib.Path(args.directory)
    if args.fit:
        table = build_table(date.fromisoformat(args.fit), directory)
        table.to_csv(TABLE, index=False, lineterminator="\n")
        print(f"wrote {len(table)} stations to {TABLE}")
        return
    entry = get_system("nexrad")
    table = entry.stations.load()
    days = [date.fromisoformat(day) for day in args.days]
    print("siting")
    siting = test_siting(table, days, directory)
    for metric in siting:
        print(f"  {metric['name']}: {metric['value']}")
    results = {}
    for key in args.configuration or list(entry.configurations):
        configuration = entry.configurations[key]
        print(key)
        metrics = list(siting)
        if configuration.grade.value == "coverage":
            metrics += test_coverage(
                key, configuration, table, days, args.samples, directory
            )
        for metric in metrics[len(siting) :]:
            print(f"  {metric['name']}: {metric['value']}")
        results[key] = {
            "digest": configuration.digest,
            "standard_version": SYSTEM_STANDARD_VERSION,
            "tatc_version": tatc.__version__,
            "source": source_of(__file__),
            "reference": REFERENCE,
            "period": f"{days[0]} to {days[-1]}",
            "platforms": [entry.id],
            "metrics": metrics,
        }
    if args.write:
        print(f"wrote {write_record(entry, results)}")


if __name__ == "__main__":
    main()
