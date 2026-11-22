"""
Validates the catalog's NEXRAD network (system `nexrad`) against the radars' own
Level II data on AWS (no account needed) and NCEI's table of NEXRAD stations.

* Siting: the table's positions and antenna altitudes (ground plus feedhorn)
  compared with NCEI's table, as downloaded on the day of the run.
* Membership: the WSR-88D radars of NCEI's table that recorded Level II volumes
  on the test's days, compared with the table's stations.
* Coverage: for each configuration, the elevation angles scanned by sampled
  volumes of its volume coverage patterns (`volume_coverage_patterns`), and the
  range of their reflectivity, compared with its lowest (or a station's own)
  elevation angle and maximum range, and the angles scanned above its highest
  elevation angle (volumes ended below it observe less, as its limitations
  state).

The test's days default to the week from `MIN_HOLDOUT_GAP_DAYS` after the
table's snapshot. With `--fit` and a day, it builds the table instead
(`resources/systems/tables/nexrad.csv`): every station of NCEI's table that
recorded a volume on the day, with the position, ground, and feedhorn height of
its first volume of the day, and its own lowest elevation angle if its pattern
scans below the configurations' lowest cut.

Usage:
    python validate_system_nexrad.py nexrad DIRECTORY 2026-10-15 --fit
    python validate_system_nexrad.py nexrad DIRECTORY [DAY ...] [--write]

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import random
from datetime import date, datetime, timedelta, timezone

import numpy as np
import pandas as pd

from tatc.catalog import SystemClass, SystemEntry, SystemGrade
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS
from tatc.catalog.system_schemas import RadarConfiguration
from validation_io import stations
from validation_io.cli import SystemRun, parse_system_run
from validation_io.records import (
    REPOSITORY,
    p95_metric,
    system_result,
    value_metric,
    write_system_record,
)

TABLES = REPOSITORY / "src" / "tatc" / "resources" / "systems" / "tables"
# the days of a test, from MIN_HOLDOUT_GAP_DAYS after the table's snapshot
TEST_DAYS = 7
# the rounding (degrees) of commanded elevation angles, below the configurations'
# lowest cut by more than which a station's lowest cut is its own
CUT_TOLERANCE = 0.05
REFERENCE = (
    "NEXRAD Level II volumes (unidata-nexrad-level2 on AWS) and NCEI's table of "
    "NEXRAD stations"
)


def default_days(entry: SystemEntry, fitting: bool) -> list[date]:
    """
    Gets the days of a test: the `TEST_DAYS` from `MIN_HOLDOUT_GAP_DAYS` after
    the table's snapshot.
    """
    if fitting:
        raise ValueError("give the day of the volumes from which to build the table")
    first = entry.stations.as_of + timedelta(days=MIN_HOLDOUT_GAP_DAYS)
    return [first + timedelta(days=k) for k in range(TEST_DAYS)]


def read_ncei(directory: pathlib.Path) -> pd.DataFrame:
    """Reads NCEI's table of NEXRAD stations, as downloaded today (UTC)."""
    today = datetime.now(timezone.utc).date()
    return stations.read_ncei_nexrad(directory / "homr" / f"{today}")


def build_table(entry: SystemEntry, day: date, directory) -> pd.DataFrame:
    """Builds the table of stations from the volumes of a day (see `--fit`)."""
    own_below = (
        min(c.min_elevation_angle for c in entry.configurations.values())
        - CUT_TOLERANCE
    )
    ncei = read_ncei(directory)
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
                "min_elevation_angle": lowest if lowest < own_below else np.nan,
            }
        )
        print(f"  {icao}: {rows[-1]}")
    return pd.DataFrame(rows)


def test_siting(table: pd.DataFrame, days: list[date], directory) -> list[dict]:
    """Tests the positions, antenna altitudes, and membership of the table."""
    ncei = read_ncei(directory)
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
        p95_metric("horizontal_error_p95", horizontal),
        p95_metric("antenna_altitude_error_p95", vertical),
        value_metric("membership_recall", recall),
        value_metric("membership_precision", precision),
    ]


def _sample_volume(
    configuration: RadarConfiguration,
    generator: random.Random,
    stations_ids: list[str],
    days: list[date],
    directory,
) -> tuple[str, str, dict] | None:
    """
    Samples a whole volume of a configuration's patterns, at a random station
    and time, or `None` if the sampled volume is not of its patterns or lacks
    radials or reflectivity.
    """
    icao = generator.choice(stations_ids)
    keys = stations.list_level_2_volumes(generator.choice(days), icao)
    if not keys:
        return None
    chosen = generator.choice(keys)
    try:
        head = stations.describe_level_2(
            stations.fetch_level_2(chosen, directory, records=6)
        )
        if head["vcp"] not in configuration.volume_coverage_patterns:
            return None
        volume = stations.describe_level_2(
            stations.fetch_level_2(chosen, directory, records=None)
        )
    except (ValueError, IndexError, KeyError):
        return None
    if "lowest" not in volume or "range" not in volume:
        return None
    return icao, chosen, volume


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
    random stations (with a fixed seed) at random times of the days, or returns
    no metrics if no volume of its patterns is found.
    """
    if not configuration.volume_coverage_patterns:
        raise ValueError(f"{key} has no volume coverage patterns to compare with")
    generator = random.Random(f"{key}-{days[0]}")
    indexed = table.set_index("id")
    stations_ids = sorted(indexed.index)
    lowest, highest, ranges, tried = [], [], [], 0
    while len(lowest) < samples and tried < 20 * samples:
        tried += 1
        sampled = _sample_volume(
            configuration, generator, stations_ids, days, directory
        )
        if sampled is None:
            continue
        icao, chosen, volume = sampled
        own = indexed.loc[icao, "min_elevation_angle"]
        expected = configuration.min_elevation_angle if pd.isna(own) else own
        lowest.append(volume["lowest"] - expected)
        # one-sided: a volume ended below the highest cut (as AVSET ends them)
        # observes less, as the configuration's limitations state
        highest.append(max(0.0, volume["highest"] - configuration.max_elevation_angle))
        ranges.append(
            100 * (volume["range"] - configuration.max_range) / configuration.max_range
        )
        print(
            f"  {chosen}: VCP {volume['vcp']}, lowest {volume['lowest']:.2f}, "
            f"highest {volume['highest']:.2f}, range {volume['range'] / 1e3:.0f} km"
        )
    if len(lowest) < samples:
        print(f"  warning: {len(lowest)} of {samples} volumes in {tried} tries")
    if not lowest:
        return []
    return [
        p95_metric("min_elevation_error_p95", lowest),
        p95_metric("max_elevation_excess_p95", highest),
        p95_metric("range_error_p95", ranges),
    ]


def _add_arguments(parser) -> None:
    parser.add_argument(
        "--samples", type=int, default=24, help="volumes per configuration"
    )


def fit(run: SystemRun) -> None:
    """Builds the table of stations from the volumes of the run's day."""
    if len(run.days) != 1:
        raise SystemExit("give one day of volumes from which to build the table")
    table = build_table(run.entry, run.days[0], run.directory)
    path = TABLES / run.entry.stations.file
    table.to_csv(path, index=False, lineterminator="\n")
    print(
        f"wrote {len(table)} stations to {path}; set the table's as_of and the "
        f"period of its provenance to {run.days[0]} in {run.entry.id}.yml"
    )


def main():
    run = parse_system_run(
        __doc__,
        system_classes=[SystemClass.RADAR_NETWORK],
        default_days=default_days,
        add_arguments=_add_arguments,
    )
    if run.fit:
        fit(run)
        return
    table = run.entry.stations.load()
    print("siting")
    siting = test_siting(table, run.days, run.directory)
    for metric in siting:
        print(f"  {metric['name']}: {metric['value']}")
    results = {}
    for key in run.configurations:
        configuration = run.entry.configurations[key]
        print(key)
        metrics = list(siting)
        if configuration.grade == SystemGrade.COVERAGE:
            coverage = test_coverage(
                key, configuration, table, run.days, run.samples, run.directory
            )
            for metric in coverage:
                print(f"  {metric['name']}: {metric['value']}")
            metrics += coverage
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
