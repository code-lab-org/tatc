"""
Validates the orbits of a catalog satellite system against element sets: the
planes of the orbits (mean altitude, inclination, and local time of ascending
node), the longitudes of geostationary satellites, and the ascending nodes of
orbits maintained on repeat ground tracks.

Element sets come from the registry (`elements/element_sets.yml`), NOAA's archives
of the JPSS satellites, or the files of element sets beside it. A run tests every
configuration of the system on each spacecraft's days at least
`MIN_HOLDOUT_GAP_DAYS` after the latest day to which its values were fitted: the
days with element sets (and, for spacecraft in NOAA's archives, every `--step`
days to the end), or the days given. Records state each spacecraft's days
(`platform_periods`).

With `--fit`, it prints the values fitted on the days of their provenance (or the
days given): the plane of an element set, the mean longitude of a geostationary
satellite, or the ascending node of a repeat ground track (from the element sets
of the days, without those after maneuvers).

Usage:
    python validate_system_orbits.py jpss [--configuration plane] [--write]
    python validate_system_orbits.py jpss --fit --configuration ground_track \
        --days 2026-08-01 2026-08-11 ...

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import json
import pathlib
from datetime import date, datetime, timedelta, timezone

import tatc
from tatc.catalog import SYSTEM_STANDARD_VERSION, get_system
from tatc.catalog.records import ValidationRecord, get_period_days
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS
from tatc.catalog.system_schemas import (
    GeostationarySlot,
    OrbitConfiguration,
    SystemEntry,
)
from tatc.schemas import RepeatGroundTrackOrbit
from validation_io import orbits
from validation_io.records import REPOSITORY, source_of

SYSTEM_RECORDS = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "systems"
REFERENCE = (
    "Element sets of the spacecraft (registered, from NOAA's archives of the JPSS "
    "satellites, or in files of element sets)"
)


def get_fit_days(configuration: OrbitConfiguration, spacecraft: str) -> list[str]:
    """Gets the days of the data to which a spacecraft's values were fitted."""
    days = set()
    for provenance in configuration.provenance[spacecraft].values():
        if provenance.kind == "measured":
            for start, end in get_period_days(provenance.period):
                days |= {
                    (start + timedelta(days=d)).isoformat()
                    for d in range((end - start).days + 1)
                }
    return sorted(days)


def get_test_days(
    configuration: OrbitConfiguration, spacecraft: str, step: int, end: date
) -> list[str]:
    """
    Gets the days on which to test a spacecraft: those with element sets at least
    `MIN_HOLDOUT_GAP_DAYS` after its latest fitted day, every `step` days to `end`
    for spacecraft in NOAA's archives.
    """
    fitted = get_fit_days(configuration, spacecraft)
    if fitted:
        first = date.fromisoformat(max(fitted)) + timedelta(days=MIN_HOLDOUT_GAP_DAYS)
    else:
        first = date.min
    days = {d for d in orbits.get_available_days(spacecraft) if d >= first.isoformat()}
    if orbits.has_archive(spacecraft) and first != date.min:
        day = first
        while day <= end:
            days.add(day.isoformat())
            day += timedelta(days=step)
    return sorted(days)


def format_period(days: list[str]) -> str:
    """Formats days as a period, with runs of consecutive days as ranges."""
    runs: list[list[date]] = []
    for day in sorted(date.fromisoformat(d) for d in set(days)):
        if runs and day - runs[-1][-1] == timedelta(days=1):
            runs[-1].append(day)
        else:
            runs.append([day])
    return ", ".join(
        f"{run[0]}" if len(run) == 1 else f"{run[0]} to {run[-1]}" for run in runs
    )


def load_element_sets(spacecraft: str, days: list[str], directory) -> tuple:
    """Loads the element sets of days, skipping days without one."""
    tles, found = [], []
    for day in days:
        try:
            tles += orbits.get_element_sets(spacecraft, [day], directory)
            found.append(day)
        except FileNotFoundError:
            print(f"  {spacecraft}: no element set on {day}")
    return tles, found


def fit(entry: SystemEntry, key: str, days: list[str] | None, directory) -> None:
    """Prints the values of a configuration fitted on its days (or those given)."""
    configuration = entry.configurations[key]
    for spacecraft, orbit in configuration.orbits.items():
        fit_days = days or get_fit_days(configuration, spacecraft)
        tles, found = load_element_sets(spacecraft, fit_days, directory)
        if not tles:
            continue
        if isinstance(orbit, RepeatGroundTrackOrbit):
            kept = orbits.drop_maneuvers(tles)
            node = orbits.fit_ascending_node(kept, orbit.revolutions, orbit.days)
            moment = f"{node.time:%Y-%m-%dT%H:%M:%S.%f}"[:-3] + "Z"
            print(
                f'{spacecraft}: ascending_node: {{time: "{moment}", longitude: '
                f"{node.longitude}}} ({len(kept)} of {len(tles)} element sets, "
                f"{', '.join(found)})"
            )
        elif isinstance(orbit, GeostationarySlot):
            slot = orbits.get_slot(tles[0])
            print(
                f"{spacecraft}: longitude: {slot['longitude']:.3f}, inclination: "
                f"{slot['inclination']:.2f} ({found[0]})"
            )
        else:
            plane = orbits.get_plane(tles[0])
            print(
                f"{spacecraft}: mean_altitude: {round(plane['mean_altitude'], -1):.0f}, "
                f"inclination: {plane['inclination']:.4f}, equator_crossing_time: "
                f"\"{orbits.to_time_of_day(plane['local_time']):%H:%M:%S}\" "
                f"({found[0]})"
            )


def test(
    entry: SystemEntry,
    key: str,
    days: list[str] | None,
    step: int,
    end: date,
    directory,
) -> dict | None:
    """
    Tests a configuration on the days of each spacecraft, returning its record
    (or `None`, if no spacecraft has days to test).
    """
    configuration = entry.configurations[key]
    metrics, tested, all_days, periods = [], [], set(), {}
    for spacecraft, orbit in configuration.orbits.items():
        test_days = days or get_test_days(configuration, spacecraft, step, end)
        tles, found = load_element_sets(spacecraft, test_days, directory)
        if not tles:
            print(f"  {spacecraft}: no days to test")
            continue
        if isinstance(orbit, GeostationarySlot):
            results = orbits.compare_slot(orbit, tles, spacecraft)
        else:
            results = orbits.compare_plane(orbit, tles, spacecraft)
            if isinstance(orbit, RepeatGroundTrackOrbit):
                delta_time, delta_longitude = orbits.compare_nodes(
                    orbit, orbits.drop_maneuvers(tles)
                )
                results += [
                    orbits._metric(
                        "node_longitude_error_p95", delta_longitude, spacecraft
                    ),
                    orbits._metric("node_time_error_p95", delta_time, spacecraft),
                ]
        for result in results:
            print(f"  {spacecraft} {result['name']}: {result['value']}")
        metrics += results
        tested.append(spacecraft)
        all_days |= set(found)
        periods[spacecraft] = format_period(found)
    if not metrics:
        return None
    return {
        "digest": configuration.digest,
        "standard_version": SYSTEM_STANDARD_VERSION,
        "tatc_version": tatc.__version__,
        "source": source_of(__file__),
        "reference": REFERENCE,
        "period": format_period(sorted(all_days)),
        "platforms": tested,
        "platform_periods": periods,
        "metrics": metrics,
    }


def write_record(entry: SystemEntry, results: dict[str, dict]) -> pathlib.Path:
    """
    Writes the records of configurations to the system's record, keeping those of
    its other configurations, after checking it against the schemas.
    """
    path = SYSTEM_RECORDS / f"{entry.id}.json"
    if path.exists():
        record = json.loads(path.read_text(encoding="utf-8"))
    else:
        record = {"entry_id": entry.id, "configurations": {}}
    record["configurations"].update(results)
    record["configurations"] = dict(sorted(record["configurations"].items()))
    ValidationRecord(**record)
    SystemEntry(**{**entry.model_dump(mode="json"), "record": record})
    SYSTEM_RECORDS.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("system", help="The identifier of the system.")
    parser.add_argument("--configuration", nargs="*", help="The configurations.")
    parser.add_argument("--days", nargs="*", help="The days (ISO 8601) of the data.")
    parser.add_argument("--fit", action="store_true", help="Print fitted values.")
    parser.add_argument("--write", action="store_true", help="Write the record.")
    parser.add_argument(
        "--step", type=int, default=7, help="Days between days of NOAA's archives."
    )
    parser.add_argument(
        "--end",
        default=datetime.now(timezone.utc).date().isoformat(),
        help="The last day of NOAA's archives to test.",
    )
    parser.add_argument(
        "--directory",
        default=str(REPOSITORY / "validation" / "data" / "elements"),
        help="Where to cache element sets downloaded from NOAA's archives.",
    )
    args = parser.parse_args()
    entry = get_system(args.system)
    if not entry.system_class.is_satellite():
        raise SystemExit(f"{entry.id} is not a satellite system")
    directory = pathlib.Path(args.directory)
    directory.mkdir(parents=True, exist_ok=True)
    keys = args.configuration or list(entry.configurations)
    if args.fit:
        for key in keys:
            print(f"{entry.id}:{key}")
            fit(entry, key, args.days, directory)
        return
    results = {}
    for key in keys:
        print(f"{entry.id}:{key}")
        result = test(
            entry, key, args.days, args.step, date.fromisoformat(args.end), directory
        )
        if result is not None:
            results[key] = result
    if args.write and results:
        print(f"wrote {write_record(entry, results)}")


if __name__ == "__main__":
    main()
