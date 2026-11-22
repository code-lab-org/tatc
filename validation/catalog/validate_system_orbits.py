"""
Validates the orbits of a catalog satellite system against element sets: the
planes of the orbits (mean altitude, inclination, and local time of ascending
node), the longitudes of geostationary satellites, and the ascending nodes of
orbits maintained on repeat ground tracks.

Element sets come from the registry (`elements/element_sets.yml`), NOAA's archives
of the JPSS satellites (cached in the directory), or the files of element sets
beside it, each standing for the day of its own epoch. A run tests every
configuration of the system on each spacecraft's days at least
`MIN_HOLDOUT_GAP_DAYS` after the latest day to which its values were fitted, to
`--end` (or a retired spacecraft's last day): the days with element sets (and,
for spacecraft in NOAA's archives, every `--step` days), or the days given.
Records state each spacecraft's days (`platform_periods`).

With `--fit`, it prints the values fitted on the days of their provenance (or the
day given): the plane of an element set, the mean longitude of a geostationary
satellite, or the ascending node of a repeat ground track (from the element sets
of the days, without those after maneuvers).

Usage:
    python validate_system_orbits.py jpss DIRECTORY [--configuration plane] [--write]
    python validate_system_orbits.py jpss DIRECTORY 2026-08-24 2026-08-25 ... \
        --fit --configuration ground_track

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import date, datetime, timedelta, timezone

from tatc.catalog import SystemClass, get_spacecraft
from tatc.catalog.records import ConfigurationRecord, format_period, get_period_days
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS
from tatc.catalog.system_schemas import OrbitConfiguration
from tatc.schemas import RepeatGroundTrackOrbit
from validation_io import orbits
from validation_io.archives import ARCHIVES
from validation_io.cli import SystemRun, parse_day, parse_system_run
from validation_io.records import (
    system_result,
    value_metric,
    write_system_record,
)

REFERENCE = (
    "Element sets of the spacecraft (registered, from NOAA's archives of the JPSS "
    "satellites, or in files of element sets)"
)


def get_fit_days(
    configuration: OrbitConfiguration, spacecraft: str, own: bool = False
) -> list[str]:
    """
    Gets the days of the data to which a spacecraft's values were fitted (with
    `own`, only those fitted to its own data, not its values relative to
    another spacecraft's).
    """
    days = set()
    for provenance in configuration.provenance[spacecraft].values():
        if provenance.kind == "measured" and (
            not own or list(provenance.platforms) == [spacecraft]
        ):
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
    `MIN_HOLDOUT_GAP_DAYS` after its latest fitted day, to `end` (or the day of
    its status, if retired), every `step` days for spacecraft in NOAA's archives.
    """
    status = get_spacecraft(spacecraft)
    if status.status == "retired":
        end = min(end, status.as_of)
    fitted = get_fit_days(configuration, spacecraft)
    first = (
        date.fromisoformat(max(fitted)) + timedelta(days=MIN_HOLDOUT_GAP_DAYS)
        if fitted
        else None
    )
    days = {
        d
        for d in orbits.get_available_days(spacecraft)
        if (first is None or first.isoformat() <= d) and d <= end.isoformat()
    }
    if spacecraft in ARCHIVES and first is not None:
        day = first
        while day <= end:
            days.add(day.isoformat())
            day += timedelta(days=step)
    return sorted(days)


def load_element_sets(
    run: SystemRun, spacecraft: str, days: list[str]
) -> dict[str, list[str]]:
    """
    Loads the element sets of days, by the days of their epochs (see
    `orbits.get_element_sets`), reporting the days without a set of their own.
    """
    sets = orbits.get_element_sets(spacecraft, days, run.directory)
    missing = sorted(set(days) - set(sets))
    if missing:
        print(f"  {spacecraft}: no element set of its own on {', '.join(missing)}")
    return sets


def fit(run: SystemRun, key: str) -> None:
    """Prints the values of a configuration fitted on its days (or the run's)."""
    configuration = run.entry.configurations[key]
    given = [day.isoformat() for day in run.days]
    fitted = {}
    for spacecraft, orbit in configuration.orbits.items():
        days = given or get_fit_days(configuration, spacecraft, own=True)
        sets = load_element_sets(run, spacecraft, days)
        fitted[spacecraft] = sets
        tles, found = list(sets.values()), list(sets)
        if not tles:
            continue
        if isinstance(orbit, RepeatGroundTrackOrbit):
            kept = orbits.drop_maneuvers(tles)
            node = orbits.fit_ascending_node(kept, orbit.revolutions, orbit.days)
            moment = f"{node.time:%Y-%m-%dT%H:%M:%S.%f}"[:-3] + "Z"
            print(
                f'{spacecraft}: ascending_node: {{time: "{moment}", longitude: '
                f"{node.longitude}}} ({len(kept)} of {len(tles)} element sets, "
                f"{format_period(found)})"
            )
            continue
        if len(tles) > 1:
            print(f"{spacecraft}: give one day to fit a plane or a slot")
            continue
        orbit_of_set = orbits.to_orbit(tles[0])
        if orbit.type == "geosynchronous":
            slot = orbits.get_mean_slot(orbit_of_set)
            print(
                f"{spacecraft}: longitude: {slot['longitude']:.3f}, inclination: "
                f"{slot['inclination']:.2f} ({found[0]})"
            )
        else:
            plane = orbits.get_plane(orbit_of_set)
            print(
                f"{spacecraft}: mean_altitude: {round(plane['mean_altitude'], -1):.0f}, "
                f"inclination: {plane['inclination']:.4f}, equator_crossing_time: "
                f"\"{orbits.to_time_of_day(plane['local_time']):%H:%M:%S}\", "
                "reference_time (if it drifts): "
                f'"{orbit_of_set.get_epoch():%Y-%m-%dT%H:%M:%S}Z" ({found[0]})'
            )
    if configuration.get_relative_values():
        fit_relative(fitted)


def fit_relative(fitted: dict[str, dict[str, list[str]]]) -> None:
    """
    Prints the node and phase of each spacecraft relative to the first's, at
    12:00 UTC on the latest day of their element sets (one set each).
    """
    sets = {sc: list(by_day.values())[0] for sc, by_day in fitted.items() if by_day}
    if not sets:
        return
    latest = max(day for by_day in fitted.values() for day in by_day)
    moment = datetime.fromisoformat(latest).replace(hour=12, tzinfo=timezone.utc)
    first = next(iter(sets))
    phase_1, node_1 = orbits.get_relative_angles(sets[first], moment)
    print(f"relative to {first} at {moment:%Y-%m-%dT%H:%MZ}:")
    for spacecraft, tle in sets.items():
        phase, node = orbits.get_relative_angles(tle, moment)
        print(
            f"{spacecraft}: right_ascension_ascending_node: "
            f"{round((node - node_1) % 360, 2) % 360:g}, phase: "
            f"{round((phase - phase_1) % 360, 2) % 360:g}"
        )


def test(run: SystemRun, key: str) -> ConfigurationRecord | None:
    """
    Tests a configuration on the days of each spacecraft, returning its record
    (or `None`, if no spacecraft has days to test).
    """
    configuration = run.entry.configurations[key]
    given = [day.isoformat() for day in run.days]
    metrics, periods, all_sets = [], {}, {}
    for spacecraft, orbit in configuration.orbits.items():
        days = given or get_test_days(configuration, spacecraft, run.step, run.end)
        sets = load_element_sets(run, spacecraft, days)
        if not sets:
            print(f"  {spacecraft}: no days to test")
            continue
        metrics += orbits.compare(orbit, list(sets.values()), spacecraft)
        span = date.fromisoformat(max(sets)) - date.fromisoformat(min(sets))
        metrics += [
            value_metric("test_sets", len(sets), spacecraft),
            value_metric("test_span_days", span.days, spacecraft),
        ]
        periods[spacecraft] = format_period(sets)
        all_sets[spacecraft] = sets
    metrics += orbits.compare_relative(configuration, all_sets)
    for metric in metrics:
        print(f"  {metric['scope']} {metric['name']}: {metric['value']}")
    all_days = {day for sets in all_sets.values() for day in sets}
    if not metrics:
        return None
    period = format_period(all_days)
    if run.write:
        run.warn_holdout(key, period, periods)
    return system_result(
        run.entry,
        key,
        source=run.source,
        reference=REFERENCE,
        period=period,
        platforms=list(periods),
        metrics=metrics,
        platform_periods=periods,
    )


def _add_arguments(parser) -> None:
    yesterday = datetime.now(timezone.utc).date() - timedelta(days=1)
    parser.add_argument(
        "--step", type=int, default=7, help="days between days of NOAA's archives"
    )
    parser.add_argument(
        "--end",
        type=lambda text: parse_day(text).date(),
        default=yesterday,
        help="the last day to test (default: yesterday)",
    )


def main():
    run = parse_system_run(
        __doc__,
        system_classes=[c for c in SystemClass if c.is_satellite()],
        add_arguments=_add_arguments,
    )
    if run.fit:
        for key in run.configurations:
            print(f"{run.entry.id}:{key}")
            fit(run, key)
        return
    results = {}
    for key in run.configurations:
        print(f"{run.entry.id}:{key}")
        result = test(run, key)
        if result is not None:
            results[key] = result
    if run.write and results:
        write_system_record(run.entry, results)


if __name__ == "__main__":
    main()
