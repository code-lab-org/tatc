"""
The orbits of the catalog's satellite systems compared with element sets: the
mean elements of the planes of the orbits, the longitudes of geostationary
satellites, and the ascending nodes of orbits maintained on repeat ground tracks,
each computed with public TAT-C methods.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import functools
import pathlib
import warnings
from datetime import datetime, time, timedelta, timezone

import numpy as np
import pandas as pd

from tatc.catalog import get_spacecraft
from tatc.catalog.system_schemas import GeostationarySlot
from tatc.constants import EARTH_EQUATORIAL_RADIUS, timescale
from tatc.schemas import GeneralPerturbationsOrbit, RepeatGroundTrackOrbit
from tatc.schemas.orbit import AscendingNode

from .archives import (
    ELEMENT_SET_FOLDER,
    _parse_element_sets,
    _read_registry,
    get_element_set,
)
from .records import p95_metric

# the time from an element set's epoch within which its nodes are compared, as
# errors of element sets grow with propagation
NODE_WINDOW = timedelta(hours=6)
# a jump of the semimajor axis (meters) between consecutive element sets above
# which a maneuver is presumed, after which sets are dropped for MANEUVER_SETTLING
MANEUVER_JUMP = 50
MANEUVER_SETTLING = timedelta(days=2)
# the length (km) of a degree of longitude along the equator
KM_PER_DEGREE = np.pi / 180 * EARTH_EQUATORIAL_RADIUS / 1e3
# the iterations within which the fit of an ascending node converges
NODE_FIT_ITERATIONS = 5


def to_orbit(tle: list[str]) -> GeneralPerturbationsOrbit:
    """Builds the orbit of an element set, propagated with its own drag."""
    return GeneralPerturbationsOrbit.from_tle(tle)


def get_epoch(tle: list[str]) -> datetime:
    """Gets the epoch of an element set."""
    return to_orbit(tle).get_epoch()


def get_local_time(orbit: GeneralPerturbationsOrbit) -> float:
    """
    Gets the local mean solar time (hours) of the ascending node of an orbit's
    first element at its epoch: the universal time (UT1) plus the node's
    longitude, 15 degrees per hour, from its right ascension (of the true equator
    and mean equinox of date, as SGP4's) less the Greenwich mean sidereal time,
    as `SunSynchronousOrbit` defines it.
    """
    epoch = timescale.from_datetime(orbit.get_epoch())
    universal_hours = np.mod(epoch.ut1 + 0.5, 1) * 24
    longitude = orbit.get_right_ascension_ascending_node() - 15 * epoch.gmst
    return float(np.mod(universal_hours + longitude / 15, 24))


def to_time_of_day(hours: float) -> time:
    """Converts hours to a time of day, to the second."""
    seconds = int(round(hours * 3600)) % 86400
    return time(seconds // 3600, seconds // 60 % 60, seconds % 60)


def get_plane(orbit: GeneralPerturbationsOrbit) -> dict:
    """
    Gets the mean elements of the plane of an element set's orbit: its mean
    altitude (meters, as TAT-C defines it), inclination (degrees), local time of
    ascending node (hours), and whether it is sun-synchronous.
    """
    return {
        "mean_altitude": orbit.get_mean_altitude(),
        "inclination": orbit.get_inclination(),
        "local_time": get_local_time(orbit),
        "sun_synchronous": orbit.elements[0].is_sun_synchronous(),
    }


def get_mean_slot(orbit: GeneralPerturbationsOrbit) -> dict:
    """
    Gets the mean longitude (degrees) of an orbit's subpoint over the day after
    its epoch, and its inclination (degrees).
    """
    epoch = orbit.get_epoch()
    times = [epoch + timedelta(minutes=m) for m in range(0, 24 * 60, 10)]
    longitude = orbit.get_geographic_position(times).longitude.degrees
    mean = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(longitude)))))
    return {"longitude": float(mean), "inclination": orbit.get_inclination()}


def _wrap(degrees):
    """Wraps angles (degrees) to [-180, 180)."""
    return (np.asarray(degrees) + 180) % 360 - 180


def compare(model, tles: list[list[str]], scope: str) -> list[dict]:
    """
    Compares a cataloged orbit with element sets: a geostationary slot (see
    `compare_slot`), a plane (see `compare_plane`), or a repeat ground track
    (see `compare_track`).

    Args:
        model: The cataloged orbit.
        tles (list[list[str]]): The element sets.
        scope (str): The scope of the metrics (the spacecraft).

    Returns:
        list[dict]: The metrics, as records state them.
    """
    if isinstance(model, GeostationarySlot):
        return compare_slot(model, tles, scope)
    if isinstance(model, RepeatGroundTrackOrbit):
        return compare_track(model, tles, scope)
    return compare_plane(model, tles, scope)


def compare_plane(model, tles: list[list[str]], scope: str) -> list[dict]:
    """
    Compares the plane of a cataloged orbit (built at each set's epoch) with the
    mean elements of element sets: the errors of mean altitude (km), inclination
    (degrees), and, for a sun-synchronous orbit, local time of ascending node
    (minutes), as 95th percentiles of their absolute values.

    Args:
        model: The cataloged orbit (with a `build(epoch)` method, or a
            `RepeatGroundTrackOrbit`).
        tles (list[list[str]]): The element sets.
        scope (str): The scope of the metrics (the spacecraft).

    Returns:
        list[dict]: The metrics, as records state them.
    """
    if isinstance(model, RepeatGroundTrackOrbit):
        # pinned to its node, so the same at every set's epoch
        pinned = (model.to_gp_orbit(), model.get_local_time_ascending_node())

        def modeled(_):
            return pinned

    else:

        def modeled(epoch):
            gp = model.build(epoch).to_gp_orbit()
            return gp, get_local_time(gp)

    altitude, inclination, local_time = [], [], []
    for tle in tles:
        orbit = to_orbit(tle)
        plane = get_plane(orbit)
        gp, model_local_time = modeled(orbit.get_epoch())
        altitude.append((plane["mean_altitude"] - gp.get_mean_altitude()) / 1e3)
        inclination.append(plane["inclination"] - gp.get_inclination())
        if plane["sun_synchronous"]:
            local_time.append(
                _wrap((plane["local_time"] - model_local_time) * 15) / 15 * 60
            )
    metrics = [
        p95_metric("altitude_error_p95", altitude, scope),
        p95_metric("inclination_error_p95", inclination, scope),
    ]
    if local_time:
        metrics.append(p95_metric("local_time_error_p95", local_time, scope))
    return metrics


def compare_slot(
    model: GeostationarySlot, tles: list[list[str]], scope: str
) -> list[dict]:
    """
    Compares a cataloged geostationary slot with element sets: the errors of the
    mean subpoint longitude (degrees, of the set's orbit less the cataloged
    orbit's, over the day after the set's epoch) and inclination (degrees).
    """
    longitude, inclination = [], []
    for tle in tles:
        orbit = to_orbit(tle)
        slot = get_mean_slot(orbit)
        # the cataloged orbit's own mean longitude over the same day
        modeled = get_mean_slot(model.build(orbit.get_epoch()).to_gp_orbit())
        longitude.append(_wrap(slot["longitude"] - modeled["longitude"]))
        inclination.append(slot["inclination"] - model.inclination)
    return [
        p95_metric("slot_longitude_error_p95", longitude, scope),
        p95_metric("inclination_error_p95", inclination, scope),
    ]


def compare_track(
    model: RepeatGroundTrackOrbit, tles: list[list[str]], scope: str
) -> list[dict]:
    """
    Compares an orbit maintained on a repeat ground track with element sets: the
    errors of its plane (see `compare_plane`) and of its ascending nodes (see
    `compare_nodes`, without the sets after maneuvers), as 95th percentiles of
    their absolute values.
    """
    delta_time, delta_longitude = compare_nodes(model, drop_maneuvers(tles))
    return [
        *compare_plane(model, tles, scope),
        p95_metric("node_longitude_error_p95", delta_longitude, scope),
        p95_metric("node_time_error_p95", delta_time, scope),
    ]


def get_nodes(orbit: GeneralPerturbationsOrbit) -> pd.DataFrame:
    """
    Gets the ascending nodes of an element set's orbit, propagated with its own
    drag, within `NODE_WINDOW` of its epoch.
    """
    epoch = orbit.get_epoch()
    return orbit.get_ascending_nodes(epoch - NODE_WINDOW, epoch + NODE_WINDOW)


def compare_nodes(
    model: RepeatGroundTrackOrbit, tles: list[list[str]]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compares the ascending nodes of element sets with the nearest ascending nodes
    of a cataloged orbit maintained on a repeat ground track.

    Args:
        model (RepeatGroundTrackOrbit): The cataloged orbit.
        tles (list[list[str]]): The element sets.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: The errors of time (seconds, the
        set's less the model's) and of longitude along the equator (km, the
        set's less the model's, positive east).
    """
    orbit = model.to_gp_orbit()
    period = 86400 * model.days / model.revolutions
    delta_time, delta_longitude = [], []
    for tle in tles:
        nodes = get_nodes(to_orbit(tle))
        modeled = orbit.get_ascending_nodes(
            nodes.time.min().to_pydatetime() - timedelta(seconds=period),
            nodes.time.max().to_pydatetime() + timedelta(seconds=period),
        )
        for node in nodes.itertuples():
            offsets = (modeled.time - node.time).dt.total_seconds().to_numpy()
            nearest = int(np.argmin(np.abs(offsets)))
            delta_time.append(-offsets[nearest])
            delta_longitude.append(
                float(_wrap(node.longitude - modeled.longitude.iloc[nearest]))
                * KM_PER_DEGREE
            )
    return np.asarray(delta_time), np.asarray(delta_longitude)


def fit_ascending_node(
    tles: list[list[str]], revolutions: int, days: int
) -> AscendingNode:
    """
    Fits the ascending node of an orbit maintained on a repeat ground track to the
    nodes of element sets: starting from a node of the first set, it moves the
    node by the median errors of time and longitude of every set's nodes (the
    middle of the band within which the orbit is maintained), until they vanish
    (within `NODE_FIT_ITERATIONS`, or else with a warning).

    Args:
        tles (list[list[str]]): The element sets of the fitting period.
        revolutions (int): The orbits in the repeat cycle.
        days (int): The days in the repeat cycle.

    Returns:
        AscendingNode: The node, rounded to a millisecond and 1e-6 degree.
    """
    first = get_nodes(to_orbit(tles[0])).iloc[0]
    node = AscendingNode(
        time=first.time.to_pydatetime(), longitude=float(first.longitude)
    )
    for _ in range(NODE_FIT_ITERATIONS):
        model = RepeatGroundTrackOrbit(
            revolutions=revolutions, days=days, ascending_node=node
        )
        delta_time, delta_longitude = compare_nodes(model, tles)
        shift_time = float(np.median(delta_time))
        shift_longitude = float(np.median(delta_longitude)) / KM_PER_DEGREE
        moment = node.time + timedelta(seconds=round(shift_time, 3))
        node = AscendingNode(
            time=moment.replace(microsecond=moment.microsecond // 1000 * 1000),
            longitude=round(float(_wrap(node.longitude + shift_longitude)), 6),
        )
        if abs(shift_time) < 0.01 and abs(shift_longitude) < 1e-5:
            return node
    warnings.warn(
        f"the ascending node did not converge in {NODE_FIT_ITERATIONS} iterations "
        f"(its last shift was {shift_time:.3f} s and {shift_longitude:.2e} deg)",
        stacklevel=2,
    )
    return node


def drop_maneuvers(tles: list[list[str]]) -> list[list[str]]:
    """
    Drops element sets within `MANEUVER_SETTLING` after a presumed maneuver: a
    jump of the semimajor axis of more than `MANEUVER_JUMP` between consecutive
    sets (by epoch).
    """
    orbits = sorted(
        ((to_orbit(tle), tle) for tle in tles), key=lambda item: item[0].get_epoch()
    )
    kept, last_maneuver, previous = [], None, None
    for orbit, tle in orbits:
        epoch = orbit.get_epoch()
        if (
            previous is not None
            and abs(orbit.get_semimajor_axis() - previous.get_semimajor_axis())
            > MANEUVER_JUMP
        ):
            last_maneuver = epoch
        if last_maneuver is None or epoch - last_maneuver > MANEUVER_SETTLING:
            kept.append(tle)
        previous = orbit
    return kept


@functools.cache
def _file_days(norad_id: int) -> dict[str, pathlib.Path]:
    """Gets the days (UTC) of a satellite's element sets in the files."""
    days = {}
    for path in sorted(ELEMENT_SET_FOLDER.glob("*.3le")):
        for epoch, _ in _parse_element_sets(path, norad_id):
            days.setdefault(f"{epoch:%Y-%m-%d}", path)
    return days


@functools.cache
def _registered_days(spacecraft_id: str) -> dict[str, list[str]]:
    """
    Gets a spacecraft's registered element sets by the days (UTC) of their own
    epochs (the registry lists them by the days of the data that instrument
    validations used them for, which may differ).
    """
    sets = {}
    for tle in _read_registry().get(spacecraft_id, {}).values():
        sets.setdefault(f"{get_epoch(tle):%Y-%m-%d}", list(tle))
    return sets


def get_available_days(spacecraft_id: str) -> list[str]:
    """
    Gets the days of the epochs of a spacecraft's element sets, registered or in
    the files of element sets (NOAA's archives of the JPSS satellites have every
    day, so they are not listed; see `ARCHIVES`).
    """
    norad_id = get_spacecraft(spacecraft_id).norad_id
    return sorted(set(_registered_days(spacecraft_id)) | set(_file_days(norad_id)))


def get_element_sets(
    spacecraft_id: str, days: list[str], directory: pathlib.Path
) -> dict[str, list[str]]:
    """
    Gets the element sets of a spacecraft for days: the registered set of each
    day's epoch, or else the one nearest the middle of the day (from NOAA's
    archives or the files), by the day (UTC) of its own epoch, without repeats
    (as the nearest set of several days may be the same) or sets whose epochs
    are outside the days' span.

    Args:
        spacecraft_id (str): The identifier of the spacecraft.
        days (list[str]): The days (ISO 8601).
        directory (pathlib.Path): Where to cache element sets downloaded from
            NOAA's archives.

    Returns:
        dict[str, list[str]]: The element sets, by the days of their epochs, in
        order.
    """
    if not days:
        return {}
    first, last = min(days), max(days)
    sets = {}
    for day in days:
        moment = datetime.fromisoformat(day).replace(tzinfo=timezone.utc)
        tle = _registered_days(spacecraft_id).get(day)
        if tle is None:
            try:
                tle = get_element_set(spacecraft_id, moment, directory=directory)
            except FileNotFoundError:
                continue
        epoch = get_epoch(tle)
        if first <= f"{epoch:%Y-%m-%d}" <= last:
            sets.setdefault(epoch, tle)
    return {f"{epoch:%Y-%m-%d}": sets[epoch] for epoch in sorted(sets)}
