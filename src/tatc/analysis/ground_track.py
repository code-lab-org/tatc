"""
Methods to collect instrument ground tracks.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time
from skyfield.toposlib import GeographicPosition

from .. import config
from ..constants import de421
from ..schemas import (
    AllInstruments,
    GeneralPerturbationsOrbit,
    Instrument,
    PointedInstrument,
    Satellite,
    Trajectory,
)
from ..utils.geometry import split_polygon
from ..utils.glint import compute_glint_angle
from ..utils.time import _index_orbit_track, _to_time, _to_time_from_offsets
from .check import (
    _check_satellites,
    _check_times,
    _combine_results,
    _dissolve_footprints,
    _get_instrument_indices,
    _is_single,
)
from .region_sampling import compute_region_access_periods


def _get_angle_columns(platform: str | None, solar_altaz: bool) -> list[str]:
    """
    Gets the names of the angle columns of a track (see `_add_altaz_columns`).

    Args:
        platform (str | None): The prefix of the platform's columns, or
            `None` to omit them.
        solar_altaz (bool): `True` to include the solar and glint angles.

    Returns:
        list[str]: the angle columns
    """
    columns = [] if platform is None else [f"{platform}_alt", f"{platform}_az"]
    if solar_altaz:
        columns += ["solar_alt", "solar_az", "glint_angle"]
    return columns


def _get_empty_track(
    column: str = "satellite",
    angles: list[str] | tuple[str, ...] = (),
    timed: bool = True,
    extra: dict[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for the track of a platform's instruments.

    Args:
        column (str): The column that names the platform (see
            `tatc.analysis.check._PLATFORM_COLUMNS`).
        angles (list[str] | tuple[str, ...]): The angle columns requested
            (see `_get_angle_columns`).
        timed (bool): `True` to include a `time` column.
        extra (dict[str, str] | None): Other columns (and their types),
            after `instrument`.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {}
    if timed:
        columns["time"] = pd.Series([], dtype="datetime64[ns, utc]")
    columns[column] = pd.Series([], dtype="str")
    columns["instrument"] = pd.Series([], dtype="str")
    for name, dtype in (extra or {}).items():
        columns[name] = pd.Series([], dtype=dtype)
    columns["valid_obs"] = pd.Series([], dtype="bool")
    columns["geometry"] = pd.Series([], dtype="object")
    for name in angles:
        columns[name] = pd.Series([], dtype="float")
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _get_mask_geometry(
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries,
) -> Polygon | MultiPolygon:
    """
    Gets a mask as a single geometry, split along the anti-meridian and
    poles (see `tatc.utils.geometry.split_polygon`) so that it compares
    correctly with footprints in the standard longitude range.

    Args:
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries):
                The mask, in WGS84 (lon/lat) coordinates.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: the mask geometry
    """
    if isinstance(mask, (Polygon, MultiPolygon)):
        geometry = mask
    elif isinstance(mask, gpd.GeoDataFrame):
        geometry = mask.dissolve().iloc[0].geometry
    else:
        geometry = mask.union_all()
    return split_polygon(geometry)


_BLIND_TOLERANCE = 0.1  # seconds, to which swept footprints are trimmed


def _get_half_spans(times: list[datetime]) -> npt.NDArray:
    """
    Gets the half span (seconds) of the time step of each time, over which a
    swept instrument's footprint sweeps (see `PointedInstrument.is_swept`):
    half the smaller of the gaps to the neighboring times (in time order), so
    that a gap between samples does not stretch the footprints beside it.

    Args:
        times (list[datetime.datetime]): The times.

    Returns:
        numpy.typing.NDArray: The half spans (seconds), zero for a single
            time or a repeated time.
    """
    seconds = np.array([(time - times[0]).total_seconds() for time in times])
    order = np.argsort(seconds, kind="stable")
    gaps = np.diff(seconds[order])
    half = np.minimum(np.append(np.inf, gaps), np.append(gaps, np.inf)) / 2
    half[~np.isfinite(half)] = 0
    result = np.empty(len(times))
    result[order] = half
    return result


def _get_observing_runs(
    orbit: GeneralPerturbationsOrbit | Trajectory,
    instrument: PointedInstrument,
    offsets: list[npt.NDArray],
    observing: list[npt.NDArray],
    times: list[datetime],
) -> list[list[npt.NDArray]]:
    """
    Gets the runs of each time step in which an instrument observes (outside
    the blind intervals of its attitude), as offsets (seconds) from the step's
    time, with the bounds of each run refined by bisection to within
    `_BLIND_TOLERANCE`, on the side on which the instrument observes.

    Args:
        orbit (GeneralPerturbationsOrbit | Trajectory): The orbit, or the
            trajectory of an aircraft.
        instrument (PointedInstrument): The instrument.
        offsets (list[numpy.typing.NDArray]): The offsets of the sub-times of
            each step.
        observing (list[numpy.typing.NDArray]): Whether the instrument
            observes at each sub-time of each step.
        times (list[datetime.datetime]): The time of each step.

    Returns:
        list[list[numpy.typing.NDArray]]: The offsets of the times of each
            run of each step.
    """
    # the transitions between consecutive sub-times, refined all at once
    steps, lows, highs, rises = [], [], [], []
    for k, (offset, observe) in enumerate(zip(offsets, observing)):
        for j in np.flatnonzero(observe[1:] != observe[:-1]):
            steps.append(k)
            lows.append(offset[j])
            highs.append(offset[j + 1])
            rises.append(bool(observe[j + 1]))
    low, high = np.array(lows, dtype=float), np.array(highs, dtype=float)
    rise = np.array(rises, dtype=bool)
    while len(low) > 0 and np.max(high - low) > _BLIND_TOLERANCE:
        middle = (low + high) / 2
        track = orbit.get_orbit_track_at_time(
            _to_time(
                [times[k] + timedelta(seconds=float(m)) for k, m in zip(steps, middle)]
            )
        )
        observe = np.reshape(instrument.is_observing(track), -1)
        # keep the observing side of each transition at the bound it rises or
        # falls past
        moves_low = observe != rise
        low = np.where(moves_low, middle, low)
        high = np.where(moves_low, high, middle)
    # the bounds of each step's runs, on their observing sides
    bounds: list[list[tuple[int, float]]] = [[] for _ in offsets]
    for k, lo, hi, up in zip(steps, low, high, rise):
        bounds[k].append((1, hi) if up else (-1, lo))
    runs = []
    for k, (offset, observe) in enumerate(zip(offsets, observing)):
        step_runs, start = [], offset[0] if observe[0] else None
        transitions = iter(bounds[k])
        for j in range(len(offset) - 1):
            if observe[j + 1] != observe[j]:
                kind, bound = next(transitions)
                if kind == 1:
                    start = bound
                else:
                    inner = offset[(offset > start) & (offset < bound)]
                    step_runs.append(np.concatenate([[start], inner, [bound]]))
                    start = None
        if start is not None:
            inner = offset[(offset > start) & (offset < offset[-1])]
            step_runs.append(np.concatenate([[start], inner, [offset[-1]]]))
        runs.append([run for run in step_runs if run[-1] - run[0] > 0])
    return runs


def _compute_step_footprints(
    orbit: GeneralPerturbationsOrbit | Trajectory,
    instrument: AllInstruments,
    times: list[datetime],
    half_spans: npt.NDArray,
    elevation: float,
    orbit_track: Geocentric | None = None,
) -> tuple[list[Polygon | MultiPolygon], list[datetime | None]]:
    """
    Computes an instrument's footprint over the time step of each time: its
    view at the time or, for a swept instrument (see
    `PointedInstrument.is_swept`), the strip its cross-track line sweeps over
    the time's half span on either side, trimmed to the times at which it
    observes (outside the blind intervals of its attitude). Also gets the time
    at which to evaluate each footprint's validity: the time itself or, for a
    swept footprint whose time falls in a blind interval, the middle of its
    longest run of observing times.

    Args:
        orbit (GeneralPerturbationsOrbit | Trajectory): The orbit, or the
            trajectory of an aircraft.
        instrument (AllInstruments): The instrument, bound to its attitude.
        times (list[datetime.datetime]): The times.
        half_spans (numpy.typing.NDArray): The half span (seconds) of each
            time's step (see `_get_half_spans`).
        elevation (float): The elevation (meters) at which to project the
            footprints.
        orbit_track (skyfield.positionlib.Geocentric | None): The orbit track
            at the times, if already propagated.

    Returns:
        tuple[list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon], list[datetime.datetime | None]]:
            the footprints and, for each, the time at which to evaluate its
            validity if not its own (otherwise `None`).
    """
    swept = (
        np.asarray(half_spans) > 0
        if isinstance(instrument, PointedInstrument) and instrument.is_swept
        else np.zeros(len(times), dtype=bool)
    )
    footprints = np.empty(len(times), dtype=object)
    represent: list[datetime | None] = [None] * len(times)
    instant = np.flatnonzero(~swept)
    if len(instant) > 0:
        track = (
            orbit.get_orbit_track_at_time(_to_time([times[i] for i in instant]))
            if orbit_track is None
            else _index_orbit_track(orbit_track, instant)
        )
        footprints[instant] = instrument.compute_footprint(track, elevation=elevation)
    index = np.flatnonzero(swept)
    if len(index) == 0:
        return list(footprints), represent
    max_substep = config.get_rc().swept_footprint_max_substep
    # the offsets (seconds) of the evenly spaced sub-times of each step, from
    # the start to the end of its span, for all steps at once
    half = np.asarray(half_spans, dtype=float)[index]
    counts = np.ceil(2 * half / max_substep).astype(int) + 1
    step = np.repeat(np.arange(len(index)), counts)
    position = np.arange(np.sum(counts)) - np.repeat(np.cumsum(counts) - counts, counts)
    offset = half[step] * (2 * position / (counts[step] - 1) - 1)
    # the sub-times, as offsets from the earliest time (without building a
    # datetime for each)
    reference = min(times[i] for i in index)
    seconds = np.array([(times[i] - reference).total_seconds() for i in index])
    track = orbit.get_orbit_track_at_time(
        _to_time_from_offsets(reference, seconds[step] + offset)
    )
    observe = np.reshape(instrument.is_observing(track), -1)
    if np.all(observe):
        footprints[index] = instrument.compute_swept_footprints(
            track, counts, elevation=elevation
        )
        return list(footprints), represent
    # trim the steps to their runs of observing times
    offsets = np.split(offset, np.cumsum(counts)[:-1])
    observing = np.split(observe, np.cumsum(counts)[:-1])
    runs = _get_observing_runs(
        orbit, instrument, offsets, observing, [times[i] for i in index]
    )
    run_times, run_counts, owners = [], [], []
    for k, (i, step_runs) in enumerate(zip(index, runs)):
        for run in step_runs:
            run_times.extend(times[i] + timedelta(seconds=float(s)) for s in run)
            run_counts.append(len(run))
            owners.append(k)
        if step_runs and not any(run[0] <= 0 <= run[-1] for run in step_runs):
            # evaluate validity in the middle of the longest run
            longest = max(step_runs, key=lambda run: run[-1] - run[0])
            represent[i] = times[i] + timedelta(
                seconds=float((longest[0] + longest[-1]) / 2)
            )
    pieces = (
        instrument.compute_swept_footprints(
            orbit.get_orbit_track_at_time(_to_time(run_times)),
            run_counts,
            elevation=elevation,
        )
        if run_counts
        else []
    )
    owners = np.array(owners, dtype=int)
    for k, i in enumerate(index):
        parts = [pieces[j] for j in np.flatnonzero(owners == k)]
        footprints[i] = (
            parts[0]
            if len(parts) == 1
            else (shapely.union_all(parts) if parts else Polygon())
        )
    return list(footprints), represent


def _cull_orbit_track(
    satellite: Satellite,
    instrument: AllInstruments,
    times: list[datetime],
    half_spans: npt.NDArray,
    mask: Polygon | MultiPolygon,
    elevation: float,
) -> tuple[Geocentric, list[Polygon | MultiPolygon], list[datetime | None]] | None:
    """
    Propagates the orbit track only at the times when the instrument's
    footprint intersects a mask. The times are first culled to those whose
    step (the time, or its span for a swept instrument; see
    `_get_half_spans`) overlaps a period when the instrument's field of regard
    may observe any part of the mask (see
    `tatc.analysis.region_sampling.compute_region_access_periods`), so that
    the orbit is not propagated far from the mask, and then to those whose
    footprint intersects the mask.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (AllInstruments): The observing instrument.
        times (list[datetime.datetime]): The times to sample.
        half_spans (numpy.typing.NDArray): The half span (seconds) of each
            time's step.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The mask, in WGS84 (lon/lat) coordinates (see `_get_mask_geometry`).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
                at which to project footprints.

    Returns:
        tuple[skyfield.positionlib.Geocentric, list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon], list[datetime.datetime | None]] | None:
            the culled orbit track, its footprints, and the times at which to
            evaluate their validity if not their own (see
            `_compute_step_footprints`), or None if no footprint intersects
            the mask
    """
    periods = compute_region_access_periods(
        mask,
        satellite,
        min(times),
        max(times),
        instrument.field_of_regard,
        elevation,
    )
    if periods.empty:
        return None
    # select the times whose steps overlap a period (the periods are sorted
    # and disjoint)
    sample = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    span = pd.to_timedelta(np.asarray(half_spans, dtype=float), unit="s")
    left = pd.DatetimeIndex([period.left for period in periods])
    right = pd.DatetimeIndex([period.right for period in periods])
    index = np.searchsorted(left, sample + span, side="right") - 1
    selected = (index >= 0) & (sample - span <= right[np.maximum(index, 0)])
    if not any(selected):
        return None
    chosen = np.flatnonzero(selected)
    orbit = satellite.orbit.to_gp_orbit()
    orbit_track = orbit.get_orbit_track([times[i] for i in chosen])
    # cull orbit track to observable footprint
    footprint, represent = _compute_step_footprints(
        orbit,
        instrument,
        [times[i] for i in chosen],
        np.asarray(half_spans)[chosen],
        elevation,
        orbit_track,
    )
    keep = np.flatnonzero([mask.intersects(f) for f in footprint])
    if len(keep) == 0:
        return None
    # keep the Earth orientation quantities cached on the times
    return (
        _index_orbit_track(orbit_track, keep),
        [footprint[i] for i in keep],
        [represent[i] for i in keep],
    )


def collect_ground_track(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    sat_altaz: bool = False,
    solar_altaz: bool = False,
    time_step: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's viewable ground footprint at each requested
    time, projected to a specified elevation by computing the exact
    ray/WGS 84 ellipsoid intersection (see
    `tatc.utils.projection.compute_footprint`). The footprint of a swept
    instrument (see `PointedInstrument.is_swept`) is the strip its
    cross-track line sweeps over the time step centered on each time, trimmed
    to the times at which it observes (see `tatc.schemas.Attitude.blind_intervals`).

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        instrument_index (int | None): The index of the observing instrument
                in each satellite, or `None` for every instrument of each
                satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results. Providing a mask limits
                the propagated orbit to the periods when the instrument may
                observe the mask (see `compute_region_access_periods`),
                improving performance for small areas.
        sat_altaz (bool): `True` to include the altitude/azimuth angles of
                the satellite (`sat_alt` and `sat_az`) seen from the
                footprint center.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at the
                footprint center. The angles are `NaN` for an empty
                footprint.
        time_step (datetime.timedelta | None): The time step that each time
                represents, over which a swept instrument's footprint sweeps,
                or `None` to take half the smaller gap to the neighboring
                times on either side of each time (so that a gap between
                times does not stretch the footprints beside it, and a single
                time is a view at an instant).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    single = _is_single(satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    _check_times(times)
    angles = _get_angle_columns("sat" if sat_altaz else None, solar_altaz)
    if len(times) == 0:
        return _get_empty_track("satellite", angles)
    if mask is not None:
        mask = _get_mask_geometry(mask)
    # the times, shared by every satellite (with their Earth orientation)
    t = _to_time(times)
    # the span of each time's step, over which a swept instrument's footprint
    # sweeps
    half_spans = _get_step_half_spans(times, time_step)
    tracks = [
        _collect_ground_track(
            satellite,
            times,
            t,
            half_spans,
            index,
            elevation,
            mask,
            sat_altaz,
            solar_altaz,
        )
        for satellite in satellites
        for index in _get_instrument_indices(satellite, instrument_index)
    ]
    return _combine_results(
        tracks, single, "time", lambda: _get_empty_track("satellite", angles)
    )


def _get_step_half_spans(
    times: list[datetime], time_step: timedelta | None
) -> npt.NDArray:
    """
    Gets the half span (seconds) of the time step of each time: half the
    given time step, or, if `None`, half the smaller gap to the neighboring
    times (see `_get_half_spans`).

    Args:
        times (list[datetime.datetime]): The times.
        time_step (datetime.timedelta | None): The time step.

    Returns:
        numpy.typing.NDArray: The half spans (seconds).
    """
    if time_step is None:
        return _get_half_spans(times)
    return np.full(len(times), time_step.total_seconds() / 2)


def _collect_ground_track(
    satellite: Satellite,
    times: list[datetime],
    t: Time,
    half_spans: npt.NDArray,
    instrument_index: int,
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    sat_altaz: bool,
    solar_altaz: bool,
) -> gpd.GeoDataFrame:
    """
    Collect a satellite's ground track at times, given as `datetime`s and as
    Skyfield times (see `collect_ground_track`), with a mask split along the anti-meridian
    and poles (see `_get_mask_geometry`).
    """
    # select the observing instrument
    instrument = satellite.get_instrument(instrument_index)
    orbit = satellite.orbit.to_gp_orbit()
    platform_altaz = "sat" if sat_altaz else None
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask,
        # reusing the footprints computed to cull it
        culled = _cull_orbit_track(
            satellite, instrument, times, half_spans, mask, elevation
        )
        if culled is None:
            return _get_empty_track(
                "satellite", _get_angle_columns(platform_altaz, solar_altaz)
            )
        orbit_track, geometries, represent = culled
    else:
        # propagate orbit
        orbit_track = orbit.get_orbit_track_at_time(t)
        # compute footprints (exact ray/WGS 84 ellipsoid intersection)
        geometries, represent = _compute_step_footprints(
            orbit, instrument, times, half_spans, elevation, orbit_track
        )
    return _assemble_track(
        "satellite",
        satellite.name,
        instrument,
        orbit_track.t.utc_datetime(),  # type: ignore
        orbit_track,
        geometries,
        elevation,
        mask,
        platform_altaz,
        solar_altaz,
        represent=represent,
        path=orbit,
    )


def _assemble_track(
    column: str,
    name: str,
    instrument: Instrument,
    times: list[datetime] | pd.DatetimeIndex,
    track: Geocentric,
    geometries: list[Polygon | MultiPolygon],
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    platform_altaz: str | None,
    solar_altaz: bool,
    represent: list[datetime | None] | None = None,
    path: GeneralPerturbationsOrbit | Trajectory | None = None,
    ground_elevation: float = 0,
    extra: dict[str, object] | None = None,
) -> gpd.GeoDataFrame:
    """
    Assembles the track of a platform's instrument (a satellite, an aircraft,
    or a fixed surface platform) from its footprints at times: one row per
    time, with the validity of each observation (see
    `Instrument.is_valid_observation`, at the footprint center), invalid for
    an empty footprint, and the angles requested (see `_add_altaz_columns`),
    clipped to a mask and sorted by time (see `_clip_track`).

    Args:
        column (str): The column that names the platform (see
            `tatc.analysis.check._PLATFORM_COLUMNS`).
        name (str): The name of the platform.
        instrument (Instrument): The instrument, bound to the platform's
            attitude.
        times (list[datetime.datetime] | pandas.DatetimeIndex): The times.
        track (skyfield.positionlib.Geocentric): The platform's track at the
            times.
        geometries (list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]):
            The footprint at each time.
        elevation (float): The elevation (meters) of the footprints.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            The mask (see `_get_mask_geometry`), if any.
        platform_altaz (str | None): The prefix of the platform's angle
            columns, or `None` to omit them.
        solar_altaz (bool): `True` to include the solar and glint angles.
        represent (list[datetime.datetime | None] | None): For each footprint,
            the time at which to evaluate its validity if not its own (see
            `_compute_step_footprints`).
        path (GeneralPerturbationsOrbit | Trajectory | None): The orbit or
            trajectory of the platform, to evaluate validity at other times.
        ground_elevation (float): The elevation (meters) of the ground, which
            blocks rays toward an elevation above the instrument.
        extra (dict[str, object] | None): Other columns, after `instrument`.

    Returns:
        geopandas.GeoDataFrame: the track
    """
    # compute targets
    target = instrument.compute_footprint_center(track, elevation, ground_elevation)
    # determine observation validity
    valid_obs = np.atleast_1d(instrument.is_valid_observation(track, target))
    moved = [i for i, time in enumerate(represent or []) if time is not None]
    if moved and path is not None:
        # swept footprints whose time falls in a blind interval are valid as
        # at the middle of their longest run of observing times
        moved_track = path.get_orbit_track_at_time(
            _to_time([represent[i] for i in moved])  # type: ignore
        )
        valid_obs = valid_obs.copy()
        valid_obs[moved] = np.atleast_1d(
            instrument.is_valid_observation(
                moved_track,
                instrument.compute_footprint_center(
                    moved_track, elevation, ground_elevation
                ),
            )
        )
    # empty footprints (as of a view that misses the elevation, or a swept
    # footprint entirely within blind intervals) are not valid
    valid_obs = valid_obs & ~shapely.is_empty(np.array(geometries, dtype=object))
    frame = gpd.GeoDataFrame(
        {
            "time": list(times),
            column: name,
            "instrument": instrument.name,
            **(extra or {}),
            "valid_obs": valid_obs,
        },
        geometry=list(geometries),
        crs="EPSG:4326",
    )
    _add_altaz_columns(frame, track, target, platform_altaz, solar_altaz)
    return _clip_track(frame, mask)


def _assemble_pixels(
    column: str,
    name: str,
    instrument: Instrument,
    times: list[datetime] | pd.DatetimeIndex | None,
    track: Geocentric,
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    platform_altaz: str | None,
    solar_altaz: bool,
    ground_elevation: float = 0,
    extra: dict[str, float] | None = None,
    region: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Assembles the pixels of a platform's instrument (a satellite, an
    aircraft, or a fixed surface platform) at times: one row per pixel
    center (see `PointedInstrument.compute_footprint_pixel_array`), with the
    validity of each observation (see `Instrument.is_valid_observation`, at
    the pixel) and the angles requested (see `_add_altaz_columns`), seen
    from the pixel, clipped to a mask and sorted by time (see `_clip_track`).

    Args:
        column (str): The column that names the platform.
        name (str): The name of the platform.
        instrument (Instrument): The instrument, bound to the platform's
            attitude: a rectangular `PointedInstrument`.
        times (list[datetime.datetime] | pandas.DatetimeIndex | None): The
            times, or `None` for the pixels of a fixed platform without
            times (all valid, without evaluating time-varying requirements).
        track (skyfield.positionlib.Geocentric): The platform's track at the
            times (at one time if `None`).
        elevation (float): The elevation (meters) of the pixels.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            The mask (see `_get_mask_geometry`), if any.
        platform_altaz (str | None): The prefix of the platform's angle
            columns, or `None` to omit them.
        solar_altaz (bool): `True` to include the solar and glint angles.
        ground_elevation (float): The elevation (meters) of the ground, which
            blocks rays toward an elevation above the instrument.
        extra (dict[str, float] | None): Other (numeric) columns, after
            `instrument`.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            A region (such as that visible over a terrain mask) outside which
            pixels are omitted.

    Returns:
        geopandas.GeoDataFrame: the pixels
    """
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )
    arrays = instrument.compute_footprint_pixel_array(
        track, elevation, ground_elevation=ground_elevation
    )
    points, index = shapely.get_parts(np.array(arrays, dtype=object), return_index=True)
    if region is not None:
        inside = shapely.intersects(region, shapely.force_2d(points))
        points, index = points[inside], index[inside]
    if len(points) == 0:
        return _get_empty_track(
            column,
            _get_angle_columns(platform_altaz, solar_altaz),
            timed=times is not None,
            extra={key: "float" for key in extra or {}},
        )
    height = shapely.get_z(points)
    target = wgs84.latlon(
        shapely.get_y(points),
        shapely.get_x(points),
        np.where(np.isnan(height), elevation, height),
    )
    rows = _index_orbit_track(track, index)
    valid_obs = (
        np.atleast_1d(instrument.is_valid_observation(rows, target))
        if times is not None
        else np.ones(len(points), dtype=bool)
    )
    timed = {} if times is None else {"time": [list(times)[i] for i in index]}
    frame = gpd.GeoDataFrame(
        {
            **timed,
            column: name,
            "instrument": instrument.name,
            **(extra or {}),
            "valid_obs": valid_obs,
        },
        geometry=list(points),
        crs="EPSG:4326",
    )
    _add_altaz_columns(frame, rows, target, platform_altaz, solar_altaz)
    return _clip_track(frame, mask)


def _add_altaz_columns(
    track: gpd.GeoDataFrame,
    orbit_track: Geocentric,
    target: GeographicPosition,
    platform: str | None,
    solar_altaz: bool,
) -> None:
    """
    Appends the altitude/azimuth angles (degrees) of an observing platform
    and of the sun, seen from the footprint centers, to a ground track (one
    row per time): `NaN` for an empty footprint (as of a view that misses
    the elevation), which has no center. An observation that is not valid
    for another reason (such as a target that is not sunlit) keeps its
    angles.

    Args:
        track (geopandas.GeoDataFrame): The ground track, modified in place.
        orbit_track (skyfield.positionlib.Geocentric): The platform's track.
        target (skyfield.toposlib.GeographicPosition): The footprint centers.
        platform (str | None): The prefix of the platform's columns
            (`{platform}_alt` and `{platform}_az`), or `None` to omit them.
        solar_altaz (bool): `True` to append the solar altitude/azimuth
            angles (`solar_alt` and `solar_az`) and the glint angle
            (`glint_angle`, see `tatc.utils.compute_glint_angle`).
    """
    angles = {}
    if platform is not None:
        alt, az, _ = (orbit_track - target.at(orbit_track.t)).altaz()
        angles[f"{platform}_alt"] = alt.degrees  # type: ignore
        angles[f"{platform}_az"] = az.degrees  # type: ignore
    if solar_altaz:
        alt, az, _ = (
            (de421["earth"] + target)
            .at(orbit_track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        angles["solar_alt"] = alt.degrees  # type: ignore
        angles["solar_az"] = az.degrees  # type: ignore
        angles["glint_angle"] = compute_glint_angle(orbit_track, target)
    empty = track.geometry.isna().to_numpy() | track.geometry.is_empty.to_numpy()
    for name, values in angles.items():
        track[name] = np.where(empty, np.nan, np.atleast_1d(values))


def _clip_track(
    track: gpd.GeoDataFrame, mask: Polygon | MultiPolygon | None
) -> gpd.GeoDataFrame:
    """
    Clips a ground track to a mask, if any, and sorts its rows by time, if it
    has times (clipping does not preserve their order), keeping the original
    order of rows at the same time.

    Args:
        track (geopandas.GeoDataFrame): The ground track.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            The mask (see `_get_mask_geometry`).

    Returns:
        geopandas.GeoDataFrame: the clipped and sorted ground track
    """
    if mask is not None:
        track = gpd.clip(track, mask)
    track = track.sort_index()
    if "time" in track.columns:
        track = track.sort_values("time", kind="stable")
    return track.reset_index(drop=True)


def _get_empty_orbit_coverage() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for the coverage of each orbit (see
    `compute_ground_track`).

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "orbit_id": pd.Series([], dtype="int"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _compute_orbit_coverage(
    satellite: Satellite, track: gpd.GeoDataFrame, first: datetime
) -> gpd.GeoDataFrame:
    """
    Dissolves the valid footprints of a satellite's instrument within each
    orbit (see `compute_ground_track`): the number of periods, of the element
    that propagates the orbit at each time (see `get_element_index`), since
    the first time.

    Args:
        satellite (Satellite): The satellite.
        track (geopandas.GeoDataFrame): The valid footprints of an instrument.
        first (datetime.datetime): The first time.

    Returns:
        geopandas.GeoDataFrame: the coverage of each orbit
    """
    orbit = satellite.orbit.to_gp_orbit()
    orbit_id = [
        int(
            (time - first)
            // orbit.elements[orbit.get_element_index(time)].get_orbit_period()
        )
        for time in track["time"]
    ]
    track = track.assign(orbit_id=orbit_id)
    times = track.groupby("orbit_id")["time"]
    dissolved = track[["orbit_id", "geometry"]].dissolve(by="orbit_id")
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": track["instrument"].iloc[0],
            "orbit_id": dissolved.index.to_numpy(),
            "start": times.min().to_numpy(),
            "end": times.max().to_numpy(),
        },
        geometry=dissolved.geometry.to_numpy(),
        crs=track.crs,
    )


def compute_ground_track(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    dissolve_orbits: bool = True,
    time_step: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on satellites:
    unlike `collect_ground_track` (one footprint per time), this dissolves
    the footprints of valid observations into one geometry (for example,
    the area observed over the times) or, if not `dissolve_orbits`, into one
    geometry for each orbit of each satellite's instrument.

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (list[datetime.datetime]): The times to sample.
        instrument_index (int | None): The index of the observing instrument
                on each satellite, or `None` for every instrument of each
                satellite.
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
                at which to compute footprints.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        dissolve_orbits (bool): `True` to dissolve the coverage of every orbit,
                satellite, and instrument into one geometry; `False` for
                one geometry for each orbit of each satellite's instrument.
        time_step (datetime.timedelta | None): The time step that each time
                represents, over which a swept instrument's footprint sweeps
                (see `collect_ground_track`), or `None` to infer it from the
                neighboring times.

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated ground
            track: if `dissolve_orbits`, one row with the observed area (or
            none if no observation is valid); otherwise, one row for each
            orbit of each satellite's instrument with a valid observation,
            with columns `satellite`, `instrument`, `orbit_id` (the number of
            orbital periods since the first time), `start` and `end` (the
            first and last times of its valid footprints), and `geometry`,
            sorted by `start`.
    """
    satellites = _check_satellites(satellites)
    _check_times(times)
    if dissolve_orbits:
        return _dissolve_footprints(
            collect_ground_track(
                satellites,
                times,
                instrument_index,
                elevation,
                mask,
                time_step=time_step,
            )
        )
    orbits = []
    for satellite in satellites:
        for index in _get_instrument_indices(satellite, instrument_index):
            track = collect_ground_track(
                satellite, times, index, elevation, mask, time_step=time_step
            )
            track = track[track["valid_obs"]]
            if not track.empty:
                orbits.append(_compute_orbit_coverage(satellite, track, min(times)))
    if len(orbits) == 0:
        return _get_empty_orbit_coverage()
    return (
        pd.concat(orbits, ignore_index=True)
        .sort_values("start", kind="stable")
        .reset_index(drop=True)
    )


def collect_ground_pixels(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    sat_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's individual ground sample point (pixel)
    locations at each requested time, rather than an aggregate observable
    geometry (see `collect_ground_track`). Only supported for a rectangular
    `PointedInstrument`, whose `cross_track_pixels` x `along_track_pixels`
    grid defines the pixel array. For the pixels of aircraft and fixed
    surface platforms, see `tatc.analysis.footprints.collect_footprints`.

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        instrument_index (int | None): The index of the observing instrument
                in each satellite, or `None` for every instrument of each
                satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground pixels.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results. Providing a mask limits
                the propagated orbit to the periods when the instrument may
                observe the mask (see `compute_region_access_periods`),
                improving performance for small areas.
        sat_altaz (bool): `True` to include the altitude/azimuth angles of
                the satellite (`sat_alt` and `sat_az`) seen from each pixel.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at
                each pixel.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground pixels results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    single = _is_single(satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    _check_times(times)
    angles = _get_angle_columns("sat" if sat_altaz else None, solar_altaz)
    if len(times) == 0:
        return _get_empty_track("satellite", angles)
    if mask is not None:
        mask = _get_mask_geometry(mask)
    # the times, shared by every satellite (with their Earth orientation)
    t = _to_time(times)
    tracks = [
        _collect_ground_pixels(
            satellite, times, t, index, elevation, mask, sat_altaz, solar_altaz
        )
        for satellite in satellites
        for index in _get_instrument_indices(satellite, instrument_index)
    ]
    return _combine_results(
        tracks, single, "time", lambda: _get_empty_track("satellite", angles)
    )


def _collect_ground_pixels(
    satellite: Satellite,
    times: list[datetime],
    t: Time,
    instrument_index: int,
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    sat_altaz: bool,
    solar_altaz: bool,
) -> gpd.GeoDataFrame:
    """
    Collect a satellite's ground pixels at times, given as `datetime`s and as
    Skyfield times (see `collect_ground_pixels`), with a mask split along the anti-meridian
    and poles (see `_get_mask_geometry`).
    """
    # select the observing instrument
    instrument = satellite.get_instrument(instrument_index)
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )
    platform_altaz = "sat" if sat_altaz else None
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask
        # pixel arrays are views at an instant, so cull by the views at the
        # times (with steps of no span)
        culled = _cull_orbit_track(
            satellite, instrument, times, np.zeros(len(times)), mask, elevation
        )
        if culled is None:
            return _get_empty_track(
                "satellite", _get_angle_columns(platform_altaz, solar_altaz)
            )
        orbit_track, _, _ = culled
    else:
        # propagate orbit
        orbit_track = satellite.get_track(t)
    return _assemble_pixels(
        "satellite",
        satellite.name,
        instrument,
        orbit_track.t.utc_datetime(),  # type: ignore
        orbit_track,
        elevation,
        mask,
        platform_altaz,
        solar_altaz,
    )
