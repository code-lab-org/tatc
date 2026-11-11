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

from .. import config
from ..constants import de421
from ..schemas import (
    AllInstruments,
    GeneralPerturbationsOrbit,
    PointedInstrument,
    Satellite,
)
from ..utils.geometry import split_polygon
from ..utils.glint import compute_glint_angle
from ..utils.time import _index_orbit_track, _to_time
from .check import (
    _check_satellite,
    _check_satellites,
    _combine_results,
    _get_instrument_indices,
    _is_single,
)
from .region_sampling import compute_region_access_periods


def _get_empty_ground_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for ground track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "valid_obs": pd.Series([], dtype="bool"),
        "geometry": pd.Series([], dtype="object"),
    }
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
    orbit: GeneralPerturbationsOrbit,
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
        orbit (GeneralPerturbationsOrbit): The orbit.
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
    orbit: GeneralPerturbationsOrbit,
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
        orbit (GeneralPerturbationsOrbit): The orbit.
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
    offsets = [
        np.linspace(
            -half_spans[i],
            half_spans[i],
            int(np.ceil(2 * half_spans[i] / max_substep)) + 1,
        )
        for i in index
    ]
    sub_times = [
        times[i] + timedelta(seconds=float(s))
        for i, offset in zip(index, offsets)
        for s in offset
    ]
    track = orbit.get_orbit_track_at_time(_to_time(sub_times))
    observe = np.reshape(instrument.is_observing(track), -1)
    counts = [len(offset) for offset in offsets]
    if np.all(observe):
        footprints[index] = instrument.compute_swept_footprints(
            track, counts, elevation=elevation
        )
        return list(footprints), represent
    # trim the steps to their runs of observing times
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
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's viewable ground footprint at each requested
    time, projected to a specified elevation by computing the exact
    ray/WGS-84-geoid intersection (see
    `tatc.utils.projection.compute_footprint`).

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
                the propagated orbit to the (buffered) region of interest,
                improving performance for small areas.
        sat_altaz (bool): `True` to include satellite altitude/azimuth angles
                for the sub-satellite point.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles
                and the glint angle (see `tatc.utils.compute_glint_angle`)
                for the sub-satellite point.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    single = _is_single(satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    if len(times) == 0:
        return _get_empty_ground_track()
    if mask is not None:
        mask = _get_mask_geometry(mask)
    # the times, shared by every satellite (with their Earth orientation)
    t = _to_time(times)
    tracks = [
        _collect_ground_track(
            satellite, times, t, index, elevation, mask, sat_altaz, solar_altaz
        )
        for satellite in satellites
        for index in _get_instrument_indices(satellite, instrument_index)
    ]
    return _combine_results(tracks, single, "time", _get_empty_ground_track)


def _collect_ground_track(
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
    Collect a satellite's ground track at times, given as `datetime`s and as
    Skyfield times (see `collect_ground_track`), with a mask split along the anti-meridian
    and poles (see `_get_mask_geometry`).
    """
    # select the observing instrument
    instrument = satellite.get_instrument(instrument_index)
    orbit = satellite.orbit.to_gp_orbit()
    # the span of each time's step, over which a swept instrument's footprint
    # sweeps
    half_spans = _get_half_spans(times)
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask,
        # reusing the footprints computed to cull it
        culled = _cull_orbit_track(
            satellite, instrument, times, half_spans, mask, elevation
        )
        if culled is None:
            return _get_empty_ground_track()
        orbit_track, geometries, represent = culled
    else:
        # propagate orbit
        orbit_track = orbit.get_orbit_track_at_time(t)
        # compute footprints (exact ray/WGS-84-geoid intersection)
        geometries, represent = _compute_step_footprints(
            orbit, instrument, times, half_spans, elevation, orbit_track
        )
    # compute targets
    target = instrument.compute_footprint_center(orbit_track, elevation)
    # determine observation validity
    valid_obs = np.atleast_1d(instrument.is_valid_observation(orbit_track, target))
    moved = [i for i, time in enumerate(represent) if time is not None]
    if moved:
        # swept footprints whose time falls in a blind interval are valid as
        # at the middle of their longest run of observing times
        moved_track = orbit.get_orbit_track_at_time(
            _to_time([represent[i] for i in moved])
        )
        valid_obs = valid_obs.copy()
        valid_obs[moved] = np.atleast_1d(
            instrument.is_valid_observation(
                moved_track,
                instrument.compute_footprint_center(moved_track, elevation),
            )
        )
    if isinstance(instrument, PointedInstrument) and instrument.is_swept:
        # swept footprints entirely within blind intervals are not valid
        valid_obs = valid_obs & ~shapely.is_empty(np.array(geometries, dtype=object))
    records = [
        {
            "time": time,
            "satellite": satellite.name,
            "instrument": instrument.name,
            "valid_obs": valid_obs[i],
            "geometry": geometries[i],
        }
        for i, time in enumerate(orbit_track.t.utc_datetime())  # type: ignore
    ]
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if sat_altaz:
        # append satellite altitude/azimuth columns
        sat_altaz_data = (orbit_track - target.at(orbit_track.t)).altaz()
        track["sat_alt"] = sat_altaz_data[0].degrees  # type: ignore
        track["sat_az"] = sat_altaz_data[1].degrees  # type: ignore
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz_data = (
            (de421["earth"] + target)
            .at(orbit_track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        track["solar_alt"] = solar_altaz_data[0].degrees  # type: ignore
        track["solar_az"] = solar_altaz_data[1].degrees  # type: ignore
        # append glint angle column (at the footprint center)
        track["glint_angle"] = np.atleast_1d(compute_glint_angle(orbit_track, target))

    if mask is not None:
        track = gpd.clip(track, mask)
    # sort rows by time (clipping does not preserve their order), keeping
    # the original order of rows at the same time
    return track.sort_index().sort_values("time", kind="stable").reset_index(drop=True)


def compute_ground_track(
    satellite: Satellite,
    times: list[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    dissolve_orbits: bool = True,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated ground track for a satellite of interest: unlike
    `collect_ground_track` (one footprint polygon per time step), this
    dissolves the valid-observation footprints within each orbit into a
    single geometry per orbit, and optionally dissolves across orbits into
    one geometry for the entire `times` range.

    Args:
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        dissolve_orbits (bool): True, to aggregate multiple orbits in one output.
    Returns:
        GeoDataFrame: The data frame of aggregated ground track results.
    """
    _check_satellite(satellite)
    track = collect_ground_track(satellite, times, instrument_index, elevation, mask)
    if not track.empty:
        # assign orbit identifier, from the period of the element that
        # propagates the orbit at each time (see `get_element_index`)
        orbit = satellite.orbit.to_gp_orbit()
        track["orbit_id"] = [
            (time - times[0])
            // orbit.elements[orbit.get_element_index(time)].get_orbit_period()
            for time in track.time
        ]
        # filter to valid observations and dissolve
        track = track[track.valid_obs].dissolve(by="orbit_id").reset_index(drop=True)
    if dissolve_orbits:
        track = track.dissolve()
    return track


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
    grid defines the pixel array.

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
                the propagated orbit to the (buffered) region of interest,
                improving performance for small areas.
        sat_altaz (bool): `True` to include satellite altitude/azimuth angles.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles
                and the glint angle (see `tatc.utils.compute_glint_angle`).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground pixels results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    single = _is_single(satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    if len(times) == 0:
        return _get_empty_ground_track()
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
    return _combine_results(tracks, single, "time", _get_empty_ground_track)


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
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask
        # pixel arrays are views at an instant, so cull by the views at the
        # times (with steps of no span)
        culled = _cull_orbit_track(
            satellite, instrument, times, np.zeros(len(times)), mask, elevation
        )
        if culled is None:
            return _get_empty_ground_track()
        orbit_track, _, _ = culled
    else:
        # propagate orbit
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    # compute the footprint pixel array
    geometries = instrument.compute_footprint_pixel_array(
        orbit_track,
        elevation,
    )
    # construct results as a list of (time index, record) tuples, keeping the
    # source orbit_track index alongside each record so satellite/solar altaz
    # below can index directly into orbit_track rather than re-deriving it
    indexed_records = [
        (
            i,
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "valid_obs": instrument.is_valid_observation(
                    orbit_track[i], wgs84.latlon(point.y, point.x, point.z)
                ).all(),
                "geometry": point,
            },
        )
        for i, time in enumerate(orbit_track.t.utc_datetime())  # type: ignore
        for point in geometries[i].geoms
    ]
    records = [record for _, record in indexed_records]
    # build geodataframe
    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if sat_altaz:
        # append satellite altitude/azimuth columns
        sat_altaz_data = [
            (
                orbit_track[i]
                - wgs84.latlon(
                    record["geometry"].y, record["geometry"].x, record["geometry"].z
                ).at(
                    orbit_track.t[i]
                )  # type: ignore
            ).altaz()
            for i, record in indexed_records
        ]
        gdf["sat_alt"] = [altaz[0].degrees for altaz in sat_altaz_data]  # type: ignore
        gdf["sat_az"] = [altaz[1].degrees for altaz in sat_altaz_data]  # type: ignore
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz_data = [
            (
                de421["earth"]
                + wgs84.latlon(
                    record["geometry"].y, record["geometry"].x, record["geometry"].z
                )
            )
            .at(orbit_track.t[i])  # type: ignore
            .observe(de421["sun"])
            .apparent()
            .altaz()
            for i, record in indexed_records
        ]
        gdf["solar_alt"] = [altaz[0].degrees for altaz in solar_altaz_data]  # type: ignore
        gdf["solar_az"] = [altaz[1].degrees for altaz in solar_altaz_data]  # type: ignore
        # append glint angle column (at each pixel)
        gdf["glint_angle"] = [
            float(
                compute_glint_angle(
                    orbit_track[i],
                    wgs84.latlon(
                        record["geometry"].y, record["geometry"].x, record["geometry"].z
                    ),
                )
            )
            for i, record in indexed_records
        ]

    if mask is not None:
        gdf = gpd.clip(gdf, mask)
    # sort rows by time (clipping does not preserve their order), keeping
    # the original order of rows at the same time
    return gdf.sort_index().sort_values("time", kind="stable").reset_index(drop=True)
