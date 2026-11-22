"""
Methods shared by the tracks of the instruments of any observing platform
(satellites, aircraft, and fixed surface platforms): projecting their
footprints over time steps, and assembling, clipping, and combining their
tracks and pixels.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import cached_property
from typing import Any, Protocol

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
from .._warnings import warn
from ..constants import de421
from ..schemas import AllInstruments, Instrument, PointedInstrument
from ..utils.geometry import split_polygon
from ..utils.glint import compute_glint_angle
from ..utils.time import _index_orbit_track, _to_time, _to_time_from_offsets
from .check import (
    _check_distinct_names,
    _check_platforms,
    _check_times,
    _combine_results,
    _get_instrument_indices,
    _get_platform_kind,
    _is_single,
    _PlatformKind,
)

Mask = Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries
"""A mask, in WGS84 (lon/lat) coordinates, to constrain the results of an analysis."""


class _Platform(Protocol):
    """
    An observing platform of any kind (`Satellite`, `Aircraft`,
    `SurfacePlatform`, or `RadarStation`): named in results, with a track
    (its position and velocity at times) from which its instruments project
    their views.
    """

    name: str

    def get_track(self, times: Any) -> Geocentric:
        """Gets the platform's track at times."""
        ...  # pylint: disable=unnecessary-ellipsis


class _InstrumentedPlatform(_Platform, Protocol):
    """
    An observing platform that carries instruments (a `Satellite`,
    `Aircraft`, or `SurfacePlatform`).
    """

    def get_instrument(self, index: int = 0) -> Instrument:
        """Gets an instrument, bound to the platform's attitude."""
        ...  # pylint: disable=unnecessary-ellipsis


def _get_platform_altaz(platform_altaz: bool, sat_altaz: bool | None) -> bool:
    """
    Gets whether to include a platform's angles, given the deprecated
    `sat_altaz` argument of the satellite analyses.

    Args:
        platform_altaz (bool): `True` to include the platform's angles.
        sat_altaz (bool | None): The deprecated name of `platform_altaz`, or
            `None` if not given.

    Returns:
        bool: `True` to include the platform's angles.
    """
    if sat_altaz is None:
        return platform_altaz
    warn(
        "sat_altaz is deprecated and will be removed in a future release: use "
        "platform_altaz (whose angles are named platform_alt and platform_az, "
        "formerly sat_alt and sat_az)",
        DeprecationWarning,
    )
    return platform_altaz or sat_altaz


def _get_angle_columns(platform_altaz: bool, solar_altaz: bool) -> list[str]:
    """
    Gets the names of the angle columns of a track (see `_add_altaz_columns`).

    Args:
        platform_altaz (bool): `True` to include the platform's angles
            (`platform_alt` and `platform_az`).
        solar_altaz (bool): `True` to include the solar and glint angles.

    Returns:
        list[str]: the angle columns
    """
    columns = ["platform_alt", "platform_az"] if platform_altaz else []
    if solar_altaz:
        columns += ["solar_alt", "solar_az", "glint_angle"]
    return columns


def _get_empty_track(
    kind: _PlatformKind,
    angles: list[str] | tuple[str, ...] = (),
    timed: bool = True,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for the track of a platform's instruments.

    Args:
        kind (_PlatformKind): The kind of platform, which names the column
            of the platform.
        angles (list[str] | tuple[str, ...]): The angle columns requested
            (see `_get_angle_columns`).
        timed (bool): `True` to include a `time` column.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {}
    if timed:
        columns["time"] = pd.Series([], dtype="datetime64[ns, utc]")
    columns[kind.column] = pd.Series([], dtype="str")
    columns["instrument"] = pd.Series([], dtype="str")
    columns["valid_obs"] = pd.Series([], dtype="bool")
    columns["geometry"] = pd.Series([], dtype="object")
    for name in angles:
        columns[name] = pd.Series([], dtype="float")
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _get_mask_geometry(mask: Mask) -> Polygon | MultiPolygon:
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


@dataclass(frozen=True)
class _FootprintRequest:
    """
    The arguments of a collection of footprints (or pixels), shared by the
    platforms of every kind (see `_collect_tracks`).

    Attributes:
        times (list[datetime.datetime] | None): The times, or `None` for the
            times of each aircraft's trajectory or the static footprints of
            fixed platforms.
        instrument_index (int | None): The index of the observing instrument
            on each platform, or `None` for every instrument of each platform.
        elevation (float | None): The elevation (meters) above the WGS 84
            ellipsoid at which to compute footprints, or `None` for the
            default of each kind (see `get_elevation`).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            The mask, split along the anti-meridian and poles (see
            `_get_mask_geometry`), if any.
        platform_altaz (bool): `True` to include the platforms' angles.
        solar_altaz (bool): `True` to include the solar and glint angles.
        time_step (datetime.timedelta | None): The time step that each time
            represents (see `_get_step_half_spans`).
        pixels (bool): `True` to collect pixel centers rather than footprints.
    """

    times: list[datetime] | None
    instrument_index: int | None = 0
    elevation: float | None = None
    mask: Polygon | MultiPolygon | None = None
    platform_altaz: bool = False
    solar_altaz: bool = False
    time_step: timedelta | None = None
    pixels: bool = False

    @classmethod
    def create(
        cls,
        times: list[datetime] | None,
        instrument_index: int | None = 0,
        elevation: float | None = None,
        mask: Mask | None = None,
        platform_altaz: bool = False,
        solar_altaz: bool = False,
        time_step: timedelta | None = None,
        pixels: bool = False,
    ) -> _FootprintRequest:
        """
        Creates a request, checking its times (see `_check_times`) and
        splitting its mask (see `_get_mask_geometry`).

        Args:
            times (list[datetime.datetime] | None): The times.
            instrument_index (int | None): The instrument index.
            elevation (float | None): The elevation (meters).
            mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                The mask, in WGS84 (lon/lat) coordinates, if any.
            platform_altaz (bool): `True` to include the platforms' angles.
            solar_altaz (bool): `True` to include the solar and glint angles.
            time_step (datetime.timedelta | None): The time step.
            pixels (bool): `True` to collect pixel centers.

        Returns:
            _FootprintRequest: the request
        """
        _check_times(times)
        return cls(
            times,
            instrument_index,
            elevation,
            None if mask is None else _get_mask_geometry(mask),
            platform_altaz,
            solar_altaz,
            time_step,
            pixels,
        )

    @property
    def angles(self) -> list[str]:
        """The angle columns requested (see `_get_angle_columns`)."""
        return _get_angle_columns(self.platform_altaz, self.solar_altaz)

    def get_elevation(self, platform: _Platform) -> float:
        """
        Gets the elevation at which to compute a platform's footprints: the
        requested elevation or, if `None`, the default of the platform's
        kind (see `tatc.analysis.check._PlatformKind.default_elevation`):
        the WGS 84 ellipsoid for satellites and aircraft, and the ground
        elevation of a surface platform. A radar station has no default.

        Args:
            platform (_Platform): The platform.

        Returns:
            float: the elevation (meters) above the WGS 84 ellipsoid
        """
        if self.elevation is not None:
            return self.elevation
        kind = _get_platform_kind(platform)
        if kind is None or kind.default_elevation is None:
            label = "these platforms" if kind is None else f"{kind.label}s"
            raise ValueError(
                f"elevation is required to collect the footprints of {label}"
            )
        return kind.default_elevation(platform)

    @cached_property
    def skyfield_times(self) -> Time:
        """
        The times as Skyfield times, shared (with their Earth orientation)
        by every platform.
        """
        return _to_time(self.times)

    @cached_property
    def half_spans(self) -> npt.NDArray:
        """
        The half span (seconds) of the step of each time (see
        `_get_step_half_spans`).
        """
        return _get_step_half_spans(self.times, self.time_step)  # type: ignore


_Collector = Callable[[Any, Any, _FootprintRequest], Any]
"""
Collects the track (or pixels) of a platform's instrument (by index, or
`None` for a platform without instruments), returning a
`geopandas.GeoDataFrame`, or `None` if it is empty.
"""


def _collect_tracks(
    platforms: object,
    kind: _PlatformKind,
    request: _FootprintRequest,
    collect: _Collector,
    empty: Callable[[bool], gpd.GeoDataFrame] | None = None,
) -> gpd.GeoDataFrame:
    """
    Collects the tracks (or pixels) of the instruments of platforms of a
    kind: checks the platforms (see `_check_platforms`, with distinct names,
    see `_check_distinct_names`), collects the track of each platform's
    instruments (or of each platform without instruments), and combines
    them (see `_combine_results`): for a single platform and instrument, its
    track; otherwise, those of each platform and instrument, concatenated
    and sorted by time (or, without times, in the order of the platforms and
    their instruments).

    Args:
        platforms (object): The platform(s).
        kind (_PlatformKind): The kind of platform.
        request (_FootprintRequest): The request, of which no times give an
            empty track.
        collect (_Collector): Collects the track of a platform's instrument.
        empty (Callable[[bool], geopandas.GeoDataFrame] | None): Gets an
            empty track, with a `time` column if timed, by default that of
            `_get_empty_track` with the angles requested.

    Returns:
        geopandas.GeoDataFrame: the tracks
    """
    single = _is_single(platforms, instrument_index=request.instrument_index)
    checked = _check_platforms(platforms, kind)
    _check_distinct_names(checked, kind, kind.argument)
    # fixed platforms without times have static footprints
    timed = not kind.fixed or request.times is not None

    def get_empty() -> gpd.GeoDataFrame:
        if empty is not None:
            return empty(timed)
        return _get_empty_track(kind, request.angles, timed=timed)

    if request.times is not None and len(request.times) == 0:
        return get_empty()
    tracks = []
    for platform in checked:
        indices = (
            _get_instrument_indices(platform, request.instrument_index)
            if kind.instrumented
            else [None]
        )
        for index in indices:
            track = collect(platform, index, request)
            tracks.append(get_empty() if track is None else track)
    if timed:
        return _combine_results(tracks, single, "time", get_empty)
    tracks = [track for track in tracks if not track.empty]
    if len(tracks) == 0:
        return get_empty()
    return pd.concat(tracks, ignore_index=True)


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


def _get_observing_runs(
    platform: _Platform,
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
        platform (Satellite | Aircraft): The platform, whose track (see
            `get_track`) is evaluated at the refined times.
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
        track = platform.get_track(
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
    platform: _Platform,
    instrument: AllInstruments,
    times: list[datetime],
    half_spans: npt.NDArray,
    elevation: float,
    track: Geocentric | None = None,
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
        platform (Satellite | Aircraft): The platform, whose track (see
            `get_track`) is evaluated at the times of each step.
        instrument (AllInstruments): The instrument, bound to its attitude.
        times (list[datetime.datetime]): The times.
        half_spans (numpy.typing.NDArray): The half span (seconds) of each
            time's step (see `_get_half_spans`).
        elevation (float): The elevation (meters) at which to project the
            footprints.
        track (skyfield.positionlib.Geocentric | None): The platform's track
            at the times, if already computed.

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
        instant_track = (
            platform.get_track(_to_time([times[i] for i in instant]))
            if track is None
            else _index_orbit_track(track, instant)
        )
        footprints[instant] = instrument.compute_footprint(
            instant_track, elevation=elevation
        )
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
    sub_track = platform.get_track(
        _to_time_from_offsets(reference, seconds[step] + offset)
    )
    observe = np.reshape(instrument.is_observing(sub_track), -1)
    if np.all(observe):
        footprints[index] = instrument.compute_swept_footprints(
            sub_track, counts, elevation=elevation
        )
        return list(footprints), represent
    # trim the steps to their runs of observing times
    offsets = np.split(offset, np.cumsum(counts)[:-1])
    observing = np.split(observe, np.cumsum(counts)[:-1])
    runs = _get_observing_runs(
        platform, instrument, offsets, observing, [times[i] for i in index]
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
            platform.get_track(_to_time(run_times)),
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


def _assemble_track(
    platform: _Platform,
    instrument: Instrument,
    times: list[datetime] | pd.DatetimeIndex | None,
    track: Geocentric,
    geometries: list[Polygon | MultiPolygon],
    elevation: float,
    request: _FootprintRequest,
    represent: list[datetime | None] | None = None,
    ground_elevation: float = 0,
) -> gpd.GeoDataFrame:
    """
    Assembles the track of a platform's instrument (a satellite, an aircraft,
    or a fixed surface platform) from its footprints at times: one row per
    time, with the validity of each observation (see
    `Instrument.is_valid_observation`, at the footprint center), invalid for
    an empty footprint, and the angles requested (see `_add_altaz_columns`),
    clipped to a mask and sorted by time (see `_clip_track`).

    Args:
        platform (Satellite | Aircraft | SurfacePlatform): The platform, which
            names the track (in the column of its kind), and whose track (see
            `get_track`) evaluates validity at other times.
        instrument (Instrument): The instrument, bound to the platform's
            attitude.
        times (list[datetime.datetime] | pandas.DatetimeIndex | None): The
            times, or `None` for the static footprints of a fixed platform
            without times (valid unless empty, without evaluating
            time-varying requirements).
        track (skyfield.positionlib.Geocentric): The platform's track at the
            times (at one time if `None`).
        geometries (list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]):
            The footprint at each time.
        elevation (float): The elevation (meters) of the footprints.
        request (_FootprintRequest): The request, of the mask (if any) and
            the angles to include.
        represent (list[datetime.datetime | None] | None): For each footprint,
            the time at which to evaluate its validity if not its own (see
            `_compute_step_footprints`).
        ground_elevation (float): The elevation (meters) of the ground, which
            blocks rays toward an elevation above the instrument.

    Returns:
        geopandas.GeoDataFrame: the track
    """
    kind = _get_platform_kind(platform)
    # compute targets
    target = instrument.compute_footprint_center(track, elevation, ground_elevation)
    # determine observation validity
    if times is None:
        valid_obs = np.ones(len(geometries), dtype=bool)
    else:
        valid_obs = np.atleast_1d(instrument.is_valid_observation(track, target))
    moved = [i for i, time in enumerate(represent or []) if time is not None]
    if moved:
        # swept footprints whose time falls in a blind interval are valid as
        # at the middle of their longest run of observing times
        moved_track = platform.get_track(
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
            **({} if times is None else {"time": list(times)}),
            kind.column: platform.name,
            "instrument": instrument.name,
            "valid_obs": valid_obs,
        },
        geometry=list(geometries),
        crs="EPSG:4326",
    )
    _add_altaz_columns(
        frame, track, target, request.platform_altaz, request.solar_altaz
    )
    return _clip_track(frame, request.mask)


def _assemble_pixels(
    platform: _Platform,
    instrument: Instrument,
    times: list[datetime] | pd.DatetimeIndex | None,
    track: Geocentric,
    elevation: float,
    request: _FootprintRequest,
    ground_elevation: float = 0,
    region: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame | None:
    """
    Assembles the pixels of a platform's instrument (a satellite, an
    aircraft, or a fixed surface platform) at times: one row per pixel
    center (see `PointedInstrument.compute_footprint_pixel_array`), with the
    validity of each observation (see `Instrument.is_valid_observation`, at
    the pixel) and the angles requested (see `_add_altaz_columns`), seen
    from the pixel, clipped to a mask and sorted by time (see `_clip_track`).

    Args:
        platform (Satellite | Aircraft | SurfacePlatform): The platform,
            which names the pixels (in the column of its kind).
        instrument (Instrument): The instrument, bound to the platform's
            attitude: a rectangular `PointedInstrument`.
        times (list[datetime.datetime] | pandas.DatetimeIndex | None): The
            times, or `None` for the pixels of a fixed platform without
            times (all valid, without evaluating time-varying requirements).
        track (skyfield.positionlib.Geocentric): The platform's track at the
            times (at one time if `None`).
        elevation (float): The elevation (meters) of the pixels.
        request (_FootprintRequest): The request, of the mask (if any) and
            the angles to include.
        ground_elevation (float): The elevation (meters) of the ground, which
            blocks rays toward an elevation above the instrument.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            A region (such as that visible over a terrain mask) outside which
            pixels are omitted.

    Returns:
        geopandas.GeoDataFrame | None: the pixels, or `None` if there are none
    """
    _check_pixel_instrument(instrument)
    kind = _get_platform_kind(platform)
    arrays = instrument.compute_footprint_pixel_array(
        track, elevation, ground_elevation=ground_elevation
    )
    points, index = shapely.get_parts(np.array(arrays, dtype=object), return_index=True)
    if region is not None:
        inside = shapely.intersects(region, shapely.force_2d(points))
        points, index = points[inside], index[inside]
    if len(points) == 0:
        return None
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
    if times is None:
        timed = {}
    else:
        time_list = list(times)
        timed = {"time": [time_list[i] for i in index]}
    frame = gpd.GeoDataFrame(
        {
            **timed,
            kind.column: platform.name,
            "instrument": instrument.name,
            "valid_obs": valid_obs,
        },
        geometry=list(points),
        crs="EPSG:4326",
    )
    _add_altaz_columns(frame, rows, target, request.platform_altaz, request.solar_altaz)
    return _clip_track(frame, request.mask)


def _assemble(
    platform: _Platform,
    instrument: Instrument,
    times: list[datetime] | pd.DatetimeIndex | None,
    track: Geocentric,
    elevation: float,
    request: _FootprintRequest,
    geometries: list[Polygon | MultiPolygon] | None = None,
    represent: list[datetime | None] | None = None,
    ground_elevation: float = 0,
    region: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame | None:
    """
    Assembles the pixels (see `_assemble_pixels`) of a platform's instrument
    if the request is of pixels, or otherwise its track from its footprints
    (see `_assemble_track`).

    Args:
        platform (_Platform): The platform.
        instrument (Instrument): The instrument, bound to the platform's
            attitude.
        times (list[datetime.datetime] | pandas.DatetimeIndex | None): The
            times, or `None` for a fixed platform without times.
        track (skyfield.positionlib.Geocentric): The platform's track at the
            times.
        elevation (float): The elevation (meters) of the footprints or
            pixels.
        request (_FootprintRequest): The request.
        geometries (list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon] | None):
            The footprint at each time (unused for pixels).
        represent (list[datetime.datetime | None] | None): For each
            footprint, the time at which to evaluate its validity if not its
            own (see `_compute_step_footprints`).
        ground_elevation (float): The elevation (meters) of the ground, which
            blocks rays toward an elevation above the instrument.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            A region outside which pixels are omitted (footprints are
            clipped to it beforehand).

    Returns:
        geopandas.GeoDataFrame | None: the track, or the pixels (`None` if
            there are none)
    """
    if request.pixels:
        return _assemble_pixels(
            platform,
            instrument,
            times,
            track,
            elevation,
            request,
            ground_elevation=ground_elevation,
            region=region,
        )
    return _assemble_track(
        platform,
        instrument,
        times,
        track,
        geometries,  # type: ignore
        elevation,
        request,
        represent=represent,
        ground_elevation=ground_elevation,
    )


def _check_pixel_instrument(instrument: Instrument) -> None:
    """
    Checks that an instrument has pixels: a rectangular `PointedInstrument`,
    whose `cross_track_pixels` x `along_track_pixels` grid defines its pixel
    array, raising a `ValueError` otherwise.

    Args:
        instrument (Instrument): The instrument.
    """
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )


def _get_viewing_geometry(
    track: Geocentric, target: GeographicPosition
) -> dict[str, npt.NDArray]:
    """
    Gets the viewing geometry of a platform from targets, for any kind of
    platform: its elevation angle (`platform_alt`) and azimuth
    (`platform_az`, degrees) seen from each target, its slant range
    (`platform_range`), and its height above the WGS 84 ellipsoid
    (`platform_height`, meters).

    Args:
        track (skyfield.positionlib.Geocentric): The platform's track (see
            `get_track`) at the time of each target.
        target (skyfield.toposlib.GeographicPosition): The targets.

    Returns:
        dict[str, numpy.typing.NDArray]: the geometry, by column
    """
    altitude, azimuth, distance = (track - target.at(track.t)).altaz()
    return {
        "platform_alt": np.atleast_1d(altitude.degrees),  # type: ignore
        "platform_az": np.atleast_1d(azimuth.degrees),  # type: ignore
        "platform_range": np.atleast_1d(distance.m),  # type: ignore
        "platform_height": np.atleast_1d(wgs84.height_of(track).m),  # type: ignore
    }


def _add_altaz_columns(
    track: gpd.GeoDataFrame,
    orbit_track: Geocentric,
    target: GeographicPosition,
    platform_altaz: bool,
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
        platform_altaz (bool): `True` to append the platform's
            altitude/azimuth angles (`platform_alt` and `platform_az`).
        solar_altaz (bool): `True` to append the solar altitude/azimuth
            angles (`solar_alt` and `solar_az`) and the glint angle
            (`glint_angle`, see `tatc.utils.compute_glint_angle`).
    """
    angles = {}
    if platform_altaz:
        geometry = _get_viewing_geometry(orbit_track, target)
        angles["platform_alt"] = geometry["platform_alt"]
        angles["platform_az"] = geometry["platform_az"]
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
