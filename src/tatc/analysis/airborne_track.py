"""
Methods to collect the observed footprints of instruments on aircraft along
their trajectories.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from ..schemas import Aircraft
from .check import (
    _check_platforms,
    _check_times,
    _combine_results,
    _dissolve_footprints,
    _get_instrument_indices,
    _is_single,
)
from .ground_track import (
    _assemble_pixels,
    _assemble_track,
    _compute_step_footprints,
    _get_angle_columns,
    _get_empty_track,
    _get_mask_geometry,
    _get_step_half_spans,
)


def collect_airborne_ground_track(
    aircraft: Aircraft | list[Aircraft],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    aircraft_altaz: bool = False,
    solar_altaz: bool = False,
    time_step: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observed footprints of instruments on aircraft at times
    along their trajectories, projected to an elevation by computing the
    exact ray/WGS 84 ellipsoid intersection, as for a satellite's ground
    track (see `collect_ground_track`). The views of a `PointedInstrument`
    are rotated by the aircraft's attitude, if its trajectory gives one (see
    `tatc.schemas.Aircraft`). The footprint of a swept instrument (see
    `PointedInstrument.is_swept`, as for a cross-track scanner or a pushbroom
    imager) is the strip its cross-track line sweeps over the time step
    centered on each time, trimmed (to within 0.1 seconds) to the times that
    the trajectory contains (see `tatc.schemas.Trajectory.contains`): it does
    not extend before the first or after the last state, or into a gap longer
    than `max_gap`.

    The footprint is the geometric overlap of an instrument's field of view
    with the elevation (for example, the ground or a layer of a storm):
    whether the instrument detects anything there (its sensitivity, or
    attenuation along the beam) is not modeled.

    Args:
        aircraft (Aircraft | list[Aircraft]): The observing aircraft.
        times (list[datetime.datetime] | None): The times to sample, of
            which those that the trajectory of an aircraft does not contain
            (before its first or after its last state, or within a gap
            longer than its `max_gap`, see `tatc.schemas.Trajectory.contains`)
            are omitted for it. Defaults to the times of each aircraft's
            trajectory.
        instrument_index (int | None): The index of the observing instrument
            on each aircraft, or `None` for every instrument of each aircraft.
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
            at which to compute footprints, in the same reference as the
            aircraft's altitudes.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
            An optional mask, always interpreted in WGS84 (lon/lat)
            coordinates, to constrain results.
        aircraft_altaz (bool): `True` to include the altitude/azimuth angles
            of the aircraft seen from the footprint center (`aircraft_alt`
            and `aircraft_az`), as `sat_altaz` does for a satellite.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
            angles (`solar_alt` and `solar_az`) and the glint angle
            (`glint_angle`, see `tatc.utils.compute_glint_angle`) at the
            footprint center. The angles are `NaN` for an empty footprint
            (as of a level at or above the aircraft).
        time_step (datetime.timedelta | None): The time step that each time
            represents, over which a swept instrument's footprint sweeps, or
            `None` to take half the smaller gap to the neighboring times on
            either side of each time (so that a gap between times, as between
            flight legs, does not stretch the footprints beside it).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected airborne track
            results, with columns `time`, `aircraft`, `instrument`,
            `valid_obs`, and `geometry` (and the angles requested): for a
            single aircraft and instrument (an integer `instrument_index`),
            sorted by time; otherwise, those of each aircraft and
            instrument, concatenated and sorted by time.
    """
    return _collect_airborne(
        aircraft,
        times,
        instrument_index,
        elevation,
        mask,
        aircraft_altaz,
        solar_altaz,
        time_step,
    )


def _collect_airborne(
    aircraft: Aircraft | list[Aircraft],
    times: list[datetime] | None,
    instrument_index: int | None,
    elevation: float,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None,
    aircraft_altaz: bool,
    solar_altaz: bool,
    time_step: timedelta | None,
    pixels: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the footprints (see `collect_airborne_ground_track`) or, if
    `pixels`, the pixels (see `tatc.analysis.footprints.collect_footprints`)
    of instruments on aircraft.
    """
    single = _is_single(aircraft, instrument_index=instrument_index)
    aircraft = _check_platforms(aircraft, Aircraft, "aircraft")
    _check_times(times)
    if mask is not None:
        mask = _get_mask_geometry(mask)
    angles = _get_angle_columns("aircraft" if aircraft_altaz else None, solar_altaz)
    tracks = [
        _collect_airborne_ground_track(
            vehicle,
            times,
            index,
            elevation,
            mask,
            aircraft_altaz,
            solar_altaz,
            time_step,
            pixels,
        )
        for vehicle in aircraft
        for index in _get_instrument_indices(vehicle, instrument_index, "aircraft")
    ]
    return _combine_results(
        tracks, single, "time", lambda: _get_empty_track("aircraft", angles)
    )


def _collect_airborne_ground_track(
    aircraft: Aircraft,
    times: list[datetime] | None,
    instrument_index: int,
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    aircraft_altaz: bool,
    solar_altaz: bool,
    time_step: timedelta | None,
    pixels: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect an aircraft instrument's footprints (see
    `collect_airborne_ground_track`) or pixels, with a mask split along the
    anti-meridian and poles (see `_get_mask_geometry`).
    """
    platform_altaz = "aircraft" if aircraft_altaz else None
    trajectory = aircraft.trajectory
    if times is None:
        times = list(trajectory.time)
    elif len(times) > 0:
        times = [
            time for time, inside in zip(times, trajectory.contains(times)) if inside
        ]
    if len(times) == 0:
        return _get_empty_track(
            "aircraft", _get_angle_columns(platform_altaz, solar_altaz)
        )
    instrument = aircraft.get_instrument(instrument_index)
    track = aircraft.get_track(times)
    if pixels:
        return _assemble_pixels(
            "aircraft",
            aircraft.name,
            instrument,
            pd.to_datetime(times, utc=True),
            track,
            elevation,
            mask,
            platform_altaz,
            solar_altaz,
        )
    half_spans = _get_step_half_spans(times, time_step)
    # the instrument observes at each time (which the trajectory contains), so
    # each footprint's validity is evaluated at its own time
    geometries, represent = _compute_step_footprints(
        trajectory, instrument, times, half_spans, elevation, track
    )
    return _assemble_track(
        "aircraft",
        aircraft.name,
        instrument,
        pd.to_datetime(times, utc=True),
        track,
        geometries,
        elevation,
        mask,
        platform_altaz,
        solar_altaz,
        represent=represent,
        path=trajectory,
    )


def compute_airborne_ground_track(
    aircraft: Aircraft | list[Aircraft],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_step: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on aircraft:
    unlike `collect_airborne_ground_track` (one footprint per time), this
    dissolves the footprints of valid observations into a single geometry
    (for example, the area observed over a flight).

    Args:
        aircraft (Aircraft | list[Aircraft]): The observing aircraft.
        times (list[datetime.datetime] | None): The times to sample (see
            `collect_airborne_ground_track`). Defaults to the times of each
            aircraft's trajectory.
        instrument_index (int | None): The index of the observing instrument
            on each aircraft, or `None` for every instrument of each aircraft.
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
            at which to compute footprints.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
            An optional mask, always interpreted in WGS84 (lon/lat)
            coordinates, to constrain results.
        time_step (datetime.timedelta | None): The time step that each time
            represents (see `collect_airborne_ground_track`).

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated airborne
            track: one row with the observed area, or none if no observation
            is valid.
    """
    return _dissolve_footprints(
        collect_airborne_ground_track(
            aircraft, times, instrument_index, elevation, mask, time_step=time_step
        )
    )
