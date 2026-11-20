"""
Methods to collect the observed footprints of instruments on aircraft along
their trajectories.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon

from ..schemas import Aircraft
from .check import (
    _check_platforms,
    _combine_results,
    _get_instrument_indices,
    _is_single,
)
from .ground_track import (
    _add_altaz_columns,
    _clip_track,
    _compute_step_footprints,
    _get_half_spans,
    _get_mask_geometry,
)


def _get_empty_airborne_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for airborne track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "aircraft": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "valid_obs": pd.Series([], dtype="bool"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


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
    centered on each time.

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
            footprint center.
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
    single = _is_single(aircraft, instrument_index=instrument_index)
    aircraft = _check_platforms(aircraft, Aircraft, "aircraft")
    if mask is not None:
        mask = _get_mask_geometry(mask)
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
        )
        for vehicle in aircraft
        for index in _get_instrument_indices(vehicle, instrument_index, "aircraft")
    ]
    return _combine_results(tracks, single, "time", _get_empty_airborne_track)


def _collect_airborne_ground_track(
    aircraft: Aircraft,
    times: list[datetime] | None,
    instrument_index: int,
    elevation: float,
    mask: Polygon | MultiPolygon | None,
    aircraft_altaz: bool,
    solar_altaz: bool,
    time_step: timedelta | None,
) -> gpd.GeoDataFrame:
    """
    Collect an aircraft instrument's footprints (see
    `collect_airborne_ground_track`), with a mask split along the
    anti-meridian and poles (see `_get_mask_geometry`).
    """
    trajectory = aircraft.trajectory
    if times is None:
        times = list(trajectory.time)
    elif len(times) > 0:
        times = [
            time for time, inside in zip(times, trajectory.contains(times)) if inside
        ]
    if len(times) == 0:
        return _get_empty_airborne_track()
    instrument = aircraft.get_instrument(instrument_index)
    half_spans = (
        _get_half_spans(times)
        if time_step is None
        else np.full(len(times), time_step.total_seconds() / 2)
    )
    track = aircraft.get_track(times)
    # the aircraft never stops observing, so each footprint's validity is
    # evaluated at its own time
    geometries, _ = _compute_step_footprints(
        trajectory, instrument, times, half_spans, elevation, track
    )
    target = instrument.compute_footprint_center(track, elevation)
    valid_obs = np.atleast_1d(
        instrument.is_valid_observation(track, target)
    ) & ~shapely.is_empty(np.array(geometries, dtype=object))
    results = gpd.GeoDataFrame(
        {
            "time": pd.to_datetime(times, utc=True),
            "aircraft": aircraft.name,
            "instrument": instrument.name,
            "valid_obs": valid_obs,
        },
        geometry=geometries,
        crs="EPSG:4326",
    )
    _add_altaz_columns(
        results, track, target, "aircraft" if aircraft_altaz else None, solar_altaz
    )
    return _clip_track(results, mask)


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
    track = collect_airborne_ground_track(
        aircraft, times, instrument_index, elevation, mask, time_step=time_step
    )
    track = track[track["valid_obs"]]
    if track.empty:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    return track[["geometry"]].dissolve().reset_index(drop=True)
