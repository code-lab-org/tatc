"""
Methods to collect satellite orbit tracks.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.functions import angle_between
from skyfield.timelib import Time

from ..constants import de421
from ..schemas import AllInstruments, ConicalInstrument, Satellite
from ..utils.observation import field_of_regard_to_swath_width
from ..utils.time import _index_orbit_track, _to_time
from .check import (
    _check_satellites,
    _combine_results,
    _get_instrument_indices,
    _is_single,
)


def _swath_width(
    instrument: AllInstruments, altitude: float, elevation: float
) -> float:
    """
    Gets an instrument's swath width: from its field of regard or, for a
    conical instrument, from its cone angle and scan sector.

    Args:
        instrument (AllInstruments): The observing instrument.
        altitude (float): The satellite altitude (meters).
        elevation (float): The elevation (meters) at which to project the swath.

    Returns:
        float: The swath width (meters).
    """
    if isinstance(instrument, ConicalInstrument):
        return instrument.get_swath_width(altitude - elevation)
    return field_of_regard_to_swath_width(
        altitude, instrument.field_of_regard, elevation
    )


def _get_empty_orbit_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for orbit track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "swath_width": pd.Series([], dtype="float"),
        "valid_obs": pd.Series([], dtype="bool"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


class OrbitCoordinate(str, Enum):
    """
    Enumeration of different orbit track coordinate systems.
    """

    WGS84 = "wgs84"
    ECEF = "ecef"
    ECI = "eci"


class OrbitOutput(str, Enum):
    """
    Enumeration of different orbit output options.
    """

    POSITION = "position"
    POSITION_VELOCITY = "velocity"


def collect_orbit_track(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    coordinates: OrbitCoordinate = OrbitCoordinate.WGS84,
    orbit_output: OrbitOutput = OrbitOutput.POSITION,
    sat_sunlit: bool = False,
    solar_altaz: bool = False,
    solar_beta: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the satellites' own positions (and, optionally, velocities) at
    each requested time, in a specified coordinate frame. Note this reports
    the satellite's location, not the zero-elevation ground point beneath it;
    use `collect_ground_track` for footprint/ground-projected results.

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (typing.List[datetime.datetime]): The list of times to sample.
        instrument_index (int | None): The index of the observing instrument
                in each satellite, or `None` for every instrument of each
                satellite.
        elevation (float): The elevation (meters) above the WGS 84 datum for
                which to project the instrument's field of regard into a
                swath width (`swath_width` output column). Does not affect
                the reported position itself, which is always the satellite's
                true altitude.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results to points whose
                sub-satellite longitude/latitude falls within the mask. This
                filter is applied consistently regardless of the requested
                output `coordinates`.
        coordinates (OrbitCoordinate): The coordinate system of orbit track
                points: `wgs84` (geodetic longitude/latitude/altitude, output
                CRS `EPSG:4326`), `ecef` (Earth-fixed Cartesian meters, output
                CRS `EPSG:4978`), or `eci` (inertial GCRS Cartesian meters,
                no fixed CRS since the frame is time-varying).
        orbit_output (OrbitOutput): `position` for position only, or
                `velocity` to also include a `velocity` column. Velocity is
                expressed in the same frame as `coordinates`, except `wgs84`
                velocity is given as local East/North/Up components (m/s)
                rather than a rate of change of longitude/latitude/altitude.
        sat_sunlit (bool): `True` to include whether the satellite is sunlit.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles as seen from the satellite's own position.
        solar_beta (bool): `True` to include solar beta angles.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected orbit track
            results: for a single satellite and instrument (an integer
            `instrument_index`), in the order of `times`; otherwise, those of
            each satellite and instrument, concatenated and sorted by time.
    """
    single = _is_single(satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    if len(times) == 0:
        return _get_empty_orbit_track()
    # the times, shared by every satellite (with their Earth orientation)
    t = _to_time(times)
    tracks = [
        _collect_orbit_track(
            satellite,
            t,
            index,
            elevation,
            mask,
            coordinates,
            orbit_output,
            sat_sunlit,
            solar_altaz,
            solar_beta,
        )
        for satellite in satellites
        for index in _get_instrument_indices(satellite, instrument_index)
    ]
    return _combine_results(tracks, single, "time", _get_empty_orbit_track)


def _collect_orbit_track(
    satellite: Satellite,
    t: Time,
    instrument_index: int,
    elevation: float,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None,
    coordinates: OrbitCoordinate,
    orbit_output: OrbitOutput,
    sat_sunlit: bool,
    solar_altaz: bool,
    solar_beta: bool,
) -> gpd.GeoDataFrame:
    """
    Collect a satellite's own position at Skyfield times (see
    `collect_orbit_track`).
    """
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    # propagate orbit
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    # geodetic (WGS84) position of the satellite itself (not the zero-elevation
    # subpoint below it -- elevation here is the satellite's own altitude)
    sat_pos = wgs84.geographic_position_of(orbit_track)
    if mask is not None:
        # trim orbit track to provided mask; mask is always interpreted in
        # WGS84 (lon/lat) coordinates regardless of the requested output
        # `coordinates`, so this filter applies consistently to all of them
        longitude = np.atleast_1d(sat_pos.longitude.degrees)
        latitude = np.atleast_1d(sat_pos.latitude.degrees)
        mask_contains_sat_pos = np.zeros(len(longitude), dtype=bool)
        for geometry in (
            mask.geometry
            if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
            else [mask]
        ):
            mask_contains_sat_pos |= shapely.contains_xy(geometry, longitude, latitude)
        if not any(mask_contains_sat_pos):
            return _get_empty_orbit_track()
        # keep the Earth orientation quantities cached on the times
        orbit_track = _index_orbit_track(
            orbit_track, np.flatnonzero(mask_contains_sat_pos)
        )
        # recompute sat_pos
        sat_pos = wgs84.geographic_position_of(orbit_track)
    altitude = np.atleast_1d(sat_pos.elevation.m)
    # create shapely points in proper coordinate system
    if coordinates == OrbitCoordinate.WGS84:
        points = shapely.points(
            np.column_stack(
                [
                    np.atleast_1d(sat_pos.longitude.degrees),
                    np.atleast_1d(sat_pos.latitude.degrees),
                    altitude,
                ]
            )
        )
    elif coordinates == OrbitCoordinate.ECEF:
        points = shapely.points(np.reshape(sat_pos.itrs_xyz.m, (3, -1)).T)
    else:
        points = shapely.points(np.reshape(orbit_track.xyz.m, (3, -1)).T)
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track)
    columns = {
        "time": list(np.atleast_1d(orbit_track.t.utc_datetime())),  # type: ignore
        "satellite": satellite.name,
        "instrument": instrument.name,
        "swath_width": [_swath_width(instrument, h, elevation) for h in altitude],
        "valid_obs": valid_obs,
        "geometry": points,
    }
    # create velocity points if needed
    if orbit_output != OrbitOutput.POSITION:
        # compute satellite velocity
        if coordinates == OrbitCoordinate.ECI:
            velocity = np.reshape(orbit_track.velocity.m_per_s, (3, -1))
        elif coordinates == OrbitCoordinate.ECEF:
            velocity = np.reshape(
                orbit_track.frame_xyz_and_velocity(itrs)[1].m_per_s, (3, -1)
            )
        else:
            # rotate ECEF velocity into local East/North/Up components at the
            # satellite's geodetic longitude/latitude
            ecef_velocity = np.reshape(
                orbit_track.frame_xyz_and_velocity(itrs)[1].m_per_s, (3, -1)
            )
            lon = np.radians(np.atleast_1d(sat_pos.longitude.degrees))
            lat = np.radians(np.atleast_1d(sat_pos.latitude.degrees))
            east = -np.sin(lon) * ecef_velocity[0] + np.cos(lon) * ecef_velocity[1]
            north = (
                -np.sin(lat) * np.cos(lon) * ecef_velocity[0]
                - np.sin(lat) * np.sin(lon) * ecef_velocity[1]
                + np.cos(lat) * ecef_velocity[2]
            )
            up = (
                np.cos(lat) * np.cos(lon) * ecef_velocity[0]
                + np.cos(lat) * np.sin(lon) * ecef_velocity[1]
                + np.sin(lat) * ecef_velocity[2]
            )
            velocity = np.array([east, north, up])
        columns["velocity"] = shapely.points(velocity.T)

    # tag the CRS to match the requested output coordinates: WGS84 is
    # geographic degrees, ECEF is geocentric meters, and ECI (GCRS) is an
    # inertial, time-varying frame with no fixed EPSG code
    track_crs = (
        "EPSG:4326"
        if coordinates == OrbitCoordinate.WGS84
        else "EPSG:4978" if coordinates == OrbitCoordinate.ECEF else None
    )
    track = gpd.GeoDataFrame(columns, crs=track_crs)
    if sat_sunlit:
        # append sat_sunlit column
        track["sat_sunlit"] = orbit_track.is_sunlit(de421)
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz_data = (
            (de421["earth"] + sat_pos)
            .at(orbit_track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        track["solar_alt"] = solar_altaz_data[0].degrees
        track["solar_az"] = solar_altaz_data[1].degrees
    if solar_beta:
        # append solar beta column
        # based on https://github.com/skyfielders/python-skyfield/issues/1054
        plane_normal = np.cross(
            np.array(orbit_track.position.m),
            np.array(orbit_track.velocity.m_per_s),
            axis=0,
        )
        sun = de421["earth"].at(orbit_track.t).observe(de421["sun"]).position.m  # type: ignore
        beta = np.pi / 2 - angle_between(plane_normal, sun)
        track["solar_beta"] = np.degrees(beta)
    # note: no further mask-based clip is needed here -- points outside the
    # (WGS84-only) mask were already dropped above, before points were
    # projected into the requested `coordinates` frame
    return track
