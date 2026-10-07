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
from shapely.geometry import (
    MultiPolygon,
    Point,
    Polygon,
)
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.functions import angle_between

from ..constants import de421
from ..schemas import AllInstruments, ConicalInstrument, Satellite
from ..utils.observation import field_of_regard_to_swath_width
from .validation import _check_satellite


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
    satellite: Satellite,
    times: list[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    coordinates: OrbitCoordinate = OrbitCoordinate.WGS84,
    orbit_output: OrbitOutput = OrbitOutput.POSITION,
    sat_sunlit: bool = False,
    solar_altaz: bool = False,
    solar_beta: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the satellite's own position (and, optionally, velocity) at each
    requested time, in a specified coordinate frame. Note this reports the
    satellite's location, not the zero-elevation ground point beneath it; use
    `collect_ground_track` for footprint/ground-projected results.

    Args:
        satellite (Satellite): The observing satellite.
        times (typing.List[datetime.datetime]): The list of times to sample.
        instrument_index (int): The index of the observing instrument in satellite.
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
        geopandas.GeoDataFrame: The data frame of collected orbit track results.
    """
    _check_satellite(satellite)
    if len(times) == 0:
        return _get_empty_orbit_track()
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    # propagate orbit
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
    # geodetic (WGS84) position of the satellite itself (not the zero-elevation
    # subpoint below it -- elevation here is the satellite's own altitude)
    sat_pos = wgs84.geographic_position_of(orbit_track)
    if mask is not None:
        # trim orbit track to provided mask; mask is always interpreted in
        # WGS84 (lon/lat) coordinates regardless of the requested output
        # `coordinates`, so this filter applies consistently to all of them
        mask_contains_sat_pos = [
            (
                any(mask.contains(Point(longitude, latitude)))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.contains(Point(longitude, latitude))
            )
            for (longitude, latitude) in zip(
                np.array(sat_pos.longitude.degrees), np.array(sat_pos.latitude.degrees)
            )
        ]
        if not any(mask_contains_sat_pos):
            return _get_empty_orbit_track()
        orbit_track = orbit_track[mask_contains_sat_pos]
        # recompute sat_pos
        sat_pos = wgs84.geographic_position_of(orbit_track)
    # create shapely points in proper coordinate system
    if coordinates == OrbitCoordinate.WGS84:
        points = [
            Point(longitude, latitude, elevation)
            for (longitude, latitude, elevation) in zip(
                np.array(sat_pos.longitude.degrees),
                np.array(sat_pos.latitude.degrees),
                np.array(sat_pos.elevation.m),
            )
        ]
    elif coordinates == OrbitCoordinate.ECEF:
        points = [
            Point(position[0], position[1], position[2])
            for position in np.array(sat_pos.itrs_xyz.m).T
        ]
    else:
        points = [
            Point(position[0], position[1], position[2])
            for position in np.array(orbit_track.xyz.m).T
        ]
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track)
    # create velocity points if needed
    if orbit_output == OrbitOutput.POSITION:
        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": _swath_width(
                    instrument, np.array(sat_pos.elevation.m)[i], elevation
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
            }
            for i, time in enumerate(orbit_track.t.utc_datetime())  # type: ignore
        ]
    else:
        # compute satellite velocity
        if coordinates == OrbitCoordinate.ECI:
            velocities = [
                Point(velocity[0], velocity[1], velocity[2])
                for velocity in np.array(orbit_track.velocity.m_per_s).T
            ]
        elif coordinates == OrbitCoordinate.ECEF:
            velocities = [
                Point(velocity[0], velocity[1], velocity[2])
                for velocity in np.array(
                    orbit_track.frame_xyz_and_velocity(itrs)[1].m_per_s
                ).T
            ]
        else:
            # rotate ECEF velocity into local East/North/Up components at the
            # satellite's geodetic longitude/latitude
            ecef_velocity = np.array(
                orbit_track.frame_xyz_and_velocity(itrs)[1].m_per_s
            )
            lon = np.radians(np.array(sat_pos.longitude.degrees))
            lat = np.radians(np.array(sat_pos.latitude.degrees))
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
            velocities = [Point(e, n, u) for e, n, u in zip(east, north, up)]

        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": _swath_width(
                    instrument, np.array(sat_pos.elevation.m)[i], elevation
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
                "velocity": velocities[i],
            }
            for i, time in enumerate(orbit_track.t.utc_datetime())  # type: ignore
        ]

    # tag the CRS to match the requested output coordinates: WGS84 is
    # geographic degrees, ECEF is geocentric meters, and ECI (GCRS) is an
    # inertial, time-varying frame with no fixed EPSG code
    track_crs = (
        "EPSG:4326"
        if coordinates == OrbitCoordinate.WGS84
        else "EPSG:4978" if coordinates == OrbitCoordinate.ECEF else None
    )
    track = gpd.GeoDataFrame(records, crs=track_crs)
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
