"""
Methods to generate coverage statistics.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
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
from ..schemas import PointedInstrument, Satellite
from ..utils.observation import field_of_regard_to_swath_width
from ..utils.projection import buffer_target


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
                "swath_width": field_of_regard_to_swath_width(
                    np.array(sat_pos.elevation.m)[i],
                    instrument.field_of_regard,
                    elevation,
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
                "swath_width": field_of_regard_to_swath_width(
                    np.array(sat_pos.elevation.m)[i],
                    instrument.field_of_regard,
                    elevation,
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


def collect_ground_track(
    satellite: Satellite,
    times: list[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    sat_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's viewable ground footprint at each requested
    time, projected to a specified elevation, using SPICE to compute the
    exact ray/WGS-84-geoid intersection (see
    `tatc.utils.projection.compute_footprint`).

    Args:
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
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
                for the sub-satellite point.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results.
    """

    if len(times) == 0:
        return _get_empty_ground_track()
    # propagate orbit
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if mask is not None and len(times) > 1:
        # use the (possibly repeat-cycle-corrected) geodetic position for this
        # rough, buffer-tolerant culling step only; final geometry/validity
        # below still uses the true orbit_track for consistency.
        sat_pos = satellite.orbit.to_gp_orbit().get_geographic_position(times)
        if isinstance(mask, (Polygon, MultiPolygon)):
            geometry = mask
        elif isinstance(mask, gpd.GeoDataFrame):
            geometry = mask.dissolve().iloc[0].geometry
        else:
            geometry = mask.union_all()
        buffered_mask = buffer_target(
            geometry=geometry,
            altitude=satellite.orbit.to_gp_orbit().get_mean_altitude(),
            inclination=satellite.orbit.to_gp_orbit().get_inclination(),
            field_of_regard=instrument.field_of_regard,
            time_step=np.diff(np.array(times)).mean() / timedelta(seconds=1),
        )
        # cull orbit track with buffered mask
        buffered_mask_contains_sat_pos = [
            buffered_mask.contains(Point(longitude, latitude))
            for (longitude, latitude) in zip(
                sat_pos.longitude.degrees, sat_pos.latitude.degrees
            )
        ]
        if not any(buffered_mask_contains_sat_pos):
            return _get_empty_ground_track()
        orbit_track = orbit_track[buffered_mask_contains_sat_pos]
        # compute footprint for culling
        footprint = instrument.compute_footprint(orbit_track, elevation=elevation)
        # cull orbit track to observable footprint
        mask_intersects_footprint = [
            (
                any(mask.intersects(f))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.intersects(f)
            )
            for f in footprint
        ]
        if not any(mask_intersects_footprint):
            return _get_empty_ground_track()
        orbit_track = orbit_track[mask_intersects_footprint]
    # compute targets
    target = instrument.compute_footprint_center(orbit_track, elevation)
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track, target)
    # compute footprints via SPICE (exact ray/WGS-84-geoid intersection)
    geometries = instrument.compute_footprint(
        orbit_track,
        None,
        elevation,
    )
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

    if mask is not None:
        track = gpd.clip(track, mask).reset_index(drop=True)
    return track


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
    track = collect_ground_track(satellite, times, instrument_index, elevation, mask)
    if not track.empty:
        # assign orbit identifier
        track["orbit_id"] = [
            (time - times[0])
            // satellite.orbit.to_gp_orbit()
            .get_closest_element(time)
            .get_orbit_period()
            for time in track.time
        ]
        # filter to valid observations and dissolve
        track = track[track.valid_obs].dissolve(by="orbit_id").reset_index(drop=True)
    if dissolve_orbits:
        track = track.dissolve()
    return track


def collect_ground_pixels(
    satellite: Satellite,
    times: list[datetime],
    instrument_index: int = 0,
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
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground pixels.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results. Providing a mask limits
                the propagated orbit to the (buffered) region of interest,
                improving performance for small areas.
        sat_altaz (bool): `True` to include satellite altitude/azimuth angles.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground pixels results.
    """

    if len(times) == 0:
        return _get_empty_ground_track()
    # propagate orbit
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )
    if mask is not None and len(times) > 1:
        # use the (possibly repeat-cycle-corrected) geodetic position for this
        # rough, buffer-tolerant culling step only; final geometry/validity
        # below still uses the true orbit_track for consistency.
        sat_pos = satellite.orbit.to_gp_orbit().get_geographic_position(times)
        if isinstance(mask, (Polygon, MultiPolygon)):
            geometry = mask
        elif isinstance(mask, gpd.GeoDataFrame):
            geometry = mask.dissolve().iloc[0].geometry
        else:
            geometry = mask.union_all()
        buffered_mask = buffer_target(
            geometry=geometry,
            altitude=satellite.orbit.to_gp_orbit().get_mean_altitude(),
            inclination=satellite.orbit.to_gp_orbit().get_inclination(),
            field_of_regard=instrument.field_of_regard,
            time_step=np.diff(np.array(times)).mean() / timedelta(seconds=1),
        )
        # cull orbit track with buffered mask
        buffered_mask_contains_sat_pos = [
            buffered_mask.contains(Point(longitude, latitude))
            for (longitude, latitude) in zip(
                sat_pos.longitude.degrees, sat_pos.latitude.degrees
            )
        ]
        if not any(buffered_mask_contains_sat_pos):
            return _get_empty_ground_track()
        orbit_track = orbit_track[buffered_mask_contains_sat_pos]
        # compute footprint for culling
        footprint = instrument.compute_footprint(orbit_track, elevation=elevation)
        # cull orbit track to observable footprint
        mask_intersects_footprint = [
            (
                any(mask.intersects(f))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.intersects(f)
            )
            for f in footprint
        ]
        if not any(mask_intersects_footprint):
            return _get_empty_ground_track()
        orbit_track = orbit_track[mask_intersects_footprint]
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

    if mask is not None:
        gdf = gpd.clip(gdf, mask).reset_index(drop=True)
    return gdf
