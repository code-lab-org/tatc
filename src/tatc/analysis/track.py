# -*- coding: utf-8 -*-
"""
Methods to generate coverage statistics.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Union, Optional
from datetime import datetime, timedelta
from enum import Enum

import pandas as pd
import numpy as np
import geopandas as gpd
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.functions import angle_between

from shapely.geometry import (
    Polygon,
    MultiPolygon,
    Point,
    LineString,
)
from shapely.ops import clip_by_rect, split
import pyproj

from ..schemas.instrument import PointedInstrument
from ..schemas.satellite import Satellite
from ..utils import (
    buffer_footprint,
    buffer_target,
    field_of_regard_to_swath_width,
)
from ..constants import de421, EARTH_MEAN_RADIUS, timescale


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
    times: List[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Optional[
        Union[
            Polygon,
            MultiPolygon,
            gpd.GeoDataFrame,
            gpd.GeoSeries,
        ]
    ] = None,
    coordinates: OrbitCoordinate = OrbitCoordinate.WGS84,
    orbit_output: OrbitOutput = OrbitOutput.POSITION,
    sat_sunlit: bool = False,
    solar_altaz: bool = False,
    solar_beta: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect orbit track points for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        times (typing.List[datetime.datetime]): The list of times to sample.
        instrument_index (int): The index of the observing instrument in satellite.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate swath width.
        mask (Polygon, MultiPolygon, geopandas.GeoDataFrame, geopandas.GeoSeries):
                An optional mask to constrain results.
        coordinates (OrbitCoordinate): The coordinate system of orbit track points.
        orbit_output (OrbitOutput): The output option.
        sat_sunlit (bool): `True` to include whether the satellite is sunlit.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles.
        solar_beta (bool): `True` to include solar beta angles.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected orbit track results.
    """
    if len(times) == 0:
        return _get_empty_orbit_track()
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    # propagate orbit
    orbit_track = satellite.orbit.to_tle().get_orbit_track(times)
    ssp = wgs84.geographic_position_of(orbit_track)
    if mask is not None:
        # trim orbit track to provided mask
        mask_contains_ssp = [
            (
                any(mask.contains(Point(longitude, latitude)))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.contains(Point(longitude, latitude))
            )
            for (longitude, latitude) in zip(
                ssp.longitude.degrees, ssp.latitude.degrees
            )
        ]
        if not any(mask_contains_ssp):
            return _get_empty_orbit_track()
        orbit_track = orbit_track[mask_contains_ssp]
        # recompute ssp
        ssp = wgs84.geographic_position_of(orbit_track)
    # create shapely points in proper coordinate system
    if coordinates == OrbitCoordinate.WGS84:
        points = [
            Point(longitude, latitude, elevation)
            for (longitude, latitude, elevation) in zip(
                ssp.longitude.degrees,
                ssp.latitude.degrees,
                ssp.elevation.m,
            )
        ]
    elif coordinates == OrbitCoordinate.ECEF:
        points = [
            Point(position[0], position[1], position[2])
            for position in ssp.itrs_xyz.m.T
        ]
    else:
        points = [
            Point(position[0], position[1], position[2])
            for position in orbit_track.xyz.m.T
        ]
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track)
    if len(orbit_track.t) == 1:
        # transform scalar to vector results
        valid_obs = np.array([valid_obs])
    # create velocity points if needed
    if orbit_output == OrbitOutput.POSITION:
        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": field_of_regard_to_swath_width(
                    ssp.elevation.m[i],
                    instrument.field_of_regard,
                    elevation,
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
            }
            for i, time in enumerate(orbit_track.t.utc_datetime())
        ]
    else:
        # compute satellite velocity
        if coordinates == OrbitCoordinate.ECI:
            velocities = [
                Point(velocity[0], velocity[1], velocity[2])
                for velocity in orbit_track.velocity.m_per_s.T
            ]
        else:
            velocities = [
                Point(velocity[0], velocity[1], velocity[2])
                for velocity in orbit_track.frame_xyz_and_velocity(itrs)[1].m_per_s.T
            ]

        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": field_of_regard_to_swath_width(
                    ssp.elevation.m[i],
                    instrument.field_of_regard,
                    elevation,
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
                "velocity": velocities[i],
            }
            for i, time in enumerate(orbit_track.t.utc_datetime())
        ]

    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if sat_sunlit:
        # append sat_sunlit column
        track["sat_sunlit"] = orbit_track.is_sunlit(de421)
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz = (
            (de421["earth"] + wgs84.geographic_position_of(orbit_track))
            .at(orbit_track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        track["solar_alt"] = solar_altaz[0].degrees
        track["solar_az"] = solar_altaz[1].degrees
    if solar_beta:
        # append solar beta column
        # based on https://github.com/skyfielders/python-skyfield/issues/1054
        plane_normal = np.cross(orbit_track.position.m, orbit_track.velocity.m_per_s, axis=0)
        sun = de421["earth"].at(orbit_track.t).observe(de421["sun"]).position.m
        beta = np.pi/2 - angle_between(plane_normal, sun)
        track["solar_beta"] = np.degrees(beta)
    if mask is not None:
        track = gpd.clip(track, mask).reset_index(drop=True)
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


def _get_utm_epsg_code(point: Point, swath_width: float) -> str:
    """
    Get the Universal Transverse Mercator (UTM) EPSG code for a ground track.

    Args:
        point (Point): the geodetic sub-satellite point
        swath_width (float): the observation swath width (meters)

    Returns:
        str: the EPSG code
    """
    # approximate footprint
    polygon = point.buffer(np.degrees(swath_width / 2 / EARTH_MEAN_RADIUS))
    results = pyproj.database.query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=pyproj.aoi.AreaOfInterest(
            *clip_by_rect(polygon, -180, -90, 180, 90).bounds
        ),
    )
    # return 5041 for UPS North; 5042 for UPS South; UTM zone, or 4087 for default
    return (
        "EPSG:5041"
        if polygon.bounds[3] > 84
        else (
            "EPSG:5042"
            if polygon.bounds[1] < -80
            else "EPSG:" + results[0].code if len(results) > 0 else "EPSG:4087"
        )
    )


def collect_ground_track(
    satellite: Satellite,
    times: List[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Optional[
        Union[
            Polygon,
            MultiPolygon,
            gpd.GeoDataFrame,
            gpd.GeoSeries,
        ]
    ] = None,
    crs: str = "EPSG:4087",
    sat_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect ground track polygons for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (Polygon, MultiPolygon, geopandas.GeoDataFrame, geopandas.GeoSeries):
                An optional mask to constrain results.
        crs (str): The coordinate reference system (CRS) in which to compute
                distance (default: World Equidistant Cylindrical `"EPSG:4087"`).
                Selecting `crs="utm"` uses Universal Transverse Mercator (UTM)
                zones for non-polar regions, and Universal Polar Stereographic
                (UPS) systems for polar regions. Selecting `crs="spice"` uses
                SPICE to compute observation footprints.
        sat_altaz (bool): `True` to include satellite altitude/azimuth angles.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results.
    """

    if len(times) == 0:
        return _get_empty_ground_track()
    # propagate orbit
    orbit_track = satellite.orbit.to_tle().get_orbit_track(times)
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if mask is not None and len(times) > 1:
        ssp = wgs84.geographic_position_of(orbit_track)
        buffered_mask = buffer_target(
            geometry=(
                mask
                if isinstance(mask, (Polygon, MultiPolygon))
                else (
                    mask.dissolve().iloc[0].geometry
                    if isinstance(mask, gpd.GeoDataFrame)
                    else mask.iloc[0]
                )
            ),
            altitude=satellite.orbit.to_tle().get_altitude(),
            inclination=satellite.orbit.to_tle().get_inclination(),
            field_of_regard=instrument.field_of_regard,
            time_step=np.diff(times).mean() / timedelta(seconds=1),
        )
        # cull orbit track with buffered mask
        buffered_mask_contains_ssp = [
            (
                any(buffered_mask.contains(Point(longitude, latitude)))
                if isinstance(buffered_mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else buffered_mask.contains(Point(longitude, latitude))
            )
            for (longitude, latitude) in zip(
                ssp.longitude.degrees, ssp.latitude.degrees
            )
        ]
        if not any(buffered_mask_contains_ssp):
            return _get_empty_ground_track()
        orbit_track = orbit_track[buffered_mask_contains_ssp]
        # compute footprint for culling
        footprint = instrument.compute_footprint(orbit_track, elevation=elevation)
        # cull orbit track to observable footprint
        mask_intersects_footprint = [
            (
                any(mask.intersects(f))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.intersects(f)
            )
            for f in (footprint if isinstance(footprint, list) else [footprint])
        ]
        if not any(mask_intersects_footprint):
            return _get_empty_ground_track()
        orbit_track = orbit_track[mask_intersects_footprint]
    # compute targets
    target = instrument.compute_footprint_center(orbit_track, elevation)
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track, target)
    if len(orbit_track.t) == 1:
        # transform scalar to vector results
        valid_obs = np.array([valid_obs])
    if crs == "spice":
        geometries = instrument.compute_footprint(
            orbit_track,
            None,
            elevation,
        )
        if len(orbit_track.t) == 1:
            geometries = [geometries]
    else:
        # compute the orbit track of the satellite
        gdf = collect_orbit_track(
            satellite, orbit_track.t.utc_datetime(), instrument_index, elevation
        )
        if crs == "utm":
            gdf["utm_crs"] = gdf.apply(
                lambda r: _get_utm_epsg_code(r.geometry, r.swath_width), axis=1
            )
            # preload transformers
            to_crs = {}
            from_crs = {}
            for code in gdf.utm_crs.unique():
                to_crs[code] = pyproj.Transformer.from_crs(
                    gdf.crs, pyproj.CRS(code), always_xy=True
                )
                from_crs[code] = pyproj.Transformer.from_crs(
                    pyproj.CRS(code), gdf.crs, always_xy=True
                )
            geometries = gdf.apply(
                lambda r: buffer_footprint(
                    r.geometry,
                    to_crs[r.utm_crs],
                    from_crs[r.utm_crs],
                    r.swath_width,
                    elevation,
                ),
                axis=1,
            ).values
        else:
            to_crs = pyproj.Transformer.from_crs(
                gdf.crs, pyproj.CRS(crs), always_xy=True
            )
            from_crs = pyproj.Transformer.from_crs(
                pyproj.CRS(crs), gdf.crs, always_xy=True
            )
            # construct polygons based on visible extent of instrument
            # project to specified elevation
            geometries = gdf.apply(
                lambda r: buffer_footprint(
                    r.geometry, to_crs, from_crs, r.swath_width, elevation
                ),
                axis=1,
            ).values
    records = [
        {
            "time": time,
            "satellite": satellite.name,
            "instrument": instrument.name,
            "valid_obs": valid_obs[i],
            "geometry": geometries[i],
        }
        for i, time in enumerate(orbit_track.t.utc_datetime())
    ]
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if sat_altaz:
        # append satellite altitude/azimuth columns
        sat_altaz = (orbit_track - target.at(orbit_track.t)).altaz()
        track["sat_alt"] = sat_altaz[0].degrees
        track["sat_az"] = sat_altaz[1].degrees
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz = (
            (de421["earth"] + target)
            .at(orbit_track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        track["solar_alt"] = solar_altaz[0].degrees
        track["solar_az"] = solar_altaz[1].degrees

    if mask is not None:
        track = gpd.clip(track, mask).reset_index(drop=True)
    return track


def compute_ground_track(
    satellite: Satellite,
    times: List[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Optional[
        Union[
            Polygon,
            MultiPolygon,
            gpd.GeoDataFrame,
            gpd.GeoSeries,
        ]
    ] = None,
    crs: str = "EPSG:4087",
    method: str = "point",
    dissolve_orbits: bool = True,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated ground track for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (Polygon, MultiPolygon, geopandas.GeoDataFrame, geopandas.GeoSeries):
                An optional mask to constrain results.
        crs (str): The coordinate reference system (CRS) in which to compute
                distance (default: World Equidistant Cylindrical `"EPSG:4087"`).
                Selecting `crs="utm"` uses Universal Transverse Mercator (UTM)
                zones for non-polar regions, and Universal Polar Stereographic
                (UPS) systems for polar regions. Selecting `crs="spice"` uses
                SPICE to compute footprints.
        method (str): The method for computing ground track: `"point"` buffers
                individual points and `"line"` buffers a line of points. The line
                method is not compatible with the `crs="spice"` option above.
        dissolve_orbits (bool): True, to aggregate multiple orbits in one output.
    Returns:
        GeoDataFrame: The data frame of aggregated ground track results.
    """
    if method not in ["point", "line"]:
        raise ValueError("Invalid method: " + str(method))
    if method == "point":
        track = collect_ground_track(
            satellite, times, instrument_index, elevation, mask, crs
        )
        if not track.empty:
            # assign orbit identifier
            track["orbit_id"] = [
                (time - times[0]) // satellite.orbit.to_tle().get_orbit_period()
                for time in track.time
            ]
            # filter to valid observations and dissolve
            track = (
                track[track.valid_obs].dissolve(by="orbit_id").reset_index(drop=True)
            )
        if dissolve_orbits:
            track = track.dissolve()
        return track
    if method == "line":
        if crs == "spice":
            raise ValueError("The line method is not compatible with spice")
        track = collect_orbit_track(satellite, times, instrument_index, elevation, None)
        # assign orbit identifier
        track["orbit_id"] = [
            (time - times[0]) // satellite.orbit.to_tle().get_orbit_period()
            for time in track.time
        ]
        # assign track identifiers to group contiguous observation periods
        track["track_id"] = (
            (track.valid_obs != track.valid_obs.shift()).astype("int").cumsum()
        )
        # filter to valid observations
        track = track[track.valid_obs].reset_index(drop=True)
        segments = []
        swath_widths = []
        for _, sub_track in track.groupby(["orbit_id", "track_id"]):
            # project points to specified elevation
            points = sub_track.geometry.apply(lambda p: Point(p.x, p.y, elevation))
            # extract longitudes
            lon = sub_track.geometry.apply(lambda p: p.x)
            # extract average swath width
            swath_widths.append(sub_track.swath_width.mean())
            # no anti-meridian crossings if all absolute longitude differences
            # are less than 180 deg
            if np.all(np.abs(np.diff(lon)) < 180):
                segments.append(LineString(points))
            else:
                # find anti-meridian crossings and calculate shift direction
                # coords from W -> E (shift < 0) will add 360 degrees to E component
                # coords from E -> W (shift > 0) will subtract 360 degrees from W component
                shift = np.insert(np.cumsum(np.around(np.diff(lon) / 360)), 0, 0)
                points = [
                    Point(p.x - 360 * shift[i], p.y, p.z) for i, p in enumerate(points)
                ]
                # split along the anti-meridian (-180 for shift > 0; 180 for shift < 0)
                shift_dir = -180 if shift.max() >= 1 else 180
                collection = split(
                    LineString(points),
                    LineString([(shift_dir, -180), (shift_dir, 180)]),
                )
                # map longitudes from (-540, -180] to (-180, 180] and from [180, 540) to [-180, 180)
                segments.extend(
                    [
                        (
                            LineString([(c[0] + 360, c[1], c[2]) for c in line.coords])
                            if np.all([c[0] <= -180 for c in line.coords])
                            else (
                                LineString(
                                    [(c[0] - 360, c[1], c[2]) for c in line.coords]
                                )
                                if np.all([c[0] >= 180 for c in line.coords])
                                else LineString(
                                    [(c[0], c[1], c[2]) for c in line.coords]
                                )
                            )
                        )
                        for line in collection.geoms
                    ]
                )
        to_crs = pyproj.Transformer.from_crs(track.crs, pyproj.CRS(crs), always_xy=True)
        from_crs = pyproj.Transformer.from_crs(
            pyproj.CRS(crs), track.crs, always_xy=True
        )
        polygons = [
            buffer_footprint(segment, to_crs, from_crs, swath_width, elevation)
            for segment, swath_width in zip(segments, swath_widths)
        ]
        # dissolve the original track
        track = track.dissolve(by=["orbit_id", "track_id"]).reset_index(drop=True)
        # and replace the geometry with the union of computed polygons
        track.geometry = polygons
        if mask is not None:
            track = gpd.clip(track, mask).reset_index(drop=True)
        if dissolve_orbits:
            track = track.dissolve()
        return track


def collect_ground_pixels(
    satellite: Satellite,
    times: List[datetime],
    instrument_index: int = 0,
    elevation: float = 0,
    mask: Optional[
        Union[
            Polygon,
            MultiPolygon,
            gpd.GeoDataFrame,
            gpd.GeoSeries,
        ]
    ] = None,
    sat_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect ground pixels for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument_index (int): The index of the observing instrument in satellite.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground pixels.
        mask (Polygon, MultiPolygon, geopandas.GeoDataFrame, geopandas.GeoSeries):
                An optional mask to constrain results.
        sat_altaz (bool): `True` to include satellite altitude/azimuth angles.
        solar_altaz (bool): `True` to include solar altitude/azimuth angles.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground pixels results.
    """

    if len(times) == 0:
        return _get_empty_ground_track()
    # propagate orbit
    orbit_track = satellite.orbit.to_tle().get_orbit_track(times)
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )
    if mask is not None and len(times) > 1:
        ssp = wgs84.geographic_position_of(orbit_track)
        buffered_mask = buffer_target(
            geometry=(
                mask
                if isinstance(mask, (Polygon, MultiPolygon))
                else (
                    mask.dissolve().iloc[0].geometry
                    if isinstance(mask, gpd.GeoDataFrame)
                    else mask.iloc[0]
                )
            ),
            altitude=satellite.orbit.to_tle().get_altitude(),
            inclination=satellite.orbit.to_tle().get_inclination(),
            field_of_regard=instrument.field_of_regard,
            time_step=np.diff(times).mean() / timedelta(seconds=1),
        )
        # cull orbit track with buffered mask
        buffered_mask_contains_ssp = [
            (
                any(buffered_mask.contains(Point(longitude, latitude)))
                if isinstance(buffered_mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else buffered_mask.contains(Point(longitude, latitude))
            )
            for (longitude, latitude) in zip(
                ssp.longitude.degrees, ssp.latitude.degrees
            )
        ]
        if not any(buffered_mask_contains_ssp):
            return _get_empty_ground_track()
        orbit_track = orbit_track[buffered_mask_contains_ssp]
        # compute footprint for culling
        footprint = instrument.compute_footprint(orbit_track, elevation=elevation)
        # cull orbit track to observable footprint
        mask_intersects_footprint = [
            (
                any(mask.intersects(f))
                if isinstance(mask, (gpd.GeoDataFrame, gpd.GeoSeries))
                else mask.intersects(f)
            )
            for f in (footprint if isinstance(footprint, list) else [footprint])
        ]
        if not any(mask_intersects_footprint):
            return _get_empty_ground_track()
        orbit_track = orbit_track[mask_intersects_footprint]
    # compute the footprint pixel array
    geometries = instrument.compute_footprint_pixel_array(
        orbit_track,
        elevation,
    )
    # construct results as a list of dictionaries
    records = [
        {
            "time": time,
            "satellite": satellite.name,
            "instrument": instrument.name,
            "valid_obs": instrument.is_valid_observation(
                orbit_track[i], wgs84.latlon(point.y, point.x, point.z)
            ),
            "geometry": point,
        }
        for i, time in enumerate(orbit_track.t.utc_datetime())
        for point in (
            geometries[i].geoms if len(orbit_track.t) > 1 else geometries.geoms
        )
    ]
    # build geodataframe
    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if sat_altaz:
        # append satellite altitude/azimuth columns
        sat_altaz = [
            (
                orbit_track[list(orbit_track.t.utc_datetime()).index(record["time"])]
                - wgs84.latlon(
                    record["geometry"].y, record["geometry"].x, record["geometry"].z
                ).at(timescale.from_datetime(record["time"]))
            ).altaz()
            for record in records
        ]
        gdf["sat_alt"] = list(map(lambda altaz: altaz[0].degrees, sat_altaz))
        gdf["sat_az"] = list(map(lambda altaz: altaz[1].degrees, sat_altaz))
    if solar_altaz:
        # append solar altitude/azimuth columns
        solar_altaz = [
            (
                de421["earth"]
                + wgs84.latlon(
                    record["geometry"].y, record["geometry"].x, record["geometry"].z
                )
            )
            .at(timescale.from_datetime(record["time"]))
            .observe(de421["sun"])
            .apparent()
            .altaz()
            for record in records
        ]
        gdf["solar_alt"] = list(map(lambda altaz: altaz[0].degrees, solar_altaz))
        gdf["solar_az"] = list(map(lambda altaz: altaz[1].degrees, solar_altaz))

    if mask is not None:
        gdf = gpd.clip(gdf, mask).reset_index(drop=True)
    return gdf
