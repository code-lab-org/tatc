# -*- coding: utf-8 -*-
"""
Methods to generate coverage statistics.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Union, Optional
from datetime import datetime
from enum import Enum

import pandas as pd
import numpy as np
import geopandas as gpd
from skyfield.api import wgs84, EarthSatellite
from skyfield.framelib import itrs

from shapely.geometry import (
    Polygon,
    MultiPolygon,
    Point,
    LineString,
)
from shapely.ops import split, transform, unary_union
import pyproj

from ..schemas.satellite import Satellite
from ..schemas.instrument import Instrument
from ..utils import (
    split_polygon,
    field_of_regard_to_swath_width,
)
from ..constants import timescale


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
    instrument: Instrument,
    times: List[datetime],
    elevation: float = 0,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    coordinates: OrbitCoordinate = OrbitCoordinate.WGS84,
    orbit_output: OrbitOutput = OrbitOutput.POSITION,
) -> gpd.GeoDataFrame:
    """
    Collect orbit track points for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        times (typing.List[datetime.datetime]): The list of times to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate swath width.
        mask (Polygon or MultiPolygon): An optional mask to constrain results.
        coordinates (OrbitCoordinate): The coordinate system of orbit track points.
        orbit_output (OrbitOutput): The output option.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected orbit track results.
    """
    if len(times) == 0:
        return _get_empty_orbit_track()
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # create skyfield time series
    ts_times = timescale.from_datetimes(times)
    # compute satellite positions
    positions = sat.at(ts_times)
    geo_positions = [wgs84.geographic_position_of(position) for position in positions]
    # create shapely points in proper coordinate system
    if coordinates == OrbitCoordinate.WGS84:
        points = [
            Point(
                position.longitude.degrees,
                position.latitude.degrees,
                position.elevation.m,
            )
            for position in geo_positions
        ]
    elif coordinates == OrbitCoordinate.ECEF:
        points = [
            Point(
                position.itrs_xyz.m[0], position.itrs_xyz.m[1], position.itrs_xyz.m[2]
            )
            for position in geo_positions
        ]
    else:
        points = [
            Point(position.xyz.m[0], position.xyz.m[1], position.xyz.m[2])
            for position in positions
        ]
    # determine observation validity
    valid_obs = instrument.is_valid_observation(sat, ts_times)
    # create velocity points if needed
    if orbit_output == OrbitOutput.POSITION:
        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": field_of_regard_to_swath_width(
                    geo_positions[i].elevation.m,
                    instrument.field_of_regard,
                    elevation,
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
            }
            for i, time in enumerate(times)
        ]
    else:
        # compute satellite velocity
        if coordinates == OrbitCoordinate.ECI:
            eci_velocity = positions.velocity.m_per_s
            velocities = [
                Point(eci_velocity[0][i], eci_velocity[1][i], eci_velocity[2][i])
                for i in range(len(eci_velocity[0]))
            ]
        else:
            velocity = positions.frame_xyz_and_velocity(itrs)[1].m_per_s
            velocities = [
                Point(velocity[0][i], velocity[1][i], velocity[2][i])
                for i in range(len(velocity[0]))
            ]

        records = [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": field_of_regard_to_swath_width(
                    geo_positions[i].elevation.m,
                    instrument.field_of_regard,
                    elevation,
                ),
                "valid_obs": valid_obs[i],
                "geometry": points[i],
                "velocity": velocities[i],
            }
            for i, time in enumerate(times)
        ]

    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)


def _get_utm_epsg_code(point: Point) -> str:
    """
    Get the Universal Transverse Mercator (UTM) EPSG code for a geodetic point.

    Args:
        point (Point): the geodetic point

    Returns:
        str: the UTM EPSG code
    """
    results = pyproj.database.query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=pyproj.aoi.AreaOfInterest(point.x, point.y, point.x, point.y),
    )
    # return the first code if exists; otherwise return a default value
    # 5041 is UPS North; 5042 is UPS South; 4087 is the default
    return (
        "EPSG:" + results[0].code
        if len(results) > 0
        else (
            "EPSG:5041"
            if point.y > 84
            else "EPSG:5042" if point.y < -80 else "EPSG:4087"
        )
    )


def collect_ground_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    elevation: float = 0,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    crs: str = "EPSG:4087",
) -> gpd.GeoDataFrame:
    """
    Collect ground track polygons for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (Polygon or MultiPolygon): An optional mask to constrain results.
        crs (str): The coordinate reference system (CRS) in which to compute
                distance (default: World Equidistant Cylindrical `"EPSG:4087"`).
                Selecting `crs="utm"` uses Universal Transverse Mercator (UTM)
                zones for non-polar regions, and Universal Polar Stereographic
                (UPS) systems for polar regions.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results.
    """
    # first, compute the orbit track of the satellite
    gdf = collect_orbit_track(satellite, instrument, times, elevation, mask)
    if gdf.empty:
        return gdf
    # project points to specified elevation
    gdf.geometry = gdf.geometry.apply(lambda p: Point(p.x, p.y, elevation))
    # at each point, draw a buffer equivalent to the swath radius
    if crs == "utm":
        # do the swath projection in the matching utm zone
        utm_crs = gdf.geometry.apply(_get_utm_epsg_code)
        for code in utm_crs.unique():
            to_crs = pyproj.Transformer.from_crs(
                gdf.crs, pyproj.CRS(code), always_xy=True
            ).transform
            from_crs = pyproj.Transformer.from_crs(
                pyproj.CRS(code), gdf.crs, always_xy=True
            ).transform
            if code in ("EPSG:5041", "EPSG:5042"):
                # keep polygons 500km away from UPS poles and apply zero buffer
                # to mitigate invalid polygons near poles in EPSG:4087
                gdf.loc[utm_crs == code, "geometry"] = gdf[utm_crs == code].apply(
                    lambda r: transform(
                        from_crs,
                        transform(to_crs, r.geometry)
                        .buffer(r.swath_width / 2)
                        .difference(Point(2e6, 2e6, 0).buffer(5e5)),
                    ).buffer(0),
                    axis=1,
                )
            else:
                # apply zero buffer to mitigate invalid polygons in EPSG:4087
                gdf.loc[utm_crs == code, "geometry"] = gdf[utm_crs == code].apply(
                    lambda r: transform(
                        from_crs,
                        transform(to_crs, r.geometry).buffer(r.swath_width / 2),
                    ).buffer(0),
                    axis=1,
                )
    else:
        # do the swath projection in the specified coordinate reference system
        to_crs = pyproj.Transformer.from_crs(
            gdf.crs, pyproj.CRS(crs), always_xy=True
        ).transform
        from_crs = pyproj.Transformer.from_crs(
            pyproj.CRS(crs), gdf.crs, always_xy=True
        ).transform
        gdf.geometry = gdf.apply(
            lambda r: transform(
                from_crs,
                transform(to_crs, r.geometry).buffer(r.swath_width / 2),
            ),
            axis=1,
        )
    # add elevation to all polygon coordinates (otherwise lost during buffering)
    gdf.geometry = gdf.geometry.apply(
        lambda g: (
            Polygon(
                [(p[0], p[1], elevation) for p in g.exterior.coords],
                [[(p[0], p[1], elevation) for p in i.coords] for i in g.interiors],
            )
            if isinstance(g, Polygon)
            else MultiPolygon(
                [
                    Polygon(
                        [(p[0], p[1], elevation) for p in n.exterior.coords],
                        [
                            [(p[0], p[1], elevation) for p in i.coords]
                            for i in n.interiors
                        ],
                    )
                    for n in g.geoms
                ]
            )
        )
    )
    # split polygons to wrap over the anti-meridian and poles
    gdf.geometry = gdf.apply(lambda r: split_polygon(r.geometry), axis=1)

    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)


def compute_ground_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    elevation: float = 0,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    crs: str = "EPSG:4087",
    method: str = "point",
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated ground track for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        elevation (float): The elevation (meters) above the datum in the
                WGS 84 coordinate system for which to calculate ground track.
        mask (Polygon or MultiPolygon): An optional mask to constrain results.
        crs (str): The coordinate reference system (CRS) in which to compute
                distance (default: World Equidistant Cylindrical `"EPSG:4087"`).
                Selecting `crs="utm"` uses Universal Transverse Mercator (UTM)
                zones for non-polar regions, and Universal Polar Stereographic
                (UPS) systems for polar regions.
        method (str): The method for computing ground track: `"point"` buffers
                individual points while `"line"` buffers a line of points. Note
                that the `"line"` method assumes continuous observations and only
                supports ground tracks that span LESS than 360 degrees longitude.

    Returns:
        GeoDataFrame: The data frame of aggregated ground track results.
    """
    if method == "point":
        track = collect_ground_track(satellite, instrument, times, elevation, mask, crs)
        # filter to valid observations and dissolve
        return track[track.valid_obs].dissolve()
    if method == "line":
        track = collect_orbit_track(satellite, instrument, times, elevation, mask)
        # assign track identifiers to group contiguous observation periods
        track["track_id"] = (
            (track.valid_obs != track.valid_obs.shift()).astype("int").cumsum()
        )
        # filter to valid observations
        track = track[track.valid_obs].reset_index(drop=True)
        segments = []
        swath_widths = []
        for track_id in track.track_id.unique():
            # project points to specified elevation
            points = track[track.track_id == track_id].geometry.apply(
                lambda p: Point(p.x, p.y, elevation)
            )
            # extract longitudes and latitudes
            lon = track[track.track_id == track_id].geometry.apply(lambda p: p.x)
            # extract average swath width
            swath_widths.append(track[track.track_id == track_id].swath_width.mean())
            # no anti-meridian crossings if all absolute longitude differences
            # are less than 180 deg for non-polar points (mean absolute latitude < 80 deg)
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
        to_crs = pyproj.Transformer.from_crs(
            track.crs, pyproj.CRS(crs), always_xy=True
        ).transform
        from_crs = pyproj.Transformer.from_crs(
            pyproj.CRS(crs), track.crs, always_xy=True
        ).transform
        # apply the coordinate reference system transformation
        polygons = [
            transform(
                from_crs,
                transform(to_crs, segment).buffer(swath_width / 2),
            )
            for segment, swath_width in zip(segments, swath_widths)
        ]
        # add elevation to all polygon coordinates (otherwise lost during buffering)
        polygons = [
            Polygon(
                [(p[0], p[1], elevation) for p in g.exterior.coords],
                [[(p[0], p[1], elevation) for p in i.coords] for i in g.interiors],
            )
            for g in polygons
        ]
        # split polygons if necessary
        polygons = list(map(split_polygon, polygons))
        # dissolve the original track
        track = track.dissolve()
        # and replace the geometry with the union of computed polygons
        track.geometry = [unary_union(polygons)]
        return track
    raise ValueError("Invalid method: " + str(method))
