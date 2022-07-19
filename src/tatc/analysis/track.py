# -*- coding: utf-8 -*-
"""
Methods to generate coverage statistics.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from typing import List, Union, Optional
import pandas as pd
import numpy as np
import geopandas as gpd
from datetime import datetime, timedelta
from skyfield.api import load, wgs84, EarthSatellite
from shapely.geometry import (
    Polygon,
    MultiPolygon,
    MultiPoint,
    Point,
    MultiLineString,
    LineString,
)
from shapely.ops import transform, clip_by_rect, unary_union
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


def collect_orbit_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
) -> gpd.GeoDataFrame:
    """
    Collect orbit track points for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        times (typing.List[datetime.datetime]): The list of times to sample.
        mask (Polygon or MultiPolygon): An optional mask to constrain results.

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
    # project to geographic positions
    subpoints = [wgs84.geographic_position_of(position) for position in positions]
    # create shapely points
    points = [
        Point(
            subpoint.longitude.degrees, subpoint.latitude.degrees, subpoint.elevation.m
        )
        for subpoint in subpoints
    ]
    valid_obs = instrument.is_valid_observation(sat, ts_times)

    records = [
        {
            "time": time,
            "satellite": satellite.name,
            "instrument": instrument.name,
            "swath_width": field_of_regard_to_swath_width(
                subpoints[i].elevation.m,
                instrument.field_of_regard,
            ),
            "valid_obs": valid_obs[i],
            "geometry": points[i],
        }
        for i, time in enumerate(times)
    ]
    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)


def collect_ground_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    crs: str = "EPSG:4087",
) -> gpd.GeoDataFrame:
    """
    Collect ground track polygons for a satellite of interest.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
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
    gdf = collect_orbit_track(satellite, instrument, times, mask)
    if gdf.empty:
        return gdf
    # project points to zero elevation
    gdf.geometry = gdf.geometry.apply(lambda p: Point(p.x, p.y))
    # at each point, draw a buffer equivalent to the swath radius
    if crs == "utm":
        # do the swath projection in the matching utm zone
        def _get_utm_epsg_code(p):
            results = pyproj.database.query_utm_crs_info(
                datum_name="WGS 84",
                area_of_interest=pyproj.aoi.AreaOfInterest(p.x, p.y, p.x, p.y),
            )
            # return the first code if exists; otherwise return a default value
            # 5041 is UPS North; 5042 is UPS South; 4087 is the default
            return (
                "EPSG:" + results[0].code
                if len(results) > 0
                else "EPSG:5041"
                if p.y > 84
                else "EPSG:5042"
                if p.y < -80
                else "EPSG:4087"
            )

        utm_crs = gdf.geometry.apply(_get_utm_epsg_code)
        for code in utm_crs.unique():
            to_crs = pyproj.Transformer.from_crs(
                gdf.crs, pyproj.CRS(code), always_xy=True
            ).transform
            from_crs = pyproj.Transformer.from_crs(
                pyproj.CRS(code), gdf.crs, always_xy=True
            ).transform
            gdf.loc[utm_crs == code, "geometry"] = gdf[utm_crs == code].apply(
                lambda r: transform(
                    from_crs,
                    transform(to_crs, r.geometry).buffer(r.swath_width / 2),
                ),
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
    # split polygons to wrap over the anti-meridian and poles
    gdf.geometry = gdf.apply(lambda r: split_polygon(r.geometry), axis=1)

    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)


def compute_ground_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
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
        mask (Polygon or MultiPolygon): An optional mask to constrain results.
        crs (str): The coordinate reference system (CRS) in which to compute
                distance (default: World Equidistant Cylindrical `"EPSG:4087"`).
                Selecting `crs="utm"` uses Universal Transverse Mercator (UTM)
                zones for non-polar regions, and Universal Polar Stereographic
                (UPS) systems for polar regions.
        method (str): The method for computing ground track: `"point"` buffers
                individual points while `"line"` buffers a line of points.

    Returns:
        GeoDataFrame: The data frame of aggregated ground track results.
    """
    if method == "point":
        track = collect_ground_track(satellite, instrument, times, mask, crs)
        # filter to valid observations and dissolve
        return track[track.valid_obs].dissolve()
    elif method == "line":
        track = collect_orbit_track(satellite, instrument, times, mask)
        # filter to valid observations
        track = track[track.valid_obs].reset_index(drop=True)
        # project points to zero elevation
        points = track.geometry.apply(lambda p: Point(p.x, p.y))
        # extract longitudes
        lons = track.geometry.apply(lambda p: p.x)
        if np.any(np.diff(np.sign(lons))):
            # split lines when crossing meridian or anti-meridian
            points_list = np.split(
                points,
                1 + np.where(np.diff(np.sign(lons)))[0],
            )
            segments = [
                LineString(pts.tolist()) if len(pts.index) > 1 else pts[0]
                for pts in points_list
            ]
        else:
            segments = [LineString(points)]
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
                transform(to_crs, segment).buffer(track.swath_width.mean() / 2),
            )
            for segment in segments
        ]
        # split polygons if necessary
        polygons = list(map(split_polygon, polygons))
        # dissolve the original track
        track = track.dissolve()
        # and replace the geometry with the union of computed polygons
        track.geometry = [unary_union(polygons)]
        return track
