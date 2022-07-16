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
    Get an empty ground track data frame.
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

    :param satellite: The observing satellite
    :type satellite: :class:`tatc.schemas.satellite.Satellite`
    :param instrument: The instrument used to make observations
    :type instrument: :class:`tatc.schemas.instrument.Instrument`
    :param times: The list of times to sample.
    :type times: list
    :param mask: A mask to constrain ground track.
    :type mask: :class:`shapely.Polygon` or :class:`shapely.MultiPolygon`, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded points.
    :rtype: :class:`geopandas.GeodataFrame`
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
    resolution: int = 4,
    split_polygons: bool = True,
) -> gpd.GeoDataFrame:
    """
    Model the ground track swath for a satellite of interest.

    :param satellite: The observing satellite
    :type satellite: :class:`tatc.schemas.satellite.Satellite`
    :param instrument: The observing instrument
    :type instrument: :class:`tatc.schemas.instrument.Instrument`
    :param times: The list of times to sample
    :type times: list
    :param mask: A mask to constrain ground track.
    :type mask: :class:`shapely.Polygon` or :class:`shapely.MultiPolygon`, optional
    :param crs: Coordinate reference system for projecting swath distance from
            sub-satellite points (default: "EPSG:4087"). Note that "utm" will
            use the corresponding Unified Transverse Mercator (UTM) zone which
            is accurate but slow.
    :type crs: str, optional
    :param resolution: Shapely buffer resolution of the projected swath
    :type resolution: int, optional
    :param split_polygons: True, if Polygons should be split into MultiPolygons
            when crossing the anti-meridian or exceeding polar boundaries.
    :type split_polygons: bool, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded polygons.
    :rtype: :class:`geopandas.GeoDataFrame`
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
            return results[0].code if len(results) > 0 else "4087"

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
                    transform(to_crs, r.geometry).buffer(
                        r.swath_width / 2, resolution=resolution
                    ),
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
                transform(to_crs, r.geometry).buffer(
                    r.swath_width / 2, resolution=resolution
                ),
            ),
            axis=1,
        )

    if split_polygons:
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
    resolution: int = 4,
    split_polygons: bool = True,
    method: str = "point",
    valid_obs: Optional[bool] = None,
) -> gpd.GeoDataFrame:
    """
    Compute the ground track swath for a satellite of interest.

    :param satellite: The observing satellite
    :type satellite: :class:`tatc.schemas.satellite.Satellite`
    :param instrument: The observing instrument
    :type instrument: :class:`tatc.schemas.instrument.Instrument`
    :param times: The list of times to sample
    :type times: list
    :param mask: A mask to constrain ground track.
    :type mask: :class:`shapely.Polygon` or :class:`shapely.MultiPolygon`, optional
    :param crs: Coordinate reference system for projecting swath distance from
            sub-satellite points (default: "EPSG:4087"). Note that "utm" will
            use the corresponding Unified Transverse Mercator (UTM) zone which
            is accurate but slow.
    :type crs: str, optional
    :param resolution: Shapely buffer resolution of the projected swath
    :type resolution: int, optional
    :param split_polygons: True, if Polygons should be split into MultiPolygons
            when crossing the anti-meridian or exceeding polar boundaries.
    :type split_polygons: bool, optional
    :param method: Method for computing the ground track: "point" buffers
            individual points while "line" buffers a line string.
    :type method: str (default: point)
    :param valid_obs: True (False), if ground track is restricted to valid (invalid) observations.
    :type valid_obs: bool, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded polygons.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    if method == "point":
        track = collect_ground_track(
            satellite, instrument, times, mask, crs, resolution, split_polygon
        )
        # filter ground track samples, if necessary
        if valid_obs is not None:
            return track[track.valid_obs == valid_obs].dissolve()
        else:
            return track.dissolve()
    elif method == "line":
        track = collect_orbit_track(satellite, instrument, times, mask)
        # filter ground track samples, if necessary
        if valid_obs is not None:
            track = track[track.valid_obs == valid_obs]
        # project points to zero elevation
        points = MultiPoint(track.geometry.apply(lambda p: Point(p.x, p.y)))
        # extract longitudes
        lons = np.array([points.geoms[i].x for i in range(len(points.geoms))])
        # split lines when crossing meridian or anti-meridian
        lines = [
            LineString(points)
            for points in np.split(
                points,
                1
                + np.where(
                    np.logical_or(np.diff(np.sign(lons)), np.diff(np.sign(lons)))
                )[0],
            )
        ]
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
                transform(to_crs, line).buffer(
                    track.swath_width.mean() / 2, resolution=resolution
                ),
            )
            for line in lines
        ]
        for p, polygon in enumerate(polygons):
            if np.median([e[0] for e in polygon.exterior.coords]) > 0:
                # fix longitudes that have been wrapped to the negative domain
                polygons[p] = Polygon(
                    [
                        [c[0] + 360 if c[0] < -90 else c[0], c[1]]
                        for c in polygon.exterior.coords
                    ],
                    [
                        [[c[0] + 360 if c[0] < -90 else c[0], c[1]] for c in i.coords]
                        for i in polygon.interiors
                    ],
                )
            else:
                # fix longitudes that have been wrapped to the positive domain
                polygons[p] = Polygon(
                    [
                        [c[0] - 360 if c[0] > 90 else c[0], c[1]]
                        for c in polygon.exterior.coords
                    ],
                    [
                        [[c[0] - 360 if c[0] > 90 else c[0], c[1]] for c in i.coords]
                        for i in polygon.interiors
                    ],
                )
            # clip the resulting polygon to the fixed domain
            polygons[p] = clip_by_rect(polygons[p], -180, -90, 180, 90)
        track = track.dissolve()
        # replace the geometry
        track.geometry = [unary_union(polygons)]
        return track
