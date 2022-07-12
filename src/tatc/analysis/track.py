# -*- coding: utf-8 -*-
"""
Methods to generate coverage statistics.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from typing import List, Union, Optional
import pandas as pd
import geopandas as gpd
from datetime import datetime, timedelta
from skyfield.api import load, wgs84, EarthSatellite
from shapely.geometry import Point, Polygon, MultiPolygon
from pyproj.database import query_utm_crs_info
from pyproj.aoi import AreaOfInterest

from ..schemas.satellite import Satellite
from ..schemas.instrument import Instrument
from ..utils import (
    normalize_geometry,
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
    fast: bool = True,
    resolution: int = 4,
    fix_polygon_longitudes: bool = True,
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
    :param fast: Whether to use a fast (True) or accurate (False)
            projection. Fast uses EPSG:4087 (World Equidistant Cylindrical)
            to project distances, accurate finds the appropriate Unified
            Transverse Mercator (UTM) zone for each point.
    :type fast: bool, option
    :param resolution: Shapely buffer resolution of the projected swath
    :type resolution: int, optional
    :param fix_polygon_longitudes: True, if polygon longitudes crossing the
            anti-meridian (180 deg longitude) should be transformed from
            the [-180, 180] to the [0,360] longitude domain to aid plotting.
    :type fix_polygon_longitudes: bool, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded polygons.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    # first, compute the orbit track of the satellite
    gdf = collect_orbit_track(satellite, instrument, times, mask)
    if gdf.empty:
        return gdf
    # project points to zero elevation
    gdf["geometry"] = gdf["geometry"].apply(lambda p: Point(p.x, p.y))
    # at each point, draw a buffer equivalent to the swath radius
    if fast:
        # note: uses EPSG:4087 (World Equidistant Cylindrical) to approximate cartesian coordinates
        gdf["geometry"] = gpd.GeoSeries(
            gdf.to_crs("EPSG:4087").apply(
                lambda r: r["geometry"].buffer(
                    r["swath_width"] / 2, resolution=resolution
                ),
                axis=1,
            ),
            crs="EPSG:4087",
        ).to_crs("EPSG:4326")
        # previously used EPSG:3857 (Pseudo-Mercator) which has poor accuracy near the poles
    else:
        # note: uses UTM zones but is 10x slower
        def _get_utm_epsg_code(p):
            results = query_utm_crs_info(
                datum_name="WGS 84", area_of_interest=AreaOfInterest(p.x, p.y, p.x, p.y)
            )
            # return the first code if exists; otherwise return a default value
            return results[0].code if len(results) > 0 else "4087"

        epsg = gdf.geometry.apply(_get_utm_epsg_code)
        for code in epsg.unique():
            gdf.loc[epsg == code, "geometry"] = gpd.GeoSeries(
                gdf[epsg == code]
                .to_crs("EPSG:" + code)
                .apply(
                    lambda r: r["geometry"].buffer(
                        r["swath_width"] / 2, resolution=resolution
                    ),
                    axis=1,
                ),
                crs="EPSG:" + code,
            ).to_crs("EPSG:4326")

    if fix_polygon_longitudes:
        gdf = normalize_geometry(gdf)

    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)
