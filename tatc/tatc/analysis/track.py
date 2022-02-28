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
from ..schemas.instrument import Instrument, DutyCycleScheme
from ..utils import (
    normalize_geometry,
    field_of_regard_to_swath_width,
    compute_orbit_period,
)


def collect_ground_track(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    target: Optional[Union[Polygon, MultiPolygon]] = None,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
) -> gpd.GeoDataFrame:
    """
    Collect ground track points for a satellite of interest.

    :param satellite: The observing satellite
    :type satellite: class:`tatc.schemas.satellite.Satellite`
    :param instrument: The instrument used to make observations
    :type instrument: class:`tatc.schemas.instrument.Instrument`
    :param times: The list of times to sample.
    :type times: list
    :param: target: A target region to prioritize operations.
    :type target: class:`shapely.Polygon` or :class:`shapely.MultyPolygon`, optional
    :param mask: A mask to constrain ground track.
    :type mask: class:`shapely.Polygon` or :class:`shapely.MultiPolygon`, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded points.
    :rtype: class:`geopandas.GeodataFrame`
    """

    # load the timescale and define starting and ending points
    ts = load.timescale()
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # load the ephemerides
    eph = load("de421.bsp")

    ts_times = [ts.from_datetime(time) for time in times]
    positions = [sat.at(time) for time in ts_times]
    subpoints = [wgs84.subpoint(position) for position in positions]
    points = [
        Point(
            subpoint.longitude.degrees, subpoint.latitude.degrees, subpoint.elevation.m
        )
        for subpoint in subpoints
    ]
    valid_obs = [
        instrument.is_valid_observation(eph, ts_times[i], positions[i])
        for i, time in enumerate(times)
    ]
    ops_intervals = instrument.generate_ops_intervals(eph, ts, sat, times, target)

    df = pd.DataFrame(
        [
            {
                "time": time,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "swath_width": field_of_regard_to_swath_width(
                    subpoints[i].elevation.m,
                    instrument.field_of_regard,
                ),
                "valid_obs": valid_obs[i]
                and (
                    instrument.duty_cycle >= 1
                    or any(ops_intervals.apply(lambda j: time in j))
                ),
                "geometry": points[i],
            }
            for i, time in enumerate(times)
        ],
        columns=[
            "time",
            "satellite",
            "instrument",
            "swath_width",
            "valid_obs",
            "geometry",
        ],
    )
    gdf = gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")
    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)


def collect_ground_track_swath(
    satellite: Satellite,
    instrument: Instrument,
    times: List[datetime],
    target: Optional[Union[Polygon, MultiPolygon]] = None,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    fast: bool = True,
    resolution: int = 4,
) -> gpd.GeoDataFrame:
    """
    Model the ground track swath for a satellite of interest.

    :param satellite: The observing satellite
    :type satellite: class:`tatc.schemas.satellite.Satellite`
    :param instrument: The observing instrument
    :type instrument: class:`tatc.schemas.instrument.Instrument`
    :param times: The list of times to sample
    :type times: list
    :param: target: A target region to prioritize operations.
    :type target: class:`shapely.Polygon` or :class:`shapely.MultyPolygon`, optional
    :param mask: A mask to constrain ground track.
    :type mask: class:`shapely.Polygon` or :class:`shapely.MultiPolygon`, optional
    :param target: An optionatarget region to prioritize operations.
    :param fast: Whether to use a fast (True) or accurate (False)
            projection. Fast uses EPSG:4087 (World Equidistant Cylindrical)
            to project distances, accurate finds the appropriate Unified
            Transverse Mercator (UTM) zone for each point.
    :type fast: bool, option
    :param resolution: Resolution of the projected swath
    :type resolution: int, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` with all recorded polygons.
    :rtype: class:`geopandas.GeoDataFrame`
    """
    # first, compute the ground track of the satellite
    gdf = collect_ground_track(satellite, instrument, times, target, mask)
    if len(gdf.index) == 0:
        return gdf
    # project points to ground
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
        def get_utm_epsg_code(p):
            results = query_utm_crs_info(
                datum_name="WGS 84", area_of_interest=AreaOfInterest(p.x, p.y, p.x, p.y)
            )
            # return the first code if exists; otherwise return a default value
            return results[0].code if len(results) > 0 else "4087"

        epsg = gdf.geometry.apply(get_utm_epsg_code)
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
    # hack to get plotting to work
    gdf = normalize_geometry(gdf)

    if mask is None:
        return gdf
    return gpd.clip(gdf, mask).reset_index(drop=True)
