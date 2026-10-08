"""
Methods to collect instrument ground tracks.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric

from ..constants import de421
from ..schemas import AllInstruments, PointedInstrument, Satellite
from ..utils.geometry import split_polygon
from .check import _check_satellite
from .region_sampling import compute_region_access_periods


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


def _get_mask_geometry(
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries,
) -> Polygon | MultiPolygon:
    """
    Gets a mask as a single geometry, split along the anti-meridian and
    poles (see `tatc.utils.geometry.split_polygon`) so that it compares
    correctly with footprints in the standard longitude range.

    Args:
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries):
                The mask, in WGS84 (lon/lat) coordinates.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: the mask geometry
    """
    if isinstance(mask, (Polygon, MultiPolygon)):
        geometry = mask
    elif isinstance(mask, gpd.GeoDataFrame):
        geometry = mask.dissolve().iloc[0].geometry
    else:
        geometry = mask.union_all()
    return split_polygon(geometry)


def _cull_orbit_track(
    satellite: Satellite,
    instrument: AllInstruments,
    times: list[datetime],
    mask: Polygon | MultiPolygon,
    elevation: float,
) -> tuple[Geocentric, list[Polygon | MultiPolygon]] | None:
    """
    Propagates the orbit track only at the times when the instrument's
    footprint intersects a mask. The times are first culled to the periods
    when the instrument's field of regard may observe any part of the mask
    (see `tatc.analysis.region_sampling.compute_region_access_periods`), so
    that the orbit is not propagated far from the mask, and then to those
    whose footprint intersects the mask.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (AllInstruments): The observing instrument.
        times (list[datetime.datetime]): The times to sample.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The mask, in WGS84 (lon/lat) coordinates (see `_get_mask_geometry`).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
                at which to project footprints.

    Returns:
        tuple[skyfield.positionlib.Geocentric, list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]] | None:
            the culled orbit track and its footprints, or None if no
            footprint intersects the mask
    """
    periods = compute_region_access_periods(
        mask,
        satellite,
        min(times),
        max(times),
        instrument.field_of_regard,
        elevation,
    )
    if periods.empty:
        return None
    # select the times within a period (which are sorted and disjoint)
    sample = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    left = pd.DatetimeIndex([period.left for period in periods])
    right = pd.DatetimeIndex([period.right for period in periods])
    index = np.searchsorted(left, sample, side="right") - 1
    selected = (index >= 0) & (sample <= right[np.maximum(index, 0)])
    if not any(selected):
        return None
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(
        [time for time, keep in zip(times, selected) if keep]
    )
    # cull orbit track to observable footprint
    footprint = instrument.compute_footprint(orbit_track, elevation=elevation)
    mask_intersects_footprint = [mask.intersects(f) for f in footprint]
    if not any(mask_intersects_footprint):
        return None
    return orbit_track[mask_intersects_footprint], [
        f for f, keep in zip(footprint, mask_intersects_footprint) if keep
    ]


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
    time, projected to a specified elevation by computing the exact
    ray/WGS-84-geoid intersection (see
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

    _check_satellite(satellite)
    if len(times) == 0:
        return _get_empty_ground_track()
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if mask is not None:
        mask = _get_mask_geometry(mask)
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask,
        # reusing the footprints computed to cull it
        culled = _cull_orbit_track(satellite, instrument, times, mask, elevation)
        if culled is None:
            return _get_empty_ground_track()
        orbit_track, geometries = culled
    else:
        # propagate orbit
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
        # compute footprints (exact ray/WGS-84-geoid intersection)
        geometries = instrument.compute_footprint(orbit_track, None, elevation)
    # compute targets
    target = instrument.compute_footprint_center(orbit_track, elevation)
    # determine observation validity
    valid_obs = instrument.is_valid_observation(orbit_track, target)
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
        track = gpd.clip(track, mask)
    # sort rows by time (clipping does not preserve their order), keeping
    # the original order of rows at the same time
    return track.sort_index().sort_values("time", kind="stable").reset_index(drop=True)


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
    _check_satellite(satellite)
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

    _check_satellite(satellite)
    if len(times) == 0:
        return _get_empty_ground_track()
    # select the observing instrument
    instrument = satellite.instruments[instrument_index]
    if not isinstance(instrument, PointedInstrument) or not instrument.is_rectangular:
        raise ValueError(
            "Ground pixels are only compatible with rectangular PointedInstrument instances"
        )
    if mask is not None:
        mask = _get_mask_geometry(mask)
    if mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask
        culled = _cull_orbit_track(satellite, instrument, times, mask, elevation)
        if culled is None:
            return _get_empty_ground_track()
        orbit_track, _ = culled
    else:
        # propagate orbit
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
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
        gdf = gpd.clip(gdf, mask)
    # sort rows by time (clipping does not preserve their order), keeping
    # the original order of rows at the same time
    return gdf.sort_index().sort_values("time", kind="stable").reset_index(drop=True)
