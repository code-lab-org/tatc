"""
Methods to sample features: regions that exist for a period of time, such as
the storms or cloud systems of a nature run in an observing system
simulation experiment (OSSE), by the footprints of satellites' instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Geod
from skyfield.api import wgs84

from ..schemas import Satellite
from ..utils.time import _to_time
from .check import _check_satellites
from .ground_track import collect_ground_track

_GEOD = Geod(ellps="WGS84")


def _get_empty_feature_observations() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for feature observation results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "feature": pd.Series([], dtype="int"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "viewed_fraction": pd.Series([], dtype="float"),
        "viewed_area": pd.Series([], dtype="float"),
        "sat_alt": pd.Series([], dtype="float"),
        "sat_az": pd.Series([], dtype="float"),
        "sat_range": pd.Series([], dtype="float"),
        "sat_height": pd.Series([], dtype="float"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _check_features(features: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Checks features: a data frame of geometries with the `start` and `end`
    of the period in which each exists.

    Args:
        features (geopandas.GeoDataFrame): The features.

    Returns:
        geopandas.GeoDataFrame: the features, in longitude/latitude coordinates.
    """
    if not isinstance(features, gpd.GeoDataFrame):
        raise TypeError(
            f"features must be a GeoDataFrame, not a {type(features).__name__}"
        )
    for column in ("start", "end"):
        if column not in features.columns:
            raise ValueError(f"features must have a '{column}' column")
        if not isinstance(features[column].dtype, pd.DatetimeTZDtype):
            raise ValueError(f"features '{column}' must be timezone-aware times")
    if features.crs is not None and not features.crs.equals("EPSG:4326"):
        features = features.to_crs("EPSG:4326")
    return features


def collect_feature_observations(
    features: gpd.GeoDataFrame,
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observations of features: regions that exist for a period
    of time (between the `start` and `end` of each feature), such as the
    storms or cloud systems of a nature run in an observing system
    simulation experiment (OSSE). The footprint of each instrument at each
    time (see `tatc.analysis.ground_track.collect_ground_track`) samples the
    features that exist at that time (from their `start`, inclusive, to
    their `end`, exclusive). The footprints should tile each swath, for
    example with an along-track field of view tailored to the time step
    (see `tatc.utils.compute_along_track_field_of_view`).

    Each observation is one satellite's (and instrument's) view of a
    feature: the part of the feature within the valid footprints of all
    times at which the feature exists, which is only part of it if the
    feature extends beyond the swath. Its time is that of the footprint
    that covers most of the feature, and the viewing geometry (the
    satellite's elevation angle `sat_alt` and azimuth `sat_az`, seen from
    the centroid of the viewed part, its slant range `sat_range`, and its
    height above the WGS 84 ellipsoid `sat_height`) is that at this time.
    For example, the size of a footprint at the viewed part (from its size
    at nadir) grows with `sat_range / sat_height` along the line of sight
    and, across it, by the secant of the zenith angle, `90 - sat_alt`.

    Args:
        features (geopandas.GeoDataFrame): The features: their geometries,
                and the `start` and `end` (timezone-aware) of the period in
                which each exists.
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (list[datetime.datetime]): The times at which to sample.
        instrument_index (int | None): The index of the observing instrument
                in each satellite, or `None` (default) for every instrument
                of each satellite.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations: the feature
            (its index label in `features`), satellite, instrument, time, the
            viewed part (geometry), its fraction of the feature's area
            (`viewed_fraction`), its area on the WGS 84 ellipsoid (square
            meters, `viewed_area`), and the viewing geometry (degrees and
            meters), sorted by feature and time.
    """
    features = _check_features(features)
    satellites = _check_satellites(satellites)
    if len(features) == 0 or len(times) == 0:
        return _get_empty_feature_observations()
    # the footprints, clipped to the bounds of the features
    mask = shapely.box(*features.total_bounds)
    track = collect_ground_track(
        satellites, times, instrument_index=instrument_index, mask=mask
    )
    track = track[track.valid_obs].reset_index(drop=True)
    if len(track) == 0:
        return _get_empty_feature_observations()
    footprints = shapely.force_2d(track.geometry.to_numpy())
    shapes = features.geometry.to_numpy()
    # the pieces: each feature's overlap with each footprint while it exists
    feature, footprint = shapely.STRtree(footprints).query(
        shapes, predicate="intersects"
    )
    time = pd.DatetimeIndex(track["time"])
    exists = np.asarray(
        (pd.DatetimeIndex(features["start"])[feature] <= time[footprint])
        & (time[footprint] < pd.DatetimeIndex(features["end"])[feature])
    )
    feature, footprint = feature[exists], footprint[exists]
    pieces = pd.DataFrame(
        {
            "feature": feature,
            "satellite": track["satellite"].to_numpy()[footprint],
            "instrument": track["instrument"].to_numpy()[footprint],
            "time": time[footprint],
            "geometry": shapely.intersection(shapes[feature], footprints[footprint]),
        }
    )
    pieces["area"] = shapely.area(pieces["geometry"].to_numpy())
    pieces = pieces[pieces["area"] > 0]
    if len(pieces) == 0:
        return _get_empty_feature_observations()
    # each observation: the union of the pieces of a feature viewed by a
    # satellite's instrument, at the time of the piece with the largest area
    observations = (
        pieces.sort_values("area", ascending=False, kind="stable")
        .groupby(["feature", "satellite", "instrument"], sort=False)
        .agg(time=("time", "first"), geometry=("geometry", shapely.union_all))
        .reset_index()
    )
    viewed = observations["geometry"].to_numpy()
    shape = shapes[observations["feature"].to_numpy()]
    observations["viewed_fraction"] = np.minimum(
        1, shapely.area(viewed) / shapely.area(shape)
    )
    observations["viewed_area"] = [
        abs(_GEOD.geometry_area_perimeter(geometry)[0]) for geometry in viewed
    ]
    # the viewing geometry from the centroid of the viewed part
    centroid = shapely.centroid(viewed)
    longitude, latitude = shapely.get_x(centroid), shapely.get_y(centroid)
    for column in ("sat_alt", "sat_az", "sat_range", "sat_height"):
        observations[column] = np.nan
    orbits = {satellite.name: satellite.orbit.to_gp_orbit() for satellite in satellites}
    for name, group in observations.groupby("satellite"):
        index = group.index.to_numpy()
        orbit_track = orbits[name].get_orbit_track_at_time(
            _to_time(list(group["time"]))
        )
        target = wgs84.latlon(latitude[index], longitude[index])
        altitude, azimuth, distance = (orbit_track - target.at(orbit_track.t)).altaz()
        observations.loc[index, "sat_alt"] = np.atleast_1d(altitude.degrees)
        observations.loc[index, "sat_az"] = np.atleast_1d(azimuth.degrees)
        observations.loc[index, "sat_range"] = np.atleast_1d(distance.m)
        observations.loc[index, "sat_height"] = np.atleast_1d(
            wgs84.height_of(orbit_track).m
        )
    # the feature's index label in the features
    observations["feature"] = features.index.to_numpy()[observations["feature"]]
    observations["time"] = pd.DatetimeIndex(observations["time"]).as_unit("ns")
    observations = gpd.GeoDataFrame(
        observations[list(_get_empty_feature_observations().columns)],
        geometry="geometry",
        crs="EPSG:4326",
    )
    return observations.sort_values(["feature", "time"], kind="stable").reset_index(
        drop=True
    )
