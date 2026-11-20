"""
Methods to sample features: regions that exist for a period of time, such as
the storms or cloud systems of a nature run in an observing system
simulation experiment (OSSE), by the footprints of the instruments of
observing platforms.

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

from ..schemas import AllPlatforms
from .check import (
    SATELLITE,
    _PlatformKind,
    _check_times,
    _get_instruments,
    _get_platform_key,
    _get_platform_kinds,
    _get_platforms,
    _group_platforms,
)
from .footprints import collect_footprints

_GEOD = Geod(ellps="WGS84")

_VIEWING_GEOMETRY = ("platform_alt", "platform_az", "platform_range", "platform_height")
"""The columns of the viewing geometry of a feature observation."""


def _get_empty_feature_observations(
    kinds: list[_PlatformKind] | None = None,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for feature observation results.

    Args:
        kinds (list[_PlatformKind] | None): The kinds of observing
            platforms, whose columns name the platforms, by default
            satellites.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    kinds = kinds or [SATELLITE]
    columns = {
        "feature": pd.Series([], dtype="int"),
        **{kind.column: pd.Series([], dtype="str") for kind in kinds},
        "instrument": pd.Series([], dtype="str"),
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "viewed_fraction": pd.Series([], dtype="float"),
        "viewed_area": pd.Series([], dtype="float"),
        **{column: pd.Series([], dtype="float") for column in _VIEWING_GEOMETRY},
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
    platforms: AllPlatforms | list[AllPlatforms],
    times: list[datetime],
    instrument_index: int | None = None,
    elevation: float | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observations of features: regions that exist for a period
    of time (between the `start` and `end` of each feature), such as the
    storms or cloud systems of a nature run in an observing system
    simulation experiment (OSSE). The footprint of each instrument at each
    time (see `tatc.analysis.footprints.collect_footprints`, for satellites,
    aircraft, fixed surface platforms, and radar stations) samples the
    features that exist at that time (from their `start`, inclusive, to
    their `end`, exclusive). The footprints should tile each swath, for
    example with an along-track field of view tailored to the time step
    (see `tatc.utils.compute_along_track_field_of_view`).

    Each observation is one platform's (and instrument's) view of a
    feature: the part of the feature within the valid footprints of all
    times at which the feature exists, which is only part of it if the
    feature extends beyond the swath. Its time is that of the footprint
    that covers most of the feature (the earliest of those that cover as
    much, as the footprints of a fixed platform at each time do), and the
    viewing geometry is that at
    this time: the platform's elevation angle (`platform_alt`) and azimuth
    (`platform_az`), seen from the centroid of the viewed part, its slant
    range (`platform_range`), and its height above the WGS 84 ellipsoid
    (`platform_height`), for any kind of platform, as for the angles of its
    footprints. For example, the size of a satellite's footprint at the
    viewed part (from its size at nadir) grows with
    `platform_range / platform_height` along the line of sight and, across
    it, by the secant of the zenith angle, `90 - platform_alt`.

    Args:
        features (geopandas.GeoDataFrame): The features: their geometries,
                and the `start` and `end` (timezone-aware) of the period in
                which each exists.
        platforms (AllPlatforms | list[AllPlatforms]): The observing
                platform(s), of any kinds (see `tatc.schemas.AllPlatforms`).
        times (list[datetime.datetime]): The times at which to sample.
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` (default) for every instrument
                of each platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid at which to compute footprints (see
                `collect_footprints`): by default, the ellipsoid for
                satellites and aircraft, and the ground elevation of each
                surface platform; required for radar stations.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations: the feature
            (its index label in `features`), the platform (named in the
            column of its kind, such as `satellite`), instrument (`None` for
            a radar station), time, the viewed part (geometry), its fraction
            of the feature's area (`viewed_fraction`), its area on the WGS 84
            ellipsoid (square meters, `viewed_area`), and the viewing
            geometry (`platform_alt`, `platform_az`, `platform_range`, and
            `platform_height`, degrees and meters), sorted by feature and
            time.
    """
    features = _check_features(features)
    groups = _group_platforms(platforms)
    _check_times(times)
    kinds = [kind for kind, members in groups.items() if members]
    # each platform, by its key in the footprints (platforms of a kind have
    # distinct names, see `_group_platforms`)
    lookup = {
        _get_platform_key(kind.column, platform.name): platform
        for kind, members in groups.items()
        for platform in members
    }
    if len(features) == 0 or len(times) == 0 or len(lookup) == 0:
        return _get_empty_feature_observations(kinds)
    # the footprints, clipped to the bounds of the features
    mask = shapely.box(*features.total_bounds)
    track = collect_footprints(
        list(lookup.values()), times, instrument_index, elevation, mask
    )
    track = track[track["valid_obs"].eq(True).to_numpy()].reset_index(drop=True)
    if len(track) == 0:
        return _get_empty_feature_observations(kinds)
    # the elevation of each footprint: the height of its (non-empty)
    # geometry, projected at a constant elevation
    coordinates, row = shapely.get_coordinates(
        track.geometry.to_numpy(), include_z=True, return_index=True
    )
    rows, first = np.unique(row, return_index=True)
    heights = np.full(len(track), np.nan)
    heights[rows] = coordinates[first, 2]
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
    instruments = _get_instruments(track)
    pieces = pd.DataFrame(
        {
            "feature": feature,
            "platform": _get_platforms(track)[footprint],
            "kind": _get_platform_kinds(track)[footprint],
            "instrument": instruments[footprint],
            "time": time[footprint],
            "height": heights[footprint],
            "geometry": shapely.intersection(shapes[feature], footprints[footprint]),
        }
    )
    pieces["area"] = shapely.area(pieces["geometry"].to_numpy())
    pieces = pieces[pieces["area"] > 0]
    if len(pieces) == 0:
        return _get_empty_feature_observations(kinds)
    # each observation: the union of the pieces of a feature viewed by a
    # platform's instrument, at the time of the piece with the largest area
    # (the earliest of those as large)
    observations = (
        pieces.sort_values(["area", "time"], ascending=[False, True], kind="stable")
        .groupby(["feature", "platform", "instrument"], sort=False, dropna=False)
        .agg(
            kind=("kind", "first"),
            time=("time", "first"),
            height=("height", "first"),
            geometry=("geometry", shapely.union_all),
        )
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
    # each platform, named in the column of its kind
    for kind in kinds:
        observations[kind.column] = [
            lookup[platform].name if column == kind.column else None
            for platform, column in zip(observations["platform"], observations["kind"])
        ]
    for column in _VIEWING_GEOMETRY:
        observations[column] = np.nan
    # the viewing geometry from the centroid of the viewed part, at the
    # elevation of the footprints
    centroid = shapely.centroid(viewed)
    longitude, latitude = shapely.get_x(centroid), shapely.get_y(centroid)
    for platform, group in observations.groupby("platform"):
        index = group.index.to_numpy()
        observer = lookup[platform]
        platform_track = observer.get_track(list(group["time"]))
        target = wgs84.latlon(
            latitude[index], longitude[index], group["height"].to_numpy()
        )
        altitude, azimuth, distance = (
            platform_track - target.at(platform_track.t)
        ).altaz()
        observations.loc[index, "platform_alt"] = np.atleast_1d(altitude.degrees)
        observations.loc[index, "platform_az"] = np.atleast_1d(azimuth.degrees)
        observations.loc[index, "platform_range"] = np.atleast_1d(distance.m)
        observations.loc[index, "platform_height"] = np.atleast_1d(
            wgs84.height_of(platform_track).m
        )
    # the feature's index label in the features
    observations["feature"] = features.index.to_numpy()[observations["feature"]]
    observations["time"] = pd.DatetimeIndex(observations["time"]).as_unit("ns")
    observations = gpd.GeoDataFrame(
        observations[list(_get_empty_feature_observations(kinds).columns)],
        geometry="geometry",
        crs="EPSG:4326",
    )
    return observations.sort_values(["feature", "time"], kind="stable").reset_index(
        drop=True
    )
