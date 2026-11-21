"""
Methods to sample features: regions that exist for a period of time, such as
the storms or cloud systems of a nature run in an observing system
simulation experiment (OSSE), by the footprints of the instruments of
observing platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Geod
from skyfield.api import wgs84

from .._warnings import warn
from ..schemas import AllPlatforms
from ._tracks import _get_viewing_geometry
from .check import (
    SATELLITE,
    _check_times,
    _get_instruments,
    _get_platform_kinds,
    _get_platform_lookup,
    _get_platforms,
    _group_platforms,
    _PlatformKind,
)
from .footprints import (
    _ELEVATION_TOLERANCE,
    _get_footprint_elevations,
    _get_kind_elevations,
    _get_kind_without_elevation,
    collect_footprints,
)

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


def _get_feature_elevations(features: gpd.GeoDataFrame) -> np.ndarray:
    """
    Gets the elevation of each feature: the z-coordinate of its vertices,
    which must be constant (within a meter), or `NaN` for a feature without
    z-coordinates.

    Args:
        features (geopandas.GeoDataFrame): The features.

    Returns:
        numpy.ndarray: the elevation (meters) of each feature
    """
    shapes = features.geometry.to_numpy()
    elevations = np.full(len(shapes), np.nan)
    rows = np.flatnonzero(shapely.has_z(shapes))
    if len(rows) == 0:
        return elevations
    coordinates, index = shapely.get_coordinates(
        shapes[rows], include_z=True, return_index=True
    )
    low = np.full(len(rows), np.inf)
    high = np.full(len(rows), -np.inf)
    np.minimum.at(low, index, coordinates[:, 2])
    np.maximum.at(high, index, coordinates[:, 2])
    spans = high - low > _ELEVATION_TOLERANCE
    if np.any(spans):
        i = np.flatnonzero(spans)[0]
        raise ValueError(
            "each feature must be at a constant elevation (its z-coordinate), at "
            "which its footprints are computed, not feature "
            f"{features.index[rows[i]]!r} from {low[i]:.0f} to {high[i]:.0f} m"
        )
    elevations[rows] = _get_footprint_elevations(shapes[rows])
    return elevations


def _get_group_elevation(
    elevation: float | Mapping[str, float | None] | None, height: float
) -> float | Mapping[str, float | None] | None:
    """
    Gets the elevation at which to compute the footprints of features at a
    height: that given for each kind of platform, or else the height (if
    not `NaN`).

    Args:
        elevation (float | Mapping[str, float | None] | None): The
            elevation (meters), or that of each kind of platform by its column.
        height (float): The elevation (meters) of the features, or `NaN`.

    Returns:
        float | Mapping[str, float | None] | None: the elevation of each kind
    """
    if np.isnan(height):
        return elevation
    return {
        kind.column: float(height) if value is None else value
        for kind, value in _get_kind_elevations(elevation).items()
    }


def _get_feature_pieces(
    rows: np.ndarray,
    shapes: np.ndarray,
    starts: pd.DatetimeIndex,
    ends: pd.DatetimeIndex,
    platforms: list[AllPlatforms],
    times: list[datetime],
    instrument_index: int | None,
    elevation: float | Mapping[str, float | None] | None,
) -> pd.DataFrame | None:
    """
    Gets the pieces of features: each feature's overlap (of positive area)
    with each valid footprint (at an elevation) while the feature exists.

    Args:
        rows (numpy.ndarray): The positions of the features.
        shapes (numpy.ndarray): The (two-dimensional) geometry of each feature.
        starts (pandas.DatetimeIndex): The start of each feature.
        ends (pandas.DatetimeIndex): The end of each feature.
        platforms (list[AllPlatforms]): The observing platforms.
        times (list[datetime.datetime]): The times at which to sample.
        instrument_index (int | None): The index of the observing instrument.
        elevation (float | Mapping[str, float | None] | None): The
            elevation of the footprints (see `collect_footprints`).

    Returns:
        pandas.DataFrame | None: the pieces: the feature (its position),
            platform (key), kind, instrument, time, footprint elevation
            (`height`), geometry, and area (square degrees), or `None`
            without any
    """
    # the footprints, clipped to the bounds of the features
    mask = shapely.box(*shapely.total_bounds(shapes[rows]))
    track = collect_footprints(platforms, times, instrument_index, elevation, mask)
    track = track[track["valid_obs"].eq(True).to_numpy()].reset_index(drop=True)
    if len(track) == 0:
        return None
    # the elevation of each footprint, projected at a constant elevation
    heights = _get_footprint_elevations(track.geometry.to_numpy())
    footprints = shapely.force_2d(track.geometry.to_numpy())
    # each feature's overlap with each footprint while it exists
    feature, footprint = shapely.STRtree(footprints).query(
        shapes[rows], predicate="intersects"
    )
    feature = rows[feature]
    time = pd.DatetimeIndex(track["time"])
    exists = np.asarray(
        (starts[feature] <= time[footprint]) & (time[footprint] < ends[feature])
    )
    feature, footprint = feature[exists], footprint[exists]
    if len(feature) == 0:
        return None
    pieces = pd.DataFrame(
        {
            "feature": feature,
            "platform": _get_platforms(track)[footprint],
            "kind": _get_platform_kinds(track)[footprint],
            "instrument": _get_instruments(track)[footprint],
            "time": time[footprint],
            "height": heights[footprint],
            "geometry": shapely.intersection(shapes[feature], footprints[footprint]),
        }
    )
    pieces["area"] = shapely.area(pieces["geometry"].to_numpy())
    pieces = pieces[pieces["area"] > 0]
    return pieces if len(pieces) > 0 else None


def _warn_feature_elevations(features: np.ndarray, footprints: np.ndarray) -> None:
    """
    Warns if features overlap footprints projected to other elevations.

    Args:
        features (numpy.ndarray): The elevation (meters) of each feature
            (`NaN` if unknown).
        footprints (numpy.ndarray): The elevation (meters) of the footprint
            that overlaps each feature (`NaN` if unknown).
    """
    differs = np.abs(features - footprints) > _ELEVATION_TOLERANCE
    if not np.any(differs):
        return
    i = np.flatnonzero(differs)[0]
    warn(
        f"features overlap footprints at other elevations (such as a feature at "
        f"{features[i]:.0f} m and a footprint at {footprints[i]:.0f} m), which are "
        "offset from them along the line of sight: leave elevation unset (or "
        "unset for the kind of platform) to compute the footprints at the "
        "elevation of each feature"
    )


def collect_feature_observations(
    features: gpd.GeoDataFrame,
    platforms: AllPlatforms | list[AllPlatforms],
    times: list[datetime],
    instrument_index: int | None = None,
    elevation: float | Mapping[str, float | None] | None = None,
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

    A feature with z-coordinates (constant over its vertices, such as a cloud
    top at its height above the WGS 84 ellipsoid) is sampled by footprints
    computed at its elevation, for each kind of platform without a given
    `elevation`: the footprints are computed once for each distinct
    elevation of the features (and once for those without z-coordinates),
    so features at many elevations (such as the continuous cloud-top heights
    of a nature run) are faster to sample if rounded to a few levels. A
    feature at another elevation than its footprints warns, as they are
    offset from it along the line of sight.

    Args:
        features (geopandas.GeoDataFrame): The features: their geometries
                (optionally at an elevation, their z-coordinate), and the
                `start` and `end` (timezone-aware) of the period in which
                each exists.
        platforms (AllPlatforms | list[AllPlatforms]): The observing
                platform(s), of any kinds (see `tatc.schemas.AllPlatforms`).
        times (list[datetime.datetime]): The times at which to sample.
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` (default) for every instrument
                of each platform.
        elevation (float | Mapping[str, float | None] | None): The
                elevation (meters) above the WGS 84 ellipsoid at which to
                compute footprints, or that of each kind of platform by its
                column (see `collect_footprints`): by default, the elevation
                of each feature (its z-coordinate), or, for a feature without
                one, the ellipsoid for satellites and aircraft and the ground
                elevation of each surface platform; required for radar
                stations unless each feature has an elevation.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations: the feature
            (its index label in `features`), the platform (named in the
            column of its kind, such as `satellite`), instrument (`None` for
            a radar station), time, the viewed part (geometry, whose
            z-coordinate is the elevation of the footprints), its fraction
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
    lookup = _get_platform_lookup(platforms)
    if len(features) == 0 or len(times) == 0 or len(lookup) == 0:
        return _get_empty_feature_observations(kinds)
    # the elevation of each feature (NaN without z-coordinates)
    feature_heights = _get_feature_elevations(features)
    missing = _get_kind_without_elevation(groups, _get_kind_elevations(elevation))
    if missing is not None and np.isnan(feature_heights).any():
        raise ValueError(
            f"elevation is required to collect the footprints of {missing.label}s, "
            "unless each feature has an elevation (its z-coordinate)"
        )
    shapes = shapely.force_2d(features.geometry.to_numpy())
    starts = pd.DatetimeIndex(features["start"])
    ends = pd.DatetimeIndex(features["end"])
    # the pieces of the features at each elevation (or without one), whose
    # footprints are computed at that elevation
    groups_of_pieces = [
        _get_feature_pieces(
            np.flatnonzero(
                np.isnan(feature_heights)
                if np.isnan(height)
                else feature_heights == height
            ),
            shapes,
            starts,
            ends,
            list(lookup.values()),
            times,
            instrument_index,
            _get_group_elevation(elevation, height),
        )
        for height in pd.unique(feature_heights)
    ]
    groups_of_pieces = [group for group in groups_of_pieces if group is not None]
    if len(groups_of_pieces) == 0:
        return _get_empty_feature_observations(kinds)
    pieces = pd.concat(groups_of_pieces, ignore_index=True)
    _warn_feature_elevations(
        feature_heights[pieces["feature"].to_numpy()], pieces["height"].to_numpy()
    )
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
        geometry = _get_viewing_geometry(
            lookup[platform].get_track(list(group["time"])),
            wgs84.latlon(latitude[index], longitude[index], group["height"].to_numpy()),
        )
        for column in _VIEWING_GEOMETRY:
            observations.loc[index, column] = geometry[column]
    # the feature's index label in the features
    observations["feature"] = features.index.to_numpy()[observations["feature"]]
    observations["time"] = pd.DatetimeIndex(observations["time"]).as_unit("ns")
    # the viewed part at the elevation of the footprints, as they are
    observations["geometry"] = shapely.force_3d(
        shapely.force_2d(viewed), np.nan_to_num(observations["height"].to_numpy())
    )
    # the names of the empty frame's dtypes (a kind without observations, or
    # radar stations without instruments, would otherwise be of objects)
    empty = _get_empty_feature_observations(kinds)
    for column in [*(kind.column for kind in kinds), "instrument"]:
        observations[column] = observations[column].astype(empty[column].dtype)
    observations = gpd.GeoDataFrame(
        observations[list(empty.columns)],
        geometry="geometry",
        crs="EPSG:4326",
    )
    return observations.sort_values(["feature", "time"], kind="stable").reset_index(
        drop=True
    )
