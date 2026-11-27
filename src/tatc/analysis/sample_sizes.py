"""
Methods to compute the sizes of the samples of observations: the extent on
the ground of an instrument's pixel at each observed target, and the
spacing between pixels, which bound the detail an observation can resolve.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from skyfield.api import wgs84

from ..schemas import (
    AllInstruments,
    AllPlatforms,
    ConicalInstrument,
    LimbSounder,
    OccultationReceiver,
    PointedInstrument,
)
from ..utils.projection import VelocityFrame, compute_view_angles
from ._tracks import _get_viewing_geometry
from .check import (
    _check_platform_lookup,
    _get_instruments,
    _get_observation_times,
    _get_platform_kind,
    _get_platform_lookup,
    _get_platforms,
)
from .footprints import _get_footprint_elevations

_COLUMNS = [
    "platform_range",
    "scan_angle",
    "sample_along",
    "sample_cross",
    "spacing_along",
    "spacing_cross",
    "sample_azimuth",
]


def _get_sized_instrument(
    platform: AllPlatforms,
    name: str | None,
    instrument: AllInstruments | Mapping[str, AllInstruments] | None,
) -> AllInstruments | None:
    """
    Gets the instrument whose samples to size at the observations of a
    platform's instrument (by name): the given instrument (or that given for
    the name), or else the observing instrument, either bound to the
    platform's attitude (see `get_instrument`), or `None` for a platform
    without instruments (a radar station) or observations without them (such
    as limb soundings, without an `instrument` column).

    Args:
        platform (AllPlatforms): The observing platform.
        name (str | None): The name of the observing instrument.
        instrument (AllInstruments | Mapping[str, AllInstruments] | None): The
            instrument to size, or one for each observing instrument's name.

    Returns:
        AllInstruments | None: the instrument, bound to the platform's attitude
    """
    kind = _get_platform_kind(platform)
    if kind is None or not kind.instrumented:
        return None
    sized = instrument.get(name) if isinstance(instrument, Mapping) else instrument
    if sized is not None:
        # bound as the platform's own instrument would be
        return platform.model_copy(update={"instruments": [sized]}).get_instrument(0)
    if name is None or pd.isna(name):
        # observations without instruments (such as limb soundings)
        return None
    index = next(
        (i for i, carried in enumerate(platform.instruments) if carried.name == name),
        None,
    )
    if index is None:
        raise ValueError(
            f"{kind.label} {platform.name!r} has no instrument named {name!r}"
        )
    return platform.get_instrument(index)


def compute_sample_sizes(
    observations: gpd.GeoDataFrame,
    platforms: AllPlatforms | list[AllPlatforms],
    instrument: AllInstruments | Mapping[str, AllInstruments] | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the sizes of the samples of observations: for each observation,
    the extent on the ground of the observing instrument's pixel at the
    target, along track and across track (`sample_along` and
    `sample_cross`, meters), and the spacing between the centers of
    neighboring pixels (`spacing_along` and `spacing_cross`, meters), which
    is larger than the pixel if there are gaps between pixels (a negative
    oversampling). An observation cannot resolve detail finer than the
    larger of the two. The along-track size lies on the ground in the
    direction of the `sample_azimuth` (degrees clockwise from north, see
    `PointedInstrument.compute_sample_azimuth` and
    `ConicalInstrument.compute_sample_azimuth`), and the cross-track size
    about perpendicular to it, which orient the sample's footprint (as to
    weight a field over it, see `sample_field`). The sizes are those of a `PointedInstrument`'s
    pixels (see `PointedInstrument.compute_sample_size` and
    `PointedInstrument.compute_sample_spacing`) or of a
    `ConicalInstrument`'s beam and samples along its scanned arc (see
    `ConicalInstrument.compute_sample_size` and
    `ConicalInstrument.compute_sample_spacing`, `NaN` unless its
    `beam_width` and `scan_sample_spacing` are given), and `NaN` for other
    instruments, whose samples are not defined. Along and across track are those
    of the instrument's velocity frame: for a geostationary platform, whose
    velocity is east in the inertial frame (and nearly zero in the Earth-fixed
    frame), a `PointedInstrument` in the inertial frame gives the sizes east-west
    (`sample_along`) and north-south (`sample_cross`), as of the pixels of a fixed
    grid of scan angles.

    Also computes the platform's slant range from the target
    (`platform_range`, meters) and the target's scan angle (`scan_angle`,
    degrees): its angle across track from nadir, positive to the left of
    the direction of motion (see `tatc.utils.compute_view_angles`),
    measured from the body's axes if the platform has an attitude, which
    gives its position in the swath. A radar station has no instruments: its
    observations have a slant range only, as do observations of limb sounders
    (see `collect_limb_soundings`) and observations without an `instrument`
    column (such as those of `collect_limb_observations`), unless `instrument`
    gives one.

    The observations may be of any kind of platform (satellites, aircraft,
    fixed surface platforms, or radar stations), each named in the column of
    its kind: those of points (see `collect_observations` and
    `collect_footprint_observations`) and of features (see
    `collect_feature_observations`), whose target is the centroid of the
    viewed part, each at its time (see `epoch` or `time`) and elevation (its
    z-coordinate, as of the footprints that viewed a feature, or zero).

    With `instrument`, the samples are those of another instrument at the
    same observations, as of a sensor whose pixels differ from the swath
    that was observed (for example, a proposed sensor flown on the same
    platforms): one instrument for every observation, or one for the
    observations of each observing instrument (by name), the others sized
    as their own.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the
                platform and `instrument` that observe each.
        platforms (AllPlatforms | list[AllPlatforms]): The observing
                platforms, including those named by the observations.
        instrument (AllInstruments | Mapping[str, AllInstruments] | None):
                The instrument whose samples to size in place of the
                observing instrument, or one for each observing instrument's
                name, or `None` for the observing instruments.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sample sizes,
            spacings, azimuths, slant ranges, and scan angles.
    """
    lookup = _get_platform_lookup(platforms)
    result = observations.copy()
    for column in _COLUMNS:
        result[column] = np.nan
    if len(observations) == 0:
        return result
    times = _get_observation_times(observations).to_numpy()
    # the targets: observed points, or the centroids of observed regions
    geometry = observations.geometry.to_numpy()
    is_point = shapely.get_type_id(geometry) == shapely.GeometryType.POINT
    centroid = shapely.centroid(geometry)
    longitude, latitude = shapely.get_x(centroid), shapely.get_y(centroid)
    # at the elevation of the target: its z-coordinate (as of the footprints
    # that observed a region, see `collect_feature_observations`), or zero
    elevation = np.nan_to_num(
        np.where(
            is_point,
            shapely.get_z(np.where(is_point, geometry, None)),
            _get_footprint_elevations(geometry),
        )
    )
    keys = _get_platforms(observations)
    _check_platform_lookup(lookup, keys)
    names = _get_instruments(observations)
    position = np.arange(len(observations))
    columns = {column: result.columns.get_loc(column) for column in _COLUMNS}
    for (key, name), index in (
        pd.DataFrame({"key": keys, "name": names})
        .groupby(["key", "name"], sort=False, dropna=False)
        .indices.items()
    ):
        platform = lookup[key]
        index = position[index]
        track = platform.get_track(list(pd.DatetimeIndex(times[index], tz="UTC")))
        target = wgs84.latlon(latitude[index], longitude[index], elevation[index])
        result.iloc[index, columns["platform_range"]] = _get_viewing_geometry(
            track, target
        )["platform_range"]
        sized = _get_sized_instrument(platform, name, instrument)
        if sized is None or isinstance(sized, (LimbSounder, OccultationReceiver)):
            # without instruments, or of sounders, which view no targets
            continue
        velocity_frame = (
            sized.velocity_frame
            if isinstance(sized, (PointedInstrument, ConicalInstrument))
            else VelocityFrame.EARTH_FIXED
        )
        roll, _ = compute_view_angles(
            track,
            target,
            velocity_frame,
            sized.nadir_reference,
            attitude=sized.platform_attitude,
        )
        result.iloc[index, columns["scan_angle"]] = np.atleast_1d(roll)
        if isinstance(sized, (PointedInstrument, ConicalInstrument)):
            for prefix, method in (
                ("sample", sized.compute_sample_size),
                ("spacing", sized.compute_sample_spacing),
            ):
                along, cross = method(track, target)
                result.iloc[index, columns[f"{prefix}_along"]] = np.atleast_1d(along)
                result.iloc[index, columns[f"{prefix}_cross"]] = np.atleast_1d(cross)
            result.iloc[index, columns["sample_azimuth"]] = np.atleast_1d(
                sized.compute_sample_azimuth(track, target)
            )
    return result
